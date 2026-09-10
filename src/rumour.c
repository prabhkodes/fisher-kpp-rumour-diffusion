#include <sys/stat.h>   /* mkdir() for the output directory */
#include <petscdmda.h>
#include <petscts.h>

typedef struct {
  PetscReal D;      /* Diffusion rate (how fast people talk) */
  PetscReal alpha;  /* Persuasion rate (how "sticky" the rumor is) */
} AppCtx;

/* RHS Function: du/dt = F(t, u) */
static PetscErrorCode RHSFunction(TS ts, PetscReal t, Vec U, Vec F, void *ptr) {
  AppCtx            *ctx = (AppCtx*)ptr;
  DM                dm;
  PetscInt          i, j, xs, ys, xm, ym, mx, my;
  PetscScalar       **u, **f;
  PetscReal         hx, hy, hx2inv, hy2inv;
  Vec               Ulocal;

  PetscFunctionBeginUser;
  PetscCall(TSGetDM(ts, &dm));
  PetscCall(DMDAGetInfo(dm, 0, &mx, &my, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0));
  
  hx = 1.0 / (PetscReal)(mx - 1);
  hy = 1.0 / (PetscReal)(my - 1);
  hx2inv = 1.0 / (hx * hx);
  hy2inv = 1.0 / (hy * hy);

  /* scatter the global vector to local vectors + ghosts to calculate Laplacian */
  PetscCall(DMGetLocalVector(dm, &Ulocal));
  PetscCall(DMGlobalToLocalBegin(dm, U, INSERT_VALUES, Ulocal));
  PetscCall(DMGlobalToLocalEnd(dm, U, INSERT_VALUES, Ulocal));

  PetscCall(DMDAVecGetArrayRead(dm, Ulocal, &u));
  PetscCall(DMDAVecGetArray(dm, F, &f));
  PetscCall(DMDAGetCorners(dm, &xs, &ys, NULL, &xm, &ym, NULL));

  for (j = ys; j < ys + ym; j++) {
    for (i = xs; i < xs + xm; i++) {
      if (i == 0 || j == 0 || i == mx - 1 || j == my - 1) {
        f[j][i] = 0.0; /* boundary */
      } else {
        /* diffusion */
        PetscScalar lap = (u[j][i-1] - 2.0*u[j][i] + u[j][i+1]) * hx2inv +
                          (u[j-1][i] - 2.0*u[j][i] + u[j+1][i]) * hy2inv;
        
        /* logistic growth (social saturation) */
        PetscScalar growth = ctx->alpha * u[j][i] * (1.0 - u[j][i]);

        f[j][i] = ctx->D * lap + growth;
      }
    }
  }

  PetscCall(DMDAVecRestoreArrayRead(dm, Ulocal, &u));
  PetscCall(DMDAVecRestoreArray(dm, F, &f));
  PetscCall(DMRestoreLocalVector(dm, &Ulocal));
  PetscFunctionReturn(PETSC_SUCCESS);
}


static PetscErrorCode MyVTKMonitor(TS ts, PetscInt step, PetscReal t, Vec U, void *ctx) {
  char           filename[PETSC_MAX_PATH_LEN];
  PetscViewer    viewer;
  PetscMPIInt    rank;

  PetscFunctionBeginUser;
  PetscCallMPI(MPI_Comm_rank(PETSC_COMM_WORLD, &rank));

  if (step % 5 == 0) {
    if (rank == 0) {
      mkdir("files", 0777); 
    }
    PetscCallMPI(MPI_Barrier(PETSC_COMM_WORLD));

    PetscCall(PetscSNPrintf(filename, sizeof(filename), "files/rumour_%03d.vts", (int)step));

    PetscCall(PetscViewerVTKOpen(PETSC_COMM_WORLD, filename, FILE_MODE_WRITE, &viewer));
    PetscCall(VecView(U, viewer));
    PetscCall(PetscViewerDestroy(&viewer));

    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "Saved state to %s at time %g\n", filename, (double)t));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

int main(int argc, char **argv) {
  DM        dm;
  TS        ts;
  Vec       U;
  AppCtx    ctx;
  PetscInt  mx = 50, my = 50;

  PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));

  /* Parameters: D=Spread speed || alpha=viral strength */
  ctx.D = 0.001; ctx.alpha = 0.5;
  PetscCall(PetscOptionsGetReal(NULL, NULL, "-D", &ctx.D, NULL));
  PetscCall(PetscOptionsGetReal(NULL, NULL, "-alpha", &ctx.alpha, NULL));




  /* make grid */
  PetscCall(DMDACreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_NONE, DM_BOUNDARY_NONE,
                         DMDA_STENCIL_STAR, mx, my, PETSC_DECIDE, PETSC_DECIDE,
                         1, 1, NULL, NULL, &dm));
  PetscCall(DMSetFromOptions(dm));
  PetscCall(DMSetUp(dm));
  PetscCall(DMDASetFieldName(dm, 0, "Awareness"));

  /* put "Patient Zero" in the center/centre lol */
  PetscCall(DMCreateGlobalVector(dm, &U));
  PetscCall(VecSet(U, 0.0));
  
  /* initial condition: center point knows the rumor */
  PetscInt center_i = mx/2, center_j = my/2;
  PetscScalar **u_init;
  PetscInt xs, ys, xm, ym;
  PetscCall(DMDAGetCorners(dm, &xs, &ys, NULL, &xm, &ym, NULL));
  PetscCall(DMDAVecGetArray(dm, U, &u_init));
  if (center_i >= xs && center_i < xs+xm && center_j >= ys && center_j < ys+ym) {
    u_init[center_j][center_i] = 1.0;
  }
  PetscCall(DMDAVecRestoreArray(dm, U, &u_init));

  /* make time stepperrr */
  PetscCall(TSCreate(PETSC_COMM_WORLD, &ts));
  PetscCall(TSSetDM(ts, dm));
  PetscCall(TSSetProblemType(ts, TS_NONLINEAR));
  PetscCall(TSSetRHSFunction(ts, NULL, RHSFunction, &ctx));

  /* solver defaults */
  PetscCall(TSSetType(ts, TSRK));
  PetscCall(TSSetMaxTime(ts, 20.0));
  PetscCall(TSSetExactFinalTime(ts, TS_EXACTFINALTIME_STEPOVER));
  PetscCall(TSSetFromOptions(ts));

  PetscCall(TSMonitorSet(ts, MyVTKMonitor, NULL, NULL));
  /* run */
  PetscCall(TSSolve(ts, U));

  /* view */
  PetscCall(VecViewFromOptions(U, NULL, "-view_rumor"));

  /* cleanup */
  PetscCall(VecDestroy(&U));
  PetscCall(TSDestroy(&ts));
  PetscCall(DMDestroy(&dm));
  PetscCall(PetscFinalize());
  return 0;
}