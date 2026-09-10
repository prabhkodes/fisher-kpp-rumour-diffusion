# Build with PETSc's own makefile machinery so the include and link flags
# always match how PETSc was configured.
#   export PETSC_DIR=/path/to/petsc
#   export PETSC_ARCH=your-arch
#   make

include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules

rumour: src/rumour.c
	$(CC) -o rumour src/rumour.c $(PETSC_CC_INCLUDES) $(PETSC_CC_LIB)

clean::
	$(RM) rumour
	$(RM) -r files

.PHONY: clean
