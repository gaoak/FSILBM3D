cat src/ConstParams.f90.f90  src/SolidSolver.f90  src/Solidbody.f90  src/FlowDomian.f90  src/FlowCondition.f90  src/Util.f90  src/main.f90  > combineforcheck.f90
gfortran -ffpe-trap=invalid,zero -fbacktrace -Wall -Wextra -pedantic -Warray-bounds -fbacktrace  -fbounds-check -Wconversion -Wconversion-extra -ffree-form -ffree-line-length-none -fopenmp -fdefault-real-8 -fdefault-double-8 -fimplicit-none -finit-real=zero -std=f2003 combineforcheck.f90 src/forkthread.o -g -g3 -Og -o combineforcheck
