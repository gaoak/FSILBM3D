CMP = gcc# intel,gcc

BUILD ?=

#######CMP settings###########
ifeq ($(CMP),intel)
FC = ifort
FFLAGS = -fpp -O3 -free -qopenmp #-real-size 32 -double-size 64
else ifeq ($(CMP),gcc)
FC = gfortran
FFLAGS = -O3
ifeq ($(BUILD),debug)
FFLAGS = -cpp -g3 -Og
FFLAGS += -ffpe-trap=invalid,zero -fbacktrace -Wall -Wextra -pedantic -Warray-bounds -fbacktrace  -fbounds-check
endif
FFLAGS += -Wconversion -Wconversion-extra -ffree-form -ffree-line-length-none -fopenmp -fimplicit-none -finit-real=zero -std=f2003 #-fdefault-real-4 -fdefault-double-8
endif


SRCDIR = ./src

### List of files for the main code
OBJDECOMP = $(SRCDECOMP:%.f90=%.o)
SRC = $(SRCDIR)/Modules.f90 $(SRCDIR)/LatticeBoltzmannSolver.f90 $(SRCDIR)/PostProcessing.f90 $(SRCDIR)/StructureSolver.f90 $(SRCDIR)/Initialization.f90  $(SRCDIR)/Interaction.f90 $(SRCDIR)/main.f90  $(SRCDIR)/Util.f90
OBJ = $(SRC:%.f90=%.o)

#######OPTIONS settings###########
OPT := -I$(SRCDIR) $(FFLAGS)
LINKOPT := $(FFLAGS)


#-----------------------------------------------------------------------
# Normally no need to change anything below

all: FSILBM3D

FSILBM3D : $(OBJ)
	$(FC) -o $@ $(LINKOPT) $(OBJDECOMP) $(OBJ)


$(OBJ):$(SRCDIR)%.o : $(SRCDIR)%.f90
	$(FC) $(FFLAGS) $(OPT) $(INC) -c $<
	mv $(@F) ${SRCDIR}
	#mv *.mod ${SRCDIR}

## This %.o : %.f90 doesn't appear to be called...
%.o : %.f90
	$(FC) $(FFLAGS) $(INC) -c $<


clean:
	rm -f $(SRCDIR)/*.o $(SRCDIR)/*.mod $(SRCDIR)/*.smod
	rm -f *.o *.mod *.smod FSILBM3D