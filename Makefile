CMP = intel# intel,gcc

BUILD ?=

#######CMP settings###########
ifeq ($(CMP),intel)
FC = ifort
FFLAGS = -fpp -O3 -free -qopenmp -heap-arrays #-real-size 32 -double-size 64
else ifeq ($(CMP),gcc)
FC = gfortran
FFLAGS = -O3
ifeq ($(BUILD),debug)
FFLAGS = -cpp -g3 -Og
FFLAGS += -ffpe-trap=invalid,zero -fbacktrace -Wall -Wextra -pedantic -Warray-bounds -fbacktrace  -fbounds-check
endif
FFLAGS += -Wconversion -Wconversion-extra -ffree-form -ffree-line-length-none -fopenmp -fimplicit-none -finit-real=zero -std=f2003 #-fdefault-real-4 -fdefault-double-8
endif
CC = gcc

SRCDIR = ./src

### List of files for the main code
OBJDECOMP = $(SRCDECOMP:%.f90=%.o)
SRC = $(SRCDIR)/Modules.f90 $(SRCDIR)/SolidSolver.f90 $(SRCDIR)/Solidbody.f90 $(SRCDIR)/LatticeBoltzmannSolver.f90 (SRCDIR)/PostProcessing.f90 $$(SRCDIR)/Initialization.f90 $(SRCDIR)/main.f90  $(SRCDIR)/Util.f90
CSRC = $(SRCDIR)/forkthread.c
OBJ = $(SRC:%.f90=%.o)
COBJ = $(CSRC:%.c=%.o)

#######OPTIONS settings###########
OPT := -I$(SRCDIR) $(FFLAGS)
LINKOPT := $(FFLAGS)


#-----------------------------------------------------------------------
# Normally no need to change anything below

all: FSILBM3D

FSILBM3D : $(OBJ) $(COBJ)
	$(FC) -o $@ $(LINKOPT) $(OBJDECOMP) $(OBJ) $(COBJ)


$(OBJ):$(SRCDIR)%.o : $(SRCDIR)%.f90
	$(FC) $(FFLAGS) $(OPT) $(INC) -c $<
	mv $(@F) ${SRCDIR}

$(COBJ):$(SRCDIR)%.o : $(SRCDIR)%.c
	$(CC)  -c $<
	mv $(@F) ${SRCDIR}

clean:
	rm -f $(SRCDIR)/*.o $(SRCDIR)/*.mod $(SRCDIR)/*.smod
	rm -f *.o *.mod *.smod FSILBM3D