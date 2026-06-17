# SPDX-License-Identifier: GPL-3.0-or-later
#
# FSILBM3D
# Copyright (C) 2025-2026 Ankang Gao and contributors

CMP = gfortran# ifort,ifx,gfortran

BUILD ?=

#######CMP settings###########
ifeq ($(CMP),ifort)
FC = ifort
FFLAGS = -diag-disable=10448 -fpp -O3 -free -qopenmp -heap-arrays #-real-size 32 -double-size 64
else ifeq ($(CMP),ifx)
FC = ifx
FFLAGS = -diag-disable=10448 -fpp -O3 -free -qopenmp -heap-arrays #-real-size 32 -double-size 64
else ifeq ($(CMP),gfortran)
FC = gfortran
FFLAGS = -O3
ifeq ($(BUILD),debug)
FFLAGS = -cpp -g -O0
FFLAGS += -ffpe-trap=invalid,zero -fbacktrace -Wall -Wextra -pedantic -Warray-bounds -fbacktrace  -fbounds-check
endif
FFLAGS += -Wconversion -Wconversion-extra -ffree-form -ffree-line-length-none -fopenmp -fimplicit-none -finit-real=zero -std=f2003 #-fdefault-real-4 -fdefault-double-8
endif
CC = cc
CPP = c++

SRCDIR = ./src

### List of files for the main code
OBJDECOMP = $(SRCDECOMP:%.f90=%.o)
SRC = $(SRCDIR)/ConstParams.f90 $(SRCDIR)/FlowCondition.f90 $(SRCDIR)/SolidSolver.f90 $(SRCDIR)/Solidbody.f90 $(SRCDIR)/FluidDomain.f90 $(SRCDIR)/LBMBlockComm.f90 $(SRCDIR)/Util.f90 $(SRCDIR)/main.f90
CSRC = $(SRCDIR)/forkthread.c
CPPSRC = $(SRCDIR)/velocitymap.cpp
OBJ = $(SRC:%.f90=%.o)
COBJ = $(CSRC:%.c=%.o)
CPPOBJ = $(CPPSRC:%.cpp=%.o)

#######OPTIONS settings###########
OPT := -I$(SRCDIR) $(FFLAGS)
LINKOPT := $(FFLAGS)


#-----------------------------------------------------------------------
# Normally no need to change anything below

all: FSILBM3D

FSILBM3D : $(OBJ) $(COBJ) #$(CPPOBJ)
	$(FC) -o $@ $(LINKOPT) $(OBJDECOMP) $(OBJ) $(COBJ)
#$(CPPOBJ)  -lstdc++


$(OBJ):$(SRCDIR)%.o : $(SRCDIR)%.f90
	$(FC) $(FFLAGS) $(OPT) $(INC) -c $<
	mv $(@F) ${SRCDIR}

$(COBJ):$(SRCDIR)%.o : $(SRCDIR)%.c
	$(CC)  -c -O3 $<
	mv $(@F) ${SRCDIR}

$(CPPOBJ):$(SRCDIR)%.o : $(SRCDIR)%.cpp
	$(CPP)  -c -O3 $<
	mv $(@F) ${SRCDIR}

clean:
	rm -f $(SRCDIR)/*.o $(SRCDIR)/*.mod $(SRCDIR)/*.smod
	rm -f *.o *.mod *.smod FSILBM3D
