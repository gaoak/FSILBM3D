# SPDX-License-Identifier: GPL-3.0-or-later
#
# FSILBM3D
# Copyright (C) 2025-2026 Ankang Gao and contributors

CMP = gfortran# nvfortran,ifort,ifx,gfortran

BUILD ?=
# CUDA is the default build: plain "make" compiles every source file,
# including src/cuda_lbm.cuf. Use "make GPU=0" for the CPU-only version.
GPU ?= 1

# STRICT controls FMA in nvfortran builds; it is not a solver tolerance.
# CPU baseline: make clean && make GPU=0 (compiler selected by CMP).
# STRICT=0 and STRICT=1 give identical flags for gfortran/ifort/ifx.
# The gfortran path uses generic x86-64 flags that do not enable FMA;
# this does not describe ifort/ifx, whose FMA behavior depends on their flags.
# Match the master's compiler, flags, inputs and thread count for CPU regressions.
# GPU comparison build: make clean && make GPU=1 STRICT=1
# This adds -Mnofma for nvfortran host code and nofma for device code.
# Defaults are GPU=1 and STRICT=0, so plain "make" builds the GPU version
# with FMA allowed (unless overridden on the command line or in the environment).
# For a clean default GPU build: make clean && make
# Explicit equivalent: make clean && make GPU=1 STRICT=0
# FMA changes rounding and can affect sensitive FSI results; disabling it
# reduces one source of CPU/GPU differences, but does not ensure bitwise equality.
# STRICT also affects GPU=0 CMP=nvfortran; it does not change gfortran/ifort/ifx flags.
STRICT ?= 0
GPUARCH = cc80
FMAOPT =
ifeq ($(STRICT),1)
GPUARCH = cc80,nofma
FMAOPT = -Mnofma
endif

#######CMP settings###########
ifeq ($(CMP),nvfortran)
FC = nvfortran
FFLAGS = -O3 -Mfree -Mpreprocess -mp -cuda -gpu=$(GPUARCH) $(FMAOPT)
else ifeq ($(CMP),ifort)
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
FFLAGS += -cpp -Wconversion -Wconversion-extra -ffree-form -ffree-line-length-none -fopenmp -fimplicit-none -finit-real=zero -std=f2003 #-fdefault-real-4 -fdefault-double-8
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

ifeq ($(GPU),1)
FC = nvfortran
FFLAGS = -O3 -Mfree -Mpreprocess -mp -cuda -gpu=$(GPUARCH) $(FMAOPT)
FFLAGS += -DUSE_CUDA
CUDAOBJ = $(SRCDIR)/cuda_lbm.o
else
CUDAOBJ =
endif

#######OPTIONS settings###########
OPT := -I$(SRCDIR) $(FFLAGS)
LINKOPT := $(FFLAGS)


#-----------------------------------------------------------------------
# Normally no need to change anything below

all: FSILBM3D

FSILBM3D : $(OBJ) $(COBJ) $(CUDAOBJ) #$(CPPOBJ)
	$(FC) -o $@ $(LINKOPT) $(OBJDECOMP) $(OBJ) $(COBJ) $(CUDAOBJ)
#$(CPPOBJ)  -lstdc++


$(OBJ):$(SRCDIR)%.o : $(SRCDIR)%.f90
	$(FC) $(FFLAGS) $(OPT) $(INC) -c $<
	mv $(@F) ${SRCDIR}

$(COBJ):$(SRCDIR)%.o : $(SRCDIR)%.c
	$(CC)  -c -O3 $<
	mv $(@F) ${SRCDIR}

$(CUDAOBJ): $(SRCDIR)/cuda_lbm.cuf $(SRCDIR)/ConstParams.o
	$(FC) $(FFLAGS) -c $< -o $@

ifeq ($(GPU),1)
# these modules import cuda_lbm_backend, so the backend module must exist first
$(SRCDIR)/Solidbody.o: $(CUDAOBJ)
$(SRCDIR)/FluidDomain.o: $(CUDAOBJ)
$(SRCDIR)/LBMBlockComm.o: $(CUDAOBJ)
endif

$(CPPOBJ):$(SRCDIR)%.o : $(SRCDIR)%.cpp
	$(CPP)  -c -O3 $<
	mv $(@F) ${SRCDIR}

clean:
	rm -f $(SRCDIR)/*.o $(SRCDIR)/*.mod $(SRCDIR)/*.smod
	rm -f *.o *.mod *.smod FSILBM3D
