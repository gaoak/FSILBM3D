# SPDX-License-Identifier: GPL-3.0-or-later
#
# FSILBM3D
# Copyright (C) 2025-2026 Ankang Gao and contributors

CMP = gfortran# ifort,ifx,gfortran

BUILD ?=
MPI ?= 0
MPIFC ?= mpifort

#######CMP settings###########
ifeq ($(CMP),ifort)
FC = ifort
FFLAGS = -diag-disable=10448 -fpp -O3 -free -qopenmp -heap-arrays #-real-size 32 -double-size 64
else ifeq ($(CMP),ifx)
FC = ifx
FFLAGS = -diag-disable=10448 -fpp -O3 -free -qopenmp -heap-arrays #-real-size 32 -double-size 64
else ifeq ($(CMP),gfortran)
FC = gfortran
FFLAGS = -cpp -O3
ifeq ($(BUILD),debug)
FFLAGS = -cpp -g -O0
FFLAGS += -ffpe-trap=invalid,zero -fbacktrace -Wall -Wextra -pedantic -Warray-bounds -fbacktrace  -fbounds-check
endif
FFLAGS += -Wconversion -Wconversion-extra -ffree-form -ffree-line-length-none -fopenmp -fimplicit-none -finit-real=zero -std=f2003 #-fdefault-real-4 -fdefault-double-8
endif
CC = cc

SRCDIR = ./src
ifeq ($(MPI),1)
FC = $(MPIFC)
FFLAGS += -DFSILBM3D_MPI
TARGET = FSILBM3D_mpi
else
TARGET = FSILBM3D
endif
OBJDIR = build-$(CMP)-mpi$(MPI)$(if $(BUILD),-$(BUILD))

### List of files for the main code
OBJ = $(addprefix $(OBJDIR)/,ConstParams.o ParallelDomain.o FlowCondition.o SolidSolver.o ImplicitIBM.o Solidbody.o FluidDomain.o LBMBlockComm.o Util.o main.o forkthread.o)

#-----------------------------------------------------------------------
# Normally no need to change anything below

.PHONY: all clean FORCE
all: $(TARGET)

$(TARGET): $(OBJ) FORCE
	$(FC) -o $@ $(FFLAGS) $(OBJ)

$(OBJDIR):
	mkdir -p $@

$(OBJDIR)/%.o: $(SRCDIR)/%.f90 | $(OBJDIR)
	cd $(OBJDIR) && $(FC) $(FFLAGS) $(INC) -c ../$< -o $(notdir $@)

$(OBJDIR)/%.o: $(SRCDIR)/%.c | $(OBJDIR)
	$(CC) -c -O3 $< -o $@

$(OBJDIR)/FlowCondition.o $(OBJDIR)/SolidSolver.o: $(OBJDIR)/ConstParams.o
$(OBJDIR)/ImplicitIBM.o: $(OBJDIR)/ParallelDomain.o
$(OBJDIR)/Solidbody.o: $(OBJDIR)/SolidSolver.o $(OBJDIR)/FlowCondition.o $(OBJDIR)/ParallelDomain.o $(OBJDIR)/ImplicitIBM.o
$(OBJDIR)/FluidDomain.o: $(OBJDIR)/Solidbody.o $(OBJDIR)/FlowCondition.o $(OBJDIR)/ParallelDomain.o
$(OBJDIR)/LBMBlockComm.o: $(OBJDIR)/FluidDomain.o
$(OBJDIR)/main.o $(OBJDIR)/Util.o: $(OBJDIR)/LBMBlockComm.o

clean:
	rm -rf $(OBJDIR)
	rm -f $(SRCDIR)/*.o $(SRCDIR)/*.mod $(SRCDIR)/*.smod
	rm -f *.o *.mod *.smod $(TARGET)
