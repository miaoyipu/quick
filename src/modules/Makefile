#!/bin/sh

#include "../config.h"

# --- Makefile for Quick Program ---
#				- v 2.00 2010/10/25 Yipu Miao
#				- v 1.18 2009/09/16 John Faver Exp $ 
#				- Makefile created by mkmf.pl $Id:
#	--------
#	 INDEX
#	--------
#	A. Compiler Setting			! Intel Fortran 9.0+ or GNU Fortran is recommended for single CPU Version
#								! mpif90 is recommended for MPI Multi-CPU Version
#	B. Make Object Files		! Source files --> Object files
#	C. Make Executed files		! Object files --> Executed files
#	D. Self-defined Option		! Make option

#************************************************************************
#                  A. Compiler Settings
# 
#   FC specifies f90 compiler
#   FFLAGS are compliation options
#   LFLAGS are linking flags
#
#************************************************************************

# Tested on LINUX machines and Mac machines
#
#---------------------
# a. MPI Version
#---------------------
#FC = mpif90
#FC = /usr/local/bin/mpif90
#MPI_INCLUDES = /usr/local/include
#MPI_LIBS = /usr/local/lib

#---------------------
# b. Single CPU Version
#---------------------
FC=gfortran
#FC=ifort

#---------------------
# other compiler options
#---------------------
#FFLAGS = -g -O3 -xW -pg -traceback
#FFLAGS = -O3 -xW -ipo
#FFLAGS = -O3
#FFLAGS = -i4 -O3 -auto -assume byterecl -w95 -cm 
#FFLAGS = -g -O3 -traceback
FFLAGS = -std=gnu -W
LD = $(FC)
LDFLAGS = $(FFLAGS)
#LDFLAGS = $(FFLAGS)
#LDFLAGS = $(FFLAGS) -static -L/opt/intel/ict/2.0/cmkl/8.0.1/lib/32
TMPFILES = *.mod *.stb

# CPP Compiler
FPP = cpp -traditional -P  -DBINTRAJ 
CPP = gcc -O3 -DCUBLAS_USE_THUNKING

# CUDA Compiler
CUDAC = /usr/local/cuda/bin/nvcc
CUDA_LIBPATH = -L/usr/local/gfortran/lib 
CUDA_FLAGS= -m64 -O3

# G++ Compiler
CXX = g++
CFLAGS = -L/usr/local/cuda/lib -lcuda -lcudart -lm -lgfortran $(CUDA_LIBPATH)    

#----------------------
# exe file location
#----------------------
exefile = ../bin/quick

#----------------------
# library file location
#----------------------
libfolder = ../lib

#----------------------
# cuda files
#----------------------
cudafolder = ./cuda
cudaobj = gpu_info.o gpu_write_info.o gpu_startup.o

#----------------------
# quick modules and object files
#----------------------
modfolder = ./modules

mpimod  = quick_mpi_module.f90
mpimod_o= quick_mpi_module.o

modfile0 = quick_constants_module.f90
modobj0  = quick_constants_module.o
modfile1 = quick_method_module.f90 quick_molspec_module.f90 quick_gaussian_class_module.f90 
modobj1  = quick_method_module.o quick_molspec_module.o quick_gaussian_class_module.o
modfile2 = quick_size_module.f90 quick_amber_interface_module.f90 quick_basis_module.f90 \
		quick_calculated_module.f90 quick_divcon_module.f90 \
		quick_ecp_module.f90 quick_electrondensity_module.f90 quick_files_module.f90 \
		quick_gridpoints_module.f90 \
		quick_mfcc_module.f90 quick_params_module.f90 quick_pb_module.f90 \
		quick_scratch_module.f90 quick_timer_module.f90 quick_all_module.f90
modobj2  = quick_size_module.o quick_amber_interface_module.o quick_basis_module.o \
		quick_calculated_module.o quick_divcon_module.o \
		quick_ecp_module.o quick_electrondensity_module.o quick_files_module.o \
		quick_gridpoints_module.o \
		quick_mfcc_module.o quick_params_module.o quick_pb_module.o \
		quick_scratch_module.o quick_timer_module.o quick_all_module.o
#.DEFAULT:
#	-touch $@

all: quick_modules

#************************************************************************
# 
#                 B. Make Object Files
# 
#************************************************************************

#================= quick module library =================================
quick_modules:
	$(FPP) quick_basis_module_.f90 > quick_basis_module.f90
	$(FPP) quick_method_module_.f90 > quick_method_module.f90
	$(FPP) quick_molspec_module_.f90 > quick_molspec_module.f90
	$(FPP) quick_calculated_module_.f90 > quick_calculated_module.f90
	$(FPP) quick_timer_module_.f90 > quick_timer_module.f90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c  ${mpimod}
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c  ${modfile0}
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c  ${modfile1}
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c  ${modfile2}
		
#************************************************************************
# 
#                 C. Self-defined Option
# 
#************************************************************************

# - 1. Clean object files
clean: neat
	-rm -f .cppdefs *.mod *.o subs/*.o quick
	
neat:
	-rm -f $(TMPFILES)

# - 2. Make tags for source files
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)