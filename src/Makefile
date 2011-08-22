#!/bin/sh

#include "config.h"

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

all: quick

#************************************************************************
# 
#                 B. Make Object Files
# 
#************************************************************************

#================= common subroutine library ============================
quick_subs:
	cp *.mod ./subs
	cd ./subs && $(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c *.f90
	ar -r $(libfolder)/quicklib.a ./subs/*.o
	rm ./subs/*.mod

#================= quick module library =================================
quick_modules:
	cd $(modfolder) && make
	mv ${modfolder}/*.mod ${modfolder}/*.o ./
#=========== targets for cuda =========================================
quick_cuda:
	cd $(cudafolder) && $(CUDAC) $(CUDA_FLAGS) -c *.cu
	cd $(cudafolder) && $(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c *.f90
	mv $(cudafolder)/*.o ./
	
#================= quick core subroutines ===============================
main.o: main.f90 quick_modules
	$(FPP) main.f90 > _main.f90
	$(FC) -o main.o $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       _main.f90
initialize.o: initialize.f90 quick_modules
	$(FPP) initialize.f90 > _initialize.f90
	$(FC) -o initialize.o $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       _initialize.f90
read_job_and_atom.o: read_job_and_atom.f90 quick_modules
	$(FPP) read_job_and_atom.f90 > _read_job_and_atom.f90
	$(FC) -o read_job_and_atom.o $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       _read_job_and_atom.f90
shell.o: shell.f90 quick_modules
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       shell.f90
getMolSad.o: getMolSad.f90 quick_modules
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       getMolSad.f90
getMol.o: getMol.f90 quick_modules
	$(FPP) getMol.f90 > _getMol.f90
	$(FC) -o getMol.o $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       _getMol.f90
schwarz.o: schwarz.f90 quick_modules
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       schwarz.f90
quick_one_electron_integral.o: quick_one_electron_integral.f90 quick_modules
	$(FPP) quick_one_electron_integral.f90 > _quick_one_electron_integral.f90
	$(FC) -o quick_one_electron_integral.o $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       _quick_one_electron_integral.f90
getEnergy.o: getEnergy.f90 quick_modules
	$(FPP) getEnergy.f90 > _getEnergy.f90
	$(FC) -o getEnergy.o $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       _getEnergy.f90
inidivcon.o: inidivcon.f90 quick_modules
	$(FPP) inidivcon.f90 > _inidivcon.f90
	$(FC) -o inidivcon.o $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       _inidivcon.f90
ecp.o: ecp.f90 quick_modules
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       ecp.f90
uelectdii.o: uelectdii.f90 quick_modules
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       uelectdii.f90
hfoperator.o: hfoperator.f90 quick_modules
	$(FPP) hfoperator.f90 >_hfoperator.f90
	$(FC) -o hfoperator.o $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       _hfoperator.f90
nuclear.o: nuclear.f90 quick_modules
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       nuclear.f90
fmm.o: fmm.f90 quick_modules
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       fmm.f90
dft.o: dft.f90 quick_modules
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       dft.f90
sedftoperator.o: sedftoperator.f90 quick_modules
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       sedftoperator.f90
scf.o: scf.f90 quick_modules
	$(FPP) scf.f90 > _scf.f90
	$(FC) -o scf.o $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c		      _scf.f90
uscf.o: uscf.f90 quick_modules
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       uscf.f90
finalize.o: finalize.f90 quick_modules
	$(FPP) finalize.f90 > _finalize.f90
	$(FC) -o finalize.o $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       _finalize.f90
uhfoperator.o: uhfoperator.f90 quick_modules
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       uhfoperator.f90
udftoperator.o: udftoperator.f90 quick_modules
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       udftoperator.f90
usedftoperator.o: usedftoperator.f90 quick_modules
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       usedftoperator.f90
quick_debug.o: quick_debug.f90 quick_modules
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       quick_debug.f90
mpi_setup.o: mpi_setup.f90 quick_modules
	$(FPP) mpi_setup.f90 > _mpi_setup.f90
	$(FC) -o mpi_setup.o $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c      _mpi_setup.f90
calMP2.o: calMP2.f90 quick_modules
	$(FPP) calMP2.f90 > _calMP2.f90
	$(FC) -o calMP2.o $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c         _calMP2.f90
optimize.o: optimize.f90 quick_modules
	$(FPP) optimize.f90 > _optimize.f90
	$(FC) -o optimize.o $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       _optimize.f90
gradient.o: gradient.f90 quick_modules
	$(FPP) gradient.f90 > _gradient.f90
	$(FC) -o gradient.o $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       _gradient.f90
hessian.o: hessian.f90 quick_modules
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       hessian.f90
CPHF.o: CPHF.f90 quick_modules
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       CPHF.f90
frequency.o: frequency.f90 quick_modules
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       frequency.f90
MFCC.o: MFCC.f90 quick_modules
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       MFCC.f90
basis.o: basis.f90 quick_modules
	$(FPP) basis.f90 > _basis.f90
	$(FC) -o basis.o $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c          _basis.f90

#=========== targets for amber-quick interface ========================
# This is a fake amber-quick interface
fake_amber_interface.o: fake_amber_interface.f90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       fake_amber_interface.f90

# This is for amber-quick interface
amber_interface.o: amber_interface.f90 quick_modules qmmm_module.mod
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       amber_interface.f90

									
OBJ =   $(mpimod_o) main.o ${modobj0} ${modobj1} ${modobj2} \
		initialize.o read_job_and_atom.o fmm.o \
		getMolSad.o getMol.o shell.o schwarz.o quick_one_electron_integral.o \
		getEnergy.o inidivcon.o ecp.o hfoperator.o nuclear.o \
		dft.o sedftoperator.o \
		scf.o uscf.o finalize.o uhfoperator.o udftoperator.o usedftoperator.o \
		uelectdii.o mpi_setup.o quick_debug.o calMP2.o optimize.o \
		gradient.o hessian.o CPHF.o frequency.o MFCC.o basis.o 

#************************************************************************
# 
#                 C. Make Executed Files
# 
#************************************************************************

quick: quick_modules quick_subs $(OBJ) fake_amber_interface.o
	$(FC) -o ${exefile} $(OBJ) $(libfolder)/quicklib.a fake_amber_interface.o $(LDFLAGS)

quick.cuda: quick_cuda quick_modules quick_subs $(OBJ) fake_amber_interface.o
	$(CXX) -o ${exefile}.cuda $(cudaobj) $(OBJ) $(libfolder)/quicklib.a fake_amber_interface.o $(CFLAGS) 


quick_lib:$(OBJ) ambermod amber_interface.o

ambermod:
	cd ../../../AmberTools/src/sqm && $(MAKE) qmmm_module.o
	cp ../../../AmberTools/src/sqm/qmmm_module.mod .
	cp ../../../AmberTools/src/sqm/qmmm_vsolv_module.mod .
	cp ../../../AmberTools/src/sqm/qmmm_struct_module.mod .
	cp ../../../AmberTools/src/sqm/qmmm_nml_module.mod .
	cp ../../../AmberTools/src/sqm/qmmm_module.o .

		
#************************************************************************
# 
#                 D. Self-defined Option
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
