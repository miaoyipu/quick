#!/bin/sh

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
FC=/usr/local/bin/mpif90
mpimod  = quick_mpi_module.f90
mpimod_o= quick_mpi_module.o
#---------------------
# b. Single CPU Version
#---------------------
#FC=gfortran
#FC=ifort

#---------------------
# other compiler options
#---------------------
#FFLAGS = -g -O3 -xW -pg -traceback
#FFLAGS = -O3 -xW -ipo
#FFLAGS = -O3
#FFLAGS = -i4 -O3 -auto -assume byterecl -w95 -cm 
#FFLAGS = -g -O3 -traceback
LD = $(FC)
#LDFLAGS = $(FFLAGS) -static
#LDFLAGS = $(FFLAGS)
#LDFLAGS = $(FFLAGS) -static -L/opt/intel/ict/2.0/cmkl/8.0.1/lib/32
TMPFILES = *.mod *.stb


#----------------------
# exe file location
#----------------------
exefile = ../bin/quick

#----------------------
# library file location
#----------------------
libfolder = ../lib/


#----------------------
# quick modules and object files
#----------------------
modfolder = ./modules

modfile1 = quick_method_module.f90
modobj1  = quick_method_module.o
modfile2 = quick_size_module.f90 quick_amber_interface_module.f90 quick_basis_module.f90 \
		quick_calculated_module.f90 quick_constants_module.f90 quick_divcon_module.f90 \
		quick_ecp_module.f90 quick_electrondensity_module.f90 quick_files_module.f90 \
		quick_gaussian_class_module.f90 quick_geocnverg_module.f90 quick_gridpoints_module.f90 \
		quick_mfcc_module.f90 quick_molspec_module.f90 quick_params_module.f90 quick_pb_module.f90 \
		quick_scratch_module.f90 quick_timer_module.f90 quick_all_module.f90
modobj2  = quick_size_module.o quick_amber_interface_module.o quick_basis_module.o \
		quick_calculated_module.o quick_constants_module.o quick_divcon_module.o \
		quick_ecp_module.o quick_electrondensity_module.o quick_files_module.o \
		quick_gaussian_class_module.o quick_geocnverg_module.o quick_gridpoints_module.o \
		quick_mfcc_module.o quick_molspec_module.o quick_params_module.o quick_pb_module.o \
		quick_scratch_module.o quick_timer_module.o quick_all_module.o
#.DEFAULT:
#	-touch $@

all: quick

#************************************************************************
# 
#                 B. Make Object Files
# 
#************************************************************************


# if you want to do single CPU calculation, please uncomment following sentence
# FC=ifort
#mpimod  =mpi_fake.f90
#mpimod_o=mpi_fake.o
#	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c	mpi_fake.f90
#

quick_modules.o:
	cd ${modfolder} && $(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c  ${mpimod}
	cd ${modfolder} && $(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c  ${modfile1}
	cd ${modfolder} && $(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c  ${modfile2}
	mv ${modfolder}/*.mod ${modfolder}/*.o ./
	
main.o: main.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       main.f90
initialize.o: initialize.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       initialize.f90
readInput.o: readInput.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       readInput.f90
shell.o: shell.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       shell.f90
subs.o: subs.f90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       subs.f90
getMolSad.o: getMolSad.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       getMolSad.f90
getMol.o: getMol.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       getMol.f90
spdfgh.o: spdfgh.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       spdfgh.f90
hrr.o: hrr.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       hrr.f90
schwarz.o: schwarz.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       schwarz.f90
fullx.o: fullx.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       fullx.f90
getEnergy.o: getEnergy.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       getEnergy.f90
inidivcon.o: inidivcon.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       inidivcon.f90
ecp.o: ecp.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       ecp.f90
electdii.o: electdii.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       electdii.f90
uelectdii.o: uelectdii.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       uelectdii.f90
hfoperator.o: hfoperator.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       hfoperator.f90
nuclear.o: nuclear.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       nuclear.f90
fmm.o: fmm.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       fmm.f90
overlap.o: overlap.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       overlap.f90
dftoperator.o: dftoperator.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       dftoperator.f90
shelldft.o: shelldft.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       shelldft.f90
dft.o: dft.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       dft.f90
gridform.o: gridform.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       gridform.f90
pteval.o: pteval.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       pteval.f90
denspt.o: denspt.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       denspt.f90
sedftoperator.o: sedftoperator.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       sedftoperator.f90
repulsion.o: repulsion.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       repulsion.f90
attraction.o: attraction.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       attraction.f90
scf.o: scf.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       scf.f90
uscf.o: uscf.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       uscf.f90
finalize.o: finalize.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       finalize.f90
uhfoperator.o: uhfoperator.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       uhfoperator.f90
udftoperator.o: udftoperator.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       udftoperator.f90
usedftoperator.o: usedftoperator.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       usedftoperator.f90
debug_quick.o: debug_quick.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       debug_quick.f90
mpi_setup.o: mpi_setup.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       mpi_setup.f90
dipole.o: dipole.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       dipole.f90
zmake.o: zmake.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       zmake.f90
calMP2.o: calMP2.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       calMP2.f90
optimize.o: optimize.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       optimize.f90
gradient.o: gradient.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       gradient.f90
lbfgs.o: lbfgs.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       lbfgs.f90
hessian.o: hessian.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       hessian.f90
CPHF.o: CPHF.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       CPHF.f90
frequency.o: frequency.f90 quick_modules.o
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       frequency.f90

# This is a fake amber-quick interface
fake_amber_interface.o: fake_amber_interface.f90
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       fake_amber_interface.f90

# This is for amber-quick interface
amber_interface.o: amber_interface.f90 quick_modules.o qmmm_module.mod
	$(FC) $(CPPDEFS) $(CPPFLAGS) $(FFLAGS) -c       amber_interface.f90
									
OBJ = $(mpimod_o) main.o ${modobj1} ${modobj2} initialize.o readInput.o fmm.o overlap.o \
		subs.o getMolSad.o getMol.o shell.o spdfgh.o hrr.o schwarz.o fullx.o \
		getEnergy.o inidivcon.o ecp.o electdii.o hfoperator.o nuclear.o dftoperator.o \
		shelldft.o dft.o gridform.o pteval.o denspt.o sedftoperator.o repulsion.o \
		attraction.o scf.o uscf.o finalize.o uhfoperator.o udftoperator.o usedftoperator.o \
		uelectdii.o mpi_setup.o debug_quick.o dipole.o zmake.o calMP2.o optimize.o \
		gradient.o lbfgs.o hessian.o CPHF.o frequency.o

#************************************************************************
# 
#                 C. Make Executed Files
# 
#************************************************************************

quick: quick_modules.o $(OBJ) fake_amber_interface.o
	$(LD) $(OBJ) fake_amber_interface.o -o ${exefile} $(LDFLAGS) 

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
	-rm -f .cppdefs $(OBJ) quick
	
neat:
	-rm -f $(TMPFILES)

# - 2. Make tags for source files
TAGS: $(SRC)
	etags $(SRC)
tags: $(SRC)
	ctags $(SRC)
