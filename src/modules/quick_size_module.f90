!
!	size_module.f90
!	new_quick
!
!	Created by Yipu Miao on 2/17/11.
!	Copyright 2011 University of Florida. All rights reserved.
!

!--------------------------------------------------------
!  SIZE Module
!--------------------------------------------------------
!
!
! First, some sizing parameters:
    module quick_size_module
    implicit none

! Maximum number of atoms:
    integer, parameter :: MAXATM = 1000

! Maximum number of basis functions (total):
! integer, parameter :: MAXBASIS = 110
! integer, parameter :: MAXBASIS3 = 3*MAXBASIS

! Maximum contraction of each basis function:
! integer, parameter :: MAXCONTRACT = 3

! Maximum number of angular and radial grid points for quadrature:
    integer, parameter :: MAXANGGRID = 1400
    integer, parameter :: MAXRADGRID = 400

! Maximum number of DIIS error vectors for scf convergence.
!    integer, parameter :: MAXDIISSCF = 5
    integer :: MAXDIISSCF,NCYC

! M value of the lbfgs optimizer.
    integer, parameter :: MLBFGS = 5

    end module quick_size_module
!
!********************************************************
