!
!	quick_calculated_module.f90
!	new_quick
!
!	Created by Yipu Miao on 2/18/11.
!	Copyright 2011 University of Florida. All rights reserved.
!

!********************************************************
! Calculated Module, storing calculated information
!--------------------------------------------------------
!
! For the above list, Smatrix is the overlap matrix, X is the
! orthogonalization matrix, O is the operator matrix, CO and COB are
! the molecular orbital coefficients, Vec is a mtrix to hold eigenvectors
! after a diagonalization, DENSE and DENSEB are the denisty matrices,
! V2 is a scratch space for diagonalization, IDEGEN is a matrix of
! orbital degeneracies, E and EB are the molecular orbital energies,
! and Eel, Ecore, Etot. alec and Belec are self explanatory. GRADIENT
! is the energy gradient with respect to atomic motion.
!
    module quick_calculated_module
    implicit none

! Xiao HE 01/13/2007

    double precision, dimension(:,:), allocatable :: Smatrix, &
    X,O,CO,COB,VEC,DENSE,DENSEB,V2,DENSEOLD,DENSESAVE,Osave,Uxiao,Osavedft
    integer, dimension(:), allocatable :: idegen
    double precision, dimension(:), allocatable :: E,EB
    double precision, dimension(:), allocatable :: gradient
    double precision, dimension(:,:), allocatable :: hessian,cphfb, &
    cphfa
    double precision, dimension(:), allocatable :: B0,BU
    double precision :: Eel,Ecore,Etot,EMP2,aelec,belec,Eelvac,Eelsol,Eelpb,Ecoremm,Gsolexp

    end module quick_calculated_module

!********************************************************
