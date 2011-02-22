!
!	quick_basis_module.f90
!	new_quick
!
!	Created by Yipu Miao on 2/18/11.
!	Copyright 2011 University of Florida. All rights reserved.
!

!********************************************************
!  Basis Set module
!--------------------------------------------------------
!
! A quick explanation :  aexp(i,j) is the alpha
! exponent of the ith primitive of the jth function, and d(i,j) is
! is multiplicative coefficient.  ncenter(i) is the center upon which
! the ith basis function is located.  ncontract(i) is the degree of
! contraction of the ith basis function, and itype(i) is the ith basis
! function's angular momentum.  The ifirst and ilast arrays denot the
! first and last basis function number on a given center, i.e. if center
! number 8 has basis functions 7-21, ifirst(8)=7 and ilast(8) = 21.
! If the calculation is closed shell, nelec is the total number of electrons,
! otherwise it is the number of alpha electrons and nelecb is the
! number of basis elctrons.  The xyz and ichrg arrays denote atomic
! position and charge, while iattype is the atomic number of the atom, i.e.
! 6 for carbon, etc. natom and nbasis are the number of atoms and
! basis functions. ibasis denotes the basis set chosen for.
!
    module quick_basis_module
    implicit none

    integer :: nshell,nbasis,nprim,jshell,jbasis
    integer :: IJKLtype,III,JJJ,KKK,LLL,IJtype,KLtype
    double precision :: Y,dnmax,Yaa(3),Ybb(3),Ycc(3)
    
    integer, allocatable, dimension(:) :: kstart,katom,ktype,kprim,kmin,kmax,ktypecp,atombasis

    double precision, allocatable, dimension(:) :: aex,gcs,gcp,gcd,gcf,gcg,eta,cons

    integer, allocatable, dimension(:) :: Qnumber,Ksumtype,Qstart,Qfinal 
    integer, allocatable, dimension(:,:) :: KLMN,Qsbasis,Qfbasis
    double precision, allocatable, dimension(:) :: gcexpomin
    double precision, allocatable, dimension(:,:) :: Apri,Kpri,gccoeff,gcexpo,Ycutoff, &
    debug1,debug2,cutmatrix,cutprim

    double precision, allocatable, dimension(:,:,:) :: Ppri,atomdens
    double precision, allocatable, dimension(:,:,:,:) :: Xcoeff,Yxiaoprim
    double precision, allocatable, dimension(:,:,:) :: orbmp2dcsub
    double precision, allocatable, dimension(:,:) :: orbmp2
    double precision, allocatable, dimension(:,:,:,:,:) :: orbmp2i331
    double precision, allocatable, dimension(:,:,:,:,:) :: orbmp2j331
    double precision, allocatable, dimension(:,:,:,:) :: orbmp2k331
    double precision, allocatable, dimension(:,:,:) :: orbmp2k331dcsub
!    double precision, allocatable, dimension(:,:,:) :: mem
!    integer, allocatable, dimension(:,:,:) :: CPmem
!    double precision, allocatable, dimension(:,:,:) :: Yxiao,Yxiaotemp 
    double precision, allocatable, dimension(:,:,:) :: Yxiao,Yxiaotemp,attraxiao,allerror,alloperator
    double precision, allocatable, dimension(:,:,:,:) :: attraxiaoopt
    double precision, allocatable, dimension(:) :: phixiao,dphidxxiao,dphidyxiao,dphidzxiao,Mulliken,Lowdin
!    dimension allerror(MAXDIISSCF,nbasis,nbasis)
!    dimension alloperator(MAXDIISSCF,nbasis,nbasis)

    ! MPI
    integer,allocatable:: mpi_jshelln(:)    ! shell no. calculated this node
    integer,allocatable:: mpi_jshell(:,:)   ! shell series no. calcluated in this node
 
    integer,allocatable:: mpi_nbasisn(:)    ! basis no. calculated this node
    integer,allocatable:: mpi_nbasis(:,:)   ! basis series no. calcluated in this node
 
    end module quick_basis_module
!********************************************************
