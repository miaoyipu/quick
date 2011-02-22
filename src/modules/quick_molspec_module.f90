!
!	quick_molspec_module.f90
!	new_quick
!
!	Created by Yipu Miao on 2/18/11.
!	Copyright 2011 University of Florida. All rights reserved.
!

!********************************************************
! molecule specification Module
!--------------------------------------------------------
!
    module quick_molspec_module
    use quick_gaussian_class_module
    implicit none
! molchg :  total molecular charge
! iscf   :  MAX SCF CYCLES
! iopt   :  MAX OPT CYCLES
! pmaxrms:  DM Maximum RMS for convenrgency
! tol    :  DM cutoff
! chg    :  Atom Charge
! extchg :  External Charges
! nextatom: No. of External Charges
! nelec  :  No. of Electrons or Alpha Electrons
! nelecb :  No. of Beta Electrons
! Tdiag  :  Diag Time
! Tdcdiag:  DC Diag Time
! nNonHATOM: No. of Non-H Atoms
! NHAtom :  No. of H Atoms
! natom  :  No. of totoal Atoms
! xyz    :  Coordinates
! extxyz :  Coordinates of external charges
! AtomDistance: Atom Distance Matrix
! 

    integer :: maxcontract
    double precision, dimension(:,:), allocatable :: aexp,dcoeff
    double precision, dimension(:), allocatable :: distnbor
    double precision, dimension(:,:), allocatable :: xyz, extxyz, AtomDistance
    type (gaussian), dimension(:), allocatable :: gauss
    integer, dimension(:,:), allocatable :: itype
    integer, dimension(:), allocatable :: iattype,ncenter, &
    ncontract,ifirst,ilast,ishellfirst,ishelllast
    double precision, dimension(:), allocatable :: chg,extchg
    double precision :: dshift,tol,pmaxrms,acutoff,signif
    integer :: molchg
    integer :: natom,nelec,nelecb,iscf,iopt,nextatom,imult,nNonHAtom,nHAtom
    integer :: iatomtype

    end module quick_molspec_module
!********************************************************