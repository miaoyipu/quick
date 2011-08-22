! Ed Brothers. 11/26/01
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP
!********************************************************
! This file list and define varibles used as modules. 
! Modules List:
!    mfccmod        : MFCC Calculation
!    sizes          : Size Parameters
!    method         : Define Method Used
!    basis          : Basis Set Parameters
!    params         : Parameters for SEDFT
!    molspec        : Molecule Specifiction
!    gridpoints     : Grid points
!    calculated     : Overlap, Transformation matrix and energy
!    geocnverg      : Geometry Optimization parameters
!    SCRATCH        
!    files          : I/O file Specifiction
!    constants      : Constant
!    ecpmod         : ECP parameters
!    quickdc        : D&C Module
!    electrondensity: Electron Density
!    divpb_interface: Div_PB modules 
!    divpb_private  :

!********************************************************
!
!
!********************************************************
!  SIZE Module
!--------------------------------------------------------
!
!
! First, some sizing parameters:
    module sizes
    implicit none

! Maximum number of atoms:
! integer, parameter :: MAXATM = 50

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

    end module sizes
!
!********************************************************
!

!********************************************************
! Method Module
!--------------------------------------------------------
!
    module method
    implicit none

    logical :: optXiao = .false.   ! Optimization
    logical :: HF =  .false.       ! HF
    logical :: DFT =  .false.      ! DFT
    logical :: MP2 =  .false.      ! MP2
    logical :: B3LYP = .false.     ! B3LYP
    logical :: BLYP = .false.      ! BLYP
    logical :: BPW91 = .false.     ! BPW91
    logical :: MPW91LYP = .false.  ! MPW91LYP
    logical :: MPW91PW91 = .false. ! MPW91PW91
    logical :: PBSOL = .false.     ! PB Solvent
    logical :: UNRST =  .false.    ! Unrestricted
    logical :: NDDO =  .false.     !

    logical :: debug =  .false.    ! debug mode
    logical :: readDMX =  .false.  ! flag to read density matrix
    logical :: diisSCF =  .false.  ! DIIS SCF
    logical :: prtGap =  .false.   ! flag to print HOMO-LUMO gap
    logical :: opt =  .false.      ! optimization
    logical :: analGrad =  .false. ! Analytical Gradient
    logical :: analHess =  .false. ! Analytical Hessian Matrix
    logical :: fDiff1 =  .false.   ! 
    logical :: fDiff2 =  .false.   !
    logical :: diisOpt =  .false.  ! DIIS Optimization
    logical :: core =  .false.     !
    logical :: annil =  .false.    !
    logical :: freq =  .false.     ! Frenquency calculation
    logical :: SEDFT = .false.     ! Semi-Empirical DFT
    logical :: zmat = .false.      ! Z-matrix
!
    logical :: ecp,custECP         ! ECP
    logical :: energyXiao          ! Print Energy each cycle
    logical :: fFunXiao

!  
    logical :: calcDens = .false.    ! calculate density
    logical :: calcDensLap = .false. ! calculate density lap
    logical :: writePMat = .false.   ! Output density matrix
 
    logical :: extCharges = .false.  ! external charge
	logical :: bPDB=.false.          ! PDB input

    end module method
!********************************************************
!

!********************************************************
! Parameter Module
!--------------------------------------------------------
!
    module params
    implicit none

    double precision :: EK1prm(0:2,0:2,0:2,0:83), &
    At1prm(0:2,0:2,0:2,0:83), &
    At2prm(0:2,0:2,0:2,0:83), &
    Bndprm(0:2,0:2,0:2,0:83)
    double precision :: param7 = 1.d0
    double precision :: param8 = 1.d0
    double precision :: param9 = 1.d0
    double precision :: param10 = 1.d0

! SEDFT parameters for H.

    DATA EK1prm(0,0,0,1) / 2.946272d0 /
    DATA At1prm(0,0,0,1) / 1.885943d0 /

! SEDFT parameters for C.

    DATA EK1prm(0,0,0,6) / 0.599653d0 /
    DATA EK1prm(1,0,0,6) / 1.258270d0 /
    DATA EK1prm(0,1,0,6) / 1.258270d0 /
    DATA EK1prm(0,0,1,6) / 1.258270d0 /

    DATA At1prm(0,0,0,6) / 0.781378d0 /
    DATA At1prm(1,0,0,6) / 1.123674d0 /
    DATA At1prm(0,1,0,6) / 1.123674d0 /
    DATA At1prm(0,0,1,6) / 1.123674d0 /

! SEDFT parameters for N.

    DATA EK1prm(0,0,0,7) / 1.123192d0 /
    DATA EK1prm(1,0,0,7) / 1.050190d0 /
    DATA EK1prm(0,1,0,7) / 1.050190d0 /
    DATA EK1prm(0,0,1,7) / 1.050190d0 /

    DATA At1prm(0,0,0,7) / 0.866199d0 /
    DATA At1prm(1,0,0,7) / 1.027186d0 /
    DATA At1prm(0,1,0,7) / 1.027186d0 /
    DATA At1prm(0,0,1,7) / 1.027186d0 /

! SEDFT parameters for O
! TEMPORARY (THESE REALLY SUCK DONKEYS.)

    DATA EK1prm(0,0,0,8) / 0.467239d0 /
    DATA EK1prm(1,0,0,8) / 1.055385d0 /
    DATA EK1prm(0,1,0,8) / 1.055385d0 /
    DATA EK1prm(0,0,1,8) / 1.055385d0 /

    DATA At1prm(0,0,0,8) / 0.723604d0 /
    DATA At1prm(1,0,0,8) / 1.062795d0 /
    DATA At1prm(0,1,0,8) / 1.062795d0 /
    DATA At1prm(0,0,1,8) / 1.062795d0 /

! SEDFT parameters for F.
! TEMPORARY

    DATA EK1prm(0,0,0,9) / 0.617189d0 /
    DATA EK1prm(1,0,0,9) / 1.175433d0 /
    DATA EK1prm(0,1,0,9) / 1.175433d0 /
    DATA EK1prm(0,0,1,9) / 1.175433d0 /

    DATA At1prm(0,0,0,9) / 0.711736d0 /
    DATA At1prm(1,0,0,9) / 1.075532d0 /
    DATA At1prm(0,1,0,9) / 1.075532d0 /
    DATA At1prm(0,0,1,9) / 1.075532d0 /

    Integer :: Sumindex(-2:7)
    DATA Sumindex /0,0,1,4,10,20,35,56,84,120 /

    Integer :: trans(0:7,0:7,0:7)

    DATA trans(0,0,0) / 1 /
    DATA trans(0,0,1) / 4 /
    DATA trans(0,0,2) / 10 /
    DATA trans(0,0,3) / 20 /
    DATA trans(0,0,4) / 35 /
    DATA trans(0,1,0) / 3 /
    DATA trans(0,1,1) / 6 /
    DATA trans(0,1,2) / 17 /
    DATA trans(0,1,3) / 32 /
    DATA trans(0,2,0) / 9 /
    DATA trans(0,2,1) / 16 /
    DATA trans(0,2,2) / 23 /
    DATA trans(0,3,0) / 19 /
    DATA trans(0,3,1) / 31 /
    DATA trans(0,4,0) / 34 /
    DATA trans(1,0,0) / 2 /
    DATA trans(1,0,1) / 7 /
    DATA trans(1,0,2) / 15 /
    DATA trans(1,0,3) / 28 /
    DATA trans(1,1,0) / 5 /
    DATA trans(1,1,1) / 11 /
    DATA trans(1,1,2) / 26 /
    DATA trans(1,2,0) / 13 /
    DATA trans(1,2,1) / 25 /
    DATA trans(1,3,0) / 30 /
    DATA trans(2,0,0) / 8 /
    DATA trans(2,0,1) / 14 /
    DATA trans(2,0,2) / 22 /
    DATA trans(2,1,0) / 12 /
    DATA trans(2,1,1) / 24 /
    DATA trans(2,2,0) / 21 /
    DATA trans(3,0,0) / 18 /
    DATA trans(3,0,1) / 27 /
    DATA trans(3,1,0) / 29 /
    DATA trans(4,0,0) / 33 /

    DATA trans(1,2,2) / 36 /
    DATA trans(2,1,2) / 37 /
    DATA trans(2,2,1) / 38 /
    DATA trans(3,1,1) / 39 /
    DATA trans(1,3,1) / 40 /
    DATA trans(1,1,3) / 41 /
    DATA trans(0,2,3) / 42 /
    DATA trans(0,3,2) / 43 /
    DATA trans(2,0,3) / 44 /
    DATA trans(3,0,2) / 45 /
    DATA trans(2,3,0) / 46 /
    DATA trans(3,2,0) / 47 /
    DATA trans(0,1,4) / 48 /
    DATA trans(0,4,1) / 49 /
    DATA trans(1,0,4) / 50 /
    DATA trans(4,0,1) / 51 /
    DATA trans(1,4,0) / 52 /
    DATA trans(4,1,0) / 53 /
    DATA trans(5,0,0) / 54 /
    DATA trans(0,5,0) / 55 /
    DATA trans(0,0,5) / 56 /

    DATA trans(4,1,1) / 57 /
    DATA trans(1,4,1) / 58 /
    DATA trans(1,1,4) / 59 /
    DATA trans(1,2,3) / 60 /
    DATA trans(1,3,2) / 61 /
    DATA trans(2,1,3) / 62 /
    DATA trans(3,1,2) / 63 /
    DATA trans(2,3,1) / 64 /
    DATA trans(3,2,1) / 65 /
    DATA trans(2,2,2) / 66 /
    DATA trans(0,1,5) / 67 /
    DATA trans(0,5,1) / 68 /
    DATA trans(1,0,5) / 69 /
    DATA trans(5,0,1) / 70 /
    DATA trans(1,5,0) / 71 /
    DATA trans(5,1,0) / 72 /
    DATA trans(0,2,4) / 73 /
    DATA trans(0,4,2) / 74 /
    DATA trans(2,0,4) / 75 /
    DATA trans(4,0,2) / 76 /
    DATA trans(2,4,0) / 77 /
    DATA trans(4,2,0) / 78 /
    DATA trans(0,3,3) / 79 /
    DATA trans(3,0,3) / 80 /
    DATA trans(3,3,0) / 81 /
    DATA trans(6,0,0) / 82 /
    DATA trans(0,6,0) / 83 /
    DATA trans(0,0,6) / 84 /

    DATA trans(5,1,1) / 85 /
    DATA trans(1,5,1) / 86 /
    DATA trans(1,1,5) / 87 /
    DATA trans(1,2,4) / 88 /
    DATA trans(1,4,2) / 89 /
    DATA trans(2,1,4) / 90 /
    DATA trans(4,1,2) / 91 /
    DATA trans(2,4,1) / 92 /
    DATA trans(4,2,1) / 93 /
    DATA trans(1,3,3) / 94 /
    DATA trans(3,1,3) / 95 /
    DATA trans(3,3,1) / 96 /
    DATA trans(3,2,2) / 97 /
    DATA trans(2,3,2) / 98 /
    DATA trans(2,2,3) / 99 /
    DATA trans(0,1,6) / 100 /
    DATA trans(0,6,1) / 101 /
    DATA trans(1,0,6) / 102 /
    DATA trans(6,0,1) / 103 /
    DATA trans(1,6,0) / 104 /
    DATA trans(6,1,0) / 105 /
    DATA trans(0,2,5) / 106 /
    DATA trans(0,5,2) / 107 /
    DATA trans(2,0,5) / 108 /
    DATA trans(5,0,2) / 109 /
    DATA trans(2,5,0) / 110 /
    DATA trans(5,2,0) / 111 /
    DATA trans(0,3,4) / 112 /
    DATA trans(0,4,3) / 113 /
    DATA trans(3,0,4) / 114 /
    DATA trans(4,0,3) / 115 /
    DATA trans(3,4,0) / 116 /
    DATA trans(4,3,0) / 117 /
    DATA trans(7,0,0) / 118 /
    DATA trans(0,7,0) / 119 /
    DATA trans(0,0,7) / 120 /

    integer :: Mcal(3,120)
    DATA Mcal(1:3,1) /0,0,0/
    DATA Mcal(1:3,2) /1,0,0/
    DATA Mcal(1:3,3) /0,1,0/
    DATA Mcal(1:3,4) /0,0,1/
    DATA Mcal(1:3,5) /1,1,0/
    DATA Mcal(1:3,6) /0,1,1/
    DATA Mcal(1:3,7) /1,0,1/
    DATA Mcal(1:3,8) /2,0,0/
    DATA Mcal(1:3,9) /0,2,0/
    DATA Mcal(1:3,10) /0,0,2/
    DATA Mcal(1:3,11) /1,1,1/
    DATA Mcal(1:3,12) /2,1,0/
    DATA Mcal(1:3,13) /1,2,0/
    DATA Mcal(1:3,14) /2,0,1/
    DATA Mcal(1:3,15) /1,0,2/
    DATA Mcal(1:3,16) /0,2,1/
    DATA Mcal(1:3,17) /0,1,2/
    DATA Mcal(1:3,18) /3,0,0/
    DATA Mcal(1:3,19) /0,3,0/
    DATA Mcal(1:3,20) /0,0,3/
    DATA Mcal(1:3,21) /2,2,0/
    DATA Mcal(1:3,22) /2,0,2/
    DATA Mcal(1:3,23) /0,2,2/
    DATA Mcal(1:3,24) /2,1,1/
    DATA Mcal(1:3,25) /1,2,1/
    DATA Mcal(1:3,26) /1,1,2/
    DATA Mcal(1:3,27) /3,0,1/
    DATA Mcal(1:3,28) /1,0,3/
    DATA Mcal(1:3,29) /3,1,0/
    DATA Mcal(1:3,30) /1,3,0/
    DATA Mcal(1:3,31) /0,3,1/
    DATA Mcal(1:3,32) /0,1,3/
    DATA Mcal(1:3,33) /4,0,0/
    DATA Mcal(1:3,34) /0,4,0/
    DATA Mcal(1:3,35) /0,0,4/

    DATA Mcal(1:3,36) /1,2,2/
    DATA Mcal(1:3,37) /2,1,2/
    DATA Mcal(1:3,38) /2,2,1/
    DATA Mcal(1:3,39) /3,1,1/
    DATA Mcal(1:3,40) /1,3,1/
    DATA Mcal(1:3,41) /1,1,3/
    DATA Mcal(1:3,42) /0,2,3/
    DATA Mcal(1:3,43) /0,3,2/
    DATA Mcal(1:3,44) /2,0,3/
    DATA Mcal(1:3,45) /3,0,2/
    DATA Mcal(1:3,46) /2,3,0/
    DATA Mcal(1:3,47) /3,2,0/
    DATA Mcal(1:3,48) /0,1,4/
    DATA Mcal(1:3,49) /0,4,1/
    DATA Mcal(1:3,50) /1,0,4/
    DATA Mcal(1:3,51) /4,0,1/
    DATA Mcal(1:3,52) /1,4,0/
    DATA Mcal(1:3,53) /4,1,0/
    DATA Mcal(1:3,54) /5,0,0/
    DATA Mcal(1:3,55) /0,5,0/
    DATA Mcal(1:3,56) /0,0,5/

    DATA Mcal(1:3,57) /4,1,1/
    DATA Mcal(1:3,58) /1,4,1/
    DATA Mcal(1:3,59) /1,1,4/
    DATA Mcal(1:3,60) /1,2,3/
    DATA Mcal(1:3,61) /1,3,2/
    DATA Mcal(1:3,62) /2,1,3/
    DATA Mcal(1:3,63) /3,1,2/
    DATA Mcal(1:3,64) /2,3,1/
    DATA Mcal(1:3,65) /3,2,1/
    DATA Mcal(1:3,66) /2,2,2/
    DATA Mcal(1:3,67) /0,1,5/
    DATA Mcal(1:3,68) /0,5,1/
    DATA Mcal(1:3,69) /1,0,5/
    DATA Mcal(1:3,70) /5,0,1/
    DATA Mcal(1:3,71) /1,5,0/
    DATA Mcal(1:3,72) /5,1,0/
    DATA Mcal(1:3,73) /0,2,4/
    DATA Mcal(1:3,74) /0,4,2/
    DATA Mcal(1:3,75) /2,0,4/
    DATA Mcal(1:3,76) /4,0,2/
    DATA Mcal(1:3,77) /2,4,0/
    DATA Mcal(1:3,78) /4,2,0/
    DATA Mcal(1:3,79) /0,3,3/
    DATA Mcal(1:3,80) /3,0,3/
    DATA Mcal(1:3,81) /3,3,0/
    DATA Mcal(1:3,82) /6,0,0/
    DATA Mcal(1:3,83) /0,6,0/
    DATA Mcal(1:3,84) /0,0,6/

    DATA Mcal(1:3,85) /5,1,1/
    DATA Mcal(1:3,86) /1,5,1/
    DATA Mcal(1:3,87) /1,1,5/
    DATA Mcal(1:3,88) /1,2,4/
    DATA Mcal(1:3,89) /1,4,2/
    DATA Mcal(1:3,90) /2,1,4/
    DATA Mcal(1:3,91) /4,1,2/
    DATA Mcal(1:3,92) /2,4,1/
    DATA Mcal(1:3,93) /4,2,1/
    DATA Mcal(1:3,94) /1,3,3/
    DATA Mcal(1:3,95) /3,1,3/
    DATA Mcal(1:3,96) /3,3,1/
    DATA Mcal(1:3,97) /3,2,2/
    DATA Mcal(1:3,98) /2,3,2/
    DATA Mcal(1:3,99) /2,2,3/
    DATA Mcal(1:3,100) /0,1,6/
    DATA Mcal(1:3,101) /0,6,1/
    DATA Mcal(1:3,102) /1,0,6/
    DATA Mcal(1:3,103) /6,0,1/
    DATA Mcal(1:3,104) /1,6,0/
    DATA Mcal(1:3,105) /6,1,0/
    DATA Mcal(1:3,106) /0,2,5/
    DATA Mcal(1:3,107) /0,5,2/
    DATA Mcal(1:3,108) /2,0,5/
    DATA Mcal(1:3,109) /5,0,2/
    DATA Mcal(1:3,110) /2,5,0/
    DATA Mcal(1:3,111) /5,2,0/
    DATA Mcal(1:3,112) /0,3,4/
    DATA Mcal(1:3,113) /0,4,3/
    DATA Mcal(1:3,114) /3,0,4/
    DATA Mcal(1:3,115) /4,0,3/
    DATA Mcal(1:3,116) /3,4,0/
    DATA Mcal(1:3,117) /4,3,0/
    DATA Mcal(1:3,118) /7,0,0/
    DATA Mcal(1:3,119) /0,7,0/
    DATA Mcal(1:3,120) /0,0,7/

    end module params
!********************************************************
!

!********************************************************
! Gaussian Module
!--------------------------------------------------------
!
! Kenneth Ayers 1/23/04
! This file defines the data types needed for a higher abstraction level
! of gaussian orbitals.
    module gaussian_class
    implicit none

    type gaussian
    integer :: ncontract
    integer, dimension(3) :: itype
    double precision, pointer, dimension(:) :: aexp,dcoeff
    end type gaussian

    end module gaussian_class
!********************************************************


!********************************************************
! molecule specification Module
!--------------------------------------------------------
!
    module molspec
    use gaussian_class
    implicit none
! molchg :  total molecular charge
! iscf   :  MAX SCF CYCLES
! pmaxrms:  DM Maximum RMS for convenrgency
! tol    :  DM cutoff
! 

    integer :: maxcontract
    double precision, dimension(:,:), allocatable :: aexp,dcoeff
    double precision, dimension(:), allocatable :: distnbor
    double precision, dimension(:,:), allocatable :: xyz, extxyz, AtomDistance
    type (gaussian), dimension(:), allocatable :: gauss
    integer, dimension(:,:), allocatable :: itype
    integer, dimension(:), allocatable :: iattype,ncenter, &
    ncontract,ifirst,ilast
    double precision, dimension(:), allocatable :: chg,extchg
    double precision :: dshift,tol,pmaxrms,acutoff,Tdiag,Tdcdiag,molchg,signif
    integer :: natom,nelec,nelecb,iscf,iopt,nextatom,imult,nNonHAtom,nHAtom

    end module molspec
!********************************************************

!********************************************************
!  ECP module
!--------------------------------------------------------
! Alessandro GENONI 03/12/2006
!
    module ecpmod
     implicit none
!
     integer, parameter :: mxproj=5,mxang=3
!
! Derived Parameters
!
     integer, parameter :: mxnnn=max(2*mxang+1,mxang+mxproj+1),&
                           mxprim=(mxang+1)*(mxang+2)/2,&
                           mxgout=mxprim*mxprim,&
                           lmax1=max(1,mxang+max(mxang,mxproj)),&
                           lfdim=lmax1+1,&                            
                           lmfdim=lfdim**2,&
                           lmxdim=(lmax1*(lmax1+2)*(lmax1+4)/3 *  (lmax1+3) +&
                                  (lmax1+2)**2 * (lmax1+4))/16,&
                           mc1dim=2*mxproj-1,&
                           len_dfac=3*lmax1+3,&
                           len_fac=mxproj*mxproj
!
     integer :: necprim,nbf12,itolecp
     double precision :: tolecp     

     integer, dimension(:), allocatable   :: nelecp,lmaxecp,nlp,kvett
     double precision, dimension (:), allocatable :: clp,zlp,ecp_int,gout

     integer, dimension(:,:), allocatable :: kfirst,klast
!
     integer, dimension(:), allocatable   :: lf,lmf,lml,lmx,lmy,lmz
     integer, dimension(:,:), allocatable :: mc,mr 

     double precision, dimension(:), allocatable   :: zlm,dfac,dfaci,factorial
     double precision, dimension(:,:), allocatable :: flmtx,fprod

    end module ecpmod 
!********************************************************
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
    module basis
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
 
    end module basis
!********************************************************

!********************************************************
!  Grid Points Module
!--------------------------------------------------------
!
! The gridpoints arrays are fairly simple: XANG, YANG, and ZANG hold
! the angular grid points location for a unit sphere, and WTANG
! is the weights of those points.  RGRID and RWEIGHT are the positions
! and weights of the radial grid points, which in use are sclaed by the
! radii and radii^3 of the atoms.
!
    module gridpoints
    use sizes
    implicit none

    double precision ::  XANG(MAXANGGRID),YANG(MAXANGGRID), &
    ZANG(MAXANGGRID),WTANG(MAXANGGRID),RGRID(MAXRADGRID), &
    RWT(MAXRADGRID)
    double precision,  dimension(:), allocatable :: sigrad2
    integer :: iradial(0:10), &
    iangular(10),iregion

    end module gridpoints
!********************************************************

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
    module calculated
    implicit none

! Xiao HE 01/13/2007
    integer ISG

    double precision, dimension(:,:), allocatable :: Smatrix, &
    X,O,CO,COB,VEC,DENSE,DENSEB,V2,DENSEOLD,DENSESAVE,Osave,Uxiao,Osavedft
    integer, dimension(:), allocatable :: idegen
    double precision, dimension(:), allocatable :: E,EB
    double precision, dimension(:), allocatable :: gradient
    double precision, dimension(:,:), allocatable :: hessian,cphfb, &
    cphfa
    double precision, dimension(:), allocatable :: B0,BU
    double precision :: Eel,Ecore,Etot,EMP2,aelec,belec,Eelvac,Eelsol,Eelpb,Ecoremm,Gsolexp

    end module calculated

!********************************************************

!********************************************************
!   Geometry Convengency module
!--------------------------------------------------------
!
    module geocnverg
    implicit none

    double precision :: stepMax,geoTest,grmsTest,gradTest, &
    gnormTest,eTest

    end module geocnverg
!
!********************************************************

!
!********************************************************
! SCRATCH module.
!--------------------------------------------------------
!
    module SCRATCH
    implicit none

    double precision, dimension(:,:), allocatable :: HOLD,HOLD2

    end module SCRATCH
!
!********************************************************
!

!
!********************************************************
! File module.
!--------------------------------------------------------
!
    module files
    implicit none

    character(len=80) :: inFileName,outFileName,dmxFileName,rstFileName, &
    CPHFFileName,basisDir,basisFileName,ECPDir,ECPFileName,basisCustName, &
	PDBFileName
    integer :: inFile = 15            ! input file
    integer :: iOutFile = 16          ! output file
    integer :: iDmxFile = 17          ! density matrix file
    integer :: iRstFile = 18          ! Restricted file
    integer :: iCPHFFile = 19         ! CPHF file
    integer :: iBasisFile = 20        ! basis set file
    integer :: iECPFile = 21          ! ECP file
    integer :: iBasisCustFile = 22
	integer :: iPDBFile=23            ! PDB file for D&C

    end module files
!
!********************************************************

!********************************************************
! Contants module
!--------------------------------------------------------
!
    module constants
    implicit none

    double precision :: PI,X0,X00,xiaoCUTOFF,xiaoCUTOFF1,primLimit,gradCutoff,PIto3half
 
    double precision, parameter :: bohr = 0.5291772083d0

    character(len=2), dimension(0:92) :: symbol = &
   & (/'XX','H ','HE','LI','BE','B ','C ','N ','O ','F ','NE', &
   & 'NA','MG','AL','SI','P ','S ','CL','AR','K ','CA', &
   & 'SC','TI','V ','CR','MN','FE','CO','NI','CU','ZN', &
   & 'GA','GE','AS','SE','BR','KR','RB','SR','Y ','ZR', &
   & 'NB','MO','TC','RU','RH','PD','AG','CD','IN','SN', &
   & 'SB','TE','I ','XE','CS','BA','LA','CE','PR','ND', &
   & 'PM','SM','EU','GD','TB','DY','HO','ER','TM','YB', &
   & 'LU','HF','TA','W ','RE','OS','IR','PT','AU','HG', &
   & 'TL','PB','BI','PO','AT','RN','FR','RA','AC','TH', &
   & 'PA','U '/)

    character(len=2), dimension(1:10) :: atomXiao

    double precision, dimension(0:83) :: eMass

    data emass &
    /0.d0, 1.0079d0, 4.0026d0, 6.941d0, 9.01218d0, &
    10.81d0, 12.011d0,14.0067d0, 15.99994d0, 18.99840d0, &
    20.179d0, 22.9898d0, 24.305d0, 26.98154d0, 28.0855d0, &
    30.97376d0, 32.06d0, 35.453d0, 39.948d0,39.0983d0, &
    40.08d0, 44.9559d0, 47.90d0, 50.9415d0, 51.996d0, &
    54.938d0, 55.847d0, 58.9332d0, 58.71d0, 63.546d0, &
    65.38d0, 69.737d0, 72.59d0, 74.9216d0, 78.96d0, 79.904d0, &
    83.80d0, 85.4678d0, 87.62d0, 88.9059d0, 91.22d0, 92.9064d0, &
    95.94d0, 98.9062d0, 101.07d0, 102.9055d0, 106.4d0, &
    107.868d0, 112.41d0, 114.82d0, 118.69d0, 121.75d0, 127.60d0, &
    126.9045d0, 131.30d0, 132.9054d0, 137.33d0, 15*0.000d0, &
    178.49d0, 180.9479d0, 183.850d0, 186.207d0, 190.20d0, &
    192.220d0, 195.090d0, 196.9665d0, 200.590d0, 204.370d0, &
    207.200d0, 208.9804d0/


    double precision, dimension(-2:30) :: fact = &
    (/   0.d0,0.d0,1.d0,1.d0,  2.000000000000000D0, &
    6.000000000000000D0,24.00000000000000D0   , &
    120.0000000000000D0,720.0000000000000D0   , &
    5040.000000000000D0,40320.00000000000D0   , &
    362880.0000000000D0,3628800.000000000D0   , &
    39916800.00000000D0,479001600.0000000D0   , &
    6227020800.000000D0,87178291200.00000D0   , &
    1307674368000.000D0,20922789888000.00D0   , &
    355687428096000.0D0,6402373705728000.D0   , &
    1.2164510040883200D+17,2.4329020081766400D+18, &
    5.1090942171709440D+19,1.1240007277776077D+21, &
    2.5852016738884978D+22,6.2044840173323941D+23, &
    1.5511210043330986D+25,4.0329146112660565D+26, &
    1.0888869450418352D+28,3.0488834461171384D+29, &
    8.8417619937397008D+30,2.6525285981219103D+32/)

    double precision, dimension(0:83) :: radii
    double precision, dimension(0:83) :: radii2

    data radii &
    /0.d0,1.d0,0.5882d0,3.0769d0,2.0513d0,1.5385d0, &
    1.2308d0,1.0256d0,0.8791d0,0.7692d0,0.6838d0, &
    4.0909d0,3.1579d0,2.5714d0,2.1687d0,1.8750d0, &
    1.6514d0,1.4754d0,1.3333d0,65*2.25d0/

! Xiao HE 02/11/2007    
    data radii2 &
    /0.d0,1.30d0,0.0d0,1.95d0,2.20d0,1.45d0, &
    1.20d0,1.10d0,1.10d0,1.20d0,0.0d0, &
    2.30d0,2.20d0,2.10d0,1.30d0,1.30d0, &
    1.10d0,1.45d0,0.0d0,65*2.25d0/

    end module constants
!********************************************************

!
!********************************************************
! electron densitry modules
!--------------------------------------------------------
!
! John F. 11/25/2008
    module electrondensity
    implicit none
   
    double precision, dimension(:,:,:), allocatable :: elecDense,elecDenseLap
    double precision :: xStart,yStart,zStart,gridSpacing,lapGridSpacing
    integer :: numPointsX,numPointsY,numPointsZ
    character(len=80) :: dxFileName,xyzFileName
     
    end module electrondensity 
!********************************************************

!********************************************************
! D&C module
!--------------------------------------------------------
!
    module quickdc
    implicit none

    integer, dimension(:,:), allocatable :: DCCore,DCBuffer1,DCBuffer2,DCSub
    integer, dimension(:), allocatable :: DCCoren,DCBuffer1n,DCBuffer2n,DCSubn,nBasisDC, &
                                          nElecDCSub,selectNN,nElecMP2Sub
    integer, dimension(:,:), allocatable :: DCOverlap,DCConnect
    integer, dimension(0:92)  :: kShell
    integer, dimension(:), allocatable :: kShellS,kShellF
    integer, dimension(:,:,:), allocatable :: DCLogic
    double precision, dimension(:,:), allocatable :: invDCOverlap
    double precision, dimension(:,:,:), allocatable :: ODCSub,PDCSub,XDCSub,SMatrixDCSub, &
                                                       coDCSub,PDCSubtran,coDCSubtran
    double precision, dimension(:,:), allocatable :: ODCSubtemp,VECtemp
    double precision, dimension(:,:), allocatable :: Vtemp,EVEC1temp,eValDCSub
    double precision, dimension(:), allocatable :: EVAL1temp,IDEGEN1temp
    logical, dimension(:,:), allocatable :: disDivMFCC
    logical, dimension(:), allocatable :: mp2Shell
    integer np,number

	character*6,allocatable:: sn(:)
    real*8,allocatable::coord(:,:)
    integer,allocatable::class(:),ttnumber(:)
    character*4,allocatable::atomname(:)
    character*3,allocatable::residue(:)
    integer,allocatable::selectC(:),charge(:),spin(:)
    integer,allocatable::selectN(:),selectCA(:),ccharge(:)
    character*200 cmdstr
	
	integer ifragbasis ! =2.residue basis,=1.atom basis,=3 non-h atom basis
	
    end module quickdc
!********************************************************


!********************************************************
! MFCC Module
!--------------------------------------------------------
    module mfccmod
    implicit none

!    integer, allocatable, dimension(:) :: mfccatom,mfcccharge
    integer :: mfccatom(50),mfcccharge(50),npmfcc,IMFCC,kxiaoconnect,IFMM,IDIVCON
    double precision :: mfcccord(3,100,50)
	integer ::Ftmp(300)
	character(len=100)::linetmp
    character(len=2) :: mfccatomxiao(100,50)
    integer :: mfccstart(50),mfccfinal(50),mfccbases(50),mfccbasef(50)
    integer :: matomstart(50),matomfinal(50),matombases(50),matombasef(50)
!    integer, dimension(:), allocatable :: matomstart,matomfinal,matombases &
!    ,matombasef

    integer :: mfccatomcap(50),mfccchargecap(50)
    double precision :: mfcccordcap(3,100,50)
    character(len=2) :: mfccatomxiaocap(100,50)
    integer :: mfccstartcap(50),mfccfinalcap(50),mfccbasescap(50),mfccbasefcap(50)
    integer :: matomstartcap(50),matomfinalcap(50),matombasescap(50) &
    ,matombasefcap(50)
!    integer, dimension(:), allocatable :: matomstartcap,matomfinalcap,matombasescap &
!    ,matombasefcap

    integer :: mfccatomcon(50),mfccchargecon(50)
    double precision :: mfcccordcon(3,100,50)
    character(len=2) :: mfccatomxiaocon(100,50)
    integer :: mfccstartcon(50),mfccfinalcon(50),mfccbasescon(50),mfccbasefcon(50)
    integer :: matomstartcon(50),matomfinalcon(50),matombasescon(50) &
    ,matombasefcon(50)
!    integer, dimension(:), allocatable :: matomstartcap,matomfinalcap,matombasescap &
!    ,matombasefcap

    integer :: mfccatomcon2(50),mfccchargecon2(50)
    double precision :: mfcccordcon2(3,100,50)
    character(len=2) :: mfccatomxiaocon2(100,50)
    integer :: mfccstartcon2(50),mfccfinalcon2(50),mfccbasescon2(50),mfccbasefcon2(50)
    integer :: matomstartcon2(50),matomfinalcon2(50),matombasescon2(50) &
    ,matombasefcon2(50)

    integer :: mfccatomconi(50),mfccchargeconi(50)
    double precision :: mfcccordconi(3,100,50)
    character(len=2) :: mfccatomxiaoconi(100,50)
    integer :: mfccstartconi(50),mfccfinalconi(50),mfccbasesconi(50),mfccbasefconi(50)
    integer :: matomstartconi(50),matomfinalconi(50),matombasesconi(50) &
    ,matombasefconi(50)

    integer :: mfccatomconj(50),mfccchargeconj(50)
    double precision :: mfcccordconj(3,100,50)
    character(len=2) :: mfccatomxiaoconj(100,50)
    integer :: mfccstartconj(50),mfccfinalconj(50),mfccbasesconj(50),mfccbasefconj(50)
    integer :: matomstartconj(50),matomfinalconj(50),matombasesconj(50) &
    ,matombasefconj(50)

    double precision, allocatable, dimension(:,:,:) :: mfccdens,mfccdenscap,mfccdenscon &
                                            ,mfccdenscon2,mfccdensconi,mfccdensconj

    end module mfccmod
!********************************************************

!********************************************************
! Following modules contains all global variables used in divpb.              
! -- Ning Liao 05/15/2004
! Add to modules.f90
! -- Yipu Miao 05/07/2010
! div PB modules includes:
! divpb_interface: 
! divpb_private:
!    contains: initialize_divpbVars()
!              deallocate_divpbVars(ierror)
!----------------------------------------------------------------------
!
!----------------------------------------------------------------------
! Interface modules
!----------------------------------------------------------------------
!
  module divpb_interface
! 
!comment
!
! Interface for DIVPB, any subroute outside of DIVPB need to use this  
! module: "use divpb_interface"                                        
!
!comment_end

  implicit none
  save

  ! original divcon variables, different name
  integer:: gnIscr, gnIout
  logical:: gbWater, gbOctanol
  logical:: gbScreen
  integer:: iiixiao

  ! the return of PB calculation - Reaction field energy & nonpolar energy
  real:: grGrf, grGnp,grGnpxiao

  ! the return of PB calculation - suface chargs & their coordinates
  ! will be used in other divcon modules (mostly pbfock)
  real, dimension(:,:), allocatable:: grpSchrgPos
  real, dimension(:), allocatable:: grpSchrg
  ! real, dimension(:,:), pointer:: grpSchrgPos
  ! real, dimension(:), pointer:: grpSchrg
  integer:: gnSchrgNum
  ! common/divpbSchrg/ grpSchrgPos, grpSchrg, gnSchrgNum

  ! keywords read in by divcon's rdkeys.F
  real :: grDielecIn, grDielecOut
  real :: grProbeRadius, grGridPerAng, grPerFill, grIonStrength
  logical :: gbUseDivpb
  logical :: gbFineSurf
  integer :: gnPBDump
  real :: grUsrCentX, grUsrCentY, grUsrCentZ
  integer :: gnUsrGridNum
  real :: grESPCentX, grESPCentY, grESPCentZ, grESPEdge
  logical:: gbPBInstant

  ! these are used to determine whether assigned value in divcon input or need
  ! to set default values in setdef.F
  data gbUseDivpb /.true./
  data gbFineSurf /.false./
  data gbPBInstant /.false./

  data grDielecIn /-1.0/
  data grDielecOut /-1.0/
  data grProbeRadius /-1.0/
  data grGridPerAng /-1.0/
  data grPerFill /-1.0/
  data grIonStrength /-1.0/
  data gnPBDump /0/

  data grUsrCentX/12.345/
  data grUsrCentY/12.345/
  data grUsrCentZ/12.345/
  data gnUsrGridNum/12345/

  data grESPCentX/12.345/
  data grESPCentY/12.345/
  data grESPCentZ/12.345/
  data grESPEdge/12.345/
end module divpb_interface
!
!----------------------------------------------------------------------
! Private modules
!----------------------------------------------------------------------

module divpb_private

  !comment
  !
  ! Private global data for DIVPB use only          
  !
  !comment_end


  implicit none
  save

include "dpparameter.h"


  ! Molecular information
  character(len=15), dimension(:), pointer::gspMolInfo
  real, dimension(:,:), pointer:: grpMolPos
  real, dimension(:), pointer:: grpMolChrg, grpMolRadii
  real:: graMolSpanMid(3)
  integer :: gnMolAtomNum
  ! common/divpbMol/ gspMolInfo,grpMolPos,grpMolChrg,grpMolRadii,
  ! .              graMolSpanMid,gnMolAtomNum

  ! Grid information
  ! grpGridPhi(1:gnGridNum, 1:gnGridNum, 1:gnGridNum)
  ! gnpGridkappa(1:gnGridNum, 1:gnGridNum, 1:gnGridNum)
  ! grpGridEpsilon(1:gnGridNum, 1:gnGridNum, 1:gnGridNum, 1:3)
  ! grpGridMol(1:3, 1:gnMolAtomNum)
  ! gnpGridBndry(1:3, 1:gnGridBndryNum)
  ! grpGridChrg(1:gnGridChrgNum)
  ! gnpGridChrgXYZ(1:3, 1:gnGridChrgNum)
  real, dimension(:,:,:), pointer:: grpGridPhi
  integer, dimension(:,:,:), pointer:: gnpGridKappa
  real, dimension(:,:,:,:), pointer:: grpGridEpsilon
  real, dimension(:,:), pointer:: grpGridMol
  integer, dimension(:,:), pointer:: gnpGridBndry
  integer, dimension(:,:), pointer:: gnpGridChrgXYZ
  real, dimension(:), pointer:: grpGridChrg
  integer:: gnGridNum, gnGridBndryNum, gnGridChrgNum
  ! common/divpbGrid/ grpGridPhi,gnpGridKappa,
  ! .               grpGridEpsilon,grpGridMol,gnpGridBndry,
  ! .               grpGridChrg, gnpGridChrgXYZ,
  ! .               gnGridNum,gnGridBndryNum, gnGridChrgNum

  ! Surface information
  ! grpSurfExtVdw(1:3, 1:gnSurfExtVdwNum) coordinates of extended vdw surface pnts
  ! grpSurfMS(1:3, 1:gnSurfMSNum) coordinates of molecular surface pnts
  real, dimension(:,:), pointer:: grpSurfExtVdw, grpSurfMS
  real, dimension(:), pointer:: grpSurfMSArea
  integer:: gnSurfExtVdwNum, gnSurfMSNum, gnSurfType
  logical:: gbSurfFirstSpt
  ! common/divpbSurf/ grpSurfExtVdw, grpSurfMS, grpSurfMSArea,
  ! .                  gnSurfExtVdwNum, gnSurfMSNum, gnSurfType,
  ! .                  gbSurfFirstSpt


  ! arrays to accelerate relaxation of PB equation.
  ! grpAccEpsRes(0:6, 1:gnGridBndryNum) the residue (due to the un-unified epsilon)
  ! add up to the laplace version of PB.
  ! grpAccChrgRes(1:gnGridChrgNum) the residue due to the charge
  ! grpAccKappaLap(1:gnGridNum,1:gnGridNum,1:gnGridNum) the laplace equation counting
  ! in the kappa item in the salt solution.
  real, dimension(:,:,:), pointer:: grpAccEpsRes
  real, dimension(:,:), pointer:: grpAccChrgRes
  real, dimension(:), pointer:: grpAccPhiOdd,grpAccKappaOdd
  real, dimension(:), pointer:: grpAccPhiEven,grpAccKappaEven
  real, dimension(:,:), pointer:: grpAccBoxBndry
  integer, dimension(:,:), pointer:: gnpAccEpsIndx
  integer, dimension(:,:), pointer:: gnpAccChrgIndx
  integer, dimension(:,:), pointer:: gnpAccBoxBndryIndx
  integer, dimension(:,:), pointer:: gnpAccBndryChrgIndx
  integer:: gnAccBndryChrgNum
  ! common/divpbAcc/ grpAccEpsRes, grpAccChrgRes,
  ! .                grpAccKappaOdd, grpAccKappaEven,
  ! .                grpAccPhiOdd, grpAccPhiEven,grpAccBoxBndry,
  ! .                gnpAccEpsIndx, gnpAccChrgIndx,
  ! .                gnpAccBoxBndryIndx, gnpAccBndryChrgIndx,
  ! .                gnAccBndryChrgNum

  ! arrays to use finer grid at the dielectric boundary
  logical :: gbFSNow
  integer :: gnDblGridChrgNum
  real, dimension(:,:), pointer:: grpDblGridMol
  real, dimension(:), pointer:: grpDblGridChrg
  integer, dimension(:,:), pointer:: gnpDblGridChrgXYZ
  real, dimension(:,:,:), pointer:: grpDblPhi
  integer, dimension(:,:), pointer:: gnpFSLocalPnts
  integer, dimension(:,:), pointer:: gnpFSLocalChrgs

  real, dimension(:,:,:), pointer:: grpFSAccEpsRes
  integer, dimension(:,:), pointer:: gnpFSAccEpsIndx
  real, dimension(:,:), pointer:: grpFSAccChrgRes
  integer, dimension(:,:), pointer:: gnpFSAccChrgIndx
  real, dimension(:), pointer:: grpFSAccKappaOdd, grpFSAccKappaEven
  integer, dimension(:,:), pointer:: gnpFSAccBoxBndryIndx
  real, dimension(:,:), pointer:: grpFSAccBoxBndry
  integer :: gnFSAccBndryChrgNum
  integer, dimension(:,:), pointer:: gnpFSAccBndryChrgIndx
  real, dimension(:), pointer:: grpFSAccPhiOdd, grpFSAccPhiEven

  ! arrays for radii information
  character(len=2):: gsaRadAtomName(MAX_RADII)
  real:: graRadAtomRadii(MAX_RADII)
  integer:: gnRadAtomNum

  ! radii information for atoms used in SCRF caculation
  ! hhk : moved to the routine, initialize_divpbVars()

  ! hhk Below is added to control variables which have 'save' attribute and
  ! need to be cleaned for each Divcon run job.

  ! Used only in pbsolver.F
  logical :: bFirst
  integer :: count

  ! Used only in divpb.F
  logical :: bFirstRun_divpb
  real :: rSpectral_divpb

  ! Used only in dpchargegrid.F
  logical :: bFirstRun_dpchrg
  real :: rLastChrg

  ! Used only in dpsetepsilon.F
  integer, dimension(:,:), pointer :: npBndry1, npBndry2
  logical :: bUseBndry1
  integer :: nAllocated

  ! Used only in pb_spt.F
  real, dimension(:,:), pointer :: rpSurfVdw1, rpSurfVdw2
  real, dimension(:,:), pointer :: rpSurfMS1, rpSurfMS2
  real, dimension(:), pointer :: rpSurfMSArea1, rpSurfMSArea2
  logical :: bUseVdw1, bUseMS1
  integer :: nVdwAllocated, nMSAllocated

  ! Used only in pb_putpnt.F
  logical :: first_pb_putpnt


CONTAINS

  !========================================================================

  subroutine initialize_divpbVars()

  !comment
  !
  ! Initialize variables (used in DivPB calculation) which need to be
  ! fresh ones for multiple runs of divcon jobs (NOT divpb jobs).
  !
  !comment_end


    gnMolAtomNum = 0

    ! Used only in pbsolver.F
    bFirst = .TRUE.
    count = 0

    ! Used only in divpb.F
    bFirstRun_divpb = .TRUE.
    rSpectral_divpb = 0.0

    ! Used only in dpchargegrid.F
    bFirstRun_dpchrg = .TRUE.
    rLastChrg = 0.0

    ! Used only in dpsetepsilon.F
    bUseBndry1 = .TRUE.
    nAllocated = 0

    ! Used only in pb_spt.F
    bUseVdw1 = .FALSE.
    bUseMS1 = .FALSE.
    nVdwAllocated = 0
    nMSAllocated = 0

    ! Used only in pb_putpnt.F
    first_pb_putpnt = .TRUE.


    ! radii information for atoms used in SCRF caculation
    gnRadAtomNum = 83

    gsaRadAtomName(1) = 'H '
    graRadAtomRadii(1) = 1.150

    gsaRadAtomName(2) = 'HE'
    graRadAtomRadii(2) = 1.181

    gsaRadAtomName(3) = 'LI'
    graRadAtomRadii(3) = 1.226

    gsaRadAtomName(4) = 'BE'
    graRadAtomRadii(4) = 1.373

    gsaRadAtomName(5) = 'B '
    graRadAtomRadii(5) = 2.042

    gsaRadAtomName(6) = 'C '
    graRadAtomRadii(6) = 1.900

    gsaRadAtomName(7) = 'N '
    graRadAtomRadii(7) = 1.600

    gsaRadAtomName(8) = 'O '
    graRadAtomRadii(8) = 1.600

    gsaRadAtomName(9) = 'F '
    graRadAtomRadii(9) = 1.682

    gsaRadAtomName(10) = 'NE'
    graRadAtomRadii(10) = 1.621

    gsaRadAtomName(11) = 'NA'
    graRadAtomRadii(11) = 1.491

    gsaRadAtomName(12) = 'MG'
    graRadAtomRadii(12) = 1.510

    gsaRadAtomName(13) = 'AL'
    graRadAtomRadii(13) = 2.249

    gsaRadAtomName(14) = 'SI'
    graRadAtomRadii(14) = 2.147

    gsaRadAtomName(15) = 'P '
    graRadAtomRadii(15) = 2.074

    gsaRadAtomName(16) = 'S '
    graRadAtomRadii(16) = 1.900

    gsaRadAtomName(17) = 'CL'
    graRadAtomRadii(17) = 1.974

    gsaRadAtomName(18) = 'AR'
    graRadAtomRadii(18) = 1.934

    gsaRadAtomName(19) = 'K '
    graRadAtomRadii(19) = 1.906

    gsaRadAtomName(20) = 'CA'
    graRadAtomRadii(20) = 1.700

    gsaRadAtomName(21) = 'SC'
    graRadAtomRadii(21) = 1.647

    gsaRadAtomName(22) = 'TI'
    graRadAtomRadii(22) = 1.587

    gsaRadAtomName(23) = 'V '
    graRadAtomRadii(23) = 1.572

    gsaRadAtomName(24) = 'CR'
    graRadAtomRadii(24) = 1.511

    gsaRadAtomName(25) = 'MN'
    graRadAtomRadii(25) = 1.480

    gsaRadAtomName(26) = 'FE'
    graRadAtomRadii(26) = 1.456

    gsaRadAtomName(27) = 'CO'
    graRadAtomRadii(27) = 1.436

    gsaRadAtomName(28) = 'NI'
    graRadAtomRadii(28) = 1.417

    gsaRadAtomName(29) = 'CU'
    graRadAtomRadii(29) = 1.748

    gsaRadAtomName(30) = 'ZN'
    graRadAtomRadii(30) = 1.381

    gsaRadAtomName(31) = 'GA'
    graRadAtomRadii(31) = 2.192

    gsaRadAtomName(32) = 'GE'
    graRadAtomRadii(32) = 2.140

    gsaRadAtomName(33) = 'AS'
    graRadAtomRadii(33) = 2.115

    gsaRadAtomName(34) = 'SE'
    graRadAtomRadii(34) = 2.103

    gsaRadAtomName(35) = 'BR'
    graRadAtomRadii(35) = 2.095

    gsaRadAtomName(36) = 'KR'
    graRadAtomRadii(36) = 2.071

    gsaRadAtomName(37) = 'RB'
    graRadAtomRadii(37) = 2.057

    gsaRadAtomName(38) = 'SR'
    graRadAtomRadii(38) = 1.821

    gsaRadAtomName(39) = 'Y '
    graRadAtomRadii(39) = 1.673

    gsaRadAtomName(40) = 'ZR'
    graRadAtomRadii(40) = 1.562

    gsaRadAtomName(41) = 'NB'
    graRadAtomRadii(41) = 1.583

    gsaRadAtomName(42) = 'MO'
    graRadAtomRadii(42) = 1.526

    gsaRadAtomName(43) = 'TC'
    graRadAtomRadii(43) = 1.499

    gsaRadAtomName(44) = 'RU'
    graRadAtomRadii(44) = 1.481

    gsaRadAtomName(45) = 'RH'
    graRadAtomRadii(45) = 1.464

    gsaRadAtomName(46) = 'PD'
    graRadAtomRadii(46) = 1.450

    gsaRadAtomName(47) = 'AG'
    graRadAtomRadii(47) = 1.574

    gsaRadAtomName(48) = 'CD'
    graRadAtomRadii(48) = 1.424

    gsaRadAtomName(49) = 'IN'
    graRadAtomRadii(49) = 2.232

    gsaRadAtomName(50) = 'SN'
    graRadAtomRadii(50) = 2.196

    gsaRadAtomName(51) = 'SB'
    graRadAtomRadii(51) = 2.210

    gsaRadAtomName(52) = 'TE'
    graRadAtomRadii(52) = 2.235

    gsaRadAtomName(53) = 'I '
    graRadAtomRadii(53) = 2.250

    gsaRadAtomName(54) = 'XE'
    graRadAtomRadii(54) = 2.202

    gsaRadAtomName(55) = 'CS'
    graRadAtomRadii(55) = 2.259

    gsaRadAtomName(56) = 'BA'
    graRadAtomRadii(56) = 1.851

    gsaRadAtomName(57) = 'LA'
    graRadAtomRadii(57) = 1.761

    gsaRadAtomName(72) = 'HF'
    graRadAtomRadii(72) = 1.570

    gsaRadAtomName(73) = 'TA'
    graRadAtomRadii(73) = 1.585

    gsaRadAtomName(74) = 'W '
    graRadAtomRadii(74) = 1.534

    gsaRadAtomName(75) = 'RE'
    graRadAtomRadii(75) = 1.477

    gsaRadAtomName(76) = 'OS'
    graRadAtomRadii(76) = 1.560

    gsaRadAtomName(77) = 'IR'
    graRadAtomRadii(77) = 1.420

    gsaRadAtomName(78) = 'PT'
    graRadAtomRadii(78) = 1.377

    gsaRadAtomName(79) = 'AU'
    graRadAtomRadii(79) = 1.647

    gsaRadAtomName(80) = 'HG'
    graRadAtomRadii(80) = 1.353

    gsaRadAtomName(81) = 'TL'
    graRadAtomRadii(81) = 2.174

    gsaRadAtomName(82) = 'PB'
    graRadAtomRadii(82) = 2.148

    gsaRadAtomName(83) = 'BI'
    graRadAtomRadii(83) = 2.185

  end subroutine initialize_divpbVars

  !========================================================================

  subroutine deallocate_divpbVars(ierror)

  !comment
  !
  ! Deallocate dynamic arrays used in DivPB calculation.
  !
  !comment_end

    use divpb_interface

    IMPLICIT NONE

    INTEGER :: iDeallocateErr, ierror

    ierror = 0

    if (allocated(grpSchrgPos)) then
       deallocate(grpSchrgPos, stat=iDeallocateErr)
       if (iDeallocateErr /= 0) then
          ierror = 1
          return
       endif
    endif
    if (allocated(grpSchrg)) then
       deallocate(grpSchrg, stat=iDeallocateErr)
       if (iDeallocateErr /= 0) then
          ierror = 1
          return
       endif
    endif
    if (associated(gspMolInfo)) then
       deallocate(gspMolInfo, stat=iDeallocateErr)
       if (iDeallocateErr /= 0) then
          ierror = 1
          return
       endif
    endif
    if (associated(grpMolPos)) then
       deallocate(grpMolPos, stat=iDeallocateErr)
       if (iDeallocateErr /= 0) then
          ierror = 1
          return
       endif
    endif
    if (associated(grpMolChrg)) then
       deallocate(grpMolChrg, stat=iDeallocateErr)
       if (iDeallocateErr /= 0) then
          ierror = 1
          return
       endif
    endif
    if (associated(grpMolRadii)) then
       deallocate(grpMolRadii, stat=iDeallocateErr)
       if (iDeallocateErr /= 0) then
          ierror = 1
          return
       endif
    endif
    if (associated(grpGridPhi)) then
       deallocate(grpGridPhi, stat=iDeallocateErr)
       if (iDeallocateErr /= 0) then
          ierror = 1
          return
       endif
    endif
    if (associated(gnpGridKappa)) then
       deallocate(gnpGridKappa, stat=iDeallocateErr)
       if (iDeallocateErr /= 0) then
          ierror = 1
          return
       endif
    endif
    if (associated(grpGridEpsilon)) then
       deallocate(grpGridEpsilon, stat=iDeallocateErr)
       if (iDeallocateErr /= 0) then
          ierror = 1
          return
       endif
    endif
    if (associated(grpGridMol)) then
       deallocate(grpGridMol, stat=iDeallocateErr)
       if (iDeallocateErr /= 0) then
          ierror = 1
          return
       endif
    endif
    if (associated(gnpGridChrgXYZ)) then
       deallocate(gnpGridChrgXYZ, stat=iDeallocateErr)
       if (iDeallocateErr /= 0) then
          ierror = 1
          return
       endif
    endif
    if (associated(grpGridChrg)) then
       deallocate(grpGridChrg, stat=iDeallocateErr)
       if (iDeallocateErr /= 0) then
          ierror = 1
          return
       endif
    endif

    ! These were not allocated but pointed to other allocated pointers.
    if (associated(gnpGridBndry)) then
       nullify(gnpGridBndry)
    endif
    if (associated(grpSurfExtVdw)) then
       nullify(grpSurfExtVdw)
    endif
    if (associated(grpSurfMS)) then
       nullify(grpSurfMS)
    endif
    if (associated(grpSurfMSArea)) then
       nullify(grpSurfMSArea)
    endif

    if (associated(grpAccEpsRes)) then
       deallocate(grpAccEpsRes, stat=iDeallocateErr)
       if (iDeallocateErr /= 0) then
          ierror = 1
          return
       endif
    endif
    if (associated(grpAccChrgRes)) then
       deallocate(grpAccChrgRes, stat=iDeallocateErr)
       if (iDeallocateErr /= 0) then
          ierror = 1
          return
       endif
    endif
    if (associated(grpAccPhiOdd)) then
       deallocate(grpAccPhiOdd, stat=iDeallocateErr)
       if (iDeallocateErr /= 0) then
          ierror = 1
          return
       endif
    endif
    if (associated(grpAccKappaOdd)) then
       deallocate(grpAccKappaOdd, stat=iDeallocateErr)
       if (iDeallocateErr /= 0) then
          ierror = 1
          return
       endif
    endif
    if (associated(grpAccPhiEven)) then
       deallocate(grpAccPhiEven, stat=iDeallocateErr)
       if (iDeallocateErr /= 0) then
          ierror = 1
          return
       endif
    endif
    if (associated(grpAccKappaEven)) then
       deallocate(grpAccKappaEven, stat=iDeallocateErr)
       if (iDeallocateErr /= 0) then
          ierror = 1
          return
       endif
    endif
    if (associated(grpAccBoxBndry)) then
       deallocate(grpAccBoxBndry, stat=iDeallocateErr)
       if (iDeallocateErr /= 0) then
          ierror = 1
          return
       endif
    endif
    if (associated(gnpAccEpsIndx)) then
       deallocate(gnpAccEpsIndx, stat=iDeallocateErr)
       if (iDeallocateErr /= 0) then
          ierror = 1
          return
       endif
    endif
    if (associated(gnpAccChrgIndx)) then
       deallocate(gnpAccChrgIndx, stat=iDeallocateErr)
       if (iDeallocateErr /= 0) then
          ierror = 1
          return
       endif
    endif
    if (associated(gnpAccBoxBndryIndx)) then
       deallocate(gnpAccBoxBndryIndx, stat=iDeallocateErr)
       if (iDeallocateErr /= 0) then
          ierror = 1
          return
       endif
    endif
    if (associated(gnpAccBndryChrgIndx)) then
       deallocate(gnpAccBndryChrgIndx, stat=iDeallocateErr)
       if (iDeallocateErr /= 0) then
          ierror = 1
          return
       endif
    endif

    ! Used only in dpsetepsilon.F

    if (associated(npBndry1)) then
       DEALLOCATE(npBndry1, STAT=iDeallocateErr)
       if (iDeallocateErr /= 0) then
          ierror = 1
          return
       endif
    endif
    if (associated(npBndry2)) then
       DEALLOCATE(npBndry2, STAT=iDeallocateErr)
       if (iDeallocateErr /= 0) then
          ierror = 1
          return
       endif
    endif

    ! Used only in pb_spt.F
    if (associated(rpSurfVdw1)) then
       DEALLOCATE(rpSurfVdw1, STAT=iDeallocateErr)
       if (iDeallocateErr /= 0) then
          ierror = 1
          return
       endif
    endif
    if (associated(rpSurfVdw2)) then
       DEALLOCATE(rpSurfVdw2, STAT=iDeallocateErr)
       if (iDeallocateErr /= 0) then
          ierror = 1
          return
       endif
    endif
    if (associated(rpSurfMS1)) then
       DEALLOCATE(rpSurfMS1, STAT=iDeallocateErr)
       if (iDeallocateErr /= 0) then
          ierror = 1
          return
       endif
    endif
    if (associated(rpSurfMS2)) then
       DEALLOCATE(rpSurfMS2, STAT=iDeallocateErr)
       if (iDeallocateErr /= 0) then
          ierror = 1
          return
       endif
    endif
    if (associated(rpSurfMSArea1)) then
       DEALLOCATE(rpSurfMSArea1, STAT=iDeallocateErr)
       if (iDeallocateErr /= 0) then
          ierror = 1
          return
       endif
    endif
    if (associated(rpSurfMSArea2)) then
       DEALLOCATE(rpSurfMSArea2, STAT=iDeallocateErr)
       if (iDeallocateErr /= 0) then
          ierror = 1
          return
       endif
    endif

    return

  end subroutine deallocate_divpbVars

end module divpb_private
!********************************************************


!********************************************************
!  Use all modules (not for DivPB)
!--------------------------------------------------------
! Stupid way of dealing with modules, should be replaced and each file
! should have individual module lists

    module allmod
    use mfccmod
    use sizes
    use method
    use basis
    use params
    use molspec
    use gridpoints
    use calculated
    use geocnverg
    use SCRATCH
    use files
    use constants
    use ecpmod
    use quickdc
    use electrondensity 
    implicit none
    end module allmod
!********************************************************