!
!	quick_method_module.f90
!	new_quick
!
!	Created by Yipu Miao on 2/18/11.
!	Copyright 2011 University of Florida. All rights reserved.
!
    module quick_method_module
    implicit none
    
    type quick_method_type

! the first some elements of namelist is for the QM method that is going to use 
        logical :: HF =  .false.       ! HF
        logical :: DFT =  .false.      ! DFT
        logical :: MP2 =  .false.      ! MP2
        logical :: B3LYP = .false.     ! B3LYP
        logical :: BLYP = .false.      ! BLYP
        logical :: BPW91 = .false.     ! BPW91
        logical :: MPW91LYP = .false.  ! MPW91LYP
        logical :: MPW91PW91 = .false. ! MPW91PW91
        logical :: SEDFT = .false.     ! Semi-Empirical DFT
        logical :: PBSOL = .false.     ! PB Solvent
        logical :: UNRST =  .false.    ! Unrestricted

        logical :: debug =  .false.    ! debug mode
        logical :: readDMX =  .false.  ! flag to read density matrix
        logical :: diisSCF =  .false.  ! DIIS SCF
        logical :: prtGap =  .false.   ! flag to print HOMO-LUMO gap
        logical :: opt =  .false.      ! optimization
        logical :: analGrad =  .false. ! Analytical Gradient
        logical :: analHess =  .false. ! Analytical Hessian Matrix

        logical :: diisOpt =  .false.  ! DIIS Optimization
        logical :: core =  .false.     ! Add core
        logical :: annil =  .false.    ! Annil Spin Contamination
        logical :: freq =  .false.     ! Frenquency calculation
        logical :: zmat = .false.      ! Z-matrix
!
        logical :: ecp                 ! ECP
        logical :: custECP             ! Custom ECP
        logical :: printEnergy = .true.! Print Energy each cycle
        logical :: fFunXiao            ! If f orbitial is contained

!  
        logical :: calcDens = .false.    ! calculate density
        logical :: calcDensLap = .false. ! calculate density lap
        logical :: writePMat = .false.   ! Output density matrix
 
        logical :: extCharges = .false.  ! external charge
        logical :: PDB = .false.         ! PDB input
        logical :: SAD = .true.          ! SAD initial guess
        logical :: FMM = .false.         ! Fast Multipole
        logical :: DIVCON = .false.      ! Div&Con
        
        integer :: ifragbasis = 1        ! =2.residue basis,=1.atom basis(DEFUALT),=3 non-h atom basis
        integer :: iSG = 1               ! =0. SG0, =1. SG1(DEFAULT)
        logical :: MFCC = .false.        ! MFCC
    end type quick_method_type
    
    type (quick_method_type) quick_method
    
    contains
        
        !------------------------
        ! Broadcast quick_method
        !------------------------
        subroutine broadcast_quick_method(self)
            use quick_MPI_module
            implicit none
        
            type(quick_method_type) self
        
            include 'mpif.h'

            call MPI_BARRIER(MPI_COMM_WORLD,mpierror)
            call MPI_BCAST(self%HF,1,mpi_logical,0,MPI_COMM_WORLD,mpierror)
            call MPI_BCAST(self%DFT,1,mpi_logical,0,MPI_COMM_WORLD,mpierror) 
            call MPI_BCAST(self%MP2,1,mpi_logical,0,MPI_COMM_WORLD,mpierror) 
            call MPI_BCAST(self%B3LYP,1,mpi_logical,0,MPI_COMM_WORLD,mpierror)
            call MPI_BCAST(self%BLYP,1,mpi_logical,0,MPI_COMM_WORLD,mpierror)
            call MPI_BCAST(self%BPW91,1,mpi_logical,0,MPI_COMM_WORLD,mpierror)
            call MPI_BCAST(self%MPW91LYP,1,mpi_logical,0,MPI_COMM_WORLD,mpierror)
            call MPI_BCAST(self%MPW91PW91,1,mpi_logical,0,MPI_COMM_WORLD,mpierror)
            call MPI_BCAST(self%PBSOL,1,mpi_logical,0,MPI_COMM_WORLD,mpierror)
            call MPI_BCAST(self%UNRST,1,mpi_logical,0,MPI_COMM_WORLD,mpierror)
            call MPI_BCAST(self%debug,1,mpi_logical,0,MPI_COMM_WORLD,mpierror)
            call MPI_BCAST(self%readDMX,1,mpi_logical,0,MPI_COMM_WORLD,mpierror)
            call MPI_BCAST(self%diisSCF,1,mpi_logical,0,MPI_COMM_WORLD,mpierror)
            call MPI_BCAST(self%prtGap,1,mpi_logical,0,MPI_COMM_WORLD,mpierror)
            call MPI_BCAST(self%opt,1,mpi_logical,0,MPI_COMM_WORLD,mpierror)
            call MPI_BCAST(self%analGrad,1,mpi_logical,0,MPI_COMM_WORLD,mpierror)
            call MPI_BCAST(self%analHess,1,mpi_logical,0,MPI_COMM_WORLD,mpierror)
            call MPI_BCAST(self%diisOpt,1,mpi_logical,0,MPI_COMM_WORLD,mpierror)
            call MPI_BCAST(self%core,1,mpi_logical,0,MPI_COMM_WORLD,mpierror)
            call MPI_BCAST(self%annil,1,mpi_logical,0,MPI_COMM_WORLD,mpierror)
            call MPI_BCAST(self%freq,1,mpi_logical,0,MPI_COMM_WORLD,mpierror)
            call MPI_BCAST(self%SEDFT,1,mpi_logical,0,MPI_COMM_WORLD,mpierror)
            call MPI_BCAST(self%Zmat,1,mpi_logical,0,MPI_COMM_WORLD,mpierror)
            call MPI_BCAST(self%ecp,1,mpi_logical,0,MPI_COMM_WORLD,mpierror)
            call MPI_BCAST(self%custECP,1,mpi_logical,0,MPI_COMM_WORLD,mpierror)
            call MPI_BCAST(self%printEnergy,1,mpi_logical,0,MPI_COMM_WORLD,mpierror)
            call MPI_BCAST(self%fFunXiao,1,mpi_logical,0,MPI_COMM_WORLD,mpierror)
            call MPI_BCAST(self%calcDens,1,mpi_logical,0,MPI_COMM_WORLD,mpierror)
            call MPI_BCAST(self%calcDensLap,1,mpi_logical,0,MPI_COMM_WORLD,mpierror)
            call MPI_BCAST(self%writePMat,1,mpi_logical,0,MPI_COMM_WORLD,mpierror)
            call MPI_BCAST(self%extCharges,1,mpi_logical,0,MPI_COMM_WORLD,mpierror)
            call MPI_BCAST(self%PDB,1,mpi_logical,0,MPI_COMM_WORLD,mpierror)
            call MPI_BCAST(self%SAD,1,mpi_logical,0,MPI_COMM_WORLD,mpierror)
            call MPI_BCAST(self%FMM,1,mpi_logical,0,MPI_COMM_WORLD,mpierror)
            call MPI_BCAST(self%DIVCON,1,mpi_logical,0,MPI_COMM_WORLD,mpierror)
            call MPI_BCAST(self%MFCC,1,mpi_logical,0,MPI_COMM_WORLD,mpierror)
            call MPI_BCAST(self%ifragbasis,1,mpi_integer,0,MPI_COMM_WORLD,mpierror)
            call MPI_BCAST(self%iSG,1,mpi_integer,0,MPI_COMM_WORLD,mpierror)
        
        end subroutine broadcast_quick_method
        
        
        
        !------------------------
        ! print quick_method
        !------------------------
        subroutine print_quick_method(self,io)
            implicit none
            integer io
            type(quick_method_type) self
                
            write(io,'(" ====== JOB CARD =====")')
            if (self%HF) then
                write(io,'("METHOD = HATREE FOCK")')
            else if (self%MP2) then
                write(io,'("METHOD = SECOND ORDER PERTURBATION THEORY")')
            else if (self%DFT) then
                write(io,'("METHOD = DENSTITY FUNCTIONAL THEORY")')
            else if (self%SEDFT) then
                write(io,'("METHOD = SEMI-EMPIRICAL DENSTITY FUNCTIONAL THEORY")')
            endif
        
            if (self%PDB)       write(io,'("PDB INPUT")')
            if (self%MFCC)      write(io,'("MFCC INITIAL GUESS")')
            if (self%SAD)       write(io,'("SAD INITAL GUESS")')
            if (self%DIVCON)    write(io,'("DIV & CON METHOD")')
            if (self%FMM)       write(io,'("USING FAST MULTIPOLE METHOD")')
            if (self%UNRST)     write(io,'("UNRESTRICTED SYSTEM")')
            if (self%annil)     write(io,'("ANNIHILATE SPIN CONTAMINAT")')
            if (self%PBSOL)     write(io,'("SOLVATION MODEL = PB")')
            if (self%printEnergy) write(io,'("PRINT ENERGY EVERY CYCLE")')
            if (self%readDMX)   write(io,'("READ DENSITY MATRIX FROM FILE")')
            if (self%diisSCF)   write(io,'("USE DIIS SCF")')
            if (self%prtGap)    write(io,'("PRINT HOMO-LUMO GAP")')
            if (self%opt)       write(io,'("GEOMETRY OPTIMIZATION")')
            if (self%freq)       write(io,'("FREQENCY CALCULATION")')
            if (self%zmat)      write(io,'("Z-MATRIX CONSTRUCTION")')
            if (self%ecp)       write(io,'("ECP BASIS SET")')
            if (self%custECP)   write(io,'("CUSTOM ECP BASIS SET")')
            if (self%writePMat) write(io,'("WRITE DENSITY MATRIX TO FILE")')
            if (self%extCharges)write(io,'("EXTERNAL CHARGES")')
            if (self%core)      write(io,'("SUM INNER ELECTRONS INTO CORE")')
            if (self%diisOpt)   write(io,'("USE DIIS FOR GEOMETRY OPTIMIZATION")')
            if (self%debug)     write(io,'("DEBUG MODE")') 
        
            if (self%DFT) then
                if (self%iSG .eq. 0) write(io,'("STANDARD GRID = SG0")')
                if (self%iSG .eq. 1) write(io,'("STANDARD GRID = SG1")')
            endif
               
            if (self%opt) then         
                if (self%analGrad)  then
                    write(io,'("ANALYTICAL GRADIENT")')
                else
                    write(io,'("NUMERICAL GRADIENT")')
                endif    

                if (self%analHess)  then
                    write(io,'("ANALYTICAL HESSIAN MATRIX")')
                else
                    write(io,'("NUMERICAL HESSIAN MATRIX")')
                endif
            endif
        
            if (self%DFT) then
                if (self%B3LYP) then
                    write(io,'("DENSITY FUNCTIONAL = B3LYP")')
                elseif(self%BLYP) then
                    write(io,'("DENSITY FUNCTIONAL = B3YP")')
                elseif(self%BPW91) then
                    write(io,'("DENSITY FUNCTIONAL = BPW91")')
                elseif(self%MPW91LYP) then
                    write(io,'("DENSITY FUNCTIONAL = MPW91LYP")')
                elseif(self%MPW91PW91) then
                    write(io,'("DENSITY FUNCTIONAL = MPW91PW91")')
                endif
            endif
        
            if (self%DIVCON) then
                if (self%ifragbasis .eq. 1) then
                    write(io,'("DIV AND CON ON ATOM BASIS")')
                elseif (self%ifragbasis .eq. 2) then
                    write(iO,'("DIV AND CON ON RESIDUE BASIS")')
                elseif (self%ifragbasis .eq. 3) then
                    write(io,'("DIV AND CON ON NON HYDROGEN ATOM BASIS")')
                else
                    write(io,'("DIV AND CON ON ATOM BASIS (BY DEFAULT)")')
                endif
            endif
            
        end subroutine print_quick_method
        
        
        
        !------------------------
        ! read quick_method
        !------------------------
        subroutine read_quick_method(self,keywd)
            implicit none
            character(len=200) :: keyWD
            type (quick_method_type) self
        
            call upcase(keyWD,200)
            if (index(keyWD,'PDB').ne. 0)       self%PDB=.true.
            if (index(keyWD,'MFCC').ne.0)       self%MFCC=.true.
            if (index(keyWD,'FMM').ne.0)        self%FMM=.true.
            if (index(keyWD,'MP2').ne.0)        self%MP2=.true. 
            if (index(keyWD,'HF').ne.0)         self%HF=.true.    
            if (index(keyWD,'DFT').ne.0)        self%DFT=.true.
            if (index(keyWD,'SEDFT').ne.0)      self%SEDFT=.true.    
            if (index(keyWD,'PBSOL').ne.0)      self%PBSOL=.true.
            if (index(keyWD,'ANNIHILATE').ne.0) self%annil=.true.
            if (index(keyWD,'B3LYP').ne.0)      self%B3LYP=.true.
            if (index(keyWD,'BLYP').ne.0)       self%BLYP=.true.
            if (index(keyWD,'BPW91').ne.0)      self%BPW91=.true.
            if (index(keyWD,'MPW91LYP').ne.0)   self%MPW91LYP=.true.
            if (index(keyWD,'MPW91PW91').ne.0)  self%MPW91PW91=.true.
            if (index(keyWD,'CORE').ne.0)       self%CORE=.true.
            if (index(keyWD,'OPT').ne.0)        self%opt=.true.
            if (index(keyWD,'DIIS-OPTIMIZE').ne.0)self%diisOpt=.true.
            if (index(keyWD,'GAP').ne.0)        self%prtGap=.true.
            if (index(keyWD,'GRAD').ne.0)       self%analGrad=.true.
            if (index(keyWD,'HESSIAN').ne.0)    self%analHess=.true.
            if (index(keyWD,'FREQ').ne.0)       self%freq=.true.
            if (index(keywd,'DEBUG').ne.0)      self%debug=.true.
            if (index(keyWD,'READ').ne.0)       self%readDMX=.true.
            if (index(keyWD,'ZMAKE').ne.0)      self%zmat=.true.
            if (index(keyWD,'WRITE').ne.0)      self%writePMat=.true.
            if (index(keyWD,'EXTCHARGES').ne.0) self%EXTCHARGES=.true.
        
        
            if (index(keywd,'DIVCON') .ne. 0) then
                self%divcon = .true.
                if (index(keywd,'ATOMBASIS') /= 0) then
                    self%ifragbasis=1
                else if (index(keywd,'RESIDUEBASIS') /= 0) then
                    self%ifragbasis=2
                else if (index(keywd,'NHAB') /= 0) then
                    self%ifragbasis=3
                else
                    self%ifragbasis=1
                endif
            endif
        
            if (self%DFT) then
                if (index(keyWD,'SG0').ne.0) then
                    self%iSG=0
                else
                    self%iSG=1
                endif
            endif
        
            self%printEnergy=.true.
            self%sad=.true.
            self%diisSCF=.true.
        
            if (index(keyWD,'ECP').ne.0)  then
                self%ECP=.true.
                if (index(keyWD,'CUSTOM').ne.0)  self%custECP=.true.
            endif
        
            if (index(keyWD,'USEDFT').ne.0) then
                self%SEDFT=.true.    
                self%UNRST=.true.
            endif
        
            if (index(keyWD,'UHF').ne.0) then
                self%HF=.true.
                self%UNRST=.true.
            endif
        
            if (index(keyWD,'UDFT').ne.0) then
                self%DFT=.true.
                self%UNRST=.true.
            endif
        end subroutine read_quick_method
        
        subroutine set_quick_method(self)
        
            implicit none
            type(quick_method_type) self        
            self%HF =  .false.       ! HF
            self%DFT =  .false.      ! DFT
            self%MP2 =  .false.      ! MP2
            self%B3LYP = .false.     ! B3LYP
            self%BLYP = .false.      ! BLYP
            self%BPW91 = .false.     ! BPW91
            self%MPW91LYP = .false.  ! MPW91LYP
            self%MPW91PW91 = .false. ! MPW91PW91
            self%SEDFT = .false.     ! Semi-Empirical DFT
            self%PBSOL = .false.     ! PB Solvent
            self%UNRST =  .false.    ! Unrestricted

            self%debug =  .false.    ! debug mode
            self%readDMX =  .false.  ! flag to read density matrix
            self%diisSCF =  .false.  ! DIIS SCF
            self%prtGap =  .false.   ! flag to print HOMO-LUMO gap
            self%opt =  .false.      ! optimization
            self%analGrad =  .false. ! Analytical Gradient
            self%analHess =  .false. ! Analytical Hessian Matrix

            self%diisOpt =  .false.  ! DIIS Optimization
            self%core =  .false.     !
            self%annil =  .false.    !
            self%freq =  .false.     ! Frenquency calculation
            self%zmat = .false.      ! Z-matrix
            self%ecp = .false.       ! ECP
            self%custECP = .false.   ! Custom ECP
            self%printEnergy = .true.! Print Energy each cycle
            self%fFunXiao = .false.            ! If f orbitial is contained 
            self%calcDens = .false.    ! calculate density
            self%calcDensLap = .false. ! calculate density lap
            self%writePMat = .false.   ! Output density matrix
            self%extCharges = .false.  ! external charge
            self%PDB = .false.         ! PDB input
            self%SAD = .true.          ! SAD initial guess
            self%FMM = .false.         ! Fast Multipole
            self%DIVCON = .false.      ! Div&Con
        
            self%ifragbasis = 1        ! =2.residue basis,=1.atom basis(DEFUALT),=3 non-h atom basis
            self%iSG = 1               ! =0. SG0, =1. SG1(DEFAULT)
            self%MFCC = .false.        ! MFCC
        
        end subroutine set_quick_method
        
        subroutine check_quick_method(self,io)
            implicit none
            type(quick_method_type) self
            integer io
            
            ! If MP2, then set HF as default
            if (self%MP2) then
                self%HF = .true.
                self%DFT = .false.
            endif
            
            
            ! OPT not available for MP2
            if (self%MP2 .and. self%OPT) then
                write(io,*) " WARN: GEOMETRY OPTIMIZAION IS NOT AVAILABLE WITH MP2, WILL DO MP2 SINGLE POINT ONLE"
                self%OPT = .false.
            endif

            ! OPT not available for other DFT except BLYP            
            if(self%DFT.and. self%OPT .and. (.NOT. self%BLYP))then
                write(iO,*) "ERROR: GEOMETRY OPTIMIZATION is only available with HF, DFT/BLYP"
                self%OPT = .false.
            endif

        end subroutine check_quick_method
        
        
    end module quick_method_module