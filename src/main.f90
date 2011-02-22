! 
!************************************************************************
!                              QUICK                                   **
!                                                                      **
!                        Copyright (c) 2010                            **
!                Regents of the University of Florida                  **
!                       All Rights Reserved.                           **
!                                                                      **
!  This software provided pursuant to a license agreement containing   **
!  restrictions on its disclosure, duplication, and use. This software **
!  contains confidential and proprietary information, and may not be   **
!  extracted or distributed, in whole or in part, for any purpose      **
!  whatsoever, without the express written permission of the authors.  **
!  This notice, and the associated author list, must be attached to    **
!  all copies, or extracts, of this software. Any additional           **
!  restrictions set forth in the license agreement also apply to this  **
!  software.                                                           **
!************************************************************************
!
!  Cite this work as:
!  Miao,Y.: He, X.: Ayers,K; Brothers, E.: Merz,K. M. QUICK
!  University of Florida, Gainesville, FL, 2010
!************************************************************************
!
! Update Log:
! -2011-2-16  YIPU MIAO Add method module and try to convert the progs into OO.
! -2010-12-1  YIPU MIAO Add freq and hessian matrix calculation. But is broken from last author.
!                       not support MPI
! -2010-11-10 YIPU MIAO Complete MPI for MP2, only valid for HF MP2,not valid for Divcon MP2.
!                       Optimize code
! -2010-11-14 YIPU MIAO Add paralle option for HF and Divcon, valid for 1e,2e and 
!                       diag(divcon only) part
! -2010-7-27  YIPU MIAO add elimination step, improve accuracy. kill some bugs.
! -2010-7-23  YIPU MIAO first final version. Divcon works. Some tests are done. 
!                       But still buggy, especially the convergence problem.
! -2010-7-14  YIPU MIAO add HF SCF calculation. Divcon is partly work, but under test.
! -2010-5-24  YIPU MIAO add keyword "atombasis“ and "residuebasis" to specify fragment method. 
!                       if it's residue based, pdb file is needed, while atom based, no pdb file is needed
! -2010-5-21  YIPU MIAO new quick can now read pdb file as input. But hasn't reach the calculation part.
! -2010.5.15  YIPU MIAO Reorganize quick program. Ready to eliminate Xiao's personal variables.
! 
! -           XIOA HE   FOR + function , tighten cutoff
!                       MULTIPOLE one-electron only upto P orbital
!                       ifMM should be earlier than MFCC
! -           XIAO HE   Clean up main.f90 and right order to call scf.f90 and electdiisdc.f90 and 
!                       xdivided.f90, weird thing when I comment the read(9999)
! -           XIAO HE   be careful of HF=.true. and DFT=.true.
! -           Xiao HE   change these 3 subroutines for high efficiency of gradient
!                       ifort   -g -O3 -pg -traceback -c        hfgrad.f90
!                       ifort   -g -O3 -pg -traceback -c       2eshellopt.f90
!                       ifort   -g -O3 -pg -traceback -c       hrrsubopt.f90
! -           Xiao HE   change allocate(Xcoeff(jbasis,jbasis,0:3,0:3))
!**************************************************************************
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!   --- QUICK ---
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Main program of QUICK
! 
    program quick
    
    use allMod
    use divPB_Private, only: initialize_DivPBVars
    implicit none
    
    include 'mpif.h'
    
    logical :: failed = .false.         ! flag to indicates SCF fail or OPT fail 
    integer :: ierr                     ! return error info
    integer :: natomSave                ! save atom for globle use
    double precision, allocatable:: xyzSave(:,:)
                                        ! save cooridnates
                                        
    integer :: i                        ! integer temp varible


!*****************************************************************
! 1. The first thing that must be done is to initialize and prepare files
!------------------------------------------------------------------!
!
    ! Initial neccessary variables 
    call initialize1(ierr)
    
    !--------------------MPI/MASTER----------------------------------
    masterwork_initial: if (master) then ! Master Work: to initial the work
    !----------------End MPI/MASTER----------------------------------

      ! Read enviromental variables: QUICK_BASIS and ECPs
      ! those can be defined in ~/.bashrc
      call getenv("QUICK_BASIS",basisdir)
      call getenv("ECPs",ecpdir)

    !--------------------MPI/MASTER-----------------------------------
    endif masterwork_initial ! if (master)
    !--------------- End MPI/MASTER-----------------------------------


    !-------------------MPI/MASTER---------------------------------------
    masterwork_readInput: if (master) then

      ! Read argument, which is input file, usually ".in" file and prepare files:
      ! .out: output file
      ! .dmx: density matrix file
      ! .rst: coordinates file
      ! .pdb: PDB file (can be input if use PDB keyword)
      ! .cphf: CPHF file
      
      call getarg(1,inFileName)
      i = index(inFileName,'.')
      if(i .eq. 0) i = index(inFileName,' ')
      
      outFileName=inFileName(1:i-1)//'.out'
      dmxFileName=inFileName(1:i-1)//'.dmx'
      rstFileName=inFileName(1:i-1)//'.rst'
      CPHFFileName=inFileName(1:i-1)//'.cphf'
      pdbFileName=inFileName(1:i-1)//'.pdb'

      open(iOutFile,file=outFileName,status='unknown')
      
      ! At the beginning of output file, copyright information will be output first 
      call outputCopyright(ierr)
      
      ! Then output file information
      call PrtDate(iOutFile,'Task starts on:')
      write (iOutFile,'(" INPUT FILE:   ",A30)') inFileName
      write (iOutFile,'(" OUTPUT FILE:  ",A30)') outFileName

      ! Output MPI Information  
      if (bMPI .and. mpisize.eq.1) then
        bMPI=.false.
        write(iOutFile,'("WARNING: NODE=1, TURN OFF MPI)')
      endif
      
      if (bMPI) then
        write (iOutFile,*)
        write (iOutFile,'(" - MPI Enabled -")')
        write (iOutFile,'(" TOTAL PROCESSOR = ",i5)') mpisize
        write (iOutFile,'(" MASTER NAME     = ",A30)') pname
      endif
    
    endif masterwork_readInput ! master
    !--------------------End MPI/MASTER----------------------------------


!*****************************************************************
! 2. Next step is to read job and initial guess
!------------------------------------------------------------------
!
    !read job spec
    call readJob
    
    if (master) then
      ! Obtain Atom information
      if (quick_method%PDB) call readPDB(infile)
      
      call getAtoms()
      
      nAtomSave=natom
    endif

    !-------------------MPI/ALL NODES------------------------------------    
    if (bMPI) then
      call MPI_BCAST(natomsave,1,mpi_integer,0,MPI_COMM_WORLD,mpierror)
      call MPI_BCAST(iatomtype,1,mpi_integer,0,MPI_COMM_WORLD,mpierror)
      call MPI_BCAST(natom,1,mpi_integer,0,MPI_COMM_WORLD,mpierror)
      call MPI_BARRIER(MPI_COMM_WORLD,mpierror)
    endif
    !------------------END MPI/ALL NODES---------------------------------    

    !allocate essential variables
    call allocateAtoms()
    
    ! Then do inital guess
    call cpu_time(timer_begin%TIniGuess)
    
    ! a. SAD intial guess
    if (quick_method%SAD) call getMolSad()

    ! b. MFCC initial guess
    if (quick_method%MFCC) then
!       call mfcc
!       call getmolmfcc
    endif
    
    call cpu_time(timer_end%TIniGuess)

!*****************************************************************
! 3. Read Molecule Structure
!-----------------------------------------------------------------
!
    call getMol(nAtomSave)

!*****************************************************************
! 4. SCF single point calculation. DFT if wanted. If it is OPT job
!    ignore this part and go to opt part. We will get variationally determined Energy.
!-----------------------------------------------------------------
!

    ! if it is div&con method, begin fragmetation step, initial and setup
    ! div&con varibles
    if (quick_method%DIVCON) call inidivcon(nAtomSave)
    

    ! if it is not opt job, begin single point calculation
    if(.not.quick_method%opt)then
!      if(.NOT.PBSOL)then
!        call getEnergy(failed)
!      else
!        HF=.true.
!        DFT=.false.
!        call getEnergy(failed)
!      endif
!   else
        call g2eshell
        call schwarzoff
        call getEnergy(failed)
    endif

    if (failed) call quick_exit(iOutFile,1)


!*****************************************************************
! 5. OPT Geometry if wanted
!-----------------------------------------------------------------

    ! Geometry optimization. Currently, only cartesian version is 
    ! available. A improvement is in optimzenew, which is based on 
    ! internal coordinates, but is under coding.    
    if (quick_method%opt)  call optimize(failed)     ! Cartesian 
    if (failed) call quick_exit(iOutFile,1)          ! If geometry optimization fails

    ! Now at this point we have an energy and a geometry.  If this is
    ! an optimization job, we now have the optimized geometry.


!*****************************************************************
! 6. Other job option
!-----------------------------------------------------------------
    
    ! 6.a PB Solvant Model
    ! 11/03/2010 Blocked by Yiao Miao
!   if (PBSOL) then
!       call initialize_DivPBVars()
!       call pPBDriver(ierror)
!   endif

    ! 6.b MP2,2nd order Møller–Plesset perturbation theory
    if(quick_method%MP2) then
        if(.not. quick_method%DIVCON) then
           if (bMPI) then
             call mpi_calmp2    ! MPI-MP2
           else
             call calmp2()      ! none-MPI MP2
           endif
        else
            call calmp2divcon   ! DIV&CON MP2
        endif
    endif   !(quick_method%MP2)

    ! 6.c Freqency calculation and mode analysis
    ! note the analytical calculation is broken and needs to be fixed
    IF (quick_method%freq) THEN
        call calcHessian(failed)
        IF (failed) call quick_exit(iOutFile,1)     ! If Hessian matrix fails
        call frequency
    ENDIF

    ! 6.d clean spin for unrestricted calculation
    ! If this is an unrestricted calculation, check out the S^2 value to
    ! see if this is a reasonable wave function.  If not, modify it.
    
!    IF (quick_method%unrst) THEN
!        IF (quick_method%debug) call debugCleanSpin
!        IF (quick_method%unrst) call spinclean
!        IF (quick_method%debug) call debugCleanSpin
!    ENDIF

    if (master) then
        
        ! Convert Cartesian coordinator to internal coordinator
        if (quick_method%zmat) call zmake
        
        ! Calculate Dipole Moment
        call dipole
        
    endif
    
    ! Now at this point we have an energy and a geometry.  If this is
    ! an optimization job, we now have the optimized geometry.

!*****************************************************************
! 8.The final job is to output energy and many other infos
!-----------------------------------------------------------------
!
    call finalize(iOutFile,0)
    

    end program quick
