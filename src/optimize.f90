! IOPT to control the cycles
! Ed Brothers. August 18,2002.
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

subroutine optimize(failed)
  use allmod
  implicit double precision(a-h,o-z)

  logical :: done,diagco,failed
  character(len=1) cartsym(3)
  dimension W(3*natom*(2*MLBFGS+1)+2*MLBFGS)
  dimension coordsnew(natom*3),hdiag(natom*3),iprint(2)
  EXTERNAL LB2
  COMMON /LB3/MP,LP,GTOL,STPMIN,STPMAX

  logical lsearch,diis
  integer IMCSRCH,nstor,ndiis
  double precision gnorm,dnorm,diagter,safeDX,gntest,gtest,sqnpar,accls,oldGrad(3*natom),coordsold(natom*3)

  include "mpif.h"


  !---------------------------------------------------------
  ! This subroutine optimizes the geometry of the molecule. It has a
  ! variety of options that are enumerated in the text.  Please note
  ! that all of the methods in this subroutine presuppose the use of
  ! cartesian space for optimization.
  !---------------------------------------------------------

  cartsym(1) = 'X'
  cartsym(2) = 'Y'
  cartsym(3) = 'Z'
  done=.false.      ! flag to show opt is done
  diagco=.false.
  iprint(1)=-1
  iprint(2)=0
  EPS=1.d-6
  XTOL=1.d-11

  ! For right now, there is no way to adjust these and only analytical gradients
  ! are available.  This should be changed later.
  quick_method%analgrad=.true.

  ! Some varibles to determine the geometry optimization
!  geotest = .0054d0
    geotest=.0018d0
!  grmstest = .0036d0
    grmstest=.0012d0
!  gradtest =0.00135d0
    gradtest=0.00045d0
!  gnrmtest = .0009d0
    gnrmtest=0.00030d0

  stepmax=.1d0/0.529177249d0    ! Max change of one step
  IFLAG=0
  I=0

  do j=1,natom
    do k=1,3
       GRADIENT((j-1)*3+K)=0d0
    enddo
  enddo
    
  !------------- MPI/MASTER --------------------------------
  if (master) then
     call PrtAct(ioutfile,"Begin Optimization Job")  

     ! At the start of this routine we have a converged density matrix.
     ! Check to be sure you should be here.
     IF (natom == 1) THEN
        write (ioutfile,'(" ONE ATOM = NO OPTIMIZATION ")')
        return
     ENDIF

     IF (iopt < 0) THEN
        error=0.d0
        DO I=1,natom*3
           temp = GRADIENT(I)/5.D-4
           error = temp*temp+error
        ENDDO
        Write (ioutfile,'(" GRADIENT BASED ERROR =",F20.10)') error
     ENDIF
  endif

  !------------- END MPI/MASTER ----------------------------


  DO WHILE (I.le.iopt.and..not.done)
     I=I+1

     if (master) then
        call PrtAct(ioutfile,"Optimize for New Step")

        ! First let's print input geometry
        write(ioutfile,*)
        write(ioutfile,'(12("="))',advance="no")
        write(ioutfile,'(2x,"GEOMETRY FOR OPTIMIZATION STEP",I4," OUT OF ",I4,2x)',advance="no") I,iopt
        write(ioutfile,'(12("="))')
        write(ioutfile,*)
        write(ioutfile,'("GEOMETRY INPUT")')
        write(ioutfile,'("ELEMENT",6x,"X",9x,"Y",9x,"Z")')
        DO J=1,natom
           Write (ioutfile,'(2x,A2,6x,F7.4,3x,F7.4,3x,F7.4)') &
                symbol(iattype(J)),xyz(1,J)*0.529177249d0, &
                xyz(2,J)*0.529177249d0,xyz(3,J)*0.529177249d0
        ENDDO

!        Block temperorly, this is integralcutoff for opt step, 
!        should be larger than single point to save time
!        integralCutoff=1.0d0/(10.0d0**6.0d0)
!        Primlimit=1.0d0/(10.0d0**6.0d0)


        ! Save grad and coordinate for further use
        do j=1,natom
           do k=1,3
              oldGrad((j-1)*3+K)=GRADIENT((j-1)*3+K)
              coordsold((j-1)*3+K)=xyz(k,j)
           enddo
        enddo
     endif

     ! calculate energy first
     call g2eshell
     call schwarzoff
     call getEnergy(failed)
        
        ! Now we have several scheme to obtain gradient. For now, 
        ! only analytical gradient is available

        ! 11/19/2010 YIPU MIAO BLOCKED SOME SUBS.
        IF (quick_method%analgrad) THEN
           !            IF (quick_method%UNRST) THEN
           !                IF (quick_method%HF) call uhfgrad
           !                IF (quick_method%DFT) call uDFTgrad
           !                IF (quick_method%SEDFT) call uSEDFTgrad
           !            ELSE
           IF (quick_method%HF) then
                if (bMPI) then
                    call mpi_hfgrad
                else
                    call hfgrad
                endif
           endif
           !                IF (quick_method%DFT) call DFTgrad
           !                IF (quick_method%SEDFT) call SEDFTgrad
           !            ENDIF

        endif
        
    if (master) then

        IF (failed) return


        !-----------------------------------------------------------------------
        ! Copy current geometry into coordsnew. Fill the rest of the array w/ zeros.
        ! Also fill the rest of gradient with zeros.
        !-----------------------------------------------------------------------
        DO J=1,natom
           DO K=1,3
              coordsnew((J-1)*3 + K) = xyz(K,J)
           ENDDO
        ENDDO
        ! Now let's call LBFGS.
        call LBFGS(natom*3,MLBFGS,coordsnew,Etot,Gradient,DIAGCO,HDIAG,IPRINT,EPS,XTOL,W,IFLAG)

        lsearch=.false.
        diis=.true.
        imcsrch=0
        nstor=10
        dnorm=1.0d0
        diagterm=1.0e-4
        safeDX=0.10d0

        sqnpar=dsqrt(natom*3.0d0)
        gtest=0.5d0
        gntest = max(gtest*sqnpar*0.25d0,gtest)
        accls=0.0d0

        !-----------------------------------------------------------------------
        ! We have a new set of coordinates, copy it onto the xyz array,
        ! and get a new energy and set of gradients. Be sure to check stepsize.
        ! Also, assemble some of the test criterion.
        !-----------------------------------------------------------------------

        geomax = -1.d0
        georms = 0.d0
        DO J=1,natom
           DO K=1,3
              tempgeo =dabs(xyz(K,J)- coordsnew((J-1)*3 + K))

              ! If the change is too much, then we have to have small change to avoid error
              IF (tempgeo > stepmax) THEN
                 xyz(K,J) =  xyz(K,J)+(coordsnew((J-1)*3+K)-xyz(K,J))*stepmax/tempgeo
                 tempgeo = stepmax*0.529177249d0
              ELSE
                 tempgeo = tempgeo*0.529177249d0
                 xyz(K,J) = coordsnew((J-1)*3 + K)
              ENDIF

              ! Max geometry change
              geomax = max(geomax,tempgeo)
              ! geometry RMS
              georms = georms + tempgeo**2.d0
           ENDDO
        ENDDO

        gradmax = -1.d0
        gradnorm = 0.d0
        write (ioutfile,'(/," ANALYTICAL GRADIENT: ")')
        write (ioutfile,'(76("-"))')
        write (ioutfile,'(" VARIBLES",4x,"OLD_X",12x,"OLD_GRAD",8x,"NEW_GRAD",10x,"NEW_X")')
        write (ioutfile,'(76("-"))')
        DO Iatm=1,natom
           DO Imomentum=1,3
              ! Max gradient change
              gradmax = max(gradmax,dabs(GRADIENT((Iatm-1)*3+Imomentum)))
              ! Grad change normalization
              gradnorm = gradnorm + GRADIENT((Iatm-1)*3+Imomentum)**2.d0
              write (ioutfile,'(I5,A1,3x,F14.10,3x,F14.10,3x,F14.10,3x,F14.10)')Iatm,cartsym(imomentum), &
                   coordsold((Iatm-1)*3+Imomentum)*0.529177249d0,oldGrad((Iatm-1)*3+Imomentum), &
                   GRADIENT((Iatm-1)*3+Imomentum),xyz(Imomentum,Iatm)*0.529177249d0
           ENDDO
        ENDDO
        write(ioutfile,'(76("-"))')
        gradnorm = (gradnorm/dble(natom*3))**.5d0

        ! geometry RMS
        georms = (georms/dble(natom*3))**.5d0
        
        if (i.gt.1) then
            Write (ioutfile,'(" OPTIMZATION STATISTICS:")')
            Write (ioutfile,'(" ENERGY CHANGE           =",E20.10)') Etot-Elast
            Write (ioutfile,'(" MAXIMUM GEOMETRY CHANGE =",E20.10,"(REQUEST=",E12.5")")') geomax,geotest
            Write (ioutfile,'(" GEOMETRY CHANGE RMS     =",E20.10,"(REQUEST=",E12.5")")') georms,grmstest
            Write (ioutfile,'(" MAXIMUM GRADIENT ELEMENT=",E20.10,"(REQUEST=",E12.5")")') gradmax,gradtest
            Write (ioutfile,'(" GRADIENT NORM           =",E20.10,"(REQUEST=",E12.5")")') gradnorm,gnrmtest

            done = geotest.gt.geomax
            done = done.and.grmstest.gt.georms
            done = done.and.gradtest.gt.gradmax
            done = done.and.gnrmtest.gt.gradnorm
        else
            Write (ioutfile,'(" OPTIMZATION STATISTICS:")')
            Write (ioutfile,'(" MAXIMUM GRADIENT ELEMENT=",E20.10,"(REQUEST=",E20.10")")') gradmax,gradtest
            done = gradtest.gt.gradmax
            if (done) then
                Write (ioutfile,'(" NO SIGNAFICANT CHANGE, NO NEED TO OPTIMIZE. USE INITIAL GEOMETRY.")')
                do j=1,natom
                    do k=1,3
                        xyz(k,j)=coordsold((j-1)*3+K)
                    enddo
                enddo
            endif
        endif

        IF (done)  Write (ioutfile,'(/" GEOMETRY OPTIMIZED AFTER",i5," CYCLES")') i
        call PrtAct(ioutfile,"Finish Optimization for This Step")
        Elast = Etot
        ! If read is on, write out a restart file.
        IF (quick_method%readdmx) call wrtrestart
     endif
     
     !-------------- END MPI/MASTER --------------------

     ! we now have new geometry, and let other nodes know the new geometry
     if (bMPI)call MPI_BCAST(xyz,natom*3,mpi_double_precision,0,MPI_COMM_WORLD,mpierror)

     ! Notify every nodes if opt is done
     if (bMPI)call MPI_BCAST(done,1,mpi_logical,0,MPI_COMM_WORLD,mpierror)
  ENDDO

    
  if (master) then
     if (done) then
        Write (ioutfile,'("================ OPTIMIZED GEOMETRY INFORMATION ==============")')
     else
        write (ioutfile,*) "WARNING: REACH MAX OPT CYCLE BUT CANNOT GET OPTIMIZED GEOMETRY"
        write (ioutfile,*) "         BUT STILL OUTPUT THE GEOMETRY WE GET."
        Write (ioutfile,'("============= GEOMETRY INFORMATION (NOT OPTIMIZED) ===========")')
     endif
     write (ioutfile,*)
     Write (ioutfile,'(" OPTIMIZED GEOMETRY IN CARTESIAN")')
     write (ioutfile,'(" ELEMENT",6x,"X",9x,"Y",9x,"Z")')
     
     DO I=1,natom
        Write (ioutfile,'(2x,A2,6x,F7.4,3x,F7.4,3x,F7.4)') &
             symbol(iattype(I)),xyz(1,I)*0.529177249d0, &
             xyz(2,I)*0.529177249d0,xyz(3,I)*0.529177249d0
     ENDDO
     
     write(ioutfile,*)
     write (ioutfile,'(" FORCE")')
     write (ioutfile,'(" ELEMENT",6x, "X",9x,"Y",9x,"Z")')
     do i=1,natom
        write(ioutfile,'(2x,A2,6x,F7.4,3x,F7.4,3x,F7.4)') &
            symbol(iattype(I)),-GRADIENT((i-1)*3+1)*0.529177249d0, &
            -GRADIENT((i-1)*3+2)*0.529177249d0,-GRADIENT((i-1)*3+3)*0.529177249d0
     enddo
    
     write (ioutfile,*)
     write (ioutfile,'(" MINIMIZED ENERGY=",F15.10)') Etot
     Write (ioutfile,'("===============================================================")')
    

     call PrtAct(ioutfile,"Finish Optimization Job")
  endif

  return
end subroutine optimize


