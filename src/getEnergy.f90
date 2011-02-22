!*******************************************************
! getEnergy(failed)
!-------------------------------------------------------
! Ed Brothers. 08/15/02
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP
!-------------------------------------------------------
! This subroutine calculates and ouptus the energy.
!
    subroutine getEnergy(failed)
    use allMod
    implicit double precision(a-h,o-z)

    include "mpif.h"
    logical :: failed
    
    !-----------------------------------------------------------------
    ! Build a transformation matrix X.
    !-----------------------------------------------------------------
    if (master) then
        call PrtAct(ioutfile,"Begin Energy Calculation")
        call fullX
        ! Construct Div & Con matrices, Overlap,X, and PDC
        if (quick_method%DivCon) then
            call DivideS
            call DivideX
            call PDCDivided
        endif   

        !-----------------------------------------------------------------
        !Classical Nuclear-Nuclear interaction energy
        !-----------------------------------------------------------------
        Ecore=0.d0      ! atom-extcharge and atom-atom replusion
        Echarge=0d0     ! extcharge-extcharge interaction
        if(.not.quick_method% extcharges) nextatom=0
        if (natom > 1) then
            do I=1,natom+nextatom
                do J=I+1,natom+nextatom
                    if(i<=natom .and. j<=natom)then
                        distance =((xyz(1,I)-xyz(1,J))**2.d0 &
                        + (xyz(2,I)-xyz(2,J))**2.d0 &
                        + (xyz(3,I)-xyz(3,J))**2.d0)**.5d0
                        Ecore=chg(I)*chg(J)/distance+Ecore
                    elseif(i<=natom .and. j>natom)then
                        distance =((xyz(1,I)-extxyz(1,J-natom))**2.d0 &
                        + (xyz(2,I)-extxyz(2,J-natom))**2.d0 &
                        + (xyz(3,I)-extxyz(3,J-natom))**2.d0)**.5d0
                        Ecore=chg(I)*extchg(J-natom)/distance+Ecore
                    elseif(i>natom .and. j>natom)then
                        distance =((extxyz(1,I-natom)-extxyz(1,J-natom))**2.d0 &
                        + (extxyz(2,I-natom)-extxyz(2,J-natom))**2.d0 &
                        + (extxyz(3,I-natom)-extxyz(3,J-natom))**2.d0)**.5d0
                        Echarge=extchg(I-natom)*extchg(J-natom)/distance+Echarge
                    endif 
                enddo
            enddo
        endif
        
    endif
    !-----------------------------------------------------------------
    ! Converge the density matrix.
    !-----------------------------------------------------------------

    if (bMPI) then
        call MPI_BCAST(Smatrix,nbasis*nbasis,mpi_double_precision,0,MPI_COMM_WORLD,mpierror)
        call MPI_BCAST(X,nbasis*nbasis,mpi_double_precision,0,MPI_COMM_WORLD,mpierror)
        call MPI_BCAST(Ecore,1,mpi_double_precision,0,MPI_COMM_WORLD,mpierror)
    endif

    if (quick_method%UNRST) then
        call uscf(failed)
    else
        call scf(failed)
    endif

    !--------------- MPI/MASTER --------------------------
    if (master) then
        !-----------------------------------------------------------------
        ! If we are saving the density matrix, do so now.
        !-----------------------------------------------------------------
        if (quick_method%readDMX) then
            call readDMXfile()
        endif
        !---------------------------------------------------------------
        ! Now that we have a converged density matrix, it is time to
        ! calculate the energy.  First it is necessary to obtain the
        ! interatomic repulsions.
        !---------------------------------------------------------------
    
        !-----------------------------------------------------------------
        ! Now calculate the energy for PB Sol
        !-----------------------------------------------------------------
        !
        ! Blocked by Yipu Miao
        !
        if(quick_method%PBSOL)then
            if (quick_method%UNRST) then
!       if (quick_method%HF) call UHFEnergy
!       if (quick_method%DFT) call uDFTEnergy
!        if (quick_method%SEDFT) call uSEDFTEnergy
            else
!        if (quick_method%HF) call HFEnergy
!        if (quick_method%DFT) call DFTenergy
!        if (quick_method%SEDFT) call SEDFTenergy
            endif
        endif

        !-----------------------------------------------------------------
        ! Output Energy
        !-----------------------------------------------------------------

        write (ioutfile,'("ELECTRONIC ENERGY  =",F16.9)') Eel
        Eelvac=Eel
        write (ioutfile,'("CORE_CORE REPULSION=",F16.9)') Ecore
        Etot = Eel + Ecore
        if (quick_method%extcharges) then
            write (ioutfile,'("EXTERNAL CHARGE REPULSION=",F16.9)') Echarge
            Etot = Etot + Echarge
        endif
        write (ioutfile,'("TOTAL ENERGY       =",F16.9)') Etot
        
        call prtact(ioutfile,"End Energy calculation")

        call flush(ioutfile)
    endif
    !--------------- END MPI/MASTER ----------------------
    
    end subroutine getenergy

!*******************************************************
! readDMXfile()
!-------------------------------------------------------
! This subroutine is to read density matrix from DMX file
!
    subroutine readDMXfile()
    use allMod
    implicit double precision(a-h,o-z)
    integer iRead,jRead
    
    open(idmxFile,file=dmxFileName,status='unknown')
    do iRead=1,nbasis
        do jRead=iRead,nbasis
            write (idmxFile,'(I4,I4,3x,E15.8)') &
                iRead,jRead,DENSE(jRead,iRead)
        enddo
    enddo
    
    if (quick_method%unrst) then
        write (idmxfile,'("BETA")')
        do iRead=1,nbasis
            do jRead=iRead,nbasis
                write (idmxfile,'(I4,I4,3x,E15.8)') &
                    iRead,jRead,DENSEB(jRead,iRead)
            enddo
        enddo
    endif
    write (idmxFile,'("END")')
    close(idmxFile)

    end subroutine readDMXfile
