#include "config.h"
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

#ifdef MPI
    include "mpif.h"
#endif

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
        quick_qm_struct%Ecore=0.d0      ! atom-extcharge and atom-atom replusion
        quick_qm_struct%ECharge=0d0     ! extcharge-extcharge interaction

        if (natom > 1) then
            do I=1,natom+quick_molspec%nextatom
                do J=I+1,natom+quick_molspec%nextatom
                    if(i<=natom .and. j<=natom)then
                        distance =((xyz(1,I)-xyz(1,J))**2.d0 &
                                 + (xyz(2,I)-xyz(2,J))**2.d0 &
                                 + (xyz(3,I)-xyz(3,J))**2.d0)**.5d0
                        quick_qm_struct%Ecore=quick_molspec%chg(I)*quick_molspec%chg(J)/distance+quick_qm_struct%Ecore
                    elseif(i<=natom .and. j>natom)then
                        distance =((xyz(1,I)-quick_molspec%extxyz(1,J-natom))**2.d0 &
                                 + (xyz(2,I)-quick_molspec%extxyz(2,J-natom))**2.d0 &
                                 + (xyz(3,I)-quick_molspec%extxyz(3,J-natom))**2.d0)**.5d0
                        quick_qm_struct%Ecore=quick_molspec%chg(I)*quick_molspec%extchg(J-natom)/distance+quick_qm_struct%Ecore
                    elseif(i>natom .and. j>natom)then
                        distance =((quick_molspec%extxyz(1,I-natom)-quick_molspec%extxyz(1,J-natom))**2.d0 &
                                 + (quick_molspec%extxyz(2,I-natom)-quick_molspec%extxyz(2,J-natom))**2.d0 &
                                 + (quick_molspec%extxyz(3,I-natom)-quick_molspec%extxyz(3,J-natom))**2.d0)**.5d0
                        quick_qm_struct%ECharge=quick_qm_struct%ECharge + &
                                  quick_molspec%extchg(I-natom)*quick_molspec%extchg(J-natom)/distance
                    endif 
                enddo
            enddo
        endif
        
    endif
    !-----------------------------------------------------------------
    ! Converge the density matrix.
    !-----------------------------------------------------------------
#ifdef MPI
    !-------------- MPI / ALL NODES ----------------------------------
    if (bMPI) then
        call MPI_BCAST(quick_qm_struct%s,nbasis*nbasis,mpi_double_precision,0,MPI_COMM_WORLD,mpierror)
        call MPI_BCAST(quick_qm_struct%x,nbasis*nbasis,mpi_double_precision,0,MPI_COMM_WORLD,mpierror)
        call MPI_BCAST(quick_qm_struct%Ecore,1,mpi_double_precision,0,MPI_COMM_WORLD,mpierror)
    endif
    !-------------- END MPI / ALL NODES ------------------------------
#endif
    if (quick_method%UNRST) then
        call uscf(failed)
    else
        call scf(failed)
    endif

    !--------------- MPI/MASTER --------------------------
    if (master) then
    
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
        quick_qm_struct%Eelvac=quick_qm_struct%Eel
        if (quick_method%extcharges) then
            quick_qm_struct%Etot = quick_qm_struct%Etot + quick_qm_struct%ECharge
        endif
        quick_qm_struct%Etot = quick_qm_struct%Eel + quick_qm_struct%Ecore

        write (ioutfile,'("ELECTRONIC ENERGY    =",F16.9)') quick_qm_struct%Eel
        write (ioutfile,'("CORE_CORE REPULSION  =",F16.9)') quick_qm_struct%Ecore
        if (quick_method%extcharges) then
            write (ioutfile,'("EXT CHARGE REPULSION =",F16.9)') quick_qm_struct%ECharge
        endif
        write (ioutfile,'("TOTAL ENERGY         =",F16.9)') quick_qm_struct%Etot        
        call prtact(ioutfile,"End Energy calculation")

        call flush(ioutfile)
    endif
    !--------------- END MPI/MASTER ----------------------
    
    end subroutine getenergy
