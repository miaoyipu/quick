!*******************************************************
! scf(failed)
!-------------------------------------------------------
! Ed Brothers. November 27, 2001
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP
!-------------------------------------------------------
! this subroutine is to do scf job for restricted system
!
    subroutine scf(failed)
    use allmod
    implicit double precision(a-h,o-z)

    logical :: done,failed
    done=.false.

    !-----------------------------------------------------------------
    ! The purpose of this subroutine is to perform scf cycles.  At this
    ! point, X has been formed. The remaining steps are:
    ! 1)  Form operator matrix.
    ! 2)  Calculate O' = Transpose[X] O X
    ! 3)  Diagonalize O' to obtain C' and eigenvalues.
    ! 4)  Calculate C = XC'
    ! 5)  Form new density matrix.
    ! 6)  Check for convergence.
    !-----------------------------------------------------------------
    
    !-----------------------------------------------------------------
    ! Each location in the code that the step is occurring will be marked.
    ! The cycles stop when prms  is less than pmaxrms or when the maximum
    ! number of scfcycles has been reached.
    !-----------------------------------------------------------------
    
    jscf=0
    PRMS=1.d30

    !-----------------------------------------------------------------
    ! Alessandro GENONI 03/21/2007
    ! ECP integrals computation exploiting Alexander V. Mitin Subroutine
    ! Note: the integrals are stored in the array ecp_int that corresponds
    !       to the lower triangular matrix of the ECP integrals
    !-----------------------------------------------------------------  
    if (quick_method%ecp) then
      call ecpint
    end if

        if (quick_method%diisscf .and. .not. quick_method%divcon) call electdiis(jscf,PRMS)       ! normal scf
        if (quick_method%diisscf .and. quick_method%divcon) call electdiisdc(jscf,PRMS)     ! div & con scf
        jscf=jscf+1
        if (quick_method%debug)  call debug_SCF()
    return
    end subroutine scf