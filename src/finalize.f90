! Xiao HE deallocate 07/17/2007
! Ken Ayers 05/26/04
! Subroutines for allocation of the various matricies in quick
! These routines are not the ideal way to deal with allocation of
! variables.  Large sized arrays should only be allocated when
! they are needed.  Eventually someone will deal with this.

subroutine deallocateall
  use allmod
  implicit double precision(a-h,o-z)

  integer i

  deallocate(Yxiao)
  deallocate(Yxiaotemp)
  deallocate(Yxiaoprim)
  deallocate(attraxiao)
  deallocate(attraxiaoopt)
  deallocate(Ycutoff)
  deallocate(cutmatrix)
  deallocate(allerror)
  deallocate(alloperator)
  !  deallocate(debug1)
  !  deallocate(debug2)
  deallocate(kstart)
  deallocate(katom)
  deallocate(ktype)
  deallocate(kprim)
  ! print*,"problem1"
  deallocate(Qnumber)
  ! print*,"here"
  deallocate(Qstart)
  deallocate(Qfinal)
  ! print*,"middle"
  deallocate(Qsbasis)
  deallocate(Qfbasis)
  deallocate(ksumtype)
  ! print*,"probelm2" 
  deallocate(KLMN)
  deallocate(cons)
  deallocate(gccoeff)
  deallocate(gcexpo)
  deallocate(gcexpomin)
  deallocate(aex)
  deallocate(gcs)
  deallocate(gcp)
  deallocate(gcd)
  deallocate(gcf)
  deallocate(gcg)

  deallocate(itype)
  deallocate(ncenter)
  deallocate(ncontract)

  deallocate(sigrad2)

  deallocate(Smatrix)
  deallocate(X)
  deallocate(O)
  deallocate(CO)
  deallocate(COB)
  deallocate(VEC)
  deallocate(DENSE)
  deallocate(DENSEB)
  deallocate(DENSEOLD)
  deallocate(DENSESAVE)
  deallocate(Osave)
  deallocate(Osavedft)
  deallocate(V2)
  deallocate(E)
  deallocate(EB)
  deallocate(idegen)
  deallocate(Uxiao)

  deallocate(hold)
  deallocate(hold2)


  deallocate(aexp)
  deallocate(dcoeff)
  deallocate(gauss)
  ! do i=1,nbasis
  !     deallocate(gauss)
  !     deallocate(gauss)
  ! enddo

  ! do i=1,nbasis
  !     deallocate(gauss%aexp(3))
  !     deallocate(gauss%dcoeff(3))
  ! enddo

  ! deallocate(Apri)
  ! deallocate(Kpri)
  ! deallocate(cutprim)
  ! deallocate(Ppri)
  ! deallocate(Xcoeff)

end subroutine deallocateall


!----------------------
! Finialize programs
!----------------------
subroutine finalize(io,status)
    use allmod
    implicit none
    integer io      !output final info and close this unit
    integer status  !exit status: 1-error 0-normal
    
    ! Deallocate all variables
    call deallocateall
    
    ! stop timer and output them
    call cpu_time(timer_end%TTotal)
    call timer_output(io)
    
    !-------------------MPI/MASTER---------------------------------------
    master_finalize:if (master) then
        if (status /=0) then
            call PrtDate(io,'Error Occured. Task Finished on:')
        else
            call PrtDate(io,'Normal Terminatation. Task Finished on:')
        endif
    endif master_finalize
    !-------------------- End MPI/MASTER ---------------------------------
    
    !-------------------- MPI/ALL NODES ----------------------------------
    call MPI_FINALIZE(mpierror)
    !-------------------- End MPI/ALL NODES-------------------------------

    close(io)
    
end subroutine finalize


!-----------------------
! Fatal exit subroutine 
!-----------------------
subroutine quick_exit(io, status)
   
   !  quick_exit: exit procedure, designed to return an gentle way to exitm.
   use allmod
   implicit none
   integer io           ! close this unit if greater than zero
   integer status       ! exit status; 1-error 0-normal


   include 'mpif.h'
   integer ierr
   
   if (status /= 0) then
      call flush(io)
      call mpi_abort(MPI_COMM_WORLD, status, ierr)
   else
      call mpi_finalize(ierr)
   end if

   call finalize(io,1)
   
   call exit(status)

end subroutine quick_exit