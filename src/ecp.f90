      subroutine ecpint
!     Contact A. V. Mitin for the ECP integral package
      return
      end
      
      Subroutine ecpoperator
!
! Alessandro GENONI 03/12/2007
! Suibroutine that adds the ECP integrals to the Fock operator
! for ECP calculations
! 
      use quick_basis_module
      use quick_ecp_module
      use quick_calculated_module
      use quick_files_module
!
      implicit double precision (a-h,o-z)
!
! Add the proper ECP integral to the corresponding Fock matrix element
!
      do i=1,nbasis
       do j=1,nbasis
        ind=kp2(i,j)
!     write(ioutfile,*) i,j,'---->',ind
        O(i,j)=O(i,j)+ecp_int(ind)
       end do
      end do
      return
      end



 integer function kp2(ii,jj) result(ind)
  use quick_ecp_module
  implicit none
  integer, intent(in)  :: ii,jj
!
  ind=kvett(max(ii,jj))+(ii+jj-max(ii,jj))
  return
 end function kp2
