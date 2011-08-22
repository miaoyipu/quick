!
!	xnorm.f90
!	new_quick
!
!	Created by Yipu Miao on 2/23/11.
!	Copyright 2011 University of Florida. All rights reserved.
!

!-----------------------------------------------------------
! xnorm(a,i,j,k)
!-----------------------------------------------------------
double precision function xnorm(a,i,j,k)
  use quick_constants_module
  implicit double precision(a-h,o-z)
  !-----------------------------------------------------------
  ! Ed Brothers. October 3, 2001
  ! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP
  ! The purpose of this function is to calculate the normalization
  ! constant of a gtf with exponent a and x,y,z exponents of i,j,k
  ! According to Mathematica, the simplied normalization constant is:
  ! 3/4 + (i + j + k)/2  3/4 + (i + j + k)/2
  ! 2                    a
  ! ut[3]= --------------------------------------------------------
  ! 1                  1                  1
  ! Sqrt[Gamma[- + i]] Sqrt[Gamma[- + j]] Sqrt[Gamma[- + k]]
  ! 2                  2                  2

  ! Constants:

  !    pi=3.1415926535897932385
  !    pito3half=5.568327996831707845284817982118835702014
  !-----------------------------------------------------------

  xnorm =   (2.d0*a)**(dble(i+j+k)/2.d0+.75d0)

  ! Before the Gamma function code, a quick note. All gamma functions are
  ! of the form:
  ! 1
  ! Gamma[- + integer].  Now since Gamma[z] = (z-1)Gamma(z-1)
  ! 2

  ! We can say Gamma(0.5 + i) = (i-1+.5)(i-2+.5)...(i-i+.5)Gamma(0.5)
  ! and Gamma(.5) is Sqrt(Pi).  Thus to calculate the three gaussians
  ! just requires a loop and multiplying by Pi^3/2

  gamma=1.d0
  do L=1,i
     gamma = gamma * (dble(i - L) + .5d0)
  enddo
  do L=1,j
     gamma = gamma * (dble(j - L) + .5d0)
  enddo
  do L=1,k
     gamma = gamma * (dble(k - L) + .5d0)
  enddo

  xnorm = xnorm*PI**(-0.75d0)*gamma**(-.5d0)

  return
END function xnorm



!-----------------------------------------------------------
! xnewnorm(a,i,j,k)
!--------------------------------------------------------
! Ed Brothers. October 3, 2001
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP
!-----------------------------------------------------------
double precision function xnewnorm(i,j,k,nprim1,ddcoe,aaexp)
  use allmod
  implicit double precision(a-h,o-z)
  real*8 ddcoe(nprim1),aaexp(nprim)

  pi32 = pi**(1.5)
  xnewnorm= 0.d0
  nxponent=-i-j-k
  xponent=-1.5d0+dble(nxponent)
  do Icon1 = 1,nprim1
     do Icon2 = 1,nprim1
        xnewnorm = xnewnorm + ddcoe(Icon1)*ddcoe(Icon2) &
             *(aaexp(Icon1)+aaexp(Icon2))**xponent
     enddo
  enddo

  gamma=1.d0
  do L=1,i
     gamma = gamma * (dble(i - L) + .5d0)
  enddo
  do L=1,j
     gamma = gamma * (dble(j - L) + .5d0)
  enddo
  do L=1,k
     gamma = gamma * (dble(k - L) + .5d0)
  enddo
  xnewnorm = (xnewnorm*gamma*pi32)**(-.5d0)

  return
END function xnewnorm
