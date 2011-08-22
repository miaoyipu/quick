!
!	matComp.f90
!	new_quick
!
!	Created by Yipu Miao on 4/18/11.
!	Copyright 2011 University of Florida. All rights reserved.
!

function rms(mat1,mat2,n)
    implicit none
    double precision rms
    integer n,i,j
    double precision mat1(n,n),mat2(n,n)
    
    rms=0.0d0
    do i=1,n
        do j=1,n
            rms=(mat1(i,j)-mat2(i,j))**2.0d0+rms
        enddo
    enddo
    rms = (rms/n**2)**0.5d0
end function rms

function maxDiff(mat1,mat2,n)
    implicit none
    double precision maxDiff
    integer n,i,j
    double precision mat1(n,n),mat2(n,n)
    
    maxDiff=0.d0
    do i=1,n
        do j=1,n
            maxDiff=max(maxDiff,dabs(mat1(i,j)-mat2(i,j)))
        enddo
    enddo
end function maxDiff