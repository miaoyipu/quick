program intdiff

integer i,j,k,l,m,n,i1,i2,i3,i4,i5,i6,i7,i8
integer j1,j2,j3,j4,j5,j6,j7,j8
double precision d1,d2
logical same
open (unit=2, file="TEST2")
open (unit=3, file="TEST")

!read(3,*) j1,j2,j3,j4,d1

do i = 1, 100000000
	same = .false.
	read(2,*) i1,i2,i3,i4,d2
	read(3,*) j1,j2,j3,j4,d1
!	if ((j1.eq.i1).and.(j2.eq.i2).and.(j3.eq.i3).and.(j4.eq.i4)) then
		if(abs(d1-d2).gt.1e-12) then
				write(*,*) "diff1",d1,d2,i1,i2,i3,i4,j1,j2,j3,j4
		else
				same = .true.
!				write(*,*) "same",i, d1,d2, i1,i2,i3,i4, j1,j2,j3,j4
		endif
		!read(3,*) j1,j2,j3,j4,d1
!	endif
!	if ((.not.same).and.abs(d2).gt.1e-10) then
!		write(*,*) "diff2",i1,i2,i3,i4,d2
!	endif
enddo


end program
