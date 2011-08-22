program intdiff

integer i,j,k,l,m,n,i1,i2,i3,i4,i5,i6,i7,i8
integer j1,j2,j3,j4,j5,j6,j7,j8
double precision d1,d2
logical same
open (unit=2, file="TEST2")
open (unit=3, file="TEST")

read(3,*) j,k,d2,j1,j2,j3,j4,j5,j6,j7,j8

do i = 1, 100000000
	same = .false.
	read(2,*) j,k,d1,i1,i2,i3,i4,i5,i6,i7,i8
	if ((j1.eq.i1).and.(j2.eq.i2).and.(j3.eq.i3).and.(j4.eq.i4) &
		.and.(j5.eq.i5).and.(j6.eq.i6).and.(j7.eq.i7).and.(j8.eq.i8)) then
		if(abs(d1-d1).gt.1e-5) then
				write(*,*) "diff1",j,k,d1,d2,i1,i2,i3,i4,i5,i6,i7,i8
		else
				same = .true.
!				write(*,*) "same",d1,d2, i1,i2,i3,i4,i5,i6,i7,i8, j1,j2,j3,j4,j5,j6,j7,j8
		endif
		read(3,*) j,k,d2,j1,j2,j3,j4,j5,j6,j7,j8
	endif
	if ((.not.same).and.abs(d1).gt.1e-10) then
		write(*,*) "diff2",j,k,d1,i1,i2,i3,i4,i5,i6,i7,i8
	endif
enddo


end program
