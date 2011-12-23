program writecode
do i = 1,3
	do j= 10,19
		write(*,'(a,i2,a,i2,a,I0,a,i0,a)') "LOC2(store,",i,",",j,", STOREDIM, STOREDIM) += x_",i,"_",j,"_0;"
	enddo
enddo
end
