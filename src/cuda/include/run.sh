for i in {0..6}
do
	for j in {0..6}
	do
	#a='h_'$i'_'$j'.h'
	#b='h_'$i'_'$j'.cu'
        c='f_'$i'_'$j'.h'
		c2='f2_'$i'_'$j'.h'
	#	cp $c $c2
	sed -i '' 's/h_/h2_/1' $c2 
#echo $b "\\"
#cp $b $c
#echo "#include " "\"./"$c"\"" 
	#rm $b
	#echo "#include" "\""../gpu_common.h"\"" >> $b
        #echo "#include" "\""./h_all_files.h"\"" >> $b
	#cat $a >> $b
	#	#echo "__device__ __noinline__ void h_$a( QUICKDouble* YVerticalTemp, QUICKDouble* store,QUICKDouble Ptempx, QUICKDouble Ptempy, QUICKDouble Ptempz, QUICKDouble WPtempx,QUICKDouble WPtempy,QUICKDouble WPtempz, QUICKDouble Qtempx, QUICKDouble Qtempy, QUICKDouble Qtempz,  QUICKDouble WQtempx,QUICKDouble WQtempy,QUICKDouble WQtempz, QUICKDouble ABCDtemp,QUICKDouble ABtemp, QUICKDouble CDtemp, QUICKDouble ABcom, QUICKDouble CDcom );" >> h_all_files.h
	done
done
