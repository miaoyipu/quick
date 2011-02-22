!********************************************************
! debug.f90
!********************************************************
! This file contains debug subroutines
! --Yipu Miao 10/09/2010
!********************************************************
! Subroutine List:
! debugSCF()                    : SCF debug
! debugElecdii()                : Elecdii debug
!

!*******************************************************
! debugSCF()
!-------------------------------------------------------
! this subroutine is to output some infos in debug mode
!
    subroutine debugSCF()
    use allmod
    implicit double precision(a-h,o-z)
    
    ioutfile=6
    !Densitry matrix
    write(ioutfile,'(//,"DENSITY MATRIX AT START OF", &
          &  " CYCLE",I4)') jscf
    do I=1,nbasis
        do J=1,nbasis
            write (ioutfile,'("DENSE[",I4,",",I4,"]=",F18.10)') &
          &          J,I,DENSE(J,I)
        enddo
    enddo
    
    ! Operator matrix
    write(ioutfile,'(//,"OPERATOR MATRIX FOR CYCLE",I4)') jscf
    do I=1,nbasis
        do J=1,nbasis
            write (ioutfile,'("o[",I4,",",I4,"]=",F18.10)') &
                    J,I,O(J,I)
        enddo
    enddo
    
    
    ! Transformed operator matrix
    write(ioutfile,'(//,"TRANSFORMED OPERATOR MATRIX", &
          &  " FOR CYCLE",I4)') jscf
    do I=1,nbasis
        do J=1,nbasis
            write (ioutfile,'("ot[",I4,",",I4,"]=",F18.10)') &
                J,I,O(J,I)
                debug1(J,I) = O(J,I)
        enddo
    enddo
    
    ! C' Matrix
    write(ioutfile,'(//,"C''  MATRIX", &
          &  " FOR CYCLE",I4)') jscf
    do I=1,nbasis
        do J=1,nbasis
            write (ioutfile,'("cprime[",I4,",",I4,"]=",F18.10)') &
                    J,I,VEC(J,I)
        enddo
    enddo
    
    write(ioutfile,'(//,"CHECK THE C'' MATRIX", &
          &  " FOR CYCLE",I4)') jscf
    write(ioutfile,'(/,"CHECK #1 IS FOR ORTHOGONALITY OF COLUMNS")')

    do I=1,nbasis
        do J=1,nbasis
            debug2IJ=0.d0
            do K=1,nbasis
                debug2IJ=debug2IJ+VEC(K,I)*VEC(K,J)
            enddo
            write (ioutfile,'(I4,I4,F18.10)') J,I,debug2IJ
        enddo
    enddo

    write(ioutfile,'(/,"CHECK #2 IS FOR ORBITAL ENERGIES")')
    do I=1,nbasis
        do J=1,nbasis
            debug2IJ=0.d0
            do K=1,nbasis
                debug2IJ=debug2IJ+VEC(K,I)*debug1(K,J)
            enddo
            debug2(I,J) = debug2IJ
        enddo
    enddo
    
    do I=1,nbasis
        do J=1,nbasis
            debug1IJ=0.d0
            do K=1,nbasis
                debug1IJ=debug1IJ+debug2(I,K)*VEC(K,J)
            enddo
            debug1(I,J) = debug1IJ
            write (ioutfile,'(I4,I4,E18.10)') J,I,debug1(I,J)
        enddo
    enddo
    
    ! Molecular orbital energy
    write(ioutfile,'(//,"MOLECULAR ORBITAL ENERGIES", &
          &  " FOR CYCLE",I4)') jscf
    do I=1,nbasis
        write (ioutfile,'(I4,F18.10)') I,E(I)
    enddo
    
    ! C matrix
    write(ioutfile,'(//,"C MATRIX", &
          &  " FOR CYCLE",I4)') jscf
    do I=1,nbasis
        do J=1,nbasis
            write (ioutfile,'("c[",I4,",",I4,"]=",F18.10)') &
                J,I,CO(J,I)
        enddo
    enddo
    
    
    ! Densitry matrix at end of this cycle  
    write(ioutfile,'(//,"DENSITY MATRIX AT end OF", &
          &  " CYCLE",I4)') jscf
    do I=1,nbasis
        do J=1,nbasis
            write (ioutfile,'("DENSE[",I4,",",I4,"]=",F18.10)') &
                J,I,DENSE(J,I)
        enddo
    enddo
    
    total=0.d0
    do I=1,nbasis
        do J=1,nbasis
            do K=1,ncontract(I)
                do L=1,ncontract(J)
                    total = total + DENSE(I,J)* &
                            dcoeff(K,I)*dcoeff(L,J) &
                            *overlap(aexp(L,J),aexp(K,I), &
                            itype(1,J),itype(2,J),itype(3,J), &
                            itype(1,I),itype(2,I),itype(3,I), &
                            xyz(1,ncenter(J)),xyz(2,ncenter(J)), &
                            xyz(3,ncenter(J)),xyz(1,ncenter(I)), &
                            xyz(2,ncenter(I)),xyz(3,ncenter(I)))
                enddo
            enddo
        enddo
    enddo
    
    write(ioutfile,'(//,"TOTAL NUMBER OF ELECTRONS=",F18.10)')total
            
    end subroutine debugSCF
    
!*******************************************************
! debugElecdii()
!-------------------------------------------------------
! this subroutine is to output some infos in debug mode
!
    subroutine debugElecdii()
    use allmod
    implicit double precision(a-h,o-z)
    
    write(iOutFile,*) "O OPERATOR for cycle ",jscf
    call PriSym(iOutFile,nbasis,O,'f14.8')
        
    write(iOutFile,*) "DENSITY MATRIX for cycle ", jscf
    call PriSym(iOutFile,nbasis,dense,'f14.8')
        
        
    end subroutine debugElecdii

!*******************************************************
! debugElecdii()
!-------------------------------------------------------
! this subroutine is to output normalization info for divcon
!   
    subroutine debugDivconNorm()
    use allmod
    
    implicit double precision(a-h,o-z)
    
    ! Xiaotemp inidicates normalization infos for system
    Xiaotemp=0.0d0
    do i=1,nbasis
        do j=1,nbasis
            Xiaotemp=Xiaotemp+DENSE(j,i)*Smatrix(j,i)
        enddo
    enddo

    ! Write Normalization for both full system and subsystem
    write(ioutfile,*)'ELECTION NORMALIZATION'
    write(ioutfile,*)'-------------------------------'
    write(ioutfile,*)'FULL=',Xiaotemp
    write(ioutfile,*)'SUBSYSTEM     NORMALIZATION'
    do itt=1,np
        Xiaotemp=0.0d0
        do i=1,nbasisdc(itt)
            do j=1,nbasisdc(itt)
                Xiaotemp=Xiaotemp+Pdcsub(itt,j,i)*smatrixdcsub(itt,j,i)
            enddo
        enddo
        write(ioutfile,*) itt,Xiaotemp
    enddo
    write(ioutfile,*)'-------------------------------'
    
    end subroutine debugDivconNorm
        