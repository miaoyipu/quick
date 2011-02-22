!*******************************************************
! getAtoms(iAtomType)
!-------------------------------------------------------
! Ed Brothers. November 13, 2001
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP
!-------------------------------------------------------
! The purpose of this subroutine is to calculate the transformation
! matrix X.  The first step is forming the overlap matrix (Smatrix).
!
    subroutine fullx
    use allmod
    implicit double precision(a-h,o-z)

    dimension Sminhalf(nbasis),V(3,nbasis),IDEGEN1(nbasis)

    do Ibas=1,nbasis
        do Jbas=Ibas,nbasis
            SJI =0.d0
            do Icon=1,ncontract(ibas)
                do Jcon=1,ncontract(jbas)
                    SJI =SJI + &
                    dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                    *overlap(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                    itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                    itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                    xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                    xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                    xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                enddo
            enddo
            Smatrix(Jbas,Ibas) = SJI
        enddo
    enddo
    
    if (quick_method%debug) then
        write(ioutfile,'(/"THE OVERLAP MATRIX")')
        call PriSym(iOutFile,nbasis,Smatrix,'F18.10')
                call flush(iOutFile)
    endif

    do Ibas=1,nbasis
        do Jbas=Ibas,nbasis
            HOLD(Jbas,Ibas) = Smatrix(Jbas,Ibas)
            Smatrix(Ibas,Jbas) = Smatrix(Jbas,Ibas)
        enddo
    enddo

! Now diagonalize HOLD to generate the eigenvectors and eigenvalues.
  
    call DIAG(NBASIS,HOLD,NBASIS,TOL,V,Sminhalf,IDEGEN1,Uxiao,IERROR)

! Consider the following:

! X = U * s^(-.5) * transpose(U)

! s^-.5 is a diagonal matrix filled with the eigenvalues of S taken to
! to the 1/square root.  If we define an intermediate matrix A for the
! purposes of this discussion:

! A   = U * s^(-.5)
! or Aij = Sum(k=1,m) Uik * s^(-.5)kj

! s^(-.5)kj = 0 unless k=j so

! Aij = Uij * s^(-.5)jj

! X   = A * transpose(U)
! Xij = Sum(k=1,m) Aik * transpose(U)kj
! Xij = Sum(k=1,m) Uik * s^(-.5)kk * transpose(U)kj
! Xij = Sum(k=1,m) Uik * s^(-.5)kk * Ujk

! Similarly:
! Xji = Sum(k=1,m) Ajk * transpose(U)ki
! Xji = Sum(k=1,m) Ujk * s^(-.5)kk * transpose(U)ki
! Xji = Sum(k=1,m) Ujk * s^(-.5)kk * Uik

! This aggravating little demonstration contains two points:
! 1)  X can be calculated without crossing columns in the array
! which adds to speed.
! 2)  X has to be symmetric. Thus we only have to fill the bottom
! half. (Lower Diagonal)

    do I=1,nbasis
        Sminhalf(I) = Sminhalf(I)**(-.5d0)
    enddo

! Transpose U onto X then copy on to U.  Now U contains U transpose.

    do I = 1,nbasis
        do J = 1,nbasis
            X(I,J) = Uxiao(J,I)
        enddo
    enddo
    do I = 1,nbasis
        do J = 1,nbasis
            Uxiao(J,I) = X(J,I)
        enddo
    enddo

! Now calculate X.
! Xij = Sum(k=1,m) Transpose(U)kj * s^(-.5)kk * Transpose(U)ki

    do I = 1,nbasis
        do J=I,nbasis
            sum = 0.d0
            do K=1,nbasis
                sum = Uxiao(K,I)*Uxiao(K,J)*Sminhalf(K)+sum
            enddo
            X(I,J) = sum
            X(J,I) = X(I,J)
        enddo
    enddo

    if (quick_method%debug) then
        write(ioutfile,'(/"THE X MATRIX")')
        call PriSym(iOutFile,nbasis,X,'F18.10')
    endif

! At this point we have the transformation matrix (X) which is necessary
! to orthogonalize the operator matrix, and the overlap matrix (S) which
! is used in the DIIS-SCF procedure.

    return
    end subroutine fullx
