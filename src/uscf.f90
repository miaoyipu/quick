! Ed Brothers. December 20, 2001
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

    subroutine uscf(failed)
    use allmod
    implicit double precision(a-h,o-z)

    logical :: done,failed
    done = .false.

! The purpose of this subroutine is to perform open shell
! scf cycles.  This is very close to the subroutine 'SCF' At this
! point, X has been formed. The remaining steps are:
! 1)  Form alpha operator matrix.
! 2)  Calculate O' = Transpose[X] O X
! 3)  Diagonalize O' to obtain C' and eigenvalues.
! 4)  Calculate C = XC'
! 5)  Form new alpha density matrix.
! 6)  Repeat steps 1-5 for beta density.
! 7)  Check for convergence.

! Each location in the code that the step is occurring will be marked.
! The cycles stop when prms is less than pmaxrms or when the maximum
! number of scfcycles has been reached.

! NOTE:  Dense and nelec are used in the USCF case to denote the alpha
! density and number of alpha electrons, respectively.

    sval = 1.d0
    jscf=0
    PRMS =1.D30

!
! Alessandro GENONI 03/21/2007
! ECP integrals computation exploiting Alexander V. Mitin Subroutine
! Note: the integrals are stored in the array ecp_int that corresponds
!       to the lower triangular matrix of the ECP integrals
!
    if (quick_method%ecp) then
      call ecpint
!
! ECP integrals DEBUG
!
      write(ioutfile,*) '  '
      write(ioutfile,*) 'ECP INTEGRALS DEBUG'
      do i=1,nbf12
        write(ioutfile,*) i,'.  ',ecp_int(i)
      end do
    end if

    done=jscf.ge.iscf

    do WHILE (.not.done)
        if (quick_method%diisscf .AND. PRMS < 1.D-1) call uelectdiis(jscf,PRMS)
        call cpu_time(t1)
        jscf=jscf+1

        if (quick_method%debug) then
            write(ioutfile,'(//,"ALPHA DENSITY MATRIX AT START OF", &
            & " CYCLE",I4)') jscf
            do I=1,nbasis
                do J=1,nbasis
                    write (ioutfile,'("DENSEA[",I4,",",I4,"]=",F18.10)') &
                    J,I,DENSE(J,I)
                enddo
            enddo
            write(ioutfile,'(//,"BETA DENSITY MATRIX AT START OF", &
            & " CYCLE",I4)') jscf
            do I=1,nbasis
                do J=1,nbasis
                    write (ioutfile,'("DENSEB[",I4,",",I4,"]=",F18.10)') &
                    J,I,DENSEB(J,I)
                enddo
            enddo
        endif


    ! 1)  Form alpha operator matrix.

        if (quick_method%HF) call uhfoperatora
        if (quick_method%DFT) call udftoperatora
        if (quick_method%SEDFT) call usedftoperatora

       ! call shift(sval,.true.,jscf)

        if (quick_method%debug) then
            write(ioutfile,'(//,"ALPHA OPERATOR MATRIX AT START OF", &
            & " CYCLE",I4)') jscf
            do I=1,nbasis
                do J=1,nbasis
                    write (ioutfile,'("OPERATORA[",I4,",",I4,"]=",F18.10)') &
                    J,I,O(J,I)
                enddo
            enddo
        endif

    ! 2)  Calculate O' = Transpose[X] O X.
    
    ! X and O are both symmetric, so O' = XOX. As X is symmetric,
    ! calculate OX first, store in HOLD, and calculate X (OX).
    
    ! The matrix multiplier comes from Steve Dixon. It calculates
    ! C = Transpose(A) B.  For symmetric matrices this is fine, but
    ! for some calculations you first have to transpose the A matrix.


        do I=1,nbasis
            do J=1,nbasis
                HOLDIJ = 0.0D0
                do K=1,nbasis
                    HOLDIJ = HOLDIJ + O(K,I)*X(K,J)
                enddo
                HOLD(I,J) = HOLDIJ
            enddo
        enddo

        do I=1,nbasis
            do J=1,nbasis
                OIJ = 0.0D0
                do K=1,nbasis
                    OIJ = OIJ + X(K,J)*HOLD(K,I)
                enddo
                O(I,J) = OIJ
            enddo
        enddo
    
    ! 3)  Diagonalize O' to obtain C' and alpha eigenvalues.
    ! In this case VEC is now C'.
    ! There is a problem if HOLD is used for evec1.
       
       
        CALL DIAG(nbasis,O,nbasis,TOL,V2,E,IDEGEN,VEC, &
        IERROR)
        

    ! 4)  Calculate C = XC'
    ! Note here too that we are actuall calculating C = Transpose[X] C',
    ! which is the same as C=XC' as X is symmetric.

        do I=1,nbasis
            do J=1,nbasis
                CIJ = 0.0D0
                do K=1,nbasis
                    CIJ = CIJ + X(K,I)*VEC(K,J)
                enddo
                CO(I,J) = CIJ
            enddo
        enddo
        

    ! 5)  Form new alpha density matrix.
    ! Save the old density matrix in HOLD for convergence comparison.
    ! Also, calculate the spin contamination.

        do I=1,nbasis
            do J=1,nbasis
                HOLD(J,I) = DENSE(J,I)
            enddo
        enddo

        do I=1,nbasis
            do J=1,nbasis
                DENSE(J,I) = 0.d0
                do K=1,nelec
                    DENSE(J,I) = DENSE(J,I) + (CO(J,K)*CO(I,K))
                enddo
            enddo
        enddo

    ! 6)  Check for change in the alpha matrix.

        OLDPRMS=PRMS
        PRMS=0.d0
        PCHANGE=0.d0
        do I=1,nbasis
            do J=1,nbasis
                PRMS=PRMS+(DENSE(J,I)-HOLD(J,I))**2.d0
                PCHANGE=MAX(PCHANGE,ABS(DENSE(J,I)-HOLD(J,I)))
            enddo
        enddo
        PRMS = (PRMS/nbasis**2.d0)**0.5d0

    ! 1)  Form beta operator matrix.

        if (quick_method%HF) call uhfoperatorb
        if (quick_method%DFT) call udftoperatorb
        if (quick_method%SEDFT) call usedftoperatorb

       ! call shift(sval,.false.,jscf)

    ! 2)  Calculate O' = Transpose[X] O X.
    
    ! X and O are both symmetric, so O' = XOX. As X is symmetric,
    ! calculate OX first, store in HOLD, and calculate X (OX).
    
    ! The matrix multiplier comes from Steve Dixon. It calculates
    ! C = Transpose(A) B.  For symmetric matrices this is fine, but
    ! for some calculations you first have to transpose the A matrix.


        do I=1,nbasis
            do J=1,nbasis
                HOLDIJ = 0.0D0
                do K=1,nbasis
                    HOLDIJ = HOLDIJ + O(K,I)*X(K,J)
                enddo
                HOLD2(I,J) = HOLDIJ
            enddo
        enddo

        do I=1,nbasis
            do J=1,nbasis
                OIJ = 0.0D0
                do K=1,nbasis
                    OIJ = OIJ + X(K,J)*HOLD2(K,I)
                enddo
                O(I,J) = OIJ
            enddo
        enddo

    
    ! 3)  Diagonalize O' to obtain C' and BETA eigenvalues.
    ! In this case VEC is now C'.
    ! There is a problem if HOLD is used for evec1.

        CALL DIAG(nbasis,O,nbasis,TOL,V2,EB,IDEGEN,VEC, &
        IERROR)

    
    ! 4)  Calculate C = XC'
    ! Note here too that we are actuall calculating C = Transpose[X] C',
    ! which is the same as C=XC' as X is symmetric.

        do I=1,nbasis
            do J=1,nbasis
                CIJ = 0.0D0
                do K=1,nbasis
                    CIJ = CIJ + X(K,I)*VEC(K,J)
                enddo
                COB(I,J) = CIJ
            enddo
        enddo

    ! 5)  Form new beta density matrix.
    ! Save the old density matrix in HOLD2 for convergence comparison.
    ! Also, check spin contamination.

        do I=1,nbasis
            do J=1,nbasis
                HOLD2(J,I) = DENSEB(J,I)
            enddo
        enddo

        do I=1,nbasis
            do J=1,nbasis
                DENSEB(J,I) = 0.d0
                do K=1,nelecb
                    DENSEB(J,I) = DENSEB(J,I) + (COB(J,K)*COB(I,K))
                enddo
            enddo
        enddo

    ! 6)  Check for change in the beta matrix.

        prms2=0.d0
        do I=1,nbasis
            do J=1,nbasis
                PRMS2=PRMS2+(DENSEB(J,I)-HOLD2(J,I))**2.d0
                PCHANGE=MAX(PCHANGE,ABS(DENSEB(J,I)-HOLD2(J,I)))
            enddo
        enddo
        PRMS2 = (PRMS2/nbasis**2.d0)**0.5d0

        prms = max(prms,prms2)

        call cpu_time(t2)

        write (ioutfile,'(/,"SCF CYCLE      = ",I8, &
        & "      TIME      = ",F8.2)') &
        jscf,T2-T1
        write (ioutfile,'("RMS CHANGE     = ",E12.6, &
        & "  MAX CHANGE= ",E12.6)') &
        PRMS,PCHANGE
        if (quick_method%DFT .OR. quick_method%SEDFT) then
            write (ioutfile,'("ALPHA ELECTRON DENSITY    =",F16.10)') &
            aelec
            write (ioutfile,'("BETA ELECTRON DENSITY     =",F16.10)') &
            belec
        endif

        if (quick_method%prtgap) write (ioutfile,'("ALPHA HOMO-LUMO GAP (EV) =", &
        & 11x,F12.6)') (E(nelec+1) - E(nelec))*27.2116d0
        if (quick_method%prtgap) write (ioutfile,'("BETA HOMO-LUMO GAP (EV)  =", &
        & 11x,F12.6)') (EB(nelecb+1) - EB(nelecb))*27.2116d0
        if (PRMS < PMAXRMS .AND. pchange < PMAXRMS*100.d0)then
            write (ioutfile,'("CONVERGED!!!!!")')
            done=.true.
        elseif(jscf >= iscf) then
            write (ioutfile,'("RAN OUT OF CYCLES.  NO CONVERGENCE.")')
            done=.true.
            failed=.true.
        elseif (prms >= oldprms) then
            if (mod(dble(jscf),2.d0) == 0.d0) then
                write (ioutfile,'("NOT IMPROVING.  ", &
                & "TRY MODifYING ALPHA DENSITY MATRIX.")')
                do Ibas=1,nbasis
                    do Jbas=1,nbasis
                        DENSE(Jbas,Ibas) =.7d0*HOLD(Jbas,Ibas) &
                        +.3d0*DENSE(Jbas,Ibas)
                    enddo
                enddo
            else
                write (ioutfile,'("NOT IMPROVING.  ", &
                & "TRY MODifYING BETA DENSITY MATRIX.")')
                do Ibas=1,nbasis
                    do Jbas=1,nbasis
                        DENSEB(Jbas,Ibas) = .7d0*HOLD2(Jbas,Ibas) &
                        +.3d0*DENSEB(Jbas,Ibas)
                    enddo
                enddo
            endif
            prms = 1.D30
        else
            continue
        endif
    enddo

    END subroutine uscf