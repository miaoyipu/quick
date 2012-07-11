#include "config.h"
! ED BROTHERS. NOVEMBER 27, 2001
! XIAO HE. SEPTEMBER 14,2008
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

SUBROUTINE CALMP2
   USE ALLMOD
   USE QUICK_GAUSSIAN_CLASS_MODULE
   IMPLICIT NONE
#ifndef CUDA
   DOUBLE PRECISION CUTOFFTEST,TESTTMP,TESTCUTOFF
   INTEGER II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2
   COMMON /HRRSTORE/II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2
   INTEGER :: NELEC,NELECB

   INTEGER :: IOCC, IVIR               ! NUMBER OF OCCUPIED AND VIRTUAL ORBITALS
   DOUBLE PRECISION :: CUTOFFMP2           ! MP2 CUTOFF
   DOUBLE PRECISION :: EMEMORYSUM         ! MEMOERY USAGE

   INTEGER :: STEP_LENGTH, NSTEP

   INTEGER :: KMAXVAL                      ! KMAXVAL: for d = 6, for f = 10


   INTEGER :: TEMP, I, J, K, L, TOTALSPACE
   INTEGER :: ICYCLE, II111, II112, IIINEW, J3, J33, J33NEW, JJ111, JJ112, JJJNEW
   INTEGER :: I3, K3, K33, KK111, KK112, L3, L3NEW, LL111, LL112, NII1, NII2, NJJ1, NJJ2, NKK1, NKK2, NLL1, NLL2
   INTEGER :: CURRENTSTEPLENGTH, i3mp2, i3mp2new
   INTEGER :: ISTEP, NSTEPS, NSTEPF, NINT, IINEW,JJNEW,KKKnew
   INTEGER :: III1, III2, JJJ1,JJJ2,KKK1,KKK2,LLL1,LLL2

   DOUBLE PRECISION :: MAX_CUTOFF, COMAX, ATEMP, ATEMP2, BTEMP
   DOUBLE PRECISION, ALLOCATABLE :: TEMP2D1(:,:), TEMP2D2(:,:), orbmp2j33(:,:,:)
   LOGICAL :: PASSCH
   DOUBLE PRECISION :: SINT, TINTS, TINTF
   DOUBLE PRECISION :: S1QT, T1QTS, T1QTF
   DOUBLE PRECISION :: T3, T4, T5, T6, T7, T8, T9
   DOUBLE PRECISION ::  S2, S3, S4, S5

   SINT=0.0
   S1QT=0.0
   S3=0.0
   S4=0.0
   S5=0.0
   NELEC = QUICK_MOLSPEC%NELEC
   NELECB = QUICK_MOLSPEC%NELECB

   CALL PRTACT(IOUTFILE,"BEGIN MP2 CALCULATION")

   CUTOFFMP2= MIN(1.0D-8, QUICK_METHOD%INTEGRALCUTOFF)          ! INTEGRAL CUTOFF CRITERIA
   QUICK_METHOD%PRIMLIMIT= MIN(1.0D-8, QUICK_METHOD%PRIMLIMIT)  ! PRIMTIVE CUTOFF
   QUICK_QM_STRUCT%EMP2=0.0D0

   ! OCCUPIED AND VIRTUAL ORBITALS NUMBER
   IOCC=NELEC/2
   IVIR=NBASIS-NELEC/2
   WRITE(IOUTFILE, '(" NUMBER OF OCCUPIED ORBITALS = ", I5)') IOCC
   WRITE(IOUTFILE, '(" NUMBER OF VIRTUAL  ORBITALS = ", I5)') IVIR
   WRITE(IOUTFILE, '(" NUMBER OF ELECTRONS         = ", I5)') NELEC
   WRITE(IOUTFILE, '(" INTEGRAL CUTOFF             = ", E12.6)') QUICK_METHOD%INTEGRALCUTOFF
   WRITE(IOUTFILE, '(" PRIM CUTOFF                 = ", E12.6)') QUICK_METHOD%PRIMLIMIT


   ! ACTUALLY NSTEP IS STEP LENGTH
   STEP_LENGTH=IOCC !MIN(INT(15D0/EMEMORYSUM),NELEC/2)

   ! IF WITH F ORBITAL
   IF(QUICK_METHOD%FFUNXIAO)THEN
      KMAXVAL=6
   ELSE
      KMAXVAL=10
   ENDIF

   ! ALLOCATE SOME VARIABLES
   !   ALLOCATE(MP2SHELL(NBASIS))
   ALLOCATE(ORBMP2(IVIR,IVIR))                                  ! O(V^2)
   ALLOCATE(ORBMP2I331(STEP_LENGTH,NBASIS,KMAXVAL,KMAXVAL,2))   ! O(IN)
   ALLOCATE(ORBMP2J331(STEP_LENGTH,IVIR,KMAXVAL,KMAXVAL,2))     ! O(IV)
   !ALLOCATE(orbmp2j33(STEP_LENGTH,IVIR,2))     ! O(ON)
   ALLOCATE(ORBMP2K331(STEP_LENGTH,IOCC,IVIR,NBASIS))           ! O(IOVN)
   ALLOCATE(AO(NBASIS, KMAXVAL, KMAXVAL, KMAXVAL))              ! O(N)
   ALLOCATE(QUICK_QM_STRUCT%COOCC(IOCC,NBASIS))
   ALLOCATE(QUICK_QM_STRUCT%COVIR(IVIR,NBASIS))

   DO I = 1, NBASIS
      DO J = 1, IOCC
         QUICK_QM_STRUCT%COOCC( J, I) = QUICK_QM_STRUCT%CO( I, J)
      ENDDO
      DO J = 1, IVIR
         QUICK_QM_STRUCT%COVIR( J, I) = QUICK_QM_STRUCT%CO( I, J+IOCC)
      ENDDO
   ENDDO

   TOTALSPACE = 0
   TEMP = NBASIS*KMAXVAL*KMAXVAL*KMAXVAL
   EMEMORYSUM=REAL(TEMP*8.0D0/1024.0D0/1024.0D0) !/1024.0D0)
   WRITE(IOUTFILE,'(" AO STORAGE SPACE            = ",F8.2," M")') EMEMORYSUM
   TOTALSPACE = TOTALSPACE + TEMP

   TEMP = STEP_LENGTH*NBASIS*KMAXVAL*KMAXVAL*2
   EMEMORYSUM=REAL(TEMP*8.0D0/1024.0D0/1024.0D0) !/1024.0D0)
   WRITE(IOUTFILE,'(" 1ST TRANSFORMATION SPACE    = ",F8.2," M")') EMEMORYSUM
   TOTALSPACE = TOTALSPACE + TEMP

   TEMP = STEP_LENGTH*IVIR*KMAXVAL*KMAXVAL*2
   EMEMORYSUM=REAL(TEMP*8.0D0/1024.0D0/1024.0D0) !/1024.0D0)
   WRITE(IOUTFILE,'(" 2ND TRANSFORMATION SPACE    = ",F8.2," M")') EMEMORYSUM
   TOTALSPACE = TOTALSPACE + TEMP

   TEMP = STEP_LENGTH*IOCC*IVIR*NBASIS
   EMEMORYSUM=REAL(TEMP*8.0D0/1024.0D0/1024.0D0) !/1024.0D0)
   WRITE(IOUTFILE,'(" 3RD TRANSFORMATION SPACE    = ",F8.2," M")') EMEMORYSUM
   TOTALSPACE = TOTALSPACE + TEMP

   TEMP = IVIR * IVIR
   EMEMORYSUM=REAL(TEMP*8.0D0/1024.0D0/1024.0D0) !/1024.0D0)
   WRITE(IOUTFILE,'(" 4TH TRANSFORMATION SPACE    = ",F8.2," M")') EMEMORYSUM
   TOTALSPACE = TOTALSPACE + TEMP

   EMEMORYSUM=REAL(TOTALSPACE*8.0D0/1024.0D0/1024.0D0) !/1024.0D0)
   WRITE(IOUTFILE,'(" TOTAL MEMORY USAGE          = ",F8.2," M")') EMEMORYSUM
   WRITE(IOUTFILE,'(" STEP LENGTH                 = [",I3,"] * 2 ELECTRONS")') STEP_LENGTH

   ! WITH STEP_LENGTH(ACUTALLY, IT REPRESETNS STEP LENGHT), WE CAN
   ! HAVE NO. OF STEPS FOR MP2 CALCULATION
   NSTEP=NELEC/2/STEP_LENGTH
   NSTEP=NSTEP+1
   IF(STEP_LENGTH*(NSTEP-1).EQ.NELEC/2)THEN
      NSTEP=NSTEP-1
   ENDIF

   WRITE(*,'(" TOTAL STEPS                 =",I6)') NSTEP

   ! PRE-STEP FOR DENSITY CUTOFF
   CALL DENSITYCUTOFF

   ! SAVE COEFFECIENT FIRST.
   DO I=1,NBASIS
      DO J=1,NBASIS
         QUICK_SCRATCH%HOLD(I,J)=QUICK_QM_STRUCT%CO(J,I)
      ENDDO
   ENDDO

   MAX_CUTOFF=MAXVAL(YCUTOFF) ! MAX VALUE OF YCUTOFF

   DO ISTEP=1,NSTEP               ! STEP COUNTER

      CALL CPU_TIME(TIMER_BEGIN%TMP2)

      ! STEP LENGTH, FROM NSTEPS TO NSTEPF
      NINT=0                         ! INTEGER COUNTER
      NSTEPS=(ISTEP-1)*STEP_LENGTH+1    ! STEP START N
      NSTEPF=ISTEP*STEP_LENGTH          ! STEP END N
      IF(ISTEP.EQ.NSTEP) NSTEPF=NELEC/2
      CURRENTSTEPLENGTH=NSTEPF-NSTEPS+1

      ALLOCATE(TEMP2D1(NBASIS,CURRENTSTEPLENGTH))
      ALLOCATE(TEMP2D2(IVIR, CURRENTSTEPLENGTH))

      WRITE(IOUTFILE,'(" STEP ",I5,": FROM ORBITAL ",I5," TO ",I5)') ISTEP, NSTEPS, NSTEPF

      ! INITIAL ORBMP2K331
      CALL INITIALORBMP2K331(ORBMP2K331,STEP_LENGTH,NBASIS,IVIR,IOCC,CURRENTSTEPLENGTH)

      ! --- BEGIN I LOOP ---
      DO II=1,JSHELL

         II111=QUICK_BASIS%KSUMTYPE(II)+QUICK_BASIS%QSBASIS(II,QUICK_BASIS%QSTART(II))
         II112=QUICK_BASIS%KSUMTYPE(II)+QUICK_BASIS%QFBASIS(II,QUICK_BASIS%QFINAL(II))


         IF ( MOD( II, 25) .EQ. 0)  WRITE(IOUTFILE,'(" EVALUATING SHELL ", I5, "")') II

         ! --- BEGIN J LOOP ---
         DO JJ=II,JSHELL
            IF(ABS(YCUTOFF(II,JJ) * MAX_CUTOFF) .GT.CUTOFFMP2)THEN

               JJ111=QUICK_BASIS%KSUMTYPE(JJ)+QUICK_BASIS%QSBASIS(JJ,QUICK_BASIS%QSTART(JJ))
               JJ112=QUICK_BASIS%KSUMTYPE(JJ)+QUICK_BASIS%QFBASIS(JJ,QUICK_BASIS%QFINAL(JJ))

               CALL INITIALORBMP2IJ(ORBMP2I331,STEP_LENGTH,CURRENTSTEPLENGTH,NBASIS,KMAXVAL,KMAXVAL)
               !CALL INITIALORBMP2J(ORBMP2J33, CURRENTSTEPLENGTH, IVIR)
               CALL INITIALORBMP2IJ(ORBMP2J331,STEP_LENGTH,CURRENTSTEPLENGTH,IVIR,KMAXVAL,KMAXVAL)
               PASSCH = .FALSE.

               ! --- BEGIN K LOOP ---
               DO KK=1,JSHELL
                  KK111=QUICK_BASIS%KSUMTYPE(KK)+QUICK_BASIS%QSBASIS(KK,QUICK_BASIS%QSTART(KK))
                  KK112=QUICK_BASIS%KSUMTYPE(KK)+QUICK_BASIS%QFBASIS(KK,QUICK_BASIS%QFINAL(KK))

                  ! --- BEGIN L LOOP ---
                  DO LL=KK,JSHELL

                     LL111=QUICK_BASIS%KSUMTYPE(LL)+QUICK_BASIS%QSBASIS(LL,QUICK_BASIS%QSTART(LL))
                     LL112=QUICK_BASIS%KSUMTYPE(LL)+QUICK_BASIS%QFBASIS(LL,QUICK_BASIS%QFINAL(LL))

                     ! SCHWARTS CUTOFF IS IMPLEMENTED HERE
                     COMAX=0.D0
                     TESTCUTOFF = YCUTOFF(II,JJ)*YCUTOFF(KK,LL)
                     DO III = II111, II112
                        DO JJJ = MAX(III,JJ111), JJ112
                           DO KKK = KK111, KK112
                              DO LLL = MAX(KKK, LL111), LL112
                                 AO(LLL, III-II111+1, JJJ-JJ111+1, KKK-KK111+1) = 0
                              ENDDO
                           ENDDO
                        ENDDO
                     ENDDO

                     IF(DABS(TESTCUTOFF).GT.CUTOFFMP2)THEN

                        ! --- BEGIN LOOP FOR OVER I BATCH
                        DO ICYCLE=1,CURRENTSTEPLENGTH
                           TEMP=NSTEPS+ICYCLE-1
                           DO KKK=KK111,KK112
                              DO LLL=MAX(KKK,LL111),LL112
                                 COMAX=MAX(COMAX,DABS(QUICK_QM_STRUCT%CO(KKK,TEMP)))
                                 COMAX=MAX(COMAX,DABS(QUICK_QM_STRUCT%CO(LLL,TEMP)))
                              ENDDO
                           ENDDO
                        ENDDO
                        TESTCUTOFF=DABS(TESTCUTOFF*COMAX)

                        ! -- SHELLMP2 INCLUDES ERI EVALUATION
                        IF(TESTCUTOFF.GT.CUTOFFMP2)THEN
                           PASSCH = .TRUE.
                           DNMAX=COMAX
                           NINT=NINT+1

                           CALL CPU_TIME(TINTS)
                           CALL SHELLMP2(NSTEPS,CURRENTSTEPLENGTH)
                           CALL CPU_TIME(TINTF)
                           SINT=SINT+TINTF-TINTS
                        ENDIF

                     ENDIF

                  ENDDO

                  ! --- BEGIN 1ST TRANSFORMATION ---
                  ! WE USE I, J INDICATES OCCUPPIED OR CORRELATED ORBITALS
                  ! WHILE USE A, B INDICTAES VIRTUAL ORBITALS
                  ! K,L,M,N INDICATES BASIS SETS
                  ! IN PROGRAM: LL -> K
                  !             KK -> L
                  !             II -> N
                  !             JJ -> M
                  ! ONCE ERI (KL|MN) IS GENERATED, IT WILL BE USED TO OBTAIN
                  !  (IL|MN) = SUM (KL|MN)*C(I, K)
                  !  (IL|MN) AND (IL|NM), M<->N SYMMETRY IS EXPOILED
                  ! (IL|MN) COORESPONDS TO ORBMP2I331( I, L, IPRIM, LPRIM, 1)
                  ! (IL|NM) CORRESPONDS TO ORBMP2I331( I, L, IPRIM, LPRIM, 2)
                  ! --- BEGIN L LOOP ---
                  DO LL=KK,JSHELL
                     LL111=QUICK_BASIS%KSUMTYPE(LL)+QUICK_BASIS%QSBASIS(LL,QUICK_BASIS%QSTART(LL))
                     LL112=QUICK_BASIS%KSUMTYPE(LL)+QUICK_BASIS%QFBASIS(LL,QUICK_BASIS%QFINAL(LL))

                     CALL  CPU_TIME(T1QTS)


                     TESTCUTOFF = YCUTOFF(II,JJ)*YCUTOFF(KK,LL)

                     IF(DABS(TESTCUTOFF).GT.CUTOFFMP2)THEN
                        DO III = II111, II112
                           IIINEW = III-II111+1
                           DO JJJ = MAX(III,JJ111), JJ112
                              JJJNEW = JJJ-JJ111+1
                              DO KKK = KK111, KK112
                                 KKKNEW = KKK-KK111+1

                                 DO LLL = MAX(KKK,LL111), LL112

                                    Y = AO(LLL, IIINEW, JJJNEW, KKKNEW)

                                    IF (DABS(Y).GT.QUICK_METHOD%INTEGRALCUTOFF) THEN
                                       DO I3MP2=1,CURRENTSTEPLENGTH
                                          I=NSTEPS+I3MP2-1
                                          ATEMP=QUICK_QM_STRUCT%COOCC(I,KKK)*Y
                                          BTEMP=QUICK_QM_STRUCT%COOCC(I,LLL)*Y


                                          ORBMP2I331(I3MP2,LLL,IIINEW,JJJNEW,1)= &
                                                ORBMP2I331(I3MP2,LLL,IIINEW,JJJNEW,1)+ATEMP
                                          IF(JJJ.NE.III)THEN
                                             ORBMP2I331(I3MP2,LLL,JJJNEW,IIINEW,2)= &
                                                   ORBMP2I331(I3MP2,LLL,JJJNEW,IIINEW,2)+ATEMP
                                          ENDIF
                                          IF(KKK.NE.LLL)THEN
                                             ORBMP2I331(I3MP2,KKK,IIINEW,JJJNEW,1)= &
                                                   ORBMP2I331(I3MP2,KKK,IIINEW,JJJNEW,1)+BTEMP
                                             IF(III.NE.JJJ)THEN
                                                ORBMP2I331(I3MP2,KKK,JJJNEW,IIINEW,2)= &
                                                      ORBMP2I331(I3MP2,KKK,JJJNEW,IIINEW,2)+BTEMP
                                             ENDIF
                                          ENDIF
                                       ENDDO
                                    ENDIF
                                 ENDDO
                              ENDDO
                           ENDDO
                        ENDDO

                     ENDIF
                     CALL CPU_TIME(T1QTF)
                     S1QT=S1QT+T1QTF-T1QTS

                  ENDDO

                  ! --- END L LOOP
               ENDDO
               ! --- END K LOOP

               IF(PASSCH) THEN
                  DO III=II111,II112
                     IIINEW=III-II111+1

                     DO JJJ=MAX(III,JJ111),JJ112
                        JJJNEW=JJJ-JJ111+1

                        CALL CPU_TIME(T3)

                        ! --- BEGIN 2ND TRANSFORMATION
                        ! AFTER GENERATING (IL|MN) AND (IL|NM), 2ND QUARTER TRANSFORMATION IS
                        ! TO GET (IA|MN) AND (IA|MN)
                        !  (IA|MN) = SUM (IL|MN) * C(A,L)
                        ! (IA|MN) COORESPONDS TO ORBMP2J331( I, A, IPRIM, APRIM, 1)
                        ! (IA|NM) COORESPONDS TO ORBMP2J331( I, A, IPRIM, APRIM, 2)
                        DO I = 1, CURRENTSTEPLENGTH
                           DO J = 1, NBASIS
                              TEMP2D1(J , I) = ORBMP2I331( I, J, IIINEW, JJJNEW, 1)
                           ENDDO
                        ENDDO


                        call DGEMM('n', 'n', IVIR, CURRENTSTEPLENGTH, NBASIS, 1.0d0, QUICK_QM_STRUCT%COVIR, &
                              IVIR, TEMP2D1, NBASIS, 0.0d0, TEMP2D2, IVIR)
                        !TEMP2D2 = MATMUL(QUICK_QM_STRUCT%COVIR,TEMP2D1)

                        DO I = 1, CURRENTSTEPLENGTH
                           DO J = 1, IVIR
                              ORBMP2J331( I, J, IIINEW, JJJNEW, 1) =  TEMP2D2( J, I)
                           ENDDO
                        ENDDO

                        IF (III.NE.JJJ) THEN
                           DO I = 1, CURRENTSTEPLENGTH
                              DO J = 1, NBASIS
                                 TEMP2D1(J , I) = ORBMP2I331( I, J, JJJNEW, IIINEW, 2)
                              ENDDO
                           ENDDO
                           !TEMP2D2 = MATMUL( QUICK_QM_STRUCT%COVIR, TEMP2D1)
                           call DGEMM('n', 'n', IVIR, CURRENTSTEPLENGTH, NBASIS, 1.0d0, QUICK_QM_STRUCT%COVIR, &
                                 IVIR, TEMP2D1, NBASIS, 0.0d0, TEMP2D2, IVIR)

                           DO I = 1, CURRENTSTEPLENGTH
                              DO J = 1, IVIR
                                 orbmp2j331( I, J, IIINEW, JJJNEW, 2) =  TEMP2D2( J, I)
                              ENDDO
                           ENDDO
                        ENDIF

                        CALL CPU_TIME(T4)
                        S2=S2+T4-T3

                        CALL CPU_TIME(T5)
                        ! -- BEGIN 3RD TRANSFORMATION
                        ! (IA|JN) += (IA|MN)*C(J,M)
                        !         = SUM ( ORBMP2J331(I, A, IPRIM, APRIM,1) * C(J, M) + ORBMP2J331(I, A, IRPIM, APRIM,2) * C(J, M))
                        ! (IA|JM) += (IA|NM)*C(J,N)
                        ! (IA|JN) CORESPONDS TO ORBMP2K331(I, J, A, N)
                        DO K33=1,NELEC/2

                           ATEMP = QUICK_QM_STRUCT%COOCC(K33, III)
                           ATEMP2= QUICK_QM_STRUCT%COOCC(K33, JJJ)

                           DO J33=1,IVIR
                              DO ICYCLE=1,CURRENTSTEPLENGTH
                                 ORBMP2K331(ICYCLE,K33,J33,JJJ) = ORBMP2K331(ICYCLE,K33,J33,JJJ)+ &
                                       orbmp2j331(ICYCLE,J33, IIINEW, JJJNEW, 1)*ATEMP
                                 IF (III.NE.JJJ) THEN
                                    ORBMP2K331(ICYCLE,K33,J33,III) = ORBMP2K331(ICYCLE,K33,J33,III) + &
                                          orbmp2j331(ICYCLE,J33, IIINEW, JJJNEW, 2)*ATEMP2
                                 ENDIF
                              ENDDO
                           ENDDO
                        ENDDO
                        CALL CPU_TIME(T6)
                        S3=S3+T6-T5
                     ENDDO
                  ENDDO

               ENDIF

            ENDIF

         ENDDO
      ENDDO

      WRITE (IOUTFILE,'(" EFFECT INTEGRALS    =",I8)') NINT


      DO ICYCLE=1,CURRENTSTEPLENGTH
         I3=NSTEPS+ICYCLE-1
         DO K3=I3,NELEC/2

            CALL CPU_TIME(T7)
            DO J3=1,NBASIS-NELEC/2
               DO L3=1,NBASIS-NELEC/2
                  ORBMP2(L3,J3)=0.0D0
                  L3NEW=L3+NELEC/2
                  DO LLL=1,NBASIS
                     ORBMP2(L3,J3)=ORBMP2(L3,J3)+ORBMP2K331(ICYCLE,K3,J3,LLL)*QUICK_QM_STRUCT%CO(LLL,L3NEW)
                  ENDDO
               ENDDO
            ENDDO
            CALL CPU_TIME(T8)
            S4=S4+T8-T7

            !     E(2) = SUM (J)  SUM(AIB) (AI/BJ) * ( 2*(AI/BJ)-(BI/AJ) ) /
            !                              (E(A)+E(B)-E(I)-E(J))
            DO J3=1,NBASIS-NELEC/2
               DO L3=1,NBASIS-NELEC/2
                  IF(K3.GT.I3)THEN
                     QUICK_QM_STRUCT%EMP2 = QUICK_QM_STRUCT%EMP2 + &
                           2.0D0/(QUICK_QM_STRUCT%E(I3)+QUICK_QM_STRUCT%E(K3)-QUICK_QM_STRUCT%E(J3+NELEC/2) &
                           -QUICK_QM_STRUCT%E(L3+NELEC/2)) &
                           *ORBMP2(J3,L3)*(2.0D0*ORBMP2(J3,L3)-ORBMP2(L3,J3))
                  ENDIF
                  IF(K3.EQ.I3)THEN
                     QUICK_QM_STRUCT%EMP2=QUICK_QM_STRUCT%EMP2+1.0D0/(QUICK_QM_STRUCT%E(I3)+QUICK_QM_STRUCT%E(K3) &
                           -QUICK_QM_STRUCT%E(J3+NELEC/2)-QUICK_QM_STRUCT%E(L3+NELEC/2)) &
                           *ORBMP2(J3,L3)*(2.0D0*ORBMP2(J3,L3)-ORBMP2(L3,J3))
                  ENDIF

               ENDDO
            ENDDO
            CALL CPU_TIME(T9)
            S5=S5+T9-T7
         ENDDO
      ENDDO

      CALL CPU_TIME(TIMER_END%TMP2)
      TIMER_CUMER%TMP2=TIMER_END%TMP2-TIMER_BEGIN%TMP2+TIMER_CUMER%TMP2
      WRITE(*,*) " ERI EVALUATION TIME =", SINT
      WRITE(*,*) " 1ST TRANSFORMATION  =", S1QT
      WRITE(*,*) " 2ND TRANSFORMATION  =", S2
      WRITE(*,*) " 3RD TRANSFORMATION  =", S3
      WRITE(*,*) " 4TH TRANSFORMATION  =", S4
      WRITE(*,*) " ENERGY              =", S5

      DEALLOCATE(TEMP2D1, TEMP2D2)
   ENDDO

   WRITE (IOUTFILE,'(" SECOND ORDER ENERGY =",F16.9)') QUICK_QM_STRUCT%EMP2

   WRITE (*,'(" SECOND ORDER ENERGY =",F16.9)') QUICK_QM_STRUCT%EMP2
   WRITE (IOUTFILE,'(" EMP2                =",F16.9)') QUICK_QM_STRUCT%ETOT+QUICK_QM_STRUCT%EMP2
   CALL PRTACT(IOUTFILE,"END MP2 CALCULATION")
   RETURN
#else
    CALL PRTACT(IOUTFILE,"BEGIN MP2 CALCULATION")
    call gpu_upload_calculated(quick_qm_struct%o,quick_qm_struct%co, &
    quick_qm_struct%vec,quick_qm_struct%dense)
    call gpu_upload_energy(quick_qm_struct%E)
    call gpu_upload_cutoff(cutmatrix, quick_method%integralCutoff,quick_method%primLimit,quick_method%DMCutoff)
    call gpu_MP2(QUICK_QM_STRUCT%EMP2)
    CALL PRTACT(IOUTFILE,"END MP2 CALCULATION")
#endif

END SUBROUTINE CALMP2


#ifdef MPI

! ED BROTHERS. NOVEMBER 27, 2001
! XIAO HE. SEPTEMBER 14,2008
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP
SUBROUTINE MPI_CALMP2
   USE ALLMOD
   USE QUICK_GAUSSIAN_CLASS_MODULE
   IMPLICIT DOUBLE PRECISION(A-H,O-Z)

   INCLUDE "MPIF.H"

   DOUBLE PRECISION CUTOFFTEST,TESTTMP,TESTCUTOFF
   DOUBLE PRECISION, ALLOCATABLE:: TEMP4D(:,:,:,:)
   INTEGER II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2
   COMMON /HRRSTORE/II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2

   IF (BMPI) THEN
      !     CALL MPI_BCAST(DENSE,NBASIS*NBASIS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,MPIERROR)
      !     CALL MPI_BCAST(CO,NBASIS*NBASIS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,MPIERROR)
      !     CALL MPI_BCAST(E,NBASIS,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,MPIERROR)
   ENDIF


   CUTOFFMP2=1.0D-8  ! CUTOFF CRITERIA
   QUICK_METHOD%PRIMLIMIT=1.0D-8
   QUICK_QM_STRUCT%EMP2=0.0D0

   ! OCCUPIED AND VIRTUAL ORBITALS NUMBER
   IOCC=NELEC/2
   IVIR=NBASIS-NELEC/2

   ! CALCULATE MEMORY USAGE AND DETERMINE STEPS
   EMEMORYSUM=REAL(IOCC*IVIR*NBASIS*8.0D0/1024.0D0/1024.0D0/1024.0D0)
   IF (MASTER) THEN
      CALL PRTACT(IOUTFILE,"BEGIN MP2 CALCULATION")
      WRITE(IOUTFILE,'("CURRENT MEMORY USAGE=",E12.6,"M")') EMEMORYSUM
   ENDIF

   ! ACTUALLY NSTEP IS STEP LENGTH
   NSTEP=MIN(INT(1.5D0/EMEMORYSUM),NELEC/2)

   ! IF WITH F ORBITAL
   IF(QUICK_METHOD%FFUNXIAO)THEN
      NBASISTEMP=6
   ELSE
      NBASISTEMP=10
   ENDIF

   ! ALLOCATE SOME VARIABLES
   ALLOCATE(MP2SHELL(NBASIS))
   ALLOCATE(ORBMP2(IVIR,IVIR))
   ALLOCATE(ORBMP2I331(NSTEP,NBASIS,NBASISTEMP,NBASISTEMP,2))
   ALLOCATE(ORBMP2J331(NSTEP,IVIR,NBASISTEMP,NBASISTEMP,2))
   ALLOCATE(ORBMP2K331(NSTEP,IOCC,IVIR,NBASIS))
   ALLOCATE(TEMP4D(NSTEP,IOCC,IVIR,NBASIS))

   ! WITH NSTEP(ACUTALLY, IT REPRESETNS STEP LENGHT), WE CAN
   ! HAVE NO. OF STEPS FOR MP2 CALCULATION
   NSTEPMP2=NELEC/2/NSTEP
   NSTEPMP2=NSTEPMP2+1
   IF(NSTEP*(NSTEPMP2-1).EQ.NELEC/2)THEN
      NSTEPMP2=NSTEPMP2-1
   ENDIF
   WRITE(IOUTFILE,'("TOTAL STEP          =",I6)') NSTEPMP2
   ! PRE-STEP FOR DENSITY CUTOFF
   CALL DENSITYCUTOFF

   ! FIRST SAVE COEFFECIENT.
   DO I=1,NBASIS
      DO J=1,NBASIS
         QUICK_SCRATCH%HOLD(I,J)=QUICK_QM_STRUCT%CO(J,I)
      ENDDO
   ENDDO

   TTT=MAXVAL(YCUTOFF) ! MAX VALUE OF YCUTOFF

   DO I3NEW=1,NSTEPMP2               ! STEP COUNTER
      CALL CPU_TIME(TIMER_BEGIN%TMP2)
      NTEMP=0    ! INTEGER COUNTER
      NSTEPMP2S=(I3NEW-1)*NSTEP+1    ! STEP START N
      NSTEPMP2F=I3NEW*NSTEP          ! STEP END N

      IF(I3NEW.EQ.NSTEPMP2)NSTEPMP2F=NELEC/2
      NSTEPLENGTH=NSTEPMP2F-NSTEPMP2S+1  ! STEP LENGH, FROM NSTEPMP2S TO NSTEPMP2F

      ! INITIAL ORBMP2K331
      CALL INITIALORBMP2K331(ORBMP2K331,NSTEP,NBASIS,IVIR,IOCC,NSTEPLENGTH)

      !---------------- MPI/ ALL NODES -----------------------------
      DO I=1,MPI_JSHELLN(MPIRANK)
         II=MPI_JSHELL(MPIRANK,I)
         !     DO II=1,JSHELL
         DO JJ=II,JSHELL

            ! FIRST WE DO INTEGRAL AND SUM THEM PROPERLY
            IF(YCUTOFF(II,JJ).GT.CUTOFFMP2/TTT)THEN
               CALL INITIALORBMP2IJ(ORBMP2I331,NSTEP,NSTEPLENGTH,NBASIS,NBASISTEMP,NBASISTEMP)
               CALL INITIALORBMP2IJ(ORBMP2J331,NSTEP,NSTEPLENGTH,IVIR,NBASISTEMP,NBASISTEMP)

               DO KK=1,JSHELL
                  DO LL=KK,JSHELL

                     ! SCHWARTS CUTOFF IS IMPLEMENTED HERE
                     COMAX=0.D0
                     TESTCUTOFF = YCUTOFF(II,JJ)*YCUTOFF(KK,LL)
                     IF(TESTCUTOFF.GT.CUTOFFMP2)THEN

                        NKK1=QUICK_BASIS%QSTART(KK)
                        NKK2=QUICK_BASIS%QFINAL(KK)
                        NLL1=QUICK_BASIS%QSTART(LL)
                        NLL2=QUICK_BASIS%QFINAL(LL)

                        NBK1=QUICK_BASIS%QSBASIS(KK,NKK1)
                        NBK2=QUICK_BASIS%QFBASIS(KK,NKK2)
                        NBL1=QUICK_BASIS%QSBASIS(LL,NLL1)
                        NBL2=QUICK_BASIS%QFBASIS(LL,NLL2)

                        KK111=QUICK_BASIS%KSUMTYPE(KK)+NBK1
                        KK112=QUICK_BASIS%KSUMTYPE(KK)+NBK2
                        LL111=QUICK_BASIS%KSUMTYPE(LL)+NBL1
                        LL112=QUICK_BASIS%KSUMTYPE(LL)+NBL2

                        ! FIND THE CO-MAX VALUE
                        DO ICYCLE=1,NSTEPLENGTH
                           I3=NSTEPMP2S+ICYCLE-1
                           DO KKK=KK111,KK112
                              DO LLL=MAX(KKK,LL111),LL112
                                 COMAX=MAX(COMAX,DABS(QUICK_QM_STRUCT%CO(KKK,I3)))
                                 COMAX=MAX(COMAX,DABS(QUICK_QM_STRUCT%CO(LLL,I3)))
                              ENDDO
                           ENDDO
                        ENDDO



                        TESTCUTOFF=TESTCUTOFF*COMAX
                        IF(TESTCUTOFF.GT.CUTOFFMP2)THEN
                           DNMAX=COMAX
                           NTEMP=NTEMP+1
                           CALL SHELLMP2(NSTEPMP2S,NSTEPLENGTH)
                        ENDIF

                     ENDIF

                  ENDDO
               ENDDO


               ! NEXT STEP IS TO FOLDING INTEGERS. WITHOUT FOLDING, THE SCALING WILL
               ! BE N^8, AND AFTER FOLDING, SCALING IS REDUCED TO 4N^5

               NII1=QUICK_BASIS%QSTART(II)
               NII2=QUICK_BASIS%QFINAL(II)
               NJJ1=QUICK_BASIS%QSTART(JJ)
               NJJ2=QUICK_BASIS%QFINAL(JJ)

               NBI1=QUICK_BASIS%QSBASIS(II,NII1)
               NBI2=QUICK_BASIS%QFBASIS(II,NII2)
               NBJ1=QUICK_BASIS%QSBASIS(JJ,NJJ1)
               NBJ2=QUICK_BASIS%QFBASIS(JJ,NJJ2)

               II111=QUICK_BASIS%KSUMTYPE(II)+NBI1
               II112=QUICK_BASIS%KSUMTYPE(II)+NBI2
               JJ111=QUICK_BASIS%KSUMTYPE(JJ)+NBJ1
               JJ112=QUICK_BASIS%KSUMTYPE(JJ)+NBJ2

               DO III=II111,II112
                  DO JJJ=MAX(III,JJ111),JJ112

                     IIINEW=III-II111+1
                     JJJNEW=JJJ-JJ111+1

                     ! FOLDING 4 INDICES INTEGERS INTO 3 INDICES
                     DO LLL=1,NBASIS
                        DO J33=1,IVIR
                           J33NEW=J33+NELEC/2
                           ATEMP=QUICK_QM_STRUCT%CO(LLL,J33NEW)
                           DO ICYCLE=1,NSTEPLENGTH
                              ORBMP2J331(ICYCLE,J33,IIINEW,JJJNEW,1)=ORBMP2J331(ICYCLE,J33,IIINEW,JJJNEW,1) + &
                                    ORBMP2I331(ICYCLE,LLL,IIINEW,JJJNEW,1)*ATEMP
                              IF(III.NE.JJJ)THEN
                                 ORBMP2J331(ICYCLE,J33,JJJNEW,IIINEW,2)=ORBMP2J331(ICYCLE,J33,JJJNEW,IIINEW,2) + &
                                       ORBMP2I331(ICYCLE,LLL,JJJNEW,IIINEW,2)*ATEMP
                              ENDIF
                           ENDDO
                        ENDDO
                     ENDDO

                     DO J33=1,IVIR
                        DO K33=1,NELEC/2
                           ATEMP=QUICK_SCRATCH%HOLD(K33,III)
                           ATEMP2=QUICK_SCRATCH%HOLD(K33,JJJ)
                           DO ICYCLE=1,NSTEPLENGTH
                              ORBMP2K331(ICYCLE,K33,J33,JJJ)=ORBMP2K331(ICYCLE,K33,J33,JJJ)+ &
                                    ORBMP2J331(ICYCLE,J33,IIINEW,JJJNEW,1)*ATEMP
                              IF(III.NE.JJJ)THEN
                                 ORBMP2K331(ICYCLE,K33,J33,III)=ORBMP2K331(ICYCLE,K33,J33,III)+ &
                                       ORBMP2J331(ICYCLE,J33,JJJNEW,IIINEW,2)*ATEMP2
                              ENDIF
                           ENDDO
                        ENDDO
                     ENDDO
                  ENDDO
               ENDDO
            ENDIF

         ENDDO
      ENDDO

      CALL CPU_TIME(TIMER_END%TMP2)
      TIMER_CUMER%TMP2=TIMER_END%TMP2-TIMER_BEGIN%TMP2+TIMER_CUMER%TMP2

      ! SEND THE INTEGRAL PACKAGE AFTER FOLDING TO MASTER NODE
      ! AND MASTER NODE WILL RECEIVE THEM AND DO NEXT TWO FOLDING STEPS

      ! SLAVE NODE WILL SEND INFOS
      IF(.NOT.MASTER) THEN
         TEMPE=QUICK_QM_STRUCT%EMP2
         DO J1=1,NBASIS
            DO K1=1,IVIR
               DO I1=1,IOCC
                  DO ICYCLE=1,NSTEPLENGTH
                     TEMP4D(ICYCLE,I1,K1,J1)= ORBMP2K331(ICYCLE,I1,K1,J1)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
         ! SEND 3 INDICES INTEGRALS TO MASTER NODE
         CALL MPI_SEND(TEMP4D,NBASIS*IVIR*IOCC*NSTEPLENGTH,MPI_DOUBLE_PRECISION,0,MPIRANK,MPI_COMM_WORLD,IERROR)
         ! MASTER NODE WILL RECEIVE INFOS FROM EVERY NODES
      ELSE
         DO I=1,MPISIZE-1
            ! RECEIVE INTEGRALS FROM SLAVE NODES
            CALL MPI_RECV(TEMP4D,NBASIS*IVIR*IOCC*NSTEPLENGTH,MPI_DOUBLE_PRECISION,I,I,MPI_COMM_WORLD,MPI_STATUS,IERROR)
            ! AND SUM THEM INTO OPERATOR
            DO J1=1,NBASIS
               DO K1=1,IVIR
                  DO I1=1,IOCC
                     DO ICYCLE=1,NSTEPLENGTH
                        ORBMP2K331(ICYCLE,I1,K1,J1)= ORBMP2K331(ICYCLE,I1,K1,J1) +TEMP4D(ICYCLE,I1,K1,J1)
                     ENDDO
                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDIF

      !---------------- END MPI/ ALL NODES ------------------------

      !---------------- MPI/ MASTER -------------------------------
      IF (MASTER) THEN


         DO ICYCLE=1,NSTEPLENGTH
            I3=NSTEPMP2S+ICYCLE-1
            DO K3=I3,NELEC/2
               ! FOLD 3 INDICES INTEGRAL INTO 2 INDICES
               DO J3=1,NBASIS-NELEC/2
                  DO L3=1,NBASIS-NELEC/2
                     ORBMP2(L3,J3)=0.0D0
                     L3NEW=L3+NELEC/2
                     DO LLL=1,NBASIS
                        ORBMP2(L3,J3)=ORBMP2(L3,J3)+ORBMP2K331(ICYCLE,K3,J3,LLL)*QUICK_QM_STRUCT%CO(LLL,L3NEW)
                     ENDDO
                  ENDDO
               ENDDO

               ! NOW WE CAN GET ENERGY
               DO J3=1,NBASIS-NELEC/2
                  DO L3=1,NBASIS-NELEC/2
                     IF(K3.GT.I3)THEN
                        QUICK_QM_STRUCT%EMP2=QUICK_QM_STRUCT%EMP2+2.0D0/(QUICK_QM_STRUCT%E(I3)+ &
                              QUICK_QM_STRUCT%E(K3)-QUICK_QM_STRUCT%E(J3+NELEC/2)-QUICK_QM_STRUCT%E(L3+NELEC/2)) &
                              *ORBMP2(J3,L3)*(2.0D0*ORBMP2(J3,L3)-ORBMP2(L3,J3))
                     ENDIF
                     IF(K3.EQ.I3)THEN
                        QUICK_QM_STRUCT%EMP2=QUICK_QM_STRUCT%EMP2+1.0D0/(QUICK_QM_STRUCT%E(I3)+ &
                              QUICK_QM_STRUCT%E(K3)-QUICK_QM_STRUCT%E(J3+NELEC/2)-QUICK_QM_STRUCT%E(L3+NELEC/2)) &
                              *ORBMP2(J3,L3)*(2.0D0*ORBMP2(J3,L3)-ORBMP2(L3,J3))
                     ENDIF

                  ENDDO
               ENDDO
            ENDDO
         ENDDO
      ENDIF
      !---------------- ALL MPI/ MASTER ---------------------------

      ! SYNC ALL NODES
      CALL MPI_BARRIER(MPI_COMM_WORLD,MPIERROR)
   ENDDO

   IF (MASTER) THEN
      WRITE (IOUTFILE,'("SECOND ORDER ENERGY=",F16.9)') QUICK_QM_STRUCT%EMP2
      WRITE (IOUTFILE,'("QUICK_QM_STRUCT%EMP2               =",F16.9)') QUICK_QM_STRUCT%ETOT+QUICK_QM_STRUCT%EMP2
      CALL PRTACT(IOUTFILE,"END MP2 CALCULATION")
   ENDIF
   RETURN
END SUBROUTINE MPI_CALMP2
#endif

! ED BROTHERS. NOVEMBER 27, 2001
! XIAO HE. SEPTEMBER 14,2008
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

SUBROUTINE CALMP2DIVCON
   USE ALLMOD
   USE QUICK_GAUSSIAN_CLASS_MODULE
   IMPLICIT DOUBLE PRECISION(A-H,O-Z)

   LOGICAL LOCALLOG1,LOCALLOG2

   DOUBLE PRECISION XIAOTEST,TESTTMP
   INTEGER II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2
   COMMON /HRRSTORE/II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2

   IS = 0
   QUICK_BASIS%FIRST_SHELL_BASIS_FUNCTION(1) = 1
   DO II=1,NATOM-1
      IS=IS+QUICK_BASIS%KSHELL(QUICK_MOLSPEC%IATTYPE(II))
      QUICK_BASIS%LAST_SHELL_BASIS_FUNCTION(II)=IS
      QUICK_BASIS%FIRST_SHELL_BASIS_FUNCTION(II+1) = IS+1
   ENDDO
   QUICK_BASIS%LAST_SHELL_BASIS_FUNCTION(NATOM) = NSHELL

   ALLOCATE(WTOSPOINT(NP,NBASIS))
   CALL WTOSCORR

   XIAOCUTOFFMP2=1.0D-7

   QUICK_QM_STRUCT%EMP2=0.0D0
   EMP2TEMP=0.0D0
   NSTEP=1


   DO ITT=1,NP

      DO I=1,NBASISDC(ITT)
         DO J=1,NBASISDC(ITT)
            QUICK_QM_STRUCT%CO(I,J)=CODCSUB(I,J,ITT)
         ENDDO
      ENDDO

      DO I=1,NBASISDC(ITT)
         DO J=1,NBASISDC(ITT)
            QUICK_SCRATCH%HOLD(I,J)=QUICK_QM_STRUCT%CO(J,I)
         ENDDO
      ENDDO


      DO K3=1,2
         DO J3=1,IVIR
            DO L3=1,IVIR
               L3NEW=J3+IOCC
               DO LLL=1,NBASISDC(ITT)
                  WRITE(10,*) K3,J3,LLL,L3NEW,ORBMP2K331(1,K3,J3,LLL),QUICK_QM_STRUCT%CO(LLL,L3NEW)
               ENDDO
            ENDDO
         ENDDO
      ENDDO

      TTT=0.0D0

      DO IIAT=1,DCSUBN(ITT)
         IIATOM=DCSUB(ITT,IIAT)
         DO II=QUICK_BASIS%FIRST_SHELL_BASIS_FUNCTION(IIATOM),QUICK_BASIS%LAST_SHELL_BASIS_FUNCTION(IIATOM)
            DO JJAT=IIAT,DCSUBN(ITT)
               JJATOM=DCSUB(ITT,JJAT)
               DO JJ=MAX(QUICK_BASIS%FIRST_SHELL_BASIS_FUNCTION(JJATOM),II),QUICK_BASIS%LAST_SHELL_BASIS_FUNCTION(JJATOM)
                  TESTTMP=YCUTOFF(II,JJ)
                  TTT=MAX(TTT,TESTTMP)
               ENDDO
            ENDDO
         ENDDO
      ENDDO


      ! DETERMINE THE ELECTRONS FOR SUBSYSTEM
      IF(MOD(NELECMP2SUB(ITT),2).EQ.1)THEN
         NELECMP2SUB(ITT)=NELECMP2SUB(ITT)+1
      ENDIF

      ! DETERMINE THE OCCUPIED AND VIRTUAL ORBITALS
      IF(MOD(NELECMP2SUB(ITT),2).EQ.1)THEN
         IOCC=NELECMP2SUB(ITT)/2+1
         IVIR=NBASISDC(ITT)-IOCC+1
      ELSE
         IOCC=NELECMP2SUB(ITT)/2
         IVIR=NBASISDC(ITT)-IOCC
      ENDIF
      WRITE(*,*) IOCC,IVIR
      ! WITH F ORBITAL
      IF (QUICK_METHOD%FFUNXIAO) THEN
         NBASISTEMP=6
      ELSE
         NBASISTEMP=10
      ENDIF

      ! ALLOCATE VARIBLES
      ALLOCATE(ORBMP2(100,100))
      ALLOCATE(ORBMP2DCSUB(IOCC,IVIR,IVIR))
      ALLOCATE(ORBMP2I331(NSTEP,NBASISDC(ITT),NBASISTEMP,NBASISTEMP,2))
      ALLOCATE(ORBMP2J331(NSTEP,IVIR,NBASISTEMP,NBASISTEMP,2))
      ALLOCATE(ORBMP2K331(NSTEP,IOCC,IVIR,NBASISDC(ITT)))
      ALLOCATE(ORBMP2K331DCSUB(IOCC,IVIR,NBASISDC(ITT)))

      ! SCHWARTZ CUTOFF IS IMPLEMENTED HERE. (AB|CD)**2<=(AB|AB)*(CD|CD)
      ! REFERENCE: STROUT DL AND SCUSERIA JCP 102(1995),8448.

      DO I3=1,IOCC
         DO L1=1,IVIR
            DO K1=1,IVIR
               ORBMP2(J1,K1)=0.0D0
            ENDDO
         ENDDO

         CALL INITIALORBMP2K331(ORBMP2K331,NSTEP,NBASISDC(ITT),IVIR,IOCC,NSTEP)

         DO L1=1,IVIR
            DO K1=1,IVIR
               DO J1=1,IOCC
                  ORBMP2DCSUB(J1,K1,L1)=0.0D0
               ENDDO
            ENDDO
         ENDDO

         DO J1=1,NBASISDC(ITT)
            DO K1=1,IVIR
               DO I1=1,IOCC
                  ORBMP2K331DCSUB(I1,K1,J1)=0.0D0
               ENDDO
            ENDDO
         ENDDO

         NTEMP=0

         DO IIAT=1,DCSUBN(ITT)
            IIATOM=DCSUB(ITT,IIAT)

            DO JJAT=IIAT,DCSUBN(ITT)
               JJATOM=DCSUB(ITT,JJAT)

               ! IIATOM AND JJATOM IS THE ATOM OF THE SUBSYSTEM
               ! FIRST, WE NEED TO FIGURE OUT WHICH SHELL SHOULD WE CALCLULATE
               ! WHICH IS IISTART1 AND IISTART2 FOR I AND JJ IS FROM JJSTART1 TO JJSTART2
               IF(IIATOM.LE.JJATOM)THEN
                  IISTART1=QUICK_BASIS%FIRST_SHELL_BASIS_FUNCTION(IIATOM)
                  IISTART2=QUICK_BASIS%LAST_SHELL_BASIS_FUNCTION(IIATOM)
                  JJSTART1=QUICK_BASIS%FIRST_SHELL_BASIS_FUNCTION(JJATOM)
                  JJSTART2=QUICK_BASIS%LAST_SHELL_BASIS_FUNCTION(JJATOM)
               ENDIF
               IF(IIATOM.GT.JJATOM)THEN
                  JJSTART1=QUICK_BASIS%FIRST_SHELL_BASIS_FUNCTION(IIATOM)
                  JJSTART2=QUICK_BASIS%LAST_SHELL_BASIS_FUNCTION(IIATOM)
                  IISTART1=QUICK_BASIS%FIRST_SHELL_BASIS_FUNCTION(JJATOM)
                  IISTART2=QUICK_BASIS%LAST_SHELL_BASIS_FUNCTION(JJATOM)
               ENDIF

               DO II=IISTART1,IISTART2
                  DO JJ=MAX(JJSTART1,II),JJSTART2

                     TESTTMP=YCUTOFF(II,JJ)
                     IF(TESTTMP.GT.XIAOCUTOFFMP2/TTT)THEN

                        CALL INITIALORBMP2IJ(ORBMP2I331,NSTEP,NSTEP,NBASISDC(ITT),NBASISTEMP,NBASISTEMP)
                        CALL INITIALORBMP2IJ(ORBMP2J331,NSTEP,NSTEP,IVIR,NBASISTEMP,NBASISTEMP)

                        ! NOW WE WILL DETERMINE K SHELL AND L SHELL, THE LAST TWO INDICES
                        DO KKAT=1,DCSUBN(ITT)
                           KKATOM=DCSUB(ITT,KKAT)
                           DO LLAT=KKAT,DCSUBN(ITT)
                              LLATOM=DCSUB(ITT,LLAT)
                              IF(KKATOM.LE.LLATOM)THEN
                                 KKSTART1=QUICK_BASIS%FIRST_SHELL_BASIS_FUNCTION(KKATOM)
                                 KKSTART2=QUICK_BASIS%LAST_SHELL_BASIS_FUNCTION(KKATOM)
                                 LLSTART1=QUICK_BASIS%FIRST_SHELL_BASIS_FUNCTION(LLATOM)
                                 LLSTART2=QUICK_BASIS%LAST_SHELL_BASIS_FUNCTION(LLATOM)
                              ENDIF
                              IF(KKATOM.GT.LLATOM)THEN
                                 LLSTART1=QUICK_BASIS%FIRST_SHELL_BASIS_FUNCTION(KKATOM)
                                 LLSTART2=QUICK_BASIS%LAST_SHELL_BASIS_FUNCTION(KKATOM)
                                 KKSTART1=QUICK_BASIS%FIRST_SHELL_BASIS_FUNCTION(LLATOM)
                                 KKSTART2=QUICK_BASIS%LAST_SHELL_BASIS_FUNCTION(LLATOM)
                              ENDIF
                              DO KK=KKSTART1,KKSTART2
                                 DO LL=MAX(LLSTART1,KK),LLSTART2
                                    COMAX=0.D0
                                    XIAOTEST1 = YCUTOFF(II,JJ)*YCUTOFF(KK,LL)
                                    IF(XIAOTEST1.GT.XIAOCUTOFFMP2)THEN

                                       NKK1=QUICK_BASIS%QSTART(KK)
                                       NKK2=QUICK_BASIS%QFINAL(KK)
                                       NLL1=QUICK_BASIS%QSTART(LL)
                                       NLL2=QUICK_BASIS%QFINAL(LL)

                                       NBK1=QUICK_BASIS%QSBASIS(KK,NKK1)
                                       NBK2=QUICK_BASIS%QFBASIS(KK,NKK2)
                                       NBL1=QUICK_BASIS%QSBASIS(LL,NLL1)
                                       NBL2=QUICK_BASIS%QFBASIS(LL,NLL2)

                                       KK111=QUICK_BASIS%KSUMTYPE(KK)+NBK1
                                       KK112=QUICK_BASIS%KSUMTYPE(KK)+NBK2
                                       LL111=QUICK_BASIS%KSUMTYPE(LL)+NBL1
                                       LL112=QUICK_BASIS%KSUMTYPE(LL)+NBL2

                                       DO KKK=KK111,KK112
                                          DO LLL=MAX(KKK,LL111),LL112

                                             COMAX=MAX(COMAX,DABS(QUICK_QM_STRUCT%CO(WTOSPOINT(ITT,KKK),I3)))
                                             COMAX=MAX(COMAX,DABS(QUICK_QM_STRUCT%CO(WTOSPOINT(ITT,LLL),I3)))

                                          ENDDO
                                       ENDDO

                                       XIAOTEST=XIAOTEST1*COMAX
                                       IF(XIAOTEST.GT.XIAOCUTOFFMP2)THEN
                                          NTEMP=NTEMP+1
                                          CALL SHELLMP2DIVCON(I3,ITT)
                                       ENDIF

                                    ENDIF

                                 ENDDO
                              ENDDO

                           ENDDO
                        ENDDO


                        NII1=QUICK_BASIS%QSTART(II)
                        NII2=QUICK_BASIS%QFINAL(II)
                        NJJ1=QUICK_BASIS%QSTART(JJ)
                        NJJ2=QUICK_BASIS%QFINAL(JJ)

                        NBI1=QUICK_BASIS%QSBASIS(II,NII1)
                        NBI2=QUICK_BASIS%QFBASIS(II,NII2)
                        NBJ1=QUICK_BASIS%QSBASIS(JJ,NJJ1)
                        NBJ2=QUICK_BASIS%QFBASIS(JJ,NJJ2)

                        II111=QUICK_BASIS%KSUMTYPE(II)+NBI1
                        II112=QUICK_BASIS%KSUMTYPE(II)+NBI2
                        JJ111=QUICK_BASIS%KSUMTYPE(JJ)+NBJ1
                        JJ112=QUICK_BASIS%KSUMTYPE(JJ)+NBJ2

                        DO III=II111,II112
                           DO JJJ=MAX(III,JJ111),JJ112

                              IIINEW=III-II111+1
                              JJJNEW=JJJ-JJ111+1

                              DO LLL=1,NBASISDC(ITT)
                                 DO J33=1,IVIR
                                    J33NEW=J33+IOCC
                                    IF(MOD(NELECMP2SUB(ITT),2).EQ.1)J33NEW=J33+IOCC-1
                                    ATEMP=QUICK_QM_STRUCT%CO(LLL,J33NEW)
                                    ORBMP2J331(NSTEP,J33,IIINEW,JJJNEW,1)=ORBMP2J331(NSTEP,J33,IIINEW,JJJNEW,1) + &
                                          ORBMP2I331(NSTEP,LLL,IIINEW,JJJNEW,1)*ATEMP
                                    IF(III.NE.JJJ)THEN
                                       ORBMP2J331(NSTEP,J33,JJJNEW,IIINEW,2)=ORBMP2J331(NSTEP,J33,JJJNEW,IIINEW,2) + &
                                             ORBMP2I331(NSTEP,LLL,JJJNEW,IIINEW,2)*ATEMP
                                    ENDIF
                                 ENDDO
                              ENDDO

                              DO J33=1,IVIR
                                 DO K33=1,IOCC
                                    ORBMP2K331(NSTEP,K33,J33,WTOSPOINT(ITT,JJJ))=ORBMP2K331(NSTEP,K33,J33,WTOSPOINT(ITT,JJJ))+ &
                                          ORBMP2J331(NSTEP,J33,IIINEW,JJJNEW,1)*QUICK_SCRATCH%HOLD(K33,WTOSPOINT(ITT,III))
                                    IF(III.NE.JJJ)THEN
                                       ORBMP2K331(NSTEP,K33,J33,WTOSPOINT(ITT,III))=ORBMP2K331(NSTEP,K33,J33,WTOSPOINT(ITT,III))+ &
                                             ORBMP2J331(NSTEP,J33,JJJNEW,IIINEW,2)*QUICK_SCRATCH%HOLD(K33,WTOSPOINT(ITT,JJJ))
                                    ENDIF
                                 ENDDO
                              ENDDO

                              LOCALLOG1=.FALSE.
                              LOCALLOG2=.FALSE.

                              DO IIATDC=1,DCCOREN(ITT)
                                 IIATOMDC=DCCORE(ITT,IIATDC)
                                 DO IINBASISDC=QUICK_BASIS%FIRST_BASIS_FUNCTION(IIATOMDC),QUICK_BASIS%LAST_BASIS_FUNCTION(IIATOMDC)
                                    IF(III.EQ.IINBASISDC)LOCALLOG1=.TRUE.
                                    IF(JJJ.EQ.IINBASISDC)LOCALLOG2=.TRUE.
                                 ENDDO
                              ENDDO

                              IF(LOCALLOG1)THEN
                                 DO J33=1,IVIR
                                    DO K33=1,IOCC
                                       ORBMP2K331DCSUB(K33,J33,WTOSPOINT(ITT,JJJ))=ORBMP2K331DCSUB(K33,J33,WTOSPOINT(ITT,JJJ))+ &
                                             ORBMP2J331(NSTEP,J33,IIINEW,JJJNEW,1)*QUICK_SCRATCH%HOLD(K33,WTOSPOINT(ITT,III))
                                    ENDDO
                                 ENDDO
                              ENDIF

                              IF(LOCALLOG2.AND.III.NE.JJJ)THEN
                                 DO J33=1,IVIR
                                    DO K33=1,IOCC
                                       ORBMP2K331DCSUB(K33,J33,WTOSPOINT(ITT,III))=ORBMP2K331DCSUB(K33,J33,WTOSPOINT(ITT,III))+ &
                                             ORBMP2J331(NSTEP,J33,JJJNEW,IIINEW,2)*QUICK_SCRATCH%HOLD(K33,WTOSPOINT(ITT,JJJ))
                                    ENDDO
                                 ENDDO
                              ENDIF

                           ENDDO
                        ENDDO

                     ENDIF

                  ENDDO
               ENDDO

            ENDDO
         ENDDO
         WRITE (IOUTFILE,*)"NTEMP=",NTEMP

         DO K3=1,2
            DO J3=1,IVIR
               DO L3=1,IVIR
                  L3NEW=J3+IOCC
                  DO LLL=1,NBASISDC(ITT)
                     WRITE(10,*) K3,J3,LLL,L3NEW,ORBMP2K331(1,K3,J3,LLL),QUICK_QM_STRUCT%CO(LLL,L3NEW)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

         DO LLL=1,NBASISDC(ITT)
            DO J3=1,IVIR
               DO L3=1,IVIR
                  L3NEW=L3+IOCC
                  ORBMP2(L3,J3)=0.0D0
                  IF(MOD(NELECMP2SUB(ITT),2).EQ.1)L3NEW=L3+IOCC-1
                  DO K3=1,IOCC
                     ORBMP2(L3,J3)=ORBMP2(L3,J3)+ORBMP2K331(NSTEP,K3,J3,LLL)*QUICK_SCRATCH%HOLD(L3NEW,LLL)
                     ORBMP2DCSUB(K3,L3,J3)=ORBMP2DCSUB(K3,L3,J3)+ORBMP2K331DCSUB(K3,J3,LLL)*QUICK_SCRATCH%HOLD(L3NEW,LLL)
                  ENDDO
               ENDDO
            ENDDO
         ENDDO

         IF(MOD(NELECMP2SUB(ITT),2).EQ.0)THEN
            DO L=1,IVIR
               DO K=1,IOCC
                  DO J=1,IVIR
                     QUICK_QM_STRUCT%EMP2=QUICK_QM_STRUCT%EMP2+1.0D0/(EVALDCSUB(ITT,I3)+EVALDCSUB(ITT,K) &
                           -EVALDCSUB(ITT,J+NELECMP2SUB(ITT)/2)-EVALDCSUB(ITT,L+NELECMP2SUB(ITT)/2)) &
                           *ORBMP2DCSUB(K,J,L)*(2.0D0*ORBMP2(J,L)-ORBMP2(L,J))
                  ENDDO
               ENDDO
            ENDDO
         ENDIF

         IF(MOD(NELECMP2SUB(ITT),2).EQ.1)THEN
            DO L=1,IVIR
               DO K=1,IOCC
                  DO J=1,IVIR
                     QUICK_QM_STRUCT%EMP2=QUICK_QM_STRUCT%EMP2+1.0D0/(EVALDCSUB(ITT,I3)+EVALDCSUB(ITT,K) &
                           -EVALDCSUB(ITT,J+NELECMP2SUB(ITT)/2)-EVALDCSUB(ITT,L+NELECMP2SUB(ITT)/2)) &
                           *ORBMP2DCSUB(K,J,L)*(2.0D0*ORBMP2(J,L)-ORBMP2(L,J))
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ENDDO
      WRITE(IOUTFILE,*) ITT,QUICK_QM_STRUCT%EMP2,QUICK_QM_STRUCT%EMP2-EMP2TEMP

      EMP2TEMP=QUICK_QM_STRUCT%EMP2

      !     DEALLOCATE(MP2SHELL)
      DEALLOCATE(ORBMP2)
      DEALLOCATE(ORBMP2I331)
      DEALLOCATE(ORBMP2J331)
      DEALLOCATE(ORBMP2K331)
      DEALLOCATE(ORBMP2DCSUB)
      DEALLOCATE(ORBMP2K331DCSUB)

   ENDDO

   999 RETURN
END SUBROUTINE CALMP2DIVCON


! INITIAL ORBMP2K331
SUBROUTINE  INITIALORBMP2K331(ORBMP2K331,NSTEP,NBASIS,IVIR,IOCC,NSTEPLENGTH)
   INTEGER NBASIS,IVIR,IOCC,NSTEPLENGTH,I1,JK1,J1,ICYCLE,NSTEP
   DOUBLE PRECISION ORBMP2K331(NSTEP,IOCC,IVIR,NBASIS)
   DO J1=1,NBASIS
      DO K1=1,IVIR
         DO I1=1,IOCC
            DO ICYCLE=1,NSTEPLENGTH
               ORBMP2K331(ICYCLE,I1,K1,J1)=0.0D0
            ENDDO
         ENDDO
      ENDDO
   ENDDO
END SUBROUTINE INITIALORBMP2K331

! INITIAL ORBMP2I331 AND ORBMP2J331
SUBROUTINE INITIALORBMP2IJ(ORBMP2I331,NSTEP,NSTEPLENGTH,NBASIS,NBASISTEMP,NBASISTEMP2)
   INTEGER NSTEP,NSTEPLENGTH,NBASIS,NBASISTEMP,NBASISTEMP2
   INTEGER L1,J1,I1,K1,ICYCLE
   DOUBLE PRECISION ORBMP2I331(NSTEP,NBASIS,NBASISTEMP,NBASISTEMP2,2)

   DO L1=1,2
      DO J1=1,NBASISTEMP
         DO I1=1,NBASISTEMP2
            DO K1=1,NBASIS
               DO ICYCLE=1,NSTEPLENGTH
                  ORBMP2I331(ICYCLE,K1,I1,J1,L1)=0.0D0
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDDO
END SUBROUTINE INITIALORBMP2IJ



SUBROUTINE INITIALORBMP2J(ORBMP2J33, IVIR, CURRENTSTEPLENGTH)
   IMPLICIT NONE
   INTEGER I,J, IVIR, CURRENTSTEPLENGTH
   DOUBLE PRECISION ORBMP2J33(CURRENTSTEPLENGTH, IVIR, 2)
   DO I = 1, IVIR
      DO J = 1, CURRENTSTEPLENGTH
         ORBMP2J33(J,I,1)=0.0
         ORBMP2J33(J,I,2)=0.0
      ENDDO
   ENDDO
END SUBROUTINE INITIALORBMP2J