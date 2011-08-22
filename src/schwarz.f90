! change call g2eshell in hfgrad.f? Answer:NO
! change shreshold!!!!!!!
!!!!! change hfoperatordelta.f hfenergy2eshell.f hfenergyshell.f
! change schwartz.f hfoperator.f 2eshell.f hfgrad.f 2eshellopt.f
! do the cutoff every OPT step
! Xiao HE 07/07/07
! Schwartz cutoff is implemented here. (ab|cd)**2<=(ab|ab)*(cd|cd)
! Reference: Strout DL and Scuseria JCP 102(1995),8448.

subroutine schwarzoff
  use allmod

  Implicit real*8(a-h,o-z)
  do II=1,nshell
     do JJ=II,nshell
        !     Ymaxtemp=0.0d0
        call shellcutoff(II,JJ,Ymaxtemp)
        !        print*,Ymaxtemp
        Ycutoff(II,JJ)=dsqrt(Ymaxtemp)
        Ycutoff(JJ,II)=dsqrt(Ymaxtemp)
     enddo
  enddo

end subroutine schwarzoff



subroutine shellcutoff(II,JJ,Ymax)
  use allmod

  Implicit real*8(a-h,o-z)
  real*8 P(3),Q(3),W(3),KAB,KCD
  Parameter(NN=13)
  real*8 FM(0:13)
  real*8 RA(3),RB(3),RC(3),RD(3)

  real*8 Qtemp(3),WQtemp(3),CDtemp,ABcom,Ptemp(3),WPtemp(3),ABtemp,CDcom,ABCDtemp
  ! integer II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2
  ! common /hrrstore/II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2

  COMMON /VRRcom/Qtemp,WQtemp,CDtemp,ABcom,Ptemp,WPtemp,ABtemp,CDcom,ABCDtemp

  COMMON /COM1/RA,RB,RC,RD

  KK=II
  LL=JJ

  Ymax=0.0d0

  do M=1,3
     RA(M)=xyz(M,quick_basis%katom(II))
     RB(M)=xyz(M,quick_basis%katom(JJ))
     RC(M)=xyz(M,quick_basis%katom(KK))
     RD(M)=xyz(M,quick_basis%katom(LL))
  enddo

  NII1=quick_basis%Qstart(II)
  NII2=quick_basis%Qfinal(II)
  NJJ1=quick_basis%Qstart(JJ)
  NJJ2=quick_basis%Qfinal(JJ)
  NKK1=quick_basis%Qstart(KK)
  NKK2=quick_basis%Qfinal(KK)
  NLL1=quick_basis%Qstart(LL)
  NLL2=quick_basis%Qfinal(LL)

  NNAB=(NII2+NJJ2)
  NNCD=(NKK2+NLL2)

  NABCDTYPE=NNAB*10+NNCD

  NNAB=sumindex(NNAB)
  NNCD=sumindex(NNCD)

  NNA=Sumindex(NII1-1)+1

  NNC=Sumindex(NKK1-1)+1

  NABCD=NII2+NJJ2+NKK2+NLL2
  ITT=0
  !ITTprim=0

  do JJJ=1,quick_basis%kprim(JJ)
     Nprij=quick_basis%kstart(JJ)+JJJ-1
     do III=1,quick_basis%kprim(II)
        Nprii=quick_basis%kstart(II)+III-1
        AB=Apri(Nprii,Nprij)
        ABtemp=0.5d0/AB
        do M=1,3
           P(M)=Ppri(M,Nprii,Nprij)
           Ptemp(M)=P(M)-RA(M)
        enddo
        !            KAB=Kpri(Nprii,Nprij)
        do LLL=1,quick_basis%kprim(LL)
           Npril=quick_basis%kstart(LL)+LLL-1
           do KKK=1,quick_basis%kprim(KK)
              Nprik=quick_basis%kstart(KK)+KKK-1
              CD=Apri(Nprik,Npril)
              ABCD=AB+CD
              ROU=AB*CD/ABCD
              RPQ=0.0d0
              ABCDxiao=dsqrt(ABCD)    

              CDtemp=0.5d0/CD
              ABcom=AB/ABCD
              CDcom=CD/ABCD
              ABCDtemp=0.5d0/ABCD

              do M=1,3
                 Q(M)=Ppri(M,Nprik,Npril)
                 W(M)=(P(M)*AB+Q(M)*CD)/ABCD
                 XXXtemp=P(M)-Q(M)
                 RPQ=RPQ+XXXtemp*XXXtemp
                 Qtemp(M)=Q(M)-RC(M)
                 WQtemp(M)=W(M)-Q(M)
                 WPtemp(M)=W(M)-P(M)
              enddo
              !                         KCD=Kpri(Nprik,Npril)

              T=RPQ*ROU

              call FmT(NABCD,T,FM)
              do iitemp=0,NABCD
                 Yxiaotemp(1,1,iitemp)=FM(iitemp)/ABCDxiao
              enddo

              ITT=ITT+1

              call vertical(NABCDTYPE)

              do I2=NNC,NNCD
                 do I1=NNA,NNAB
                    Yxiao(ITT,I1,I2)=Yxiaotemp(I1,I2,0)
                 enddo
              enddo

              If(KKK.eq.III.and.JJJ.eq.LLL)then

                 !                             ITTprim=ITTprim+1

                 do I2=NNC,NNCD
                    do I1=NNA,NNAB
                       Yxiaoprim(III,JJJ,I1,I2)=Yxiaotemp(I1,I2,0)
                    enddo
                 enddo

              endif
           enddo
        enddo
     enddo
  enddo


  !      do JJJ=1,quick_basis%kprim(JJ)
  !        Nprij=quick_basis%kstart(JJ)+JJJ-1
  ! do III=1,quick_basis%kprim(II)
  !   Nprii=quick_basis%kstart(II)+III-1

  do IIxiao=1,quick_basis%kprim(II)
     Nprii=quick_basis%kstart(II)+IIxiao-1
     do JJxiao=1,quick_basis%kprim(JJ)
        Nprij=quick_basis%kstart(JJ)+JJxiao-1

        Ymaxprim=0.0d0
        do I=NII1,NII2
           if(I.eq.0)then
              NNA=1
           else
              NNA=Sumindex(I-1)+1
           endif
           do J=NJJ1,NJJ2
              NNAB=SumINDEX(I+J)
              K=I
              L=J
              NNC=NNA
              NNCD=SumIndex(K+L)
              call classprim(I,J,K,L,II,JJ,KK,LL,NNA,NNC,NNAB,NNCD,Ymaxprim,IIxiao,JJxiao)
           enddo
        enddo

        cutprim(Nprii,Nprij)=dsqrt(Ymaxprim)

     enddo
  enddo


  do I=NII1,NII2
     if(I.eq.0)then
        NNA=1
     else
        NNA=Sumindex(I-1)+1
     endif
     

     do J=NJJ1,NJJ2
        NNAB=SumINDEX(I+J)
        !       do K=NKK1,NKK2
        !            if(K.eq.0)then
        !            NNC=1
        !            else
        !            NNC=Sumindex(k-1)+1
        !            endif
        !         do L=NLL1,NLL2
        NNC=NNA
        K=I
        L=J
        NNCD=SumIndex(K+L)
        call classcutoff(I,J,K,L,II,JJ,KK,LL,NNA,NNC,NNAB,NNCD,Ymax)
        !                 call classprim(I,J,K,L,II,JJ,KK,LL,NNA,NNC,NNAB,NNCD,Ymaxprim,ITT)
        !          enddo
        !        enddo
     enddo


  enddo

end subroutine shellcutoff



subroutine classcutoff(I,J,K,L,II,JJ,KK,LL,NNA,NNC,NNAB,NNCD,Ymax)
  use allmod

  Implicit real*8(A-H,O-Z)
  real*8 store(120,120)
  INTEGER NA(3),NB(3),NC(3),ND(3)
  real*8 P(3),Q(3),W(3),KAB,KCD
  Parameter(NN=13)
  real*8 FM(0:13)
  real*8 RA(3),RB(3),RC(3),RD(3)
  real*8 X44(100000)

  real*8 coefangxiaoL(20),coefangxiaoR(20)
  integer angxiaoL(20),angxiaoR(20),numangularL,numangularR

  COMMON /COM1/RA,RB,RC,RD
  COMMON /COM2/AA,BB,CC,DD,AB,CD,ROU,ABCD
  COMMON /COM4/P,Q,W
  COMMON /COM5/FM

  ! integer II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2
  common /xiaostore/store
  ! common /hrrstore/II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2

  ITT=0
  do JJJ=1,quick_basis%kprim(JJ)
     Nprij=quick_basis%kstart(JJ)+JJJ-1
     do III=1,quick_basis%kprim(II)
        Nprii=quick_basis%kstart(II)+III-1
            
        X2=X0*XCoeff(Nprii,Nprij,I,J)
        do LLL=1,quick_basis%kprim(LL)
           Npril=quick_basis%kstart(LL)+LLL-1
           do KKK=1,quick_basis%kprim(KK)
              Nprik=quick_basis%kstart(KK)+KKK-1
              ITT=ITT+1
              X44(ITT)=X2*XCoeff(Nprik,Npril,K,L)
            enddo
        enddo
     enddo
  enddo
  do MM2=NNC,NNCD
     do MM1=NNA,NNAB
        Ytemp=0.0d0
        do itemp=1,ITT
           Ytemp=Ytemp+X44(itemp)*Yxiao(itemp,MM1,MM2)
           !                        print*,'***',X4,Yxiao(MM1,MM2,ITT),store(MM1,MM2)
        enddo
        store(MM1,MM2)=Ytemp
     enddo
  enddo

  NBI1=quick_basis%Qsbasis(II,I)
  NBI2=quick_basis%Qfbasis(II,I)
  NBJ1=quick_basis%Qsbasis(JJ,J)
  NBJ2=quick_basis%Qfbasis(JJ,J)
  NBK1=quick_basis%Qsbasis(KK,K)
  NBK2=quick_basis%Qfbasis(KK,K)
  NBL1=quick_basis%Qsbasis(LL,L)
  NBL2=quick_basis%Qfbasis(LL,L)

  do III=quick_basis%ksumtype(II)+NBI1,quick_basis%ksumtype(II)+NBI2
     do JJJ=quick_basis%ksumtype(JJ)+NBJ1,quick_basis%ksumtype(JJ)+NBJ2
        !            do KKK=quick_basis%ksumtype(KK)+NBK1,quick_basis%ksumtype(KK)+NBK2
        !              do LLL=quick_basis%ksumtype(LL)+NBL1,quick_basis%ksumtype(LL)+NBL2
        KKK=III
        LLL=JJJ

        If((I.eq.0.and.J.eq.0.and.K.eq.0.and.L.eq.0).or. &
             (I.eq.1.and.J.eq.0.and.K.eq.1.and.L.eq.0))then

           do M=1,3
              NA(M)=KLMN(M,III)
              !                     NB(M)=KLMN(M,JJJ)
              NC(M)=KLMN(M,KKK)
              !                     ND(M)=KLMN(M,LLL) 
           enddo

           M1=trans(NA(1),NA(2),NA(3))
           M3=trans(NC(1),NC(2),NC(3))
           Y=store(M1,M3)

        elseif(I.eq.0.and.J.eq.1.and.K.eq.0.and.L.eq.1)then

           do M=1,3
              !                     NA(M)=KLMN(M,III)
              NB(M)=KLMN(M,JJJ)
              !                     NC(M)=KLMN(M,KKK)
              ND(M)=KLMN(M,LLL)
           enddo

           !                    Y=HRR(NA,NB,NC,ND)
           !                    Y=Y*cons(III)*cons(JJJ)*cons(KKK)*cons(LLL)

           M1=trans(NB(1),NB(2),NB(3))
           M3=trans(ND(1),ND(2),ND(3))

           do itemp=1,3
              if(ND(itemp).ne.0)then
                 ctemp=(RC(itemp)-RD(itemp))
                 Y1=store(M1,M3)+ctemp*store(M1,1)
                 Y2=store(1,M3)+ctemp*store(1,1)
                 goto 117
              endif
           enddo
117        continue

           do jtemp=1,3
              if(NB(jtemp).ne.0)then
                 Y=Y1+(RA(jtemp)-RB(jtemp))*Y2
                 goto 118
              endif
           enddo
118        continue

        elseif(I.eq.1.and.J.eq.1.and.K.eq.1.and.L.eq.1)then

           do M=1,3
              NA(M)=KLMN(M,III)
              NB(M)=KLMN(M,JJJ)
              NC(M)=KLMN(M,KKK)
              ND(M)=KLMN(M,LLL)
           enddo

           !                    Y=HRR(NA,NB,NC,ND)
           !                    Y=Y*cons(III)*cons(JJJ)*cons(KKK)*cons(LLL)

           MA=trans(NA(1),NA(2),NA(3))
           MAB=trans(NA(1)+NB(1),NA(2)+NB(2),NA(3)+NB(3))
           MCX=trans(NC(1),NC(2),NC(3))
           MCD=trans(NC(1)+ND(1),NC(2)+ND(2),NC(3)+ND(3))

           do itemp=1,3
              if(ND(itemp).ne.0)then
                 ctemp=(RC(itemp)-RD(itemp))
                 Y1=store(MAB,MCD)+ctemp*store(MAB,MCX)
                 Y2=store(MA,MCD)+ctemp*store(MA,MCX)
                 goto 141
              endif
           enddo
141        continue

           do jtemp=1,3
              if(NB(jtemp).ne.0)then
                 Y=Y1+(RA(jtemp)-RB(jtemp))*Y2
                 goto 142
              endif
           enddo
142        continue

        elseif((I.eq.2.and.J.eq.0.and.K.eq.2.and.L.eq.0).or. &
             (I.eq.0.and.J.eq.2.and.K.eq.0.and.L.eq.2).or. &
             (I.eq.2.and.J.eq.1.and.K.eq.2.and.L.eq.1).or. &
             (I.eq.1.and.J.eq.2.and.K.eq.1.and.L.eq.2).or. &
             (I.eq.2.and.J.eq.2.and.K.eq.2.and.L.eq.2).or. &
             (I.eq.3.and.J.eq.0.and.K.eq.3.and.L.eq.0).or. &
             (I.eq.0.and.J.eq.3.and.K.eq.0.and.L.eq.3).or. &
             (I.eq.3.and.J.eq.1.and.K.eq.3.and.L.eq.1).or. &
             (I.eq.1.and.J.eq.3.and.K.eq.1.and.L.eq.3).or. &
             (I.eq.3.and.J.eq.2.and.K.eq.3.and.L.eq.2).or. &
             (I.eq.2.and.J.eq.3.and.K.eq.2.and.L.eq.3).or. &
             (I.eq.3.and.J.eq.3.and.K.eq.3.and.L.eq.3))then

           IJtype=10*I+J
           KLtype=10*K+L

           call lefthrr(RA,RB,KLMN(1:3,III),KLMN(1:3,JJJ),IJtype,coefangxiaoL,angxiaoL,numangularL)
           call lefthrr(RC,RD,KLMN(1:3,KKK),KLMN(1:3,LLL),KLtype,coefangxiaoR,angxiaoR,numangularR)

           Y=0.0d0
           do ixiao=1,numangularL
              do jxiao=1,numangularR
                 Y=Y+coefangxiaoL(ixiao)*coefangxiaoR(jxiao)* &
                      store(angxiaoL(ixiao),angxiaoR(jxiao))
              enddo
           enddo

           Y=Y*cons(III)*cons(JJJ)*cons(KKK)*cons(LLL)

        endif

        Ytemp=dabs(Y)
        If(dabs(Ytemp).gt.Ymax)Ymax=Ytemp

     enddo
  enddo

End subroutine classcutoff

subroutine classprim(I,J,K,L,II,JJ,KK,LL,NNA,NNC,NNAB,NNCD,Ymax1,IIIxiao,JJJxiao)
  use allmod

  Implicit real*8(A-H,O-Z)
  real*8 store(120,120)
  INTEGER NA(3),NB(3),NC(3),ND(3)
  real*8 P(3),Q(3),W(3),KAB,KCD
  Parameter(NN=13)
  real*8 FM(0:13)
  real*8 RA(3),RB(3),RC(3),RD(3)
  real*8 X44(1296)
  real*8 X4444(6,6)

  real*8 coefangxiaoL(20),coefangxiaoR(20)
  integer angxiaoL(20),angxiaoR(20),numangularL,numangularR

  COMMON /COM1/RA,RB,RC,RD
  COMMON /COM2/AA,BB,CC,DD,AB,CD,ROU,ABCD
  COMMON /COM4/P,Q,W
  COMMON /COM5/FM

  ! integer II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2
  common /xiaostore/store
  ! common /hrrstore/II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2

  ! ITT=0
  !

  ITT=0
  do JJJ=1,quick_basis%kprim(JJ)
     Nprij=quick_basis%kstart(JJ)+JJJ-1
     do III=1,quick_basis%kprim(II)
        Nprii=quick_basis%kstart(II)+III-1
        X2=X0*XCoeff(Nprii,Nprij,I,J)
        X4444(III,JJJ)=X2
     enddo
  enddo

  do MM2=NNC,NNCD
     do MM1=NNA,NNAB
        Ytemp=Yxiaoprim(IIIxiao,JJJxiao,MM1,MM2)*X4444(IIIxiao,JJJxiao)**2.0d0
        store(MM1,MM2)=Ytemp
        !                        print*,'***',X4,Yxiao(MM1,MM2,ITT),store(MM1,MM2)
     enddo
  enddo

  NBI1=quick_basis%Qsbasis(II,I)
  NBI2=quick_basis%Qfbasis(II,I)
  NBJ1=quick_basis%Qsbasis(JJ,J)
  NBJ2=quick_basis%Qfbasis(JJ,J)
  NBK1=quick_basis%Qsbasis(KK,K)
  NBK2=quick_basis%Qfbasis(KK,K)
  NBL1=quick_basis%Qsbasis(LL,L)
  NBL2=quick_basis%Qfbasis(LL,L)


  do III=quick_basis%ksumtype(II)+NBI1,quick_basis%ksumtype(II)+NBI2
     do JJJ=quick_basis%ksumtype(JJ)+NBJ1,quick_basis%ksumtype(JJ)+NBJ2
        !            do KKK=quick_basis%ksumtype(KK)+NBK1,quick_basis%ksumtype(KK)+NBK2
        !              do LLL=quick_basis%ksumtype(LL)+NBL1,quick_basis%ksumtype(LL)+NBL2
        KKK=III
        LLL=JJJ

        If((I.eq.0.and.J.eq.0.and.K.eq.0.and.L.eq.0).or. &
             (I.eq.1.and.J.eq.0.and.K.eq.1.and.L.eq.0))then

           do M=1,3
              NA(M)=KLMN(M,III)
              !                     NB(M)=KLMN(M,JJJ)
              NC(M)=KLMN(M,KKK)
              !                     ND(M)=KLMN(M,LLL) 
           enddo

           M1=trans(NA(1),NA(2),NA(3))
           M3=trans(NC(1),NC(2),NC(3))
           Y=store(M1,M3)

        elseif(I.eq.0.and.J.eq.1.and.K.eq.0.and.L.eq.1)then

           do M=1,3
              !                     NA(M)=KLMN(M,III)
              NB(M)=KLMN(M,JJJ)
              !                     NC(M)=KLMN(M,KKK)
              ND(M)=KLMN(M,LLL)
           enddo

           !                    Y=HRR(NA,NB,NC,ND)
           !                    Y=Y*cons(III)*cons(JJJ)*cons(KKK)*cons(LLL)

           M1=trans(NB(1),NB(2),NB(3))
           M3=trans(ND(1),ND(2),ND(3))

           do itemp=1,3
              if(ND(itemp).ne.0)then
                 ctemp=(RC(itemp)-RD(itemp))
                 Y1=store(M1,M3)+ctemp*store(M1,1)
                 Y2=store(1,M3)+ctemp*store(1,1)
                 goto 117
              endif
           enddo
117        continue

           do jtemp=1,3
              if(NB(jtemp).ne.0)then
                 Y=Y1+(RA(jtemp)-RB(jtemp))*Y2
                 goto 118
              endif
           enddo
118        continue

        elseif(I.eq.1.and.J.eq.1.and.K.eq.1.and.L.eq.1)then

           do M=1,3
              NA(M)=KLMN(M,III)
              NB(M)=KLMN(M,JJJ)
              NC(M)=KLMN(M,KKK)
              ND(M)=KLMN(M,LLL)
           enddo

           !                    Y=HRR(NA,NB,NC,ND)
           !                    Y=Y*cons(III)*cons(JJJ)*cons(KKK)*cons(LLL)

           MA=trans(NA(1),NA(2),NA(3))
           MAB=trans(NA(1)+NB(1),NA(2)+NB(2),NA(3)+NB(3))
           MCX=trans(NC(1),NC(2),NC(3))
           MCD=trans(NC(1)+ND(1),NC(2)+ND(2),NC(3)+ND(3))

           do itemp=1,3
              if(ND(itemp).ne.0)then
                 ctemp=(RC(itemp)-RD(itemp))
                 Y1=store(MAB,MCD)+ctemp*store(MAB,MCX)
                 Y2=store(MA,MCD)+ctemp*store(MA,MCX)
                 goto 141
              endif
           enddo
141        continue

           do jtemp=1,3
              if(NB(jtemp).ne.0)then
                 Y=Y1+(RA(jtemp)-RB(jtemp))*Y2
                 goto 142
              endif
           enddo
142        continue

        elseif((I.eq.2.and.J.eq.0.and.K.eq.2.and.L.eq.0).or. &
             (I.eq.0.and.J.eq.2.and.K.eq.0.and.L.eq.2).or. &
             (I.eq.2.and.J.eq.1.and.K.eq.2.and.L.eq.1).or. &
             (I.eq.1.and.J.eq.2.and.K.eq.1.and.L.eq.2).or. &
             (I.eq.2.and.J.eq.2.and.K.eq.2.and.L.eq.2).or. &
             (I.eq.3.and.J.eq.0.and.K.eq.3.and.L.eq.0).or. &
             (I.eq.0.and.J.eq.3.and.K.eq.0.and.L.eq.3).or. &
             (I.eq.3.and.J.eq.1.and.K.eq.3.and.L.eq.1).or. &
             (I.eq.1.and.J.eq.3.and.K.eq.1.and.L.eq.3).or. &
             (I.eq.3.and.J.eq.2.and.K.eq.3.and.L.eq.2).or. &
             (I.eq.2.and.J.eq.3.and.K.eq.2.and.L.eq.3).or. &
             (I.eq.3.and.J.eq.3.and.K.eq.3.and.L.eq.3))then

           IJtype=10*I+J
           KLtype=10*K+L

           call lefthrr(RA,RB,KLMN(1:3,III),KLMN(1:3,JJJ),IJtype,coefangxiaoL,angxiaoL,numangularL)
           call lefthrr(RC,RD,KLMN(1:3,KKK),KLMN(1:3,LLL),KLtype,coefangxiaoR,angxiaoR,numangularR)

           Y=0.0d0
           do ixiao=1,numangularL
              do jxiao=1,numangularR
                 Y=Y+coefangxiaoL(ixiao)*coefangxiaoR(jxiao)* &
                      store(angxiaoL(ixiao),angxiaoR(jxiao))
              enddo
           enddo

           Y=Y*cons(III)*cons(JJJ)*cons(KKK)*cons(LLL)

        endif

        Ytemp=dabs(Y)
        If(dabs(Ytemp).gt.Ymax1)Ymax1=Ytemp

     enddo
  enddo

End subroutine classprim

subroutine DNscreen(II,JJ,DNmax1)
  use allmod

  Implicit real*8(a-h,o-z)

  NII1=quick_basis%Qstart(II)
  NII2=quick_basis%Qfinal(II)
  NJJ1=quick_basis%Qstart(JJ)
  NJJ2=quick_basis%Qfinal(JJ)

  NBI1=quick_basis%Qsbasis(II,NII1)
  NBI2=quick_basis%Qfbasis(II,NII2)
  NBJ1=quick_basis%Qsbasis(JJ,NJJ1)
  NBJ2=quick_basis%Qfbasis(JJ,NJJ2)

  II111=quick_basis%ksumtype(II)+NBI1
  II112=quick_basis%ksumtype(II)+NBI2
  JJ111=quick_basis%ksumtype(JJ)+NBJ1
  JJ112=quick_basis%ksumtype(JJ)+NBJ2

  do III=II111,II112
     !         do JJJ=max(III,JJ111),JJ112
     do JJJ=JJ111,JJ112
        DENSEJI=dabs(quick_qm_struct%dense(JJJ,III))
        If(DENSEJI.gt.DNmax1)DNmax1=DENSEJI
     enddo
  enddo

End subroutine DNscreen

