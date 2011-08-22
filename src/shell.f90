subroutine g2eshell
  !*******************************************************
  ! g2eshell
  !--------------------------------------------------------
  ! BE careful of FmT.f90!!!!!!!!!!!!
  ! BE careful of IJKLtype!!!!!!!!!!!!!!!!!
  ! BE careful of size FM,store
  ! Xiao HE Prepare pair-shell constant 07/07/07 version
  !--------------------------------------------------------
  ! This subroutine is to Use the shell structure as initial guess
  ! to save the computational time
  !--------------------------------------------------------
  !
  use allmod
  !
  implicit double precision(a-h,o-z)
  real*8 Kab 
  double precision temp(3)     

  call zeroVec(temp,3)

  do ics=1,jshell                       ! ics is the shell no.
     do ips=1,quick_basis%kprim(ics)                ! ips is prim no. for certain shell
        Nprii=quick_basis%kstart(ics)+IPS-1         ! we can get globle prim no. for ips
        AA=gcexpo(ips,quick_basis%ksumtype(ics))    ! so we have the exponent part for ics shell ips prim
        do jcs=1,jshell
           do jps=1,quick_basis%kprim(jcs)
              Nprij=quick_basis%kstart(jcs)+jps-1
              BB=gcexpo(jps,quick_basis%ksumtype(jcs))
              Apri(Nprii,Nprij)=AA+BB
              do I=1,3
                 Ppri(i,Nprii,Nprij)=(xyz(i,quick_basis%katom(ics))*AA+ &
                      xyz(i,quick_basis%katom(jcs))*BB)/(AA+BB)
              enddo
              ! K(A,B)=1/(a+b)*exp(-ab/(a+b)*(A-B)^2)
              Kpri(Nprii,Nprij)=dexp(-AA*BB/(AA+BB)* &
                (quick_molspec%atomdistance(quick_basis%katom(ics),quick_basis%katom(jcs)))**2)/(AA+BB)
              X1temp=ONE*Kpri(Nprii,Nprij)
              do  itemp=quick_basis%Qstart(ics),quick_basis%Qfinal(ics)
                 X1=X1temp*gccoeff(ips,quick_basis%ksumtype(ics)+Itemp)
                 do itemp2=quick_basis%Qstart(jcs),quick_basis%Qfinal(jcs)
                    Xcoeff(Nprii,Nprij,itemp,itemp2)=X1*gccoeff(jps,quick_basis%ksumtype(jcs)+Itemp2)
                 enddo
              enddo
           enddo
        enddo
     enddo
  enddo
end subroutine g2eshell

! Vertical Recursion by Xiao HE 07/07/07 version
subroutine shell
  use allmod

  Implicit real*8(a-h,o-z)
  real*8 P(3),Q(3),W(3),KAB,KCD,AAtemp(3)
  Parameter(NN=13)
  real*8 FM(0:13)
  real*8 RA(3),RB(3),RC(3),RD(3)

  real*8 Qtemp(3),WQtemp(3),CDtemp,ABcom,Ptemp(3),WPtemp(3),ABtemp,CDcom,ABCDtemp
  integer II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2
  common /hrrstore/II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2
  COMMON /VRRcom/Qtemp,WQtemp,CDtemp,ABcom,Ptemp,WPtemp,ABtemp,CDcom,ABCDtemp
  COMMON /COM1/RA,RB,RC,RD


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


  do JJJ=1,quick_basis%kprim(JJ)
     Nprij=quick_basis%kstart(JJ)+JJJ-1
     do III=1,quick_basis%kprim(II)
        Nprii=quick_basis%kstart(II)+III-1
        AB=Apri(Nprii,Nprij)
        ABtemp=0.5d0/AB
        cutoffprim1=dnmax*cutprim(Nprii,Nprij)

        do M=1,3
           P(M)=Ppri(M,Nprii,Nprij)
           AAtemp(M)=P(M)*AB
           Ptemp(M)=P(M)-RA(M)
        enddo

        do LLL=1,quick_basis%kprim(LL)
           Npril=quick_basis%kstart(LL)+LLL-1
           do KKK=1,quick_basis%kprim(KK)
              Nprik=quick_basis%kstart(KK)+KKK-1
              cutoffprim=cutoffprim1*cutprim(Nprik,Npril)
              If(cutoffprim.gt.quick_method%primLimit)then
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
                    W(M)=(AAtemp(M)+Q(M)*CD)/ABCD
                    XXXtemp=P(M)-Q(M)
                    RPQ=RPQ+XXXtemp*XXXtemp
                    Qtemp(M)=Q(M)-RC(M)
                    WQtemp(M)=W(M)-Q(M)
                    WPtemp(M)=W(M)-P(M)
                 enddo
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
              endif
           enddo
        enddo
     enddo
  enddo


  do I=NII1,NII2
     NNA=Sumindex(I-1)+1
     do J=NJJ1,NJJ2
        NNAB=SumINDEX(I+J)
        do K=NKK1,NKK2
           NNC=Sumindex(k-1)+1
           do L=NLL1,NLL2
              NNCD=SumIndex(K+L)
              call iclass(I,J,K,L,NNA,NNC,NNAB,NNCD)
           enddo
        enddo
     enddo
  enddo
201 return
end subroutine shell

! Horrizontal recursion and Fock matrix builder by Xiao HE 07/07/07 version
subroutine iclass(I,J,K,L,NNA,NNC,NNAB,NNCD)
  use allmod

  Implicit real*8(A-H,O-Z)
  real*8 store(120,120)
  INTEGER NA(3),NB(3),NC(3),ND(3)
  real*8 P(3),Q(3),W(3),KAB,KCD
  Parameter(NN=13)
  real*8 FM(0:13)
  real*8 RA(3),RB(3),RC(3),RD(3)
  real*8 X44(129600)

  COMMON /COM1/RA,RB,RC,RD
  COMMON /COM2/AA,BB,CC,DD,AB,CD,ROU,ABCD
  COMMON /COM4/P,Q,W
  COMMON /COM5/FM

  integer II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2
  common /xiaostore/store
  common /hrrstore/II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2

  ITT=0
  do JJJ=1,quick_basis%kprim(JJ)
     Nprij=quick_basis%kstart(JJ)+JJJ-1
     do III=1,quick_basis%kprim(II)
        Nprii=quick_basis%kstart(II)+III-1
        X2=X0*XCoeff(Nprii,Nprij,I,J)
        cutoffprim1=dnmax*cutprim(Nprii,Nprij)
        do LLL=1,quick_basis%kprim(LL)
           Npril=quick_basis%kstart(LL)+LLL-1
           do KKK=1,quick_basis%kprim(KK)
              Nprik=quick_basis%kstart(KK)+KKK-1
              cutoffprim=cutoffprim1*cutprim(Nprik,Npril)
              If(cutoffprim.gt.quick_method%primLimit)then
                 ITT=ITT+1
                 X44(ITT)=X2*XCoeff(Nprik,Npril,K,L)
              endif
           enddo
        enddo
     enddo
  enddo

  do MM2=NNC,NNCD
     do MM1=NNA,NNAB
        Ytemp=0.0d0
        do itemp=1,ITT
           Ytemp=Ytemp+X44(itemp)*Yxiao(itemp,MM1,MM2)
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


  IJtype=10*I+J
  KLtype=10*K+L
  IJKLtype=100*IJtype+KLtype


  if((max(I,J,K,L).eq.2.and.(J.ne.0.or.L.ne.0)).or.(max(I,J,K,L).ge.3))IJKLtype=999

  III1=quick_basis%ksumtype(II)+NBI1
  III2=quick_basis%ksumtype(II)+NBI2
  JJJ1=quick_basis%ksumtype(JJ)+NBJ1
  JJJ2=quick_basis%ksumtype(JJ)+NBJ2
  KKK1=quick_basis%ksumtype(KK)+NBK1
  KKK2=quick_basis%ksumtype(KK)+NBK2
  LLL1=quick_basis%ksumtype(LL)+NBL1
  LLL2=quick_basis%ksumtype(LL)+NBL2

  If(II.lt.JJ.and.II.lt.KK.and.KK.lt.LL)then
     do III=III1,III2
        do JJJ=JJJ1,JJJ2
           do KKK=KKK1,KKK2
              do LLL=LLL1,LLL2
                 call hrrwhole
                 DENSEKI=quick_qm_struct%dense(KKK,III)
                 DENSEKJ=quick_qm_struct%dense(KKK,JJJ)
                 DENSELJ=quick_qm_struct%dense(LLL,JJJ)
                 DENSELI=quick_qm_struct%dense(LLL,III)
                 DENSELK=quick_qm_struct%dense(LLL,KKK)
                 DENSEJI=quick_qm_struct%dense(JJJ,III)

                 ! Find the (ij|kl) integrals where j>i,k>i,l>k. Note that k and j
                 ! can be equal.
                 quick_qm_struct%o(JJJ,III) = quick_qm_struct%o(JJJ,III)+2.d0*DENSELK*Y
                 quick_qm_struct%o(LLL,KKK) = quick_qm_struct%o(LLL,KKK)+2.d0*DENSEJI*Y
                 quick_qm_struct%o(KKK,III) = quick_qm_struct%o(KKK,III)-.5d0*DENSELJ*Y
                 quick_qm_struct%o(LLL,III) = quick_qm_struct%o(LLL,III)-.5d0*DENSEKJ*Y
                 quick_qm_struct%o(JJJ,KKK) = quick_qm_struct%o(JJJ,KKK)-.5d0*DENSELI*Y
                 quick_qm_struct%o(JJJ,LLL) = quick_qm_struct%o(JJJ,LLL)-.5d0*DENSEKI*Y
                 quick_qm_struct%o(KKK,JJJ) = quick_qm_struct%o(KKK,JJJ)-.5d0*DENSELI*Y
                 quick_qm_struct%o(LLL,JJJ) = quick_qm_struct%o(LLL,JJJ)-.5d0*DENSEKI*Y

              enddo
           enddo
        enddo
     enddo
  else
     do III=III1,III2                          
        do JJJ=max(III,JJJ1),JJJ2                          
           do KKK=max(III,KKK1),KKK2                          
              do LLL=max(KKK,LLL1),LLL2                          
                 If(III.LT.KKK)then
                    call hrrwhole
                    If(III.lt.JJJ.and.KKK.lt.LLL)then
                       DENSEKI=quick_qm_struct%dense(KKK,III)
                       DENSEKJ=quick_qm_struct%dense(KKK,JJJ)
                       DENSELJ=quick_qm_struct%dense(LLL,JJJ)
                       DENSELI=quick_qm_struct%dense(LLL,III)
                       DENSELK=quick_qm_struct%dense(LLL,KKK)
                       DENSEJI=quick_qm_struct%dense(JJJ,III)

                       ! Find the (ij|kl) integrals where j>i,k>i,l>k. Note that k and j
                       ! can be equal.

                       quick_qm_struct%o(JJJ,III) = quick_qm_struct%o(JJJ,III)+2.d0*DENSELK*Y
                       quick_qm_struct%o(LLL,KKK) = quick_qm_struct%o(LLL,KKK)+2.d0*DENSEJI*Y
                       quick_qm_struct%o(KKK,III) = quick_qm_struct%o(KKK,III)-.5d0*DENSELJ*Y
                       quick_qm_struct%o(LLL,III) = quick_qm_struct%o(LLL,III)-.5d0*DENSEKJ*Y
                       quick_qm_struct%o(JJJ,KKK) = quick_qm_struct%o(JJJ,KKK)-.5d0*DENSELI*Y
                       quick_qm_struct%o(JJJ,LLL) = quick_qm_struct%o(JJJ,LLL)-.5d0*DENSEKI*Y
                       quick_qm_struct%o(KKK,JJJ) = quick_qm_struct%o(KKK,JJJ)-.5d0*DENSELI*Y
                       quick_qm_struct%o(LLL,JJJ) = quick_qm_struct%o(LLL,JJJ)-.5d0*DENSEKI*Y
                    elseif(III.eq.JJJ.and.KKK.eq.LLL)then
                       DENSEJI=quick_qm_struct%dense(KKK,III)
                       DENSEJJ=quick_qm_struct%dense(KKK,KKK)
                       DENSEII=quick_qm_struct%dense(III,III)

                       ! Find  all the (ii|jj) integrals.
                       quick_qm_struct%o(III,III) = quick_qm_struct%o(III,III)+DENSEJJ*Y
                       quick_qm_struct%o(KKK,KKK) = quick_qm_struct%o(KKK,KKK)+DENSEII*Y
                       quick_qm_struct%o(KKK,III) = quick_qm_struct%o(KKK,III)-.5d0*DENSEJI*Y

                    elseif(JJJ.eq.KKK.and.JJJ.eq.LLL)then
                       DENSEJI=quick_qm_struct%dense(JJJ,III)
                       DENSEJJ=quick_qm_struct%dense(JJJ,JJJ)

                       ! Find  all the (ij|jj) integrals.
                       quick_qm_struct%o(JJJ,III) = quick_qm_struct%o(JJJ,III)+.5d0*DENSEJJ*Y
                       quick_qm_struct%o(JJJ,JJJ) = quick_qm_struct%o(JJJ,JJJ)+DENSEJI*Y

                       !        ! Find  all the (ii|ij) integrals.
                       !        ! Find all the (ij|ij) integrals

                       ! Find all the (ij|ik) integrals where j>i,k>j
                    elseif(KKK.eq.LLL.and.III.lt.JJJ.and.JJJ.ne.KKK)then
                       DENSEKI=quick_qm_struct%dense(KKK,III)
                       DENSEKJ=quick_qm_struct%dense(KKK,JJJ)
                       DENSEKK=quick_qm_struct%dense(KKK,KKK)
                       DENSEJI=quick_qm_struct%dense(JJJ,III)

                       ! Find all the (ij|kk) integrals where j>i, k>j.
                       quick_qm_struct%o(JJJ,III) = quick_qm_struct%o(JJJ,III)+DENSEKK*Y
                       quick_qm_struct%o(KKK,KKK) = quick_qm_struct%o(KKK,KKK)+2.d0*DENSEJI*Y
                       quick_qm_struct%o(KKK,III) = quick_qm_struct%o(KKK,III)-.5d0*DENSEKJ*Y
                       quick_qm_struct%o(KKK,JJJ) = quick_qm_struct%o(KKK,JJJ)-.5d0*DENSEKI*Y
                       quick_qm_struct%o(JJJ,KKK) = quick_qm_struct%o(JJJ,KKK)-.5d0*DENSEKI*Y

                       !        ! Find all the (ik|jj) integrals where j>i, k>j.
                    elseif(III.eq.JJJ.and.KKK.lt.LLL)then
                       DENSEII=quick_qm_struct%dense(III,III)
                       DENSEJI=quick_qm_struct%dense(KKK,III)
                       DENSEKI=quick_qm_struct%dense(LLL,III)
                       DENSEKJ=quick_qm_struct%dense(LLL,KKK)

                       ! Find all the (ii|jk) integrals where j>i, k>j.
                       quick_qm_struct%o(LLL,KKK) = quick_qm_struct%o(LLL,KKK)+DENSEII*Y
                       quick_qm_struct%o(III,III) = quick_qm_struct%o(III,III)+2.d0*DENSEKJ*Y
                       quick_qm_struct%o(KKK,III) = quick_qm_struct%o(KKK,III)-.5d0*DENSEKI*Y
                       quick_qm_struct%o(LLL,III) = quick_qm_struct%o(LLL,III)-.5d0*DENSEJI*Y

                    endif

                 else
                    If(JJJ.LE.LLL)then
                       call hrrwhole
                       if(III.eq.JJJ.and.III.eq.KKK.and.III.eq.LLL)then
                          DENSEII=quick_qm_struct%dense(III,III)

                          ! do all the (ii|ii) integrals.
                          quick_qm_struct%o(III,III) = quick_qm_struct%o(III,III)+.5d0*DENSEII*Y
                       elseif(III.eq.JJJ.and.III.eq.KKK.and.III.lt.LLL)then
                          DENSEJI=quick_qm_struct%dense(LLL,III)
                          DENSEII=quick_qm_struct%dense(III,III)

                          ! Find  all the (ii|ij) integrals.
                          quick_qm_struct%o(LLL,III) = quick_qm_struct%o(LLL,III)+.5d0*DENSEII*Y
                          quick_qm_struct%o(III,III) = quick_qm_struct%o(III,III)+DENSEJI*Y

                       elseif(III.eq.KKK.and.JJJ.eq.LLL.and.III.lt.JJJ)then
                          DENSEJI=quick_qm_struct%dense(JJJ,III)
                          DENSEJJ=quick_qm_struct%dense(JJJ,JJJ)
                          DENSEII=quick_qm_struct%dense(III,III)

                          ! Find all the (ij|ij) integrals
                          quick_qm_struct%o(JJJ,III) = quick_qm_struct%o(JJJ,III)+1.50*DENSEJI*Y
                          quick_qm_struct%o(JJJ,JJJ) = quick_qm_struct%o(JJJ,JJJ)-.5d0*DENSEII*Y
                          quick_qm_struct%o(III,III) = quick_qm_struct%o(III,III)-.5d0*DENSEJJ*Y

                       elseif(III.eq.KKK.and.III.lt.JJJ.and.JJJ.lt.LLL)then
                          DENSEKI=quick_qm_struct%dense(LLL,III)
                          DENSEKJ=quick_qm_struct%dense(LLL,JJJ)
                          DENSEII=quick_qm_struct%dense(III,III)
                          DENSEJI=quick_qm_struct%dense(JJJ,III)

                          ! Find all the (ij|ik) integrals where j>i,k>j
                          quick_qm_struct%o(JJJ,III) = quick_qm_struct%o(JJJ,III)+1.5d0*DENSEKI*Y
                          quick_qm_struct%o(LLL,III) = quick_qm_struct%o(LLL,III)+1.5d0*DENSEJI*Y
                          quick_qm_struct%o(III,III) = quick_qm_struct%o(III,III)-1.d0*DENSEKJ*Y
                          quick_qm_struct%o(LLL,JJJ) = quick_qm_struct%o(LLL,JJJ)-.5d0*DENSEII*Y
                       endif
                    endif
                 endif
              enddo
           enddo
        enddo
     enddo
  endif

201 return    
End subroutine iclass


subroutine vertical(NABCDTYPE)
  use allMod

  implicit double precision(a-h,o-z)

  select case (NABCDTYPE)

  case(0)
     !                              goto 123
  case(10)
     call PSSS(0)
  case(1)
     call SSPS(0)
  case(11)
     call PSSS(0)

     call SSPS(0)
     call SSPS(1)
     call PSPS(0)
  case(20)
     call PSSS(0)
     call PSSS(1)
     call DSSS(0)
  case(2)
     call SSPS(0)
     call SSPS(1)
     call SSDS(0)
  case(21)
     call SSPS(0)

     call PSSS(0)
     call SSPS(1)
     call PSPS(0)

     call PSSS(1)
     call DSSS(0)

     call PSSS(2)
     call DSSS(1)
     call DSPS(0)
  case(12)
     call SSPS(0)

     call PSSS(0)
     call SSPS(1)
     call PSPS(0)

     call SSDS(0)

     call SSPS(2) 
     call SSDS(1)
     call PSDS(0)
  case(22)
     call SSPS(0)

     call PSSS(0)
     call SSPS(1)
     call PSPS(0)

     call SSDS(0)
     call SSPS(2) 
     call SSDS(1)
     call PSDS(0)

     call PSSS(1)
     call DSSS(0)
     call PSSS(2)
     call DSSS(1)
     call DSPS(0)

     !                             call SSPS(2)
     call SSPS(3)
     call SSDS(2)
     call PSDS(1)

     call PSPS(1)
     call DSDS(0)

  case(30)
     call PSSS(0)
     call PSSS(1)
     call DSSS(0)

     call PSSS(2)
     call DSSS(1)

     call FSSS(0)

  case(3)
     call SSPS(0)
     call SSPS(1)
     call SSDS(0)

     call SSPS(2)
     call SSDS(1)

     call SSFS(0)

  case(40)
     call PSSS(0)
     call PSSS(1)
     call DSSS(0)

     call PSSS(2)
     call DSSS(1)

     call FSSS(0)

     call PSSS(3)
     call DSSS(2)

     call FSSS(1)

     call GSSS(0)

  case(4)
     call SSPS(0)
     call SSPS(1)
     call SSDS(0)

     call SSPS(2)
     call SSDS(1)

     call SSFS(0)

     call SSPS(3)
     call SSDS(2)

     call SSFS(1)

     call SSGS(0)

  case(31)
     call PSSS(0)
     call PSSS(1)
     call DSSS(0)

     call PSSS(2)
     call DSSS(1)

     call FSSS(0)

     call PSSS(3)
     call DSSS(2)

     call FSSS(1)

     call SSPS(0)
     call SSPS(1)
     call PSPS(0)

     call DSPS(0)

     call FSPS(0)

  case(13)
     call SSPS(0)
     call SSPS(1)
     call SSDS(0)

     call SSPS(2)
     call SSDS(1)

     call SSFS(0)

     call SSPS(3)
     call SSDS(2)

     call SSFS(1)

     call PSSS(0)
     !                             call PSSS(1)
     call PSPS(0)

     call PSDS(0)

     call PSFS(0)

  case(41)
     call PSSS(0)
     call PSSS(1)
     call DSSS(0)

     call PSSS(2)
     call DSSS(1)

     call FSSS(0)

     call PSSS(3)
     call DSSS(2)

     call FSSS(1)

     call GSSS(0)

     call PSSS(4)
     call DSSS(3)

     call FSSS(2)

     call GSSS(1)

     call SSPS(0)
     call SSPS(1)
     call PSPS(0)

     !                             call PSSS(1)
     call DSPS(0)

     call FSPS(0)

     call GSPS(0)

  case(14)
     call SSPS(0)
     call SSPS(1)
     call SSDS(0)

     call SSPS(2)
     call SSDS(1)

     call SSFS(0)

     call SSPS(3)
     call SSDS(2)

     call SSFS(1)

     call PSSS(0)
     call PSSS(1)
     call PSPS(0)

     call SSGS(0)

     call SSPS(4)
     call SSDS(3)

     call SSFS(2)

     call SSGS(1)

     !                             call PSSS(1)
     call PSDS(0)

     call PSFS(0)

     call PSGS(0)

  case(32)
     call PSSS(0)
     call PSSS(1)
     call DSSS(0)

     call PSSS(2)
     call DSSS(1)

     call FSSS(0)

     call PSSS(3)
     call DSSS(2)

     call FSSS(1)

     call SSPS(0)
     call SSPS(1)
     call PSPS(0)

     call DSPS(0)

     call FSPS(0)

     call PSSS(4)
     call DSSS(3)

     call FSSS(2)

     call FSPS(1)

     call DSPS(1)

     call SSDS(0)

     !                             call SSDS(0)

     call SSPS(2)
     call SSDS(1)
     call PSDS(0)    

     call SSPS(3)
     call SSDS(2)
     call PSDS(1)

     !                             call SSPS(1)
     !                             call SSDS(1)

     !                             call PSSS(1)
     call PSPS(1)
     call DSDS(0)

     call FSDS(0)                                                      

  case(23)
     !                           if(NNAB.eq.2.and.NNCD.eq.3)then
     call SSPS(0)
     call SSPS(1)
     call SSDS(0)

     call SSPS(2)
     call SSDS(1)

     call SSFS(0)

     call SSPS(3)
     call SSDS(2)

     call SSFS(1)

     call PSSS(0)
     call PSSS(1)
     call PSPS(0)

     call PSDS(0)

     call PSFS(0)

     call SSPS(4)
     call SSDS(3)

     call SSFS(2)

     call PSFS(1)

     call PSDS(1)

     call DSSS(0)

     call PSSS(2)
     call DSSS(1)
     call DSPS(0)

     !                             call PSSS(3)
     !                             call DSSS(2)
     !                             call PSDS(1)

     !                             call SSPS(1)
     !                             call SSDS(1)

     !                             call PSSS(1)
     call PSPS(1)
     call DSDS(0)

     call DSFS(0)                                        

  case(42)
     !                           if(NNAB.eq.4.and.NNCD.eq.2)then
     call PSSS(0)
     call PSSS(1)
     call DSSS(0)

     call PSSS(2)
     call DSSS(1)

     call FSSS(0)

     call PSSS(3)
     call DSSS(2)

     call FSSS(1)

     call GSSS(0)

     call PSSS(4)
     call DSSS(3)

     call FSSS(2)

     call GSSS(1)

     !                             call PSSS(1)
     call DSPS(0)

     call FSPS(0)

     call GSPS(0)

     call PSSS(5)
     call DSSS(4)
     call FSSS(3)

     call GSSS(2)
     call DSPS(1)

     call FSPS(1)

     call GSPS(1)

     !********************DSDS

     call SSPS(0)

     !                             call PSSS(0)
     call SSPS(1)
     call PSPS(0)

     call SSDS(0)
     call SSPS(2)
     call SSDS(1)
     call PSDS(0)

     !                             call PSSS(1)
     !                             call DSSS(0)
     !                             call PSSS(2)
     !                             call DSSS(1)
     !                             call DSPS(0)

     call SSPS(3)
     call SSDS(2)
     call PSDS(1)

     call PSPS(1)
     call DSDS(0)

     !*******************FSDS
     !                             call FSPS(1)

     !                              call DSPS(1)

     !                             call SSDS(0)

     !                             call SSPS(2)
     !                             call SSDS(1)
     !                             call PSDS(0)

     !                             call SSPS(3)
     !                             call SSDS(2)
     !                             call PSDS(1)

     !                             call PSPS(1)
     !                             call DSDS(0)

     call FSDS(0)

     call GSDS(0)


  case(24)
     !                           if(NNAB.eq.2.and.NNCD.eq.4)then
     call SSPS(0)
     call SSPS(1)
     call SSDS(0)

     call SSPS(2)
     call SSDS(1)

     call SSFS(0)

     call SSPS(3)
     call SSDS(2)

     call SSFS(1)

     call SSGS(0)

     call SSPS(4)
     call SSDS(3)

     call SSFS(2)

     call SSGS(1)

     !                             call PSSS(1)
     call PSDS(0)

     call PSFS(0)

     call PSGS(0)

     call SSPS(5)
     call SSDS(4)

     call SSFS(3)

     call SSGS(2)
     call PSDS(1)
     call PSFS(1)

     call PSGS(1)

     !********************DSDS

     call PSSS(0)

     !                             call PSSS(0)
     call PSSS(1)
     call PSPS(0)

     call DSSS(0)
     call PSSS(2)
     call DSSS(1)
     call DSPS(0)

     !                             call PSSS(1)
     !                             call DSSS(0)
     !                             call PSSS(2)
     !                             call DSSS(1)
     !                             call DSPS(0)

     call PSSS(3)
     call DSSS(2)
     call DSPS(1)

     call PSPS(1)
     call DSDS(0)

     !*******************FSDS

     !                             call PSFS(1)

     !                             call PSDS(1)

     !                             call SSDS(0)

     !                             call SSPS(2)
     !                             call SSDS(1)
     !                             call PSDS(0)

     !                             call SSPS(3)
     !                             call SSDS(2)
     !                             call PSDS(1)

     !                             call PSPS(1)
     !                             call DSDS(0)

     call DSFS(0)
     call DSGS(0)


  case(33)
     !                           if(NNAB.eq.3.and.NNCD.eq.3)then
     do ixiao=0,5
        call PSSS(ixiao)
        !                              call SSPS(ixiao)
     enddo
     do ixiao=0,4
        call SSPS(ixiao)
     enddo
     do ixiao=0,3
        call PSPS(ixiao)
     enddo
     do ixiao=0,4
        call DSSS(ixiao)
        !                              call SSDS(ixiao)
     enddo
     do ixiao=0,3
        call SSDS(ixiao)
     enddo
     do ixiao=0,3
        call DSPS(ixiao)
        !                              call PSDS(ixiao)
     enddo
     do ixiao=0,2
        call PSDS(ixiao)
     enddo
     do ixiao=0,1
        call DSDS(ixiao)
     enddo
     do ixiao=0,3
        call FSSS(ixiao)
        !                              call SSFS(ixiao)
     enddo
     do ixiao=0,2
        call SSFS(ixiao)
     enddo
     do ixiao=0,2
        call FSPS(ixiao)
        !                              call PSFS(ixiao)
     enddo
     do ixiao=0,1
        call PSFS(ixiao)
     enddo
     call FSDS(0)
     call DSFS(0)
     call FSDS(1)
     !                              call DSFS(1)
     call FSFS(0)

  case(43)
     do ixiao=0,6
        call PSSS(ixiao)
        !                              call SSPS(ixiao)
     enddo
     do ixiao=0,5
        call SSPS(ixiao)
     enddo
     do ixiao=0,4
        call PSPS(ixiao)
     enddo
     do ixiao=0,5
        call DSSS(ixiao)
        !                             call SSDS(ixiao)
     enddo
     do ixiao=0,4
        call SSDS(ixiao)
     enddo
     do ixiao=0,4
        call DSPS(ixiao)
        !                              call PSDS(ixiao)
     enddo
     do ixiao=0,3
        call PSDS(ixiao)
     enddo
     do ixiao=0,2
        call DSDS(ixiao)
     enddo
     do ixiao=0,4
        call FSSS(ixiao)
        !                              call SSFS(ixiao)
     enddo
     do ixiao=0,3
        call SSFS(ixiao)
     enddo
     do ixiao=0,3
        call FSPS(ixiao)
        !                              call PSFS(ixiao)
     enddo
     do ixiao=0,2
        call PSFS(ixiao)
     enddo
     call FSDS(0)
     call DSFS(0)
     call FSDS(1)
     call DSFS(1)
     call FSDS(2)
     !                              call DSFS(2)
     call FSFS(0)
     !                              call FSFS(1)
     do ixiao=0,3
        call GSSS(ixiao)
     enddo
     do ixiao=0,2
        call GSPS(ixiao)
     enddo
     do ixiao=0,1
        call GSDS(ixiao)
     enddo
     call GSFS(0)

  case(34)
     !                           if(NNAB.eq.3.and.NNCD.eq.4)then
     do ixiao=0,6
        !                              call PSSS(ixiao)
        call SSPS(ixiao)
     enddo
     do ixiao=0,5
        call PSSS(ixiao)
     enddo
     do ixiao=0,4
        call PSPS(ixiao)
     enddo
     do ixiao=0,5
        !                              call DSSS(ixiao)
        call SSDS(ixiao)
     enddo
     do ixiao=0,4
        call DSSS(ixiao)
     enddo
     do ixiao=0,4
        !                              call DSPS(ixiao)
        call PSDS(ixiao)
     enddo
     do ixiao=0,3
        call DSPS(ixiao)
     enddo
     do ixiao=0,2
        call DSDS(ixiao)
     enddo
     do ixiao=0,4
        !                              call FSSS(ixiao)
        call SSFS(ixiao)
     enddo
     do ixiao=0,3
        call FSSS(ixiao)
     enddo
     do ixiao=0,3
        !                              call FSPS(ixiao)
        call PSFS(ixiao)
     enddo
     do ixiao=0,2
        call FSPS(ixiao)
     enddo
     call FSDS(0)
     call DSFS(0)
     call FSDS(1)
     call DSFS(1)
     !                              call FSDS(2)
     call DSFS(2)
     call FSFS(0)
     !                              call FSFS(1)
     do ixiao=0,3
        call SSGS(ixiao)
     enddo
     do ixiao=0,2
        call PSGS(ixiao)
     enddo
     do ixiao=0,1
        call DSGS(ixiao)
     enddo
     call FSGS(0)

  case(44)
     !                           if(NNAB.eq.4.and.NNCD.eq.4)then
     do ixiao=0,7
        call PSSS(ixiao)
        !                              call SSPS(ixiao)
     enddo
     do ixiao=0,6
        call SSPS(ixiao)
     enddo
     do ixiao=0,5
        call PSPS(ixiao)
     enddo
     do ixiao=0,6
        call DSSS(ixiao)
        !                              call SSDS(ixiao)
     enddo
     do ixiao=0,5
        call SSDS(ixiao)
     enddo
     do ixiao=0,5
        call DSPS(ixiao)
        !                              call PSDS(ixiao)
     enddo
     do ixiao=0,4
        call PSDS(ixiao)
     enddo
     do ixiao=0,3
        call DSDS(ixiao)
     enddo
     do ixiao=0,5
        call FSSS(ixiao)
        !                              call SSFS(ixiao)
     enddo
     do ixiao=0,4
        call SSFS(ixiao)
     enddo
     do ixiao=0,4
        call FSPS(ixiao)
        !                              call PSFS(ixiao)
     enddo
     do ixiao=0,3
        call PSFS(ixiao)
     enddo
     do ixiao=0,3
        call FSDS(ixiao)
        !                              call DSFS(ixiao)
     enddo
     do ixiao=0,2
        call DSFS(ixiao)
     enddo
     do ixiao=0,1 
        call FSFS(ixiao)
     enddo
     do ixiao=0,4
        call GSSS(ixiao)
        !                              call SSGS(ixiao)
     enddo
     do ixiao=0,3
        call SSGS(ixiao)
     enddo
     do ixiao=0,3
        call GSPS(ixiao)
        !                              call PSGS(ixiao)
     enddo
     do ixiao=0,2
        call PSGS(ixiao)
     enddo
     do ixiao=0,2
        call GSDS(ixiao)
        !                              call DSGS(ixiao)
     enddo
     do ixiao=0,1
        call DSGS(ixiao)
     enddo
     do ixiao=0,1
        call GSFS(ixiao)
        !                              call FSGS(ixiao)
     enddo
     call FSGS(0)
     call GSGS(0)

  case(50)
     do ixiao=0,4
        call PSSS(ixiao)
     enddo
     do ixiao=0,3
        call DSSS(ixiao)
     enddo
     do ixiao=0,2
        call FSSS(ixiao)
     enddo
     do ixiao=0,1
        call GSSS(ixiao)
     enddo
     call BSLS(5,0,0)

  case(5)
     do ixiao=0,4
        call SSPS(ixiao)
     enddo
     do ixiao=0,3
        call SSDS(ixiao)
     enddo
     do ixiao=0,2
        call SSFS(ixiao)
     enddo
     do ixiao=0,1
        call SSGS(ixiao)
     enddo
     call LSBS(0,5,0)

  case(51)
     do ixiao=0,5
        call PSSS(ixiao)
        call SSPS(ixiao)
     enddo
     do ixiao=0,4
        call PSPS(ixiao)
     enddo
     do ixiao=0,4
        call DSSS(ixiao)
     enddo
     do ixiao=0,3
        call DSPS(ixiao)
        call FSSS(ixiao)
     enddo
     do ixiao=0,2
        call FSPS(ixiao)
        call GSSS(ixiao)
     enddo
     do ixiao=0,1
        call GSPS(ixiao)
        call BSLS(5,0,ixiao)
     enddo
     call BSLS(5,1,0)

  case(15)
     do ixiao=0,5
        call SSPS(ixiao)
        call PSSS(ixiao)
     enddo
     do ixiao=0,4
        call PSPS(ixiao)
        call SSDS(ixiao)
     enddo
     do ixiao=0,3
        call PSDS(ixiao)
        call SSFS(ixiao)
     enddo
     do ixiao=0,2
        call PSFS(ixiao)
        call SSGS(ixiao)
     enddo
     do ixiao=0,1
        call PSGS(ixiao)
        call LSBS(0,5,ixiao)
     enddo
     call LSBS(1,5,0)

  case(52)
     do ixiao=0,6
        call PSSS(ixiao)
        call SSPS(ixiao)
     enddo
     do ixiao=0,5
        call PSPS(ixiao)
        call DSSS(ixiao)
        call SSDS(ixiao)
     enddo
     do ixiao=0,4
        call PSDS(ixiao)
        call DSPS(ixiao)
        call FSSS(ixiao)
     enddo
     do ixiao=0,3
        call DSDS(ixiao)
        call FSPS(ixiao)
        call GSSS(ixiao)
     enddo
     do ixiao=0,2
        call FSDS(ixiao)
        call GSPS(ixiao)
        call BSLS(5,0,ixiao)
     enddo
     do ixiao=0,1
        call GSDS(ixiao) 
        call BSLS(5,1,ixiao)
     enddo
     call BSLS(5,2,0)

  case(25)
     do ixiao=0,6
        call SSPS(ixiao)
        call PSSS(ixiao)
     enddo
     do ixiao=0,5
        call PSPS(ixiao)
        call SSDS(ixiao)
        call DSSS(ixiao)
     enddo
     do ixiao=0,4
        call DSPS(ixiao)
        call PSDS(ixiao)
        call SSFS(ixiao)
     enddo
     do ixiao=0,3
        call DSDS(ixiao)
        call PSFS(ixiao)
        call SSGS(ixiao)
     enddo
     do ixiao=0,2
        call DSFS(ixiao)
        call PSGS(ixiao)
        call LSBS(0,5,ixiao)
     enddo
     do ixiao=0,1
        call DSGS(ixiao)
        call LSBS(1,5,ixiao)
     enddo
     call LSBS(2,5,0)

  case(53)
     do ixiao=0,7
        call PSSS(ixiao)
        call SSPS(ixiao)
     enddo
     do ixiao=0,6
        call PSPS(ixiao)
        call DSSS(ixiao)
        call SSDS(ixiao)
     enddo
     do ixiao=0,5
        call PSDS(ixiao)
        call DSPS(ixiao)
        call FSSS(ixiao)
        call SSFS(ixiao)
     enddo
     do ixiao=0,4
        call PSFS(ixiao)
        call DSDS(ixiao)
        call FSPS(ixiao)
        call GSSS(ixiao)
     enddo
     do ixiao=0,3
        call DSFS(ixiao)
        call FSDS(ixiao)
        call GSPS(ixiao)
        call BSLS(5,0,ixiao)
     enddo
     do ixiao=0,2
        call FSFS(ixiao)
        call GSDS(ixiao)
        call BSLS(5,1,ixiao)
     enddo
     do ixiao=0,1
        call GSFS(ixiao)
        call BSLS(5,2,ixiao)
     enddo
     call BSLS(5,3,0)

  case(35)
     do ixiao=0,7
        call SSPS(ixiao)
        call PSSS(ixiao)
     enddo
     do ixiao=0,6
        call PSPS(ixiao)
        call SSDS(ixiao)
        call DSSS(ixiao)
     enddo
     do ixiao=0,5
        call DSPS(ixiao)
        call PSDS(ixiao)
        call SSFS(ixiao)
        call FSSS(ixiao)
     enddo
     do ixiao=0,4
        call FSPS(ixiao)
        call DSDS(ixiao)
        call PSFS(ixiao)
        call SSGS(ixiao)
     enddo
     do ixiao=0,3
        call FSDS(ixiao)
        call DSFS(ixiao)
        call PSGS(ixiao)
        call LSBS(0,5,ixiao)
     enddo
     do ixiao=0,2
        call FSFS(ixiao)
        call DSGS(ixiao)
        call LSBS(1,5,ixiao)
     enddo
     do ixiao=0,1
        call FSGS(ixiao)
        call LSBS(2,5,ixiao)
     enddo
     call LSBS(3,5,0)

  case(54)
     do ixiao=0,8
        call PSSS(ixiao)
        call SSPS(ixiao)
     enddo
     do ixiao=0,7
        call PSPS(ixiao)
        call DSSS(ixiao)
        call SSDS(ixiao)
     enddo
     do ixiao=0,6
        call FSSS(ixiao)
        call PSDS(ixiao)
        call DSPS(ixiao)
        call SSFS(ixiao)
     enddo
     do ixiao=0,5
        call SSGS(ixiao)
        call PSFS(ixiao)
        call DSDS(ixiao)
        call FSPS(ixiao)
        call GSSS(ixiao)
     enddo
     do ixiao=0,4
        call PSGS(ixiao)
        call DSFS(ixiao)
        call FSDS(ixiao)
        call GSPS(ixiao)
        call BSLS(5,0,ixiao)
     enddo
     do ixiao=0,3
        call DSGS(ixiao)
        call FSFS(ixiao)
        call GSDS(ixiao)
        call BSLS(5,1,ixiao)
     enddo
     do ixiao=0,2
        call FSGS(ixiao)
        call GSFS(ixiao)
        call BSLS(5,2,ixiao)
     enddo
     do ixiao=0,1    
        call GSGS(ixiao)
        call BSLS(5,3,ixiao)
     enddo
     call BSLS(5,4,0)

  case(45)
     do ixiao=0,8
        call SSPS(ixiao)
        call PSSS(ixiao)
     enddo
     do ixiao=0,7
        call PSPS(ixiao)
        call DSSS(ixiao)
        call SSDS(ixiao)
     enddo
     do ixiao=0,6
        call FSSS(ixiao)
        call PSDS(ixiao)
        call DSPS(ixiao)
        call SSFS(ixiao)
     enddo
     do ixiao=0,5
        call SSGS(ixiao)
        call PSFS(ixiao)
        call DSDS(ixiao)
        call FSPS(ixiao)
        call GSSS(ixiao)
     enddo
     do ixiao=0,4
        call GSPS(ixiao)
        call FSDS(ixiao)
        call DSFS(ixiao)
        call PSGS(ixiao)
        call LSBS(0,5,ixiao)
     enddo
     do ixiao=0,3
        call GSDS(ixiao)
        call FSFS(ixiao)
        call DSGS(ixiao)
        call LSBS(1,5,ixiao)
     enddo
     do ixiao=0,2
        call GSFS(ixiao)
        call FSGS(ixiao)
        call LSBS(2,5,ixiao)
     enddo
     do ixiao=0,1
        call GSGS(ixiao)
        call LSBS(3,5,ixiao)
     enddo
     call LSBS(4,5,0)

  case(55)
     do ixiao=0,9
        call PSSS(ixiao)
        call SSPS(ixiao)
     enddo
     do ixiao=0,8
        call PSPS(ixiao)
        call DSSS(ixiao)
        call SSDS(ixiao)
     enddo
     do ixiao=0,7
        call FSSS(ixiao)
        call PSDS(ixiao)
        call DSPS(ixiao)
        call SSFS(ixiao)
     enddo
     do ixiao=0,6
        call SSGS(ixiao)
        call PSFS(ixiao)
        call DSDS(ixiao)
        call FSPS(ixiao)
        call GSSS(ixiao)
     enddo
     do ixiao=0,5
        call LSBS(0,5,ixiao)
        call PSGS(ixiao)
        call DSFS(ixiao)
        call FSDS(ixiao)
        call GSPS(ixiao)
        call BSLS(5,0,ixiao)
     enddo
     do ixiao=0,4
        call LSBS(1,5,ixiao)
        call DSGS(ixiao)
        call FSFS(ixiao)
        call GSDS(ixiao)
        call BSLS(5,1,ixiao)
     enddo
     do ixiao=0,3
        call LSBS(2,5,ixiao)
        call FSGS(ixiao)
        call GSFS(ixiao)
        call BSLS(5,2,ixiao)
     enddo
     do ixiao=0,2
        call LSBS(3,5,ixiao)
        call GSGS(ixiao)
        call BSLS(5,3,ixiao)
     enddo
     do ixiao=0,1
        call LSBS(4,5,ixiao)
        call BSLS(5,4,ixiao)
     enddo
     call BSLS(5,5,0)

  case(60)
     do ixiao=0,5
        call PSSS(ixiao)
     enddo
     do ixiao=0,4
        call DSSS(ixiao)
     enddo
     do ixiao=0,3
        call FSSS(ixiao)
     enddo
     do ixiao=0,2
        call GSSS(ixiao)
     enddo
     do ixiao=0,1
        call BSLS(5,0,ixiao)
     enddo
     call BSLS(6,0,0)

  case(6)
     do ixiao=0,5
        call SSPS(ixiao)
     enddo
     do ixiao=0,4
        call SSDS(ixiao)
     enddo
     do ixiao=0,3
        call SSFS(ixiao)
     enddo
     do ixiao=0,2
        call SSGS(ixiao)
     enddo
     do ixiao=0,1
        call LSBS(0,5,ixiao)
     enddo
     call LSBS(0,6,0)

  case(61)
     do ixiao=0,6
        call PSSS(ixiao)
        call SSPS(ixiao)
     enddo
     do ixiao=0,5
        call PSPS(ixiao)
        call DSSS(ixiao)
     enddo
     do ixiao=0,4
        call DSPS(ixiao)
        call FSSS(ixiao)
     enddo
     do ixiao=0,3
        call FSPS(ixiao)
        call GSSS(ixiao)
     enddo
     do ixiao=0,2
        call GSPS(ixiao)
        call BSLS(5,0,ixiao)
     enddo
     do ixiao=0,1
        call BSLS(5,1,ixiao)
        call BSLS(6,0,ixiao)
     enddo
     call BSLS(6,1,0)

  case(16)
     do ixiao=0,6
        call SSPS(ixiao)
        call PSSS(ixiao)
     enddo
     do ixiao=0,5
        call PSPS(ixiao)
        call SSDS(ixiao)
     enddo
     do ixiao=0,4
        call PSDS(ixiao)
        call SSFS(ixiao)
     enddo
     do ixiao=0,3
        call PSFS(ixiao)
        call SSGS(ixiao)
     enddo
     do ixiao=0,2
        call PSGS(ixiao)
        call LSBS(0,5,ixiao)
     enddo
     do ixiao=0,1
        call LSBS(1,5,ixiao)
        call LSBS(0,6,ixiao)
     enddo
     call LSBS(1,6,0)

  case(62)
     do ixiao=0,7
        call PSSS(ixiao)
        call SSPS(ixiao)
     enddo
     do ixiao=0,6
        call PSPS(ixiao)
        call DSSS(ixiao)
        call SSDS(ixiao)
     enddo
     do ixiao=0,5
        call PSDS(ixiao)
        call DSPS(ixiao)
        call FSSS(ixiao)
     enddo
     do ixiao=0,4
        call DSDS(ixiao)
        call FSPS(ixiao)
        call GSSS(ixiao)
     enddo
     do ixiao=0,3
        call FSDS(ixiao)
        call GSPS(ixiao)
        call BSLS(5,0,ixiao)
     enddo
     do ixiao=0,2
        call GSDS(ixiao) 
        call BSLS(5,1,ixiao)
        call BSLS(6,0,ixiao)
     enddo
     do ixiao=0,1
        call BSLS(5,2,ixiao)
        call BSLS(6,1,ixiao)
     enddo
     call BSLS(6,2,0)

  case(26)
     do ixiao=0,7
        call SSPS(ixiao)
        call PSSS(ixiao)
     enddo
     do ixiao=0,6
        call PSPS(ixiao)
        call SSDS(ixiao)
        call DSSS(ixiao)
     enddo
     do ixiao=0,5
        call DSPS(ixiao)
        call PSDS(ixiao)
        call SSFS(ixiao)
     enddo
     do ixiao=0,4
        call DSDS(ixiao)
        call PSFS(ixiao)
        call SSGS(ixiao)
     enddo
     do ixiao=0,3
        call DSFS(ixiao)
        call PSGS(ixiao)
        call LSBS(0,5,ixiao)
     enddo
     do ixiao=0,2
        call DSGS(ixiao)
        call LSBS(1,5,ixiao)
        call LSBS(0,6,ixiao)
     enddo
     do ixiao=0,1
        call LSBS(2,5,ixiao)
        call LSBS(1,6,ixiao)
     enddo
     call LSBS(2,6,0)

  case(63)
     do ixiao=0,8
        call PSSS(ixiao)
        call SSPS(ixiao)
     enddo
     do ixiao=0,7
        call PSPS(ixiao)
        call DSSS(ixiao)
        call SSDS(ixiao)
     enddo
     do ixiao=0,6
        call PSDS(ixiao)
        call DSPS(ixiao)
        call FSSS(ixiao)
        call SSFS(ixiao)
     enddo
     do ixiao=0,5
        call PSFS(ixiao)
        call DSDS(ixiao)
        call FSPS(ixiao)
        call GSSS(ixiao)
     enddo
     do ixiao=0,4
        call DSFS(ixiao)
        call FSDS(ixiao)
        call GSPS(ixiao)
        call BSLS(5,0,ixiao)
     enddo
     do ixiao=0,3
        call FSFS(ixiao)
        call GSDS(ixiao)
        call BSLS(5,1,ixiao)
        call BSLS(6,0,ixiao)
     enddo
     do ixiao=0,2
        call GSFS(ixiao)
        call BSLS(5,2,ixiao)
        call BSLS(6,1,ixiao)
     enddo
     do ixiao=0,1
        call BSLS(5,3,ixiao)
        call BSLS(6,2,ixiao)
     enddo
     call BSLS(6,3,0)

  case(36)
     do ixiao=0,8
        call SSPS(ixiao)
        call PSSS(ixiao)
     enddo
     do ixiao=0,7
        call PSPS(ixiao)
        call SSDS(ixiao)
        call DSSS(ixiao)
     enddo
     do ixiao=0,6
        call DSPS(ixiao)
        call PSDS(ixiao)
        call SSFS(ixiao)
        call FSSS(ixiao)
     enddo
     do ixiao=0,5
        call FSPS(ixiao)
        call DSDS(ixiao)
        call PSFS(ixiao)
        call SSGS(ixiao)
     enddo
     do ixiao=0,4
        call FSDS(ixiao)
        call DSFS(ixiao)
        call PSGS(ixiao)
        call LSBS(0,5,ixiao)
     enddo
     do ixiao=0,3
        call FSFS(ixiao)
        call DSGS(ixiao)
        call LSBS(1,5,ixiao)
        call LSBS(0,6,ixiao)
     enddo
     do ixiao=0,2
        call FSGS(ixiao)
        call LSBS(2,5,ixiao)
        call LSBS(1,6,ixiao)
     enddo
     do ixiao=0,1
        call LSBS(3,5,ixiao)
        call LSBS(2,6,ixiao)
     enddo
     call LSBS(3,6,0)

  case(64)
     do ixiao=0,9
        call PSSS(ixiao)
        call SSPS(ixiao)
     enddo
     do ixiao=0,8
        call PSPS(ixiao)
        call DSSS(ixiao)
        call SSDS(ixiao)
     enddo
     do ixiao=0,7
        call FSSS(ixiao)
        call PSDS(ixiao)
        call DSPS(ixiao)
        call SSFS(ixiao)
     enddo
     do ixiao=0,6
        call SSGS(ixiao)
        call PSFS(ixiao)
        call DSDS(ixiao)
        call FSPS(ixiao)
        call GSSS(ixiao)
     enddo
     do ixiao=0,5
        call PSGS(ixiao)
        call DSFS(ixiao)
        call FSDS(ixiao)
        call GSPS(ixiao)
        call BSLS(5,0,ixiao)
     enddo
     do ixiao=0,4
        call DSGS(ixiao)
        call FSFS(ixiao)
        call GSDS(ixiao)
        call BSLS(5,1,ixiao)
        call BSLS(6,0,ixiao)
     enddo
     do ixiao=0,3
        call FSGS(ixiao)
        call GSFS(ixiao)
        call BSLS(5,2,ixiao)
        call BSLS(6,1,ixiao)
     enddo
     do ixiao=0,2    
        call GSGS(ixiao)
        call BSLS(5,3,ixiao)
        call BSLS(6,2,ixiao)
     enddo
     do ixiao=0,1
        call BSLS(5,4,ixiao)
        call BSLS(6,3,ixiao)
     enddo
     call BSLS(6,4,0)

  case(46)
     do ixiao=0,9
        call SSPS(ixiao)
        call PSSS(ixiao)
     enddo
     do ixiao=0,8
        call PSPS(ixiao)
        call DSSS(ixiao)
        call SSDS(ixiao)
     enddo
     do ixiao=0,7
        call FSSS(ixiao)
        call PSDS(ixiao)
        call DSPS(ixiao)
        call SSFS(ixiao)
     enddo
     do ixiao=0,6
        call SSGS(ixiao)
        call PSFS(ixiao)
        call DSDS(ixiao)
        call FSPS(ixiao)
        call GSSS(ixiao)
     enddo
     do ixiao=0,5
        call GSPS(ixiao)
        call FSDS(ixiao)
        call DSFS(ixiao)
        call PSGS(ixiao)
        call LSBS(0,5,ixiao)
     enddo
     do ixiao=0,4
        call GSDS(ixiao)
        call FSFS(ixiao)
        call DSGS(ixiao)
        call LSBS(1,5,ixiao)
        call LSBS(0,6,ixiao)
     enddo
     do ixiao=0,3
        call GSFS(ixiao)
        call FSGS(ixiao)
        call LSBS(2,5,ixiao)
        call LSBS(1,6,ixiao)
     enddo
     do ixiao=0,2
        call GSGS(ixiao)
        call LSBS(3,5,ixiao)
        call LSBS(2,6,ixiao)
     enddo
     do ixiao=0,1
        call LSBS(4,5,ixiao)
        call LSBS(3,6,ixiao)
     enddo
     call LSBS(4,6,0)

  case(65)
     do ixiao=0,10
        call PSSS(ixiao)
        call SSPS(ixiao)
     enddo
     do ixiao=0,9
        call PSPS(ixiao)
        call DSSS(ixiao)
        call SSDS(ixiao)
     enddo
     do ixiao=0,8
        call FSSS(ixiao)
        call PSDS(ixiao)
        call DSPS(ixiao)
        call SSFS(ixiao)
     enddo
     do ixiao=0,7
        call SSGS(ixiao)
        call PSFS(ixiao)
        call DSDS(ixiao)
        call FSPS(ixiao)
        call GSSS(ixiao)
     enddo
     do ixiao=0,6
        call LSBS(0,5,ixiao)
        call PSGS(ixiao)
        call DSFS(ixiao)
        call FSDS(ixiao)
        call GSPS(ixiao)
        call BSLS(5,0,ixiao)
     enddo
     do ixiao=0,5
        call LSBS(1,5,ixiao)
        call DSGS(ixiao)
        call FSFS(ixiao)
        call GSDS(ixiao)
        call BSLS(5,1,ixiao)
        call BSLS(6,0,ixiao)
     enddo
     do ixiao=0,4
        call LSBS(2,5,ixiao)
        call FSGS(ixiao)
        call GSFS(ixiao)
        call BSLS(5,2,ixiao)
        call BSLS(6,1,ixiao)
     enddo
     do ixiao=0,3
        call LSBS(3,5,ixiao)
        call GSGS(ixiao)
        call BSLS(5,3,ixiao)
        call BSLS(6,2,ixiao)
     enddo
     do ixiao=0,2
        call LSBS(4,5,ixiao)
        call BSLS(5,4,ixiao)
        call BSLS(6,3,ixiao)
     enddo
     do ixiao=0,1
        call BSLS(5,5,ixiao)
        call BSLS(6,4,ixiao)
     enddo
     call BSLS(6,5,0)

  case(56)
     do ixiao=0,10
        call PSSS(ixiao)
        call SSPS(ixiao)
     enddo
     do ixiao=0,9
        call PSPS(ixiao)
        call DSSS(ixiao)
        call SSDS(ixiao)
     enddo
     do ixiao=0,8
        call FSSS(ixiao)
        call PSDS(ixiao)
        call DSPS(ixiao)
        call SSFS(ixiao)
     enddo
     do ixiao=0,7
        call SSGS(ixiao)
        call PSFS(ixiao)
        call DSDS(ixiao)
        call FSPS(ixiao)
        call GSSS(ixiao)
     enddo
     do ixiao=0,6
        call LSBS(0,5,ixiao)
        call PSGS(ixiao)
        call DSFS(ixiao)
        call FSDS(ixiao)
        call GSPS(ixiao)
        call BSLS(5,0,ixiao)
     enddo
     do ixiao=0,5
        call BSLS(5,1,ixiao)
        call GSDS(ixiao)
        call FSFS(ixiao)
        call DSGS(ixiao)
        call LSBS(1,5,ixiao)
        call LSBS(0,6,ixiao)
     enddo
     do ixiao=0,4
        call BSLS(5,2,ixiao)
        call FSGS(ixiao)
        call GSFS(ixiao)
        call LSBS(2,5,ixiao)
        call LSBS(1,6,ixiao)
     enddo
     do ixiao=0,3
        call BSLS(5,3,ixiao)
        call GSGS(ixiao)
        call LSBS(3,5,ixiao)
        call LSBS(2,6,ixiao)
     enddo
     do ixiao=0,2
        call BSLS(5,4,ixiao)
        call LSBS(4,5,ixiao)
        call LSBS(3,6,ixiao)
     enddo
     do ixiao=0,1
        call BSLS(5,5,ixiao)
        call LSBS(4,6,ixiao)
     enddo
     call LSBS(5,6,0)

  case(66)
     do ixiao=0,11
        call PSSS(ixiao)
        call SSPS(ixiao)
     enddo
     do ixiao=0,10
        call PSPS(ixiao)
        call DSSS(ixiao)
        call SSDS(ixiao)
     enddo
     do ixiao=0,9
        call FSSS(ixiao)
        call PSDS(ixiao)
        call DSPS(ixiao)
        call SSFS(ixiao)
     enddo
     do ixiao=0,8
        call SSGS(ixiao)
        call PSFS(ixiao)
        call DSDS(ixiao)
        call FSPS(ixiao)
        call GSSS(ixiao)
     enddo
     do ixiao=0,7
        call LSBS(0,5,ixiao)
        call PSGS(ixiao)
        call DSFS(ixiao)
        call FSDS(ixiao)
        call GSPS(ixiao)
        call BSLS(5,0,ixiao)
     enddo
     do ixiao=0,6
        call LSBS(0,6,ixiao)
        call LSBS(1,5,ixiao)
        call DSGS(ixiao)
        call FSFS(ixiao)
        call GSDS(ixiao)
        call BSLS(5,1,ixiao)
        call BSLS(6,0,ixiao)
     enddo
     do ixiao=0,5
        call LSBS(1,6,ixiao)
        call LSBS(2,5,ixiao)
        call FSGS(ixiao)
        call GSFS(ixiao)
        call BSLS(5,2,ixiao)
        call BSLS(6,1,ixiao)
     enddo
     do ixiao=0,4
        call LSBS(2,6,ixiao)
        call LSBS(3,5,ixiao)
        call GSGS(ixiao)
        call BSLS(5,3,ixiao)
        call BSLS(6,2,ixiao)
     enddo
     do ixiao=0,3
        call LSBS(3,6,ixiao)
        call LSBS(4,5,ixiao)
        call BSLS(5,4,ixiao)
        call BSLS(6,3,ixiao)
     enddo
     do ixiao=0,2
        call LSBS(4,6,ixiao)
        call BSLS(5,5,ixiao)
        call BSLS(6,4,ixiao)
     enddo
     do ixiao=0,1
        call LSBS(5,6,ixiao)
        call BSLS(6,5,ixiao)
     enddo
     call BSLS(6,6,0)

  case(70)
     do ixiao=0,6
        call PSSS(ixiao)
     enddo
     do ixiao=0,5
        call DSSS(ixiao)
     enddo
     do ixiao=0,4
        call FSSS(ixiao)
     enddo
     do ixiao=0,3
        call GSSS(ixiao)
     enddo
     do ixiao=0,2
        call BSLS(5,0,ixiao)
     enddo
     do ixiao=0,1
        call BSLS(6,0,ixiao)
     enddo
     call BSLS(7,0,0)

  case(7)
     do ixiao=0,6
        call SSPS(ixiao)
     enddo
     do ixiao=0,5
        call SSDS(ixiao)
     enddo
     do ixiao=0,4
        call SSFS(ixiao)
     enddo
     do ixiao=0,3
        call SSGS(ixiao)
     enddo
     do ixiao=0,2
        call LSBS(0,5,ixiao)
     enddo
     do ixiao=0,1
        call LSBS(0,6,ixiao)
     enddo
     call LSBS(0,7,0)

  case(71)
     do ixiao=0,7
        call PSSS(ixiao)
        call SSPS(ixiao)
     enddo
     do ixiao=0,6
        call PSPS(ixiao)
        call DSSS(ixiao)
     enddo
     do ixiao=0,5
        call DSPS(ixiao)
        call FSSS(ixiao)
     enddo
     do ixiao=0,4
        call FSPS(ixiao)
        call GSSS(ixiao)
     enddo
     do ixiao=0,3
        call GSPS(ixiao)
        call BSLS(5,0,ixiao)
     enddo
     do ixiao=0,2
        call BSLS(5,1,ixiao)
        call BSLS(6,0,ixiao)
     enddo
     do ixiao=0,1
        call BSLS(6,1,ixiao)
        call BSLS(7,0,ixiao)
     enddo
     call BSLS(7,1,0)

  case(17)
     do ixiao=0,7
        call SSPS(ixiao)
        call PSSS(ixiao)
     enddo
     do ixiao=0,6
        call PSPS(ixiao)
        call SSDS(ixiao)
     enddo
     do ixiao=0,5
        call PSDS(ixiao)
        call SSFS(ixiao)
     enddo
     do ixiao=0,4
        call PSFS(ixiao)
        call SSGS(ixiao)
     enddo
     do ixiao=0,3
        call PSGS(ixiao)
        call LSBS(0,5,ixiao)
     enddo
     do ixiao=0,2
        call LSBS(1,5,ixiao)
        call LSBS(0,6,ixiao)
     enddo
     do ixiao=0,1
        call LSBS(1,6,ixiao)
        call LSBS(0,7,ixiao)
     enddo
     call LSBS(1,7,0)

  case(72)
     do ixiao=0,8
        call PSSS(ixiao)
        call SSPS(ixiao)
     enddo
     do ixiao=0,7
        call PSPS(ixiao)
        call DSSS(ixiao)
        call SSDS(ixiao)
     enddo
     do ixiao=0,6
        call PSDS(ixiao)
        call DSPS(ixiao)
        call FSSS(ixiao)
     enddo
     do ixiao=0,5
        call DSDS(ixiao)
        call FSPS(ixiao)
        call GSSS(ixiao)
     enddo
     do ixiao=0,4
        call FSDS(ixiao)
        call GSPS(ixiao)
        call BSLS(5,0,ixiao)
     enddo
     do ixiao=0,3
        call GSDS(ixiao) 
        call BSLS(5,1,ixiao)
        call BSLS(6,0,ixiao)
     enddo
     do ixiao=0,2
        call BSLS(5,2,ixiao)
        call BSLS(6,1,ixiao)
        call BSLS(7,0,ixiao)
     enddo
     do ixiao=0,1
        call BSLS(6,2,ixiao)
        call BSLS(7,1,ixiao)
     enddo
     call BSLS(7,2,0)

  case(27)
     do ixiao=0,8
        call SSPS(ixiao)
        call PSSS(ixiao)
     enddo
     do ixiao=0,7
        call PSPS(ixiao)
        call SSDS(ixiao)
        call DSSS(ixiao)
     enddo
     do ixiao=0,6
        call DSPS(ixiao)
        call PSDS(ixiao)
        call SSFS(ixiao)
     enddo
     do ixiao=0,5
        call DSDS(ixiao)
        call PSFS(ixiao)
        call SSGS(ixiao)
     enddo
     do ixiao=0,4
        call DSFS(ixiao)
        call PSGS(ixiao)
        call LSBS(0,5,ixiao)
     enddo
     do ixiao=0,3
        call DSGS(ixiao)
        call LSBS(1,5,ixiao)
        call LSBS(0,6,ixiao)
     enddo
     do ixiao=0,2
        call LSBS(2,5,ixiao)
        call LSBS(1,6,ixiao)
        call LSBS(0,7,ixiao)
     enddo
     do ixiao=0,1
        call LSBS(2,6,ixiao)
        call LSBS(1,7,ixiao)
     enddo
     call LSBS(2,7,0)

  case(73)
     do ixiao=0,9
        call PSSS(ixiao)
        call SSPS(ixiao)
     enddo
     do ixiao=0,8
        call PSPS(ixiao)
        call DSSS(ixiao)
        call SSDS(ixiao)
     enddo
     do ixiao=0,7
        call PSDS(ixiao)
        call DSPS(ixiao)
        call FSSS(ixiao)
        call SSFS(ixiao)
     enddo
     do ixiao=0,6
        call PSFS(ixiao)
        call DSDS(ixiao)
        call FSPS(ixiao)
        call GSSS(ixiao)
     enddo
     do ixiao=0,5
        call DSFS(ixiao)
        call FSDS(ixiao)
        call GSPS(ixiao)
        call BSLS(5,0,ixiao)
     enddo
     do ixiao=0,4
        call FSFS(ixiao)
        call GSDS(ixiao)
        call BSLS(5,1,ixiao)
        call BSLS(6,0,ixiao)
     enddo
     do ixiao=0,3
        call GSFS(ixiao)
        call BSLS(5,2,ixiao)
        call BSLS(6,1,ixiao)
        call BSLS(7,0,ixiao)
     enddo
     do ixiao=0,2
        call BSLS(5,3,ixiao)
        call BSLS(6,2,ixiao)
        call BSLS(7,1,ixiao)
     enddo
     do ixiao=0,1
        call BSLS(6,3,ixiao)
        call BSLS(7,2,ixiao)
     enddo
     call BSLS(7,3,0)

  case(37)
     do ixiao=0,9
        call SSPS(ixiao)
        call PSSS(ixiao)
     enddo
     do ixiao=0,8
        call PSPS(ixiao)
        call SSDS(ixiao)
        call DSSS(ixiao)
     enddo
     do ixiao=0,7
        call DSPS(ixiao)
        call PSDS(ixiao)
        call SSFS(ixiao)
        call FSSS(ixiao)
     enddo
     do ixiao=0,6
        call FSPS(ixiao)
        call DSDS(ixiao)
        call PSFS(ixiao)
        call SSGS(ixiao)
     enddo
     do ixiao=0,5
        call FSDS(ixiao)
        call DSFS(ixiao)
        call PSGS(ixiao)
        call LSBS(0,5,ixiao)
     enddo
     do ixiao=0,4
        call FSFS(ixiao)
        call DSGS(ixiao)
        call LSBS(1,5,ixiao)
        call LSBS(0,6,ixiao)
     enddo
     do ixiao=0,3
        call FSGS(ixiao)
        call LSBS(2,5,ixiao)
        call LSBS(1,6,ixiao)
        call LSBS(0,7,ixiao)
     enddo
     do ixiao=0,2
        call LSBS(3,5,ixiao)
        call LSBS(2,6,ixiao)
        call LSBS(1,7,ixiao)
     enddo
     do ixiao=0,1
        call LSBS(3,6,ixiao)
        call LSBS(2,7,ixiao)
     enddo
     call LSBS(3,7,0)

  case(74)
     do ixiao=0,10
        call PSSS(ixiao)
        call SSPS(ixiao)
     enddo
     do ixiao=0,9
        call PSPS(ixiao)
        call DSSS(ixiao)
        call SSDS(ixiao)
     enddo
     do ixiao=0,8
        call FSSS(ixiao)
        call PSDS(ixiao)
        call DSPS(ixiao)
        call SSFS(ixiao)
     enddo
     do ixiao=0,7
        call SSGS(ixiao)
        call PSFS(ixiao)
        call DSDS(ixiao)
        call FSPS(ixiao)
        call GSSS(ixiao)
     enddo
     do ixiao=0,6
        call PSGS(ixiao)
        call DSFS(ixiao)
        call FSDS(ixiao)
        call GSPS(ixiao)
        call BSLS(5,0,ixiao)
     enddo
     do ixiao=0,5
        call DSGS(ixiao)
        call FSFS(ixiao)
        call GSDS(ixiao)
        call BSLS(5,1,ixiao)
        call BSLS(6,0,ixiao)
     enddo
     do ixiao=0,4
        call FSGS(ixiao)
        call GSFS(ixiao)
        call BSLS(5,2,ixiao)
        call BSLS(6,1,ixiao)
        call BSLS(7,0,ixiao)
     enddo
     do ixiao=0,3    
        call GSGS(ixiao)
        call BSLS(5,3,ixiao)
        call BSLS(6,2,ixiao)
        call BSLS(7,1,ixiao)
     enddo
     do ixiao=0,2
        call BSLS(5,4,ixiao)
        call BSLS(6,3,ixiao)
        call BSLS(7,2,ixiao)
     enddo
     do ixiao=0,1
        call BSLS(6,4,ixiao)
        call BSLS(7,3,ixiao)
     enddo
     call BSLS(7,4,0)

  case(47)
     do ixiao=0,10
        call SSPS(ixiao)
        call PSSS(ixiao)
     enddo
     do ixiao=0,9
        call PSPS(ixiao)
        call DSSS(ixiao)
        call SSDS(ixiao)
     enddo
     do ixiao=0,8
        call FSSS(ixiao)
        call PSDS(ixiao)
        call DSPS(ixiao)
        call SSFS(ixiao)
     enddo
     do ixiao=0,7
        call SSGS(ixiao)
        call PSFS(ixiao)
        call DSDS(ixiao)
        call FSPS(ixiao)
        call GSSS(ixiao)
     enddo
     do ixiao=0,6
        call GSPS(ixiao)
        call FSDS(ixiao)
        call DSFS(ixiao)
        call PSGS(ixiao)
        call LSBS(0,5,ixiao)
     enddo
     do ixiao=0,5
        call GSDS(ixiao)
        call FSFS(ixiao)
        call DSGS(ixiao)
        call LSBS(1,5,ixiao)
        call LSBS(0,6,ixiao)
     enddo
     do ixiao=0,4
        call GSFS(ixiao)
        call FSGS(ixiao)
        call LSBS(2,5,ixiao)
        call LSBS(1,6,ixiao)
        call LSBS(0,7,ixiao)
     enddo
     do ixiao=0,3
        call GSGS(ixiao)
        call LSBS(3,5,ixiao)
        call LSBS(2,6,ixiao)
        call LSBS(1,7,ixiao)
     enddo
     do ixiao=0,2
        call LSBS(4,5,ixiao)
        call LSBS(3,6,ixiao)
        call LSBS(2,7,ixiao)
     enddo
     do ixiao=0,1
        call LSBS(4,6,ixiao)
        call LSBS(3,7,ixiao)
     enddo
     call LSBS(4,7,0)

  case(75)
     do ixiao=0,11
        call PSSS(ixiao)
        call SSPS(ixiao)
     enddo
     do ixiao=0,10
        call PSPS(ixiao)
        call DSSS(ixiao)
        call SSDS(ixiao)
     enddo
     do ixiao=0,9
        call FSSS(ixiao)
        call PSDS(ixiao)
        call DSPS(ixiao)
        call SSFS(ixiao)
     enddo
     do ixiao=0,8
        call SSGS(ixiao)
        call PSFS(ixiao)
        call DSDS(ixiao)
        call FSPS(ixiao)
        call GSSS(ixiao)
     enddo
     do ixiao=0,7
        call LSBS(0,5,ixiao)
        call PSGS(ixiao)
        call DSFS(ixiao)
        call FSDS(ixiao)
        call GSPS(ixiao)
        call BSLS(5,0,ixiao)
     enddo
     do ixiao=0,6
        call LSBS(1,5,ixiao)
        call DSGS(ixiao)
        call FSFS(ixiao)
        call GSDS(ixiao)
        call BSLS(5,1,ixiao)
        call BSLS(6,0,ixiao)
     enddo
     do ixiao=0,5
        call LSBS(2,5,ixiao)
        call FSGS(ixiao)
        call GSFS(ixiao)
        call BSLS(5,2,ixiao)
        call BSLS(6,1,ixiao)
        call BSLS(7,0,ixiao)
     enddo
     do ixiao=0,4
        call LSBS(3,5,ixiao)
        call GSGS(ixiao)
        call BSLS(5,3,ixiao)
        call BSLS(6,2,ixiao)
        call BSLS(7,1,ixiao)
     enddo
     do ixiao=0,3
        call LSBS(4,5,ixiao)
        call BSLS(5,4,ixiao)
        call BSLS(6,3,ixiao)
        call BSLS(7,2,ixiao)
     enddo
     do ixiao=0,2
        call BSLS(5,5,ixiao)
        call BSLS(6,4,ixiao)
        call BSLS(7,3,ixiao)
     enddo
     do ixiao=0,1
        call BSLS(6,5,ixiao)
        call BSLS(7,4,ixiao)
     enddo
     call BSLS(7,5,0)

  case(57)
     do ixiao=0,11
        call PSSS(ixiao)
        call SSPS(ixiao)
     enddo
     do ixiao=0,10
        call PSPS(ixiao)
        call DSSS(ixiao)
        call SSDS(ixiao)
     enddo
     do ixiao=0,9
        call FSSS(ixiao)
        call PSDS(ixiao)
        call DSPS(ixiao)
        call SSFS(ixiao)
     enddo
     do ixiao=0,8
        call SSGS(ixiao)
        call PSFS(ixiao)
        call DSDS(ixiao)
        call FSPS(ixiao)
        call GSSS(ixiao)
     enddo
     do ixiao=0,7
        call LSBS(0,5,ixiao)
        call PSGS(ixiao)
        call DSFS(ixiao)
        call FSDS(ixiao)
        call GSPS(ixiao)
        call BSLS(5,0,ixiao)
     enddo
     do ixiao=0,6
        call BSLS(5,1,ixiao)
        call GSDS(ixiao)
        call FSFS(ixiao)
        call DSGS(ixiao)
        call LSBS(1,5,ixiao)
        call LSBS(0,6,ixiao)
     enddo
     do ixiao=0,5
        call BSLS(5,2,ixiao)
        call FSGS(ixiao)
        call GSFS(ixiao)
        call LSBS(2,5,ixiao)
        call LSBS(1,6,ixiao)
        call LSBS(0,7,ixiao)
     enddo
     do ixiao=0,4
        call BSLS(5,3,ixiao)
        call GSGS(ixiao)
        call LSBS(3,5,ixiao)
        call LSBS(2,6,ixiao)
        call LSBS(1,7,ixiao)
     enddo
     do ixiao=0,3
        call BSLS(5,4,ixiao)
        call LSBS(4,5,ixiao)
        call LSBS(3,6,ixiao)
        call LSBS(2,7,ixiao)
     enddo
     do ixiao=0,2
        call BSLS(5,5,ixiao)
        call LSBS(4,6,ixiao)
        call LSBS(3,7,ixiao)
     enddo
     do ixiao=0,1
        call LSBS(5,6,ixiao)
        call LSBS(4,7,ixiao)
     enddo
     call LSBS(5,7,0)

  case(76)
     do ixiao=0,12
        call PSSS(ixiao)
        call SSPS(ixiao)
     enddo
     do ixiao=0,11
        call PSPS(ixiao)
        call DSSS(ixiao)
        call SSDS(ixiao)
     enddo
     do ixiao=0,10
        call FSSS(ixiao)
        call PSDS(ixiao)
        call DSPS(ixiao)
        call SSFS(ixiao)
     enddo
     do ixiao=0,9
        call SSGS(ixiao)
        call PSFS(ixiao)
        call DSDS(ixiao)
        call FSPS(ixiao)
        call GSSS(ixiao)
     enddo
     do ixiao=0,8
        call LSBS(0,5,ixiao)
        call PSGS(ixiao)
        call DSFS(ixiao)
        call FSDS(ixiao)
        call GSPS(ixiao)
        call BSLS(5,0,ixiao)
     enddo
     do ixiao=0,7
        call LSBS(0,6,ixiao)
        call LSBS(1,5,ixiao)
        call DSGS(ixiao)
        call FSFS(ixiao)
        call GSDS(ixiao)
        call BSLS(5,1,ixiao)
        call BSLS(6,0,ixiao)
     enddo
     do ixiao=0,6
        call LSBS(1,6,ixiao)
        call LSBS(2,5,ixiao)
        call FSGS(ixiao)
        call GSFS(ixiao)
        call BSLS(5,2,ixiao)
        call BSLS(6,1,ixiao)
        call BSLS(7,0,ixiao)
     enddo
     do ixiao=0,5
        call LSBS(2,6,ixiao)
        call LSBS(3,5,ixiao)
        call GSGS(ixiao)
        call BSLS(5,3,ixiao)
        call BSLS(6,2,ixiao)
        call BSLS(7,1,ixiao)
     enddo
     do ixiao=0,4
        call LSBS(3,6,ixiao)
        call LSBS(4,5,ixiao)
        call BSLS(5,4,ixiao)
        call BSLS(6,3,ixiao)
        call BSLS(7,2,ixiao)
     enddo
     do ixiao=0,3
        call LSBS(4,6,ixiao)
        call BSLS(5,5,ixiao)
        call BSLS(6,4,ixiao)
        call BSLS(7,3,ixiao)
     enddo
     do ixiao=0,2
        call LSBS(5,6,ixiao)
        call BSLS(6,5,ixiao)
        call BSLS(7,4,ixiao)
     enddo
     do ixiao=0,1
        call BSLS(6,6,ixiao)
        call BSLS(7,5,ixiao)
     enddo
     call BSLS(7,6,0)

  case(67)
     do ixiao=0,12
        call PSSS(ixiao)
        call SSPS(ixiao)
     enddo
     do ixiao=0,11
        call PSPS(ixiao)
        call DSSS(ixiao)
        call SSDS(ixiao)
     enddo
     do ixiao=0,10
        call FSSS(ixiao)
        call PSDS(ixiao)
        call DSPS(ixiao)
        call SSFS(ixiao)
     enddo
     do ixiao=0,9
        call SSGS(ixiao)
        call PSFS(ixiao)
        call DSDS(ixiao)
        call FSPS(ixiao)
        call GSSS(ixiao)
     enddo
     do ixiao=0,8
        call LSBS(0,5,ixiao)
        call PSGS(ixiao)
        call DSFS(ixiao)
        call FSDS(ixiao)
        call GSPS(ixiao)
        call BSLS(5,0,ixiao)
     enddo
     do ixiao=0,7
        call LSBS(0,6,ixiao)
        call LSBS(1,5,ixiao)
        call DSGS(ixiao)
        call FSFS(ixiao)
        call GSDS(ixiao)
        call BSLS(5,1,ixiao)
        call BSLS(6,0,ixiao)
     enddo
     do ixiao=0,6
        call BSLS(6,1,ixiao)
        call BSLS(5,2,ixiao)
        call FSGS(ixiao)
        call GSFS(ixiao)
        call LSBS(2,5,ixiao)
        call LSBS(1,6,ixiao)
        call LSBS(0,7,ixiao)
     enddo
     do ixiao=0,5
        call LSBS(2,6,ixiao)
        call LSBS(3,5,ixiao)
        call GSGS(ixiao)
        call BSLS(5,3,ixiao)
        call BSLS(6,2,ixiao)
        call LSBS(1,7,ixiao)
     enddo
     do ixiao=0,4
        call LSBS(3,6,ixiao)
        call LSBS(4,5,ixiao)
        call BSLS(5,4,ixiao)
        call BSLS(6,3,ixiao)
        call LSBS(2,7,ixiao)
     enddo
     do ixiao=0,3
        call LSBS(4,6,ixiao)
        call BSLS(5,5,ixiao)
        call BSLS(6,4,ixiao)
        call LSBS(3,7,ixiao)
     enddo
     do ixiao=0,2
        call LSBS(5,6,ixiao)
        call BSLS(6,5,ixiao)
        call LSBS(4,7,ixiao)
     enddo
     do ixiao=0,1
        call BSLS(6,6,ixiao)
        call LSBS(5,7,ixiao)
     enddo
     call LSBS(6,7,0)

  case(77)
     do ixiao=0,13
        call PSSS(ixiao)
        call SSPS(ixiao)
     enddo
     do ixiao=0,12
        call PSPS(ixiao)
        call DSSS(ixiao)
        call SSDS(ixiao)
     enddo
     do ixiao=0,11
        call FSSS(ixiao)
        call PSDS(ixiao)
        call DSPS(ixiao)
        call SSFS(ixiao)
     enddo
     do ixiao=0,10
        call SSGS(ixiao)
        call PSFS(ixiao)
        call DSDS(ixiao)
        call FSPS(ixiao)
        call GSSS(ixiao)
     enddo
     do ixiao=0,9
        call LSBS(0,5,ixiao)
        call PSGS(ixiao)
        call DSFS(ixiao)
        call FSDS(ixiao)
        call GSPS(ixiao)
        call BSLS(5,0,ixiao)
     enddo
     do ixiao=0,8
        call LSBS(0,6,ixiao)
        call LSBS(1,5,ixiao)
        call DSGS(ixiao)
        call FSFS(ixiao)
        call GSDS(ixiao)
        call BSLS(5,1,ixiao)
        call BSLS(6,0,ixiao)
     enddo
     do ixiao=0,7
        call LSBS(0,7,ixiao)
        call LSBS(1,6,ixiao)
        call LSBS(2,5,ixiao)
        call FSGS(ixiao)
        call GSFS(ixiao)
        call BSLS(5,2,ixiao)
        call BSLS(6,1,ixiao)
        call BSLS(7,0,ixiao)
     enddo
     do ixiao=0,6
        call LSBS(1,7,ixiao)
        call LSBS(2,6,ixiao)
        call LSBS(3,5,ixiao)
        call GSGS(ixiao)
        call BSLS(5,3,ixiao)
        call BSLS(6,2,ixiao)
        call BSLS(7,1,ixiao)
     enddo
     do ixiao=0,5
        call LSBS(2,7,ixiao)
        call LSBS(3,6,ixiao)
        call LSBS(4,5,ixiao)
        call BSLS(5,4,ixiao)
        call BSLS(6,3,ixiao)
        call BSLS(7,2,ixiao)
     enddo
     do ixiao=0,4
        call LSBS(3,7,ixiao)
        call LSBS(4,6,ixiao)
        call BSLS(5,5,ixiao)
        call BSLS(6,4,ixiao)
        call BSLS(7,3,ixiao)
     enddo
     do ixiao=0,3
        call LSBS(4,7,ixiao)
        call LSBS(5,6,ixiao)
        call BSLS(6,5,ixiao)
        call BSLS(7,4,ixiao)
     enddo
     do ixiao=0,2
        call LSBS(5,7,ixiao)
        call BSLS(6,6,ixiao)
        call BSLS(7,5,ixiao)
     enddo
     call LSBS(6,7,0)
     call BSLS(7,6,0)
     call LSBS(6,7,1)
     call BSLS(7,6,1)
     call BSLS(7,7,0)

  end select
end subroutine vertical



! Ed Brothers. October 23, 2001
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

!    subroutine attrashell(a,b,i,j,k,ii,jj,kk,Ax,Ay,Az, &
!    Bx,By,Bz,Cx,Cy,Cz,Z)
subroutine attrashellenergy(IIsh,JJsh)
  use allmod
  !    use xiaoconstants
  implicit double precision(a-h,o-z)
  dimension aux(0:20)
  real*8 AA(3),BB(3),CC(3),PP(3)
  common /xiaoattra/attra,aux,AA,BB,CC,PP,g

  real*8 RA(3),RB(3),RP(3)

  ! Variables needed later:
  !    pi=3.1415926535897932385

  Ax=xyz(1,quick_basis%katom(IIsh))
  Ay=xyz(2,quick_basis%katom(IIsh))
  Az=xyz(3,quick_basis%katom(IIsh))

  Bx=xyz(1,quick_basis%katom(JJsh))
  By=xyz(2,quick_basis%katom(JJsh))
  Bz=xyz(3,quick_basis%katom(JJsh))

  !   Cx=sumx
  !   Cy=sumy
  !   Cz=sumz

  ! The purpose of this subroutine is to calculate the nuclear attraction
  ! of an electron  distributed between gtfs with orbital exponents a
  ! and b on A and B with angular momentums defined by i,j,k (a's x, y
  ! and z exponents, respectively) and ii,jj,k and kk on B with the core at
  ! (Cx,Cy,Cz) with charge Z. m is the "order" of the integral which
  ! arises from the recusion relationship.

  ! The this is taken from the recursive relation found in Obara and Saika,
  ! J. Chem. Phys. 84 (7) 1986, 3963.

  ! The first step is generating all the necessary auxillary integrals.
  ! These are (0|1/rc|0)^(m) = 2 Sqrt (g/Pi) (0||0) Fm(g(Rpc)^2)
  ! The values of m range from 0 to i+j+k+ii+jj+kk.

  NII2=quick_basis%Qfinal(IIsh)
  NJJ2=quick_basis%Qfinal(JJsh)
  Maxm=NII2+NJJ2


  do ips=1,quick_basis%kprim(IIsh)
     a=gcexpo(ips,quick_basis%ksumtype(IIsh))
     do jps=1,quick_basis%kprim(JJsh)
        b=gcexpo(jps,quick_basis%ksumtype(JJsh))       

        g = a+b
        Px = (a*Ax + b*Bx)/g
        Py = (a*Ay + b*By)/g
        Pz = (a*Az + b*Bz)/g

        constant = overlap(a,b,0,0,0,0,0,0,Ax,Ay,Az,Bx,By,Bz) &
             * 2.d0 * sqrt(g/Pi)

        do iatom=1,natom+quick_molspec%nextatom
           if(iatom<=natom)then
              Cx=xyz(1,iatom)
              Cy=xyz(2,iatom)
              Cz=xyz(3,iatom)
              Z=-1.0d0*quick_molspec%chg(iatom)
           else
              Cx=quick_molspec%extxyz(1,iatom-natom)
              Cy=quick_molspec%extxyz(2,iatom-natom)
              Cz=quick_molspec%extxyz(3,iatom-natom)
              Z=-quick_molspec%extchg(iatom-natom)
           endif

           PCsquare = (Px-Cx)**2 + (Py -Cy)**2 + (Pz -Cz)**2

           U = g* PCsquare
           !    Maxm = i+j+k+ii+jj+kk
           call FmT(Maxm,U,aux)
           do L = 0,maxm
              aux(L) = aux(L)*constant*Z
              attraxiao(1,1,L)=aux(L)
           enddo

           ! At this point all the auxillary integrals have been calculated.
           ! It is now time to decompase the attraction integral to it's
           ! auxillary integrals through the recursion scheme.  To do this we use
           ! a recursive function.

           !    attraction = attrecurse(i,j,k,ii,jj,kk,0,aux,Ax,Ay,Az,Bx,By,Bz, &
           !    Cx,Cy,Cz,Px,Py,Pz,g)
           NIJ1=10*NII2+NJJ2

           call nuclearattraenergy(ips,jps,IIsh,JJsh,NIJ1,Ax,Ay,Az,Bx,By,Bz, &
                Cx,Cy,Cz,Px,Py,Pz,iatom)

        enddo

     enddo
  enddo

  ! Xiao HE remember to multiply Z   01/12/2008
  !    attraction = attraction*(-1.d0)* Z
  return
end subroutine attrashellenergy



! Vertical Recursion by Xiao HE 07/07/07 version
subroutine shellmp2(nstepmp2s,nsteplength)
  use allmod

  Implicit real*8(a-h,o-z)
  real*8 P(3),Q(3),W(3),KAB,KCD,AAtemp(3)
  Parameter(NN=13)
  real*8 FM(0:13)
  real*8 RA(3),RB(3),RC(3),RD(3)

  real*8 Qtemp(3),WQtemp(3),CDtemp,ABcom,Ptemp(3),WPtemp(3),ABtemp,CDcom,ABCDtemp
  integer II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2
  common /hrrstore/II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2

  COMMON /VRRcom/Qtemp,WQtemp,CDtemp,ABcom,Ptemp,WPtemp,ABtemp,CDcom,ABCDtemp

  COMMON /COM1/RA,RB,RC,RD

  Do M=1,3
     RA(M)=xyz(M,quick_basis%katom(II))
     RB(M)=xyz(M,quick_basis%katom(JJ))
     RC(M)=xyz(M,quick_basis%katom(KK))
     RD(M)=xyz(M,quick_basis%katom(LL))
  Enddo

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
  Do JJJ=1,quick_basis%kprim(JJ)
     Nprij=quick_basis%kstart(JJ)+JJJ-1
     Do III=1,quick_basis%kprim(II)
        Nprii=quick_basis%kstart(II)+III-1
        AB=Apri(Nprii,Nprij)
        ABtemp=0.5d0/AB
        cutoffprim1=dnmax*cutprim(Nprii,Nprij)
        Do M=1,3
           P(M)=Ppri(M,Nprii,Nprij)
           AAtemp(M)=P(M)*AB
           Ptemp(M)=P(M)-RA(M)
        Enddo
        !            KAB=Kpri(Nprii,Nprij)
        Do LLL=1,quick_basis%kprim(LL)
           Npril=quick_basis%kstart(LL)+LLL-1
           Do KKK=1,quick_basis%kprim(KK)
              Nprik=quick_basis%kstart(KK)+KKK-1
              cutoffprim=cutoffprim1*cutprim(Nprik,Npril)
              !                       print*,cutoffprim
              If(cutoffprim.gt.quick_method%primLimit)then
                 CD=Apri(Nprik,Npril)
                 ABCD=AB+CD
                 ROU=AB*CD/ABCD
                 RPQ=0.0d0
                 ABCDxiao=dsqrt(ABCD) 

                 CDtemp=0.5d0/CD
                 ABcom=AB/ABCD
                 CDcom=CD/ABCD
                 ABCDtemp=0.5d0/ABCD

                 Do M=1,3
                    Q(M)=Ppri(M,Nprik,Npril)
                    W(M)=(AAtemp(M)+Q(M)*CD)/ABCD
                    XXXtemp=P(M)-Q(M)
                    RPQ=RPQ+XXXtemp*XXXtemp
                    Qtemp(M)=Q(M)-RC(M)
                    WQtemp(M)=W(M)-Q(M)
                    WPtemp(M)=W(M)-P(M)
                 Enddo
                 !                         KCD=Kpri(Nprik,Npril)

                 T=RPQ*ROU

                 !                         NABCD=0 
                 !                         call FmT(0,T,FM)
                 !                         do iitemp=0,0
                 !                           Yxiaotemp(1,1,iitemp)=FM(iitemp)/ABCDxiao
                 !                         enddo
                 call FmT(NABCD,T,FM)
                 do iitemp=0,NABCD
                    Yxiaotemp(1,1,iitemp)=FM(iitemp)/ABCDxiao
                 enddo
                 !                         if(II.eq.1.and.JJ.eq.4.and.KK.eq.10.and.LL.eq.16)then
                 !                          print*,III,JJJ,KKK,LLL,T,NABCD,FM(0:NABCD)
                 !                         endif 
                 !                         print*,III,JJJ,KKK,LLL,FM
                 ITT=ITT+1

                 call vertical(NABCDTYPE) 

                 Do I2=NNC,NNCD
                    Do I1=NNA,NNAB
                       Yxiao(ITT,I1,I2)=Yxiaotemp(I1,I2,0)
                    Enddo
                 Enddo
                 !                           else
                 !!                             print*,cutoffprim
                 !                             ITT=ITT+1
                 !                           Do I2=NNC,NNCD
                 !                             Do I1=NNA,NNAB
                 !                               Yxiao(ITT,I1,I2)=0.0d0
                 !                             Enddo
                 !                           Enddo
              Endif
           enddo
        enddo
     enddo
  enddo


  Do I=NII1,NII2
     NNA=Sumindex(I-1)+1
     Do J=NJJ1,NJJ2
        NNAB=SumINDEX(I+J)
        Do K=NKK1,NKK2
           NNC=Sumindex(k-1)+1
           Do L=NLL1,NLL2
              NNCD=SumIndex(K+L)
              call classmp2(I,J,K,L,NNA,NNC,NNAB,NNCD,nstepmp2s,nsteplength)
              !                   call class
           enddo
        enddo
     enddo
  enddo

end subroutine shellmp2


! Vertical Recursion by Xiao HE 07/07/07 version
subroutine shellmp2divcon(i33,ittsub)
  use allmod

  Implicit real*8(a-h,o-z)
  real*8 P(3),Q(3),W(3),KAB,KCD,AAtemp(3)
  Parameter(NN=13)
  real*8 FM(0:13)
  real*8 RA(3),RB(3),RC(3),RD(3)

  real*8 Qtemp(3),WQtemp(3),CDtemp,ABcom,Ptemp(3),WPtemp(3),ABtemp,CDcom,ABCDtemp
  integer II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2
  common /hrrstore/II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2

  COMMON /VRRcom/Qtemp,WQtemp,CDtemp,ABcom,Ptemp,WPtemp,ABtemp,CDcom,ABCDtemp

  COMMON /COM1/RA,RB,RC,RD

  Do M=1,3
     RA(M)=xyz(M,quick_basis%katom(II))
     RB(M)=xyz(M,quick_basis%katom(JJ))
     RC(M)=xyz(M,quick_basis%katom(KK))
     RD(M)=xyz(M,quick_basis%katom(LL))
  Enddo

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


  Do JJJ=1,quick_basis%kprim(JJ)
     Nprij=quick_basis%kstart(JJ)+JJJ-1
     Do III=1,quick_basis%kprim(II)
        Nprii=quick_basis%kstart(II)+III-1
        AB=Apri(Nprii,Nprij)
        ABtemp=0.5d0/AB
        cutoffprim1=dnmax*cutprim(Nprii,Nprij)
        Do M=1,3
           P(M)=Ppri(M,Nprii,Nprij)
           AAtemp(M)=P(M)*AB
           Ptemp(M)=P(M)-RA(M)
        Enddo
        Do LLL=1,quick_basis%kprim(LL)
           Npril=quick_basis%kstart(LL)+LLL-1
           Do KKK=1,quick_basis%kprim(KK)
              Nprik=quick_basis%kstart(KK)+KKK-1
              cutoffprim=cutoffprim1*cutprim(Nprik,Npril)
              If(cutoffprim.gt.quick_method%primLimit)then
                 CD=Apri(Nprik,Npril)
                 ABCD=AB+CD
                 ROU=AB*CD/ABCD
                 RPQ=0.0d0
                 ABCDxiao=dsqrt(ABCD) 

                 CDtemp=0.5d0/CD
                 ABcom=AB/ABCD
                 CDcom=CD/ABCD
                 ABCDtemp=0.5d0/ABCD

                 Do M=1,3
                    Q(M)=Ppri(M,Nprik,Npril)
                    W(M)=(AAtemp(M)+Q(M)*CD)/ABCD
                    XXXtemp=P(M)-Q(M)
                    RPQ=RPQ+XXXtemp*XXXtemp
                    Qtemp(M)=Q(M)-RC(M)
                    WQtemp(M)=W(M)-Q(M)
                    WPtemp(M)=W(M)-P(M)
                 Enddo
                 T=RPQ*ROU

                 call FmT(NABCD,T,FM)
                 do iitemp=0,NABCD
                    Yxiaotemp(1,1,iitemp)=FM(iitemp)/ABCDxiao
                 enddo
                 ITT=ITT+1

                 call vertical(NABCDTYPE) 

                 Do I2=NNC,NNCD
                    Do I1=NNA,NNAB
                       Yxiao(ITT,I1,I2)=Yxiaotemp(I1,I2,0)
                    Enddo
                 Enddo
              Endif
           enddo
        enddo
     enddo
  enddo


  Do I=NII1,NII2
     NNA=Sumindex(I-1)+1
     Do J=NJJ1,NJJ2
        NNAB=SumINDEX(I+J)
        Do K=NKK1,NKK2
           NNC=Sumindex(k-1)+1
           Do L=NLL1,NLL2
              NNCD=SumIndex(K+L)
              call classmp2divcon(I,J,K,L,NNA,NNC,NNAB,NNCD,i33,ittsub)
           enddo
        enddo
     enddo
  enddo

end subroutine shellmp2divcon

! Horrizontal recursion and Fock matrix builder by Xiao HE 07/07/07 version
 subroutine classmp2divcon(I,J,K,L,NNA,NNC,NNAB,NNCD,i33,ittsub)
   ! subroutine class
   use allmod

   Implicit real*8(A-H,O-Z)
   real*8 store(120,120)
   INTEGER NA(3),NB(3),NC(3),ND(3)
   real*8 P(3),Q(3),W(3),KAB,KCD
   Parameter(NN=13)
   real*8 FM(0:13)
   real*8 RA(3),RB(3),RC(3),RD(3)
   real*8 X44(129600)

   COMMON /COM1/RA,RB,RC,RD
   COMMON /COM2/AA,BB,CC,DD,AB,CD,ROU,ABCD
   COMMON /COM4/P,Q,W
   COMMON /COM5/FM

   integer II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2
   common /xiaostore/store
   common /hrrstore/II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2

   ITT=0

   Do JJJ=1,quick_basis%kprim(JJ)
      Nprij=quick_basis%kstart(JJ)+JJJ-1
      Do III=1,quick_basis%kprim(II)
         Nprii=quick_basis%kstart(II)+III-1
         X2=X0*XCoeff(Nprii,Nprij,I,J)
         cutoffprim1=dnmax*cutprim(Nprii,Nprij)
         Do LLL=1,quick_basis%kprim(LL)
            Npril=quick_basis%kstart(LL)+LLL-1
            Do KKK=1,quick_basis%kprim(KK)
               Nprik=quick_basis%kstart(KK)+KKK-1
               cutoffprim=cutoffprim1*cutprim(Nprik,Npril)
               If(cutoffprim.gt.quick_method%primLimit)then
                  ITT=ITT+1
                  X44(ITT)=X2*XCoeff(Nprik,Npril,K,L)
               endif
            enddo
         enddo
      enddo
   enddo

   Do MM2=NNC,NNCD
      Do MM1=NNA,NNAB
         Ytemp=0.0d0
         Do itemp=1,ITT
            Ytemp=Ytemp+X44(itemp)*Yxiao(itemp,MM1,MM2)
            !                           Ytemp=Ytemp+Yxiao(itemp,MM1,MM2)
         Enddo
         store(MM1,MM2)=Ytemp
      Enddo
   enddo


   NBI1=quick_basis%Qsbasis(II,I)
   NBI2=quick_basis%Qfbasis(II,I)
   NBJ1=quick_basis%Qsbasis(JJ,J)
   NBJ2=quick_basis%Qfbasis(JJ,J)
   NBK1=quick_basis%Qsbasis(KK,K)
   NBK2=quick_basis%Qfbasis(KK,K)
   NBL1=quick_basis%Qsbasis(LL,L)
   NBL2=quick_basis%Qfbasis(LL,L)

   !       IJKLtype=1000*I+100*J+10*K+L
   IJtype=10*I+J
   KLtype=10*K+L
   IJKLtype=100*IJtype+KLtype

   !*****       IF(max(I,J,K,L).eq.2.and.(J.ne.0.or.L.ne.0))IJKLtype=999
   IF((max(I,J,K,L).eq.2.and.(J.ne.0.or.L.ne.0)).or.(max(I,J,K,L).ge.3))IJKLtype=999
   !       IJKLtype=999
   !      If(J.eq.0.and.L.eq.0)then

   III1=quick_basis%ksumtype(II)+NBI1
   III2=quick_basis%ksumtype(II)+NBI2
   JJJ1=quick_basis%ksumtype(JJ)+NBJ1
   JJJ2=quick_basis%ksumtype(JJ)+NBJ2
   KKK1=quick_basis%ksumtype(KK)+NBK1
   KKK2=quick_basis%ksumtype(KK)+NBK2
   LLL1=quick_basis%ksumtype(LL)+NBL1
   LLL2=quick_basis%ksumtype(LL)+NBL2


   NII1=quick_basis%Qstart(II)
   NJJ1=quick_basis%Qstart(JJ)

   NBI1=quick_basis%Qsbasis(II,NII1)
   NBJ1=quick_basis%Qsbasis(JJ,NJJ1)

   II111=quick_basis%ksumtype(II)+NBI1
   JJ111=quick_basis%ksumtype(JJ)+NBJ1

   if(II.lt.JJ.and.KK.lt.LL)then

      Do III=III1,III2
         Do JJJ=JJJ1,JJJ2
            Do KKK=KKK1,KKK2
               Do LLL=LLL1,LLL2

                  call hrrwhole

                  KKKsub=wtospoint(ittsub,KKK)
                  LLLsub=wtospoint(ittsub,LLL) 

                  atemp=quick_qm_struct%co(KKKsub,i33)*Y
                  btemp=quick_qm_struct%co(LLLsub,i33)*Y

                  IIInew=III-II111+1
                  JJJnew=JJJ-JJ111+1

                  orbmp2i331(1,LLLsub,IIInew,JJJnew,1)=orbmp2i331(1,LLLsub,IIInew,JJJnew,1)+atemp
                  orbmp2i331(1,LLLsub,JJJnew,IIInew,2)=orbmp2i331(1,LLLsub,JJJnew,IIInew,2)+atemp
                  orbmp2i331(1,KKKsub,IIInew,JJJnew,1)=orbmp2i331(1,KKKsub,IIInew,JJJnew,1)+btemp
                  orbmp2i331(1,KKKsub,JJJnew,IIInew,2)=orbmp2i331(1,KKKsub,JJJnew,IIInew,2)+btemp

               Enddo
            Enddo
         Enddo
      Enddo

   else

      Do III=III1,III2
         !         if(max(III,JJJ1).le.JJJ2)then
         Do JJJ=max(III,JJJ1),JJJ2
            Do KKK=KKK1,KKK2
               !            if(max(KKK,LLL1).le.LLL2)then
               Do LLL=max(KKK,LLL1),LLL2

                  call hrrwhole

                  KKKsub=wtospoint(ittsub,KKK)
                  LLLsub=wtospoint(ittsub,LLL)

                  atemp=quick_qm_struct%co(KKKsub,i33)*Y
                  btemp=quick_qm_struct%co(LLLsub,i33)*Y

                  IIInew=III-II111+1
                  JJJnew=JJJ-JJ111+1             

                  !                 mp2shell(KKK)=.true.
                  !                 mp2shell(LLL)=.true.                

                  orbmp2i331(1,LLLsub,IIInew,JJJnew,1)=orbmp2i331(1,LLLsub,IIInew,JJJnew,1)+atemp
                  if(JJJ.ne.III)then
                     orbmp2i331(1,LLLsub,JJJnew,IIInew,2)=orbmp2i331(1,LLLsub,JJJnew,IIInew,2)+atemp
                  endif
                  if(KKK.ne.LLL)then
                     orbmp2i331(1,KKKsub,IIInew,JJJnew,1)=orbmp2i331(1,KKKsub,IIInew,JJJnew,1)+btemp
                     if(III.ne.JJJ)then
                        orbmp2i331(1,KKKsub,JJJnew,IIInew,2)=orbmp2i331(1,KKKsub,JJJnew,IIInew,2)+btemp
                     endif
                  endif

               Enddo
               !            endif
            Enddo
         Enddo
         !        endif
      Enddo

   Endif

 End subroutine classmp2divcon
      

! Horrizontal recursion and Fock matrix builder by Xiao HE 07/07/07 version
     subroutine classmp2(I,J,K,L,NNA,NNC,NNAB,NNCD,nstepmp2s,nsteplength)
       ! subroutine class
       use allmod

       Implicit real*8(A-H,O-Z)
       real*8 store(120,120)
       INTEGER NA(3),NB(3),NC(3),ND(3)
       real*8 P(3),Q(3),W(3),KAB,KCD
       Parameter(NN=13)
       real*8 FM(0:13)
       real*8 RA(3),RB(3),RC(3),RD(3)
       real*8 X44(129600)

       COMMON /COM1/RA,RB,RC,RD
       COMMON /COM2/AA,BB,CC,DD,AB,CD,ROU,ABCD
       COMMON /COM4/P,Q,W
       COMMON /COM5/FM

       integer II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2
       common /xiaostore/store
       common /hrrstore/II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2

       ITT=0
       Do JJJ=1,quick_basis%kprim(JJ)
          Nprij=quick_basis%kstart(JJ)+JJJ-1
          Do III=1,quick_basis%kprim(II)
             Nprii=quick_basis%kstart(II)+III-1
             X2=X0*XCoeff(Nprii,Nprij,I,J)
             cutoffprim1=dnmax*cutprim(Nprii,Nprij)
             Do LLL=1,quick_basis%kprim(LL)
                Npril=quick_basis%kstart(LL)+LLL-1
                Do KKK=1,quick_basis%kprim(KK)
                   Nprik=quick_basis%kstart(KK)+KKK-1
                   cutoffprim=cutoffprim1*cutprim(Nprik,Npril)
                   If(cutoffprim.gt.quick_method%primLimit)then
                      ITT=ITT+1
                      X44(ITT)=X2*XCoeff(Nprik,Npril,K,L)
                   endif
                enddo
             enddo
          enddo
       enddo

       Do MM2=NNC,NNCD
          Do MM1=NNA,NNAB
             Ytemp=0.0d0
             Do itemp=1,ITT
                Ytemp=Ytemp+X44(itemp)*Yxiao(itemp,MM1,MM2)
                !                           Ytemp=Ytemp+Yxiao(itemp,MM1,MM2)
             Enddo
             store(MM1,MM2)=Ytemp
          Enddo
       enddo


       NBI1=quick_basis%Qsbasis(II,I)
       NBI2=quick_basis%Qfbasis(II,I)
       NBJ1=quick_basis%Qsbasis(JJ,J)
       NBJ2=quick_basis%Qfbasis(JJ,J)
       NBK1=quick_basis%Qsbasis(KK,K)
       NBK2=quick_basis%Qfbasis(KK,K)
       NBL1=quick_basis%Qsbasis(LL,L)
       NBL2=quick_basis%Qfbasis(LL,L)

       !       IJKLtype=1000*I+100*J+10*K+L
       IJtype=10*I+J
       KLtype=10*K+L
       IJKLtype=100*IJtype+KLtype

       !*****       IF(max(I,J,K,L).eq.2.and.(J.ne.0.or.L.ne.0))IJKLtype=999
       IF((max(I,J,K,L).eq.2.and.(J.ne.0.or.L.ne.0)).or.(max(I,J,K,L).ge.3))IJKLtype=999
       !       IJKLtype=999
       !      If(J.eq.0.and.L.eq.0)then

       III1=quick_basis%ksumtype(II)+NBI1
       III2=quick_basis%ksumtype(II)+NBI2
       JJJ1=quick_basis%ksumtype(JJ)+NBJ1
       JJJ2=quick_basis%ksumtype(JJ)+NBJ2
       KKK1=quick_basis%ksumtype(KK)+NBK1
       KKK2=quick_basis%ksumtype(KK)+NBK2
       LLL1=quick_basis%ksumtype(LL)+NBL1
       LLL2=quick_basis%ksumtype(LL)+NBL2


       NII1=quick_basis%Qstart(II)
       NJJ1=quick_basis%Qstart(JJ)

       NBI1=quick_basis%Qsbasis(II,NII1)
       NBJ1=quick_basis%Qsbasis(JJ,NJJ1)

       II111=quick_basis%ksumtype(II)+NBI1
       JJ111=quick_basis%ksumtype(JJ)+NBJ1

       if(II.lt.JJ.and.KK.lt.LL)then

          Do III=III1,III2
             Do JJJ=JJJ1,JJJ2
                Do KKK=KKK1,KKK2
                   Do LLL=LLL1,LLL2

                      call hrrwhole

                      do i3mp2=1,nsteplength
                         i3mp2new=nstepmp2s+i3mp2-1
                         atemp=quick_qm_struct%co(KKK,i3mp2new)*Y
                         btemp=quick_qm_struct%co(LLL,i3mp2new)*Y

                         IIInew=III-II111+1
                         JJJnew=JJJ-JJ111+1

                         orbmp2i331(i3mp2,LLL,IIInew,JJJnew,1)= &
                              orbmp2i331(i3mp2,LLL,IIInew,JJJnew,1)+atemp
                         orbmp2i331(i3mp2,LLL,JJJnew,IIInew,2)= &
                              orbmp2i331(i3mp2,LLL,JJJnew,IIInew,2)+atemp
                         orbmp2i331(i3mp2,KKK,IIInew,JJJnew,1)= &
                              orbmp2i331(i3mp2,KKK,IIInew,JJJnew,1)+btemp
                         orbmp2i331(i3mp2,KKK,JJJnew,IIInew,2)= &
                              orbmp2i331(i3mp2,KKK,JJJnew,IIInew,2)+btemp
                      enddo

                   Enddo
                Enddo
             Enddo
          Enddo

       else

          Do III=III1,III2
             !         if(max(III,JJJ1).le.JJJ2)then
             Do JJJ=max(III,JJJ1),JJJ2
                Do KKK=KKK1,KKK2
                   !            if(max(KKK,LLL1).le.LLL2)then
                   Do LLL=max(KKK,LLL1),LLL2

                      call hrrwhole

                      do i3mp2=1,nsteplength
                         i3mp2new=nstepmp2s+i3mp2-1
                         atemp=quick_qm_struct%co(KKK,i3mp2new)*Y
                         btemp=quick_qm_struct%co(LLL,i3mp2new)*Y

                         IIInew=III-II111+1
                         JJJnew=JJJ-JJ111+1

                         orbmp2i331(i3mp2,LLL,IIInew,JJJnew,1)= &
                              orbmp2i331(i3mp2,LLL,IIInew,JJJnew,1)+atemp
                         if(JJJ.ne.III)then
                            orbmp2i331(i3mp2,LLL,JJJnew,IIInew,2)= &
                                 orbmp2i331(i3mp2,LLL,JJJnew,IIInew,2)+atemp
                         endif
                         if(KKK.ne.LLL)then
                            orbmp2i331(i3mp2,KKK,IIInew,JJJnew,1)= &
                                 orbmp2i331(i3mp2,KKK,IIInew,JJJnew,1)+btemp
                            if(III.ne.JJJ)then
                               orbmp2i331(i3mp2,KKK,JJJnew,IIInew,2)= &
                                    orbmp2i331(i3mp2,KKK,JJJnew,IIInew,2)+btemp
                            endif
                         endif

                      enddo

                   Enddo
                   !            endif
                Enddo
             Enddo
             !        endif
          Enddo

       Endif

     End subroutine classmp2

! Ed Brothers. October 23, 2001
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

!    subroutine attrashell(a,b,i,j,k,ii,jj,kk,Ax,Ay,Az, &
!    Bx,By,Bz,Cx,Cy,Cz,Z)
    subroutine attrashellopt(IIsh,JJsh)
    use allmod
!    use xiaoconstants
    implicit double precision(a-h,o-z)
    dimension aux(0:20)
           real*8 AA(3),BB(3),CC(3),PP(3)
           common /xiaoattra/attra,aux,AA,BB,CC,PP,g

    real*8 RA(3),RB(3),RP(3)

! Variables needed later:
!    pi=3.1415926535897932385

   Ax=xyz(1,quick_basis%katom(IIsh))
   Ay=xyz(2,quick_basis%katom(IIsh))
   Az=xyz(3,quick_basis%katom(IIsh))

   Bx=xyz(1,quick_basis%katom(JJsh))
   By=xyz(2,quick_basis%katom(JJsh))
   Bz=xyz(3,quick_basis%katom(JJsh))

!   Cx=sumx
!   Cy=sumy
!   Cz=sumz

! The purpose of this subroutine is to calculate the nuclear attraction
! of an electron  distributed between gtfs with orbital exponents a
! and b on A and B with angular momentums defined by i,j,k (a's x, y
! and z exponents, respectively) and ii,jj,k and kk on B with the core at
! (Cx,Cy,Cz) with charge Z. m is the "order" of the integral which
! arises from the recusion relationship.

! The this is taken from the recursive relation found in Obara and Saika,
! J. Chem. Phys. 84 (7) 1986, 3963.

! The first step is generating all the necessary auxillary integrals.
! These are (0|1/rc|0)^(m) = 2 Sqrt (g/Pi) (0||0) Fm(g(Rpc)^2)
! The values of m range from 0 to i+j+k+ii+jj+kk.
  
  NII2=quick_basis%Qfinal(IIsh)
  NJJ2=quick_basis%Qfinal(JJsh)
  Maxm=NII2+NJJ2+1+1

  DO ips=1,quick_basis%kprim(IIsh)
   a=gcexpo(ips,quick_basis%ksumtype(IIsh))
     Do jps=1,quick_basis%kprim(JJsh)
       b=gcexpo(jps,quick_basis%ksumtype(JJsh))       

    g = a+b
    Px = (a*Ax + b*Bx)/g
    Py = (a*Ay + b*By)/g
    Pz = (a*Az + b*Bz)/g

    constant = overlap(a,b,0,0,0,0,0,0,Ax,Ay,Az,Bx,By,Bz) &
    * 2.d0 * sqrt(g/Pi)

      Do iatom=1,natom
       if(quick_basis%katom(IIsh).eq.iatom.and.quick_basis%katom(JJsh).eq.iatom)then
         continue
       else
        Cx=xyz(1,iatom)
        Cy=xyz(2,iatom)
        Cz=xyz(3,iatom)
        Z=-1.0d0*quick_molspec%chg(iatom)
        
    PCsquare = (Px-Cx)**2 + (Py -Cy)**2 + (Pz -Cz)**2

    U = g* PCsquare
!    Maxm = i+j+k+ii+jj+kk
    call FmT(Maxm,U,aux)
    DO L = 0,maxm
        aux(L) = aux(L)*constant*Z
        attraxiao(1,1,L)=aux(L)
    ENDDO
 
    DO L = 0,maxm-1
        attraxiaoopt(1,1,1,L)=2.0d0*g*(Px-Cx)*aux(L+1)
        attraxiaoopt(2,1,1,L)=2.0d0*g*(Py-Cy)*aux(L+1)
        attraxiaoopt(3,1,1,L)=2.0d0*g*(Pz-Cz)*aux(L+1)
    ENDDO


! At this point all the auxillary integrals have been calculated.
! It is now time to decompase the attraction integral to it's
! auxillary integrals through the recursion scheme.  To do this we use
! a recursive function.

!    attraction = attrecurse(i,j,k,ii,jj,kk,0,aux,Ax,Ay,Az,Bx,By,Bz, &
!    Cx,Cy,Cz,Px,Py,Pz,g)
    NIJ1=10*NII2+NJJ2

    call nuclearattraopt(ips,jps,IIsh,JJsh,NIJ1,Ax,Ay,Az,Bx,By,Bz, &
    Cx,Cy,Cz,Px,Py,Pz,iatom)

     endif   

    enddo

   enddo
  enddo

! Xiao HE remember to multiply Z   01/12/2008
!    attraction = attraction*(-1.d0)* Z
    return
    end



!!!!!!!!BE careful of IJKLtype and Fmt.f!!!!!!!That's the difference of 6d and 10f
! vertical and hrr for gradient 10/01/2007

!Be careful of store(1,1),IASTART,size of FM, STORE, coefangxiaoL(4),coefangxiaoR(4)!!!!
! Vertical Recursion by Xiao HE 07/07/07 version
    subroutine shellopt
      use allmod

      Implicit real*8(a-h,o-z)
      real*8 P(3),Q(3),W(3),KAB,KCD
      Parameter(NN=14)
      real*8 FM(0:14)
      real*8 RA(3),RB(3),RC(3),RD(3)

      real*8 Qtemp(3),WQtemp(3),CDtemp,ABcom,Ptemp(3),WPtemp(3),ABtemp,CDcom,ABCDtemp
      integer II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2,NNABfirst,NNCDfirst
      common /hrrstore/II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2

      COMMON /VRRcom/Qtemp,WQtemp,CDtemp,ABcom,Ptemp,WPtemp,ABtemp,CDcom,ABCDtemp

      COMMON /COM1/RA,RB,RC,RD

      !    logical same
      !    same = .false.

      ! print*,II,JJ,KK,LL

      Do M=1,3
         RA(M)=xyz(M,quick_basis%katom(II))
         RB(M)=xyz(M,quick_basis%katom(JJ))
         RC(M)=xyz(M,quick_basis%katom(KK))
         RD(M)=xyz(M,quick_basis%katom(LL))
      Enddo

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

      NNABfirst=sumindex(NII2+NJJ2+1)
      NNCDfirst=sumindex(NKK2+NLL2+1)

      !     print*,'NNCDfirst=',KK,LL,NKK1,NKK2,NLL1,NLL2,NNCD,NNCDfirst

      !      NNA=Sumindex(NII1-1)+1

      !            NNC=Sumindex(NKK1-1)+1

      !       NNA=1
      !       NNC=1

      NNA=Sumindex(NII1-2)+1

      NNC=Sumindex(NKK1-2)+1

      NABCD=NII2+NJJ2+NKK2+NLL2

      ! For first derivative of nuclui motion, the total angular momentum is raised by 1
      NABCD=NABCD+1+1

      !print*,'NABCD=',NABCD

      ITT=0
      Do JJJ=1,quick_basis%kprim(JJ)
         Nprij=quick_basis%kstart(JJ)+JJJ-1
         Do III=1,quick_basis%kprim(II)
            Nprii=quick_basis%kstart(II)+III-1
            AB=Apri(Nprii,Nprij)
            ABtemp=0.5d0/AB
            cutoffprim1=dnmax*cutprim(Nprii,Nprij)
            Do M=1,3
               P(M)=Ppri(M,Nprii,Nprij)
               Ptemp(M)=P(M)-RA(M)
            Enddo
            !            KAB=Kpri(Nprii,Nprij)
            Do LLL=1,quick_basis%kprim(LL)
               Npril=quick_basis%kstart(LL)+LLL-1
               Do KKK=1,quick_basis%kprim(KK)
                  Nprik=quick_basis%kstart(KK)+KKK-1
                  cutoffprim=cutoffprim1*cutprim(Nprik,Npril)
                  If(cutoffprim.gt.quick_method%gradCutoff)then
                     CD=Apri(Nprik,Npril)
                     ABCD=AB+CD
                     ROU=AB*CD/ABCD
                     RPQ=0.0d0
                     ABCDxiao=dsqrt(ABCD) 

                     CDtemp=0.5d0/CD
                     ABcom=AB/ABCD
                     CDcom=CD/ABCD
                     ABCDtemp=0.5d0/ABCD

                     Do M=1,3
                        Q(M)=Ppri(M,Nprik,Npril)
                        W(M)=(P(M)*AB+Q(M)*CD)/ABCD
                        XXXtemp=P(M)-Q(M)
                        RPQ=RPQ+XXXtemp*XXXtemp
                        Qtemp(M)=Q(M)-RC(M)
                        WQtemp(M)=W(M)-Q(M)
                        WPtemp(M)=W(M)-P(M)
                     Enddo
                     !                         KCD=Kpri(Nprik,Npril)

                     T=RPQ*ROU

                     call FmT(NABCD,T,FM)
                     do iitemp=0,NABCD
                        Yxiaotemp(1,1,iitemp)=FM(iitemp)/ABCDxiao
                     enddo

                     ITT=ITT+1

                     call vertical(NABCDTYPE+11)

                     !                           if(NABCDTYPE.eq.44)print*,'xiao',NABCD,FM

                     Do I2=NNC,NNCDfirst
                        Do I1=NNA,NNABfirst
                           Yxiao(ITT,I1,I2)=Yxiaotemp(I1,I2,0)
                        Enddo
                     Enddo

                  Endif
               enddo
            enddo
         enddo
      enddo

      ! NNA=1
      ! NNC=1

      Do I=NII1,NII2
         NNA=Sumindex(I-2)+1
         Do J=NJJ1,NJJ2
            NNAB=SumINDEX(I+J)
            ! change for first derivative
            NNABfirst=SumINDEX(I+J+1)
            Do K=NKK1,NKK2
               NNC=Sumindex(k-2)+1
               Do L=NLL1,NLL2
                  NNCD=SumIndex(K+L)
                  NNCDfirst=SumIndex(K+L+1)
                  call classopt(I,J,K,L,NNA,NNC,NNAB,NNCD,NNABfirst,NNCDfirst)
                  !                   call class
               enddo
            enddo
         enddo
      enddo

    end subroutine shellopt

    ! Horrizontal recursion and Fock matrix builder by Xiao HE 07/07/07 version
    subroutine classopt(I,J,K,L,NNA,NNC,NNAB,NNCD,NNABfirst,NNCDfirst)
      ! subroutine class
      use allmod

      Implicit real*8(A-H,O-Z)
      real*8 store(120,120)

      real*8 storeaa(120,120)
      real*8 storebb(120,120)
      real*8 storecc(120,120) 
      real*8 storedd(120,120)

      INTEGER NA(3),NB(3),NC(3),ND(3)
      real*8 P(3),Q(3),W(3),KAB,KCD
      Parameter(NN=14)
      real*8 FM(0:13)
      real*8 RA(3),RB(3),RC(3),RD(3)
      real*8 X44(129600)

      real*8 X44aa(1296)
      real*8 X44bb(1296)
      real*8 X44cc(1296)
      real*8 X44dd(1296)
      real*8 AA,BB,CC,DD

      COMMON /COM1/RA,RB,RC,RD
      COMMON /COM2/AA,BB,CC,DD,AB,CD,ROU,ABCD
      COMMON /COM4/P,Q,W
      COMMON /COM5/FM

      integer II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2
      common /xiaostore/store
      common /xiaostoreopt/storeaa,storebb,storecc,storedd
      common /hrrstore/II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2

      ITT=0
      Do JJJ=1,quick_basis%kprim(JJ)
         Nprij=quick_basis%kstart(JJ)+JJJ-1
         BB=gcexpo(JJJ,quick_basis%ksumtype(JJ))
         Do III=1,quick_basis%kprim(II)
            Nprii=quick_basis%kstart(II)+III-1
            AA=gcexpo(III,quick_basis%ksumtype(II))

            X2=X0*XCoeff(Nprii,Nprij,I,J)
            cutoffprim1=dnmax*cutprim(Nprii,Nprij)
            Do LLL=1,quick_basis%kprim(LL)
               Npril=quick_basis%kstart(LL)+LLL-1
               DD=gcexpo(LLL,quick_basis%ksumtype(LL))
               Do KKK=1,quick_basis%kprim(KK)
                  Nprik=quick_basis%kstart(KK)+KKK-1
                  CC=gcexpo(KKK,quick_basis%ksumtype(KK))
                  cutoffprim=cutoffprim1*cutprim(Nprik,Npril)
                  If(cutoffprim.gt.quick_method%gradCutoff)then
                     ITT=ITT+1
                     X44(ITT)=X2*XCoeff(Nprik,Npril,K,L)
                     X44AA(ITT)=X44(ITT)*AA*2.0d0
                     X44BB(ITT)=X44(ITT)*BB*2.0d0
                     X44CC(ITT)=X44(ITT)*CC*2.0d0
                     !                       X44DD(ITT)=X44(ITT)*DD*2.0d0
                  endif
               enddo
            enddo
         enddo
      enddo

      !       do III=1,quick_basis%kprim(II)
      !         AA=gcexpo(III,quick_basis%ksumtype(II))


      Do MM2=NNC,NNCD
         Do MM1=NNA,NNAB
            Ytemp=0.0d0
            !                           YtempAA=0.0d0
            !                           YtempBB=0.0d0
            !                           YtempCC=0.0d0
            Do itemp=1,ITT
               Ytemp=Ytemp+X44(itemp)*Yxiao(itemp,MM1,MM2)
               !                           YtempAA=YtempAA+X44AA(itemp)*Yxiao(itemp,MM1,MM2)
               !                           YtempBB=YtempBB+X44BB(itemp)*Yxiao(itemp,MM1,MM2)
               !                           YtempCC=YtempCC+X44CC(itemp)*Yxiao(itemp,MM1,MM2)

            Enddo
            store(MM1,MM2)=Ytemp
            !                         storeAA(MM1,MM2)=YtempAA
            !                         storeBB(MM1,MM2)=YtempBB
            !                         storeCC(MM1,MM2)=YtempCC

         Enddo
      enddo


      Do MM2=NNC,NNCDfirst
         Do MM1=NNA,NNABfirst
            YtempAA=0.0d0
            YtempBB=0.0d0
            YtempCC=0.0d0
            Do itemp=1,ITT
               YtempAA=YtempAA+X44AA(itemp)*Yxiao(itemp,MM1,MM2)
               YtempBB=YtempBB+X44BB(itemp)*Yxiao(itemp,MM1,MM2)
               YtempCC=YtempCC+X44CC(itemp)*Yxiao(itemp,MM1,MM2)

               !                  if(II.eq.1.and.JJ.eq.1.and.KK.eq.13.and.LL.eq.13)then
               !                    print*,KKK-39,LLL-39,III+12,JJJ+12,itemp,X44AA(itemp),Yxiao(itemp,MM1,MM2)
               !                  endif

            Enddo
            storeAA(MM1,MM2)=YtempAA
            storeBB(MM1,MM2)=YtempBB
            storeCC(MM1,MM2)=YtempCC

         Enddo
      enddo

      NBI1=quick_basis%Qsbasis(II,I)
      NBI2=quick_basis%Qfbasis(II,I)
      NBJ1=quick_basis%Qsbasis(JJ,J)
      NBJ2=quick_basis%Qfbasis(JJ,J)
      NBK1=quick_basis%Qsbasis(KK,K)
      NBK2=quick_basis%Qfbasis(KK,K)
      NBL1=quick_basis%Qsbasis(LL,L)
      NBL2=quick_basis%Qfbasis(LL,L)

      !    iA = quick_basis%ncenter(II)
      !    iB = quick_basis%ncenter(JJ)
      !    iC = quick_basis%ncenter(KK)
      !    iD = quick_basis%ncenter(LL)

      !    same = iA.eq.iB .and. iB.eq.iC.and. iC.eq.iD

      !    print*,II,JJ,KK,LL

      !    IF (same.eq..true.) return

      !    iAstart = (iA-1)*3
      !    iBstart = (iB-1)*3
      !    iCstart = (iC-1)*3
      !    iDstart = (iD-1)*3

      Agrad1=0.d0
      Bgrad1=0.d0
      Cgrad1=0.d0
      Dgrad1=0.d0
      Agrad2=0.d0
      Bgrad2=0.d0
      Cgrad2=0.d0
      Dgrad2=0.d0
      Agrad3=0.d0
      Bgrad3=0.d0
      Cgrad3=0.d0
      Dgrad3=0.d0

      !       IJKLtype=1000*I+100*J+10*K+L
      IJtype=10*I+J
      KLtype=10*K+L
      IJKLtype=100*IJtype+KLtype

      !       IF(max(I,J,K,L).eq.2.and.(J.ne.0.or.L.ne.0))IJKLtype=999
      !       IJKLtype=999
      !      If(J.eq.0.and.L.eq.0)then

      III1=quick_basis%ksumtype(II)+NBI1
      III2=quick_basis%ksumtype(II)+NBI2
      JJJ1=quick_basis%ksumtype(JJ)+NBJ1
      JJJ2=quick_basis%ksumtype(JJ)+NBJ2
      KKK1=quick_basis%ksumtype(KK)+NBK1
      KKK2=quick_basis%ksumtype(KK)+NBK2
      LLL1=quick_basis%ksumtype(LL)+NBL1
      LLL2=quick_basis%ksumtype(LL)+NBL2

      iA = quick_basis%ncenter(III2)
      iB = quick_basis%ncenter(JJJ2)
      iC = quick_basis%ncenter(KKK2)
      iD = quick_basis%ncenter(LLL2)

      iAstart = (iA-1)*3
      iBstart = (iB-1)*3
      iCstart = (iC-1)*3
      iDstart = (iD-1)*3

      If(II.lt.JJ.and.II.lt.KK.and.KK.lt.LL)then

         !       Do III=quick_basis%ksumtype(II)+NBI1,quick_basis%ksumtype(II)+NBI2
         !         Do JJJ=quick_basis%ksumtype(JJ)+NBJ1,quick_basis%ksumtype(JJ)+NBJ2
         !            Do KKK=quick_basis%ksumtype(KK)+NBK1,quick_basis%ksumtype(KK)+NBK2
         !              Do LLL=quick_basis%ksumtype(LL)+NBL1,quick_basis%ksumtype(LL)+NBL2

         Do III=III1,III2
            Do JJJ=JJJ1,JJJ2
               Do KKK=KKK1,KKK2
                  Do LLL=LLL1,LLL2

                     !                call hrrwhole(IJKLtype,III,JJJ,KKK,LLL,Y)
                     call hrrwholeopt

                     DENSEKI=quick_qm_struct%dense(KKK,III)
                     DENSEKJ=quick_qm_struct%dense(KKK,JJJ)
                     DENSELJ=quick_qm_struct%dense(LLL,JJJ)
                     DENSELI=quick_qm_struct%dense(LLL,III)
                     DENSELK=quick_qm_struct%dense(LLL,KKK)

                     DENSEJI=quick_qm_struct%dense(JJJ,III)
                     ! Find the (ij|kl) integrals where j>i,k>i,l>k. Note that k and j
                     ! can be equal.

                     !                    O(JJJ,III) = O(JJJ,III)+2.d0*DENSELK*Y
                     !                    O(LLL,KKK) = O(LLL,KKK)+2.d0*DENSEJI*Y
                     !                    O(KKK,III) = O(KKK,III)-.5d0*DENSELJ*Y
                     !                    O(LLL,III) = O(LLL,III)-.5d0*DENSEKJ*Y
                     !                        O(JJJ,KKK) = O(JJJ,KKK)-.5d0*DENSELI*Y
                     !                        O(JJJ,LLL) = O(JJJ,LLL)-.5d0*DENSEKI*Y
                     !                        O(KKK,JJJ) = O(KKK,JJJ)-.5d0*DENSELI*Y
                     !                        O(LLL,JJJ) = O(LLL,JJJ)-.5d0*DENSEKI*Y

                     constant = (4.d0*DENSEJI*DENSELK-DENSEKI*DENSELJ &
                          -DENSELI*DENSEKJ)

                     !                    print*,'here',constant                    

                     Agrad1=Agrad1+Yaa(1)*constant
                     Agrad2=Agrad2+Yaa(2)*constant
                     Agrad3=Agrad3+Yaa(3)*constant
                     Bgrad1=Bgrad1+Ybb(1)*constant
                     Bgrad2=Bgrad2+Ybb(2)*constant
                     Bgrad3=Bgrad3+Ybb(3)*constant
                     Cgrad1=Cgrad1+Ycc(1)*constant
                     Cgrad2=Cgrad2+Ycc(2)*constant
                     Cgrad3=Cgrad3+Ycc(3)*constant

                     !                      print*,III,JJJ,KKK,LLL,Yaa(1)*constant,Ybb(1)*constant,Ycc(1)*constant, &
                     !                             Yaa(2)*constant,Ybb(2)*constant,Ycc(2)*constant, &
                     !                             Yaa(3)*constant,Ybb(3)*constant,Ycc(3)*constant

                     !                      print*,III,JJJ,KKK,LLL,Yaa(1),Ybb(1),Ycc(1), &
                     !                             Yaa(2),Ybb(2),Ycc(2), &
                     !                             Yaa(3),Ybb(3),Ycc(3)


                     !            if(III.eq.4.and.JJJ.eq.8.and.KKK.eq.8.and.LLL.eq.10)then
                     !                do ixiao=1,20
                     !                  do jxiao=1,10
                     !                    print*,ixiao,jxiao,storeAA(ixiao,jxiao),KLMN(1:3,III),KLMN(1:3,JJJ),&
                     !                           KLMN(1:3,KKK),KLMN(1:3,LLL),store(ixiao,jxiao),Yaa(1), &
                     !                           Yaa(2),Yaa(3)
                     !                  enddo
                     !                enddo
                     !            endif

                     !    iA = quick_basis%ncenter(III)        
                     !    iB = quick_basis%ncenter(JJJ)        
                     !    iC = quick_basis%ncenter(KKK)
                     !    iD = quick_basis%ncenter(LLL)
                     !
                     !    iAstart = (iA-1)*3
                     !    iBstart = (iB-1)*3
                     !    iCstart = (iC-1)*3    
                     !    iDstart = (iD-1)*3      

                  enddo
               enddo
            enddo
         enddo

      else

         !       Do III=quick_basis%ksumtype(II)+NBI1,quick_basis%ksumtype(II)+NBI2
         !         Do JJJ=max(III,quick_basis%ksumtype(JJ)+NBJ1),quick_basis%ksumtype(JJ)+NBJ2
         !            Do KKK=max(III,quick_basis%ksumtype(KK)+NBK1),quick_basis%ksumtype(KK)+NBK2
         !              Do LLL=max(KKK,quick_basis%ksumtype(LL)+NBL1),quick_basis%ksumtype(LL)+NBL2

         Do III=III1,III2                          
            Do JJJ=max(III,JJJ1),JJJ2                          
               Do KKK=max(III,KKK1),KKK2                          
                  Do LLL=max(KKK,LLL1),LLL2                          

                     !    iA = quick_basis%ncenter(III)
                     !    iB = quick_basis%ncenter(JJJ)
                     !    iC = quick_basis%ncenter(KKK)
                     !    iD = quick_basis%ncenter(LLL)
                     !
                     !    iAstart = (iA-1)*3
                     !    iBstart = (iB-1)*3
                     !    iCstart = (iC-1)*3
                     !    iDstart = (iD-1)*3

                     If(III.LT.KKK)then

                        !                call hrrwhole(IJKLtype,III,JJJ,KKK,LLL,Y)
                        call hrrwholeopt

                        If(III.lt.JJJ.and.KKK.lt.LLL)then
                           DENSEKI=quick_qm_struct%dense(KKK,III)
                           DENSEKJ=quick_qm_struct%dense(KKK,JJJ)

                           DENSELJ=quick_qm_struct%dense(LLL,JJJ)
                           DENSELI=quick_qm_struct%dense(LLL,III)
                           DENSELK=quick_qm_struct%dense(LLL,KKK)

                           DENSEJI=quick_qm_struct%dense(JJJ,III)
                           ! Find the (ij|kl) integrals where j>i,k>i,l>k. Note that k and j
                           ! can be equal.

                           !                    O(JJJ,III) = O(JJJ,III)+2.d0*DENSELK*Y
                           !                    O(LLL,KKK) = O(LLL,KKK)+2.d0*DENSEJI*Y
                           !                    O(KKK,III) = O(KKK,III)-.5d0*DENSELJ*Y
                           !                    O(LLL,III) = O(LLL,III)-.5d0*DENSEKJ*Y
                           !                        O(JJJ,KKK) = O(JJJ,KKK)-.5d0*DENSELI*Y
                           !                        O(JJJ,LLL) = O(JJJ,LLL)-.5d0*DENSEKI*Y
                           !                        O(KKK,JJJ) = O(KKK,JJJ)-.5d0*DENSELI*Y
                           !                        O(LLL,JJJ) = O(LLL,JJJ)-.5d0*DENSEKI*Y

                           constant = (4.d0*DENSEJI*DENSELK-DENSEKI*DENSELJ &
                                -DENSELI*DENSEKJ)

                           Agrad1=Agrad1+Yaa(1)*constant
                           Agrad2=Agrad2+Yaa(2)*constant
                           Agrad3=Agrad3+Yaa(3)*constant
                           Bgrad1=Bgrad1+Ybb(1)*constant
                           Bgrad2=Bgrad2+Ybb(2)*constant
                           Bgrad3=Bgrad3+Ybb(3)*constant
                           Cgrad1=Cgrad1+Ycc(1)*constant
                           Cgrad2=Cgrad2+Ycc(2)*constant
                           Cgrad3=Cgrad3+Ycc(3)*constant

                           !                      print*,III,JJJ,KKK,LLL,Yaa(1)*constant,Ybb(1)*constant,Ycc(1)*constant, &
                           !                             Yaa(2)*constant,Ybb(2)*constant,Ycc(2)*constant, &
                           !                             Yaa(3)*constant,Ybb(3)*constant,Ycc(3)*constant

                           !    ! DO all the (ii|ii) integrals.
                           !        ! Set some variables to reduce access time for some of the more
                           !        ! used quantities. (AGAIN)
                        ElseIf(III.eq.JJJ.and.KKK.eq.LLL)then
                           DENSEJI=quick_qm_struct%dense(KKK,III)
                           DENSEJJ=quick_qm_struct%dense(KKK,KKK)
                           DENSEII=quick_qm_struct%dense(III,III)

                           ! Find  all the (ii|jj) integrals.
                           !            O(III,III) = O(III,III)+DENSEJJ*Y
                           !            O(KKK,KKK) = O(KKK,KKK)+DENSEII*Y
                           !            O(KKK,III) = O(KKK,III)-.5d0*DENSEJI*Y

                           constant = (DENSEII*DENSEJJ-.5d0*DENSEJI*DENSEJI)

                           Agrad1=Agrad1+Yaa(1)*constant
                           Agrad2=Agrad2+Yaa(2)*constant
                           Agrad3=Agrad3+Yaa(3)*constant
                           Bgrad1=Bgrad1+Ybb(1)*constant
                           Bgrad2=Bgrad2+Ybb(2)*constant
                           Bgrad3=Bgrad3+Ybb(3)*constant
                           Cgrad1=Cgrad1+Ycc(1)*constant
                           Cgrad2=Cgrad2+Ycc(2)*constant
                           Cgrad3=Cgrad3+Ycc(3)*constant

                        Elseif(JJJ.eq.KKK.and.JJJ.eq.LLL)then
                           DENSEJI=quick_qm_struct%dense(JJJ,III)
                           DENSEJJ=quick_qm_struct%dense(JJJ,JJJ)

                           ! Find  all the (ij|jj) integrals.
                           !            O(JJJ,III) = O(JJJ,III)+.5d0*DENSEJJ*Y
                           !            O(JJJ,JJJ) = O(JJJ,JJJ)+DENSEJI*Y

                           constant =  DENSEJJ*DENSEJI

                           Agrad1=Agrad1+Yaa(1)*constant
                           Agrad2=Agrad2+Yaa(2)*constant
                           Agrad3=Agrad3+Yaa(3)*constant
                           Bgrad1=Bgrad1+Ybb(1)*constant
                           Bgrad2=Bgrad2+Ybb(2)*constant
                           Bgrad3=Bgrad3+Ybb(3)*constant
                           Cgrad1=Cgrad1+Ycc(1)*constant
                           Cgrad2=Cgrad2+Ycc(2)*constant
                           Cgrad3=Cgrad3+Ycc(3)*constant

                           !        ! Find  all the (ii|ij) integrals.
                           !
                           !        ! Find all the (ij|ij) integrals
                           !
                           ! Find all the (ij|ik) integrals where j>i,k>j
                        Elseif(KKK.eq.LLL.and.III.lt.JJJ.and.JJJ.ne.KKK)then
                           DENSEKI=quick_qm_struct%dense(KKK,III)
                           DENSEKJ=quick_qm_struct%dense(KKK,JJJ)
                           DENSEKK=quick_qm_struct%dense(KKK,KKK)
                           DENSEJI=quick_qm_struct%dense(JJJ,III)

                           ! Find all the (ij|kk) integrals where j>i, k>j.
                           !                O(JJJ,III) = O(JJJ,III)+DENSEKK*Y
                           !                O(KKK,KKK) = O(KKK,KKK)+2.d0*DENSEJI*Y
                           !                O(KKK,III) = O(KKK,III)-.5d0*DENSEKJ*Y
                           !                O(KKK,JJJ) = O(KKK,JJJ)-.5d0*DENSEKI*Y
                           !                O(JJJ,KKK) = O(JJJ,KKK)-.5d0*DENSEKI*Y

                           constant=(2.d0*DENSEJI*DENSEKK-DENSEKI*DENSEKJ)

                           Agrad1=Agrad1+Yaa(1)*constant
                           Agrad2=Agrad2+Yaa(2)*constant
                           Agrad3=Agrad3+Yaa(3)*constant
                           Bgrad1=Bgrad1+Ybb(1)*constant
                           Bgrad2=Bgrad2+Ybb(2)*constant
                           Bgrad3=Bgrad3+Ybb(3)*constant
                           Cgrad1=Cgrad1+Ycc(1)*constant
                           Cgrad2=Cgrad2+Ycc(2)*constant
                           Cgrad3=Cgrad3+Ycc(3)*constant

                           !            ! Find all the (ik|jj) integrals where j>i, k>j.
                        Elseif(III.eq.JJJ.and.KKK.lt.LLL)then
                           DENSEII=quick_qm_struct%dense(III,III)
                           DENSEJI=quick_qm_struct%dense(KKK,III)
                           DENSEKI=quick_qm_struct%dense(LLL,III)
                           DENSEKJ=quick_qm_struct%dense(LLL,KKK)

                           ! Find all the (ii|jk) integrals where j>i, k>j.
                           !                O(LLL,KKK) = O(LLL,KKK)+DENSEII*Y
                           !                O(III,III) = O(III,III)+2.d0*DENSEKJ*Y
                           !                O(KKK,III) = O(KKK,III)-.5d0*DENSEKI*Y
                           !                O(LLL,III) = O(LLL,III)-.5d0*DENSEJI*Y

                           constant = (2.d0*DENSEKJ*DENSEII-DENSEJI*DENSEKI)

                           Agrad1=Agrad1+Yaa(1)*constant
                           Agrad2=Agrad2+Yaa(2)*constant
                           Agrad3=Agrad3+Yaa(3)*constant
                           Bgrad1=Bgrad1+Ybb(1)*constant
                           Bgrad2=Bgrad2+Ybb(2)*constant
                           Bgrad3=Bgrad3+Ybb(3)*constant
                           Cgrad1=Cgrad1+Ycc(1)*constant
                           Cgrad2=Cgrad2+Ycc(2)*constant
                           Cgrad3=Cgrad3+Ycc(3)*constant

                        Endif

                     Else
                        If(JJJ.LE.LLL)then

                           !                call hrrwhole(IJKLtype,III,JJJ,KKK,LLL,Y)
                           call hrrwholeopt

                           if(III.eq.JJJ.and.III.eq.KKK.and.III.eq.LLL)then
                              DENSEII=quick_qm_struct%dense(III,III)

                              ! DO all the (ii|ii) integrals.
                              !        O(III,III) = O(III,III)+.5d0*DENSEII*Y

                              constant=0.0d0          

                           Elseif(III.eq.JJJ.and.III.eq.KKK.and.III.lt.LLL)then
                              DENSEJI=quick_qm_struct%dense(LLL,III)
                              DENSEII=quick_qm_struct%dense(III,III)

                              ! Find  all the (ii|ij) integrals.
                              !            O(LLL,III) = O(LLL,III)+.5d0*DENSEII*Y
                              !            O(III,III) = O(III,III)+DENSEJI*Y

                              constant= DENSEJI*DENSEII

                              Agrad1=Agrad1+Yaa(1)*constant
                              Agrad2=Agrad2+Yaa(2)*constant
                              Agrad3=Agrad3+Yaa(3)*constant
                              Bgrad1=Bgrad1+Ybb(1)*constant
                              Bgrad2=Bgrad2+Ybb(2)*constant
                              Bgrad3=Bgrad3+Ybb(3)*constant
                              Cgrad1=Cgrad1+Ycc(1)*constant
                              Cgrad2=Cgrad2+Ycc(2)*constant
                              Cgrad3=Cgrad3+Ycc(3)*constant

                           Elseif(III.eq.KKK.and.JJJ.eq.LLL.and.III.lt.JJJ)then
                              DENSEJI=quick_qm_struct%dense(JJJ,III)
                              DENSEJJ=quick_qm_struct%dense(JJJ,JJJ)
                              DENSEII=quick_qm_struct%dense(III,III)

                              ! Find all the (ij|ij) integrals
                              !            O(JJJ,III) = O(JJJ,III)+1.50*DENSEJI*Y
                              !            O(JJJ,JJJ) = O(JJJ,JJJ)-.5d0*DENSEII*Y
                              !            O(III,III) = O(III,III)-.5d0*DENSEJJ*Y

                              constant =(1.5d0*DENSEJI*DENSEJI-0.50d0*DENSEJJ*DENSEII)

                              Agrad1=Agrad1+Yaa(1)*constant
                              Agrad2=Agrad2+Yaa(2)*constant
                              Agrad3=Agrad3+Yaa(3)*constant
                              Bgrad1=Bgrad1+Ybb(1)*constant
                              Bgrad2=Bgrad2+Ybb(2)*constant
                              Bgrad3=Bgrad3+Ybb(3)*constant
                              Cgrad1=Cgrad1+Ycc(1)*constant
                              Cgrad2=Cgrad2+Ycc(2)*constant
                              Cgrad3=Cgrad3+Ycc(3)*constant

                           Elseif(III.eq.KKK.and.III.lt.JJJ.and.JJJ.lt.LLL)then
                              DENSEKI=quick_qm_struct%dense(LLL,III)
                              DENSEKJ=quick_qm_struct%dense(LLL,JJJ)
                              !                DENSEKK=quick_qm_struct%dense(LLL,LLL)
                              DENSEII=quick_qm_struct%dense(III,III)
                              DENSEJI=quick_qm_struct%dense(JJJ,III)

                              ! Find all the (ij|ik) integrals where j>i,k>j
                              !                O(JJJ,III) = O(JJJ,III)+1.5d0*DENSEKI*Y
                              !                O(LLL,III) = O(LLL,III)+1.5d0*DENSEJI*Y
                              !                O(III,III) = O(III,III)-1.d0*DENSEKJ*Y
                              !                O(LLL,JJJ) = O(LLL,JJJ)-.5d0*DENSEII*Y

                              constant = (3.0d0*DENSEJI*DENSEKI-DENSEKJ*DENSEII)

                              Agrad1=Agrad1+Yaa(1)*constant
                              Agrad2=Agrad2+Yaa(2)*constant
                              Agrad3=Agrad3+Yaa(3)*constant
                              Bgrad1=Bgrad1+Ybb(1)*constant
                              Bgrad2=Bgrad2+Ybb(2)*constant
                              Bgrad3=Bgrad3+Ybb(3)*constant
                              Cgrad1=Cgrad1+Ycc(1)*constant
                              Cgrad2=Cgrad2+Ycc(2)*constant
                              Cgrad3=Cgrad3+Ycc(3)*constant

                           Endif

                        Endif
                     Endif

                  Enddo
               Enddo
            Enddo
         Enddo
      endif

      quick_qm_struct%gradient(iASTART+1) = quick_qm_struct%gradient(iASTART+1)+ &
           AGrad1
      quick_qm_struct%gradient(iBSTART+1) = quick_qm_struct%gradient(iBSTART+1)+ &
           BGrad1
      quick_qm_struct%gradient(iCSTART+1) = quick_qm_struct%gradient(iCSTART+1)+ &
           CGrad1
      quick_qm_struct%gradient(iDSTART+1) = quick_qm_struct%gradient(iDSTART+1) &
           -AGrad1-BGrad1-CGrad1

      quick_qm_struct%gradient(iASTART+2) = quick_qm_struct%gradient(iASTART+2)+ &
           AGrad2
      quick_qm_struct%gradient(iBSTART+2) = quick_qm_struct%gradient(iBSTART+2)+ &
           BGrad2
      quick_qm_struct%gradient(iCSTART+2) = quick_qm_struct%gradient(iCSTART+2)+ &
           CGrad2
      quick_qm_struct%gradient(iDSTART+2) = quick_qm_struct%gradient(iDSTART+2) &
           -AGrad2-BGrad2-CGrad2

      quick_qm_struct%gradient(iASTART+3) = quick_qm_struct%gradient(iASTART+3)+ &
           AGrad3
      quick_qm_struct%gradient(iBSTART+3) = quick_qm_struct%gradient(iBSTART+3)+ &
           BGrad3
      quick_qm_struct%gradient(iCSTART+3) = quick_qm_struct%gradient(iCSTART+3)+ &
           CGrad3
      quick_qm_struct%gradient(iDSTART+3) = quick_qm_struct%gradient(iDSTART+3) &
           -AGrad3-BGrad3-CGrad3
      !   print*,Agrad1,BGrad1,CGrad1,AGrad2,BGrad2,CGrad2,AGrad3,BGrad3,CGrad3

      !    print*,'here' 
      !    print*,constant,Yaa1,Yaa2,Yaa3,Ybb1,Ybb2,Ybb3,Ycc1,Ycc2,Ycc3

      !      stop 
    End subroutine classopt



! Vertical Recursion by Xiao HE 07/07/07 version
 subroutine shelldft(IItemp,JJtemp,KKtemp,LLtemp)
 use allmod

 Implicit double precision(a-h,o-z)
 double precision P(3),Q(3),W(3),KAB,KCD
 Parameter(NN=13)
 double precision FM(0:13)
 double precision RA(3),RB(3),RC(3),RD(3)

 double precision Qtemp(3),WQtemp(3),CDtemp,ABcom,Ptemp(3),WPtemp(3),ABtemp,CDcom,ABCDtemp
 integer II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2
 common /hrrstore/II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2

 COMMON /VRRcom/Qtemp,WQtemp,CDtemp,ABcom,Ptemp,WPtemp,ABtemp,CDcom,ABCDtemp

 COMMON /COM1/RA,RB,RC,RD

 II=IItemp
 JJ=JJtemp
 KK=KKtemp
 LL=LLtemp

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
      do JJJ=1,quick_basis%kprim(JJ)
        Nprij=quick_basis%kstart(JJ)+JJJ-1
 do III=1,quick_basis%kprim(II)
   Nprii=quick_basis%kstart(II)+III-1
          AB=Apri(Nprii,Nprij)
          ABtemp=0.5d0/AB
          cutoffprim1=dnmax*cutprim(Nprii,Nprij)
            do M=1,3
               P(M)=Ppri(M,Nprii,Nprij)
               Ptemp(M)=P(M)-RA(M)
            enddo
!            KAB=Kpri(Nprii,Nprij)
                   do LLL=1,quick_basis%kprim(LL)
                     Npril=quick_basis%kstart(LL)+LLL-1
              do KKK=1,quick_basis%kprim(KK)
                Nprik=quick_basis%kstart(KK)+KKK-1
                       cutoffprim=cutoffprim1*cutprim(Nprik,Npril)
!                       print*,cutoffprim,quick_method%primLimit
!                       stop
                     If(cutoffprim.gt.quick_method%primLimit)then
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
                        
!                         NABCD=0 
!                         call FmT(0,T,FM)
!                         do iitemp=0,0
!                           Yxiaotemp(1,1,iitemp)=FM(iitemp)/ABCDxiao
!                         enddo
                         call FmT(NABCD,T,FM)
                         do iitemp=0,NABCD
!                           print*,iitemp,FM(iitemp),ABCDxiao,Yxiaotemp(1,1,iitemp)
                           Yxiaotemp(1,1,iitemp)=FM(iitemp)/ABCDxiao
                         enddo
!                         if(II.eq.1.and.JJ.eq.4.and.KK.eq.10.and.LL.eq.16)then
!                          print*,III,JJJ,KKK,LLL,T,NABCD,FM(0:NABCD)
!                         endif 
!                         print*,III,JJJ,KKK,LLL,FM
                           ITT=ITT+1

                           call vertical(NABCDTYPE) 
 
                           do I2=NNC,NNCD
                             do I1=NNA,NNAB
                               Yxiao(ITT,I1,I2)=Yxiaotemp(I1,I2,0)
                             enddo
                           enddo
!                           else
!!                             print*,cutoffprim
!                             ITT=ITT+1
!                           do I2=NNC,NNCD
!                             do I1=NNA,NNAB
!                               Yxiao(ITT,I1,I2)=0.0d0
!                             enddo
!                           enddo
                        endif          
                      enddo
                   enddo
         enddo
       enddo


 do I=NII1,NII2
      NNA=Sumindex(I-1)+1
   do J=NJJ1,NJJ2
     NNAB=SumINDEX(I+J)
       do K=NKK1,NKK2
            NNC=Sumindex(k-1)+1
         do L=NLL1,NLL2
           NNCD=SumIndex(K+L)
                  call classdft(I,J,K,L,NNA,NNC,NNAB,NNCD)
!                   call class
          enddo
        enddo
   enddo
 enddo
 
 end

! Horrizontal recursion and Fock matrix builder by Xiao HE 07/07/07 version
 subroutine classdft(I,J,K,L,NNA,NNC,NNAB,NNCD)
! subroutine class
 use allmod

 Implicit double precision(A-H,O-Z)
 double precision store(120,120)
 INTEGER NA(3),NB(3),NC(3),ND(3)
 double precision P(3),Q(3),W(3),KAB,KCD
 Parameter(NN=13)
 double precision FM(0:13)
 double precision RA(3),RB(3),RC(3),RD(3)
 double precision X44(1296)

 COMMON /COM1/RA,RB,RC,RD
 COMMON /COM2/AA,BB,CC,DD,AB,CD,ROU,ABCD
 COMMON /COM4/P,Q,W
 COMMON /COM5/FM

 integer II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2
 common /xiaostore/store
 common /hrrstore/II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2
 
 ITT=0
      do JJJ=1,quick_basis%kprim(JJ)
        Nprij=quick_basis%kstart(JJ)+JJJ-1
 do III=1,quick_basis%kprim(II)
    Nprii=quick_basis%kstart(II)+III-1
             X2=X0*XCoeff(Nprii,Nprij,I,J)
             cutoffprim1=dnmax*cutprim(Nprii,Nprij)
                   do LLL=1,quick_basis%kprim(LL)
                     Npril=quick_basis%kstart(LL)+LLL-1
              do KKK=1,quick_basis%kprim(KK)
                Nprik=quick_basis%kstart(KK)+KKK-1
                       cutoffprim=cutoffprim1*cutprim(Nprik,Npril)
                     If(cutoffprim.gt.quick_method%primLimit)then
                       ITT=ITT+1
                       X44(ITT)=X2*XCoeff(Nprik,Npril,K,L)
                     endif
                      enddo
                   enddo
         enddo
       enddo

                       do MM2=NNC,NNCD
                         do MM1=NNA,NNAB
                           Ytemp=0.0d0
                           do itemp=1,ITT
                           Ytemp=Ytemp+X44(itemp)*Yxiao(itemp,MM1,MM2)
!                           Ytemp=Ytemp+Yxiao(itemp,MM1,MM2)
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
   
!       IJKLtype=1000*I+100*J+10*K+L
       IJtype=10*I+J
       KLtype=10*K+L
       IJKLtype=100*IJtype+KLtype

!*****       if(max(I,J,K,L).eq.2.and.(J.ne.0.or.L.ne.0))IJKLtype=999
       if((max(I,J,K,L).eq.2.and.(J.ne.0.or.L.ne.0)).or.(max(I,J,K,L).ge.3))IJKLtype=999
!       IJKLtype=999
!      If(J.eq.0.and.L.eq.0)then

          III1=quick_basis%ksumtype(II)+NBI1
          III2=quick_basis%ksumtype(II)+NBI2
          JJJ1=quick_basis%ksumtype(JJ)+NBJ1
          JJJ2=quick_basis%ksumtype(JJ)+NBJ2
          KKK1=quick_basis%ksumtype(KK)+NBK1
          KKK2=quick_basis%ksumtype(KK)+NBK2
          LLL1=quick_basis%ksumtype(LL)+NBL1
          LLL2=quick_basis%ksumtype(LL)+NBL2
 
         If(II.lt.JJ.and.II.lt.KK.and.KK.lt.LL)then

!       do III=quick_basis%ksumtype(II)+NBI1,quick_basis%ksumtype(II)+NBI2
!         do JJJ=quick_basis%ksumtype(JJ)+NBJ1,quick_basis%ksumtype(JJ)+NBJ2
!            do KKK=quick_basis%ksumtype(KK)+NBK1,quick_basis%ksumtype(KK)+NBK2
!              do LLL=quick_basis%ksumtype(LL)+NBL1,quick_basis%ksumtype(LL)+NBL2

       do III=III1,III2
         do JJJ=JJJ1,JJJ2
            do KKK=KKK1,KKK2
              do LLL=LLL1,LLL2

!                 if(III.eq.1.and.JJJ.eq.2.and.KKK.eq.2.and.LLL.eq.30)then
!                   print*,'xiao',IJtype,KLtype,IJKLtype,KLMN(1:3,III),KLMN(1:3,JJJ),KLMN(1:3,KKK),KLMN(1:3,LLL), &
!                         store(1,18),store(1,8),store(1,2),store(1,1),RC(1)-RD(1)
!                 endif
!                 if(III.eq.1.and.JJJ.eq.4.and.KKK.eq.7.and.LLL.eq.19)then
!                   print*,'xiao',IJtype,KLtype,IJKLtype,KLMN(1:3,III),KLMN(1:3,JJJ),KLMN(1:3,KKK),KLMN(1:3,LLL), &
!                         store(3,12)*dsqrt(3.0d0),III,JJJ,KKK,LLL,'xiao11'
!                 endif

!                 if(III.eq.1.and.JJJ.eq.10.and.KKK.eq.41.and.LLL.eq.47)then
!                   print*,'xiao',IJtype,KLtype,IJKLtype,KLMN(1:3,III),KLMN(1:3,JJJ),KLMN(1:3,KKK),KLMN(1:3,LLL), &
!                         store(1,1),III,JJJ,KKK,LLL,'xiao00'
!                           do itemp=1,ITT
!                           print*,itemp,X44(itemp)*Yxiao(itemp,1,1)
!                           enddo
!                         stop
!                 endif

!                call hrrwhole(IJKLtype,III,JJJ,KKK,LLL,Y)
                call hrrwhole

                    DENSELK=quick_qm_struct%dense(LLL,KKK)
                    DENSEJI=quick_qm_struct%dense(JJJ,III)
                ! Find the (ij|kl) integrals where j>i,k>i,l>k. Note that k and j
                ! can be equal.

                    quick_qm_struct%o(JJJ,III) = quick_qm_struct%o(JJJ,III)+2.d0*DENSELK*Y
                    quick_qm_struct%o(LLL,KKK) = quick_qm_struct%o(LLL,KKK)+2.d0*DENSEJI*Y

!                      print*,III,JJJ,KKK,LLL,Y

                      enddo
                    enddo
                  enddo
                enddo

       else

!       do III=quick_basis%ksumtype(II)+NBI1,quick_basis%ksumtype(II)+NBI2
!         do JJJ=max(III,quick_basis%ksumtype(JJ)+NBJ1),quick_basis%ksumtype(JJ)+NBJ2
!            do KKK=max(III,quick_basis%ksumtype(KK)+NBK1),quick_basis%ksumtype(KK)+NBK2
!              do LLL=max(KKK,quick_basis%ksumtype(LL)+NBL1),quick_basis%ksumtype(LL)+NBL2

       do III=III1,III2                          
         do JJJ=max(III,JJJ1),JJJ2                          
            do KKK=max(III,KKK1),KKK2                          
              do LLL=max(KKK,LLL1),LLL2                          

                 If(III.LT.KKK)then

!                call hrrwhole(IJKLtype,III,JJJ,KKK,LLL,Y)
!                 if(III.eq.1.and.JJJ.eq.2.and.KKK.eq.2.and.LLL.eq.30)then
!                   print*,'xiao',store(1,18),store(1,8),store(1,2),store(1,1),RC(1)-RD(1)
!                 endif
                call hrrwhole

         If(III.lt.JJJ.and.KKK.lt.LLL)then
                    DENSELK=quick_qm_struct%dense(LLL,KKK)
                    DENSEJI=quick_qm_struct%dense(JJJ,III)
                ! Find the (ij|kl) integrals where j>i,k>i,l>k. Note that k and j
                ! can be equal.

                    quick_qm_struct%o(JJJ,III) = quick_qm_struct%o(JJJ,III)+2.d0*DENSELK*Y
                    quick_qm_struct%o(LLL,KKK) = quick_qm_struct%o(LLL,KKK)+2.d0*DENSEJI*Y

!                      print*,III,JJJ,KKK,LLL,Y

!    ! do all the (ii|ii) integrals.
!        ! Set some variables to reduce access time for some of the more
!        ! used quantities. (AGAIN)
       elseif(III.eq.JJJ.and.KKK.eq.LLL)then
            DENSEJJ=quick_qm_struct%dense(KKK,KKK)
            DENSEII=quick_qm_struct%dense(III,III)

        ! Find  all the (ii|jj) integrals.
            quick_qm_struct%o(III,III) = quick_qm_struct%o(III,III)+DENSEJJ*Y
            quick_qm_struct%o(KKK,KKK) = quick_qm_struct%o(KKK,KKK)+DENSEII*Y

       elseif(JJJ.eq.KKK.and.JJJ.eq.LLL)then
            DENSEJI=quick_qm_struct%dense(JJJ,III)
            DENSEJJ=quick_qm_struct%dense(JJJ,JJJ)

        ! Find  all the (ij|jj) integrals.
            quick_qm_struct%o(JJJ,III) = quick_qm_struct%o(JJJ,III)+DENSEJJ*Y
            quick_qm_struct%o(JJJ,JJJ) = quick_qm_struct%o(JJJ,JJJ)+2.0d0*DENSEJI*Y

!        ! Find  all the (ii|ij) integrals.
!
!        ! Find all the (ij|ij) integrals
!
            ! Find all the (ij|ik) integrals where j>i,k>j
         elseif(KKK.eq.LLL.and.III.lt.JJJ.and.JJJ.ne.KKK)then
                DENSEKK=quick_qm_struct%dense(KKK,KKK)
                DENSEJI=quick_qm_struct%dense(JJJ,III)

            ! Find all the (ij|kk) integrals where j>i, k>j.
                quick_qm_struct%o(JJJ,III) = quick_qm_struct%o(JJJ,III)+DENSEKK*Y
                quick_qm_struct%o(KKK,KKK) = quick_qm_struct%o(KKK,KKK)+2.d0*DENSEJI*Y

!            ! Find all the (ik|jj) integrals where j>i, k>j.
         elseif(III.eq.JJJ.and.KKK.lt.LLL)then
                DENSEII=quick_qm_struct%dense(III,III)
                DENSEKJ=quick_qm_struct%dense(LLL,KKK)

            ! Find all the (ii|jk) integrals where j>i, k>j.
                quick_qm_struct%o(LLL,KKK) = quick_qm_struct%o(LLL,KKK)+DENSEII*Y
                quick_qm_struct%o(III,III) = quick_qm_struct%o(III,III)+2.d0*DENSEKJ*Y

             endif

                  else
                    If(JJJ.LE.LLL)then

!                call hrrwhole(IJKLtype,III,JJJ,KKK,LLL,Y)
                call hrrwhole

      if(III.eq.JJJ.and.III.eq.KKK.and.III.eq.LLL)then
        DENSEII=quick_qm_struct%dense(III,III)

    ! do all the (ii|ii) integrals.
        quick_qm_struct%o(III,III) = quick_qm_struct%o(III,III)+DENSEII*Y
  
       elseif(III.eq.JJJ.and.III.eq.KKK.and.III.lt.LLL)then
            DENSEJI=quick_qm_struct%dense(LLL,III)
            DENSEII=quick_qm_struct%dense(III,III)

        ! Find  all the (ii|ij) integrals.
            quick_qm_struct%o(LLL,III) = quick_qm_struct%o(LLL,III)+DENSEII*Y
            quick_qm_struct%o(III,III) = quick_qm_struct%o(III,III)+2.0d0*DENSEJI*Y
        
       elseif(III.eq.KKK.and.JJJ.eq.LLL.and.III.lt.JJJ)then
            DENSEJI=quick_qm_struct%dense(JJJ,III)

        ! Find all the (ij|ij) integrals
            quick_qm_struct%o(JJJ,III) = quick_qm_struct%o(JJJ,III)+2.0d0*DENSEJI*Y

        elseif(III.eq.KKK.and.III.lt.JJJ.and.JJJ.lt.LLL)then
                DENSEKI=quick_qm_struct%dense(LLL,III)
!                DENSEKJ=quick_qm_struct%dense(LLL,JJJ)
!                DENSEKK=quick_qm_struct%dense(LLL,LLL)
                DENSEJI=quick_qm_struct%dense(JJJ,III)

            ! Find all the (ij|ik) integrals where j>i,k>j
                quick_qm_struct%o(JJJ,III) = quick_qm_struct%o(JJJ,III)+2.0d0*DENSEKI*Y
                quick_qm_struct%o(LLL,III) = quick_qm_struct%o(LLL,III)+2.0d0*DENSEJI*Y

             endif

                    endif
                  endif

             enddo
           enddo
         enddo
       enddo
     endif
    
     End
      
! Vertical Recursion by Xiao HE 07/07/07 version
 subroutine shelldftb3lyp(IItemp,JJtemp,KKtemp,LLtemp)
 use allmod

 Implicit double precision(a-h,o-z)
 double precision P(3),Q(3),W(3),KAB,KCD
 Parameter(NN=13)
 double precision FM(0:13)
 double precision RA(3),RB(3),RC(3),RD(3)

 double precision Qtemp(3),WQtemp(3),CDtemp,ABcom,Ptemp(3),WPtemp(3),ABtemp,CDcom,ABCDtemp
 integer II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2
 common /hrrstore/II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2

 COMMON /VRRcom/Qtemp,WQtemp,CDtemp,ABcom,Ptemp,WPtemp,ABtemp,CDcom,ABCDtemp

 COMMON /COM1/RA,RB,RC,RD

 II=IItemp
 JJ=JJtemp
 KK=KKtemp
 LL=LLtemp

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
      do JJJ=1,quick_basis%kprim(JJ)
        Nprij=quick_basis%kstart(JJ)+JJJ-1
 do III=1,quick_basis%kprim(II)
   Nprii=quick_basis%kstart(II)+III-1
          AB=Apri(Nprii,Nprij)
          ABtemp=0.5d0/AB
          cutoffprim1=dnmax*cutprim(Nprii,Nprij)
            do M=1,3
               P(M)=Ppri(M,Nprii,Nprij)
               Ptemp(M)=P(M)-RA(M)
            enddo
!            KAB=Kpri(Nprii,Nprij)
                   do LLL=1,quick_basis%kprim(LL)
                     Npril=quick_basis%kstart(LL)+LLL-1
              do KKK=1,quick_basis%kprim(KK)
                Nprik=quick_basis%kstart(KK)+KKK-1
                       cutoffprim=cutoffprim1*cutprim(Nprik,Npril)
!                       print*,cutoffprim,quick_method%primLimit
!                       stop
                     If(cutoffprim.gt.quick_method%primLimit)then
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
                        
!                         NABCD=0 
!                         call FmT(0,T,FM)
!                         do iitemp=0,0
!                           Yxiaotemp(1,1,iitemp)=FM(iitemp)/ABCDxiao
!                         enddo
                         call FmT(NABCD,T,FM)
                         do iitemp=0,NABCD
!                           print*,iitemp,FM(iitemp),ABCDxiao,Yxiaotemp(1,1,iitemp)
                           Yxiaotemp(1,1,iitemp)=FM(iitemp)/ABCDxiao
                         enddo
!                         if(II.eq.1.and.JJ.eq.4.and.KK.eq.10.and.LL.eq.16)then
!                          print*,III,JJJ,KKK,LLL,T,NABCD,FM(0:NABCD)
!                         endif 
!                         print*,III,JJJ,KKK,LLL,FM
                           ITT=ITT+1

                           call vertical(NABCDTYPE) 
 
                           do I2=NNC,NNCD
                             do I1=NNA,NNAB
                               Yxiao(ITT,I1,I2)=Yxiaotemp(I1,I2,0)
                             enddo
                           enddo
!                           else
!!                             print*,cutoffprim
!                             ITT=ITT+1
!                           do I2=NNC,NNCD
!                             do I1=NNA,NNAB
!                               Yxiao(ITT,I1,I2)=0.0d0
!                             enddo
!                           enddo
                        endif          
                      enddo
                   enddo
         enddo
       enddo


 do I=NII1,NII2
      NNA=Sumindex(I-1)+1
   do J=NJJ1,NJJ2
     NNAB=SumINDEX(I+J)
       do K=NKK1,NKK2
            NNC=Sumindex(k-1)+1
         do L=NLL1,NLL2
           NNCD=SumIndex(K+L)
                  call classdftb3lyp(I,J,K,L,NNA,NNC,NNAB,NNCD)
!                   call class
          enddo
        enddo
   enddo
 enddo
 
 end

! Horrizontal recursion and Fock matrix builder by Xiao HE 07/07/07 version
 subroutine classdftb3lyp(I,J,K,L,NNA,NNC,NNAB,NNCD)
! subroutine class
 use allmod

 Implicit double precision(A-H,O-Z)
 double precision store(120,120)
 INTEGER NA(3),NB(3),NC(3),ND(3)
 double precision P(3),Q(3),W(3),KAB,KCD
 Parameter(NN=13)
 double precision FM(0:13)
 double precision RA(3),RB(3),RC(3),RD(3)
 double precision X44(1296)

 COMMON /COM1/RA,RB,RC,RD
 COMMON /COM2/AA,BB,CC,DD,AB,CD,ROU,ABCD
 COMMON /COM4/P,Q,W
 COMMON /COM5/FM

 integer II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2
 common /xiaostore/store
 common /hrrstore/II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2
 
 ITT=0
      do JJJ=1,quick_basis%kprim(JJ)
        Nprij=quick_basis%kstart(JJ)+JJJ-1
 do III=1,quick_basis%kprim(II)
    Nprii=quick_basis%kstart(II)+III-1
             X2=X0*XCoeff(Nprii,Nprij,I,J)
             cutoffprim1=dnmax*cutprim(Nprii,Nprij)
                   do LLL=1,quick_basis%kprim(LL)
                     Npril=quick_basis%kstart(LL)+LLL-1
              do KKK=1,quick_basis%kprim(KK)
                Nprik=quick_basis%kstart(KK)+KKK-1
                       cutoffprim=cutoffprim1*cutprim(Nprik,Npril)
                     If(cutoffprim.gt.quick_method%primLimit)then
                       ITT=ITT+1
                       X44(ITT)=X2*XCoeff(Nprik,Npril,K,L)
                     endif
                      enddo
                   enddo
         enddo
       enddo

                       do MM2=NNC,NNCD
                         do MM1=NNA,NNAB
                           Ytemp=0.0d0
                           do itemp=1,ITT
                           Ytemp=Ytemp+X44(itemp)*Yxiao(itemp,MM1,MM2)
!                           Ytemp=Ytemp+Yxiao(itemp,MM1,MM2)
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
   
!       IJKLtype=1000*I+100*J+10*K+L
       IJtype=10*I+J
       KLtype=10*K+L
       IJKLtype=100*IJtype+KLtype

!*****       if(max(I,J,K,L).eq.2.and.(J.ne.0.or.L.ne.0))IJKLtype=999
       if((max(I,J,K,L).eq.2.and.(J.ne.0.or.L.ne.0)).or.(max(I,J,K,L).ge.3))IJKLtype=999
!       IJKLtype=999
!      If(J.eq.0.and.L.eq.0)then

          III1=quick_basis%ksumtype(II)+NBI1
          III2=quick_basis%ksumtype(II)+NBI2
          JJJ1=quick_basis%ksumtype(JJ)+NBJ1
          JJJ2=quick_basis%ksumtype(JJ)+NBJ2
          KKK1=quick_basis%ksumtype(KK)+NBK1
          KKK2=quick_basis%ksumtype(KK)+NBK2
          LLL1=quick_basis%ksumtype(LL)+NBL1
          LLL2=quick_basis%ksumtype(LL)+NBL2
 
         If(II.lt.JJ.and.II.lt.KK.and.KK.lt.LL)then

!       do III=quick_basis%ksumtype(II)+NBI1,quick_basis%ksumtype(II)+NBI2
!         do JJJ=quick_basis%ksumtype(JJ)+NBJ1,quick_basis%ksumtype(JJ)+NBJ2
!            do KKK=quick_basis%ksumtype(KK)+NBK1,quick_basis%ksumtype(KK)+NBK2
!              do LLL=quick_basis%ksumtype(LL)+NBL1,quick_basis%ksumtype(LL)+NBL2

       do III=III1,III2
         do JJJ=JJJ1,JJJ2
            do KKK=KKK1,KKK2
              do LLL=LLL1,LLL2

!                call hrrwhole(IJKLtype,III,JJJ,KKK,LLL,Y)
                call hrrwhole

                DENSEKI=quick_qm_struct%dense(KKK,III)
                DENSEKJ=quick_qm_struct%dense(KKK,JJJ)
                    DENSELJ=quick_qm_struct%dense(LLL,JJJ)
                    DENSELI=quick_qm_struct%dense(LLL,III)

                    DENSELK=quick_qm_struct%dense(LLL,KKK)
                    DENSEJI=quick_qm_struct%dense(JJJ,III)
                ! Find the (ij|kl) integrals where j>i,k>i,l>k. Note that k and j
                ! can be equal.

                    quick_qm_struct%o(JJJ,III) = quick_qm_struct%o(JJJ,III)+2.d0*DENSELK*Y
                    quick_qm_struct%o(LLL,KKK) = quick_qm_struct%o(LLL,KKK)+2.d0*DENSEJI*Y

                    quick_qm_struct%o(KKK,III) = quick_qm_struct%o(KKK,III)-.1d0*DENSELJ*Y
                    quick_qm_struct%o(LLL,III) = quick_qm_struct%o(LLL,III)-.1d0*DENSEKJ*Y
                        quick_qm_struct%o(JJJ,KKK) = quick_qm_struct%o(JJJ,KKK)-.1d0*DENSELI*Y
                        quick_qm_struct%o(JJJ,LLL) = quick_qm_struct%o(JJJ,LLL)-.1d0*DENSEKI*Y
                        quick_qm_struct%o(KKK,JJJ) = quick_qm_struct%o(KKK,JJJ)-.1d0*DENSELI*Y
                        quick_qm_struct%o(LLL,JJJ) = quick_qm_struct%o(LLL,JJJ)-.1d0*DENSEKI*Y

!                      print*,III,JJJ,KKK,LLL,Y

                      enddo
                    enddo
                  enddo
                enddo

       else

!       do III=quick_basis%ksumtype(II)+NBI1,quick_basis%ksumtype(II)+NBI2
!         do JJJ=max(III,quick_basis%ksumtype(JJ)+NBJ1),quick_basis%ksumtype(JJ)+NBJ2
!            do KKK=max(III,quick_basis%ksumtype(KK)+NBK1),quick_basis%ksumtype(KK)+NBK2
!              do LLL=max(KKK,quick_basis%ksumtype(LL)+NBL1),quick_basis%ksumtype(LL)+NBL2

       do III=III1,III2                          
         do JJJ=max(III,JJJ1),JJJ2                          
            do KKK=max(III,KKK1),KKK2                          
              do LLL=max(KKK,LLL1),LLL2                          

                 If(III.LT.KKK)then

!                call hrrwhole(IJKLtype,III,JJJ,KKK,LLL,Y)
!                 if(III.eq.1.and.JJJ.eq.2.and.KKK.eq.2.and.LLL.eq.30)then
!                   print*,'xiao',store(1,18),store(1,8),store(1,2),store(1,1),RC(1)-RD(1)
!                 endif
                call hrrwhole

         If(III.lt.JJJ.and.KKK.lt.LLL)then
                DENSEKI=quick_qm_struct%dense(KKK,III)
                DENSEKJ=quick_qm_struct%dense(KKK,JJJ)

                    DENSELJ=quick_qm_struct%dense(LLL,JJJ)
                    DENSELI=quick_qm_struct%dense(LLL,III)

                    DENSELK=quick_qm_struct%dense(LLL,KKK)
                    DENSEJI=quick_qm_struct%dense(JJJ,III)
                ! Find the (ij|kl) integrals where j>i,k>i,l>k. Note that k and j
                ! can be equal.

                    quick_qm_struct%o(JJJ,III) = quick_qm_struct%o(JJJ,III)+2.d0*DENSELK*Y
                    quick_qm_struct%o(LLL,KKK) = quick_qm_struct%o(LLL,KKK)+2.d0*DENSEJI*Y

                    quick_qm_struct%o(KKK,III) = quick_qm_struct%o(KKK,III)-.1d0*DENSELJ*Y
                    quick_qm_struct%o(LLL,III) = quick_qm_struct%o(LLL,III)-.1d0*DENSEKJ*Y
                        quick_qm_struct%o(JJJ,KKK) = quick_qm_struct%o(JJJ,KKK)-.1d0*DENSELI*Y
                        quick_qm_struct%o(JJJ,LLL) = quick_qm_struct%o(JJJ,LLL)-.1d0*DENSEKI*Y
                        quick_qm_struct%o(KKK,JJJ) = quick_qm_struct%o(KKK,JJJ)-.1d0*DENSELI*Y
                        quick_qm_struct%o(LLL,JJJ) = quick_qm_struct%o(LLL,JJJ)-.1d0*DENSEKI*Y

!                      print*,III,JJJ,KKK,LLL,Y

!    ! do all the (ii|ii) integrals.
!        ! Set some variables to reduce access time for some of the more
!        ! used quantities. (AGAIN)
       elseif(III.eq.JJJ.and.KKK.eq.LLL)then
            DENSEJI=quick_qm_struct%dense(KKK,III)
            DENSEJJ=quick_qm_struct%dense(KKK,KKK)
            DENSEII=quick_qm_struct%dense(III,III)

        ! Find  all the (ii|jj) integrals.
            quick_qm_struct%o(III,III) = quick_qm_struct%o(III,III)+DENSEJJ*Y
            quick_qm_struct%o(KKK,KKK) = quick_qm_struct%o(KKK,KKK)+DENSEII*Y
            quick_qm_struct%o(KKK,III) = quick_qm_struct%o(KKK,III)-.1d0*DENSEJI*Y

       elseif(JJJ.eq.KKK.and.JJJ.eq.LLL)then
            DENSEJI=quick_qm_struct%dense(JJJ,III)
            DENSEJJ=quick_qm_struct%dense(JJJ,JJJ)

        ! Find  all the (ij|jj) integrals.
            quick_qm_struct%o(JJJ,III) = quick_qm_struct%o(JJJ,III)+.9d0*DENSEJJ*Y
            quick_qm_struct%o(JJJ,JJJ) = quick_qm_struct%o(JJJ,JJJ)+1.8d0*DENSEJI*Y

!        ! Find  all the (ii|ij) integrals.
!
!        ! Find all the (ij|ij) integrals
!
            ! Find all the (ij|ik) integrals where j>i,k>j
         elseif(KKK.eq.LLL.and.III.lt.JJJ.and.JJJ.ne.KKK)then
                DENSEKI=quick_qm_struct%dense(KKK,III)
                DENSEKJ=quick_qm_struct%dense(KKK,JJJ)

                DENSEKK=quick_qm_struct%dense(KKK,KKK)
                DENSEJI=quick_qm_struct%dense(JJJ,III)

            ! Find all the (ij|kk) integrals where j>i, k>j.
                quick_qm_struct%o(JJJ,III) = quick_qm_struct%o(JJJ,III)+DENSEKK*Y
                quick_qm_struct%o(KKK,KKK) = quick_qm_struct%o(KKK,KKK)+2.d0*DENSEJI*Y

                quick_qm_struct%o(KKK,III) = quick_qm_struct%o(KKK,III)-.1d0*DENSEKJ*Y
                quick_qm_struct%o(KKK,JJJ) = quick_qm_struct%o(KKK,JJJ)-.1d0*DENSEKI*Y
                quick_qm_struct%o(JJJ,KKK) = quick_qm_struct%o(JJJ,KKK)-.1d0*DENSEKI*Y

!            ! Find all the (ik|jj) integrals where j>i, k>j.
         elseif(III.eq.JJJ.and.KKK.lt.LLL)then
                DENSEII=quick_qm_struct%dense(III,III)
                DENSEKJ=quick_qm_struct%dense(LLL,KKK)

                DENSEJI=quick_qm_struct%dense(KKK,III)
                DENSEKI=quick_qm_struct%dense(LLL,III)

            ! Find all the (ii|jk) integrals where j>i, k>j.
                quick_qm_struct%o(LLL,KKK) = quick_qm_struct%o(LLL,KKK)+DENSEII*Y
                quick_qm_struct%o(III,III) = quick_qm_struct%o(III,III)+2.d0*DENSEKJ*Y

                quick_qm_struct%o(KKK,III) = quick_qm_struct%o(KKK,III)-.1d0*DENSEKI*Y
                quick_qm_struct%o(LLL,III) = quick_qm_struct%o(LLL,III)-.1d0*DENSEJI*Y

             endif

                  else
                    If(JJJ.LE.LLL)then

!                call hrrwhole(IJKLtype,III,JJJ,KKK,LLL,Y)
                call hrrwhole

      if(III.eq.JJJ.and.III.eq.KKK.and.III.eq.LLL)then
        DENSEII=quick_qm_struct%dense(III,III)

    ! do all the (ii|ii) integrals.
        quick_qm_struct%o(III,III) = quick_qm_struct%o(III,III)+0.9d0*DENSEII*Y
  
       elseif(III.eq.JJJ.and.III.eq.KKK.and.III.lt.LLL)then
            DENSEJI=quick_qm_struct%dense(LLL,III)
            DENSEII=quick_qm_struct%dense(III,III)

        ! Find  all the (ii|ij) integrals.
            quick_qm_struct%o(LLL,III) = quick_qm_struct%o(LLL,III)+0.9d0*DENSEII*Y
            quick_qm_struct%o(III,III) = quick_qm_struct%o(III,III)+1.8d0*DENSEJI*Y
        
       elseif(III.eq.KKK.and.JJJ.eq.LLL.and.III.lt.JJJ)then
            DENSEJI=quick_qm_struct%dense(JJJ,III)

            DENSEJJ=quick_qm_struct%dense(JJJ,JJJ)
            DENSEII=quick_qm_struct%dense(III,III)

        ! Find all the (ij|ij) integrals
            quick_qm_struct%o(JJJ,III) = quick_qm_struct%o(JJJ,III)+1.9d0*DENSEJI*Y
            quick_qm_struct%o(JJJ,JJJ) = quick_qm_struct%o(JJJ,JJJ)-.1d0*DENSEII*Y
            quick_qm_struct%o(III,III) = quick_qm_struct%o(III,III)-.1d0*DENSEJJ*Y

        elseif(III.eq.KKK.and.III.lt.JJJ.and.JJJ.lt.LLL)then
                DENSEKI=quick_qm_struct%dense(LLL,III)
!                DENSEKJ=quick_qm_struct%dense(LLL,JJJ)
!                DENSEKK=quick_qm_struct%dense(LLL,LLL)
                DENSEJI=quick_qm_struct%dense(JJJ,III)

                DENSEKJ=quick_qm_struct%dense(LLL,JJJ)
!                DENSEKK=quick_qm_struct%dense(LLL,LLL)
                DENSEII=quick_qm_struct%dense(III,III)

            ! Find all the (ij|ik) integrals where j>i,k>j
                quick_qm_struct%o(JJJ,III) = quick_qm_struct%o(JJJ,III)+1.9d0*DENSEKI*Y
                quick_qm_struct%o(LLL,III) = quick_qm_struct%o(LLL,III)+1.9d0*DENSEJI*Y

                quick_qm_struct%o(III,III) = quick_qm_struct%o(III,III)-0.2d0*DENSEKJ*Y
                quick_qm_struct%o(LLL,JJJ) = quick_qm_struct%o(LLL,JJJ)-0.1d0*DENSEII*Y

             endif

                    endif
                  endif

             enddo
           enddo
         enddo
       enddo
     endif
    
     End
      
