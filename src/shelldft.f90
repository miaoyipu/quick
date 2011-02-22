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
   RA(M)=xyz(M,katom(II))
   RB(M)=xyz(M,katom(JJ))
   RC(M)=xyz(M,katom(KK))
   RD(M)=xyz(M,katom(LL))
 enddo

 NII1=Qstart(II)
 NII2=Qfinal(II)
 NJJ1=Qstart(JJ)
 NJJ2=Qfinal(JJ)
 NKK1=Qstart(KK)
 NKK2=Qfinal(KK)
 NLL1=Qstart(LL)
 NLL2=Qfinal(LL)

!print*,II,JJ,KK,LL,Qstart(II),Qfinal(II),Qstart(JJ),Qfinal(JJ),Qstart(KK),Qfinal(KK),Qstart(LL),Qfinal(LL)

     NNAB=(NII2+NJJ2)
     NNCD=(NKK2+NLL2)

     NABCDTYPE=NNAB*10+NNCD

     NNAB=sumindex(NNAB)
     NNCD=sumindex(NNCD)

      NNA=Sumindex(NII1-1)+1

            NNC=Sumindex(NKK1-1)+1

NABCD=NII2+NJJ2+NKK2+NLL2
ITT=0
      do JJJ=1,Kprim(JJ)
        Nprij=kstart(JJ)+JJJ-1
 do III=1,Kprim(II)
   Nprii=kstart(II)+III-1
          AB=Apri(Nprii,Nprij)
          ABtemp=0.5d0/AB
          cutoffprim1=dnmax*cutprim(Nprii,Nprij)
            do M=1,3
               P(M)=Ppri(M,Nprii,Nprij)
               Ptemp(M)=P(M)-RA(M)
            enddo
!            KAB=Kpri(Nprii,Nprij)
                   do LLL=1,Kprim(LL)
                     Npril=kstart(LL)+LLL-1
              do KKK=1,Kprim(KK)
                Nprik=kstart(KK)+KKK-1
                       cutoffprim=cutoffprim1*cutprim(Nprik,Npril)
!                       print*,cutoffprim,primlimit
!                       stop
                     If(cutoffprim.gt.primlimit)then
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
      do JJJ=1,Kprim(JJ)
        Nprij=kstart(JJ)+JJJ-1
 do III=1,Kprim(II)
    Nprii=kstart(II)+III-1
             X2=X0*XCoeff(Nprii,Nprij,I,J)
             cutoffprim1=dnmax*cutprim(Nprii,Nprij)
                   do LLL=1,Kprim(LL)
                     Npril=kstart(LL)+LLL-1
              do KKK=1,Kprim(KK)
                Nprik=kstart(KK)+KKK-1
                       cutoffprim=cutoffprim1*cutprim(Nprik,Npril)
                     If(cutoffprim.gt.primlimit)then
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


   NBI1=Qsbasis(II,I)
   NBI2=Qfbasis(II,I)
   NBJ1=Qsbasis(JJ,J)
   NBJ2=Qfbasis(JJ,J)
   NBK1=Qsbasis(KK,K)
   NBK2=Qfbasis(KK,K)
   NBL1=Qsbasis(LL,L)
   NBL2=Qfbasis(LL,L)
   
!       IJKLtype=1000*I+100*J+10*K+L
       IJtype=10*I+J
       KLtype=10*K+L
       IJKLtype=100*IJtype+KLtype

!*****       if(max(I,J,K,L).eq.2.and.(J.ne.0.or.L.ne.0))IJKLtype=999
       if((max(I,J,K,L).eq.2.and.(J.ne.0.or.L.ne.0)).or.(max(I,J,K,L).ge.3))IJKLtype=999
!       IJKLtype=999
!      If(J.eq.0.and.L.eq.0)then

          III1=Ksumtype(II)+NBI1
          III2=Ksumtype(II)+NBI2
          JJJ1=Ksumtype(JJ)+NBJ1
          JJJ2=Ksumtype(JJ)+NBJ2
          KKK1=Ksumtype(KK)+NBK1
          KKK2=Ksumtype(KK)+NBK2
          LLL1=Ksumtype(LL)+NBL1
          LLL2=Ksumtype(LL)+NBL2
 
         If(II.lt.JJ.and.II.lt.KK.and.KK.lt.LL)then

!       do III=Ksumtype(II)+NBI1,Ksumtype(II)+NBI2
!         do JJJ=Ksumtype(JJ)+NBJ1,Ksumtype(JJ)+NBJ2
!            do KKK=Ksumtype(KK)+NBK1,Ksumtype(KK)+NBK2
!              do LLL=Ksumtype(LL)+NBL1,Ksumtype(LL)+NBL2

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

                    DENSELK=DENSE(LLL,KKK)
                    DENSEJI=DENSE(JJJ,III)
                ! Find the (ij|kl) integrals where j>i,k>i,l>k. Note that k and j
                ! can be equal.

                    O(JJJ,III) = O(JJJ,III)+2.d0*DENSELK*Y
                    O(LLL,KKK) = O(LLL,KKK)+2.d0*DENSEJI*Y

!                      print*,III,JJJ,KKK,LLL,Y

                      enddo
                    enddo
                  enddo
                enddo

       else

!       do III=Ksumtype(II)+NBI1,Ksumtype(II)+NBI2
!         do JJJ=max(III,Ksumtype(JJ)+NBJ1),Ksumtype(JJ)+NBJ2
!            do KKK=max(III,Ksumtype(KK)+NBK1),Ksumtype(KK)+NBK2
!              do LLL=max(KKK,Ksumtype(LL)+NBL1),Ksumtype(LL)+NBL2

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
                    DENSELK=DENSE(LLL,KKK)
                    DENSEJI=DENSE(JJJ,III)
                ! Find the (ij|kl) integrals where j>i,k>i,l>k. Note that k and j
                ! can be equal.

                    O(JJJ,III) = O(JJJ,III)+2.d0*DENSELK*Y
                    O(LLL,KKK) = O(LLL,KKK)+2.d0*DENSEJI*Y

!                      print*,III,JJJ,KKK,LLL,Y

!    ! do all the (ii|ii) integrals.
!        ! Set some variables to reduce access time for some of the more
!        ! used quantities. (AGAIN)
       elseif(III.eq.JJJ.and.KKK.eq.LLL)then
            DENSEJJ=DENSE(KKK,KKK)
            DENSEII=DENSE(III,III)

        ! Find  all the (ii|jj) integrals.
            O(III,III) = O(III,III)+DENSEJJ*Y
            O(KKK,KKK) = O(KKK,KKK)+DENSEII*Y

       elseif(JJJ.eq.KKK.and.JJJ.eq.LLL)then
            DENSEJI=DENSE(JJJ,III)
            DENSEJJ=DENSE(JJJ,JJJ)

        ! Find  all the (ij|jj) integrals.
            O(JJJ,III) = O(JJJ,III)+DENSEJJ*Y
            O(JJJ,JJJ) = O(JJJ,JJJ)+2.0d0*DENSEJI*Y

!        ! Find  all the (ii|ij) integrals.
!
!        ! Find all the (ij|ij) integrals
!
            ! Find all the (ij|ik) integrals where j>i,k>j
         elseif(KKK.eq.LLL.and.III.lt.JJJ.and.JJJ.ne.KKK)then
                DENSEKK=DENSE(KKK,KKK)
                DENSEJI=DENSE(JJJ,III)

            ! Find all the (ij|kk) integrals where j>i, k>j.
                O(JJJ,III) = O(JJJ,III)+DENSEKK*Y
                O(KKK,KKK) = O(KKK,KKK)+2.d0*DENSEJI*Y

!            ! Find all the (ik|jj) integrals where j>i, k>j.
         elseif(III.eq.JJJ.and.KKK.lt.LLL)then
                DENSEII=DENSE(III,III)
                DENSEKJ=DENSE(LLL,KKK)

            ! Find all the (ii|jk) integrals where j>i, k>j.
                O(LLL,KKK) = O(LLL,KKK)+DENSEII*Y
                O(III,III) = O(III,III)+2.d0*DENSEKJ*Y

             endif

                  else
                    If(JJJ.LE.LLL)then

!                call hrrwhole(IJKLtype,III,JJJ,KKK,LLL,Y)
                call hrrwhole

      if(III.eq.JJJ.and.III.eq.KKK.and.III.eq.LLL)then
        DENSEII=DENSE(III,III)

    ! do all the (ii|ii) integrals.
        O(III,III) = O(III,III)+DENSEII*Y
  
       elseif(III.eq.JJJ.and.III.eq.KKK.and.III.lt.LLL)then
            DENSEJI=DENSE(LLL,III)
            DENSEII=DENSE(III,III)

        ! Find  all the (ii|ij) integrals.
            O(LLL,III) = O(LLL,III)+DENSEII*Y
            O(III,III) = O(III,III)+2.0d0*DENSEJI*Y
        
       elseif(III.eq.KKK.and.JJJ.eq.LLL.and.III.lt.JJJ)then
            DENSEJI=DENSE(JJJ,III)

        ! Find all the (ij|ij) integrals
            O(JJJ,III) = O(JJJ,III)+2.0d0*DENSEJI*Y

        elseif(III.eq.KKK.and.III.lt.JJJ.and.JJJ.lt.LLL)then
                DENSEKI=DENSE(LLL,III)
!                DENSEKJ=DENSE(LLL,JJJ)
!                DENSEKK=DENSE(LLL,LLL)
                DENSEJI=DENSE(JJJ,III)

            ! Find all the (ij|ik) integrals where j>i,k>j
                O(JJJ,III) = O(JJJ,III)+2.0d0*DENSEKI*Y
                O(LLL,III) = O(LLL,III)+2.0d0*DENSEJI*Y

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
   RA(M)=xyz(M,katom(II))
   RB(M)=xyz(M,katom(JJ))
   RC(M)=xyz(M,katom(KK))
   RD(M)=xyz(M,katom(LL))
 enddo

 NII1=Qstart(II)
 NII2=Qfinal(II)
 NJJ1=Qstart(JJ)
 NJJ2=Qfinal(JJ)
 NKK1=Qstart(KK)
 NKK2=Qfinal(KK)
 NLL1=Qstart(LL)
 NLL2=Qfinal(LL)

!print*,II,JJ,KK,LL,Qstart(II),Qfinal(II),Qstart(JJ),Qfinal(JJ),Qstart(KK),Qfinal(KK),Qstart(LL),Qfinal(LL)

     NNAB=(NII2+NJJ2)
     NNCD=(NKK2+NLL2)

     NABCDTYPE=NNAB*10+NNCD

     NNAB=sumindex(NNAB)
     NNCD=sumindex(NNCD)

      NNA=Sumindex(NII1-1)+1

            NNC=Sumindex(NKK1-1)+1

NABCD=NII2+NJJ2+NKK2+NLL2
ITT=0
      do JJJ=1,Kprim(JJ)
        Nprij=kstart(JJ)+JJJ-1
 do III=1,Kprim(II)
   Nprii=kstart(II)+III-1
          AB=Apri(Nprii,Nprij)
          ABtemp=0.5d0/AB
          cutoffprim1=dnmax*cutprim(Nprii,Nprij)
            do M=1,3
               P(M)=Ppri(M,Nprii,Nprij)
               Ptemp(M)=P(M)-RA(M)
            enddo
!            KAB=Kpri(Nprii,Nprij)
                   do LLL=1,Kprim(LL)
                     Npril=kstart(LL)+LLL-1
              do KKK=1,Kprim(KK)
                Nprik=kstart(KK)+KKK-1
                       cutoffprim=cutoffprim1*cutprim(Nprik,Npril)
!                       print*,cutoffprim,primlimit
!                       stop
                     If(cutoffprim.gt.primlimit)then
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
      do JJJ=1,Kprim(JJ)
        Nprij=kstart(JJ)+JJJ-1
 do III=1,Kprim(II)
    Nprii=kstart(II)+III-1
             X2=X0*XCoeff(Nprii,Nprij,I,J)
             cutoffprim1=dnmax*cutprim(Nprii,Nprij)
                   do LLL=1,Kprim(LL)
                     Npril=kstart(LL)+LLL-1
              do KKK=1,Kprim(KK)
                Nprik=kstart(KK)+KKK-1
                       cutoffprim=cutoffprim1*cutprim(Nprik,Npril)
                     If(cutoffprim.gt.primlimit)then
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


   NBI1=Qsbasis(II,I)
   NBI2=Qfbasis(II,I)
   NBJ1=Qsbasis(JJ,J)
   NBJ2=Qfbasis(JJ,J)
   NBK1=Qsbasis(KK,K)
   NBK2=Qfbasis(KK,K)
   NBL1=Qsbasis(LL,L)
   NBL2=Qfbasis(LL,L)
   
!       IJKLtype=1000*I+100*J+10*K+L
       IJtype=10*I+J
       KLtype=10*K+L
       IJKLtype=100*IJtype+KLtype

!*****       if(max(I,J,K,L).eq.2.and.(J.ne.0.or.L.ne.0))IJKLtype=999
       if((max(I,J,K,L).eq.2.and.(J.ne.0.or.L.ne.0)).or.(max(I,J,K,L).ge.3))IJKLtype=999
!       IJKLtype=999
!      If(J.eq.0.and.L.eq.0)then

          III1=Ksumtype(II)+NBI1
          III2=Ksumtype(II)+NBI2
          JJJ1=Ksumtype(JJ)+NBJ1
          JJJ2=Ksumtype(JJ)+NBJ2
          KKK1=Ksumtype(KK)+NBK1
          KKK2=Ksumtype(KK)+NBK2
          LLL1=Ksumtype(LL)+NBL1
          LLL2=Ksumtype(LL)+NBL2
 
         If(II.lt.JJ.and.II.lt.KK.and.KK.lt.LL)then

!       do III=Ksumtype(II)+NBI1,Ksumtype(II)+NBI2
!         do JJJ=Ksumtype(JJ)+NBJ1,Ksumtype(JJ)+NBJ2
!            do KKK=Ksumtype(KK)+NBK1,Ksumtype(KK)+NBK2
!              do LLL=Ksumtype(LL)+NBL1,Ksumtype(LL)+NBL2

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

                DENSEKI=DENSE(KKK,III)
                DENSEKJ=DENSE(KKK,JJJ)
                    DENSELJ=DENSE(LLL,JJJ)
                    DENSELI=DENSE(LLL,III)

                    DENSELK=DENSE(LLL,KKK)
                    DENSEJI=DENSE(JJJ,III)
                ! Find the (ij|kl) integrals where j>i,k>i,l>k. Note that k and j
                ! can be equal.

                    O(JJJ,III) = O(JJJ,III)+2.d0*DENSELK*Y
                    O(LLL,KKK) = O(LLL,KKK)+2.d0*DENSEJI*Y

                    O(KKK,III) = O(KKK,III)-.1d0*DENSELJ*Y
                    O(LLL,III) = O(LLL,III)-.1d0*DENSEKJ*Y
                        O(JJJ,KKK) = O(JJJ,KKK)-.1d0*DENSELI*Y
                        O(JJJ,LLL) = O(JJJ,LLL)-.1d0*DENSEKI*Y
                        O(KKK,JJJ) = O(KKK,JJJ)-.1d0*DENSELI*Y
                        O(LLL,JJJ) = O(LLL,JJJ)-.1d0*DENSEKI*Y

!                      print*,III,JJJ,KKK,LLL,Y

                      enddo
                    enddo
                  enddo
                enddo

       else

!       do III=Ksumtype(II)+NBI1,Ksumtype(II)+NBI2
!         do JJJ=max(III,Ksumtype(JJ)+NBJ1),Ksumtype(JJ)+NBJ2
!            do KKK=max(III,Ksumtype(KK)+NBK1),Ksumtype(KK)+NBK2
!              do LLL=max(KKK,Ksumtype(LL)+NBL1),Ksumtype(LL)+NBL2

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
                DENSEKI=DENSE(KKK,III)
                DENSEKJ=DENSE(KKK,JJJ)

                    DENSELJ=DENSE(LLL,JJJ)
                    DENSELI=DENSE(LLL,III)

                    DENSELK=DENSE(LLL,KKK)
                    DENSEJI=DENSE(JJJ,III)
                ! Find the (ij|kl) integrals where j>i,k>i,l>k. Note that k and j
                ! can be equal.

                    O(JJJ,III) = O(JJJ,III)+2.d0*DENSELK*Y
                    O(LLL,KKK) = O(LLL,KKK)+2.d0*DENSEJI*Y

                    O(KKK,III) = O(KKK,III)-.1d0*DENSELJ*Y
                    O(LLL,III) = O(LLL,III)-.1d0*DENSEKJ*Y
                        O(JJJ,KKK) = O(JJJ,KKK)-.1d0*DENSELI*Y
                        O(JJJ,LLL) = O(JJJ,LLL)-.1d0*DENSEKI*Y
                        O(KKK,JJJ) = O(KKK,JJJ)-.1d0*DENSELI*Y
                        O(LLL,JJJ) = O(LLL,JJJ)-.1d0*DENSEKI*Y

!                      print*,III,JJJ,KKK,LLL,Y

!    ! do all the (ii|ii) integrals.
!        ! Set some variables to reduce access time for some of the more
!        ! used quantities. (AGAIN)
       elseif(III.eq.JJJ.and.KKK.eq.LLL)then
            DENSEJI=DENSE(KKK,III)
            DENSEJJ=DENSE(KKK,KKK)
            DENSEII=DENSE(III,III)

        ! Find  all the (ii|jj) integrals.
            O(III,III) = O(III,III)+DENSEJJ*Y
            O(KKK,KKK) = O(KKK,KKK)+DENSEII*Y
            O(KKK,III) = O(KKK,III)-.1d0*DENSEJI*Y

       elseif(JJJ.eq.KKK.and.JJJ.eq.LLL)then
            DENSEJI=DENSE(JJJ,III)
            DENSEJJ=DENSE(JJJ,JJJ)

        ! Find  all the (ij|jj) integrals.
            O(JJJ,III) = O(JJJ,III)+.9d0*DENSEJJ*Y
            O(JJJ,JJJ) = O(JJJ,JJJ)+1.8d0*DENSEJI*Y

!        ! Find  all the (ii|ij) integrals.
!
!        ! Find all the (ij|ij) integrals
!
            ! Find all the (ij|ik) integrals where j>i,k>j
         elseif(KKK.eq.LLL.and.III.lt.JJJ.and.JJJ.ne.KKK)then
                DENSEKI=DENSE(KKK,III)
                DENSEKJ=DENSE(KKK,JJJ)

                DENSEKK=DENSE(KKK,KKK)
                DENSEJI=DENSE(JJJ,III)

            ! Find all the (ij|kk) integrals where j>i, k>j.
                O(JJJ,III) = O(JJJ,III)+DENSEKK*Y
                O(KKK,KKK) = O(KKK,KKK)+2.d0*DENSEJI*Y

                O(KKK,III) = O(KKK,III)-.1d0*DENSEKJ*Y
                O(KKK,JJJ) = O(KKK,JJJ)-.1d0*DENSEKI*Y
                O(JJJ,KKK) = O(JJJ,KKK)-.1d0*DENSEKI*Y

!            ! Find all the (ik|jj) integrals where j>i, k>j.
         elseif(III.eq.JJJ.and.KKK.lt.LLL)then
                DENSEII=DENSE(III,III)
                DENSEKJ=DENSE(LLL,KKK)

                DENSEJI=DENSE(KKK,III)
                DENSEKI=DENSE(LLL,III)

            ! Find all the (ii|jk) integrals where j>i, k>j.
                O(LLL,KKK) = O(LLL,KKK)+DENSEII*Y
                O(III,III) = O(III,III)+2.d0*DENSEKJ*Y

                O(KKK,III) = O(KKK,III)-.1d0*DENSEKI*Y
                O(LLL,III) = O(LLL,III)-.1d0*DENSEJI*Y

             endif

                  else
                    If(JJJ.LE.LLL)then

!                call hrrwhole(IJKLtype,III,JJJ,KKK,LLL,Y)
                call hrrwhole

      if(III.eq.JJJ.and.III.eq.KKK.and.III.eq.LLL)then
        DENSEII=DENSE(III,III)

    ! do all the (ii|ii) integrals.
        O(III,III) = O(III,III)+0.9d0*DENSEII*Y
  
       elseif(III.eq.JJJ.and.III.eq.KKK.and.III.lt.LLL)then
            DENSEJI=DENSE(LLL,III)
            DENSEII=DENSE(III,III)

        ! Find  all the (ii|ij) integrals.
            O(LLL,III) = O(LLL,III)+0.9d0*DENSEII*Y
            O(III,III) = O(III,III)+1.8d0*DENSEJI*Y
        
       elseif(III.eq.KKK.and.JJJ.eq.LLL.and.III.lt.JJJ)then
            DENSEJI=DENSE(JJJ,III)

            DENSEJJ=DENSE(JJJ,JJJ)
            DENSEII=DENSE(III,III)

        ! Find all the (ij|ij) integrals
            O(JJJ,III) = O(JJJ,III)+1.9d0*DENSEJI*Y
            O(JJJ,JJJ) = O(JJJ,JJJ)-.1d0*DENSEII*Y
            O(III,III) = O(III,III)-.1d0*DENSEJJ*Y

        elseif(III.eq.KKK.and.III.lt.JJJ.and.JJJ.lt.LLL)then
                DENSEKI=DENSE(LLL,III)
!                DENSEKJ=DENSE(LLL,JJJ)
!                DENSEKK=DENSE(LLL,LLL)
                DENSEJI=DENSE(JJJ,III)

                DENSEKJ=DENSE(LLL,JJJ)
!                DENSEKK=DENSE(LLL,LLL)
                DENSEII=DENSE(III,III)

            ! Find all the (ij|ik) integrals where j>i,k>j
                O(JJJ,III) = O(JJJ,III)+1.9d0*DENSEKI*Y
                O(LLL,III) = O(LLL,III)+1.9d0*DENSEJI*Y

                O(III,III) = O(III,III)-0.2d0*DENSEKJ*Y
                O(LLL,JJJ) = O(LLL,JJJ)-0.1d0*DENSEII*Y

             endif

                    endif
                  endif

             enddo
           enddo
         enddo
       enddo
     endif
    
     End
      
