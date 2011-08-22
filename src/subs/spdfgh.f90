!
!	spdfgh.f90
!	new_quick
!
!	Created by Yipu Miao on 3/4/11.
!	Copyright 2011 University of Florida. All rights reserved.
!
!   Vertical Recursion subroutines by hand, these parts can be optimized by MAPLE

subroutine PSSS(mtemp)
   use allmod
   Implicit double precision(a-h,o-z)
   double precision mem(35,35,0:8)
   integer CPmem(35,35,0:8)

   integer MA(3),MB(3),NA(3),NB(3),LA(3),LB(3)
   double precision RA(3),RB(3),RC(3),RD(3),P(3),Q(3),W(3)
   double precision FM(0:13)
   double precision Qtemp(3),WQtemp(3),CDtemp,ABcom,Ptemp(3),WPtemp(3),ABtemp,CDcom,ABCDtemp
   COMMON /VRRcom/Qtemp,WQtemp,CDtemp,ABcom,Ptemp,WPtemp,ABtemp,CDcom,ABCDtemp

   COMMON /COM1/RA,RB,RC,RD
   COMMON /COM2/AA,BB,CC,DD,AB,CD,ROU,ABCD
   COMMON /COM4/P,Q,W
   COMMON /COM5/FM

   do i=1,3
      Yxiaotemp(i+1,1,mtemp)=Ptemp(i)*Yxiaotemp(1,1,mtemp)+WPtemp(i)*Yxiaotemp(1,1,mtemp+1)
   enddo


end



subroutine SSPS(mtemp)
   use allmod
   Implicit double precision(a-h,o-z)
   double precision mem(35,35,0:8)
   integer CPmem(35,35,0:8)

   integer MA(3),MB(3),NA(3),NB(3),LA(3),LB(3)
   double precision RA(3),RB(3),RC(3),RD(3),P(3),Q(3),W(3)
   double precision FM(0:13)
   double precision Qtemp(3),WQtemp(3),CDtemp,ABcom,Ptemp(3),WPtemp(3),ABtemp,CDcom,ABCDtemp
   COMMON /VRRcom/Qtemp,WQtemp,CDtemp,ABcom,Ptemp,WPtemp,ABtemp,CDcom,ABCDtemp

   COMMON /COM1/RA,RB,RC,RD
   COMMON /COM2/AA,BB,CC,DD,AB,CD,ROU,ABCD
   COMMON /COM4/P,Q,W
   COMMON /COM5/FM

   do i=1,3
      Yxiaotemp(1,i+1,mtemp)=Qtemp(i)*Yxiaotemp(1,1,mtemp)+WQtemp(i)*Yxiaotemp(1,1,mtemp+1)
   enddo

end



subroutine PSPS(mtemp)
   use allmod
   Implicit double precision(a-h,o-z)
   double precision mem(35,35,0:8)
   integer CPmem(35,35,0:8)

   integer MA(3),MB(3),NA(3),NB(3),LA(3),LB(3)
   double precision RA(3),RB(3),RC(3),RD(3),P(3),Q(3),W(3)
   double precision FM(0:13)
   double precision Qtemp(3),WQtemp(3),CDtemp,ABcom,Ptemp(3),WPtemp(3),ABtemp,CDcom,ABCDtemp
   COMMON /VRRcom/Qtemp,WQtemp,CDtemp,ABcom,Ptemp,WPtemp,ABtemp,CDcom,ABCDtemp

   COMMON /COM1/RA,RB,RC,RD
   COMMON /COM2/AA,BB,CC,DD,AB,CD,ROU,ABCD
   COMMON /COM4/P,Q,W
   COMMON /COM5/FM

   do i=2,4
      do j=2,4
         Yxiaotemp(i,j,mtemp)=Ptemp(i-1)*Yxiaotemp(1,j,mtemp)+WPtemp(i-1)*Yxiaotemp(1,j,mtemp+1)
         if(i.eq.j)then
            Yxiaotemp(i,j,mtemp)=Yxiaotemp(i,j,mtemp)+ABCDtemp*Yxiaotemp(1,1,mtemp+1)
         endif
      enddo
   enddo

End


subroutine DSSS(mtemp)
   use allmod
   Implicit double precision(a-h,o-z)
   double precision mem(35,35,0:8)
   integer CPmem(35,35,0:8)

   integer MA(3),MB(3),NA(3),NB(3),LA(3),LB(3),B(3)
   double precision RA(3),RB(3),RC(3),RD(3),P(3),Q(3),W(3)
   double precision FM(0:13)
   double precision Qtemp(3),WQtemp(3),CDtemp,ABcom,Ptemp(3),WPtemp(3),ABtemp,CDcom,ABCDtemp
   COMMON /VRRcom/Qtemp,WQtemp,CDtemp,ABcom,Ptemp,WPtemp,ABtemp,CDcom,ABCDtemp

   COMMON /COM1/RA,RB,RC,RD
   COMMON /COM2/AA,BB,CC,DD,AB,CD,ROU,ABCD
   COMMON /COM4/P,Q,W
   COMMON /COM5/FM

   do i=5,10
      B(1)=Mcal(1,i)
      B(2)=Mcal(2,i)
      B(3)=Mcal(3,i)
      do j=1,3
         if(Mcal(j,i).ne.0)then
            B(j)=Mcal(j,i)-1
            itemp=trans(B(1),B(2),B(3))
            Yxiaotemp(i,1,mtemp)=Ptemp(j)*Yxiaotemp(itemp,1,mtemp)+WPtemp(j)*Yxiaotemp(itemp,1,mtemp+1)
            if(Mcal(j,i).gt.1)then
               Yxiaotemp(i,1,mtemp)=Yxiaotemp(i,1,mtemp)+ABtemp*(Yxiaotemp(1,1,mtemp) &
                     -CDcom*Yxiaotemp(1,1,mtemp+1))
            endif
            goto 111
         endif
      enddo
      111 continue
   enddo

End


subroutine SSDS(mtemp)
   use allmod
   Implicit double precision(a-h,o-z)
   double precision mem(35,35,0:8)
   integer CPmem(35,35,0:8)

   integer MA(3),MB(3),NA(3),NB(3),LA(3),LB(3),B(3)
   double precision RA(3),RB(3),RC(3),RD(3),P(3),Q(3),W(3)
   double precision FM(0:13)
   double precision Qtemp(3),WQtemp(3),CDtemp,ABcom,Ptemp(3),WPtemp(3),ABtemp,CDcom,ABCDtemp
   COMMON /VRRcom/Qtemp,WQtemp,CDtemp,ABcom,Ptemp,WPtemp,ABtemp,CDcom,ABCDtemp

   COMMON /COM1/RA,RB,RC,RD
   COMMON /COM2/AA,BB,CC,DD,AB,CD,ROU,ABCD
   COMMON /COM4/P,Q,W
   COMMON /COM5/FM

   do i=5,10
      B(1)=Mcal(1,i)
      B(2)=Mcal(2,i)
      B(3)=Mcal(3,i)
      do j=1,3
         if(Mcal(j,i).ne.0)then
            B(j)=Mcal(j,i)-1
            itemp=trans(B(1),B(2),B(3))
            Yxiaotemp(1,i,mtemp)=Qtemp(j)*Yxiaotemp(1,itemp,mtemp)+WQtemp(j)*Yxiaotemp(1,itemp,mtemp+1)
            if(Mcal(j,i).gt.1)then
               Yxiaotemp(1,i,mtemp)=Yxiaotemp(1,i,mtemp)+CDtemp*(Yxiaotemp(1,1,mtemp) &
                     -ABcom*Yxiaotemp(1,1,mtemp+1))
            endif
            goto 222
         endif
      enddo
      222 continue
   enddo

End



subroutine DSPS(mtemp)
   use allmod
   Implicit double precision(a-h,o-z)
   double precision mem(35,35,0:8)
   integer CPmem(35,35,0:8)

   integer MA(3),MB(3),NA(3),NB(3),LA(3),LB(3),B(3),Axiao(3),firstxiao,secondxiao
   double precision RA(3),RB(3),RC(3),RD(3),P(3),Q(3),W(3)
   double precision FM(0:13)
   double precision Qtemp(3),WQtemp(3),CDtemp,ABcom,Ptemp(3),WPtemp(3),ABtemp,CDcom,ABCDtemp
   COMMON /VRRcom/Qtemp,WQtemp,CDtemp,ABcom,Ptemp,WPtemp,ABtemp,CDcom,ABCDtemp

   COMMON /COM1/RA,RB,RC,RD
   COMMON /COM2/AA,BB,CC,DD,AB,CD,ROU,ABCD
   COMMON /COM4/P,Q,W
   COMMON /COM5/FM

   do i=5,10
      do jtemp=2,4
         B(1)=Mcal(1,i)
         B(2)=Mcal(2,i)
         B(3)=Mcal(3,i)
         do j=1,3
            if(Mcal(j,jtemp).ne.0)then
               Yxiaotemp(i,jtemp,mtemp)=Qtemp(j)*Yxiaotemp(i,1,mtemp)+WQtemp(j)*Yxiaotemp(i,1,mtemp+1)
               if(B(j).ne.0)then
                  B(j)=Mcal(j,i)-1
                  itemp=trans(B(1),B(2),B(3))
                  Yxiaotemp(i,jtemp,mtemp)=Yxiaotemp(i,jtemp,mtemp)+Mcal(j,i)*ABCDtemp*Yxiaotemp(itemp,1,mtemp+1)
               endif
               goto 333
            endif
         enddo
         333     continue
      enddo
   enddo

End



subroutine PSDS(mtemp)
   use allmod
   Implicit double precision(a-h,o-z)
   double precision mem(35,35,0:8)
   !     double precision fact
   integer CPmem(35,35,0:8)
   ! double precision MEM(10,10,0:4)
   ! Integer CPMEM(10,10,0:4)

   integer MA(3),MB(3),NA(3),NB(3),LA(3),LB(3),B(3),Axiao(3),firstxiao,secondxiao
   double precision RA(3),RB(3),RC(3),RD(3),P(3),Q(3),W(3)
   double precision FM(0:13)
   double precision Qtemp(3),WQtemp(3),CDtemp,ABcom,Ptemp(3),WPtemp(3),ABtemp,CDcom,ABCDtemp
   COMMON /VRRcom/Qtemp,WQtemp,CDtemp,ABcom,Ptemp,WPtemp,ABtemp,CDcom,ABCDtemp

   COMMON /COM1/RA,RB,RC,RD
   COMMON /COM2/AA,BB,CC,DD,AB,CD,ROU,ABCD
   COMMON /COM4/P,Q,W
   COMMON /COM5/FM

   !     call SSDS(mtemp)
   !     call SSDS(mtemp+1)
   !     call SSPS(mtemp+1)

   do i=5,10
      !        B(1)=Mcal(1,i)
      !        B(2)=Mcal(2,i)
      !        B(3)=Mcal(3,i)
      do jtemp=2,4
         B(1)=Mcal(1,i)
         B(2)=Mcal(2,i)
         B(3)=Mcal(3,i)
         !          Axiao(1)=Mcal(1,jtemp)
         !          Axiao(2)=Mcal(2,jtemp)
         !          Axiao(3)=Mcal(3,jtemp)
         do j=1,3
            if(Mcal(j,jtemp).ne.0)then
               Yxiaotemp(jtemp,i,mtemp)=Ptemp(j)*Yxiaotemp(1,i,mtemp)+WPtemp(j)*Yxiaotemp(1,i,mtemp+1)
               if(B(j).ne.0)then
                  B(j)=Mcal(j,i)-1
                  itemp=trans(B(1),B(2),B(3))
                  Yxiaotemp(jtemp,i,mtemp)=Yxiaotemp(jtemp,i,mtemp)+Mcal(j,i)*ABCDtemp*Yxiaotemp(1,itemp,mtemp+1)
               endif
               goto 444
            endif
         enddo
         444     continue
      enddo
   enddo

End


subroutine DSDS(mtemp)
   use allmod
   Implicit double precision(a-h,o-z)
   double precision mem(35,35,0:8)
   !     double precision fact
   integer CPmem(35,35,0:8)
   ! double precision MEM(10,10,0:4)
   ! Integer CPMEM(10,10,0:4)

   integer MA(3),MB(3),NA(3),NB(3),LA(3),LB(3),B(3),Axiao(3),firstxiao,secondxiao
   double precision RA(3),RB(3),RC(3),RD(3),P(3),Q(3),W(3)
   double precision FM(0:13)
   double precision Qtemp(3),WQtemp(3),CDtemp,ABcom,Ptemp(3),WPtemp(3),ABtemp,CDcom,ABCDtemp
   COMMON /VRRcom/Qtemp,WQtemp,CDtemp,ABcom,Ptemp,WPtemp,ABtemp,CDcom,ABCDtemp

   COMMON /COM1/RA,RB,RC,RD
   COMMON /COM2/AA,BB,CC,DD,AB,CD,ROU,ABCD
   COMMON /COM4/P,Q,W
   COMMON /COM5/FM

   !     call SSDS(mtemp)
   !     call SSDS(mtemp+1)
   !     call SSDS(mtemp+2)
   !     call SSPS(mtemp+1)
   !     call SSPS(mtemp+2)

   do i=5,10
      !        B(1)=Mcal(1,i)
      !        B(2)=Mcal(2,i)
      !        B(3)=Mcal(3,i)
      do jtemp=5,10
         B(1)=Mcal(1,i)
         B(2)=Mcal(2,i)
         B(3)=Mcal(3,i)
         Axiao(1)=Mcal(1,jtemp)
         Axiao(2)=Mcal(2,jtemp)
         Axiao(3)=Mcal(3,jtemp)
         do j=1,3
            if(Mcal(j,i).ne.0)then
               B(j)=Mcal(j,i)-1
               ixiao=trans(B(1),B(2),B(3))
               Yxiaotemp(i,jtemp,mtemp)=Ptemp(j)*Yxiaotemp(ixiao,jtemp,mtemp) &
                     +WPtemp(j)*Yxiaotemp(ixiao,jtemp,mtemp+1)
               if(Mcal(j,i).ge.2)then
                  !               B(j)=B(j)-1
                  !               itemp=trans(B(1),B(2),B(3))
                  !               Yxiaotemp(i,jtemp,mtemp)=Yxiaotemp(i,jtemp,mtemp)+ABtemp*(Yxiaotemp(itemp,jtemp,mtemp) &
                        !                                       -CDcom*Yxiaotemp(itemp,jtemp,mtemp+1))
                  Yxiaotemp(i,jtemp,mtemp)=Yxiaotemp(i,jtemp,mtemp)+ABtemp*(Yxiaotemp(1,jtemp,mtemp) &
                        -CDcom*Yxiaotemp(1,jtemp,mtemp+1))
               endif
               if(Axiao(j).ne.0)then
                  Axiao(j)=Axiao(j)-1
                  secondxiao=trans(Axiao(1),Axiao(2),Axiao(3))
                  Yxiaotemp(i,jtemp,mtemp)=Yxiaotemp(i,jtemp,mtemp) &
                        +Mcal(j,jtemp)*ABCDtemp*Yxiaotemp(ixiao,secondxiao,mtemp+1)
               endif
               goto 555
            endif
         enddo
         555     continue
      enddo
   enddo

End


subroutine FSSS(mtemp)
   use allmod
   Implicit double precision(a-h,o-z)
   double precision mem(35,35,0:8)
   !     double precision fact
   integer CPmem(35,35,0:8)
   ! double precision MEM(10,10,0:4)
   ! Integer CPMEM(10,10,0:4)

   integer MA(3),MB(3),NA(3),NB(3),LA(3),LB(3),B(3),Axiao(3),firstxiao,secondxiao
   double precision RA(3),RB(3),RC(3),RD(3),P(3),Q(3),W(3)
   double precision FM(0:13)
   double precision Qtemp(3),WQtemp(3),CDtemp,ABcom,Ptemp(3),WPtemp(3),ABtemp,CDcom,ABCDtemp
   COMMON /VRRcom/Qtemp,WQtemp,CDtemp,ABcom,Ptemp,WPtemp,ABtemp,CDcom,ABCDtemp

   COMMON /COM1/RA,RB,RC,RD
   COMMON /COM2/AA,BB,CC,DD,AB,CD,ROU,ABCD
   COMMON /COM4/P,Q,W
   COMMON /COM5/FM

   do i=11,20
      B(1)=Mcal(1,i)
      B(2)=Mcal(2,i)
      B(3)=Mcal(3,i)
      do j=1,3
         if(Mcal(j,i).ne.0)then
            B(j)=Mcal(j,i)-1
            itemp=trans(B(1),B(2),B(3))
            Yxiaotemp(i,1,mtemp)=Ptemp(j)*Yxiaotemp(itemp,1,mtemp)+WPtemp(j)*Yxiaotemp(itemp,1,mtemp+1)
            if(Mcal(j,i).gt.1)then
               B(j)=Mcal(j,i)-2
               inewtemp=trans(B(1),B(2),B(3))
               Yxiaotemp(i,1,mtemp)=Yxiaotemp(i,1,mtemp)+(B(j)+1)*ABtemp*(Yxiaotemp(inewtemp,1,mtemp) &
                     -CDcom*Yxiaotemp(inewtemp,1,mtemp+1))
            endif
            goto 111
         endif
      enddo
      111 continue
   enddo

End

subroutine SSFS(mtemp)
   use allmod
   Implicit double precision(a-h,o-z)
   double precision mem(35,35,0:8)
   !     double precision fact
   integer CPmem(35,35,0:8)
   ! double precision MEM(10,10,0:4)
   ! Integer CPMEM(10,10,0:4)

   integer MA(3),MB(3),NA(3),NB(3),LA(3),LB(3),B(3),Axiao(3),firstxiao,secondxiao
   double precision RA(3),RB(3),RC(3),RD(3),P(3),Q(3),W(3)
   double precision FM(0:13)
   double precision Qtemp(3),WQtemp(3),CDtemp,ABcom,Ptemp(3),WPtemp(3),ABtemp,CDcom,ABCDtemp
   COMMON /VRRcom/Qtemp,WQtemp,CDtemp,ABcom,Ptemp,WPtemp,ABtemp,CDcom,ABCDtemp

   COMMON /COM1/RA,RB,RC,RD
   COMMON /COM2/AA,BB,CC,DD,AB,CD,ROU,ABCD
   COMMON /COM4/P,Q,W
   COMMON /COM5/FM

   do i=11,20
      B(1)=Mcal(1,i)
      B(2)=Mcal(2,i)
      B(3)=Mcal(3,i)
      do j=1,3
         if(Mcal(j,i).ne.0)then
            B(j)=Mcal(j,i)-1
            itemp=trans(B(1),B(2),B(3))
            Yxiaotemp(1,i,mtemp)=Qtemp(j)*Yxiaotemp(1,itemp,mtemp)+WQtemp(j)*Yxiaotemp(1,itemp,mtemp+1)
            if(Mcal(j,i).gt.1)then
               B(j)=Mcal(j,i)-2
               inewtemp=trans(B(1),B(2),B(3))
               Yxiaotemp(1,i,mtemp)=Yxiaotemp(1,i,mtemp)+(B(j)+1)*CDtemp*(Yxiaotemp(1,inewtemp,mtemp) &
                     -ABcom*Yxiaotemp(1,inewtemp,mtemp+1))
            endif
            goto 111
         endif
      enddo
      111 continue
   enddo

End

subroutine GSSS(mtemp)
   use allmod
   Implicit double precision(a-h,o-z)
   double precision mem(35,35,0:8)
   !     double precision fact
   integer CPmem(35,35,0:8)
   ! double precision MEM(10,10,0:4)
   ! Integer CPMEM(10,10,0:4)

   integer MA(3),MB(3),NA(3),NB(3),LA(3),LB(3),B(3),Axiao(3),firstxiao,secondxiao
   double precision RA(3),RB(3),RC(3),RD(3),P(3),Q(3),W(3)
   double precision FM(0:13)
   double precision Qtemp(3),WQtemp(3),CDtemp,ABcom,Ptemp(3),WPtemp(3),ABtemp,CDcom,ABCDtemp
   COMMON /VRRcom/Qtemp,WQtemp,CDtemp,ABcom,Ptemp,WPtemp,ABtemp,CDcom,ABCDtemp

   COMMON /COM1/RA,RB,RC,RD
   COMMON /COM2/AA,BB,CC,DD,AB,CD,ROU,ABCD
   COMMON /COM4/P,Q,W
   COMMON /COM5/FM

   do i=21,35
      B(1)=Mcal(1,i)
      B(2)=Mcal(2,i)
      B(3)=Mcal(3,i)
      do j=1,3
         if(Mcal(j,i).ne.0)then
            B(j)=Mcal(j,i)-1
            itemp=trans(B(1),B(2),B(3))
            Yxiaotemp(i,1,mtemp)=Ptemp(j)*Yxiaotemp(itemp,1,mtemp)+WPtemp(j)*Yxiaotemp(itemp,1,mtemp+1)
            if(Mcal(j,i).gt.1)then
               B(j)=Mcal(j,i)-2
               inewtemp=trans(B(1),B(2),B(3))
               Yxiaotemp(i,1,mtemp)=Yxiaotemp(i,1,mtemp)+(B(j)+1)*ABtemp*(Yxiaotemp(inewtemp,1,mtemp) &
                     -CDcom*Yxiaotemp(inewtemp,1,mtemp+1))
            endif
            goto 111
         endif
      enddo
      111 continue
   enddo

End

subroutine SSGS(mtemp)
   use allmod
   Implicit double precision(a-h,o-z)
   double precision mem(35,35,0:8)
   !     double precision fact
   integer CPmem(35,35,0:8)
   ! double precision MEM(10,10,0:4)
   ! Integer CPMEM(10,10,0:4)

   integer MA(3),MB(3),NA(3),NB(3),LA(3),LB(3),B(3),Axiao(3),firstxiao,secondxiao
   double precision RA(3),RB(3),RC(3),RD(3),P(3),Q(3),W(3)
   double precision FM(0:13)
   double precision Qtemp(3),WQtemp(3),CDtemp,ABcom,Ptemp(3),WPtemp(3),ABtemp,CDcom,ABCDtemp
   COMMON /VRRcom/Qtemp,WQtemp,CDtemp,ABcom,Ptemp,WPtemp,ABtemp,CDcom,ABCDtemp

   COMMON /COM1/RA,RB,RC,RD
   COMMON /COM2/AA,BB,CC,DD,AB,CD,ROU,ABCD
   COMMON /COM4/P,Q,W
   COMMON /COM5/FM

   do i=21,35
      B(1)=Mcal(1,i)
      B(2)=Mcal(2,i)
      B(3)=Mcal(3,i)
      do j=1,3
         if(Mcal(j,i).ne.0)then
            B(j)=Mcal(j,i)-1
            itemp=trans(B(1),B(2),B(3))
            Yxiaotemp(1,i,mtemp)=Qtemp(j)*Yxiaotemp(1,itemp,mtemp)+WQtemp(j)*Yxiaotemp(1,itemp,mtemp+1)
            if(Mcal(j,i).gt.1)then
               B(j)=Mcal(j,i)-2
               inewtemp=trans(B(1),B(2),B(3))
               Yxiaotemp(1,i,mtemp)=Yxiaotemp(1,i,mtemp)+(B(j)+1)*CDtemp*(Yxiaotemp(1,inewtemp,mtemp) &
                     -ABcom*Yxiaotemp(1,inewtemp,mtemp+1))
            endif
            goto 111
         endif
      enddo
      111 continue
   enddo

End


subroutine FSPS(mtemp)
   use allmod
   Implicit double precision(a-h,o-z)
   double precision mem(35,35,0:8)
   !     double precision fact
   integer CPmem(35,35,0:8)
   ! double precision MEM(10,10,0:4)
   ! Integer CPMEM(10,10,0:4)

   integer MA(3),MB(3),NA(3),NB(3),LA(3),LB(3),B(3),Axiao(3),firstxiao,secondxiao
   double precision RA(3),RB(3),RC(3),RD(3),P(3),Q(3),W(3)
   double precision FM(0:13)
   double precision Qtemp(3),WQtemp(3),CDtemp,ABcom,Ptemp(3),WPtemp(3),ABtemp,CDcom,ABCDtemp
   COMMON /VRRcom/Qtemp,WQtemp,CDtemp,ABcom,Ptemp,WPtemp,ABtemp,CDcom,ABCDtemp

   COMMON /COM1/RA,RB,RC,RD
   COMMON /COM2/AA,BB,CC,DD,AB,CD,ROU,ABCD
   COMMON /COM4/P,Q,W
   COMMON /COM5/FM

   do i=11,20
      !        B(1)=Mcal(1,i)
      !        B(2)=Mcal(2,i)
      !        B(3)=Mcal(3,i)
      do jtemp=2,4
         B(1)=Mcal(1,i)
         B(2)=Mcal(2,i)
         B(3)=Mcal(3,i)
         !          Axiao(1)=Mcal(1,jtemp)
         !          Axiao(2)=Mcal(2,jtemp)
         !          Axiao(3)=Mcal(3,jtemp)
         do j=1,3
            if(Mcal(j,jtemp).ne.0)then
               Yxiaotemp(i,jtemp,mtemp)=Qtemp(j)*Yxiaotemp(i,1,mtemp)+WQtemp(j)*Yxiaotemp(i,1,mtemp+1)
               if(B(j).ne.0)then
                  B(j)=Mcal(j,i)-1
                  itemp=trans(B(1),B(2),B(3))
                  Yxiaotemp(i,jtemp,mtemp)=Yxiaotemp(i,jtemp,mtemp)+Mcal(j,i)*ABCDtemp*Yxiaotemp(itemp,1,mtemp+1)
               endif
               !             if(i.eq.5.and.jtemp.eq.3)then
               !               print*,B(1),B(2),B(3),trans(B(1),B(2),B(3)),itemp
               !               print*,'quick',Yxiaotemp(5,1,0),Yxiaotemp(5,1,1),Yxiaotemp(2,1,1),Yxiaotemp(5,3,0) &
                     !,Qtemp(2),WQtemp(2),ABCDtemp,Qtemp(2)*Yxiaotemp(5,1,0) &
                     !                               +WQtemp(2)*Yxiaotemp(5,1,1)+ABCDtemp*Yxiaotemp(2,1,1)
               !
               !             endif
               goto 333
            endif
         enddo
         333     continue
      enddo
   enddo

End

subroutine PSFS(mtemp)
   use allmod
   Implicit double precision(a-h,o-z)
   double precision mem(35,35,0:8)
   !     double precision fact
   integer CPmem(35,35,0:8)
   ! double precision MEM(10,10,0:4)
   ! Integer CPMEM(10,10,0:4)

   integer MA(3),MB(3),NA(3),NB(3),LA(3),LB(3),B(3),Axiao(3),firstxiao,secondxiao
   double precision RA(3),RB(3),RC(3),RD(3),P(3),Q(3),W(3)
   double precision FM(0:13)
   double precision Qtemp(3),WQtemp(3),CDtemp,ABcom,Ptemp(3),WPtemp(3),ABtemp,CDcom,ABCDtemp
   COMMON /VRRcom/Qtemp,WQtemp,CDtemp,ABcom,Ptemp,WPtemp,ABtemp,CDcom,ABCDtemp

   COMMON /COM1/RA,RB,RC,RD
   COMMON /COM2/AA,BB,CC,DD,AB,CD,ROU,ABCD
   COMMON /COM4/P,Q,W
   COMMON /COM5/FM

   do i=11,20
      !        B(1)=Mcal(1,i)
      !        B(2)=Mcal(2,i)
      !        B(3)=Mcal(3,i)
      do jtemp=2,4
         B(1)=Mcal(1,i)
         B(2)=Mcal(2,i)
         B(3)=Mcal(3,i)
         !          Axiao(1)=Mcal(1,jtemp)
         !          Axiao(2)=Mcal(2,jtemp)
         !          Axiao(3)=Mcal(3,jtemp)
         do j=1,3
            if(Mcal(j,jtemp).ne.0)then
               Yxiaotemp(jtemp,i,mtemp)=Ptemp(j)*Yxiaotemp(1,i,mtemp)+WPtemp(j)*Yxiaotemp(1,i,mtemp+1)
               if(B(j).ne.0)then
                  B(j)=Mcal(j,i)-1
                  itemp=trans(B(1),B(2),B(3))
                  Yxiaotemp(jtemp,i,mtemp)=Yxiaotemp(jtemp,i,mtemp)+Mcal(j,i)*ABCDtemp*Yxiaotemp(1,itemp,mtemp+1)
               endif
               goto 444
            endif
         enddo
         444     continue
      enddo
   enddo

End

subroutine GSPS(mtemp)
   use allmod
   Implicit double precision(a-h,o-z)
   double precision mem(35,35,0:8)
   !     double precision fact
   integer CPmem(35,35,0:8)
   ! double precision MEM(10,10,0:4)
   ! Integer CPMEM(10,10,0:4)

   integer MA(3),MB(3),NA(3),NB(3),LA(3),LB(3),B(3),Axiao(3),firstxiao,secondxiao
   double precision RA(3),RB(3),RC(3),RD(3),P(3),Q(3),W(3)
   double precision FM(0:13)
   double precision Qtemp(3),WQtemp(3),CDtemp,ABcom,Ptemp(3),WPtemp(3),ABtemp,CDcom,ABCDtemp
   COMMON /VRRcom/Qtemp,WQtemp,CDtemp,ABcom,Ptemp,WPtemp,ABtemp,CDcom,ABCDtemp

   COMMON /COM1/RA,RB,RC,RD
   COMMON /COM2/AA,BB,CC,DD,AB,CD,ROU,ABCD
   COMMON /COM4/P,Q,W
   COMMON /COM5/FM

   do i=21,35
      !        B(1)=Mcal(1,i)
      !        B(2)=Mcal(2,i)
      !        B(3)=Mcal(3,i)
      do jtemp=2,4
         B(1)=Mcal(1,i)
         B(2)=Mcal(2,i)
         B(3)=Mcal(3,i)
         !          Axiao(1)=Mcal(1,jtemp)
         !          Axiao(2)=Mcal(2,jtemp)
         !          Axiao(3)=Mcal(3,jtemp)
         do j=1,3
            if(Mcal(j,jtemp).ne.0)then
               Yxiaotemp(i,jtemp,mtemp)=Qtemp(j)*Yxiaotemp(i,1,mtemp)+WQtemp(j)*Yxiaotemp(i,1,mtemp+1)
               if(B(j).ne.0)then
                  B(j)=Mcal(j,i)-1
                  itemp=trans(B(1),B(2),B(3))
                  Yxiaotemp(i,jtemp,mtemp)=Yxiaotemp(i,jtemp,mtemp)+Mcal(j,i)*ABCDtemp*Yxiaotemp(itemp,1,mtemp+1)
               endif
               !             if(i.eq.5.and.jtemp.eq.3)then
               !               print*,B(1),B(2),B(3),trans(B(1),B(2),B(3)),itemp
               !               print*,'quick',Yxiaotemp(5,1,0),Yxiaotemp(5,1,1),Yxiaotemp(2,1,1),Yxiaotemp(5,3,0) &
                     !,Qtemp(2),WQtemp(2),ABCDtemp,Qtemp(2)*Yxiaotemp(5,1,0) &
                     !                               +WQtemp(2)*Yxiaotemp(5,1,1)+ABCDtemp*Yxiaotemp(2,1,1)
               !
               !             endif
               goto 333
            endif
         enddo
         333     continue
      enddo
   enddo

End

subroutine PSGS(mtemp)
   use allmod
   Implicit double precision(a-h,o-z)
   double precision mem(35,35,0:8)
   !     double precision fact
   integer CPmem(35,35,0:8)
   ! double precision MEM(10,10,0:4)
   ! Integer CPMEM(10,10,0:4)

   integer MA(3),MB(3),NA(3),NB(3),LA(3),LB(3),B(3),Axiao(3),firstxiao,secondxiao
   double precision RA(3),RB(3),RC(3),RD(3),P(3),Q(3),W(3)
   double precision FM(0:13)
   double precision Qtemp(3),WQtemp(3),CDtemp,ABcom,Ptemp(3),WPtemp(3),ABtemp,CDcom,ABCDtemp
   COMMON /VRRcom/Qtemp,WQtemp,CDtemp,ABcom,Ptemp,WPtemp,ABtemp,CDcom,ABCDtemp

   COMMON /COM1/RA,RB,RC,RD
   COMMON /COM2/AA,BB,CC,DD,AB,CD,ROU,ABCD
   COMMON /COM4/P,Q,W
   COMMON /COM5/FM

   do i=21,35
      !        B(1)=Mcal(1,i)
      !        B(2)=Mcal(2,i)
      !        B(3)=Mcal(3,i)
      do jtemp=2,4
         B(1)=Mcal(1,i)
         B(2)=Mcal(2,i)
         B(3)=Mcal(3,i)
         !          Axiao(1)=Mcal(1,jtemp)
         !          Axiao(2)=Mcal(2,jtemp)
         !          Axiao(3)=Mcal(3,jtemp)
         do j=1,3
            if(Mcal(j,jtemp).ne.0)then
               Yxiaotemp(jtemp,i,mtemp)=Ptemp(j)*Yxiaotemp(1,i,mtemp)+WPtemp(j)*Yxiaotemp(1,i,mtemp+1)
               if(B(j).ne.0)then
                  B(j)=Mcal(j,i)-1
                  itemp=trans(B(1),B(2),B(3))
                  Yxiaotemp(jtemp,i,mtemp)=Yxiaotemp(jtemp,i,mtemp)+Mcal(j,i)*ABCDtemp*Yxiaotemp(1,itemp,mtemp+1)
               endif
               goto 444
            endif
         enddo
         444     continue
      enddo
   enddo

End


subroutine FSDS(mtemp)
   use allmod
   Implicit double precision(a-h,o-z)
   double precision mem(35,35,0:8)
   !     double precision fact
   integer CPmem(35,35,0:8)
   ! double precision MEM(10,10,0:4)
   ! Integer CPMEM(10,10,0:4)

   integer MA(3),MB(3),NA(3),NB(3),LA(3),LB(3),B(3),Axiao(3),firstxiao,secondxiao
   double precision RA(3),RB(3),RC(3),RD(3),P(3),Q(3),W(3)
   double precision FM(0:13)
   double precision Qtemp(3),WQtemp(3),CDtemp,ABcom,Ptemp(3),WPtemp(3),ABtemp,CDcom,ABCDtemp
   COMMON /VRRcom/Qtemp,WQtemp,CDtemp,ABcom,Ptemp,WPtemp,ABtemp,CDcom,ABCDtemp

   COMMON /COM1/RA,RB,RC,RD
   COMMON /COM2/AA,BB,CC,DD,AB,CD,ROU,ABCD
   COMMON /COM4/P,Q,W
   COMMON /COM5/FM

   do i=5,10
      !        B(1)=Mcal(1,i)
      !        B(2)=Mcal(2,i)
      !        B(3)=Mcal(3,i)
      do jtemp=11,20
         B(1)=Mcal(1,i)
         B(2)=Mcal(2,i)
         B(3)=Mcal(3,i)
         Axiao(1)=Mcal(1,jtemp)
         Axiao(2)=Mcal(2,jtemp)
         Axiao(3)=Mcal(3,jtemp)
         do j=1,3
            if(Mcal(j,i).ne.0)then
               B(j)=Mcal(j,i)-1
               ixiao=trans(B(1),B(2),B(3))
               Yxiaotemp(jtemp,i,mtemp)=Qtemp(j)*Yxiaotemp(jtemp,ixiao,mtemp) &
                     +WQtemp(j)*Yxiaotemp(jtemp,ixiao,mtemp+1)
               if(Mcal(j,i).ge.2)then
                  !               B(j)=B(j)-1
                  !               itemp=trans(B(1),B(2),B(3))
                  !               Yxiaotemp(i,jtemp,mtemp)=Yxiaotemp(i,jtemp,mtemp)+ABtemp*(Yxiaotemp(itemp,jtemp,mtemp) &
                        !                                       -CDcom*Yxiaotemp(itemp,jtemp,mtemp+1))
                  Yxiaotemp(jtemp,i,mtemp)=Yxiaotemp(jtemp,i,mtemp)+CDtemp*(Yxiaotemp(jtemp,1,mtemp) &
                        -ABcom*Yxiaotemp(jtemp,1,mtemp+1))
               endif
               if(Axiao(j).ne.0)then
                  Axiao(j)=Axiao(j)-1
                  secondxiao=trans(Axiao(1),Axiao(2),Axiao(3))
                  Yxiaotemp(jtemp,i,mtemp)=Yxiaotemp(jtemp,i,mtemp) &
                        +Mcal(j,jtemp)*ABCDtemp*Yxiaotemp(secondxiao,ixiao,mtemp+1)
               endif
               goto 555
            endif
         enddo

         555     continue
      enddo
   enddo

end

subroutine DSFS(mtemp)
   use allmod
   Implicit double precision(a-h,o-z)
   double precision mem(35,35,0:8)
   !     double precision fact
   integer CPmem(35,35,0:8)
   ! double precision MEM(10,10,0:4)
   ! Integer CPMEM(10,10,0:4)

   integer MA(3),MB(3),NA(3),NB(3),LA(3),LB(3),B(3),Axiao(3),firstxiao,secondxiao
   double precision RA(3),RB(3),RC(3),RD(3),P(3),Q(3),W(3)
   double precision FM(0:13)
   double precision Qtemp(3),WQtemp(3),CDtemp,ABcom,Ptemp(3),WPtemp(3),ABtemp,CDcom,ABCDtemp
   COMMON /VRRcom/Qtemp,WQtemp,CDtemp,ABcom,Ptemp,WPtemp,ABtemp,CDcom,ABCDtemp

   COMMON /COM1/RA,RB,RC,RD
   COMMON /COM2/AA,BB,CC,DD,AB,CD,ROU,ABCD
   COMMON /COM4/P,Q,W
   COMMON /COM5/FM

   do i=5,10
      !        B(1)=Mcal(1,i)
      !        B(2)=Mcal(2,i)
      !        B(3)=Mcal(3,i)
      do jtemp=11,20
         B(1)=Mcal(1,i)
         B(2)=Mcal(2,i)
         B(3)=Mcal(3,i)
         Axiao(1)=Mcal(1,jtemp)
         Axiao(2)=Mcal(2,jtemp)
         Axiao(3)=Mcal(3,jtemp)
         do j=1,3
            if(Mcal(j,i).ne.0)then
               B(j)=Mcal(j,i)-1
               ixiao=trans(B(1),B(2),B(3))
               Yxiaotemp(i,jtemp,mtemp)=Ptemp(j)*Yxiaotemp(ixiao,jtemp,mtemp) &
                     +WPtemp(j)*Yxiaotemp(ixiao,jtemp,mtemp+1)
               if(Mcal(j,i).ge.2)then
                  !               B(j)=B(j)-1
                  !               itemp=trans(B(1),B(2),B(3))
                  !               Yxiaotemp(i,jtemp,mtemp)=Yxiaotemp(i,jtemp,mtemp)+ABtemp*(Yxiaotemp(itemp,jtemp,mtemp) &
                        !                                       -CDcom*Yxiaotemp(itemp,jtemp,mtemp+1))
                  Yxiaotemp(i,jtemp,mtemp)=Yxiaotemp(i,jtemp,mtemp)+ABtemp*(Yxiaotemp(1,jtemp,mtemp) &
                        -CDcom*Yxiaotemp(1,jtemp,mtemp+1))
               endif
               if(Axiao(j).ne.0)then
                  Axiao(j)=Axiao(j)-1
                  secondxiao=trans(Axiao(1),Axiao(2),Axiao(3))
                  Yxiaotemp(i,jtemp,mtemp)=Yxiaotemp(i,jtemp,mtemp) &
                        +Mcal(j,jtemp)*ABCDtemp*Yxiaotemp(ixiao,secondxiao,mtemp+1)
               endif
               goto 555
            endif
         enddo

         555     continue
      enddo
   enddo

End


subroutine GSDS(mtemp)
   use allmod
   Implicit double precision(a-h,o-z)
   double precision mem(35,35,0:8)
   !     double precision fact
   integer CPmem(35,35,0:8)
   ! double precision MEM(10,10,0:4)
   ! Integer CPMEM(10,10,0:4)

   integer MA(3),MB(3),NA(3),NB(3),LA(3),LB(3),B(3),Axiao(3),firstxiao,secondxiao
   double precision RA(3),RB(3),RC(3),RD(3),P(3),Q(3),W(3)
   double precision FM(0:13)
   double precision Qtemp(3),WQtemp(3),CDtemp,ABcom,Ptemp(3),WPtemp(3),ABtemp,CDcom,ABCDtemp
   COMMON /VRRcom/Qtemp,WQtemp,CDtemp,ABcom,Ptemp,WPtemp,ABtemp,CDcom,ABCDtemp

   COMMON /COM1/RA,RB,RC,RD
   COMMON /COM2/AA,BB,CC,DD,AB,CD,ROU,ABCD
   COMMON /COM4/P,Q,W
   COMMON /COM5/FM

   do i=5,10
      !        B(1)=Mcal(1,i)
      !        B(2)=Mcal(2,i)
      !        B(3)=Mcal(3,i)
      do jtemp=21,35
         B(1)=Mcal(1,i)
         B(2)=Mcal(2,i)
         B(3)=Mcal(3,i)
         Axiao(1)=Mcal(1,jtemp)
         Axiao(2)=Mcal(2,jtemp)
         Axiao(3)=Mcal(3,jtemp)
         do j=1,3
            if(Mcal(j,i).ne.0)then
               B(j)=Mcal(j,i)-1
               ixiao=trans(B(1),B(2),B(3))
               Yxiaotemp(jtemp,i,mtemp)=Qtemp(j)*Yxiaotemp(jtemp,ixiao,mtemp) &
                     +WQtemp(j)*Yxiaotemp(jtemp,ixiao,mtemp+1)
               if(Mcal(j,i).ge.2)then
                  !               B(j)=B(j)-1
                  !               itemp=trans(B(1),B(2),B(3))
                  !               Yxiaotemp(i,jtemp,mtemp)=Yxiaotemp(i,jtemp,mtemp)+ABtemp*(Yxiaotemp(itemp,jtemp,mtemp) &
                        !                                       -CDcom*Yxiaotemp(itemp,jtemp,mtemp+1))
                  Yxiaotemp(jtemp,i,mtemp)=Yxiaotemp(jtemp,i,mtemp)+CDtemp*(Yxiaotemp(jtemp,1,mtemp) &
                        -ABcom*Yxiaotemp(jtemp,1,mtemp+1))
               endif
               if(Axiao(j).ne.0)then
                  Axiao(j)=Axiao(j)-1
                  secondxiao=trans(Axiao(1),Axiao(2),Axiao(3))
                  Yxiaotemp(jtemp,i,mtemp)=Yxiaotemp(jtemp,i,mtemp) &
                        +Mcal(j,jtemp)*ABCDtemp*Yxiaotemp(secondxiao,ixiao,mtemp+1)
               endif
               goto 555
            endif
         enddo
         555     continue
      enddo
   enddo

end

subroutine DSGS(mtemp)
   use allmod
   Implicit double precision(a-h,o-z)
   double precision mem(35,35,0:8)
   !     double precision fact
   integer CPmem(35,35,0:8)
   ! double precision MEM(10,10,0:4)
   ! Integer CPMEM(10,10,0:4)

   integer MA(3),MB(3),NA(3),NB(3),LA(3),LB(3),B(3),Axiao(3),firstxiao,secondxiao
   double precision RA(3),RB(3),RC(3),RD(3),P(3),Q(3),W(3)
   double precision FM(0:13)
   double precision Qtemp(3),WQtemp(3),CDtemp,ABcom,Ptemp(3),WPtemp(3),ABtemp,CDcom,ABCDtemp
   COMMON /VRRcom/Qtemp,WQtemp,CDtemp,ABcom,Ptemp,WPtemp,ABtemp,CDcom,ABCDtemp

   COMMON /COM1/RA,RB,RC,RD
   COMMON /COM2/AA,BB,CC,DD,AB,CD,ROU,ABCD
   COMMON /COM4/P,Q,W
   COMMON /COM5/FM

   do i=5,10
      !        B(1)=Mcal(1,i)
      !        B(2)=Mcal(2,i)
      !        B(3)=Mcal(3,i)
      do jtemp=21,35
         B(1)=Mcal(1,i)
         B(2)=Mcal(2,i)
         B(3)=Mcal(3,i)
         Axiao(1)=Mcal(1,jtemp)
         Axiao(2)=Mcal(2,jtemp)
         Axiao(3)=Mcal(3,jtemp)
         do j=1,3
            if(Mcal(j,i).ne.0)then
               B(j)=Mcal(j,i)-1
               ixiao=trans(B(1),B(2),B(3))
               Yxiaotemp(i,jtemp,mtemp)=Ptemp(j)*Yxiaotemp(ixiao,jtemp,mtemp) &
                     +WPtemp(j)*Yxiaotemp(ixiao,jtemp,mtemp+1)
               if(Mcal(j,i).ge.2)then
                  !               B(j)=B(j)-1
                  !               itemp=trans(B(1),B(2),B(3))
                  !               Yxiaotemp(i,jtemp,mtemp)=Yxiaotemp(i,jtemp,mtemp)+ABtemp*(Yxiaotemp(itemp,jtemp,mtemp) &
                        !                                       -CDcom*Yxiaotemp(itemp,jtemp,mtemp+1))
                  Yxiaotemp(i,jtemp,mtemp)=Yxiaotemp(i,jtemp,mtemp)+ABtemp*(Yxiaotemp(1,jtemp,mtemp) &
                        -CDcom*Yxiaotemp(1,jtemp,mtemp+1))
               endif
               if(Axiao(j).ne.0)then
                  Axiao(j)=Axiao(j)-1
                  secondxiao=trans(Axiao(1),Axiao(2),Axiao(3))
                  Yxiaotemp(i,jtemp,mtemp)=Yxiaotemp(i,jtemp,mtemp) &
                        +Mcal(j,jtemp)*ABCDtemp*Yxiaotemp(ixiao,secondxiao,mtemp+1)
               endif
               goto 555
            endif
         enddo
         555     continue
      enddo
   enddo

End

subroutine FSFS(mtemp)
   use allmod
   Implicit double precision(a-h,o-z)
   double precision mem(35,35,0:8)
   !     double precision fact
   integer CPmem(35,35,0:8)
   ! double precision MEM(10,10,0:4)
   ! Integer CPMEM(10,10,0:4)

   integer MA(3),MB(3),NA(3),NB(3),LA(3),LB(3),B(3),Axiao(3),firstxiao,secondxiao
   double precision RA(3),RB(3),RC(3),RD(3),P(3),Q(3),W(3)
   double precision FM(0:13)
   double precision Qtemp(3),WQtemp(3),CDtemp,ABcom,Ptemp(3),WPtemp(3),ABtemp,CDcom,ABCDtemp
   COMMON /VRRcom/Qtemp,WQtemp,CDtemp,ABcom,Ptemp,WPtemp,ABtemp,CDcom,ABCDtemp

   COMMON /COM1/RA,RB,RC,RD
   COMMON /COM2/AA,BB,CC,DD,AB,CD,ROU,ABCD
   COMMON /COM4/P,Q,W
   COMMON /COM5/FM

   do i=11,20
      !        B(1)=Mcal(1,i)
      !        B(2)=Mcal(2,i)
      !        B(3)=Mcal(3,i)
      do jtemp=11,20
         B(1)=Mcal(1,i)
         B(2)=Mcal(2,i)
         B(3)=Mcal(3,i)
         Axiao(1)=Mcal(1,jtemp)
         Axiao(2)=Mcal(2,jtemp)
         Axiao(3)=Mcal(3,jtemp)
         do j=1,3
            if(Mcal(j,i).ne.0)then
               B(j)=Mcal(j,i)-1
               ixiao=trans(B(1),B(2),B(3))
               Yxiaotemp(jtemp,i,mtemp)=Qtemp(j)*Yxiaotemp(jtemp,ixiao,mtemp) &
                     +WQtemp(j)*Yxiaotemp(jtemp,ixiao,mtemp+1)
               if(Mcal(j,i).ge.2)then
                  !               B(j)=B(j)-1
                  !               itemp=trans(B(1),B(2),B(3))
                  !               Yxiaotemp(i,jtemp,mtemp)=Yxiaotemp(i,jtemp,mtemp)+ABtemp*(Yxiaotemp(itemp,jtemp,mtemp) &
                        !                                       -CDcom*Yxiaotemp(itemp,jtemp,mtemp+1))
                  B(j)=B(j)-1
                  ihigh=trans(B(1),B(2),B(3))
                  Yxiaotemp(jtemp,i,mtemp)=Yxiaotemp(jtemp,i,mtemp)+(B(j)+1)* &
                        CDtemp*(Yxiaotemp(jtemp,ihigh,mtemp) &
                        -ABcom*Yxiaotemp(jtemp,ihigh,mtemp+1))
               endif
               if(Axiao(j).ne.0)then
                  Axiao(j)=Axiao(j)-1
                  secondxiao=trans(Axiao(1),Axiao(2),Axiao(3))
                  Yxiaotemp(jtemp,i,mtemp)=Yxiaotemp(jtemp,i,mtemp) &
                        +Mcal(j,jtemp)*ABCDtemp*Yxiaotemp(secondxiao,ixiao,mtemp+1)
               endif
               goto 555
            endif
         enddo
         555     continue
      enddo
   enddo

end

subroutine GSFS(mtemp)
   use allmod
   Implicit double precision(a-h,o-z)
   double precision mem(35,35,0:8)
   !     double precision fact
   integer CPmem(35,35,0:8)
   ! double precision MEM(10,10,0:4)
   ! Integer CPMEM(10,10,0:4)

   integer MA(3),MB(3),NA(3),NB(3),LA(3),LB(3),B(3),Axiao(3),firstxiao,secondxiao
   double precision RA(3),RB(3),RC(3),RD(3),P(3),Q(3),W(3)
   double precision FM(0:13)
   double precision Qtemp(3),WQtemp(3),CDtemp,ABcom,Ptemp(3),WPtemp(3),ABtemp,CDcom,ABCDtemp
   COMMON /VRRcom/Qtemp,WQtemp,CDtemp,ABcom,Ptemp,WPtemp,ABtemp,CDcom,ABCDtemp

   COMMON /COM1/RA,RB,RC,RD
   COMMON /COM2/AA,BB,CC,DD,AB,CD,ROU,ABCD
   COMMON /COM4/P,Q,W
   COMMON /COM5/FM

   do i=11,20
      !        B(1)=Mcal(1,i)
      !        B(2)=Mcal(2,i)
      !        B(3)=Mcal(3,i)
      do jtemp=21,35
         B(1)=Mcal(1,i)
         B(2)=Mcal(2,i)
         B(3)=Mcal(3,i)
         Axiao(1)=Mcal(1,jtemp)
         Axiao(2)=Mcal(2,jtemp)
         Axiao(3)=Mcal(3,jtemp)
         do j=1,3
            if(Mcal(j,i).ne.0)then
               B(j)=Mcal(j,i)-1
               ixiao=trans(B(1),B(2),B(3))
               Yxiaotemp(jtemp,i,mtemp)=Qtemp(j)*Yxiaotemp(jtemp,ixiao,mtemp) &
                     +WQtemp(j)*Yxiaotemp(jtemp,ixiao,mtemp+1)
               if(Mcal(j,i).ge.2)then
                  !               B(j)=B(j)-1
                  !               itemp=trans(B(1),B(2),B(3))
                  !               Yxiaotemp(i,jtemp,mtemp)=Yxiaotemp(i,jtemp,mtemp)+ABtemp*(Yxiaotemp(itemp,jtemp,mtemp) &
                        !                                       -CDcom*Yxiaotemp(itemp,jtemp,mtemp+1))
                  B(j)=B(j)-1
                  ihigh=trans(B(1),B(2),B(3))
                  !                Yxiaotemp(jtemp,i,mtemp)=Yxiaotemp(jtemp,i,mtemp)+CDtemp*(Yxiaotemp(jtemp,1,mtemp) &
                        !                                       -ABcom*Yxiaotemp(jtemp,1,mtemp+1))
                  Yxiaotemp(jtemp,i,mtemp)=Yxiaotemp(jtemp,i,mtemp)+(B(j)+1)* &
                        CDtemp*(Yxiaotemp(jtemp,ihigh,mtemp) &
                        -ABcom*Yxiaotemp(jtemp,ihigh,mtemp+1))
               endif
               if(Axiao(j).ne.0)then
                  Axiao(j)=Axiao(j)-1
                  secondxiao=trans(Axiao(1),Axiao(2),Axiao(3))
                  Yxiaotemp(jtemp,i,mtemp)=Yxiaotemp(jtemp,i,mtemp) &
                        +Mcal(j,jtemp)*ABCDtemp*Yxiaotemp(secondxiao,ixiao,mtemp+1)
               endif
               goto 555
            endif
         enddo
         555     continue
      enddo
   enddo

end

subroutine FSGS(mtemp)
   use allmod
   Implicit double precision(a-h,o-z)
   double precision mem(35,35,0:8)
   !     double precision fact
   integer CPmem(35,35,0:8)
   ! double precision MEM(10,10,0:4)
   ! Integer CPMEM(10,10,0:4)

   integer MA(3),MB(3),NA(3),NB(3),LA(3),LB(3),B(3),Axiao(3),firstxiao,secondxiao
   double precision RA(3),RB(3),RC(3),RD(3),P(3),Q(3),W(3)
   double precision FM(0:13)
   double precision Qtemp(3),WQtemp(3),CDtemp,ABcom,Ptemp(3),WPtemp(3),ABtemp,CDcom,ABCDtemp
   COMMON /VRRcom/Qtemp,WQtemp,CDtemp,ABcom,Ptemp,WPtemp,ABtemp,CDcom,ABCDtemp

   COMMON /COM1/RA,RB,RC,RD
   COMMON /COM2/AA,BB,CC,DD,AB,CD,ROU,ABCD
   COMMON /COM4/P,Q,W
   COMMON /COM5/FM

   do i=11,20
      !        B(1)=Mcal(1,i)
      !        B(2)=Mcal(2,i)
      !        B(3)=Mcal(3,i)
      do jtemp=21,35
         B(1)=Mcal(1,i)
         B(2)=Mcal(2,i)
         B(3)=Mcal(3,i)
         Axiao(1)=Mcal(1,jtemp)
         Axiao(2)=Mcal(2,jtemp)
         Axiao(3)=Mcal(3,jtemp)
         do j=1,3
            if(Mcal(j,i).ne.0)then
               B(j)=Mcal(j,i)-1
               ixiao=trans(B(1),B(2),B(3))
               Yxiaotemp(i,jtemp,mtemp)=Ptemp(j)*Yxiaotemp(ixiao,jtemp,mtemp) &
                     +WPtemp(j)*Yxiaotemp(ixiao,jtemp,mtemp+1)
               if(Mcal(j,i).ge.2)then
                  !               B(j)=B(j)-1
                  !               itemp=trans(B(1),B(2),B(3))
                  !               Yxiaotemp(i,jtemp,mtemp)=Yxiaotemp(i,jtemp,mtemp)+ABtemp*(Yxiaotemp(itemp,jtemp,mtemp) &
                        !                                       -CDcom*Yxiaotemp(itemp,jtemp,mtemp+1))
                  B(j)=B(j)-1
                  ihigh=trans(B(1),B(2),B(3))
                  Yxiaotemp(i,jtemp,mtemp)=Yxiaotemp(i,jtemp,mtemp)+(B(j)+1)* &
                        ABtemp*(Yxiaotemp(ihigh,jtemp,mtemp) &
                        -CDcom*Yxiaotemp(ihigh,jtemp,mtemp+1))
                  !                Yxiaotemp(i,jtemp,mtemp)=Yxiaotemp(i,jtemp,mtemp)+ABtemp*(Yxiaotemp(1,jtemp,mtemp) &
                        !                                       -CDcom*Yxiaotemp(1,jtemp,mtemp+1))
               endif
               if(Axiao(j).ne.0)then
                  Axiao(j)=Axiao(j)-1
                  secondxiao=trans(Axiao(1),Axiao(2),Axiao(3))
                  Yxiaotemp(i,jtemp,mtemp)=Yxiaotemp(i,jtemp,mtemp) &
                        +Mcal(j,jtemp)*ABCDtemp*Yxiaotemp(ixiao,secondxiao,mtemp+1)
               endif
               goto 555
            endif
         enddo
         555     continue
      enddo
   enddo

End



subroutine GSGS(mtemp)
   use allmod
   Implicit double precision(a-h,o-z)
   double precision mem(35,35,0:8)
   !     double precision fact
   integer CPmem(35,35,0:8)
   ! double precision MEM(10,10,0:4)
   ! Integer CPMEM(10,10,0:4)

   integer MA(3),MB(3),NA(3),NB(3),LA(3),LB(3),B(3),Axiao(3),firstxiao,secondxiao
   double precision RA(3),RB(3),RC(3),RD(3),P(3),Q(3),W(3)
   double precision FM(0:13)
   double precision Qtemp(3),WQtemp(3),CDtemp,ABcom,Ptemp(3),WPtemp(3),ABtemp,CDcom,ABCDtemp
   COMMON /VRRcom/Qtemp,WQtemp,CDtemp,ABcom,Ptemp,WPtemp,ABtemp,CDcom,ABCDtemp

   COMMON /COM1/RA,RB,RC,RD
   COMMON /COM2/AA,BB,CC,DD,AB,CD,ROU,ABCD
   COMMON /COM4/P,Q,W
   COMMON /COM5/FM

   do i=21,35
      !        B(1)=Mcal(1,i)
      !        B(2)=Mcal(2,i)
      !        B(3)=Mcal(3,i)
      do jtemp=21,35
         B(1)=Mcal(1,i)
         B(2)=Mcal(2,i)
         B(3)=Mcal(3,i)
         Axiao(1)=Mcal(1,jtemp)
         Axiao(2)=Mcal(2,jtemp)
         Axiao(3)=Mcal(3,jtemp)
         do j=1,3
            if(Mcal(j,i).ne.0)then
               B(j)=Mcal(j,i)-1
               ixiao=trans(B(1),B(2),B(3))
               Yxiaotemp(jtemp,i,mtemp)=Qtemp(j)*Yxiaotemp(jtemp,ixiao,mtemp) &
                     +WQtemp(j)*Yxiaotemp(jtemp,ixiao,mtemp+1)
               if(Mcal(j,i).ge.2)then
                  !               B(j)=B(j)-1
                  !               itemp=trans(B(1),B(2),B(3))
                  !               Yxiaotemp(i,jtemp,mtemp)=Yxiaotemp(i,jtemp,mtemp)+ABtemp*(Yxiaotemp(itemp,jtemp,mtemp) &
                        !                                       -CDcom*Yxiaotemp(itemp,jtemp,mtemp+1))
                  B(j)=B(j)-1
                  ihigh=trans(B(1),B(2),B(3))
                  Yxiaotemp(jtemp,i,mtemp)=Yxiaotemp(jtemp,i,mtemp)+(B(j)+1)* &
                        CDtemp*(Yxiaotemp(jtemp,ihigh,mtemp) &
                        -ABcom*Yxiaotemp(jtemp,ihigh,mtemp+1))
               endif
               if(Axiao(j).ne.0)then
                  Axiao(j)=Axiao(j)-1
                  secondxiao=trans(Axiao(1),Axiao(2),Axiao(3))
                  Yxiaotemp(jtemp,i,mtemp)=Yxiaotemp(jtemp,i,mtemp) &
                        +Mcal(j,jtemp)*ABCDtemp*Yxiaotemp(secondxiao,ixiao,mtemp+1)
               endif
               goto 555
            endif
         enddo
         555     continue
      enddo
   enddo

end

subroutine BSLS(IBxiao,ILxiao,mtemp)
   use allmod
   Implicit double precision(a-h,o-z)
   double precision mem(35,35,0:8)
   !     double precision fact
   integer CPmem(35,35,0:8)
   ! double precision MEM(10,10,0:4)
   ! Integer CPMEM(10,10,0:4)
   integer IBxiao,ILxiao,mtemp

   integer MA(3),MB(3),NA(3),NB(3),LA(3),LB(3),B(3),Axiao(3),firstxiao,secondxiao
   double precision RA(3),RB(3),RC(3),RD(3),P(3),Q(3),W(3)
   double precision FM(0:13)
   double precision Qtemp(3),WQtemp(3),CDtemp,ABcom,Ptemp(3),WPtemp(3),ABtemp,CDcom,ABCDtemp
   COMMON /VRRcom/Qtemp,WQtemp,CDtemp,ABcom,Ptemp,WPtemp,ABtemp,CDcom,ABCDtemp

   COMMON /COM1/RA,RB,RC,RD
   COMMON /COM2/AA,BB,CC,DD,AB,CD,ROU,ABCD
   COMMON /COM4/P,Q,W
   COMMON /COM5/FM

   if(ILxiao.ne.0)then
      ! GSFS situation
      !     do i=11,20
      do i=Sumindex(ILxiao-1)+1,Sumindex(ILxiao)
         !        B(1)=Mcal(1,i)
         !        B(2)=Mcal(2,i)
         !        B(3)=Mcal(3,i)
         !        do jtemp=21,35
         do jtemp=Sumindex(IBxiao-1)+1,Sumindex(IBxiao)
            B(1)=Mcal(1,i)
            B(2)=Mcal(2,i)
            B(3)=Mcal(3,i)
            Axiao(1)=Mcal(1,jtemp)
            Axiao(2)=Mcal(2,jtemp)
            Axiao(3)=Mcal(3,jtemp)
            do j=1,3
               if(Mcal(j,i).ne.0)then
                  B(j)=Mcal(j,i)-1
                  ixiao=trans(B(1),B(2),B(3))
                  Yxiaotemp(jtemp,i,mtemp)=Qtemp(j)*Yxiaotemp(jtemp,ixiao,mtemp) &
                        +WQtemp(j)*Yxiaotemp(jtemp,ixiao,mtemp+1)
                  if(Mcal(j,i).ge.2)then
                     !               B(j)=B(j)-1
                     !               itemp=trans(B(1),B(2),B(3))
                     !               Yxiaotemp(i,jtemp,mtemp)=Yxiaotemp(i,jtemp,mtemp)+ABtemp*(Yxiaotemp(itemp,jtemp,mtemp) &
                           !                                       -CDcom*Yxiaotemp(itemp,jtemp,mtemp+1))
                     B(j)=B(j)-1
                     ihigh=trans(B(1),B(2),B(3))
                     !                Yxiaotemp(jtemp,i,mtemp)=Yxiaotemp(jtemp,i,mtemp)+CDtemp*(Yxiaotemp(jtemp,1,mtemp) &
                           !                                       -ABcom*Yxiaotemp(jtemp,1,mtemp+1))
                     Yxiaotemp(jtemp,i,mtemp)=Yxiaotemp(jtemp,i,mtemp)+(B(j)+1)* &
                           CDtemp*(Yxiaotemp(jtemp,ihigh,mtemp) &
                           -ABcom*Yxiaotemp(jtemp,ihigh,mtemp+1))
                  endif
                  if(Axiao(j).ne.0)then
                     Axiao(j)=Axiao(j)-1
                     secondxiao=trans(Axiao(1),Axiao(2),Axiao(3))
                     Yxiaotemp(jtemp,i,mtemp)=Yxiaotemp(jtemp,i,mtemp) &
                           +Mcal(j,jtemp)*ABCDtemp*Yxiaotemp(secondxiao,ixiao,mtemp+1)
                  endif
                  goto 555
               endif
            enddo
            555     continue
         enddo
      enddo

   else

      do i=Sumindex(IBxiao-1)+1,Sumindex(IBxiao)
         B(1)=Mcal(1,i)
         B(2)=Mcal(2,i)
         B(3)=Mcal(3,i)
         do j=1,3
            if(Mcal(j,i).ne.0)then
               B(j)=Mcal(j,i)-1
               itemp=trans(B(1),B(2),B(3))
               Yxiaotemp(i,1,mtemp)=Ptemp(j)*Yxiaotemp(itemp,1,mtemp)+WPtemp(j)*Yxiaotemp(itemp,1,mtemp+1)
               if(Mcal(j,i).gt.1)then
                  B(j)=Mcal(j,i)-2
                  inewtemp=trans(B(1),B(2),B(3))
                  Yxiaotemp(i,1,mtemp)=Yxiaotemp(i,1,mtemp)+(B(j)+1)*ABtemp*(Yxiaotemp(inewtemp,1,mtemp) &
                        -CDcom*Yxiaotemp(inewtemp,1,mtemp+1))
               endif
               goto 111
            endif
         enddo
         111 continue
      enddo

   endif

end

subroutine LSBS(ILxiao,IBxiao,mtemp)
   use allmod
   Implicit double precision(a-h,o-z)
   double precision mem(35,35,0:8)
   !     double precision fact
   integer CPmem(35,35,0:8)
   ! double precision MEM(10,10,0:4)
   ! Integer CPMEM(10,10,0:4)
   integer IBxiao,ILxiao,mtemp

   integer MA(3),MB(3),NA(3),NB(3),LA(3),LB(3),B(3),Axiao(3),firstxiao,secondxiao
   double precision RA(3),RB(3),RC(3),RD(3),P(3),Q(3),W(3)
   double precision FM(0:13)
   double precision Qtemp(3),WQtemp(3),CDtemp,ABcom,Ptemp(3),WPtemp(3),ABtemp,CDcom,ABCDtemp
   COMMON /VRRcom/Qtemp,WQtemp,CDtemp,ABcom,Ptemp,WPtemp,ABtemp,CDcom,ABCDtemp

   COMMON /COM1/RA,RB,RC,RD
   COMMON /COM2/AA,BB,CC,DD,AB,CD,ROU,ABCD
   COMMON /COM4/P,Q,W
   COMMON /COM5/FM

   if(ILxiao.ne.0)then

      do i=Sumindex(ILxiao-1)+1,Sumindex(ILxiao)
         !     do i=11,20
         !        B(1)=Mcal(1,i)
         !        B(2)=Mcal(2,i)
         !        B(3)=Mcal(3,i)
         do jtemp=Sumindex(IBxiao-1)+1,Sumindex(IBxiao)
            !        do jtemp=21,35
            B(1)=Mcal(1,i)
            B(2)=Mcal(2,i)
            B(3)=Mcal(3,i)
            Axiao(1)=Mcal(1,jtemp)
            Axiao(2)=Mcal(2,jtemp)
            Axiao(3)=Mcal(3,jtemp)
            do j=1,3
               if(Mcal(j,i).ne.0)then
                  B(j)=Mcal(j,i)-1
                  ixiao=trans(B(1),B(2),B(3))
                  Yxiaotemp(i,jtemp,mtemp)=Ptemp(j)*Yxiaotemp(ixiao,jtemp,mtemp) &
                        +WPtemp(j)*Yxiaotemp(ixiao,jtemp,mtemp+1)
                  if(Mcal(j,i).ge.2)then
                     !               B(j)=B(j)-1
                     !               itemp=trans(B(1),B(2),B(3))
                     !               Yxiaotemp(i,jtemp,mtemp)=Yxiaotemp(i,jtemp,mtemp)+ABtemp*(Yxiaotemp(itemp,jtemp,mtemp) &
                           !                                       -CDcom*Yxiaotemp(itemp,jtemp,mtemp+1))
                     B(j)=B(j)-1
                     ihigh=trans(B(1),B(2),B(3))
                     Yxiaotemp(i,jtemp,mtemp)=Yxiaotemp(i,jtemp,mtemp)+(B(j)+1)* &
                           ABtemp*(Yxiaotemp(ihigh,jtemp,mtemp) &
                           -CDcom*Yxiaotemp(ihigh,jtemp,mtemp+1))
                     !                Yxiaotemp(i,jtemp,mtemp)=Yxiaotemp(i,jtemp,mtemp)+ABtemp*(Yxiaotemp(1,jtemp,mtemp) &
                           !                                       -CDcom*Yxiaotemp(1,jtemp,mtemp+1))
                  endif
                  if(Axiao(j).ne.0)then
                     Axiao(j)=Axiao(j)-1
                     secondxiao=trans(Axiao(1),Axiao(2),Axiao(3))
                     Yxiaotemp(i,jtemp,mtemp)=Yxiaotemp(i,jtemp,mtemp) &
                           +Mcal(j,jtemp)*ABCDtemp*Yxiaotemp(ixiao,secondxiao,mtemp+1)
                  endif
                  goto 555
               endif
            enddo
            555     continue
         enddo
      enddo

   else

      do i=Sumindex(IBxiao-1)+1,Sumindex(IBxiao)
         B(1)=Mcal(1,i)
         B(2)=Mcal(2,i)
         B(3)=Mcal(3,i)
         do j=1,3
            if(Mcal(j,i).ne.0)then
               B(j)=Mcal(j,i)-1
               itemp=trans(B(1),B(2),B(3))
               Yxiaotemp(1,i,mtemp)=Qtemp(j)*Yxiaotemp(1,itemp,mtemp)+WQtemp(j)*Yxiaotemp(1,itemp,mtemp+1)
               if(Mcal(j,i).gt.1)then
                  B(j)=Mcal(j,i)-2
                  inewtemp=trans(B(1),B(2),B(3))
                  Yxiaotemp(1,i,mtemp)=Yxiaotemp(1,i,mtemp)+(B(j)+1)*CDtemp*(Yxiaotemp(1,inewtemp,mtemp) &
                        -ABcom*Yxiaotemp(1,inewtemp,mtemp+1))
               endif
               goto 111
            endif
         enddo
         111 continue
      enddo

   endif

End

