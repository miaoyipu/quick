! Be careful of (NB(itemp).eq.2.and.NB(jtemp).eq.2)
! Xiao HE 07/07/07 version
! BE careful of the array coefangxiaoL(20),allocate,main,module array,2eshell(opt),hrrsub,vertical
! subroutine hrrwhole(IJKLtype,III,JJJ,KKK,LLL,Y)
!Horrizontal Recursion subroutines by hand, these parts can be optimized by MAPLE
 subroutine hrrwhole
 use allmod

 Implicit real*8(A-H,O-Z)
 real*8 store(120,120)
 INTEGER NA(3),NB(3),NC(3),ND(3)
 real*8 RA(3),RB(3),RC(3),RD(3)
 Integer M1,M2,M3,M4

 real*8 coefangxiaoL(20),coefangxiaoR(20)
 integer angxiaoL(20),angxiaoR(20),numangularL,numangularR

 COMMON /COM1/RA,RB,RC,RD

 integer II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2
 common /xiaostore/store
 common /hrrstore/II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2

                   select case (IJKLtype)
                   
                   case (0,10,1000,1010)
!                   Do M=1,3
!                     NA(M)=KLMN(M,III)
!!                     NB(M)=KLMN(M,JJJ)
!                     NC(M)=KLMN(M,KKK)
!!                     ND(M)=KLMN(M,LLL) 
!                   Enddo

!                    Y=HRR(NA,NB,NC,ND)
!                    Y=Y*cons(III)*cons(JJJ)*cons(KKK)*cons(LLL)

!                     M1=trans(NA(1),NA(2),NA(3))
!                     M3=trans(NC(1),NC(2),NC(3))
                     M1=trans(KLMN(1,III),KLMN(2,III),KLMN(3,III))
                     M3=trans(KLMN(1,KKK),KLMN(2,KKK),KLMN(3,KKK))
                     Y=store(M1,M3)

                   case (2000,20,2010,1020,2020)
                     M1=trans(KLMN(1,III),KLMN(2,III),KLMN(3,III))
                     M3=trans(KLMN(1,KKK),KLMN(2,KKK),KLMN(3,KKK))
                     Y=store(M1,M3)

                     Y=Y*cons(III)*cons(JJJ)*cons(KKK)*cons(LLL)
 
                   case(100)
                   Do M=1,3
!                     NA(M)=KLMN(M,III)
                     NB(M)=KLMN(M,JJJ)
!                     NC(M)=KLMN(M,KKK)
!                     ND(M)=KLMN(M,LLL) 
                   Enddo

!                    Y=HRR(NA,NB,NC,ND)
!                    Y=Y*cons(III)*cons(JJJ)*cons(KKK)*cons(LLL)

                     M1=trans(NB(1),NB(2),NB(3))
!                     M3=trans(NC(1),NC(2),NC(3))
                     do itemp=1,3
                       if(NB(itemp).ne.0)then
                         Y=store(M1,1)+(RA(itemp)-RB(itemp))*store(1,1)
                         goto 111
                       endif
                     enddo

                  case(110)
                   Do M=1,3
!                     NA(M)=KLMN(M,III)
                     NB(M)=KLMN(M,JJJ)
                     NC(M)=KLMN(M,KKK)
!                     ND(M)=KLMN(M,LLL) 
                   Enddo

!                    Y=HRR(NA,NB,NC,ND)
!                    Y=Y*cons(III)*cons(JJJ)*cons(KKK)*cons(LLL)

                     M1=trans(NB(1),NB(2),NB(3))
                     M3=trans(NC(1),NC(2),NC(3))
                     do itemp=1,3
                       if(NB(itemp).ne.0)then
                         Y=store(M1,M3)+(RA(itemp)-RB(itemp))*store(1,M3)
                         goto 111
                       endif
                     enddo

                 case(101)
                   Do M=1,3
!                     NA(M)=KLMN(M,III)
                     NB(M)=KLMN(M,JJJ)
!                     NC(M)=KLMN(M,KKK)
                     ND(M)=KLMN(M,LLL)
                   Enddo

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
117 continue

                     do jtemp=1,3
                       if(NB(jtemp).ne.0)then
                         Y=Y1+(RA(jtemp)-RB(jtemp))*Y2
                         goto 111
                       endif
                     enddo

                case(111)
                   Do M=1,3
!                     NA(M)=KLMN(M,III)
                     NB(M)=KLMN(M,JJJ)
                     NC(M)=KLMN(M,KKK)
                     ND(M)=KLMN(M,LLL)
                   Enddo

!                    Y=HRR(NA,NB,NC,ND)
!                    Y=Y*cons(III)*cons(JJJ)*cons(KKK)*cons(LLL)

                     MB=trans(NB(1),NB(2),NB(3))
                     MCX=trans(NC(1),NC(2),NC(3))
                     MCD=trans(NC(1)+ND(1),NC(2)+ND(2),NC(3)+ND(3))

                     do itemp=1,3
                       if(ND(itemp).ne.0)then
                         ctemp=(RC(itemp)-RD(itemp))
                         Y1=store(MB,MCD)+ctemp*store(MB,MCX)
                         Y2=store(1,MCD)+ctemp*store(1,MCX)
                         goto 1230
                       endif
                     enddo
1230 continue

                     do jtemp=1,3
                       if(NB(jtemp).ne.0)then
                         Y=Y1+(RA(jtemp)-RB(jtemp))*Y2
                         goto 111
                       endif
                     enddo

              case(1100)
                   Do M=1,3
                     NA(M)=KLMN(M,III)
                     NB(M)=KLMN(M,JJJ)
!                     NC(M)=KLMN(M,KKK)
!                     ND(M)=KLMN(M,LLL) 
                   Enddo

!                    Y=HRR(NA,NB,NC,ND)
!                    Y=Y*cons(III)*cons(JJJ)*cons(KKK)*cons(LLL)

                     MA=trans(NA(1),NA(2),NA(3))
                     MAB=trans(NA(1)+NB(1),NA(2)+NB(2),NA(3)+NB(3))
                     do itemp=1,3
                       if(NB(itemp).ne.0)then
                         Y=store(MAB,1)+(RA(itemp)-RB(itemp))*store(MA,1)
                         goto 111
                       endif
                     enddo

             case(1110)
                   Do M=1,3
                     NA(M)=KLMN(M,III)
                     NB(M)=KLMN(M,JJJ)
                     NC(M)=KLMN(M,KKK)
!                     ND(M)=KLMN(M,LLL) 
                   Enddo

!                    Y=HRR(NA,NB,NC,ND)
!                    Y=Y*cons(III)*cons(JJJ)*cons(KKK)*cons(LLL)

                     MA=trans(NA(1),NA(2),NA(3))
                     MAB=trans(NA(1)+NB(1),NA(2)+NB(2),NA(3)+NB(3))
                     MCX=trans(NC(1),NC(2),NC(3))
                     do itemp=1,3
                       if(NB(itemp).ne.0)then
                         Y=store(MAB,MCX)+(RA(itemp)-RB(itemp))*store(MA,MCX)
                         goto 111
                       endif
                     enddo

            case(1101)
                   Do M=1,3
                     NA(M)=KLMN(M,III)
                     NB(M)=KLMN(M,JJJ)
!                     NC(M)=KLMN(M,KKK)
                     ND(M)=KLMN(M,LLL)
                   Enddo

!                    Y=HRR(NA,NB,NC,ND)
!                    Y=Y*cons(III)*cons(JJJ)*cons(KKK)*cons(LLL)

                     MA=trans(NA(1),NA(2),NA(3))
                     MAB=trans(NA(1)+NB(1),NA(2)+NB(2),NA(3)+NB(3))
                     MCX=trans(ND(1),ND(2),ND(3))

                     do itemp=1,3
                       if(ND(itemp).ne.0)then
                         ctemp=(RC(itemp)-RD(itemp))
                         Y1=store(MAB,MCX)+ctemp*store(MAB,1)
                         Y2=store(MA,MCX)+ctemp*store(MA,1)
                         goto 135
                       endif
                     enddo
135 continue

                     do jtemp=1,3
                       if(NB(jtemp).ne.0)then
                         Y=Y1+(RA(jtemp)-RB(jtemp))*Y2
                         goto 111
                       endif
                     enddo

           case(1111)
                   Do M=1,3
                     NA(M)=KLMN(M,III)
                     NB(M)=KLMN(M,JJJ)
                     NC(M)=KLMN(M,KKK)
                     ND(M)=KLMN(M,LLL)
                   Enddo

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
141 continue

                     do jtemp=1,3
                       if(NB(jtemp).ne.0)then
                         Y=Y1+(RA(jtemp)-RB(jtemp))*Y2
                         goto 111
                       endif
                     enddo

                  case(1)
                   Do M=1,3
!                     NA(M)=KLMN(M,III)
!                     NB(M)=KLMN(M,JJJ)
!                     NC(M)=KLMN(M,KKK)
                     ND(M)=KLMN(M,LLL)
                   Enddo

!                    Y=HRR(NA,NB,NC,ND)
!                    Y=Y*cons(III)*cons(JJJ)*cons(KKK)*cons(LLL)

!                     M1=trans(NB(1),NB(2),NB(3))
                     MD=trans(ND(1),ND(2),ND(3))
                     do itemp=1,3
                       if(ND(itemp).ne.0)then
                         Y=store(1,MD)+(RC(itemp)-RD(itemp))*store(1,1)
                         goto 111
                       endif
                     enddo

                   case(11)
                   Do M=1,3
!                     NA(M)=KLMN(M,III)
!                     NB(M)=KLMN(M,JJJ)
                     NC(M)=KLMN(M,KKK)
                     ND(M)=KLMN(M,LLL)
                   Enddo

!                    Y=HRR(NA,NB,NC,ND)
!                    Y=Y*cons(III)*cons(JJJ)*cons(KKK)*cons(LLL)

                     MCD=trans(NC(1)+ND(1),NC(2)+ND(2),NC(3)+ND(3))
                     MCX=trans(NC(1),NC(2),NC(3))
                     do itemp=1,3
                       if(ND(itemp).ne.0)then
                         Y=store(1,MCD)+(RC(itemp)-RD(itemp))*store(1,MCX)
                         goto 111
                       endif
                     enddo

                  case(1001)
                   Do M=1,3
                     NA(M)=KLMN(M,III)
!                     NB(M)=KLMN(M,JJJ)
!                     NC(M)=KLMN(M,KKK)
                     ND(M)=KLMN(M,LLL)
                   Enddo

!                    Y=HRR(NA,NB,NC,ND)
!                    Y=Y*cons(III)*cons(JJJ)*cons(KKK)*cons(LLL)

                     MA=trans(NA(1),NA(2),NA(3))
                     MD=trans(ND(1),ND(2),ND(3))
                     do itemp=1,3
                       if(ND(itemp).ne.0)then
                         Y=store(MA,MD)+(RC(itemp)-RD(itemp))*store(MA,1)
                         goto 111
                       endif
                     enddo

                 case(1011)
                   Do M=1,3
                     NA(M)=KLMN(M,III)
!                     NB(M)=KLMN(M,JJJ)
                     NC(M)=KLMN(M,KKK)
                     ND(M)=KLMN(M,LLL)
                   Enddo

!                    Y=HRR(NA,NB,NC,ND)
!                    Y=Y*cons(III)*cons(JJJ)*cons(KKK)*cons(LLL)

                     MA=trans(NA(1),NA(2),NA(3))
                     MCX=trans(NC(1),NC(2),NC(3))
                     MCD=trans(NC(1)+ND(1),NC(2)+ND(2),NC(3)+ND(3))
                     do itemp=1,3
                       if(ND(itemp).ne.0)then
                         Y=store(MA,MCD)+(RC(itemp)-RD(itemp))*store(MA,MCX)
                         goto 111
                       endif
                     enddo

           case(999) 

             call lefthrr(RA,RB,KLMN(1:3,III),KLMN(1:3,JJJ),IJtype,coefangxiaoL,angxiaoL,numangularL)
             call lefthrr(RC,RD,KLMN(1:3,KKK),KLMN(1:3,LLL),KLtype,coefangxiaoR,angxiaoR,numangularR)
                  
                  Y=0.0d0
                  do ixiao=1,numangularL
                    do jxiao=1,numangularR
                      Y=Y+coefangxiaoL(ixiao)*coefangxiaoR(jxiao)* &
                          store(angxiaoL(ixiao),angxiaoR(jxiao))
!                 if(III.eq.1.and.JJJ.eq.2.and.KKK.eq.2.and.LLL.eq.30)then
!                   print*,'xiao11',Y,coefangxiaoL(ixiao),coefangxiaoR(jxiao),store(angxiaoL(ixiao),angxiaoR(jxiao))
!                 endif
!                 if(III.eq.1.and.JJJ.eq.4.and.KKK.eq.7.and.LLL.eq.19)then
!                   print*,'xiao111',Y,coefangxiaoL(ixiao),coefangxiaoR(jxiao),store(angxiaoL(ixiao),angxiaoR(jxiao))
!                 endif
                    enddo
                  enddo
!                  print*,numangularL,numangularR,coefangxiaoL(numangularL),coefangxiaoR(numangularR),Y
                  
                  Y=Y*cons(III)*cons(JJJ)*cons(KKK)*cons(LLL)

                  end select

111 continue

      end


     subroutine lefthrr(RA,RB,KLMNA,KLMNB,IKnumber,coefangxiao,angxiao,numangular)

 use allmod

 Implicit real*8(A-H,O-Z)
! real*8 store(35,35)
! INTEGER NA(3),NB(3),NC(3),ND(3)
! real*8 RA(3),RB(3),RC(3),RD(3)
 INTEGER KLMNA(3),KLMNB(3),NA(3),NB(3)
 real*8 RA(3),RB(3)
 Integer M1,M2,M3,M4

 real*8 coefangxiao(20)
 integer angxiao(20),numangular,IKnumber

! integer II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2
! common /xiaostore/store
! common /hrrstore/II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2

                   select case (IKnumber)

                   case (0)
                   numangular=1
                   coefangxiao(1)=1.0d0
                   angxiao(1)=1
                   
                   case (1)
                   Do M=1,3
!                     NA(M)=KLMN(M,III)
                     NB(M)=KLMNB(M)
!                     NC(M)=KLMN(M,KKK)
!                     ND(M)=KLMN(M,LLL) 
                   Enddo

                     M1=trans(NB(1),NB(2),NB(3))

                     do itemp=1,3
                       if(NB(itemp).ne.0)then
                         Y=(RA(itemp)-RB(itemp))
                         goto 111
                       endif
                     enddo

111                numangular=2
                   coefangxiao(1)=1.0d0
                   angxiao(1)=M1
                   coefangxiao(2)=Y
                   angxiao(2)=1
                   

                   case (10)
                   Do M=1,3
                     NA(M)=KLMNA(M)
!                     NB(M)=KLMN(M,JJJ)
!                     NC(M)=KLMN(M,KKK)
!                     ND(M)=KLMN(M,LLL) 
                   Enddo
        
                     M1=trans(NA(1),NA(2),NA(3))
             
                   numangular=1
                   coefangxiao(1)=1.0d0 
                   angxiao(1)=M1

                   case (11)
                   Do M=1,3
                     NA(M)=KLMNA(M)
                     NB(M)=KLMNB(M)
!                     NC(M)=KLMN(M,KKK)
!                     ND(M)=KLMN(M,LLL) 
                   Enddo

                     M1=trans(NA(1)+NB(1),NA(2)+NB(2),NA(3)+NB(3))
                     M2=trans(NA(1),NA(2),NA(3))

                     do itemp=1,3
                       if(NB(itemp).ne.0)then
                         Y=(RA(itemp)-RB(itemp))
                         goto 222
                       endif
                     enddo

222                numangular=2
                   coefangxiao(1)=1.0d0 
                   angxiao(1)=M1
                   coefangxiao(2)=Y
                   angxiao(2)=M2

                   case (20,30,40)
                   Do M=1,3
                     NA(M)=KLMNA(M)
!                     NB(M)=KLMN(M,JJJ)
!                     NC(M)=KLMN(M,KKK)
!                     ND(M)=KLMN(M,LLL) 
                   Enddo
        
                     M1=trans(NA(1),NA(2),NA(3))
      
                   numangular=1
                   coefangxiao(1)=1.0d0 
                   angxiao(1)=M1

                   case (2)
                   Do M=1,3
!                    NA(M)=KLMN(M,III)
                     NB(M)=KLMNB(M)
!                     NC(M)=KLMN(M,KKK)
!                     ND(M)=KLMN(M,LLL) 
                   Enddo

                     M1=trans(NB(1),NB(2),NB(3))
                     coefangxiao(1)=1.0d0
                     angxiao(1)=M1

                     do itemp=1,3
                         if(NB(itemp).eq.2)then
                           numangular=3
                           Y=(RA(itemp)-RB(itemp))
                           coefangxiao(2)=2.0d0*Y
                           angxiao(2)=itemp+1
                           coefangxiao(3)=Y*Y
                           angxiao(3)=1
                           goto 333
                         endif
                         do jtemp=itemp+1,3
                           if(NB(itemp).eq.1.and.NB(jtemp).eq.1)then
                             numangular=4
                             coefangxiao(2)=(RA(itemp)-RB(itemp))
                             angxiao(2)=jtemp+1
                             coefangxiao(3)=(RA(jtemp)-RB(jtemp))
                             angxiao(3)=itemp+1
                             coefangxiao(4)=coefangxiao(2)*coefangxiao(3)
                             angxiao(4)=1
                             goto 333
                           endif
                        enddo
                     enddo

333                continue

                   case (21,31,41)
                   Do M=1,3
                     NA(M)=KLMNA(M)
                     NB(M)=KLMNB(M)
!                     NC(M)=KLMN(M,KKK)
!                     ND(M)=KLMN(M,LLL) 
                   Enddo

                     M1=trans(NA(1)+NB(1),NA(2)+NB(2),NA(3)+NB(3))
                     M2=trans(NA(1),NA(2),NA(3))
                     coefangxiao(1)=1.0d0
                     angxiao(1)=M1

                     do itemp=1,3
                         if(NB(itemp).eq.1)then
                           numangular=2
                           coefangxiao(2)=(RA(itemp)-RB(itemp))
                           angxiao(2)=M2
                           goto 444
                         endif
                     enddo

444                continue


                   case (12)
                   Do M=1,3
                     NA(M)=KLMNA(M)
                     NB(M)=KLMNB(M)
!                     NC(M)=KLMN(M,KKK)
!                     ND(M)=KLMN(M,LLL) 
                   Enddo

                     M1=trans(NA(1)+NB(1),NA(2)+NB(2),NA(3)+NB(3))
                     M2=trans(NA(1),NA(2),NA(3))
                     coefangxiao(1)=1.0d0
                     angxiao(1)=M1

                     do itemp=1,3
                         if(NB(itemp).eq.2)then
                           numangular=3
                           Y=(RA(itemp)-RB(itemp))
                           coefangxiao(2)=2.0d0*Y
                           NA(itemp)=NA(itemp)+1
                           angxiao(2)=trans(NA(1),NA(2),NA(3))
                           coefangxiao(3)=Y*Y
                           angxiao(3)=M2
                           goto 555
                         endif
                         do jtemp=itemp+1,3
                           if(NB(itemp).eq.1.and.NB(jtemp).eq.1)then
                             numangular=4
                             coefangxiao(2)=(RA(itemp)-RB(itemp))
                             NA(jtemp)=NA(jtemp)+1
                             angxiao(2)=trans(NA(1),NA(2),NA(3))
                             NA(jtemp)=NA(jtemp)-1
                             
                             coefangxiao(3)=(RA(jtemp)-RB(jtemp))
                             NA(itemp)=NA(itemp)+1
                             angxiao(3)=trans(NA(1),NA(2),NA(3))
                             coefangxiao(4)=coefangxiao(2)*coefangxiao(3)
                             angxiao(4)=M2
                             goto 555
                           endif
                        enddo
                     enddo

555                continue


                   case (22,32,42)
                   Do M=1,3
                     NA(M)=KLMNA(M)
                     NB(M)=KLMNB(M)
!                     NC(M)=KLMN(M,KKK)
!                     ND(M)=KLMN(M,LLL) 
                   Enddo

                     M1=trans(NA(1)+NB(1),NA(2)+NB(2),NA(3)+NB(3))
                     M2=trans(NA(1),NA(2),NA(3))
                     coefangxiao(1)=1.0d0
                     angxiao(1)=M1

                     do itemp=1,3
                         if(NB(itemp).eq.2)then
                           numangular=3
                           Y=(RA(itemp)-RB(itemp))
                           coefangxiao(2)=2.0d0*Y
                           NA(itemp)=NA(itemp)+1
                           angxiao(2)=trans(NA(1),NA(2),NA(3))
                           coefangxiao(3)=Y*Y
                           angxiao(3)=M2
                           goto 666
                         endif
                         do jtemp=itemp+1,3
                           if(NB(itemp).eq.1.and.NB(jtemp).eq.1)then
                             numangular=4
                             coefangxiao(2)=(RA(itemp)-RB(itemp))
                             NA(jtemp)=NA(jtemp)+1
                             angxiao(2)=trans(NA(1),NA(2),NA(3))
                             NA(jtemp)=NA(jtemp)-1

                             coefangxiao(3)=(RA(jtemp)-RB(jtemp))
                             NA(itemp)=NA(itemp)+1
                             angxiao(3)=trans(NA(1),NA(2),NA(3))
                             coefangxiao(4)=coefangxiao(2)*coefangxiao(3)
                             angxiao(4)=M2
                             goto 666
                           endif
                        enddo
                     enddo

666                continue

                   case(3,13,23,33,43)
                   Do M=1,3
                     NA(M)=KLMNA(M)
                     NB(M)=KLMNB(M)
!                     NC(M)=KLMN(M,KKK)
!                     ND(M)=KLMN(M,LLL) 
                   Enddo

                     M1=trans(NA(1)+NB(1),NA(2)+NB(2),NA(3)+NB(3))
                     M2=trans(NA(1),NA(2),NA(3))
                     coefangxiao(1)=1.0d0
                     angxiao(1)=M1

                     do itemp=1,3
                         if(NB(itemp).eq.3)then
                           numangular=4
                           Y=(RA(itemp)-RB(itemp))

                           coefangxiao(2)=3.0d0*Y
                           NA(itemp)=NA(itemp)+2
                           angxiao(2)=trans(NA(1),NA(2),NA(3))
                           NA(itemp)=NA(itemp)-2

                           coefangxiao(3)=3.0d0*Y*Y
                           NA(itemp)=NA(itemp)+1
                           angxiao(3)=trans(NA(1),NA(2),NA(3))
                           NA(itemp)=NA(itemp)-1

                           coefangxiao(4)=Y*Y*Y
                           angxiao(4)=M2
                           goto 777
                         endif

                         do jtemp=1,3
                           if(NB(itemp).eq.1.and.NB(jtemp).eq.2)then
                             numangular=6

                             Yxiaotemp1=(RA(itemp)-RB(itemp))
                             coefangxiao(2)=Yxiaotemp1
                             NA(jtemp)=NA(jtemp)+2
                             angxiao(2)=trans(NA(1),NA(2),NA(3))
                             NA(jtemp)=NA(jtemp)-2

                             Yxiaotemp2=(RA(jtemp)-RB(jtemp))
                             coefangxiao(3)=2.0d0*Yxiaotemp2
                             NA(itemp)=NA(itemp)+1
                             NA(jtemp)=NA(jtemp)+1
                             angxiao(3)=trans(NA(1),NA(2),NA(3))
                             NA(itemp)=NA(itemp)-1
                             NA(jtemp)=NA(jtemp)-1

                             coefangxiao(4)=2.0d0*Yxiaotemp1*Yxiaotemp2
                             NA(jtemp)=NA(jtemp)+1
                             angxiao(4)=trans(NA(1),NA(2),NA(3))
                             NA(jtemp)=NA(jtemp)-1

                             coefangxiao(5)=Yxiaotemp2*Yxiaotemp2
                             NA(itemp)=NA(itemp)+1
                             angxiao(5)=trans(NA(1),NA(2),NA(3))
                             NA(itemp)=NA(itemp)-1

                             coefangxiao(6)=Yxiaotemp1*Yxiaotemp2*Yxiaotemp2
                             angxiao(6)=M2

                             goto 777
                           endif

                        enddo
                     enddo

                           if(NB(1).eq.1.and.NB(2).eq.1)then
                             numangular=8

                             Yxiaotemp1=(RA(1)-RB(1))
                             Yxiaotemp2=(RA(2)-RB(2))
                             Yxiaotemp3=(RA(3)-RB(3))

                             coefangxiao(2)=Yxiaotemp1
                             NA(2)=NA(2)+1
                             NA(3)=NA(3)+1
                             angxiao(2)=trans(NA(1),NA(2),NA(3))
                             NA(2)=NA(2)-1
                             NA(3)=NA(3)-1

                             coefangxiao(3)=Yxiaotemp2
                             NA(1)=NA(1)+1
                             NA(3)=NA(3)+1
                             angxiao(3)=trans(NA(1),NA(2),NA(3))
                             NA(1)=NA(1)-1
                             NA(3)=NA(3)-1

                             coefangxiao(4)=Yxiaotemp3
                             NA(1)=NA(1)+1
                             NA(2)=NA(2)+1
                             angxiao(4)=trans(NA(1),NA(2),NA(3))
                             NA(1)=NA(1)-1
                             NA(2)=NA(2)-1

                             coefangxiao(5)=Yxiaotemp1*Yxiaotemp2
                             NA(3)=NA(3)+1
                             angxiao(5)=trans(NA(1),NA(2),NA(3))
                             NA(3)=NA(3)-1

                             coefangxiao(6)=Yxiaotemp1*Yxiaotemp3
                             NA(2)=NA(2)+1
                             angxiao(6)=trans(NA(1),NA(2),NA(3))
                             NA(2)=NA(2)-1

                             coefangxiao(7)=Yxiaotemp2*Yxiaotemp3
                             NA(1)=NA(1)+1
                             angxiao(7)=trans(NA(1),NA(2),NA(3))
                             NA(1)=NA(1)-1

                             coefangxiao(8)=Yxiaotemp1*Yxiaotemp2*Yxiaotemp3
                             angxiao(8)=M2

                             goto 777
                           endif

777                continue

                   case(4,14,24,34,44)
                   Do M=1,3
                     NA(M)=KLMNA(M)
                     NB(M)=KLMNB(M)
!                     NC(M)=KLMN(M,KKK)
!                     ND(M)=KLMN(M,LLL) 
                   Enddo

                     M1=trans(NA(1)+NB(1),NA(2)+NB(2),NA(3)+NB(3))
                     M2=trans(NA(1),NA(2),NA(3))
                     coefangxiao(1)=1.0d0
                     angxiao(1)=M1

                     do itemp=1,3
                         if(NB(itemp).eq.4)then
                           numangular=5
                           Y=(RA(itemp)-RB(itemp))

                           coefangxiao(2)=4.0d0*Y
                           NA(itemp)=NA(itemp)+3
                           angxiao(2)=trans(NA(1),NA(2),NA(3))
                           NA(itemp)=NA(itemp)-3

                           coefangxiao(3)=6.0d0*Y*Y
                           NA(itemp)=NA(itemp)+2
                           angxiao(3)=trans(NA(1),NA(2),NA(3))
                           NA(itemp)=NA(itemp)-2

                           coefangxiao(4)=4.0d0*Y*Y*Y
                           NA(itemp)=NA(itemp)+1
                           angxiao(4)=trans(NA(1),NA(2),NA(3))
                           NA(itemp)=NA(itemp)-1

                           coefangxiao(5)=Y*Y*Y*Y
                           angxiao(5)=M2

                           goto 888
                         endif

                         do jtemp=1,3
                           if(NB(itemp).eq.1.and.NB(jtemp).eq.3)then
                             numangular=8

                             Yxiaotemp1=(RA(itemp)-RB(itemp))
                             coefangxiao(2)=Yxiaotemp1
                             NA(jtemp)=NA(jtemp)+3
                             angxiao(2)=trans(NA(1),NA(2),NA(3))
                             NA(jtemp)=NA(jtemp)-3

                             Yxiaotemp2=(RA(jtemp)-RB(jtemp))
                             coefangxiao(3)=3.0d0*Yxiaotemp2
                             NA(itemp)=NA(itemp)+1
                             NA(jtemp)=NA(jtemp)+2
                             angxiao(3)=trans(NA(1),NA(2),NA(3))
                             NA(itemp)=NA(itemp)-1
                             NA(jtemp)=NA(jtemp)-2

                             coefangxiao(4)=3.0d0*Yxiaotemp1*Yxiaotemp2
                             NA(jtemp)=NA(jtemp)+2
                             angxiao(4)=trans(NA(1),NA(2),NA(3))
                             NA(jtemp)=NA(jtemp)-2

                             coefangxiao(5)=3.0d0*Yxiaotemp2*Yxiaotemp2
                             NA(itemp)=NA(itemp)+1
                             NA(jtemp)=NA(jtemp)+1
                             angxiao(5)=trans(NA(1),NA(2),NA(3))
                             NA(itemp)=NA(itemp)-1
                             NA(jtemp)=NA(jtemp)-1

                             coefangxiao(6)=3.0d0*Yxiaotemp1*Yxiaotemp2*Yxiaotemp2
                             NA(jtemp)=NA(jtemp)+1
                             angxiao(6)=trans(NA(1),NA(2),NA(3))
                             NA(jtemp)=NA(jtemp)-1

                             coefangxiao(7)=Yxiaotemp2*Yxiaotemp2*Yxiaotemp2
                             NA(itemp)=NA(itemp)+1
                             angxiao(7)=trans(NA(1),NA(2),NA(3))
                             NA(itemp)=NA(itemp)-1

                             coefangxiao(8)=Yxiaotemp1*Yxiaotemp2*Yxiaotemp2*Yxiaotemp2
                             angxiao(8)=M2

                             goto 888
                           endif

                           if(NB(itemp).eq.2.and.NB(jtemp).eq.2.and.itemp.ne.jtemp)then
                             numangular=9

                             Yxiaotemp1=(RA(itemp)-RB(itemp))
                             Yxiaotemp2=(RA(jtemp)-RB(jtemp))

                             coefangxiao(2)=2.0d0*Yxiaotemp1
                             NA(itemp)=NA(itemp)+1
                             NA(jtemp)=NA(jtemp)+2
                             angxiao(2)=trans(NA(1),NA(2),NA(3))
                             NA(itemp)=NA(itemp)-1
                             NA(jtemp)=NA(jtemp)-2

                             coefangxiao(3)=2.0d0*Yxiaotemp2
                             NA(itemp)=NA(itemp)+2
                             NA(jtemp)=NA(jtemp)+1
                             angxiao(3)=trans(NA(1),NA(2),NA(3))
                             NA(itemp)=NA(itemp)-2
                             NA(jtemp)=NA(jtemp)-1

                             coefangxiao(4)=4.0d0*Yxiaotemp1*Yxiaotemp2
                             NA(itemp)=NA(itemp)+1
                             NA(jtemp)=NA(jtemp)+1
                             angxiao(4)=trans(NA(1),NA(2),NA(3))
                             NA(itemp)=NA(itemp)-1
                             NA(jtemp)=NA(jtemp)-1

                             coefangxiao(5)=Yxiaotemp1*Yxiaotemp1
                             NA(jtemp)=NA(jtemp)+2
                             angxiao(5)=trans(NA(1),NA(2),NA(3))
                             NA(jtemp)=NA(jtemp)-2

                             coefangxiao(6)=Yxiaotemp2*Yxiaotemp2
                             NA(itemp)=NA(itemp)+2
                             angxiao(6)=trans(NA(1),NA(2),NA(3))
                             NA(itemp)=NA(itemp)-2

                             coefangxiao(7)=2.0d0*Yxiaotemp1*Yxiaotemp2*Yxiaotemp2
                             NA(itemp)=NA(itemp)+1
                             angxiao(7)=trans(NA(1),NA(2),NA(3))
                             NA(itemp)=NA(itemp)-1

                             coefangxiao(8)=2.0d0*Yxiaotemp1*Yxiaotemp1*Yxiaotemp2
                             NA(jtemp)=NA(jtemp)+1
                             angxiao(8)=trans(NA(1),NA(2),NA(3))
                             NA(jtemp)=NA(jtemp)-1

                             coefangxiao(9)=Yxiaotemp1*Yxiaotemp1*Yxiaotemp2*Yxiaotemp2
                             angxiao(9)=M2

                             goto 888
                           endif

                           do ktemp=1,3
                           if(NB(itemp).eq.1.and.NB(jtemp).eq.1.and.NB(ktemp).eq.2.and.itemp.ne.jtemp)then
                             numangular=12

                             Yxiaotemp1=(RA(itemp)-RB(itemp))
                             Yxiaotemp2=(RA(jtemp)-RB(jtemp))
                             Yxiaotemp3=(RA(ktemp)-RB(ktemp))

                             coefangxiao(2)=Yxiaotemp1
                             NA(jtemp)=NA(jtemp)+1
                             NA(ktemp)=NA(ktemp)+2
                             angxiao(2)=trans(NA(1),NA(2),NA(3))
                             NA(jtemp)=NA(jtemp)-1
                             NA(ktemp)=NA(ktemp)-2

                             coefangxiao(3)=Yxiaotemp2
                             NA(itemp)=NA(itemp)+1
                             NA(ktemp)=NA(ktemp)+2
                             angxiao(3)=trans(NA(1),NA(2),NA(3))
                             NA(itemp)=NA(itemp)-1
                             NA(ktemp)=NA(ktemp)-2

                             coefangxiao(4)=2.0d0*Yxiaotemp3
                             NA(itemp)=NA(itemp)+1
                             NA(jtemp)=NA(jtemp)+1
                             NA(ktemp)=NA(ktemp)+1
                             angxiao(4)=trans(NA(1),NA(2),NA(3))
                             NA(itemp)=NA(itemp)-1
                             NA(jtemp)=NA(jtemp)-1
                             NA(ktemp)=NA(ktemp)-1

                             coefangxiao(5)=Yxiaotemp1*Yxiaotemp2
                             NA(ktemp)=NA(ktemp)+2
                             angxiao(5)=trans(NA(1),NA(2),NA(3))
                             NA(ktemp)=NA(ktemp)-2

                             coefangxiao(6)=2.0d0*Yxiaotemp1*Yxiaotemp3
                             NA(jtemp)=NA(jtemp)+1
                             NA(ktemp)=NA(ktemp)+1
                             angxiao(6)=trans(NA(1),NA(2),NA(3))
                             NA(jtemp)=NA(jtemp)-1
                             NA(ktemp)=NA(ktemp)-1

                             coefangxiao(7)=2.0d0*Yxiaotemp2*Yxiaotemp3
                             NA(itemp)=NA(itemp)+1
                             NA(ktemp)=NA(ktemp)+1
                             angxiao(7)=trans(NA(1),NA(2),NA(3))
                             NA(itemp)=NA(itemp)-1
                             NA(ktemp)=NA(ktemp)-1

                             coefangxiao(8)=Yxiaotemp3*Yxiaotemp3
                             NA(itemp)=NA(itemp)+1
                             NA(jtemp)=NA(jtemp)+1
                             angxiao(8)=trans(NA(1),NA(2),NA(3))
                             NA(itemp)=NA(itemp)-1
                             NA(jtemp)=NA(jtemp)-1

                             coefangxiao(9)=2.0d0*Yxiaotemp1*Yxiaotemp2*Yxiaotemp3
                             NA(ktemp)=NA(ktemp)+1
                             angxiao(9)=trans(NA(1),NA(2),NA(3))
                             NA(ktemp)=NA(ktemp)-1

                             coefangxiao(10)=Yxiaotemp1*Yxiaotemp3*Yxiaotemp3
                             NA(jtemp)=NA(jtemp)+1
                             angxiao(10)=trans(NA(1),NA(2),NA(3))
                             NA(jtemp)=NA(jtemp)-1

                             coefangxiao(11)=Yxiaotemp2*Yxiaotemp3*Yxiaotemp3
                             NA(itemp)=NA(itemp)+1
                             angxiao(11)=trans(NA(1),NA(2),NA(3))
                             NA(itemp)=NA(itemp)-1

                             coefangxiao(12)=Yxiaotemp1*Yxiaotemp2*Yxiaotemp3*Yxiaotemp3
                             angxiao(12)=M2

                             goto 888
                           endif                                                          

                           enddo 
                        enddo
                     enddo

888                continue

                   end select

                   End 

! PAY ATTENTION TO THE INDICE OF IJTYPE AND KLTYPE
! Xiao HE 07/07/07 version
! subroutine hrrwhole(IJKLtype,III,JJJ,KKK,LLL,Y)
!Horrizontal Recursion subroutines by hand, these parts can be optimized by MAPLE
 subroutine hrrwholeopt
 use allmod

 Implicit real*8(A-H,O-Z)
 real*8 store(120,120)
 INTEGER NA(3),NB(3),NC(3),ND(3)
 real*8 RA(3),RB(3),RC(3),RD(3)
 Integer M1,M2,M3,M4

 real*8 storeaa(120,120)
 real*8 storebb(120,120)
 real*8 storecc(120,120)
 real*8 storedd(120,120)

 real*8 coefangxiaoL(20),coefangxiaoR(20)
 integer angxiaoL(20),angxiaoR(20),numangularL,numangularR

 real*8 coefangxiaoLnew(20),coefangxiaoRnew(20)
 integer angxiaoLnew(20),angxiaoRnew(20),numangularLnew,numangularRnew

 COMMON /COM1/RA,RB,RC,RD

 integer II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2
 common /xiaostore/store
 common /xiaostoreopt/storeaa,storebb,storecc,storedd
 common /hrrstore/II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2

 tempconstant=cons(III)*cons(JJJ)*cons(KKK)*cons(LLL)

             call lefthrr(RA,RB,KLMN(1:3,III),KLMN(1:3,JJJ),IJtype,coefangxiaoL,angxiaoL,numangularL)
             call lefthrr(RC,RD,KLMN(1:3,KKK),KLMN(1:3,LLL),KLtype,coefangxiaoR,angxiaoR,numangularR)

!                  Y=0.0d0
!                  do ixiao=1,numangularL
!                    do jxiao=1,numangularR
!                      Y=Y+coefangxiaoL(ixiao)*coefangxiaoR(jxiao)* &
!                          store(angxiaoL(ixiao),angxiaoR(jxiao))
!                    enddo
!                  enddo
!
!                  Y=Y*cons(III)*cons(JJJ)*cons(KKK)*cons(LLL)

! alpha (a+1b|cd)
                
!                IJtype=IJtype+10           

           do itemp=1,3

              do itempxiao=1,3
                NA(itempxiao)=KLMN(itempxiao,III) 
              enddo
         
              NA(itemp)=KLMN(itemp,III)+1

             call lefthrr(RA,RB,NA(1:3),KLMN(1:3,JJJ),IJtype+10,coefangxiaoLnew,angxiaoLnew,numangularLnew)
!             call lefthrr(RC,RD,KLMN(1:3,KKK),KLMN(1:3,LLL),KLtype,coefangxiaoR,angxiaoR,numangularR)

                  Yaa(itemp)=0.0d0
                  do ixiao=1,numangularLnew
                    do jxiao=1,numangularR
                      Yaa(itemp)=Yaa(itemp)+coefangxiaoLnew(ixiao)*coefangxiaoR(jxiao)* &
                          storeAA(angxiaoLnew(ixiao),angxiaoR(jxiao))
                    enddo
                  enddo
 
         
              if(KLMN(itemp,III).ge.1)then

              do itempxiao=1,3
                NA(itempxiao)=KLMN(itempxiao,III)
              enddo

                NA(itemp)=KLMN(itemp,III)-1

             call lefthrr(RA,RB,NA(1:3),KLMN(1:3,JJJ),IJtype-10,coefangxiaoLnew,angxiaoLnew,numangularLnew)
!             call lefthrr(RC,RD,KLMN(1:3,KKK),KLMN(1:3,LLL),KLtype,coefangxiaoR,angxiaoR,numangularR)

                  do ixiao=1,numangularLnew
                    do jxiao=1,numangularR
                      Yaa(itemp)=Yaa(itemp)-KLMN(itemp,III)*coefangxiaoLnew(ixiao)* &
                          coefangxiaoR(jxiao)*store(angxiaoLnew(ixiao),angxiaoR(jxiao))
                    enddo
                  enddo
              endif

                  Yaa(itemp)=Yaa(itemp)*tempconstant
!                  if(III.le.12.and.JJJ.le.12.and.KKK.ge.13.and.LLL.ge.13)then
!                    print*,III,JJJ,KKK,LLL,itemp,Yaa(itemp)
!                    print*,III,JJJ,KKK,LLL,itemp,Yaa(itemp),numangularLnew,numangularR
!                  do ixiao=1,numangularLnew
!                    do jxiao=1,numangularR
!                      print*,coefangxiaoLnew(ixiao),coefangxiaoR(jxiao), &
!                          storeAA(angxiaoLnew(ixiao),angxiaoR(jxiao)), &
!                          coefangxiaoLnew(ixiao)*coefangxiaoR(jxiao)* &
!                          storeAA(angxiaoLnew(ixiao),angxiaoR(jxiao))
!                    enddo
!                  enddo
!                  endif
!                  if(III.le.39.and.JJJ.le.39.and.KKK.ge.40.and.LLL.ge.40)then
!                    print*,KKK-39,LLL-39,III+12,JJJ+12,itemp,Yaa(itemp),numangularLnew,numangularR
!                  endif
          enddo
     
!        IJtype=IJtype-10

!(ab+1|cd) beta

!        IJtype=IJtype+1

           do itemp=1,3

              do itempxiao=1,3
                NA(itempxiao)=KLMN(itempxiao,III)
                NB(itempxiao)=KLMN(itempxiao,JJJ)
              enddo

              NB(itemp)=KLMN(itemp,JJJ)+1

             call lefthrr(RA,RB,NA(1:3),NB(1:3),IJtype+1,coefangxiaoLnew,angxiaoLnew,numangularLnew)
!             call lefthrr(RC,RD,KLMN(1:3,KKK),KLMN(1:3,LLL),KLtype,coefangxiaoR,angxiaoR,numangularR)

                  Ybb(itemp)=0.0d0
                  do ixiao=1,numangularLnew
                    do jxiao=1,numangularR
                      Ybb(itemp)=Ybb(itemp)+coefangxiaoLnew(ixiao)*coefangxiaoR(jxiao)* &
                          storeBB(angxiaoLnew(ixiao),angxiaoR(jxiao))
                    enddo
                  enddo

              if(KLMN(itemp,JJJ).ge.1)then

              do itempxiao=1,3
                NA(itempxiao)=KLMN(itempxiao,III)
                NB(itempxiao)=KLMN(itempxiao,JJJ)
              enddo

                NB(itemp)=KLMN(itemp,JJJ)-1

             call lefthrr(RA,RB,NA(1:3),NB(1:3),IJtype-1,coefangxiaoLnew,angxiaoLnew,numangularLnew)
!             call lefthrr(RC,RD,KLMN(1:3,KKK),KLMN(1:3,LLL),KLtype,coefangxiaoR,angxiaoR,numangularR)

                  do ixiao=1,numangularLnew
                    do jxiao=1,numangularR
                      Ybb(itemp)=Ybb(itemp)-KLMN(itemp,JJJ)*coefangxiaoLnew(ixiao)* &
                          coefangxiaoR(jxiao)*store(angxiaoLnew(ixiao),angxiaoR(jxiao))
                    enddo
                  enddo
              endif

                  Ybb(itemp)=Ybb(itemp)*tempconstant
          enddo

!        IJtype=IJtype-1

!(ab|c+1d) gamma

!        KLtype=KLtype+10

           do itemp=1,3

              do itempxiao=1,3
                NC(itempxiao)=KLMN(itempxiao,KKK)
              enddo

              NC(itemp)=KLMN(itemp,KKK)+1

!             call lefthrr(RA,RB,NA(1:3),NB(1:3),IJtype,coefangxiaoLnew,angxiaoLnew,numangularLnew)
             call lefthrr(RC,RD,NC(1:3),KLMN(1:3,LLL),KLtype+10,coefangxiaoRnew,angxiaoRnew,numangularRnew)

                  Ycc(itemp)=0.0d0
                  do ixiao=1,numangularL
                    do jxiao=1,numangularRnew
                      Ycc(itemp)=Ycc(itemp)+coefangxiaoL(ixiao)*coefangxiaoRnew(jxiao)* &
                          storeCC(angxiaoL(ixiao),angxiaoRnew(jxiao))
                    enddo
                  enddo

!                  if(III.ge.30.and.JJJ.ge.30.and.KKK.ge.30.and.LLL.ge.40 &
!                     .and.III.le.39.and.JJJ.le.39.and.KKK.le.39)then
!                    print*,III,JJJ,KKK,LLL,itemp,Ycc(itemp),KLtype,NC,RC,RD
!                  do ixiao=1,numangularL
!                    do jxiao=1,numangularRnew
!                      print*,ixiao,jxiao,coefangxiaoL(ixiao),coefangxiaoRnew(jxiao), &
!                          storeCC(angxiaoL(ixiao),angxiaoRnew(jxiao))
!                    enddo
!                  enddo
!                  endif

!                          if(III.eq.1.and.JJJ.eq.2.and.KKK.eq.3.and.LLL.eq.10)then
!                            if(itemp.eq.1.or.itemp.eq.2)then
!                              print*,'xiao',YCC(1),storecc(1,8),YCC(2),storecc(1,5)
!                            endif
!                          endif

              if(KLMN(itemp,KKK).ge.1)then

              do itempxiao=1,3
                NC(itempxiao)=KLMN(itempxiao,KKK)
              enddo

                NC(itemp)=KLMN(itemp,KKK)-1

!             call lefthrr(RA,RB,NA(1:3),NB(1:3),IJtype,coefangxiaoLnew,angxiaoLnew,numangularLnew)
             call lefthrr(RC,RD,NC(1:3),KLMN(1:3,LLL),KLtype-10,coefangxiaoRnew,angxiaoRnew,numangularRnew)

                  do ixiao=1,numangularL
                    do jxiao=1,numangularRnew
                      Ycc(itemp)=Ycc(itemp)-KLMN(itemp,KKK)*coefangxiaoL(ixiao)* &
                          coefangxiaoRnew(jxiao)*store(angxiaoL(ixiao),angxiaoRnew(jxiao))
                    enddo
                  enddo
              endif

!                          if(III.eq.1.and.JJJ.eq.2.and.KKK.eq.3.and.LLL.eq.10)then
!                            if(itemp.eq.1.or.itemp.eq.2)then 
!                              print*,'xiao',YCC(1),store(1,1),YCC(2),&
!NC(1:3),KLMN(1:3,LLL),KLtype-10,coefangxiaoRnew,angxiaoRnew,numangularRnew
!                            endif
!                          endif

                  Ycc(itemp)=Ycc(itemp)*tempconstant

!                  if(III.ge.42.and.JJJ.ge.42.and.KKK.ge.42.and.LLL.le.12)then
!                  if(III.ge.30.and.JJJ.ge.30.and.KKK.ge.30.and.LLL.ge.40 &
!                     .and.III.le.39.and.JJJ.le.39.and.KKK.le.39)then
!                    print*,III,JJJ,KKK,LLL,itemp,Ycc(itemp),KLtype,NC,RC,RD
!                  do ixiao=1,numangularL
!                    do jxiao=1,numangularRnew
!                      print*,ixiao,jxiao,coefangxiaoL(ixiao),coefangxiaoRnew(jxiao), &
!                          store(angxiaoL(ixiao),angxiaoRnew(jxiao))
!                    enddo
!                  enddo
!                  endif
          enddo

!      KLtype=KLtype-10

      end

