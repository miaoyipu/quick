! Ed Brothers. November 27, 2001
! Xiao HE. September 14,2008
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

subroutine calmp2
  use allmod
  use quick_gaussian_class_module
  implicit double precision(a-h,o-z)

  double precision cutoffTest,testtmp,testCutoff
  integer II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2
  common /hrrstore/II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2

  call PrtAct(ioutfile,"Begin MP2 Calculation")
  cutoffmp2=1.0d-8  ! cutoff criteria
  primlimit=1.0d-8
  emp2=0.0d0

  ! occupied and virtual orbitals number
  iocc=Nelec/2
  ivir=Nbasis-Nelec/2

  ! calculate memory usage and determine steps
  ememorysum=real(iocc*ivir*nbasis*8.0d0/1024.0d0/1024.0d0/1024.0d0)
  write(ioutfile,'("CURRENT MEMORY USAGE=",E12.6,"M")') ememorysum

  ! actually nstep is step length
  nstep=min(int(1.5d0/ememorysum),Nelec/2)

  ! if with f orbital
  if(quick_method%ffunxiao)then
     nbasistemp=6
  else
     nbasistemp=10
  endif

  ! Allocate some variables
  allocate(mp2shell(nbasis))
  allocate(orbmp2(ivir,ivir))
  allocate(orbmp2i331(nstep,nbasis,nbasistemp,nbasistemp,2))
  allocate(orbmp2j331(nstep,ivir,nbasistemp,nbasistemp,2))
  allocate(orbmp2k331(nstep,iocc,ivir,nbasis))

  ! with nstep(acutally, it represetns step lenght), we can 
  ! have no. of steps for mp2 calculation
  nstepmp2=nelec/2/nstep
  nstepmp2=nstepmp2+1
  if(nstep*(nstepmp2-1).eq.nelec/2)then
     nstepmp2=nstepmp2-1
  endif
  write(ioutfile,'("TOTAL STEP          =",I6)') nstepmp2

  ! Pre-step for density cutoff
  call densityCutoff    

  ! first save coeffecient.
  Do i=1,nbasis
     Do j=1,nbasis
        HOLD(i,j)=CO(j,i)
     ENDDO
  ENDDO

  ttt=MAXVAL(Ycutoff) ! Max Value of Ycutoff

  Do i3new=1,nstepmp2               ! Step counter

     call cpu_time(timer_begin%TMP2)
     ntemp=0    ! integer counter
     nstepmp2s=(i3new-1)*nstep+1    ! Step start n
     nstepmp2f=i3new*nstep          ! Step end n

     if(i3new.eq.nstepmp2)nstepmp2f=nelec/2
     nsteplength=nstepmp2f-nstepmp2s+1  ! Step Lengh, from nstepmp2s to nstepmp2f
     
     ! Initial orbmp2k331
     call initialOrbmp2k331(orbmp2k331,nstep,nbasis,ivir,iocc,nsteplength)
     Do II=1,jshell
        do JJ=II,jshell
           if(Ycutoff(II,JJ).gt.cutoffmp2/ttt)then
              call initialOrbmp2ij(orbmp2i331,nstep,nsteplength,nbasis,nbasistemp,nbasistemp)
              call initialOrbmp2ij(orbmp2j331,nstep,nsteplength,ivir,nbasistemp,nbasistemp)

              do KK=1,jshell
                 do LL=KK,jshell

                    ! Schwarts cutoff is implemented here
                    comax=0.d0
                    testCutoff = Ycutoff(II,JJ)*Ycutoff(KK,LL)
                    If(testCutoff.gt.cutoffmp2)then

                       NKK1=Qstart(KK)
                       NKK2=Qfinal(KK)
                       NLL1=Qstart(LL)
                       NLL2=Qfinal(LL)

                       NBK1=Qsbasis(KK,NKK1)
                       NBK2=Qfbasis(KK,NKK2)
                       NBL1=Qsbasis(LL,NLL1)
                       NBL2=Qfbasis(LL,NLL2)

                       KK111=Ksumtype(KK)+NBK1
                       KK112=Ksumtype(KK)+NBK2
                       LL111=Ksumtype(LL)+NBL1
                       LL112=Ksumtype(LL)+NBL2

                       do icycle=1,nsteplength
                            i3=nstepmp2s+icycle-1
                            Do KKK=KK111,KK112
                                Do LLL=max(KKK,LL111),LL112
                                    comax=max(comax,dabs(co(kkk,i3)))
                                    comax=max(comax,dabs(co(lll,i3)))    
                                Enddo
                            Enddo        
                       enddo
                       
                       testCutoff=testCutoff*comax
                       If(testCutoff.gt.cutoffmp2)then
                          dnmax=comax
                          ntemp=ntemp+1
                          call shellmp2(nstepmp2s,nsteplength)
                       Endif

                    endif

                 enddo
              enddo



              NII1=Qstart(II)
              NII2=Qfinal(II)
              NJJ1=Qstart(JJ)
              NJJ2=Qfinal(JJ)

              NBI1=Qsbasis(II,NII1)
              NBI2=Qfbasis(II,NII2)
              NBJ1=Qsbasis(JJ,NJJ1)
              NBJ2=Qfbasis(JJ,NJJ2)

              II111=Ksumtype(II)+NBI1
              II112=Ksumtype(II)+NBI2
              JJ111=Ksumtype(JJ)+NBJ1
              JJ112=Ksumtype(JJ)+NBJ2

              Do III=II111,II112
                 Do JJJ=max(III,JJ111),JJ112

                    IIInew=III-II111+1
                    JJJnew=JJJ-JJ111+1

                    do LLL=1,nbasis
                       do j33=1,ivir
                          j33new=j33+nelec/2
                          atemp=CO(LLL,j33new)
                          do icycle=1,nsteplength
                             orbmp2j331(icycle,j33,IIInew,JJJnew,1)=orbmp2j331(icycle,j33,IIInew,JJJnew,1) + &
                                  orbmp2i331(icycle,LLL,IIInew,JJJnew,1)*atemp
                             if(III.ne.JJJ)then
                                orbmp2j331(icycle,j33,JJJnew,IIInew,2)=orbmp2j331(icycle,j33,JJJnew,IIInew,2) + &
                                     orbmp2i331(icycle,LLL,JJJnew,IIInew,2)*atemp
                             endif
                          enddo
                       enddo
                    enddo

                    do j33=1,ivir
                       do k33=1,nelec/2
                          atemp=HOLD(k33,III)
                          atemp2=HOLD(k33,JJJ)
                          do icycle=1,nsteplength
                             orbmp2k331(icycle,k33,j33,JJJ)=orbmp2k331(icycle,k33,j33,JJJ)+ &
                                  orbmp2j331(icycle,j33,IIInew,JJJnew,1)*atemp
                             if(III.ne.JJJ)then
                                orbmp2k331(icycle,k33,j33,III)=orbmp2k331(icycle,k33,j33,III)+ &
                                     orbmp2j331(icycle,j33,JJJnew,IIInew,2)*atemp2
                             endif
                          enddo
                       enddo
                    enddo
                 Enddo
              Enddo
           endif

        enddo
     enddo

     write (ioutfile,'("EFFECT INTEGRALS    =",i8)') ntemp
     
     do icycle=1,nsteplength
        do k3=1,nelec/2
            do j3=1,nbasis-nelec/2
            do l3=1,nbasis-nelec/2
            L3new=l3+nelec/2
                do lll=1,nbasis
                 write(10,*) k3,j3,LLL,L3new,orbmp2k331(icycle,k3,j3,LLL),co(LLL,L3new)
                enddo
            enddo
            enddo
        enddo
     enddo

     do icycle=1,nsteplength
        i3=nstepmp2s+icycle-1
        Do k3=i3,nelec/2

           Do J3=1,nbasis-nelec/2
              Do L3=1,nbasis-nelec/2
                 orbmp2(L3,J3)=0.0d0
                 L3new=L3+nelec/2
                 Do LLL=1,nbasis
                    orbmp2(L3,J3)=orbmp2(L3,J3)+orbmp2k331(icycle,k3,j3,LLL)*CO(LLL,L3new)
                 Enddo
              Enddo
           Enddo
       
           Do J3=1,nbasis-nelec/2
              Do L3=1,nbasis-nelec/2
                 if(k3.gt.i3)then
                    emp2=emp2+2.0d0/(E(i3)+E(k3)-E(j3+nelec/2)-E(l3+nelec/2)) &
                         *orbmp2(j3,l3)*(2.0d0*orbmp2(j3,l3)-orbmp2(l3,j3))
                 endif
                 if(k3.eq.i3)then
                    emp2=emp2+1.0d0/(E(i3)+E(k3)-E(j3+nelec/2)-E(l3+nelec/2)) &
                         *orbmp2(j3,l3)*(2.0d0*orbmp2(j3,l3)-orbmp2(l3,j3))
                 endif

              Enddo
           Enddo
        enddo
     Enddo
        
     call cpu_time(timer_end%TMP2)
     timer_cumer%TMP2=timer_end%TMP2-timer_begin%TMP2+timer_cumer%TMP2

  ENDDO

  write (iOutFile,'("SECOND ORDER ENERGY =",F16.9)') EMP2
  write (iOutFile,'("EMP2                =",F16.9)') Etot+EMP2
  call PrtAct(ioutfile,"End MP2 Calculation")
  return
end subroutine calmp2




! Ed Brothers. November 27, 2001
! Xiao HE. September 14,2008
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP
subroutine MPI_calmp2
  use allmod
  use quick_gaussian_class_module
  implicit double precision(a-h,o-z)

  include "mpif.h"

  double precision cutoffTest,testtmp,testCutoff
  double precision, allocatable:: temp4d(:,:,:,:)
  integer II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2
  common /hrrstore/II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2

  if (bMPI) then
     call MPI_BCAST(DENSE,nbasis*nbasis,mpi_double_precision,0,MPI_COMM_WORLD,mpierror)
     call MPI_BCAST(CO,nbasis*nbasis,mpi_double_precision,0,MPI_COMM_WORLD,mpierror)
     call MPI_BCAST(E,nbasis,mpi_double_precision,0,MPI_COMM_WORLD,mpierror)
  endif


  cutoffmp2=1.0d-8  ! cutoff criteria
  primlimit=1.0d-8
  emp2=0.0d0

  ! occupied and virtual orbitals number
  iocc=Nelec/2
  ivir=Nbasis-Nelec/2

  ! calculate memory usage and determine steps
  ememorysum=real(iocc*ivir*nbasis*8.0d0/1024.0d0/1024.0d0/1024.0d0)
  if (master) then
     call PrtAct(ioutfile,"Begin MP2 Calculation")
     write(ioutfile,'("Current memory usage=",E12.6,"M")') ememorysum
  endif

  ! actually nstep is step length
  nstep=min(int(1.5d0/ememorysum),Nelec/2)

  ! if with f orbital
  if(quick_method%ffunxiao)then
     nbasistemp=6
  else
     nbasistemp=10
  endif

  ! Allocate some variables
  allocate(mp2shell(nbasis))
  allocate(orbmp2(ivir,ivir))
  allocate(orbmp2i331(nstep,nbasis,nbasistemp,nbasistemp,2))
  allocate(orbmp2j331(nstep,ivir,nbasistemp,nbasistemp,2))
  allocate(orbmp2k331(nstep,iocc,ivir,nbasis))
  allocate(temp4d(nstep,iocc,ivir,nbasis))

  ! with nstep(acutally, it represetns step lenght), we can 
  ! have no. of steps for mp2 calculation
  nstepmp2=nelec/2/nstep
  nstepmp2=nstepmp2+1
  if(nstep*(nstepmp2-1).eq.nelec/2)then
     nstepmp2=nstepmp2-1
  endif
  write(ioutfile,'("TOTAL STEP          =",I6)') nstepmp2
  ! Pre-step for density cutoff
  call densityCutoff    

  ! first save coeffecient.
  Do i=1,nbasis
     Do j=1,nbasis
        HOLD(i,j)=CO(j,i)
     ENDDO
  ENDDO

  ttt=MAXVAL(Ycutoff) ! Max Value of Ycutoff

  Do i3new=1,nstepmp2               ! Step counter
     call cpu_time(timer_begin%TMP2)
     ntemp=0    ! integer counter
     nstepmp2s=(i3new-1)*nstep+1    ! Step start n
     nstepmp2f=i3new*nstep          ! Step end n

     if(i3new.eq.nstepmp2)nstepmp2f=nelec/2
     nsteplength=nstepmp2f-nstepmp2s+1  ! Step Lengh, from nstepmp2s to nstepmp2f

     ! Initial orbmp2k331
     call initialOrbmp2k331(orbmp2k331,nstep,nbasis,ivir,iocc,nsteplength)

     !---------------- MPI/ ALL NODES -----------------------------
     do i=1,mpi_jshelln(mpirank)
        II=mpi_jshell(mpirank,i)
        !     Do II=1,jshell
        do JJ=II,jshell

           ! First we do integral and sum them properly
           if(Ycutoff(II,JJ).gt.cutoffmp2/ttt)then
              call initialOrbmp2ij(orbmp2i331,nstep,nsteplength,nbasis,nbasistemp,nbasistemp)
              call initialOrbmp2ij(orbmp2j331,nstep,nsteplength,ivir,nbasistemp,nbasistemp)

              do KK=1,jshell
                 do LL=KK,jshell

                    ! Schwarts cutoff is implemented here
                    comax=0.d0
                    testCutoff = Ycutoff(II,JJ)*Ycutoff(KK,LL)
                    If(testCutoff.gt.cutoffmp2)then

                       NKK1=Qstart(KK)
                       NKK2=Qfinal(KK)
                       NLL1=Qstart(LL)
                       NLL2=Qfinal(LL)

                       NBK1=Qsbasis(KK,NKK1)
                       NBK2=Qfbasis(KK,NKK2)
                       NBL1=Qsbasis(LL,NLL1)
                       NBL2=Qfbasis(LL,NLL2)

                       KK111=Ksumtype(KK)+NBK1
                       KK112=Ksumtype(KK)+NBK2
                       LL111=Ksumtype(LL)+NBL1
                       LL112=Ksumtype(LL)+NBL2

                       ! find the co-max value
                       do icycle=1,nsteplength
                            i3=nstepmp2s+icycle-1
                            Do KKK=KK111,KK112
                                Do LLL=max(KKK,LL111),LL112
                                    comax=max(comax,dabs(co(kkk,i3)))
                                    comax=max(comax,dabs(co(lll,i3)))    
                                Enddo
                            Enddo        
                       enddo
                       
                       

                       testCutoff=testCutoff*comax
                       If(testCutoff.gt.cutoffmp2)then
                          dnmax=comax
                          ntemp=ntemp+1
                          call shellmp2(nstepmp2s,nsteplength)
                       Endif

                    endif

                 enddo
              enddo


              ! Next step is to folding integers. Without folding, the scaling will 
              ! be n^8, and after folding, scaling is reduced to 4n^5

              NII1=Qstart(II)
              NII2=Qfinal(II)
              NJJ1=Qstart(JJ)
              NJJ2=Qfinal(JJ)

              NBI1=Qsbasis(II,NII1)
              NBI2=Qfbasis(II,NII2)
              NBJ1=Qsbasis(JJ,NJJ1)
              NBJ2=Qfbasis(JJ,NJJ2)

              II111=Ksumtype(II)+NBI1
              II112=Ksumtype(II)+NBI2
              JJ111=Ksumtype(JJ)+NBJ1
              JJ112=Ksumtype(JJ)+NBJ2

              Do III=II111,II112
                 Do JJJ=max(III,JJ111),JJ112

                    IIInew=III-II111+1
                    JJJnew=JJJ-JJ111+1

                    ! Folding 4 indices integers into 3 indices
                    do LLL=1,nbasis
                       do j33=1,ivir
                          j33new=j33+nelec/2
                          atemp=CO(LLL,j33new)
                          do icycle=1,nsteplength
                             orbmp2j331(icycle,j33,IIInew,JJJnew,1)=orbmp2j331(icycle,j33,IIInew,JJJnew,1) + &
                                  orbmp2i331(icycle,LLL,IIInew,JJJnew,1)*atemp
                             if(III.ne.JJJ)then
                                orbmp2j331(icycle,j33,JJJnew,IIInew,2)=orbmp2j331(icycle,j33,JJJnew,IIInew,2) + &
                                     orbmp2i331(icycle,LLL,JJJnew,IIInew,2)*atemp
                             endif
                          enddo
                       enddo
                    enddo

                    do j33=1,ivir
                       do k33=1,nelec/2
                          atemp=HOLD(k33,III)
                          atemp2=HOLD(k33,JJJ)
                          do icycle=1,nsteplength
                             orbmp2k331(icycle,k33,j33,JJJ)=orbmp2k331(icycle,k33,j33,JJJ)+ &
                                  orbmp2j331(icycle,j33,IIInew,JJJnew,1)*atemp
                             if(III.ne.JJJ)then
                                orbmp2k331(icycle,k33,j33,III)=orbmp2k331(icycle,k33,j33,III)+ &
                                     orbmp2j331(icycle,j33,JJJnew,IIInew,2)*atemp2
                             endif
                          enddo
                       enddo
                    enddo
                 Enddo
              Enddo
           endif

        enddo
     enddo

     call cpu_time(timer_end%TMP2)
     timer_cumer%TMP2=timer_end%TMP2-timer_begin%TMP2+timer_cumer%TMP2

     ! Send the integral package after folding to master node
     ! and master node will receive them and do next two folding steps

     ! slave node will send infos
     if(.not.master) then
        tempE=emp2
        do j1=1,nbasis
           do k1=1,ivir
              do i1=1,iocc
                 do icycle=1,nsteplength
                    temp4d(icycle,i1,k1,j1)= orbmp2k331(icycle,i1,k1,j1)  
                 enddo
              enddo
           enddo
        enddo
        ! send 3 indices integrals to master node
        call MPI_SEND(temp4d,nbasis*ivir*iocc*nsteplength,mpi_double_precision,0,mpirank,MPI_COMM_WORLD,IERROR)
        ! master node will receive infos from every nodes
     else
        do i=1,mpisize-1
           ! receive integrals from slave nodes
           call MPI_RECV(temp4d,nbasis*ivir*iocc*nsteplength,mpi_double_precision,i,i,MPI_COMM_WORLD,MPI_STATUS,IERROR)
           ! and sum them into operator
           do j1=1,nbasis
              do k1=1,ivir
                 do i1=1,iocc
                    do icycle=1,nsteplength
                       orbmp2k331(icycle,i1,k1,j1)= orbmp2k331(icycle,i1,k1,j1) +temp4d(icycle,i1,k1,j1)
                    enddo
                 enddo
              enddo
           enddo
        enddo
     endif

     !---------------- END MPI/ ALL NODES ------------------------

     !---------------- MPI/ MASTER -------------------------------
     if (master) then


        do icycle=1,nsteplength
           i3=nstepmp2s+icycle-1
           Do k3=i3,nelec/2
              ! fold 3 indices integral into 2 indices
              Do J3=1,nbasis-nelec/2
                 Do L3=1,nbasis-nelec/2
                    orbmp2(L3,J3)=0.0d0
                    L3new=L3+nelec/2
                    Do LLL=1,nbasis
                       orbmp2(L3,J3)=orbmp2(L3,J3)+orbmp2k331(icycle,k3,j3,LLL)*CO(LLL,L3new)
                    Enddo
                 Enddo
              Enddo

              ! Now we can get energy
              Do J3=1,nbasis-nelec/2
                 Do L3=1,nbasis-nelec/2
                    if(k3.gt.i3)then
                       emp2=emp2+2.0d0/(E(i3)+E(k3)-E(j3+nelec/2)-E(l3+nelec/2)) &
                            *orbmp2(j3,l3)*(2.0d0*orbmp2(j3,l3)-orbmp2(l3,j3))
                    endif
                    if(k3.eq.i3)then
                       emp2=emp2+1.0d0/(E(i3)+E(k3)-E(j3+nelec/2)-E(l3+nelec/2)) &
                            *orbmp2(j3,l3)*(2.0d0*orbmp2(j3,l3)-orbmp2(l3,j3))
                    endif

                 Enddo
              Enddo
           enddo
        Enddo
     endif
     !---------------- ALL MPI/ MASTER ---------------------------

     ! sync all nodes
     call MPI_BARRIER(MPI_COMM_WORLD,mpierror)
  enddo

  if (master) then
     write (iOutFile,'("Second order energy=",F16.9)') EMP2
     write (iOutFile,'("EMP2               =",F16.9)') Etot+EMP2
     call PrtAct(ioutfile,"End MP2 Calculation")
  endif
  return
end subroutine mpi_calmp2


! Ed Brothers. November 27, 2001
! Xiao HE. September 14,2008
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

subroutine calmp2divcon
  use allmod
  use quick_gaussian_class_module
  implicit double precision(a-h,o-z)

  logical locallog1,locallog2

  real*8 Xiaotest,testtmp
  integer II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2
  common /hrrstore/II,JJ,KK,LL,NBI1,NBI2,NBJ1,NBJ2,NBK1,NBK2,NBL1,NBL2

  is = 0
  ishellfirst(1) = 1
  do ii=1,natom-1
     is=is+kshell(iattype(ii))
     ishelllast(ii)=is
     ishellfirst(ii+1) = is+1
  enddo
  ishelllast(natom) = nshell

  allocate(wtospoint(np,nbasis))
  call wtoscorr

  xiaocutoffmp2=1.0d-7

  emp2=0.0d0
  emp2temp=0.0d0
  nstep=1
  

  do itt=1,np

     Do i=1,nbasisdc(itt)
        Do j=1,nbasisdc(itt)
           CO(i,j)=COdcsub(i,j,itt)
        ENDDO
     ENDDO

     Do i=1,nbasisdc(itt)
        Do j=1,nbasisdc(itt)
           HOLD(i,j)=CO(j,i)
        ENDDO
     ENDDO


        do k3=1,2
            do j3=1,ivir
            do l3=1,ivir
                L3new=j3+iocc
                do lll=1,nbasisdc(itt)
                 write(10,*) k3,j3,LLL,L3new,orbmp2k331(1,k3,j3,LLL),co(LLL,L3new)
                enddo
            enddo
            enddo
        enddo
        
     ttt=0.0d0

     do iiat=1,dcsubn(itt)
        iiatom=dcsub(itt,iiat)
        do II=ishellfirst(iiatom),ishelllast(iiatom)
           do jjat=iiat,dcsubn(itt)
              jjatom=dcsub(itt,jjat)
              do JJ=max(ishellfirst(jjatom),II),ishelllast(jjatom)
                 Testtmp=Ycutoff(II,JJ)
                 ttt=max(ttt,Testtmp)
              enddo
           enddo
        enddo
     Enddo


    ! determine the electrons for subsystem
     if(mod(nelecmp2sub(itt),2).eq.1)then
        nelecmp2sub(itt)=nelecmp2sub(itt)+1
     endif

    ! determine the occupied and virtual orbitals
     if(mod(nelecmp2sub(itt),2).eq.1)then
        iocc=Nelecmp2sub(itt)/2+1
        ivir=nbasisdc(itt)-iocc+1
     else
        iocc=Nelecmp2sub(itt)/2
        ivir=nbasisdc(itt)-iocc
     endif
write(*,*) iocc,ivir
     ! with f orbital
     if (quick_method%ffunxiao) then
        nbasistemp=6
     else
        nbasistemp=10
     endif
     
     ! allocate varibles
     allocate(orbmp2(100,100))
     allocate(orbmp2dcsub(iocc,ivir,ivir))
     allocate(orbmp2i331(nstep,nbasisdc(itt),nbasistemp,nbasistemp,2))
     allocate(orbmp2j331(nstep,ivir,nbasistemp,nbasistemp,2))
     allocate(orbmp2k331(nstep,iocc,ivir,nbasisdc(itt)))
     allocate(orbmp2k331dcsub(iocc,ivir,nbasisdc(itt)))

     ! Schwartz cutoff is implemented here. (ab|cd)**2<=(ab|ab)*(cd|cd)
     ! Reference: Strout DL and Scuseria JCP 102(1995),8448.

     DO i3=1,iocc
        do l1=1,ivir
           do k1=1,ivir
                 orbmp2(j1,k1)=0.0d0
           enddo
        enddo

        call initialOrbmp2k331(orbmp2k331,nstep,nbasisdc(itt),ivir,iocc,nstep)

        do l1=1,ivir
           do k1=1,ivir
              do j1=1,iocc
                 orbmp2dcsub(j1,k1,l1)=0.0d0
              enddo
           enddo
        enddo

        do j1=1,nbasisdc(itt)
           do k1=1,ivir
              do i1=1,iocc
                 orbmp2k331dcsub(i1,k1,j1)=0.0d0
              enddo
           enddo
        enddo

        ntemp=0

        do iiat=1,dcsubn(itt)
           iiatom=dcsub(itt,iiat)
           
           do jjat=iiat,dcsubn(itt)
              jjatom=dcsub(itt,jjat)
              
              ! iiatom and jjatom is the atom of the subsystem
              ! first, we need to figure out which shell should we calclulate
              ! which is IIstart1 and IIstart2 for I and JJ is from JJstart1 to JJstart2
              if(iiatom.le.jjatom)then
                 IIstart1=ishellfirst(iiatom)
                 IIstart2=ishelllast(iiatom)
                 JJstart1=ishellfirst(jjatom)
                 JJstart2=ishelllast(jjatom)
              endif
              if(iiatom.gt.jjatom)then
                 JJstart1=ishellfirst(iiatom)
                 JJstart2=ishelllast(iiatom)
                 IIstart1=ishellfirst(jjatom)
                 IIstart2=ishelllast(jjatom)
              endif
              
              do II=IIstart1,IIstart2
                 do JJ=max(JJstart1,II),JJstart2

                    Testtmp=Ycutoff(II,JJ)
                    if(Testtmp.gt.XiaoCUTOFFmp2/ttt)then

                        call initialOrbmp2ij(orbmp2i331,nstep,nstep,nbasisdc(itt),nbasistemp,nbasistemp)
                        call initialOrbmp2ij(orbmp2j331,nstep,nstep,ivir,nbasistemp,nbasistemp)
                        
                       ! Now we will determine K shell and L shell, the last two indices
                       do kkat=1,dcsubn(itt)
                          kkatom=dcsub(itt,kkat)
                          do LLat=kkat,dcsubn(itt)
                             LLatom=dcsub(itt,LLat)
                             if(KKatom.le.LLatom)then
                                KKstart1=ishellfirst(KKatom)
                                KKstart2=ishelllast(KKatom)
                                LLstart1=ishellfirst(LLatom)
                                LLstart2=ishelllast(LLatom)
                             endif
                             if(KKatom.gt.LLatom)then
                                LLstart1=ishellfirst(KKatom)
                                LLstart2=ishelllast(KKatom)
                                KKstart1=ishellfirst(LLatom)
                                KKstart2=ishelllast(LLatom)
                             endif
                             do KK=KKstart1,KKstart2
                                do LL=max(LLstart1,KK),LLstart2
                                   comax=0.d0
                                   XiaoTEST1 = Ycutoff(II,JJ)*Ycutoff(KK,LL)
                                   If(XiaoTEST1.gt.XiaoCUTOFFmp2)then

                                      NKK1=Qstart(KK)
                                      NKK2=Qfinal(KK)
                                      NLL1=Qstart(LL)
                                      NLL2=Qfinal(LL)

                                      NBK1=Qsbasis(KK,NKK1)
                                      NBK2=Qfbasis(KK,NKK2)
                                      NBL1=Qsbasis(LL,NLL1)
                                      NBL2=Qfbasis(LL,NLL2)

                                      KK111=Ksumtype(KK)+NBK1
                                      KK112=Ksumtype(KK)+NBK2
                                      LL111=Ksumtype(LL)+NBL1
                                      LL112=Ksumtype(LL)+NBL2

                                      Do KKK=KK111,KK112
                                         Do LLL=max(KKK,LL111),LL112

                                            comax=max(comax,dabs(co(wtospoint(itt,kkk),i3)))
                                            comax=max(comax,dabs(co(wtospoint(itt,lll),i3)))    

                                         Enddo
                                      Enddo

                                      Xiaotest=xiaotest1*comax
                                      If(XiaoTEST.gt.XiaoCUTOFFmp2)then
                                         ntemp=ntemp+1
                                         call shellmp2divcon(i3,itt)
                                      Endif

                                   endif

                                enddo
                             enddo

                          ENDDO
                       ENDDO


                       NII1=Qstart(II)
                       NII2=Qfinal(II)
                       NJJ1=Qstart(JJ)
                       NJJ2=Qfinal(JJ)

                       NBI1=Qsbasis(II,NII1)
                       NBI2=Qfbasis(II,NII2)
                       NBJ1=Qsbasis(JJ,NJJ1)
                       NBJ2=Qfbasis(JJ,NJJ2)

                       II111=Ksumtype(II)+NBI1
                       II112=Ksumtype(II)+NBI2
                       JJ111=Ksumtype(JJ)+NBJ1
                       JJ112=Ksumtype(JJ)+NBJ2

                       Do III=II111,II112
                          Do JJJ=max(III,JJ111),JJ112

                             IIInew=III-II111+1
                             JJJnew=JJJ-JJ111+1

                             do LLL=1,nbasisdc(itt)
                                do j33=1,ivir
                                   j33new=j33+iocc
                                   if(mod(nelecmp2sub(itt),2).eq.1)j33new=j33+iocc-1
                                   atemp=CO(LLL,j33new)
                                   orbmp2j331(nstep,j33,IIInew,JJJnew,1)=orbmp2j331(nstep,j33,IIInew,JJJnew,1) + &
                                        orbmp2i331(nstep,LLL,IIInew,JJJnew,1)*atemp
                                   if(III.ne.JJJ)then
                                      orbmp2j331(nstep,j33,JJJnew,IIInew,2)=orbmp2j331(nstep,j33,JJJnew,IIInew,2) + &
                                           orbmp2i331(nstep,LLL,JJJnew,IIInew,2)*atemp
                                   endif
                                enddo
                             enddo

                             do j33=1,ivir
                                do k33=1,iocc
                                   orbmp2k331(nstep,k33,j33,wtospoint(itt,JJJ))=orbmp2k331(nstep,k33,j33,wtospoint(itt,JJJ))+ &
                                        orbmp2j331(nstep,j33,IIInew,JJJnew,1)*HOLD(k33,wtospoint(itt,III))
                                   if(III.ne.JJJ)then
                                      orbmp2k331(nstep,k33,j33,wtospoint(itt,III))=orbmp2k331(nstep,k33,j33,wtospoint(itt,III))+ &
                                           orbmp2j331(nstep,j33,JJJnew,IIInew,2)*HOLD(k33,wtospoint(itt,JJJ))
                                   endif
                                enddo
                             enddo

                             locallog1=.false.
                             locallog2=.false.

                             do iiatdc=1,dccoren(itt)
                                iiatomdc=dccore(itt,iiatdc)
                                do IInbasisdc=ifirst(iiatomdc),ilast(iiatomdc)
                                   if(III.eq.IInbasisdc)locallog1=.true.
                                   if(JJJ.eq.IInbasisdc)locallog2=.true.
                                enddo
                             enddo

                             if(locallog1)then
                                do j33=1,ivir
                                   do k33=1,iocc
                                      orbmp2k331dcsub(k33,j33,wtospoint(itt,JJJ))=orbmp2k331dcsub(k33,j33,wtospoint(itt,JJJ))+ &
                                           orbmp2j331(nstep,j33,IIInew,JJJnew,1)*HOLD(k33,wtospoint(itt,III))
                                   enddo
                                enddo
                             endif

                             if(locallog2.and.III.ne.JJJ)then
                                do j33=1,ivir
                                   do k33=1,iocc
                                      orbmp2k331dcsub(k33,j33,wtospoint(itt,III))=orbmp2k331dcsub(k33,j33,wtospoint(itt,III))+ &
                                           orbmp2j331(nstep,j33,JJJnew,IIInew,2)*HOLD(k33,wtospoint(itt,JJJ))
                                   enddo
                                enddo
                             endif

                          Enddo
                       Enddo

                    endif

                 enddo
              enddo

           ENDDO
        ENDDO
        write (ioutfile,*)"ntemp=",ntemp
        
        do k3=1,2
            do j3=1,ivir
            do l3=1,ivir
                L3new=j3+iocc
                do lll=1,nbasisdc(itt)
                 write(10,*) k3,j3,LLL,L3new,orbmp2k331(1,k3,j3,LLL),co(LLL,L3new)
                enddo
            enddo
            enddo
        enddo
        
        Do LLL=1,nbasisdc(itt)
           Do J3=1,ivir
              Do L3=1,ivir
                 L3new=L3+iocc
                 orbmp2(L3,j3)=0.0d0
                 if(mod(nelecmp2sub(itt),2).eq.1)L3new=L3+iocc-1
                 Do k3=1,iocc
                    orbmp2(l3,j3)=orbmp2(l3,j3)+orbmp2k331(nstep,k3,j3,LLL)*HOLD(L3new,LLL)
                    orbmp2dcsub(k3,l3,j3)=orbmp2dcsub(k3,l3,j3)+orbmp2k331dcsub(k3,j3,LLL)*HOLD(L3new,LLL)
                 Enddo
              Enddo
           Enddo
        Enddo

        if(mod(nelecmp2sub(itt),2).eq.0)then
           do l=1,ivir
              do k=1,iocc
                 do j=1,ivir
                    emp2=emp2+1.0d0/(Evaldcsub(itt,i3)+Evaldcsub(itt,k) &
                         -Evaldcsub(itt,j+nelecmp2sub(itt)/2)-Evaldcsub(itt,l+nelecmp2sub(itt)/2)) &
                         *orbmp2dcsub(k,j,l)*(2.0d0*orbmp2(j,l)-orbmp2(l,j))
                 enddo
              enddo
           enddo
        endif

        if(mod(nelecmp2sub(itt),2).eq.1)then
           do l=1,ivir
              do k=1,iocc
                 do j=1,ivir
                    emp2=emp2+1.0d0/(Evaldcsub(itt,i3)+Evaldcsub(itt,k) &
                         -Evaldcsub(itt,j+nelecmp2sub(itt)/2)-Evaldcsub(itt,l+nelecmp2sub(itt)/2)) &
                         *orbmp2dcsub(k,j,l)*(2.0d0*orbmp2(j,l)-orbmp2(l,j))
                 enddo
              enddo
           enddo
        endif
     ENDDO
     write(ioutfile,*) itt,emp2,emp2-emp2temp

     emp2temp=emp2

!     deallocate(mp2shell)
     deallocate(orbmp2)
     deallocate(orbmp2i331)
     deallocate(orbmp2j331)
     deallocate(orbmp2k331)
     deallocate(orbmp2dcsub)
     deallocate(orbmp2k331dcsub)

  enddo

999 return
end subroutine calmp2divcon


! Initial Orbmp2k331
subroutine  initialOrbmp2k331(orbmp2k331,nstep,nbasis,ivir,iocc,nsteplength)
  integer nbasis,ivir,iocc,nsteplength,i1,jk1,j1,icycle,nstep
  double precision orbmp2k331(nstep,iocc,ivir,nbasis)
  do j1=1,nbasis
     do k1=1,ivir
        do i1=1,iocc
           do icycle=1,nsteplength
              orbmp2k331(icycle,i1,k1,j1)=0.0d0  
           enddo
        enddo
     enddo
  enddo
end subroutine initialOrbmp2k331

! Initial Orbmp2i331 and Orbmp2j331
subroutine initialOrbmp2ij(orbmp2i331,nstep,nsteplength,nbasis,nbasistemp,nbasistemp2)
  integer nstep,nsteplength,nbasis,nbasistemp,nbasistemp2
  integer l1,j1,i1,k1,icycle
  double precision orbmp2i331(nstep,nbasis,nbasistemp,nbasistemp2,2)

  do l1=1,2
     do j1=1,nbasistemp
        do i1=1,nbasistemp2
           do k1=1,nbasis
              do icycle=1,nsteplength
                 orbmp2i331(icycle,k1,i1,j1,l1)=0.0d0
              enddo
           enddo
        enddo
     enddo
  enddo
end subroutine initialOrbmp2ij
