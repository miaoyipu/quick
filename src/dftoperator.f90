! Ed Brothers. November 27, 2001
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

    subroutine dftoperator
    use allmod
    implicit double precision(a-h,o-z)

! The purpose of this subroutine is to form the operator matrix
! for a full Density Functional calculation, i.e. the KS matrix.  The
! KS  matrix is as follows:

! O(I,J) =  F(I,J) = KE(I,J) + IJ attraction to each atom + repulsion_prim
! with each possible basis  + Exchange/Correlation functional.

! Note that the KS operator matrix is symmetric.

! The first part is the one elctron code.

    do Ibas=1,nbasis
        do Jbas=Ibas,nbasis
            O(Jbas,Ibas) = 0.d0
            do Icon=1,ncontract(ibas)
                do Jcon=1,ncontract(jbas)

                    O(Jbas,Ibas)=O(Jbas,Ibas)+ &
                    dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas)* &
                    ekinetic(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                    itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                    itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                    xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                    xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                    xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))

!                    do iatom = 1,natom
!                        O(Jbas,Ibas)=O(Jbas,Ibas)+ &
!                        dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas)* &
!                        attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
!                        itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
!                        itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
!                        xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
!                        xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
!                        xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)), &
!                        xyz(1,iatom),xyz(2,iatom),xyz(3,iatom), &
!                        chg(iatom))
!                    enddo
                enddo
            enddo
        enddo
    enddo

    do IIsh=1,jshell
      do JJsh=IIsh,jshell
        call attrashell(IIsh,JJsh)
      enddo
    enddo

  call cpu_time(t1)

    Eelxc=0.0d0

  if(quick_method%printEnergy)then
    Eel=0.d0
    do Ibas=1,nbasis
        do Icon=1,ncontract(Ibas)
            do Jcon=1,ncontract(Ibas)

            ! Kinetic energy.

                Eel=Eel+DENSE(Ibas,Ibas)* &
                dcoeff(Jcon,Ibas)*dcoeff(Icon,Ibas)* &
                ekinetic(aexp(Jcon,Ibas),aexp(Icon,Ibas), &
                itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                xyz(1,ncenter(Ibas)),xyz(2,ncenter(Ibas)), &
                xyz(3,ncenter(Ibas)),xyz(1,ncenter(Ibas)), &
                xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))

            ! Nuclear attraction.

!                do iatom = 1,natom
!                    Eel=Eel+DENSE(Ibas,Ibas)* &
!                    dcoeff(Jcon,Ibas)*dcoeff(Icon,Ibas)* &
!                    attraction(aexp(Jcon,Ibas),aexp(Icon,Ibas), &
!                    itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
!                    itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
!                    xyz(1,ncenter(Ibas)),xyz(2,ncenter(Ibas)), &
!                    xyz(3,ncenter(Ibas)),xyz(1,ncenter(Ibas)), &
!                    xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)), &
!                    xyz(1,iatom),xyz(2,iatom),xyz(3,iatom), &
!                    chg(iatom))
!                enddo
            enddo
        enddo
    enddo

    do Ibas=1,nbasis
        do Jbas=Ibas+1,nbasis
            do Icon=1,ncontract(ibas)
                do Jcon=1,ncontract(jbas)

                ! Kinetic energy.

                    Eel=Eel+DENSE(Jbas,Ibas)* &
                    dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas)* &
                    2.d0*ekinetic(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                    itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                    itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                    xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                    xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                    xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))

                ! Nuclear attraction.

!                    do iatom = 1,natom
!                        Eel=Eel+DENSE(Jbas,Ibas)* &
!                        dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas)* &
!                        2.d0*attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
!                        itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
!                        itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
!                        xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
!                        xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
!                        xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)), &
!                        xyz(1,iatom),xyz(2,iatom),xyz(3,iatom), &
!                        chg(iatom))
!                    enddo
                enddo
            enddo
        enddo
    enddo

    do IIsh=1,jshell
      do JJsh=IIsh,jshell
        call attrashellenergy(IIsh,JJsh)
      enddo
    enddo

  endif
!    print*,'Eel=', Eel

  call cpu_time(t2)

         if(quick_method%printEnergy)then
           write (ioutfile,'("Time for one-electron energy evaluation=",F16.9)') t2-t1
         endif

!


!
! Alessandro GENONI 03/21/2007
! Sum the ECP integrals to the partial Fock matrix
!
    if (quick_method%ecp) then
      call ecpoperator
    end if

! The previous two terms are the one electron part of the Fock matrix.
! The next term defines the electron repulsion_prim.

! Delta density matrix cutoff

 call cpu_time(T1)

 do II=1,jshell
   do JJ=II,jshell
     DNtemp=0.0d0
     call DNscreen(II,JJ,DNtemp)
     Cutmatrix(II,JJ)=DNtemp
     Cutmatrix(JJ,II)=DNtemp
   enddo
 enddo

! Schwartz cutoff is implemented here. (ab|cd)**2<=(ab|ab)*(cd|cd)
! Reference: Strout DL and Scuseria JCP 102(1995),8448.

! print*,"before 2e"
if(quick_method%B3LYP)then
 do II=1,jshell
   do JJ=II,jshell
     Testtmp=Ycutoff(II,JJ)
     do KK=II,jshell
       do LL=KK,jshell
!          Nxiao1=Nxiao1+1
            cutoffTest1 = TESTtmp*Ycutoff(KK,LL)
          If(cutoffTest1.gt.quick_method%integralCutoff)then
            DNmax=max(4.0d0*cutmatrix(II,JJ),4.0d0*cutmatrix(KK,LL), &
                  cutmatrix(II,LL),cutmatrix(II,KK),cutmatrix(JJ,KK),cutmatrix(JJ,LL))
!            DNmax=max(cutmatrix(II,JJ),cutmatrix(KK,LL) &
!                  )
            cutoffTest=testCutoff*DNmax
            If(cutoffTest.gt.quick_method%integralCutoff)then
              IIxiao=II
              JJxiao=JJ
              KKxiao=KK
              LLxiao=LL
              call shelldftb3lyp(IIxiao,JJxiao,KKxiao,LLxiao)
!            Nxiao2=Nxiao2+1
            endif
!            else
!             print*,II,JJ,KK,LL,cutoffTest,testCutoff,DNmax
!            print*,'***',O(1,1)
          endif
       enddo
     enddo
   enddo
 enddo
else
 do II=1,jshell
   do JJ=II,jshell
     Testtmp=Ycutoff(II,JJ)
     do KK=II,jshell
       do LL=KK,jshell
!          Nxiao1=Nxiao1+1
            cutoffTest1 = TESTtmp*Ycutoff(KK,LL)
          If(cutoffTest1.gt.quick_method%integralCutoff)then
!            DNmax=max(4.0d0*cutmatrix(II,JJ),4.0d0*cutmatrix(KK,LL), &
!                  cutmatrix(II,LL),cutmatrix(II,KK),cutmatrix(JJ,KK),cutmatrix(JJ,LL))
            DNmax=max(cutmatrix(II,JJ),cutmatrix(KK,LL) &
                  )
            cutoffTest=testCutoff*DNmax
            If(cutoffTest.gt.quick_method%integralCutoff)then
              IIxiao=II
              JJxiao=JJ
              KKxiao=KK
              LLxiao=LL
              call shelldft(IIxiao,JJxiao,KKxiao,LLxiao)
!            Nxiao2=Nxiao2+1
            endif
!            else
!             print*,II,JJ,KK,LL,cutoffTest,testCutoff,DNmax
!            print*,'***',O(1,1)
          endif
       enddo
     enddo
   enddo
 enddo
endif

        do I=1,nbasis
          do J=1,nbasis
            Osavedft(i,j)=O(i,j)
          enddo
        enddo

        call cpu_time(t2)

        write (ioutfile,'(" TIME of evaluation integral = ",F12.2)') &
        T2-T1

! print*,'Nxiao1=',Nxiao1,'Nxiao2=',Nxiao2,integralCutOff

    do Ibas=1,nbasis
        do Jbas=Ibas+1,nbasis
            O(Ibas,Jbas) = O(Jbas,Ibas)
        enddo
    enddo

   if(quick_method%printEnergy)then
    do Ibas=1,nbasis
        do Jbas=1,nbasis
            Eel=Eel+DENSE(Ibas,Jbas)*O(Jbas,Ibas)
        enddo
    enddo

    Eel=Eel/2.0d0
   endif


! stop

! The next portion is the exchange/correlation functional.
! The angular grid code came from CCL.net.  The radial grid
! formulas (position and wieghts) is from Gill, Johnson and Pople,
! Chem. Phys. Lett. v209, n 5+6, 1993, pg 506-512.  The weighting scheme
! is from Stratmann, Scuseria, and Frisch, Chem. Phys. Lett., v 257,
! 1996, pg 213-223.

! The actual element is:
! F alpha mu nu = Integral((df/drhoa Phimu Phinu)+
! (2 df/dgaa Grad(rho a) + df/dgab Grad(rho b)) doT Grad(Phimu Phinu))

! where F alpha mu nu is the the alpha spin portion of the operator matrix
! element mu, nu,
! df/drhoa is the derivative of the functional by the alpha density,
! df/dgaa is the derivative of the functional by the alpha gradient
! invariant, i.e. the dot product of the gradient of the alpha
! density with itself.
! df/dgab is the derivative of the functional by the dot product of
! the gradient of the alpha density with the beta density.
! Grad(Phimu Phinu) is the gradient of Phimu times Phinu.

! First, find the grid point.

    aelec=0.d0
    belec=0.d0

! Xiao HE 02/09/1007 SG0 SG1 grids
!    do Ireg=1,iregion
!        call gridform(iangular(Ireg))
!        do Irad=iradial(ireg-1)+1,iradial(ireg)
!            do Iatm=1,natom
!                rad = radii(quick_molspec%iattype(iatm))
!                rad3 = rad*rad*rad
!                do Iang=1,iangular(Ireg)

     if(quick_method%B3LYP)then
         do Iatm=1,natom
           If(quick_method%ISG.eq.1)then
             Iradtemp=50
           else
             If(quick_molspec%iattype(iatm).le.10)then
               Iradtemp=23
             else
               Iradtemp=26
             endif
           endif

           do Irad=1,Iradtemp
             If(quick_method%ISG.eq.1)then
               call gridformnew(iatm,RGRID(Irad),iiangt)
               rad = radii(quick_molspec%iattype(iatm))
             else
               call gridformSG0(iatm,Iradtemp+1-Irad,iiangt,RGRID,RWT)
               rad = radii2(quick_molspec%iattype(iatm))
             endif

             rad3 = rad*rad*rad
!             print*,iiangt
               do Iang=1,iiangt
                    gridx=xyz(1,Iatm)+rad*RGRID(Irad)*XANG(Iang)
                    gridy=xyz(2,Iatm)+rad*RGRID(Irad)*YANG(Iang)
                    gridz=xyz(3,Iatm)+rad*RGRID(Irad)*ZANG(Iang)

                ! Next, calculate the weight of the grid point in the SSW scheme.  If
                ! the grid point has a zero weight, we can skip it.

                    weight=SSW(gridx,gridy,gridz,Iatm) &
                    *WTANG(Iang)*RWT(Irad)*rad3

                    if (weight < quick_method%DMCutoff ) then
                        continue
                    else

                            do Ibas=1,nbasis
                                call pteval(gridx,gridy,gridz,phi,dphidx,dphidy, &
                                dphidz,Ibas)
                                phixiao(Ibas)=phi
                                dphidxxiao(Ibas)=dphidx
                                dphidyxiao(Ibas)=dphidy
                                dphidzxiao(Ibas)=dphidz
                            enddo
                            
                    ! Next, evaluate the densities at the grid point and the gradient
                    ! at that grid point.

                        call denspt(gridx,gridy,gridz,density,densityb,gax,gay,gaz, &
                        gbx,gby,gbz)

                        if (density < quick_method%DMCutoff ) then
                            continue
                        else

                     densitysum=2.0d0*density
                     sigma=4.0d0*(gax*gax+gay*gay+gaz*gaz)
                     call b3lyp_e(densitysum,sigma,zkec)
!
!                    Eel = Eel + (param7*Ex+param8*Ec) &
                     Eelxc = Eelxc + zkec*weight


                            aelec = weight*density+aelec
                            belec = weight*densityb+belec

                        ! This allows the calculation of the derivative of the functional
                        ! with regard to the density (dfdr), with regard to the alpha-alpha
                        ! density invariant (df/dgaa), and the alpha-beta density invariant.

!                            call becke(density,gax,gay,gaz,gbx,gby,gbz, &
!                            dfdr,dfdgaa,dfdgab)
!
!                            call lyp(density,densityb,gax,gay,gaz,gbx,gby,gbz, &
!                            dfdr2,dfdgaa2,dfdgab2)

                            densitysum=2.0d0*density
                            sigma=4.0d0*(gax*gax+gay*gay+gaz*gaz)

                            call b3lypf(densitysum,sigma,dfdr,xiaodot)
                          
!                            dfdr=dfdr*2.0d0

!                            dfdr = dfdr+dfdr2
!                            dfdgaa = dfdgaa + dfdgaa2
!                            dfdgab = dfdgab + dfdgab2

                        ! Calculate the first term in the dot product shown above,i.e.:
                        ! (2 df/dgaa Grad(rho a) + df/dgab Grad(rho b)) doT Grad(Phimu Phinu))
                             xdot=xiaodot*gax
                             ydot=xiaodot*gay
                             zdot=xiaodot*gaz

!                            xdot=2.d0*dfdgaa*gax+dfdgab*gbx
!                            ydot=2.d0*dfdgaa*gay+dfdgab*gby
!                            zdot=2.d0*dfdgaa*gaz+dfdgab*gbz

                        ! Now loop over basis functions and compute the addition to the matrix
                        ! element.

                            do Ibas=1,nbasis
!                                call pteval(gridx,gridy,gridz,phi,dphidx,dphidy, &
!                                dphidz,Ibas)
        phi=phixiao(Ibas)
        dphidx=dphidxxiao(Ibas)
        dphidy=dphidyxiao(Ibas)
        dphidz=dphidzxiao(Ibas)
!                                quicktest = DABS(dphidx)+DABS(dphidy)+DABS(dphidz) &
!                                +DABS(phi)
                                quicktest = DABS(dphidx+dphidy+dphidz+ &
                                phi)
                                if (quicktest < quick_method%DMCutoff ) then
                                    continue
                                else
                                    do Jbas=Ibas,nbasis
!                                        call pteval(gridx,gridy,gridz,phi2,dphi2dx,dphi2dy, &
!                                        dphi2dz,Jbas)
        phi2=phixiao(Jbas)
        dphi2dx=dphidxxiao(Jbas)
        dphi2dy=dphidyxiao(Jbas)
        dphi2dz=dphidzxiao(Jbas)
!                                        quicktest = DABS(dphi2dx)+DABS(dphi2dy)+DABS(dphi2dz) &
!                                        +DABS(phi2)
!                                        if (quicktest < quick_method%DMCutoff ) then
!                                            continue
!                                        else
                                            temp = phi*phi2
                                            tempgx = phi*dphi2dx + phi2*dphidx
                                            tempgy = phi*dphi2dy + phi2*dphidy
                                            tempgz = phi*dphi2dz + phi2*dphidz
                                            O(Jbas,Ibas)=O(Jbas,Ibas)+(temp*dfdr+ &
                                            xdot*tempgx+ydot*tempgy+zdot*tempgz)*weight
!                                        endif
                                    enddo
                                endif
                            enddo
                        endif
                    endif
                enddo
            enddo
        enddo
!    enddo
     endif

    if(quick_method%BLYP)then
         do Iatm=1,natom
           If(quick_method%ISG.eq.1)then
             Iradtemp=50
           else
             If(quick_molspec%iattype(iatm).le.10)then
               Iradtemp=23
             else
               Iradtemp=26
             endif
           endif

           do Irad=1,Iradtemp
             If(quick_method%ISG.eq.1)then
               call gridformnew(iatm,RGRID(Irad),iiangt)
               rad = radii(quick_molspec%iattype(iatm))
             else
               call gridformSG0(iatm,Iradtemp+1-Irad,iiangt,RGRID,RWT)
               rad = radii2(quick_molspec%iattype(iatm))
             endif

             rad3 = rad*rad*rad
!             print*,iiangt
               do Iang=1,iiangt
                    gridx=xyz(1,Iatm)+rad*RGRID(Irad)*XANG(Iang)
                    gridy=xyz(2,Iatm)+rad*RGRID(Irad)*YANG(Iang)
                    gridz=xyz(3,Iatm)+rad*RGRID(Irad)*ZANG(Iang)

                ! Next, calculate the weight of the grid point in the SSW scheme.  If
                ! the grid point has a zero weight, we can skip it.

                    weight=SSW(gridx,gridy,gridz,Iatm) &
                    *WTANG(Iang)*RWT(Irad)*rad3

                    if (weight < quick_method%DMCutoff ) then
                        continue
                    else

                            do Ibas=1,nbasis
                                call pteval(gridx,gridy,gridz,phi,dphidx,dphidy, &
                                dphidz,Ibas)
                                phixiao(Ibas)=phi
                                dphidxxiao(Ibas)=dphidx
                                dphidyxiao(Ibas)=dphidy
                                dphidzxiao(Ibas)=dphidz
                            enddo
                            
                    ! Next, evaluate the densities at the grid point and the gradient
                    ! at that grid point.

                        call denspt(gridx,gridy,gridz,density,densityb,gax,gay,gaz, &
                        gbx,gby,gbz)

                        if (density < quick_method%DMCutoff ) then
                            continue
                        else

                    call becke_E(density,densityb,gax,gay,gaz,gbx,gby,gbz,Ex)

                    call lyp_e(density,densityb,gax,gay,gaz,gbx,gby,gbz,Ec)

                    Eelxc = Eelxc + (param7*Ex+param8*Ec) &
                            *weight


                            aelec = weight*density+aelec
                            belec = weight*densityb+belec

                        ! This allows the calculation of the derivative of the functional
                        ! with regard to the density (dfdr), with regard to the alpha-alpha
                        ! density invariant (df/dgaa), and the alpha-beta density invariant.

                            call becke(density,gax,gay,gaz,gbx,gby,gbz, &
                            dfdr,dfdgaa,dfdgab)

                            call lyp(density,densityb,gax,gay,gaz,gbx,gby,gbz, &
                            dfdr2,dfdgaa2,dfdgab2)

!                            dfdr = param7*dfdr+param8*dfdr2
!                            dfdgaa = param7*dfdgaa + param8*dfdgaa2
!                            dfdgab = param7*dfdgab + param8*dfdgab2

                            dfdr = dfdr+dfdr2
                            dfdgaa = dfdgaa + dfdgaa2
                            dfdgab = dfdgab + dfdgab2

                        ! Calculate the first term in the dot product shown above,i.e.:
                        ! (2 df/dgaa Grad(rho a) + df/dgab Grad(rho b)) doT Grad(Phimu Phinu))

                            xdot=2.d0*dfdgaa*gax+dfdgab*gbx
                            ydot=2.d0*dfdgaa*gay+dfdgab*gby
                            zdot=2.d0*dfdgaa*gaz+dfdgab*gbz

                        ! Now loop over basis functions and compute the addition to the matrix
                        ! element.

                            do Ibas=1,nbasis
!                                call pteval(gridx,gridy,gridz,phi,dphidx,dphidy, &
!                                dphidz,Ibas)
        phi=phixiao(Ibas)
        dphidx=dphidxxiao(Ibas)
        dphidy=dphidyxiao(Ibas)
        dphidz=dphidzxiao(Ibas)
!                                quicktest = DABS(dphidx)+DABS(dphidy)+DABS(dphidz) &
!                                +DABS(phi)
                                quicktest = DABS(dphidx+dphidy+dphidz+ &
                                phi)
                                if (quicktest < quick_method%DMCutoff ) then
                                    continue
                                else
                                    do Jbas=Ibas,nbasis
!                                        call pteval(gridx,gridy,gridz,phi2,dphi2dx,dphi2dy, &
!                                        dphi2dz,Jbas)
        phi2=phixiao(Jbas)
        dphi2dx=dphidxxiao(Jbas)
        dphi2dy=dphidyxiao(Jbas)
        dphi2dz=dphidzxiao(Jbas)
!                                        quicktest = DABS(dphi2dx)+DABS(dphi2dy)+DABS(dphi2dz) &
!                                        +DABS(phi2)
!                                        if (quicktest < quick_method%DMCutoff ) then
!                                            continue
!                                        else
                                            temp = phi*phi2
                                            tempgx = phi*dphi2dx + phi2*dphidx
                                            tempgy = phi*dphi2dy + phi2*dphidy
                                            tempgz = phi*dphi2dz + phi2*dphidz
                                            O(Jbas,Ibas)=O(Jbas,Ibas)+(temp*dfdr+ &
                                            xdot*tempgx+ydot*tempgy+zdot*tempgz)*weight
!                                        endif
                                    enddo
                                endif
                            enddo
                        endif
                    endif
                enddo
            enddo
        enddo
     endif

! Finally, copy lower diagonal to upper diagonal.

    if(quick_method%MPW91LYP)then
         do Iatm=1,natom
           If(quick_method%ISG.eq.1)then
             Iradtemp=50
           else
             If(quick_molspec%iattype(iatm).le.10)then
               Iradtemp=23
             else
               Iradtemp=26
             endif
           endif

           do Irad=1,Iradtemp
             If(quick_method%ISG.eq.1)then
               call gridformnew(iatm,RGRID(Irad),iiangt)
               rad = radii(quick_molspec%iattype(iatm))
             else
               call gridformSG0(iatm,Iradtemp+1-Irad,iiangt,RGRID,RWT)
               rad = radii2(quick_molspec%iattype(iatm))
             endif

             rad3 = rad*rad*rad
!             print*,iiangt
               do Iang=1,iiangt
                    gridx=xyz(1,Iatm)+rad*RGRID(Irad)*XANG(Iang)
                    gridy=xyz(2,Iatm)+rad*RGRID(Irad)*YANG(Iang)
                    gridz=xyz(3,Iatm)+rad*RGRID(Irad)*ZANG(Iang)

                ! Next, calculate the weight of the grid point in the SSW scheme.  If
                ! the grid point has a zero weight, we can skip it.

                    weight=SSW(gridx,gridy,gridz,Iatm) &
                    *WTANG(Iang)*RWT(Irad)*rad3

                    if (weight < quick_method%DMCutoff ) then
                        continue
                    else

                            do Ibas=1,nbasis
                                call pteval(gridx,gridy,gridz,phi,dphidx,dphidy, &
                                dphidz,Ibas)
                                phixiao(Ibas)=phi
                                dphidxxiao(Ibas)=dphidx
                                dphidyxiao(Ibas)=dphidy
                                dphidzxiao(Ibas)=dphidz
                            enddo
                            
                    ! Next, evaluate the densities at the grid point and the gradient
                    ! at that grid point.

                        call denspt(gridx,gridy,gridz,density,densityb,gax,gay,gaz, &
                        gbx,gby,gbz)

                        if (density < quick_method%DMCutoff ) then
                            continue
                        else

                    call mpw91_E(density,densityb,gax,gay,gaz,gbx,gby,gbz,Ex)

                    call lyp_e(density,densityb,gax,gay,gaz,gbx,gby,gbz,Ec)

                    Eelxc = Eelxc + (param7*Ex+param8*Ec) &
                            *weight


                            aelec = weight*density+aelec
                            belec = weight*densityb+belec

                        ! This allows the calculation of the derivative of the functional
                        ! with regard to the density (dfdr), with regard to the alpha-alpha
                        ! density invariant (df/dgaa), and the alpha-beta density invariant.

                            call mpw91(density,gax,gay,gaz,gbx,gby,gbz, &
                            dfdr,dfdgaa,dfdgab)

                            call lyp(density,densityb,gax,gay,gaz,gbx,gby,gbz, &
                            dfdr2,dfdgaa2,dfdgab2)

!                            dfdr = param7*dfdr+param8*dfdr2
!                            dfdgaa = param7*dfdgaa + param8*dfdgaa2
!                            dfdgab = param7*dfdgab + param8*dfdgab2

                            dfdr = dfdr+dfdr2
                            dfdgaa = dfdgaa + dfdgaa2
                            dfdgab = dfdgab + dfdgab2

                        ! Calculate the first term in the dot product shown above,i.e.:
                        ! (2 df/dgaa Grad(rho a) + df/dgab Grad(rho b)) doT Grad(Phimu Phinu))

                            xdot=2.d0*dfdgaa*gax+dfdgab*gbx
                            ydot=2.d0*dfdgaa*gay+dfdgab*gby
                            zdot=2.d0*dfdgaa*gaz+dfdgab*gbz

                        ! Now loop over basis functions and compute the addition to the matrix
                        ! element.

                            do Ibas=1,nbasis
!                                call pteval(gridx,gridy,gridz,phi,dphidx,dphidy, &
!                                dphidz,Ibas)
        phi=phixiao(Ibas)
        dphidx=dphidxxiao(Ibas)
        dphidy=dphidyxiao(Ibas)
        dphidz=dphidzxiao(Ibas)
!                                quicktest = DABS(dphidx)+DABS(dphidy)+DABS(dphidz) &
!                                +DABS(phi)
                                quicktest = DABS(dphidx+dphidy+dphidz+ &
                                phi)
                                if (quicktest < quick_method%DMCutoff ) then
                                    continue
                                else
                                    do Jbas=Ibas,nbasis
!                                        call pteval(gridx,gridy,gridz,phi2,dphi2dx,dphi2dy, &
!                                        dphi2dz,Jbas)
        phi2=phixiao(Jbas)
        dphi2dx=dphidxxiao(Jbas)
        dphi2dy=dphidyxiao(Jbas)
        dphi2dz=dphidzxiao(Jbas)
!                                        quicktest = DABS(dphi2dx)+DABS(dphi2dy)+DABS(dphi2dz) &
!                                        +DABS(phi2)
!                                        if (quicktest < quick_method%DMCutoff ) then
!                                            continue
!                                        else
                                            temp = phi*phi2
                                            tempgx = phi*dphi2dx + phi2*dphidx
                                            tempgy = phi*dphi2dy + phi2*dphidy
                                            tempgz = phi*dphi2dz + phi2*dphidz
                                            O(Jbas,Ibas)=O(Jbas,Ibas)+(temp*dfdr+ &
                                            xdot*tempgx+ydot*tempgy+zdot*tempgz)*weight
!                                        endif
                                    enddo
                                endif
                            enddo
                        endif
                    endif
                enddo
            enddo
        enddo
     endif

    if(quick_method%BPW91)then
         do Iatm=1,natom
           If(quick_method%ISG.eq.1)then
             Iradtemp=50
           else
             If(quick_molspec%iattype(iatm).le.10)then
               Iradtemp=23
             else
               Iradtemp=26
             endif
           endif

           do Irad=1,Iradtemp
             If(quick_method%ISG.eq.1)then
               call gridformnew(iatm,RGRID(Irad),iiangt)
               rad = radii(quick_molspec%iattype(iatm))
             else
               call gridformSG0(iatm,Iradtemp+1-Irad,iiangt,RGRID,RWT)
               rad = radii2(quick_molspec%iattype(iatm))
             endif

             rad3 = rad*rad*rad
!             print*,iiangt
               do Iang=1,iiangt
                    gridx=xyz(1,Iatm)+rad*RGRID(Irad)*XANG(Iang)
                    gridy=xyz(2,Iatm)+rad*RGRID(Irad)*YANG(Iang)
                    gridz=xyz(3,Iatm)+rad*RGRID(Irad)*ZANG(Iang)

                ! Next, calculate the weight of the grid point in the SSW scheme.  If
                ! the grid point has a zero weight, we can skip it.

                    weight=SSW(gridx,gridy,gridz,Iatm) &
                    *WTANG(Iang)*RWT(Irad)*rad3

                    if (weight < quick_method%DMCutoff ) then
                        continue
                    else

                            do Ibas=1,nbasis
                                call pteval(gridx,gridy,gridz,phi,dphidx,dphidy, &
                                dphidz,Ibas)
                                phixiao(Ibas)=phi
                                dphidxxiao(Ibas)=dphidx
                                dphidyxiao(Ibas)=dphidy
                                dphidzxiao(Ibas)=dphidz
                            enddo
                            
                    ! Next, evaluate the densities at the grid point and the gradient
                    ! at that grid point.

                        call denspt(gridx,gridy,gridz,density,densityb,gax,gay,gaz, &
                        gbx,gby,gbz)

                        if (density < quick_method%DMCutoff ) then
                            continue
                        else

                    call becke_E(density,densityb,gax,gay,gaz,gbx,gby,gbz,Ex)

                    call lyp_e(density,densityb,gax,gay,gaz,gbx,gby,gbz,Ec)

                    Eelxc = Eelxc + (param7*Ex+param8*Ec) &
                            *weight


                            aelec = weight*density+aelec
                            belec = weight*densityb+belec

                        ! This allows the calculation of the derivative of the functional
                        ! with regard to the density (dfdr), with regard to the alpha-alpha
                        ! density invariant (df/dgaa), and the alpha-beta density invariant.

                            call becke(density,gax,gay,gaz,gbx,gby,gbz, &
                            dfdr,dfdgaa,dfdgab)

                            call lyp(density,densityb,gax,gay,gaz,gbx,gby,gbz, &
                            dfdr2,dfdgaa2,dfdgab2)

!                            dfdr = param7*dfdr+param8*dfdr2
!                            dfdgaa = param7*dfdgaa + param8*dfdgaa2
!                            dfdgab = param7*dfdgab + param8*dfdgab2

                            dfdr = dfdr+dfdr2
                            dfdgaa = dfdgaa + dfdgaa2
                            dfdgab = dfdgab + dfdgab2

                        ! Calculate the first term in the dot product shown above,i.e.:
                        ! (2 df/dgaa Grad(rho a) + df/dgab Grad(rho b)) doT Grad(Phimu Phinu))

                            xdot=2.d0*dfdgaa*gax+dfdgab*gbx
                            ydot=2.d0*dfdgaa*gay+dfdgab*gby
                            zdot=2.d0*dfdgaa*gaz+dfdgab*gbz

                        ! Now loop over basis functions and compute the addition to the matrix
                        ! element.

                            do Ibas=1,nbasis
!                                call pteval(gridx,gridy,gridz,phi,dphidx,dphidy, &
!                                dphidz,Ibas)
        phi=phixiao(Ibas)
        dphidx=dphidxxiao(Ibas)
        dphidy=dphidyxiao(Ibas)
        dphidz=dphidzxiao(Ibas)
!                                quicktest = DABS(dphidx)+DABS(dphidy)+DABS(dphidz) &
!                                +DABS(phi)
                                quicktest = DABS(dphidx+dphidy+dphidz+ &
                                phi)
                                if (quicktest < quick_method%DMCutoff ) then
                                    continue
                                else
                                    do Jbas=Ibas,nbasis
!                                        call pteval(gridx,gridy,gridz,phi2,dphi2dx,dphi2dy, &
!                                        dphi2dz,Jbas)
        phi2=phixiao(Jbas)
        dphi2dx=dphidxxiao(Jbas)
        dphi2dy=dphidyxiao(Jbas)
        dphi2dz=dphidzxiao(Jbas)
!                                        quicktest = DABS(dphi2dx)+DABS(dphi2dy)+DABS(dphi2dz) &
!                                        +DABS(phi2)
!                                        if (quicktest < quick_method%DMCutoff ) then
!                                            continue
!                                        else
                                            temp = phi*phi2
                                            tempgx = phi*dphi2dx + phi2*dphidx
                                            tempgy = phi*dphi2dy + phi2*dphidy
                                            tempgz = phi*dphi2dz + phi2*dphidz
                                            O(Jbas,Ibas)=O(Jbas,Ibas)+(temp*dfdr+ &
                                            xdot*tempgx+ydot*tempgy+zdot*tempgz)*weight
!                                        endif
                                    enddo
                                endif
                            enddo
                        endif
                    endif
                enddo
            enddo
        enddo
     endif

    if(quick_method%MPW91PW91)then
         do Iatm=1,natom
           If(quick_method%ISG.eq.1)then
             Iradtemp=50
           else
             If(quick_molspec%iattype(iatm).le.10)then
               Iradtemp=23
             else
               Iradtemp=26
             endif
           endif

           do Irad=1,Iradtemp
             If(quick_method%ISG.eq.1)then
               call gridformnew(iatm,RGRID(Irad),iiangt)
               rad = radii(quick_molspec%iattype(iatm))
             else
               call gridformSG0(iatm,Iradtemp+1-Irad,iiangt,RGRID,RWT)
               rad = radii2(quick_molspec%iattype(iatm))
             endif

             rad3 = rad*rad*rad
!             print*,iiangt
               do Iang=1,iiangt
                    gridx=xyz(1,Iatm)+rad*RGRID(Irad)*XANG(Iang)
                    gridy=xyz(2,Iatm)+rad*RGRID(Irad)*YANG(Iang)
                    gridz=xyz(3,Iatm)+rad*RGRID(Irad)*ZANG(Iang)

                ! Next, calculate the weight of the grid point in the SSW scheme.  If
                ! the grid point has a zero weight, we can skip it.

                    weight=SSW(gridx,gridy,gridz,Iatm) &
                    *WTANG(Iang)*RWT(Irad)*rad3

                    if (weight < quick_method%DMCutoff ) then
                        continue
                    else

                            do Ibas=1,nbasis
                                call pteval(gridx,gridy,gridz,phi,dphidx,dphidy, &
                                dphidz,Ibas)
                                phixiao(Ibas)=phi
                                dphidxxiao(Ibas)=dphidx
                                dphidyxiao(Ibas)=dphidy
                                dphidzxiao(Ibas)=dphidz
                            enddo
                            
                    ! Next, evaluate the densities at the grid point and the gradient
                    ! at that grid point.

                        call denspt(gridx,gridy,gridz,density,densityb,gax,gay,gaz, &
                        gbx,gby,gbz)

                        if (density < quick_method%DMCutoff ) then
                            continue
                        else

                    call becke_E(density,densityb,gax,gay,gaz,gbx,gby,gbz,Ex)

                    call lyp_e(density,densityb,gax,gay,gaz,gbx,gby,gbz,Ec)

                    Eelxc = Eelxc + (param7*Ex+param8*Ec) &
                            *weight


                            aelec = weight*density+aelec
                            belec = weight*densityb+belec

                        ! This allows the calculation of the derivative of the functional
                        ! with regard to the density (dfdr), with regard to the alpha-alpha
                        ! density invariant (df/dgaa), and the alpha-beta density invariant.

                            call becke(density,gax,gay,gaz,gbx,gby,gbz, &
                            dfdr,dfdgaa,dfdgab)

                            call lyp(density,densityb,gax,gay,gaz,gbx,gby,gbz, &
                            dfdr2,dfdgaa2,dfdgab2)

!                            dfdr = param7*dfdr+param8*dfdr2
!                            dfdgaa = param7*dfdgaa + param8*dfdgaa2
!                            dfdgab = param7*dfdgab + param8*dfdgab2

                            dfdr = dfdr+dfdr2
                            dfdgaa = dfdgaa + dfdgaa2
                            dfdgab = dfdgab + dfdgab2

                        ! Calculate the first term in the dot product shown above,i.e.:
                        ! (2 df/dgaa Grad(rho a) + df/dgab Grad(rho b)) doT Grad(Phimu Phinu))

                            xdot=2.d0*dfdgaa*gax+dfdgab*gbx
                            ydot=2.d0*dfdgaa*gay+dfdgab*gby
                            zdot=2.d0*dfdgaa*gaz+dfdgab*gbz

                        ! Now loop over basis functions and compute the addition to the matrix
                        ! element.

                            do Ibas=1,nbasis
!                                call pteval(gridx,gridy,gridz,phi,dphidx,dphidy, &
!                                dphidz,Ibas)
        phi=phixiao(Ibas)
        dphidx=dphidxxiao(Ibas)
        dphidy=dphidyxiao(Ibas)
        dphidz=dphidzxiao(Ibas)
!                                quicktest = DABS(dphidx)+DABS(dphidy)+DABS(dphidz) &
!                                +DABS(phi)
                                quicktest = DABS(dphidx+dphidy+dphidz+ &
                                phi)
                                if (quicktest < quick_method%DMCutoff ) then
                                    continue
                                else
                                    do Jbas=Ibas,nbasis
!                                        call pteval(gridx,gridy,gridz,phi2,dphi2dx,dphi2dy, &
!                                        dphi2dz,Jbas)
        phi2=phixiao(Jbas)
        dphi2dx=dphidxxiao(Jbas)
        dphi2dy=dphidyxiao(Jbas)
        dphi2dz=dphidzxiao(Jbas)
!                                        quicktest = DABS(dphi2dx)+DABS(dphi2dy)+DABS(dphi2dz) &
!                                        +DABS(phi2)
!                                        if (quicktest < quick_method%DMCutoff ) then
!                                            continue
!                                        else
                                            temp = phi*phi2
                                            tempgx = phi*dphi2dx + phi2*dphidx
                                            tempgy = phi*dphi2dy + phi2*dphidy
                                            tempgz = phi*dphi2dz + phi2*dphidz
                                            O(Jbas,Ibas)=O(Jbas,Ibas)+(temp*dfdr+ &
                                            xdot*tempgx+ydot*tempgy+zdot*tempgz)*weight
!                                        endif
                                    enddo
                                endif
                            enddo
                        endif
                    endif
                enddo
            enddo
        enddo
     endif



    do Ibas=1,nbasis
        do Jbas=Ibas+1,nbasis
            O(Ibas,Jbas) = O(Jbas,Ibas)
        enddo
    enddo

        call cpu_time(t3)

        write (ioutfile,'(" TIME of evaluation numerical integral = ",F12.2)') &
        T3-T2

!   if(printEnergy)then
!    do Ibas=1,nbasis
!        do Jbas=1,nbasis
!            Eel=Eel+DENSE(Ibas,Jbas)*O(Jbas,Ibas)
!        enddo
!    enddo

    Eel=Eel+Eelxc
!   endif



    END subroutine dftoperator

! Ed Brothers. November 27, 2001
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

    subroutine dftoperatordelta
    use allmod
    implicit double precision(a-h,o-z)

! The purpose of this subroutine is to form the operator matrix
! for a full Density Functional calculation, i.e. the KS matrix.  The
! KS  matrix is as follows:

! O(I,J) =  F(I,J) = KE(I,J) + IJ attraction to each atom + repulsion_prim
! with each possible basis  + Exchange/Correlation functional.

! Note that the KS operator matrix is symmetric.

  call cpu_time(t1)

  Eelxc=0.0d0

  if(quick_method%printEnergy)then
    Eel=0.d0
    do Ibas=1,nbasis
        do Icon=1,ncontract(Ibas)
            do Jcon=1,ncontract(Ibas)

            ! Kinetic energy.

                Eel=Eel+DENSESAVE(Ibas,Ibas)* &
                dcoeff(Jcon,Ibas)*dcoeff(Icon,Ibas)* &
                ekinetic(aexp(Jcon,Ibas),aexp(Icon,Ibas), &
                itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                xyz(1,ncenter(Ibas)),xyz(2,ncenter(Ibas)), &
                xyz(3,ncenter(Ibas)),xyz(1,ncenter(Ibas)), &
                xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))

            ! Nuclear attraction.

!                do iatom = 1,natom
!                    Eel=Eel+DENSESAVE(Ibas,Ibas)* &
!                    dcoeff(Jcon,Ibas)*dcoeff(Icon,Ibas)* &
!                    attraction(aexp(Jcon,Ibas),aexp(Icon,Ibas), &
!                    itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
!                    itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
!                    xyz(1,ncenter(Ibas)),xyz(2,ncenter(Ibas)), &
!                    xyz(3,ncenter(Ibas)),xyz(1,ncenter(Ibas)), &
!                    xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)), &
!                    xyz(1,iatom),xyz(2,iatom),xyz(3,iatom), &
!                    chg(iatom))
!                enddo
            enddo
        enddo
    enddo

    do Ibas=1,nbasis
        do Jbas=Ibas+1,nbasis
            do Icon=1,ncontract(ibas)
                do Jcon=1,ncontract(jbas)

                ! Kinetic energy.

                    Eel=Eel+DENSESAVE(Jbas,Ibas)* &
                    dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas)* &
                    2.d0*ekinetic(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                    itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                    itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                    xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                    xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                    xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))

                ! Nuclear attraction.

!                    do iatom = 1,natom
!                        Eel=Eel+DENSESAVE(Jbas,Ibas)* &
!                        dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas)* &
!                        2.d0*attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
!                        itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
!                        itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
!                        xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
!                        xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
!                        xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)), &
!                        xyz(1,iatom),xyz(2,iatom),xyz(3,iatom), &
!                        chg(iatom))
!                    enddo
                enddo
            enddo
        enddo
    enddo

   do itemp1=1,nbasis
    do jtemp2=1,nbasis
      HOLD(jtemp2,itemp1)=DENSE(jtemp2,itemp1)
      DENSE(jtemp2,itemp1)=DENSEsave(jtemp2,itemp1)
    enddo
   enddo

    do IIsh=1,jshell
      do JJsh=IIsh,jshell
        call attrashellenergy(IIsh,JJsh)
      enddo
    enddo

   do itemp1=1,nbasis
    do jtemp2=1,nbasis
      DENSE(jtemp2,itemp1)=HOLD(jtemp2,itemp1)
    enddo
   enddo

  endif

  call cpu_time(t2)

!           write (ioutfile,'("TOTAL ENERGY OF CURRENT CYCLE=",F16.9)') Eel

         if(quick_method%printEnergy)then
           write (ioutfile,'("Time for one-electron energy evaluation=",F16.9)') t2-t1
         endif

!    print*,'Eel=', Eel


! The first part is the one elctron code.

        do I=1,nbasis
          do J=1,nbasis
            O(i,j)=Osavedft(i,j)
          enddo
        enddo

!
! Alessandro GENONI 03/21/2007
! Sum the ECP integrals to the partial Fock matrix
!
    if (quick_method%ecp) then
      call ecpoperator
    end if

! The previous two terms are the one electron part of the Fock matrix.
! The next term defines the electron repulsion_prim.

 call cpu_time(T1)

! Delta density matrix cutoff

 do II=1,jshell
   do JJ=II,jshell
     DNtemp=0.0d0
     call DNscreen(II,JJ,DNtemp)
     Cutmatrix(II,JJ)=DNtemp
     Cutmatrix(JJ,II)=DNtemp
   enddo
 enddo

! Schwartz cutoff is implemented here. (ab|cd)**2<=(ab|ab)*(cd|cd)
! Reference: Strout DL and Scuseria JCP 102(1995),8448.

! print*,"before 2e"
if(quick_method%B3LYP)then
 do II=1,jshell
   do JJ=II,jshell
     Testtmp=Ycutoff(II,JJ)
     do KK=II,jshell
       do LL=KK,jshell
!          Nxiao1=Nxiao1+1
            cutoffTest1 = TESTtmp*Ycutoff(KK,LL)
          If(cutoffTest1.gt.quick_method%integralCutoff)then
            DNmax=max(4.0d0*cutmatrix(II,JJ),4.0d0*cutmatrix(KK,LL), &
                  cutmatrix(II,LL),cutmatrix(II,KK),cutmatrix(JJ,KK),cutmatrix(JJ,LL))
!            DNmax=max(cutmatrix(II,JJ),cutmatrix(KK,LL) &
!                  )
            cutoffTest=testCutoff*DNmax
            If(cutoffTest.gt.quick_method%integralCutoff)then
              IIxiao=II
              JJxiao=JJ
              KKxiao=KK
              LLxiao=LL
              call shelldftb3lyp(IIxiao,JJxiao,KKxiao,LLxiao)
!            Nxiao2=Nxiao2+1
            endif
!            else
!             print*,II,JJ,KK,LL,cutoffTest,testCutoff,DNmax
!            print*,'***',O(1,1)
          endif
       enddo
     enddo
   enddo
 enddo
else
 do II=1,jshell
   do JJ=II,jshell
     Testtmp=Ycutoff(II,JJ)
     do KK=II,jshell
       do LL=KK,jshell
!          Nxiao1=Nxiao1+1
            cutoffTest1 = TESTtmp*Ycutoff(KK,LL)
          If(cutoffTest1.gt.quick_method%integralCutoff)then
!            DNmax=max(4.0d0*cutmatrix(II,JJ),4.0d0*cutmatrix(KK,LL), &
!                  cutmatrix(II,LL),cutmatrix(II,KK),cutmatrix(JJ,KK),cutmatrix(JJ,LL))
            DNmax=max(cutmatrix(II,JJ),cutmatrix(KK,LL) &
                  )
            cutoffTest=testCutoff*DNmax
            If(cutoffTest.gt.quick_method%integralCutoff)then
              IIxiao=II
              JJxiao=JJ
              KKxiao=KK
              LLxiao=LL
              call shelldft(IIxiao,JJxiao,KKxiao,LLxiao)
!            Nxiao2=Nxiao2+1
            endif
!            else
!             print*,II,JJ,KK,LL,cutoffTest,testCutoff,DNmax
!            print*,'***',O(1,1)
          endif
       enddo
     enddo
   enddo
 enddo
endif


        do I=1,nbasis
          do J=1,nbasis
            Osavedft(i,j)=O(i,j)
          enddo
        enddo

            do I=1,nbasis
              do J=1,nbasis
                DENSE(I,J)=DENSESAVE(I,J)
              enddo
           enddo

        call cpu_time(t2)

        write (ioutfile,'(" TIME of evaluation integral = ",F12.2)') &
        T2-T1

! print*,'Nxiao1=',Nxiao1,'Nxiao2=',Nxiao2,integralCutOff

    do Ibas=1,nbasis
        do Jbas=Ibas+1,nbasis
            O(Ibas,Jbas) = O(Jbas,Ibas)
        enddo
    enddo

  if(quick_method%printEnergy)then
    do Ibas=1,nbasis
        do Jbas=1,nbasis
            Eel=Eel+DENSESAVE(Ibas,Jbas)*O(Jbas,Ibas)
        enddo
    enddo

    Eel=Eel/2.0d0
 endif

!           write (ioutfile,'("TOTAL ENERGY OF CURRENT CYCLE=",F16.9)') Eel

! stop

! The next portion is the exchange/correlation functional.
! The angular grid code came from CCL.net.  The radial grid
! formulas (position and wieghts) is from Gill, Johnson and Pople,
! Chem. Phys. Lett. v209, n 5+6, 1993, pg 506-512.  The weighting scheme
! is from Stratmann, Scuseria, and Frisch, Chem. Phys. Lett., v 257,
! 1996, pg 213-223.

! The actual element is:
! F alpha mu nu = Integral((df/drhoa Phimu Phinu)+
! (2 df/dgaa Grad(rho a) + df/dgab Grad(rho b)) doT Grad(Phimu Phinu))

! where F alpha mu nu is the the alpha spin portion of the operator matrix
! element mu, nu,
! df/drhoa is the derivative of the functional by the alpha density,
! df/dgaa is the derivative of the functional by the alpha gradient
! invariant, i.e. the dot product of the gradient of the alpha
! density with itself.
! df/dgab is the derivative of the functional by the dot product of
! the gradient of the alpha density with the beta density.
! Grad(Phimu Phinu) is the gradient of Phimu times Phinu.

! First, find the grid point.

    aelec=0.d0
    belec=0.d0

! Xiao HE 02/09/1007 SG0 SG1 grids
!    do Ireg=1,iregion
!        call gridform(iangular(Ireg))
!        do Irad=iradial(ireg-1)+1,iradial(ireg)
!            do Iatm=1,natom
!                rad = radii(quick_molspec%iattype(iatm))
!                rad3 = rad*rad*rad
!                do Iang=1,iangular(Ireg)

if(quick_method%B3LYP)then
         do Iatm=1,natom
           If(quick_method%ISG.eq.1)then
             Iradtemp=50
           else
             If(quick_molspec%iattype(iatm).le.10)then
               Iradtemp=23
             else
               Iradtemp=26
             endif
           endif

           do Irad=1,Iradtemp
             If(quick_method%ISG.eq.1)then
               call gridformnew(iatm,RGRID(Irad),iiangt)
               rad = radii(quick_molspec%iattype(iatm))
             else
               call gridformSG0(iatm,Iradtemp+1-Irad,iiangt,RGRID,RWT)
               rad = radii2(quick_molspec%iattype(iatm))
             endif

             rad3 = rad*rad*rad
               do Iang=1,iiangt
                    gridx=xyz(1,Iatm)+rad*RGRID(Irad)*XANG(Iang)
                    gridy=xyz(2,Iatm)+rad*RGRID(Irad)*YANG(Iang)
                    gridz=xyz(3,Iatm)+rad*RGRID(Irad)*ZANG(Iang)

                ! Next, calculate the weight of the grid point in the SSW scheme.  If
                ! the grid point has a zero weight, we can skip it.

                    weight=SSW(gridx,gridy,gridz,Iatm) &
                    *WTANG(Iang)*RWT(Irad)*rad3

                    if (weight < quick_method%DMCutoff ) then
                        continue
                    else

                            do Ibas=1,nbasis
                                call pteval(gridx,gridy,gridz,phi,dphidx,dphidy, &
                                dphidz,Ibas)
                                phixiao(Ibas)=phi
                                dphidxxiao(Ibas)=dphidx
                                dphidyxiao(Ibas)=dphidy
                                dphidzxiao(Ibas)=dphidz
                            enddo
                            
                    ! Next, evaluate the densities at the grid point and the gradient
                    ! at that grid point.

                        call denspt(gridx,gridy,gridz,density,densityb,gax,gay,gaz, &
                        gbx,gby,gbz)

                        if (density < quick_method%DMCutoff ) then
                            continue
                        else

                     densitysum=2.0d0*density
                     sigma=4.0d0*(gax*gax+gay*gay+gaz*gaz)
                     call b3lyp_e(densitysum,sigma,zkec)
!
!                    Eel = Eel + (param7*Ex+param8*Ec) &
                     Eelxc = Eelxc + zkec*weight

                            aelec = weight*density+aelec
                            belec = weight*densityb+belec

                        ! This allows the calculation of the derivative of the functional
                        ! with regard to the density (dfdr), with regard to the alpha-alpha
                        ! density invariant (df/dgaa), and the alpha-beta density invariant.

!                            call becke(density,gax,gay,gaz,gbx,gby,gbz, &
!                            dfdr,dfdgaa,dfdgab)
!
!                            call lyp(density,densityb,gax,gay,gaz,gbx,gby,gbz, &
!                            dfdr2,dfdgaa2,dfdgab2)

                            densitysum=2.0d0*density
                            sigma=4.0d0*(gax*gax+gay*gay+gaz*gaz)

                            call b3lypf(densitysum,sigma,dfdr,xiaodot)

!                            dfdr=dfdr*2.0d0

!                            dfdr = dfdr+dfdr2
!                            dfdgaa = dfdgaa + dfdgaa2
!                            dfdgab = dfdgab + dfdgab2

                        ! Calculate the first term in the dot product shown above,i.e.:
                        ! (2 df/dgaa Grad(rho a) + df/dgab Grad(rho b)) doT Grad(Phimu Phinu))

                             xdot=xiaodot*gax
                             ydot=xiaodot*gay
                             zdot=xiaodot*gaz

!                            xdot=2.d0*dfdgaa*gax+dfdgab*gbx
!                            ydot=2.d0*dfdgaa*gay+dfdgab*gby
!                            zdot=2.d0*dfdgaa*gaz+dfdgab*gbz

                        ! Now loop over basis functions and compute the addition to the matrix
                        ! element.

                            do Ibas=1,nbasis
!                                call pteval(gridx,gridy,gridz,phi,dphidx,dphidy, &
!                                dphidz,Ibas)
        phi=phixiao(Ibas)
        dphidx=dphidxxiao(Ibas)
        dphidy=dphidyxiao(Ibas)
        dphidz=dphidzxiao(Ibas)
!                                quicktest = DABS(dphidx)+DABS(dphidy)+DABS(dphidz) &
!                                +DABS(phi)
                                quicktest = DABS(dphidx+dphidy+dphidz &
                                +phi)
                                if (quicktest < quick_method%DMCutoff ) then
                                    continue
                                else
                                    do Jbas=Ibas,nbasis
!                                        call pteval(gridx,gridy,gridz,phi2,dphi2dx,dphi2dy, &
!                                        dphi2dz,Jbas)
        phi2=phixiao(Jbas)
        dphi2dx=dphidxxiao(Jbas)
        dphi2dy=dphidyxiao(Jbas)
        dphi2dz=dphidzxiao(Jbas)
!                                        quicktest = DABS(dphi2dx)+DABS(dphi2dy)+DABS(dphi2dz) &
!                                        +DABS(phi2)
!                                        if (quicktest < quick_method%DMCutoff ) then
!                                            continue
!                                        else
                                            temp = phi*phi2
                                            tempgx = phi*dphi2dx + phi2*dphidx
                                            tempgy = phi*dphi2dy + phi2*dphidy
                                            tempgz = phi*dphi2dz + phi2*dphidz
                                            O(Jbas,Ibas)=O(Jbas,Ibas)+(temp*dfdr+ &
                                            xdot*tempgx+ydot*tempgy+zdot*tempgz)*weight
!                                        endif
                                    enddo
                                endif
                            enddo
                        endif
                    endif
                enddo
            enddo
        enddo
!    enddo
endif

    if(quick_method%BLYP)then
         do Iatm=1,natom
           If(quick_method%ISG.eq.1)then
             Iradtemp=50
           else
             If(quick_molspec%iattype(iatm).le.10)then
               Iradtemp=23
             else
               Iradtemp=26
             endif
           endif

           do Irad=1,Iradtemp
             If(quick_method%ISG.eq.1)then
               call gridformnew(iatm,RGRID(Irad),iiangt)
               rad = radii(quick_molspec%iattype(iatm))
             else
               call gridformSG0(iatm,Iradtemp+1-Irad,iiangt,RGRID,RWT)
               rad = radii2(quick_molspec%iattype(iatm))
             endif

             rad3 = rad*rad*rad
!             print*,iiangt
               do Iang=1,iiangt
                    gridx=xyz(1,Iatm)+rad*RGRID(Irad)*XANG(Iang)
                    gridy=xyz(2,Iatm)+rad*RGRID(Irad)*YANG(Iang)
                    gridz=xyz(3,Iatm)+rad*RGRID(Irad)*ZANG(Iang)

                ! Next, calculate the weight of the grid point in the SSW scheme.  If
                ! the grid point has a zero weight, we can skip it.

                    weight=SSW(gridx,gridy,gridz,Iatm) &
                    *WTANG(Iang)*RWT(Irad)*rad3

                    if (weight < quick_method%DMCutoff ) then
                        continue
                    else

                            do Ibas=1,nbasis
                                call pteval(gridx,gridy,gridz,phi,dphidx,dphidy, &
                                dphidz,Ibas)
                                phixiao(Ibas)=phi
                                dphidxxiao(Ibas)=dphidx
                                dphidyxiao(Ibas)=dphidy
                                dphidzxiao(Ibas)=dphidz
                            enddo
                            
                    ! Next, evaluate the densities at the grid point and the gradient
                    ! at that grid point.

                        call denspt(gridx,gridy,gridz,density,densityb,gax,gay,gaz, &
                        gbx,gby,gbz)

                        if (density < quick_method%DMCutoff ) then
                            continue
                        else

                    call becke_E(density,densityb,gax,gay,gaz,gbx,gby,gbz,Ex)

                    call lyp_e(density,densityb,gax,gay,gaz,gbx,gby,gbz,Ec)

                    Eelxc = Eelxc + (param7*Ex+param8*Ec) &
                    *weight


                            aelec = weight*density+aelec
                            belec = weight*densityb+belec

                        ! This allows the calculation of the derivative of the functional
                        ! with regard to the density (dfdr), with regard to the alpha-alpha
                        ! density invariant (df/dgaa), and the alpha-beta density invariant.

                            call becke(density,gax,gay,gaz,gbx,gby,gbz, &
                            dfdr,dfdgaa,dfdgab)

                            call lyp(density,densityb,gax,gay,gaz,gbx,gby,gbz, &
                            dfdr2,dfdgaa2,dfdgab2)

!                            dfdr = param7*dfdr+param8*dfdr2
!                            dfdgaa = param7*dfdgaa + param8*dfdgaa2
!                            dfdgab = param7*dfdgab + param8*dfdgab2

                            dfdr = dfdr+dfdr2
                            dfdgaa = dfdgaa + dfdgaa2
                            dfdgab = dfdgab + dfdgab2

                        ! Calculate the first term in the dot product shown above,i.e.:
                        ! (2 df/dgaa Grad(rho a) + df/dgab Grad(rho b)) doT Grad(Phimu Phinu))

                            xdot=2.d0*dfdgaa*gax+dfdgab*gbx
                            ydot=2.d0*dfdgaa*gay+dfdgab*gby
                            zdot=2.d0*dfdgaa*gaz+dfdgab*gbz

                        ! Now loop over basis functions and compute the addition to the matrix
                        ! element.

                            do Ibas=1,nbasis
!                                call pteval(gridx,gridy,gridz,phi,dphidx,dphidy, &
!                                dphidz,Ibas)
        phi=phixiao(Ibas)
        dphidx=dphidxxiao(Ibas)
        dphidy=dphidyxiao(Ibas)
        dphidz=dphidzxiao(Ibas)
!                                quicktest = DABS(dphidx)+DABS(dphidy)+DABS(dphidz) &
!                                +DABS(phi)
                                quicktest = DABS(dphidx+dphidy+dphidz+ &
                                phi)
                                if (quicktest < quick_method%DMCutoff ) then
                                    continue
                                else
                                    do Jbas=Ibas,nbasis
!                                        call pteval(gridx,gridy,gridz,phi2,dphi2dx,dphi2dy, &
!                                        dphi2dz,Jbas)
        phi2=phixiao(Jbas)
        dphi2dx=dphidxxiao(Jbas)
        dphi2dy=dphidyxiao(Jbas)
        dphi2dz=dphidzxiao(Jbas)
!                                        quicktest = DABS(dphi2dx)+DABS(dphi2dy)+DABS(dphi2dz) &
!                                        +DABS(phi2)
!                                        if (quicktest < quick_method%DMCutoff ) then
!                                            continue
!                                        else
                                            temp = phi*phi2
                                            tempgx = phi*dphi2dx + phi2*dphidx
                                            tempgy = phi*dphi2dy + phi2*dphidy
                                            tempgz = phi*dphi2dz + phi2*dphidz
                                            O(Jbas,Ibas)=O(Jbas,Ibas)+(temp*dfdr+ &
                                            xdot*tempgx+ydot*tempgy+zdot*tempgz)*weight
!                                        endif
                                    enddo
                                endif
                            enddo
                        endif
                    endif
                enddo
            enddo
        enddo
     endif

    if(quick_method%MPW91LYP)then
         do Iatm=1,natom
           If(quick_method%ISG.eq.1)then
             Iradtemp=50
           else
             If(quick_molspec%iattype(iatm).le.10)then
               Iradtemp=23
             else
               Iradtemp=26
             endif
           endif

           do Irad=1,Iradtemp
             If(quick_method%ISG.eq.1)then
               call gridformnew(iatm,RGRID(Irad),iiangt)
               rad = radii(quick_molspec%iattype(iatm))
             else
               call gridformSG0(iatm,Iradtemp+1-Irad,iiangt,RGRID,RWT)
               rad = radii2(quick_molspec%iattype(iatm))
             endif

             rad3 = rad*rad*rad
!             print*,iiangt
               do Iang=1,iiangt
                    gridx=xyz(1,Iatm)+rad*RGRID(Irad)*XANG(Iang)
                    gridy=xyz(2,Iatm)+rad*RGRID(Irad)*YANG(Iang)
                    gridz=xyz(3,Iatm)+rad*RGRID(Irad)*ZANG(Iang)

                ! Next, calculate the weight of the grid point in the SSW scheme.  If
                ! the grid point has a zero weight, we can skip it.

                    weight=SSW(gridx,gridy,gridz,Iatm) &
                    *WTANG(Iang)*RWT(Irad)*rad3

                    if (weight < quick_method%DMCutoff ) then
                        continue
                    else

                            do Ibas=1,nbasis
                                call pteval(gridx,gridy,gridz,phi,dphidx,dphidy, &
                                dphidz,Ibas)
                                phixiao(Ibas)=phi
                                dphidxxiao(Ibas)=dphidx
                                dphidyxiao(Ibas)=dphidy
                                dphidzxiao(Ibas)=dphidz
                            enddo
                            
                    ! Next, evaluate the densities at the grid point and the gradient
                    ! at that grid point.

                        call denspt(gridx,gridy,gridz,density,densityb,gax,gay,gaz, &
                        gbx,gby,gbz)

                        if (density < quick_method%DMCutoff ) then
                            continue
                        else

                    call mpw91_E(density,densityb,gax,gay,gaz,gbx,gby,gbz,Ex)

                    call lyp_e(density,densityb,gax,gay,gaz,gbx,gby,gbz,Ec)

                    Eelxc = Eelxc + (param7*Ex+param8*Ec) &
                    *weight


                            aelec = weight*density+aelec
                            belec = weight*densityb+belec

                        ! This allows the calculation of the derivative of the functional
                        ! with regard to the density (dfdr), with regard to the alpha-alpha
                        ! density invariant (df/dgaa), and the alpha-beta density invariant.

                            call mpw91(density,gax,gay,gaz,gbx,gby,gbz, &
                            dfdr,dfdgaa,dfdgab)

                            call lyp(density,densityb,gax,gay,gaz,gbx,gby,gbz, &
                            dfdr2,dfdgaa2,dfdgab2)

!                            dfdr = param7*dfdr+param8*dfdr2
!                            dfdgaa = param7*dfdgaa + param8*dfdgaa2
!                            dfdgab = param7*dfdgab + param8*dfdgab2

                            dfdr = dfdr+dfdr2
                            dfdgaa = dfdgaa + dfdgaa2
                            dfdgab = dfdgab + dfdgab2

                        ! Calculate the first term in the dot product shown above,i.e.:
                        ! (2 df/dgaa Grad(rho a) + df/dgab Grad(rho b)) doT Grad(Phimu Phinu))

                            xdot=2.d0*dfdgaa*gax+dfdgab*gbx
                            ydot=2.d0*dfdgaa*gay+dfdgab*gby
                            zdot=2.d0*dfdgaa*gaz+dfdgab*gbz

                        ! Now loop over basis functions and compute the addition to the matrix
                        ! element.

                            do Ibas=1,nbasis
!                                call pteval(gridx,gridy,gridz,phi,dphidx,dphidy, &
!                                dphidz,Ibas)
        phi=phixiao(Ibas)
        dphidx=dphidxxiao(Ibas)
        dphidy=dphidyxiao(Ibas)
        dphidz=dphidzxiao(Ibas)
!                                quicktest = DABS(dphidx)+DABS(dphidy)+DABS(dphidz) &
!                                +DABS(phi)
                                quicktest = DABS(dphidx+dphidy+dphidz+ &
                                phi)
                                if (quicktest < quick_method%DMCutoff ) then
                                    continue
                                else
                                    do Jbas=Ibas,nbasis
!                                        call pteval(gridx,gridy,gridz,phi2,dphi2dx,dphi2dy, &
!                                        dphi2dz,Jbas)
        phi2=phixiao(Jbas)
        dphi2dx=dphidxxiao(Jbas)
        dphi2dy=dphidyxiao(Jbas)
        dphi2dz=dphidzxiao(Jbas)
!                                        quicktest = DABS(dphi2dx)+DABS(dphi2dy)+DABS(dphi2dz) &
!                                        +DABS(phi2)
!                                        if (quicktest < quick_method%DMCutoff ) then
!                                            continue
!                                        else
                                            temp = phi*phi2
                                            tempgx = phi*dphi2dx + phi2*dphidx
                                            tempgy = phi*dphi2dy + phi2*dphidy
                                            tempgz = phi*dphi2dz + phi2*dphidz
                                            O(Jbas,Ibas)=O(Jbas,Ibas)+(temp*dfdr+ &
                                            xdot*tempgx+ydot*tempgy+zdot*tempgz)*weight
!                                        endif
                                    enddo
                                endif
                            enddo
                        endif
                    endif
                enddo
            enddo
        enddo
     endif

    if(quick_method%BPW91)then
         do Iatm=1,natom
           If(quick_method%ISG.eq.1)then
             Iradtemp=50
           else
             If(quick_molspec%iattype(iatm).le.10)then
               Iradtemp=23
             else
               Iradtemp=26
             endif
           endif

           do Irad=1,Iradtemp
             If(quick_method%ISG.eq.1)then
               call gridformnew(iatm,RGRID(Irad),iiangt)
               rad = radii(quick_molspec%iattype(iatm))
             else
               call gridformSG0(iatm,Iradtemp+1-Irad,iiangt,RGRID,RWT)
               rad = radii2(quick_molspec%iattype(iatm))
             endif

             rad3 = rad*rad*rad
!             print*,iiangt
               do Iang=1,iiangt
                    gridx=xyz(1,Iatm)+rad*RGRID(Irad)*XANG(Iang)
                    gridy=xyz(2,Iatm)+rad*RGRID(Irad)*YANG(Iang)
                    gridz=xyz(3,Iatm)+rad*RGRID(Irad)*ZANG(Iang)

                ! Next, calculate the weight of the grid point in the SSW scheme.  If
                ! the grid point has a zero weight, we can skip it.

                    weight=SSW(gridx,gridy,gridz,Iatm) &
                    *WTANG(Iang)*RWT(Irad)*rad3

                    if (weight < quick_method%DMCutoff ) then
                        continue
                    else

                            do Ibas=1,nbasis
                                call pteval(gridx,gridy,gridz,phi,dphidx,dphidy, &
                                dphidz,Ibas)
                                phixiao(Ibas)=phi
                                dphidxxiao(Ibas)=dphidx
                                dphidyxiao(Ibas)=dphidy
                                dphidzxiao(Ibas)=dphidz
                            enddo
                            
                    ! Next, evaluate the densities at the grid point and the gradient
                    ! at that grid point.

                        call denspt(gridx,gridy,gridz,density,densityb,gax,gay,gaz, &
                        gbx,gby,gbz)

                        if (density < quick_method%DMCutoff ) then
                            continue
                        else

                    call becke_E(density,densityb,gax,gay,gaz,gbx,gby,gbz,Ex)

                    call lyp_e(density,densityb,gax,gay,gaz,gbx,gby,gbz,Ec)

                    Eelxc = Eelxc + (param7*Ex+param8*Ec) &
                    *weight


                            aelec = weight*density+aelec
                            belec = weight*densityb+belec

                        ! This allows the calculation of the derivative of the functional
                        ! with regard to the density (dfdr), with regard to the alpha-alpha
                        ! density invariant (df/dgaa), and the alpha-beta density invariant.

                            call becke(density,gax,gay,gaz,gbx,gby,gbz, &
                            dfdr,dfdgaa,dfdgab)

                            call lyp(density,densityb,gax,gay,gaz,gbx,gby,gbz, &
                            dfdr2,dfdgaa2,dfdgab2)

!                            dfdr = param7*dfdr+param8*dfdr2
!                            dfdgaa = param7*dfdgaa + param8*dfdgaa2
!                            dfdgab = param7*dfdgab + param8*dfdgab2

                            dfdr = dfdr+dfdr2
                            dfdgaa = dfdgaa + dfdgaa2
                            dfdgab = dfdgab + dfdgab2

                        ! Calculate the first term in the dot product shown above,i.e.:
                        ! (2 df/dgaa Grad(rho a) + df/dgab Grad(rho b)) doT Grad(Phimu Phinu))

                            xdot=2.d0*dfdgaa*gax+dfdgab*gbx
                            ydot=2.d0*dfdgaa*gay+dfdgab*gby
                            zdot=2.d0*dfdgaa*gaz+dfdgab*gbz

                        ! Now loop over basis functions and compute the addition to the matrix
                        ! element.

                            do Ibas=1,nbasis
!                                call pteval(gridx,gridy,gridz,phi,dphidx,dphidy, &
!                                dphidz,Ibas)
        phi=phixiao(Ibas)
        dphidx=dphidxxiao(Ibas)
        dphidy=dphidyxiao(Ibas)
        dphidz=dphidzxiao(Ibas)
!                                quicktest = DABS(dphidx)+DABS(dphidy)+DABS(dphidz) &
!                                +DABS(phi)
                                quicktest = DABS(dphidx+dphidy+dphidz+ &
                                phi)
                                if (quicktest < quick_method%DMCutoff ) then
                                    continue
                                else
                                    do Jbas=Ibas,nbasis
!                                        call pteval(gridx,gridy,gridz,phi2,dphi2dx,dphi2dy, &
!                                        dphi2dz,Jbas)
        phi2=phixiao(Jbas)
        dphi2dx=dphidxxiao(Jbas)
        dphi2dy=dphidyxiao(Jbas)
        dphi2dz=dphidzxiao(Jbas)
!                                        quicktest = DABS(dphi2dx)+DABS(dphi2dy)+DABS(dphi2dz) &
!                                        +DABS(phi2)
!                                        if (quicktest < quick_method%DMCutoff ) then
!                                            continue
!                                        else
                                            temp = phi*phi2
                                            tempgx = phi*dphi2dx + phi2*dphidx
                                            tempgy = phi*dphi2dy + phi2*dphidy
                                            tempgz = phi*dphi2dz + phi2*dphidz
                                            O(Jbas,Ibas)=O(Jbas,Ibas)+(temp*dfdr+ &
                                            xdot*tempgx+ydot*tempgy+zdot*tempgz)*weight
!                                        endif
                                    enddo
                                endif
                            enddo
                        endif
                    endif
                enddo
            enddo
        enddo
     endif

    if(quick_method%MPW91PW91)then
         do Iatm=1,natom
           If(quick_method%ISG.eq.1)then
             Iradtemp=50
           else
             If(quick_molspec%iattype(iatm).le.10)then
               Iradtemp=23
             else
               Iradtemp=26
             endif
           endif

           do Irad=1,Iradtemp
             If(quick_method%ISG.eq.1)then
               call gridformnew(iatm,RGRID(Irad),iiangt)
               rad = radii(quick_molspec%iattype(iatm))
             else
               call gridformSG0(iatm,Iradtemp+1-Irad,iiangt,RGRID,RWT)
               rad = radii2(quick_molspec%iattype(iatm))
             endif

             rad3 = rad*rad*rad
!             print*,iiangt
               do Iang=1,iiangt
                    gridx=xyz(1,Iatm)+rad*RGRID(Irad)*XANG(Iang)
                    gridy=xyz(2,Iatm)+rad*RGRID(Irad)*YANG(Iang)
                    gridz=xyz(3,Iatm)+rad*RGRID(Irad)*ZANG(Iang)

                ! Next, calculate the weight of the grid point in the SSW scheme.  If
                ! the grid point has a zero weight, we can skip it.

                    weight=SSW(gridx,gridy,gridz,Iatm) &
                    *WTANG(Iang)*RWT(Irad)*rad3

                    if (weight < quick_method%DMCutoff ) then
                        continue
                    else

                            do Ibas=1,nbasis
                                call pteval(gridx,gridy,gridz,phi,dphidx,dphidy, &
                                dphidz,Ibas)
                                phixiao(Ibas)=phi
                                dphidxxiao(Ibas)=dphidx
                                dphidyxiao(Ibas)=dphidy
                                dphidzxiao(Ibas)=dphidz
                            enddo
                            
                    ! Next, evaluate the densities at the grid point and the gradient
                    ! at that grid point.

                        call denspt(gridx,gridy,gridz,density,densityb,gax,gay,gaz, &
                        gbx,gby,gbz)

                        if (density < quick_method%DMCutoff ) then
                            continue
                        else

                    call becke_E(density,densityb,gax,gay,gaz,gbx,gby,gbz,Ex)

                    call lyp_e(density,densityb,gax,gay,gaz,gbx,gby,gbz,Ec)

                    Eelxc = Eelxc + (param7*Ex+param8*Ec) &
                    *weight


                            aelec = weight*density+aelec
                            belec = weight*densityb+belec

                        ! This allows the calculation of the derivative of the functional
                        ! with regard to the density (dfdr), with regard to the alpha-alpha
                        ! density invariant (df/dgaa), and the alpha-beta density invariant.

                            call becke(density,gax,gay,gaz,gbx,gby,gbz, &
                            dfdr,dfdgaa,dfdgab)

                            call lyp(density,densityb,gax,gay,gaz,gbx,gby,gbz, &
                            dfdr2,dfdgaa2,dfdgab2)

!                            dfdr = param7*dfdr+param8*dfdr2
!                            dfdgaa = param7*dfdgaa + param8*dfdgaa2
!                            dfdgab = param7*dfdgab + param8*dfdgab2

                            dfdr = dfdr+dfdr2
                            dfdgaa = dfdgaa + dfdgaa2
                            dfdgab = dfdgab + dfdgab2

                        ! Calculate the first term in the dot product shown above,i.e.:
                        ! (2 df/dgaa Grad(rho a) + df/dgab Grad(rho b)) doT Grad(Phimu Phinu))

                            xdot=2.d0*dfdgaa*gax+dfdgab*gbx
                            ydot=2.d0*dfdgaa*gay+dfdgab*gby
                            zdot=2.d0*dfdgaa*gaz+dfdgab*gbz

                        ! Now loop over basis functions and compute the addition to the matrix
                        ! element.

                            do Ibas=1,nbasis
!                                call pteval(gridx,gridy,gridz,phi,dphidx,dphidy, &
!                                dphidz,Ibas)
        phi=phixiao(Ibas)
        dphidx=dphidxxiao(Ibas)
        dphidy=dphidyxiao(Ibas)
        dphidz=dphidzxiao(Ibas)
!                                quicktest = DABS(dphidx)+DABS(dphidy)+DABS(dphidz) &
!                                +DABS(phi)
                                quicktest = DABS(dphidx+dphidy+dphidz+ &
                                phi)
                                if (quicktest < quick_method%DMCutoff ) then
                                    continue
                                else
                                    do Jbas=Ibas,nbasis
!                                        call pteval(gridx,gridy,gridz,phi2,dphi2dx,dphi2dy, &
!                                        dphi2dz,Jbas)
        phi2=phixiao(Jbas)
        dphi2dx=dphidxxiao(Jbas)
        dphi2dy=dphidyxiao(Jbas)
        dphi2dz=dphidzxiao(Jbas)
!                                        quicktest = DABS(dphi2dx)+DABS(dphi2dy)+DABS(dphi2dz) &
!                                        +DABS(phi2)
!                                        if (quicktest < quick_method%DMCutoff ) then
!                                            continue
!                                        else
                                            temp = phi*phi2
                                            tempgx = phi*dphi2dx + phi2*dphidx
                                            tempgy = phi*dphi2dy + phi2*dphidy
                                            tempgz = phi*dphi2dz + phi2*dphidz
                                            O(Jbas,Ibas)=O(Jbas,Ibas)+(temp*dfdr+ &
                                            xdot*tempgx+ydot*tempgy+zdot*tempgz)*weight
!                                        endif
                                    enddo
                                endif
                            enddo
                        endif
                    endif
                enddo
            enddo
        enddo
     endif



! Finally, copy lower diagonal to upper diagonal.


    do Ibas=1,nbasis
        do Jbas=Ibas+1,nbasis
            O(Ibas,Jbas) = O(Jbas,Ibas)
        enddo
    enddo

        call cpu_time(t3)

        write (ioutfile,'(" TIME of evaluation numerical integral = ",F12.2)') &
        T3-T2

!  if(printEnergy)then
!    do Ibas=1,nbasis
!        do Jbas=1,nbasis
!            Eel=Eel+DENSESAVE(Ibas,Jbas)*O(Jbas,Ibas)
!        enddo
!    enddo

    Eel=Eel+Eelxc
! endif

!           write (ioutfile,'("TOTAL ENERGY OF CURRENT CYCLE=",F16.9)') Eel


    END subroutine dftoperatordelta

