! Ed Brothers. October 22,2002
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

subroutine calchessian(failed)
  use allmod
  implicit double precision(a-h,o-z)
  logical failed
  
  failed=.false.

  ! This subroutine calculates the second derivative of energy with respect
  ! to nuclear displacement.  It then uses the analytical Hessian to
  ! optimize geometry if this is an optimization job, and finally calculates
  ! the frequency.  Note that if this is an optimization job it should have
  ! already passed though the LBFGS optimizer before getting here, and thus
  ! requires only refinement.
  
  call PrtAct(ioutfile,"Begin Hessian calculation")

  ! First print out a warning.
  IF ( .NOT. quick_method%opt) &
       write (ioutfile,'(/" WARNING !! FREQUENCIES ONLY VALID AT ", &
       & "OPTIMIZED GEOMETRY!!!")')
       
  ! Now calculate the Hessian.
  IF (quick_method%analhess.and.quick_method%HF) then
    ! Analytical Hessian Matrix, but is broken now
    write (ioutfile,'(/"ANALYTICAL HESSIAN CURRENTLY BROKEN.")')
    IF (quick_method%HF) call HFHessian
  ELSE
    ! Numerical Hessian Matrix
    call fdhessian(failed)
  endif
  
  ! Output Hessian Matrix
  write (ioutfile,'(/"HESSIAN MATRIX ")')
  call PriHessian(ioutfile,3*natom,HESSIAN,'f12.6')
  
  call PrtAct(ioutfile,"Finish Hessian calculation")
end subroutine calchessian


! Ed Brothers. January 21, 2003.
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

subroutine fdhessian(failed)
  use allmod
  implicit double precision(a-h,o-z)

  character(len=1) cartsym(3)
  logical :: failed

  cartsym(1) = 'X'
  cartsym(2) = 'Y'
  cartsym(3) = 'Z'


  ! Finite difference hessian:  When you just can't solve the CPSCF.

  ! Now take a step in each of the cartesian directions, perform an
  ! scf cycle, calculate the gradient, and add that into the Hessian.
  ! Store everything (additional scf output) in the .cphf file. Note that this
  ! is a central difference finite difference.

  stepsize = 1.d-4
  if (master) then
    itemp = ioutfile
    open(iCPHFfile,file=CPHFfilename,status='unknown')
    ioutfile = icphffile
  endif

  DO Iatom = 1,natom
     DO Idirection = 1,3
        xyz(Idirection,Iatom) = xyz(Idirection,Iatom) + stepsize
        call getenergy(failed)
        IF (failed) return

        ! BLOCKED by YIPU MIAO
        IF (quick_method%unrst) THEN
           !                IF (quick_method%HF) call uhfgrad
           !                IF (quick_method%DFT) call udftgrad
           !                IF (quick_method%SEDFT) call usedftgrad
        ELSE
           IF (quick_method%HF) call hfgrad
           !                IF (quick_method%DFT) call dftgrad
           !               IF (quick_method%SEDFT) call sedftgrad
        ENDIF
        Idest = (Iatom-1)*3 + Idirection
        DO Iadd = 1,natom*3
           HESSIAN(Iadd,Idest) = Gradient(Iadd)
        ENDDO

        xyz(Idirection,Iatom) = xyz(Idirection,Iatom)-2.d0*stepsize
        call getenergy(failed)
        IF (failed) return
        IF (quick_method%unrst) THEN
           !                IF (quick_method%HF) call uhfgrad
           !                IF (quick_method%DFT) call udftgrad
           !                IF (quick_method%SEDFT) call usedftgrad
        ELSE
           IF (quick_method%HF) call hfgrad
           !                IF (quick_method%DFT) call dftgrad
           !                IF (quick_method%SEDFT) call sedftgrad
        ENDIF
        Idest = (Iatom-1)*3 + Idirection
        DO Iadd = 1,natom*3
           HESSIAN(Iadd,Idest) =(HESSIAN(Iadd,Idest) - Gradient(Iadd)) &
                /(2.d0*stepsize)
        ENDDO

        ! Finally, return the coordinates to the original location.

        xyz(Idirection,Iatom) = xyz(Idirection,Iatom) + stepsize
     ENDDO
  ENDDO

  if (master) ioutfile=itemp

end subroutine fdhessian




! Ed Brothers. October 22, 2002.
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

subroutine HFHessian
  use allmod
  implicit double precision(a-h,o-z)
  ! dimension W(2*(maxbasis/2)**2,2*(maxbasis/2)**2),
  dimension itype2(3,2),ielecfld(3)
  allocatable W(:,:)
  character(len=1) cartsym(3)

  cartsym(1) = 'X'
  cartsym(2) = 'Y'
  cartsym(3) = 'Z'

  ! The purpose of this subroutine is to calculate the 2nd derivative of
  ! the HF energy with respect to nuclear displacement.  The results
  ! of this are stored in Hessian, which is organized by atom and then
  ! by direction of displacement, i.e. element 1 is the gradient of the
  ! x diplacement of atom 1, element 5 is the y displacement of atom 2,
  ! and Hessian(1,5) is the d2E/dX(1)dY(2).

  ! Please note that there is only one set of code for restricted and
  ! unrestricted HF.

  ! Please also note the Hessian is symmetric.

  DO I=1,3
     ielecfld(I)=0
  ENDDO
  DO Iatm=1,natom*3
     DO Jatm=1,natom*3
        Hessian(Jatm,Iatm)=0.d0
     ENDDO
  ENDDO

  ! Note that the hessian is the sum of a series of term, each explained
  ! when we get to it.
  ! 1)  The 2nd derivative of the nuclear repulsion.

  ! We saw in the gradient code that:
  ! dVnn/dXA = ZA (Sum over B) ZB*(XB-XA) RAB^-3

  ! Using the same math, we can see that:

  ! dVnn/dXAdXA = ZA ZB ((3 (XA-XB)^2 RAB^-5) -  RAB^-3)
  ! dVnn/dXAdXB = ZA ZB ((-3 (XA-XB)^2 RAB^-5) +  RAB^-3)
  ! dVnn/dXAdYA = ZA ZB (3 (XA-XB)(YA-YB) RAB^-5)
  ! dVnn/dXAdYB = ZA ZB (-3 (XA-XB)(YA-YB) RAB^-5)

  ! An interesting fact is that d2Vnn/dXAdYB=-d2Vnn/dXAdYA.  Another
  ! intesting fact is d2Vnn/dXAdXA=d2Vnn/dXBdXB.  We use this
  ! in the next loop.

  DO Iatm=1,natom
     DO Jatm=Iatm+1,natom
        Istart = (Iatm-1)*3
        Jstart = (Jatm-1)*3
        XAminXB = xyz(1,Iatm)-xyz(1,Jatm)
        YAminYB = xyz(2,Iatm)-xyz(2,Jatm)
        ZAminZB = xyz(3,Iatm)-xyz(3,Jatm)
        temp = (XAminXB**2.d0 + YAminYB**2.d0 + ZAminZB**2.d0)
        RAB3 = temp**(-1.5d0)
        RAB5 = temp**(-2.5d0)
        ZA = chg(Iatm)
        ZB = chg(Jatm)

        temp = ZA*ZB*(3.d0*RAB5*XAminXB**2.d0-RAB3)
        Hessian(Istart+1,Istart+1) = temp
        Hessian(Jstart+1,Jstart+1) = temp
        Hessian(Jstart+1,Istart+1) = -temp

        temp = ZA*ZB*(3.d0*RAB5*YAminYB**2.d0-RAB3)
        Hessian(Istart+2,Istart+2) = temp
        Hessian(Jstart+2,Jstart+2) = temp
        Hessian(Jstart+2,Istart+2) = -temp

        temp = ZA*ZB*(3.d0*RAB5*ZAminZB**2.d0-RAB3)
        Hessian(Istart+3,Istart+3) = temp
        Hessian(Jstart+3,Jstart+3) = temp
        Hessian(Jstart+3,Istart+3) = -temp

        temp = ZA*ZB*(3.d0*RAB5*XAminXB*YAminYB)
        Hessian(Istart+2,Istart+1) = temp
        Hessian(Jstart+2,Jstart+1) = temp
        Hessian(Jstart+1,Istart+2) = -temp
        Hessian(Jstart+2,Istart+1) = -temp

        temp = ZA*ZB*(3.d0*RAB5*XAminXB*ZAminZB)
        Hessian(Istart+3,Istart+1) = temp
        Hessian(Jstart+3,Jstart+1) = temp
        Hessian(Jstart+1,Istart+3) = -temp
        Hessian(Jstart+3,Istart+1) = -temp

        temp = ZA*ZB*(3.d0*RAB5*YAminYB*ZAminZB)
        Hessian(Istart+3,Istart+2) = temp
        Hessian(Jstart+3,Jstart+2) = temp
        Hessian(Jstart+2,Istart+3) = -temp
        Hessian(Jstart+3,Istart+2) = -temp
     ENDDO
  ENDDO


  ! 2)  The negative of the energy weighted density matrix element i j
  ! with the second derivative of the ij overlap.

  ! 3)  The second derivative of the 1 electron kinetic energy term ij
  ! times the density matrix element ij.

  ! These terms are grouped together since we loop over the same terms.
  ! Also note that these are the 2-center terms.

  ! First we need to form the energy weighted density matrix.

  ! The energy weighted denisty matrix fora a restricted calculation is:
  ! Q(i,j) =2*(Sum over alpha electrons a)  E(a) C(I,a) C(J,a)
  ! Where C is the alpha or beta molecular orbital coefficients, and
  ! E is the alpha or beta molecular orbital energies.
  ! We'll store this in HOLD as we don't really need it except breifly.

  IF ( .NOT. quick_method%unrst) THEN
     DO I=1,nbasis
        DO J=1,nbasis
           HOLDJI = 0.d0
           DO K=1,nelec/2
              HOLDJI = HOLDJI + (E(K)*CO(J,K)*CO(I,K))
           ENDDO
           HOLD(J,I) = 2.d0*HOLDJI
        ENDDO
     ENDDO

     ! The energy weighted denisty matrix for unrestricted calculation is:
     ! Q(i,j) = (Sum over alpha electrons a)  E(a) C(I,a) C(J,a)
     ! +(Sum over alpha electrons b)  EB(b) CB(I,b) CB(J,b)
     ! Where C is the alpha or beta molecular orbital coefficients, and
     ! E is the alpha or beta molecular orbital energies.
     ! We'll store this in HOLD as we don't really need it (except for hessian
     ! calculations later).

  ELSE
     DO I=1,nbasis
        DO J=1,nbasis
           HOLDJI = 0.d0
           ! Alpha part
           DO K=1,nelec
              HOLDJI = HOLDJI + (E(K)*CO(J,K)*CO(I,K))
           ENDDO
           ! Beta part
           DO K=1,nelecb
              HOLDJI = HOLDJI + (EB(K)*COB(J,K)*COB(I,K))
           ENDDO
           HOLD(J,I) = HOLDJI
        ENDDO
     ENDDO
  ENDIF

  ! Before we begin this, a quick note on the second derivative of Gaussian
  ! orbitals.  If there is only one center, the second derivative is zero.
  ! (As with the first derivative.)  If there are two or more centers, there
  ! are three possibilities.

  ! d/dXA d/dXA ((x-XA)^i (y-YA)^j (z-ZA)^k e^(-ar^2))
  ! = 4a^2((x-XA)^(i+2) (y-YA)^j (z-ZA)^k e^(-ar^2))
  ! - 2a(2i+1)((x-XA)^(i) (y-YA)^j (z-ZA)^k e^(-ar^2))
  ! + (i-1)i((x-XA)^(i-2) (y-YA)^j (z-ZA)^k e^(-ar^2))

  ! d/dXA d/dYA ((x-XA)^i (y-YA)^j (z-ZA)^k e^(-ar^2))
  ! = 4a^2 ((x-XA)^(i+1)(y-YA)^(j+1)(z-ZA)^k e^(-ar^2))
  ! - 2ai ((x-XA)^(i-1)(y-YA)^(j+1)(z-ZA)^k e^(-ar^2))
  ! - 2aj ((x-XA)^(i+1)(y-YA)^(j-1)(z-ZA)^k e^(-ar^2))
  ! -  ij ((x-XA)^(i-1)(y-YA)^(j-1)(z-ZA)^k e^(-ar^2))

  ! d/dXA d/dXB gtfonA gtfonB = dgtfonA/dXA gtfonB + gtfonA dgtfonB/dXB

  ! Note the final case is explained in the gradient code, as it is just a
  ! sum of first derivatives.


  DO Ibas=1,nbasis
     ISTART = (ncenter(Ibas)-1) *3
     DO Jbas=ilast(ncenter(IBAS))+1,nbasis
        JSTART = (ncenter(Jbas)-1) *3
        DENSEJI = DENSE(Jbas,Ibas)
        IF(quick_method%unrst) DENSEJI = DENSEJI+DENSEB(Jbas,Ibas)

        ! We have selected our two basis functions.  First,calculate the
        ! d^2/dXA^2 type terms.

        DO Imomentum=1,3
           d2SI = 0.d0
           d2SJ =0.d0
           d2KEI = 0.d0
           d2KEJ = 0.d0

           ! Do the Ibas derivatives first.

           itype(Imomentum,Ibas) = itype(Imomentum,Ibas)+2
           DO Icon=1,ncontract(Ibas)
              DO Jcon=1,ncontract(Jbas)
                 d2SI = d2SI + 4.d0*aexp(Icon,Ibas)*aexp(Icon,Ibas) &
                      *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                      *overlap(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                      itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                      itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                      xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                      xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                      xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                 d2KEI = d2KEI + 4.d0*aexp(Icon,Ibas)*aexp(Icon,Ibas) &
                      *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                      *ekinetic(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                      itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                      itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                      xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                      xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                      xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
              ENDDO
           ENDDO
           itype(Imomentum,Ibas) = itype(Imomentum,Ibas)-2

           DO Icon=1,ncontract(Ibas)
              DO Jcon=1,ncontract(Jbas)
                 d2SI = d2SI - 2.d0*aexp(Icon,Ibas) &
                      *(1.d0+2.d0*dble(itype(Imomentum,Ibas))) &
                      *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                      *overlap(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                      itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                      itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                      xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                      xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                      xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                 d2KEI = d2KEI - 2.d0*aexp(Icon,Ibas) &
                      *(1.d0+2.d0*dble(itype(Imomentum,Ibas))) &
                      *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                      *ekinetic(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                      itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                      itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                      xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                      xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                      xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
              ENDDO
           ENDDO


           IF (itype(Imomentum,Ibas) >= 2) THEN
              const = dble(itype(Imomentum,Ibas)) &
                   *dble(itype(Imomentum,Ibas)-1)
              itype(Imomentum,Ibas) = itype(Imomentum,Ibas)-2
              DO Icon=1,ncontract(Ibas)
                 DO Jcon=1,ncontract(Jbas)
                    d2SI = d2SI + const* &
                         dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *overlap(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                         itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                         itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                         xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                         xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                         xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                    d2KEI = d2KEI + const* &
                         dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *ekinetic(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                         itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                         itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                         xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                         xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                         xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                 ENDDO
              ENDDO
              itype(Imomentum,Ibas) = itype(Imomentum,Ibas)+2
           ENDIF
           Hessian(ISTART+Imomentum,ISTART+Imomentum) = &
                Hessian(ISTART+Imomentum,ISTART+Imomentum) &
                -d2SI*HOLD(Jbas,Ibas)*2.d0 &
                +d2KeI*DENSEJI*2.d0

           ! Now do the Jbas derivatives.

           itype(Imomentum,Jbas) = itype(Imomentum,Jbas)+2
           DO Icon=1,ncontract(Ibas)
              DO Jcon=1,ncontract(Jbas)
                 d2SJ = d2SJ + 4.d0*aexp(Jcon,Jbas)*aexp(Jcon,Jbas) &
                      *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                      *overlap(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                      itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                      itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                      xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                      xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                      xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                 d2KEJ = d2KEJ + 4.d0*aexp(Jcon,Jbas)*aexp(Jcon,Jbas) &
                      *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                      *ekinetic(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                      itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                      itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                      xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                      xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                      xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
              ENDDO
           ENDDO
           itype(Imomentum,Jbas) = itype(Imomentum,Jbas)-2

           DO Icon=1,ncontract(Ibas)
              DO Jcon=1,ncontract(Jbas)
                 d2SJ = d2SJ - 2.d0*aexp(Jcon,Jbas) &
                      *(1.d0+2.d0*dble(itype(Imomentum,Jbas))) &
                      *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                      *overlap(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                      itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                      itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                      xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                      xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                      xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                 d2KEJ = d2KEJ - 2.d0*aexp(Jcon,Jbas) &
                      *(1.d0+2.d0*dble(itype(Imomentum,Jbas))) &
                      *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                      *ekinetic(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                      itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                      itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                      xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                      xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                      xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
              ENDDO
           ENDDO


           IF (itype(Imomentum,Jbas) >= 2) THEN
              const = dble(itype(Imomentum,Jbas)) &
                   *dble(itype(Imomentum,Jbas)-1)
              itype(Imomentum,Jbas) = itype(Imomentum,Jbas)-2
              DO Icon=1,ncontract(Ibas)
                 DO Jcon=1,ncontract(Jbas)
                    d2SJ = d2SJ + const* &
                         dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *overlap(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                         itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                         itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                         xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                         xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                         xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                    d2KEJ = d2KEJ + const* &
                         dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *ekinetic(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                         itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                         itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                         xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                         xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                         xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                 ENDDO
              ENDDO
              itype(Imomentum,Jbas) = itype(Imomentum,Jbas)+2
           ENDIF
           Hessian(JSTART+Imomentum,JSTART+Imomentum) = &
                Hessian(JSTART+Imomentum,JSTART+Imomentum) &
                -d2SJ*HOLD(Jbas,Ibas)*2.d0 &
                +d2KeJ*DENSEJI*2.d0


           ! Now we are going to do derivatives of the d2/dXAdYA type.  Note that
           ! we are still in the IMOMENTUM loop.

           DO Imomentum2=Imomentum+1,3
              d2SI = 0.d0
              d2SJ = 0.d0
              d2KEI = 0.d0
              d2KEJ = 0.d0

              ! Do the Ibas derivatives first.

              itype(Imomentum,Ibas) = itype(Imomentum,Ibas)+1
              itype(Imomentum2,Ibas) = itype(Imomentum2,Ibas)+1
              DO Icon=1,ncontract(Ibas)
                 DO Jcon=1,ncontract(Jbas)
                    d2SI = d2SI + 4.d0*aexp(Icon,Ibas)*aexp(Icon,Ibas) &
                         *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *overlap(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                         itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                         itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                         xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                         xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                         xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                    d2KEI = d2KEI + 4.d0*aexp(Icon,Ibas)*aexp(Icon,Ibas) &
                         *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *ekinetic(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                         itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                         itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                         xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                         xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                         xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                 ENDDO
              ENDDO
              itype(Imomentum,Ibas) = itype(Imomentum,Ibas)-1
              itype(Imomentum2,Ibas) = itype(Imomentum2,Ibas)-1

              IF (itype(Imomentum,Ibas) /= 0) THEN
                 const = dble(itype(Imomentum,Ibas))
                 itype(Imomentum,Ibas) = itype(Imomentum,Ibas)-1
                 itype(Imomentum2,Ibas) = itype(Imomentum2,Ibas)+1
                 DO Icon=1,ncontract(Ibas)
                    DO Jcon=1,ncontract(Jbas)
                       d2SI = d2SI - 2.d0*aexp(Icon,Ibas)*const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *overlap(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                            itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                            itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                            xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                            xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                            xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                       d2KEI = d2KEI - 2.d0*aexp(Icon,Ibas)*const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *ekinetic(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                            itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                            itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                            xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                            xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                            xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                    ENDDO
                 ENDDO
                 itype(Imomentum,Ibas) = itype(Imomentum,Ibas)+1
                 itype(Imomentum2,Ibas) = itype(Imomentum2,Ibas)-1
              ENDIF

              IF (itype(Imomentum2,Ibas) /= 0) THEN
                 const = dble(itype(Imomentum2,Ibas))
                 itype(Imomentum,Ibas) = itype(Imomentum,Ibas)+1
                 itype(Imomentum2,Ibas) = itype(Imomentum2,Ibas)-1
                 DO Icon=1,ncontract(Ibas)
                    DO Jcon=1,ncontract(Jbas)
                       d2SI = d2SI - 2.d0*aexp(Icon,Ibas)*const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *overlap(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                            itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                            itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                            xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                            xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                            xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                       d2KEI = d2KEI - 2.d0*aexp(Icon,Ibas)*const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *ekinetic(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                            itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                            itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                            xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                            xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                            xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                    ENDDO
                 ENDDO
                 itype(Imomentum,Ibas) = itype(Imomentum,Ibas)-1
                 itype(Imomentum2,Ibas) = itype(Imomentum2,Ibas)+1
              ENDIF

              IF (itype(Imomentum2,Ibas) /= 0 .AND. &
                   itype(Imomentum,Ibas) /= 0) THEN
                 const = dble(itype(Imomentum2,Ibas))* &
                      dble(itype(Imomentum,Ibas))
                 itype(Imomentum,Ibas) = itype(Imomentum,Ibas)-1
                 itype(Imomentum2,Ibas) = itype(Imomentum2,Ibas)-1
                 DO Icon=1,ncontract(Ibas)
                    DO Jcon=1,ncontract(Jbas)
                       d2SI = d2SI +const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *overlap(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                            itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                            itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                            xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                            xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                            xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                       d2KEI = d2KEI + const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *ekinetic(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                            itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                            itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                            xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                            xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                            xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                    ENDDO
                 ENDDO
                 itype(Imomentum,Ibas) = itype(Imomentum,Ibas)+1
                 itype(Imomentum2,Ibas) = itype(Imomentum2,Ibas)+1
              ENDIF

              ! Now do the Jbas derivatives.

              itype(Imomentum,Jbas) = itype(Imomentum,Jbas)+1
              itype(Imomentum2,Jbas) = itype(Imomentum2,Jbas)+1
              DO Icon=1,ncontract(Ibas)
                 DO Jcon=1,ncontract(Jbas)
                    d2SJ = d2SJ + 4.d0*aexp(Jcon,Jbas)*aexp(Jcon,Jbas) &
                         *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *overlap(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                         itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                         itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                         xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                         xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                         xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                    d2KEJ = d2KEJ + 4.d0*aexp(Jcon,Jbas)*aexp(Jcon,Jbas) &
                         *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *ekinetic(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                         itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                         itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                         xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                         xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                         xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                 ENDDO
              ENDDO
              itype(Imomentum,Jbas) = itype(Imomentum,Jbas)-1
              itype(Imomentum2,Jbas) = itype(Imomentum2,Jbas)-1

              IF (itype(Imomentum,Jbas) /= 0) THEN
                 const = dble(itype(Imomentum,Jbas))
                 itype(Imomentum,Jbas) = itype(Imomentum,Jbas)-1
                 itype(Imomentum2,Jbas) = itype(Imomentum2,Jbas)+1
                 DO Icon=1,ncontract(Ibas)
                    DO Jcon=1,ncontract(Jbas)
                       d2SJ = d2SJ - 2.d0*aexp(Jcon,Jbas)*const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *overlap(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                            itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                            itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                            xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                            xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                            xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                       d2KEJ = d2KEJ - 2.d0*aexp(Jcon,Jbas)*const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *ekinetic(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                            itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                            itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                            xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                            xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                            xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                    ENDDO
                 ENDDO
                 itype(Imomentum,Jbas) = itype(Imomentum,Jbas)+1
                 itype(Imomentum2,Jbas) = itype(Imomentum2,Jbas)-1
              ENDIF

              IF (itype(Imomentum2,Jbas) /= 0) THEN
                 const = dble(itype(Imomentum2,Jbas))
                 itype(Imomentum,Jbas) = itype(Imomentum,Jbas)+1
                 itype(Imomentum2,Jbas) = itype(Imomentum2,Jbas)-1
                 DO Icon=1,ncontract(Ibas)
                    DO Jcon=1,ncontract(Jbas)
                       d2SJ = d2SJ - 2.d0*aexp(Jcon,Jbas)*const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *overlap(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                            itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                            itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                            xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                            xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                            xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                       d2KEJ = d2KEJ - 2.d0*aexp(Jcon,Jbas)*const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *ekinetic(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                            itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                            itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                            xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                            xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                            xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                    ENDDO
                 ENDDO
                 itype(Imomentum,Jbas) = itype(Imomentum,Jbas)-1
                 itype(Imomentum2,Jbas) = itype(Imomentum2,Jbas)+1
              ENDIF

              IF (itype(Imomentum2,Jbas) /= 0 .AND. &
                   itype(Imomentum,Jbas) /= 0) THEN
                 const = dble(itype(Imomentum2,Jbas))* &
                      dble(itype(Imomentum,Jbas))
                 itype(Imomentum,Jbas) = itype(Imomentum,Jbas)-1
                 itype(Imomentum2,Jbas) = itype(Imomentum2,Jbas)-1
                 DO Icon=1,ncontract(Ibas)
                    DO Jcon=1,ncontract(Jbas)
                       d2SJ = d2SJ + const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *overlap(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                            itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                            itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                            xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                            xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                            xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                       d2KEJ = d2KEJ + const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *ekinetic(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                            itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                            itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                            xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                            xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                            xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                    ENDDO
                 ENDDO
                 itype(Imomentum,Jbas) = itype(Imomentum,Jbas)+1
                 itype(Imomentum2,Jbas) = itype(Imomentum2,Jbas)+1
              ENDIF

              ! Now add the contributions to the Hessian Array.

              Hessian(ISTART+Imomentum2,ISTART+Imomentum) = &
                   Hessian(ISTART+Imomentum2,ISTART+Imomentum) &
                   -d2SI*HOLD(Jbas,Ibas)*2.d0 &
                   +d2KeI*DENSEJI*2.d0
              Hessian(JSTART+Imomentum2,JSTART+Imomentum) = &
                   Hessian(JSTART+Imomentum2,JSTART+Imomentum) &
                   -d2SJ*HOLD(Jbas,Ibas)*2.d0 &
                   +d2KeJ*DENSEJI*2.d0
           ENDDO

           ! The last part is the d^2/dXAdYB portion.  Note that we are still
           ! inside the IMOMENTUM loop.

           DO Jmomentum=1,3
              d2SIJ = 0.d0
              d2KEIJ = 0.d0

              itype(Imomentum,Ibas) = itype(Imomentum,Ibas)+1
              itype(Jmomentum,Jbas) = itype(Jmomentum,Jbas)+1
              DO Icon=1,ncontract(Ibas)
                 DO Jcon=1,ncontract(Jbas)
                    d2SIJ = d2SIJ + 4.d0*aexp(Icon,Ibas)*aexp(Jcon,Jbas) &
                         *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *overlap(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                         itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                         itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                         xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                         xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                         xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                    d2KEIJ = d2KEIJ + 4.d0*aexp(Icon,Ibas)*aexp(Jcon,Jbas) &
                         *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *ekinetic(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                         itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                         itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                         xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                         xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                         xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                 ENDDO
              ENDDO
              itype(Imomentum,Ibas) = itype(Imomentum,Ibas)-1
              itype(Jmomentum,Jbas) = itype(Jmomentum,Jbas)-1

              IF (itype(Jmomentum,Jbas) /= 0) THEN
                 const = dble(itype(Jmomentum,Jbas))
                 itype(Imomentum,Ibas) = itype(Imomentum,Ibas)+1
                 itype(Jmomentum,Jbas) = itype(Jmomentum,Jbas)-1
                 DO Icon=1,ncontract(Ibas)
                    DO Jcon=1,ncontract(Jbas)
                       d2SIJ = d2SIJ - 2.d0*aexp(Icon,Ibas)*const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *overlap(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                            itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                            itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                            xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                            xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                            xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                       d2KEIJ = d2KEIJ - 2.d0*aexp(Icon,Ibas)*const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *ekinetic(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                            itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                            itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                            xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                            xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                            xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                    ENDDO
                 ENDDO
                 itype(Imomentum,Ibas) = itype(Imomentum,Ibas)-1
                 itype(Jmomentum,Jbas) = itype(Jmomentum,Jbas)+1
              ENDIF

              IF (itype(Imomentum,Ibas) /= 0) THEN
                 const = dble(itype(Imomentum,Ibas))
                 itype(Imomentum,Ibas) = itype(Imomentum,Ibas)-1
                 itype(Jmomentum,Jbas) = itype(Jmomentum,Jbas)+1
                 DO Icon=1,ncontract(Ibas)
                    DO Jcon=1,ncontract(Jbas)
                       d2SIJ = d2SIJ - 2.d0*aexp(Jcon,Jbas)*const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *overlap(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                            itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                            itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                            xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                            xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                            xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                       d2KEIJ = d2KEIJ - 2.d0*aexp(Jcon,Jbas)*const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *ekinetic(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                            itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                            itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                            xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                            xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                            xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                    ENDDO
                 ENDDO
                 itype(Imomentum,Ibas) = itype(Imomentum,Ibas)+1
                 itype(Jmomentum,Jbas) = itype(Jmomentum,Jbas)-1
              ENDIF

              IF (itype(Imomentum,Ibas) /= 0 .AND. &
                   itype(Jmomentum,Jbas) /= 0) THEN
                 const = dble(itype(Imomentum,Ibas))* &
                      dble(itype(Jmomentum,Jbas))
                 itype(Imomentum,Ibas) = itype(Imomentum,Ibas)-1
                 itype(Jmomentum,Jbas) = itype(Jmomentum,Jbas)-1
                 DO Icon=1,ncontract(Ibas)
                    DO Jcon=1,ncontract(Jbas)
                       d2SIJ = d2SIJ +const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *overlap(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                            itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                            itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                            xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                            xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                            xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                       d2KEIJ = d2KEIJ +const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *ekinetic(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                            itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                            itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                            xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                            xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                            xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                    ENDDO
                 ENDDO
                 itype(Imomentum,Ibas) = itype(Imomentum,Ibas)+1
                 itype(Jmomentum,Jbas) = itype(Jmomentum,Jbas)+1
              ENDIF

              ! Now we add the contribution to the Hessian array.

              Hessian(JSTART+Jmomentum,ISTART+Imomentum) = &
                   Hessian(JSTART+Jmomentum,ISTART+Imomentum) &
                   -d2SIJ*HOLD(Jbas,Ibas)*2.d0 &
                   +d2KeIJ*DENSEJI*2.d0

           ENDDO

           ! Quick note:  The three loops closed here are Imomentum,Jbas, and Ibas.
        ENDDO
     ENDDO
  ENDDO

  ! 4)  The second derivative of the 1 electron nuclear attraction term ij
  ! ij times the density matrix element ij.

  ! Please note that these are the three center terms.

  DO Ibas=1,nbasis
     iA=ncenter(Ibas)
     ISTART = (iA-1)*3
     DO I=1,3
        itype2(I,1) = itype(I,Ibas)
     ENDDO

     DO Jbas=Ibas,nbasis
        iB = ncenter(Jbas)
        JSTART = (iB-1)*3
        DO I=1,3
           itype2(I,2) = itype(I,Jbas)
        ENDDO

        DO iC = 1,natom
           iCSTART = (iC-1)*3

           ! As before, if all terms are on the same atom, they move with the
           ! atom and the derivative is zero.

           IF (iA == iC .AND. iB == iC) THEN
              continue
           ELSE
              DENSEJI=DENSE(Jbas,Ibas)
              IF (quick_method%unrst) DENSEJI = DENSEJI+DENSEB(Jbas,Ibas)
              IF (Ibas /= Jbas) DENSEJI=2.d0*DENSEJI
              chgC = chg(iC)

              ! First, take the second derivative of the center the electron is being
              ! attracted to.  These are the d2/dCXdCY and d2/dCX2 derivatives.


              DO ICmomentum=1,3
                 ielecfld(ICmomentum) =  ielecfld(ICmomentum)+1
                 DO JCmomentum=ICmomentum,3
                    ielecfld(JCmomentum) =  ielecfld(JCmomentum)+1
                    d2NAC=0.d0
                    DO Icon=1,ncontract(Ibas)
                       DO Jcon=1,ncontract(Jbas)
                          d2NAC= d2NAC + dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                               *electricfld(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                               itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                               itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                               ielecfld(1),ielecfld(2),ielecfld(3), &
                               xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                               xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                               xyz(1,iC),xyz(2,iC),xyz(3,iC), chgC)
                       ENDDO
                    ENDDO
                    Hessian(ICstart+JCmomentum,ICstart+ICmomentum) = &
                         Hessian(ICstart+JCmomentum,ICstart+ICmomentum)+ &
                         DENSEJI*d2NAC
                    ielecfld(JCmomentum) =  ielecfld(JCmomentum)-1
                 ENDDO

                 ! Now calculate the derivatives of the type d2/dCXdIbasCenterX.  This is
                 ! basically moving the attracting center and one of the centers that a
                 ! basis function is on.

                 ! This is where we begin using the itype2 array.  If Ibas=Jbas and
                 ! we adjust the angular momentum of Ibas in itype, we also inadvertently
                 ! adjust the angular momentum of Jbas.  This leads to large errors.  Thus
                 ! we use a dummy array.

                 ! Note we are still in the ICmomentum loop.  First, loop over Ibas
                 ! momentums.

                 DO IImomentum = 1,3
                    d2NAC=0.d0
                    itype2(IImomentum,1) = itype2(IImomentum,1)+1
                    DO Icon=1,ncontract(Ibas)
                       DO Jcon=1,ncontract(Jbas)
                          d2NAC= d2NAC + 2.d0*aexp(Icon,Ibas)* &
                               dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                               *electricfld(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                               itype2(1,2),itype2(2,2),itype2(3,2), &
                               itype2(1,1),itype2(2,1),itype2(3,1), &
                               ielecfld(1),ielecfld(2),ielecfld(3), &
                               xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                               xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                               xyz(1,iC),xyz(2,iC),xyz(3,iC), chgC)
                       ENDDO
                    ENDDO
                    itype2(IImomentum,1) = itype2(IImomentum,1)-1

                    IF (itype2(IImomentum,1) /= 0) THEN
                       const = dble(itype2(IImomentum,1))
                       itype2(IImomentum,1) = itype2(IImomentum,1)-1
                       DO Icon=1,ncontract(Ibas)
                          DO Jcon=1,ncontract(Jbas)
                             d2NAC= d2NAC -const* &
                                  dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                                  *electricfld(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                                  itype2(1,2),itype2(2,2),itype2(3,2), &
                                  itype2(1,1),itype2(2,1),itype2(3,1), &
                                  ielecfld(1),ielecfld(2),ielecfld(3), &
                                  xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                                  xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                                  xyz(1,iC),xyz(2,iC),xyz(3,iC), chgC)
                          ENDDO
                       ENDDO
                       itype2(IImomentum,1) = itype2(IImomentum,1)+1
                    ENDIF

                    IF (iA == iC .AND. ICmomentum == IImomentum) &
                         d2NAC=d2NAC*2.d0
                    IF (ICstart+ICmomentum >= Istart+IImomentum) THEN
                       Hessian(ICstart+ICmomentum,Istart+IImomentum) = &
                            Hessian(ICstart+ICmomentum,Istart+IImomentum)+ &
                            DENSEJI*d2NAC
                    ELSE
                       Hessian(Istart+IImomentum,ICstart+ICmomentum) = &
                            Hessian(Istart+IImomentum,ICstart+ICmomentum)+ &
                            DENSEJI*d2NAC
                    ENDIF
                 ENDDO

                 ! Now loop over Jbas momentums.

                 DO JJmomentum = 1,3
                    d2NAC=0.d0
                    itype2(JJmomentum,2) = itype2(JJmomentum,2)+1
                    DO Icon=1,ncontract(Ibas)
                       DO Jcon=1,ncontract(Jbas)
                          d2NAC= d2NAC + 2.d0*aexp(Jcon,Jbas)* &
                               dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                               *electricfld(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                               itype2(1,2),itype2(2,2),itype2(3,2), &
                               itype2(1,1),itype2(2,1),itype2(3,1), &
                               ielecfld(1),ielecfld(2),ielecfld(3), &
                               xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                               xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                               xyz(1,iC),xyz(2,iC),xyz(3,iC), chgC)
                       ENDDO
                    ENDDO
                    itype2(JJmomentum,2) = itype2(JJmomentum,2)-1

                    IF (itype2(JJmomentum,2) /= 0) THEN
                       const = dble(itype2(JJmomentum,2))
                       itype2(JJmomentum,2) = itype2(JJmomentum,2)-1
                       DO Icon=1,ncontract(Ibas)
                          DO Jcon=1,ncontract(Jbas)
                             d2NAC= d2NAC -const* &
                                  dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                                  *electricfld(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                                  itype2(1,2),itype2(2,2),itype2(3,2), &
                                  itype2(1,1),itype2(2,1),itype2(3,1), &
                                  ielecfld(1),ielecfld(2),ielecfld(3), &
                                  xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                                  xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                                  xyz(1,iC),xyz(2,iC),xyz(3,iC), chgC)
                          ENDDO
                       ENDDO
                       itype2(JJmomentum,2) = itype2(JJmomentum,2)+1
                    ENDIF

                    IF (iB == iC .AND. ICmomentum == JJmomentum) &
                         d2NAC=d2NAC*2.d0
                    IF (ICstart+ICmomentum >= Jstart+JJmomentum) THEN
                       Hessian(ICstart+ICmomentum,Jstart+JJmomentum) = &
                            Hessian(ICstart+ICmomentum,Jstart+JJmomentum)+ &
                            DENSEJI*d2NAC
                    ELSE
                       Hessian(Jstart+JJmomentum,ICstart+ICmomentum) = &
                            Hessian(Jstart+JJmomentum,ICstart+ICmomentum)+ &
                            DENSEJI*d2NAC
                    ENDIF

                 ENDDO
                 ielecfld(ICmomentum) =  ielecfld(ICmomentum)-1
              ENDDO

              ! Please note we have exited all inner loops at this point and are only
              ! inside the Ibas,Jbas,IC loop here.

              ! At this point we have found all of the elements of the hessian that
              ! involve the attractive center.  Now we perturb the atoms on which the
              ! basis functions lie.  This is exactly analogous to what was done with
              ! the two center integrals above in 2) and 3) with one exception.(d2/dXAdYB)
              ! First,calculate the d^2/dXA^2 type terms.

              DO Imomentum=1,3
                 d2AI = 0.d0
                 d2AJ =0.d0

                 ! Do the Ibas derivatives first.

                 itype2(Imomentum,1) = itype2(Imomentum,1)+2
                 DO Icon=1,ncontract(Ibas)
                    DO Jcon=1,ncontract(Jbas)
                       d2AI = d2AI + 4.d0*aexp(Icon,Ibas)*aexp(Icon,Ibas) &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                            itype2(1,2),itype2(2,2),itype2(3,2), &
                            itype2(1,1),itype2(2,1),itype2(3,1), &
                            xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                            xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                            xyz(1,iC),xyz(2,iC),xyz(3,iC), chgC)
                    ENDDO
                 ENDDO
                 itype2(Imomentum,1) = itype2(Imomentum,1)-2

                 DO Icon=1,ncontract(Ibas)
                    DO Jcon=1,ncontract(Jbas)
                       d2AI = d2AI - 2.d0*aexp(Icon,Ibas) &
                            *(1.d0+2.d0*dble(itype2(Imomentum,1))) &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                            itype2(1,2),itype2(2,2),itype2(3,2), &
                            itype2(1,1),itype2(2,1),itype2(3,1), &
                            xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                            xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                            xyz(1,iC),xyz(2,iC),xyz(3,iC), chgC)
                    ENDDO
                 ENDDO


                 IF (itype2(Imomentum,1) >= 2) THEN
                    const = dble(itype2(Imomentum,1)) &
                         *dble(itype2(Imomentum,1)-1)
                    itype2(Imomentum,1) = itype2(Imomentum,1)-2
                    DO Icon=1,ncontract(Ibas)
                       DO Jcon=1,ncontract(Jbas)
                          d2AI = d2AI + const* &
                               dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                               *attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                               itype2(1,2),itype2(2,2),itype2(3,2), &
                               itype2(1,1),itype2(2,1),itype2(3,1), &
                               xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                               xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                               xyz(1,iC),xyz(2,iC),xyz(3,iC), chgC)
                       ENDDO
                    ENDDO
                    itype2(Imomentum,1) = itype2(Imomentum,1)+2
                 ENDIF
                 Hessian(ISTART+Imomentum,ISTART+Imomentum) = &
                      Hessian(ISTART+Imomentum,ISTART+Imomentum) &
                      +d2AI*DENSEJI

                 ! Now do the Jbas derivatives.

                 itype2(Imomentum,2) = itype2(Imomentum,2)+2
                 DO Icon=1,ncontract(Ibas)
                    DO Jcon=1,ncontract(Jbas)
                       d2AJ = d2AJ + 4.d0*aexp(Jcon,Jbas)*aexp(Jcon,Jbas) &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                            itype2(1,2),itype2(2,2),itype2(3,2), &
                            itype2(1,1),itype2(2,1),itype2(3,1), &
                            xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                            xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                            xyz(1,iC),xyz(2,iC),xyz(3,iC), chgC)
                    ENDDO
                 ENDDO
                 itype2(Imomentum,2) = itype2(Imomentum,2)-2

                 DO Icon=1,ncontract(Ibas)
                    DO Jcon=1,ncontract(Jbas)
                       d2AJ = d2AJ - 2.d0*aexp(Jcon,Jbas) &
                            *(1.d0+2.d0*dble(itype2(Imomentum,2))) &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                            itype2(1,2),itype2(2,2),itype2(3,2), &
                            itype2(1,1),itype2(2,1),itype2(3,1), &
                            xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                            xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                            xyz(1,iC),xyz(2,iC),xyz(3,iC), chgC)
                    ENDDO
                 ENDDO


                 IF (itype2(Imomentum,2) >= 2) THEN
                    const = dble(itype2(Imomentum,2)) &
                         *dble(itype2(Imomentum,2)-1)
                    itype2(Imomentum,2) = itype2(Imomentum,2)-2
                    DO Icon=1,ncontract(Ibas)
                       DO Jcon=1,ncontract(Jbas)
                          d2AJ = d2AJ + const* &
                               dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                               *attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                               itype2(1,2),itype2(2,2),itype2(3,2), &
                               itype2(1,1),itype2(2,1),itype2(3,1), &
                               xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                               xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                               xyz(1,iC),xyz(2,iC),xyz(3,iC), chgC)
                       ENDDO
                    ENDDO
                    itype2(Imomentum,2) = itype2(Imomentum,2)+2
                 ENDIF
                 Hessian(JSTART+Imomentum,JSTART+Imomentum) = &
                      Hessian(JSTART+Imomentum,JSTART+Imomentum) &
                      +d2AJ*DENSEJI


                 ! Now we are going to do derivatives of the d2/dXAdYA type.  Note that
                 ! we are still in the IMOMENTUM loop.

                 DO Imomentum2=Imomentum+1,3
                    d2AI = 0.d0
                    d2AJ = 0.d0

                    ! Do the Ibas derivatives first.

                    itype2(Imomentum,1) = itype2(Imomentum,1)+1
                    itype2(Imomentum2,1) = itype2(Imomentum2,1)+1
                    DO Icon=1,ncontract(Ibas)
                       DO Jcon=1,ncontract(Jbas)
                          d2AI = d2AI + 4.d0*aexp(Icon,Ibas)*aexp(Icon,Ibas) &
                               *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                               *attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                               itype2(1,2),itype2(2,2),itype2(3,2), &
                               itype2(1,1),itype2(2,1),itype2(3,1), &
                               xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                               xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                               xyz(1,iC),xyz(2,iC),xyz(3,iC), chgC)
                       ENDDO
                    ENDDO
                    itype2(Imomentum,1) = itype2(Imomentum,1)-1
                    itype2(Imomentum2,1) = itype2(Imomentum2,1)-1

                    IF (itype2(Imomentum,1) /= 0) THEN
                       const = dble(itype2(Imomentum,1))
                       itype2(Imomentum,1) = itype2(Imomentum,1)-1
                       itype2(Imomentum2,1) = itype2(Imomentum2,1)+1
                       DO Icon=1,ncontract(Ibas)
                          DO Jcon=1,ncontract(Jbas)
                             d2AI = d2AI - 2.d0*aexp(Icon,Ibas)*const &
                                  *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                                  *attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                                  itype2(1,2),itype2(2,2),itype2(3,2), &
                                  itype2(1,1),itype2(2,1),itype2(3,1), &
                                  xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                                  xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                                  xyz(1,iC),xyz(2,iC),xyz(3,iC), chgC)
                          ENDDO
                       ENDDO
                       itype2(Imomentum,1) = itype2(Imomentum,1)+1
                       itype2(Imomentum2,1) = itype2(Imomentum2,1)-1
                    ENDIF

                    IF (itype2(Imomentum2,1) /= 0) THEN
                       const = dble(itype2(Imomentum2,1))
                       itype2(Imomentum,1) = itype2(Imomentum,1)+1
                       itype2(Imomentum2,1) = itype2(Imomentum2,1)-1
                       DO Icon=1,ncontract(Ibas)
                          DO Jcon=1,ncontract(Jbas)
                             d2AI = d2AI - 2.d0*aexp(Icon,Ibas)*const &
                                  *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                                  *attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                                  itype2(1,2),itype2(2,2),itype2(3,2), &
                                  itype2(1,1),itype2(2,1),itype2(3,1), &
                                  xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                                  xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                                  xyz(1,iC),xyz(2,iC),xyz(3,iC), chgC)
                          ENDDO
                       ENDDO
                       itype2(Imomentum,1) = itype2(Imomentum,1)-1
                       itype2(Imomentum2,1) = itype2(Imomentum2,1)+1
                    ENDIF

                    IF (itype2(Imomentum2,1) /= 0 .AND. &
                         itype2(Imomentum,1) /= 0) THEN
                       const = dble(itype2(Imomentum2,1))* &
                            dble(itype2(Imomentum,1))
                       itype2(Imomentum,1) = itype2(Imomentum,1)-1
                       itype2(Imomentum2,1) = itype2(Imomentum2,1)-1
                       DO Icon=1,ncontract(Ibas)
                          DO Jcon=1,ncontract(Jbas)
                             d2AI = d2AI +const &
                                  *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                                  *attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                                  itype2(1,2),itype2(2,2),itype2(3,2), &
                                  itype2(1,1),itype2(2,1),itype2(3,1), &
                                  xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                                  xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                                  xyz(1,iC),xyz(2,iC),xyz(3,iC), chgC)
                          ENDDO
                       ENDDO
                       itype2(Imomentum,1) = itype2(Imomentum,1)+1
                       itype2(Imomentum2,1) = itype2(Imomentum2,1)+1
                    ENDIF


                    ! Now do the Jbas derivatives.

                    itype2(Imomentum,2) = itype2(Imomentum,2)+1
                    itype2(Imomentum2,2) = itype2(Imomentum2,2)+1
                    DO Icon=1,ncontract(Ibas)
                       DO Jcon=1,ncontract(Jbas)
                          d2AJ = d2AJ + 4.d0*aexp(Jcon,Jbas)*aexp(Jcon,Jbas) &
                               *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                               *attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                               itype2(1,2),itype2(2,2),itype2(3,2), &
                               itype2(1,1),itype2(2,1),itype2(3,1), &
                               xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                               xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                               xyz(1,iC),xyz(2,iC),xyz(3,iC), chgC)
                       ENDDO
                    ENDDO
                    itype2(Imomentum,2) = itype2(Imomentum,2)-1
                    itype2(Imomentum2,2) = itype2(Imomentum2,2)-1

                    IF (itype2(Imomentum,2) /= 0) THEN
                       const = dble(itype2(Imomentum,2))
                       itype2(Imomentum,2) = itype2(Imomentum,2)-1
                       itype2(Imomentum2,2) = itype2(Imomentum2,2)+1
                       DO Icon=1,ncontract(Ibas)
                          DO Jcon=1,ncontract(Jbas)
                             d2AJ = d2AJ - 2.d0*aexp(Jcon,Jbas)*const &
                                  *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                                  *attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                                  itype2(1,2),itype2(2,2),itype2(3,2), &
                                  itype2(1,1),itype2(2,1),itype2(3,1), &
                                  xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                                  xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                                  xyz(1,iC),xyz(2,iC),xyz(3,iC), chgC)
                          ENDDO
                       ENDDO
                       itype2(Imomentum,2) = itype2(Imomentum,2)+1
                       itype2(Imomentum2,2) = itype2(Imomentum2,2)-1
                    ENDIF

                    IF (itype2(Imomentum2,2) /= 0) THEN
                       const = dble(itype2(Imomentum2,2))
                       itype2(Imomentum,2) = itype2(Imomentum,2)+1
                       itype2(Imomentum2,2) = itype2(Imomentum2,2)-1
                       DO Icon=1,ncontract(Ibas)
                          DO Jcon=1,ncontract(Jbas)
                             d2AJ = d2AJ - 2.d0*aexp(Jcon,Jbas)*const &
                                  *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                                  *attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                                  itype2(1,2),itype2(2,2),itype2(3,2), &
                                  itype2(1,1),itype2(2,1),itype2(3,1), &
                                  xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                                  xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                                  xyz(1,iC),xyz(2,iC),xyz(3,iC), chgC)
                          ENDDO
                       ENDDO
                       itype2(Imomentum,2) = itype2(Imomentum,2)-1
                       itype2(Imomentum2,2) = itype2(Imomentum2,2)+1
                    ENDIF

                    IF (itype2(Imomentum2,2) /= 0 .AND. &
                         itype2(Imomentum,2) /= 0) THEN
                       const = dble(itype2(Imomentum2,2))* &
                            dble(itype2(Imomentum,2))
                       itype2(Imomentum,2) = itype2(Imomentum,2)-1
                       itype2(Imomentum2,2) = itype2(Imomentum2,2)-1
                       DO Icon=1,ncontract(Ibas)
                          DO Jcon=1,ncontract(Jbas)
                             d2AJ = d2AJ + const &
                                  *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                                  *attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                                  itype2(1,2),itype2(2,2),itype2(3,2), &
                                  itype2(1,1),itype2(2,1),itype2(3,1), &
                                  xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                                  xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                                  xyz(1,iC),xyz(2,iC),xyz(3,iC), chgC)
                          ENDDO
                       ENDDO
                       itype2(Imomentum,2) = itype2(Imomentum,2)+1
                       itype2(Imomentum2,2) = itype2(Imomentum2,2)+1
                    ENDIF

                    ! Now add the contributions to the Hessian Array.

                    Hessian(ISTART+Imomentum2,ISTART+Imomentum) = &
                         Hessian(ISTART+Imomentum2,ISTART+Imomentum) &
                         +d2AI*DENSEJI
                    Hessian(JSTART+Imomentum2,JSTART+Imomentum) = &
                         Hessian(JSTART+Imomentum2,JSTART+Imomentum) &
                         +d2AJ*DENSEJI
                 ENDDO

                 ! Close the Imomentum2 loop.

                 ! The last part is the d^2/dXAdYB portion.  Note that we are still
                 ! inside the IMOMENTUM loop.



                 DO Jmomentum=1,3
                    d2AIJ = 0.d0

                    itype2(Imomentum,1) = itype2(Imomentum,1)+1
                    itype2(Jmomentum,2) = itype2(Jmomentum,2)+1
                    DO Icon=1,ncontract(Ibas)
                       DO Jcon=1,ncontract(Jbas)
                          d2AIJ = d2AIJ + 4.d0*aexp(Icon,Ibas)*aexp(Jcon,Jbas) &
                               *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                               *attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                               itype2(1,2),itype2(2,2),itype2(3,2), &
                               itype2(1,1),itype2(2,1),itype2(3,1), &
                               xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                               xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                               xyz(1,iC),xyz(2,iC),xyz(3,iC), chgC)
                       ENDDO
                    ENDDO
                    itype2(Imomentum,1) = itype2(Imomentum,1)-1
                    itype2(Jmomentum,2) = itype2(Jmomentum,2)-1

                    IF (itype2(Jmomentum,2) /= 0) THEN
                       const = dble(itype2(Jmomentum,2))
                       itype2(Imomentum,1) = itype2(Imomentum,1)+1
                       itype2(Jmomentum,2) = itype2(Jmomentum,2)-1
                       DO Icon=1,ncontract(Ibas)
                          DO Jcon=1,ncontract(Jbas)
                             d2AIJ = d2AIJ - 2.d0*aexp(Icon,Ibas)*const &
                                  *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                                  *attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                                  itype2(1,2),itype2(2,2),itype2(3,2), &
                                  itype2(1,1),itype2(2,1),itype2(3,1), &
                                  xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                                  xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                                  xyz(1,iC),xyz(2,iC),xyz(3,iC), chgC)
                          ENDDO
                       ENDDO
                       itype2(Imomentum,1) = itype2(Imomentum,1)-1
                       itype2(Jmomentum,2) = itype2(Jmomentum,2)+1
                    ENDIF

                    IF (itype2(Imomentum,1) /= 0) THEN
                       const = dble(itype2(Imomentum,1))
                       itype2(Imomentum,1) = itype2(Imomentum,1)-1
                       itype2(Jmomentum,2) = itype2(Jmomentum,2)+1
                       DO Icon=1,ncontract(Ibas)
                          DO Jcon=1,ncontract(Jbas)
                             d2AIJ = d2AIJ - 2.d0*aexp(Jcon,Jbas)*const &
                                  *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                                  *attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                                  itype2(1,2),itype2(2,2),itype2(3,2), &
                                  itype2(1,1),itype2(2,1),itype2(3,1), &
                                  xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                                  xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                                  xyz(1,iC),xyz(2,iC),xyz(3,iC), chgC)
                          ENDDO
                       ENDDO
                       itype2(Imomentum,1) = itype2(Imomentum,1)+1
                       itype2(Jmomentum,2) = itype2(Jmomentum,2)-1
                    ENDIF

                    IF (itype2(Imomentum,1) /= 0 .AND. &
                         itype2(Jmomentum,2) /= 0) THEN
                       const = dble(itype2(Imomentum,1))* &
                            dble(itype2(Jmomentum,2))
                       itype2(Imomentum,1) = itype2(Imomentum,1)-1
                       itype2(Jmomentum,2) = itype2(Jmomentum,2)-1
                       DO Icon=1,ncontract(Ibas)
                          DO Jcon=1,ncontract(Jbas)
                             d2AIJ = d2AIJ +const &
                                  *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                                  *attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                                  itype2(1,2),itype2(2,2),itype2(3,2), &
                                  itype2(1,1),itype2(2,1),itype2(3,1), &
                                  xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                                  xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                                  xyz(1,iC),xyz(2,iC),xyz(3,iC), chgC)
                          ENDDO
                       ENDDO
                       itype2(Imomentum,1) = itype2(Imomentum,1)+1
                       itype2(Jmomentum,2) = itype2(Jmomentum,2)+1
                    ENDIF

                    ! Now we add the contribution to the Hessian array.

                    IF (iA /= iB) THEN
                       Hessian(JSTART+Jmomentum,ISTART+Imomentum) = &
                            Hessian(JSTART+Jmomentum,ISTART+Imomentum) &
                            +d2AIJ*DENSEJI
                    ELSE
                       IF (Imomentum == Jmomentum) THEN
                          Hessian(JSTART+Jmomentum,ISTART+Imomentum) = &
                               Hessian(JSTART+Jmomentum,ISTART+Imomentum) &
                               +2.d0*d2AIJ*DENSEJI
                       ELSEIF (Jmomentum > Imomentum) THEN
                          Hessian(JSTART+Jmomentum,ISTART+Imomentum) = &
                               Hessian(JSTART+Jmomentum,ISTART+Imomentum) &
                               +d2AIJ*DENSEJI
                       ELSE
                          Hessian(ISTART+Imomentum,JSTART+Jmomentum) = &
                               Hessian(ISTART+Imomentum,JSTART+Jmomentum) &
                               +d2AIJ*DENSEJI
                       ENDIF
                    ENDIF

                    ! Now close the Jmomentum and Imomentum loops.

                 ENDDO
              ENDDO

           ENDIF
        ENDDO
     ENDDO
  ENDDO


  ! 5)  The 2nd derivative of the 4center 2e- terms with respect to X times
  ! the coefficient found in the energy. (i.e. the multiplicative
  ! constants from the density matrix that arise as these are both
  ! the exchange and correlation integrals.)


  DO I=1,nbasis
     ! Set some variables to reduce access time for some of the more
     ! used quantities.

     DENSEII=DENSE(I,I)
     IF (quick_method%unrst) DENSEII=DENSEII+DENSEB(I,I)

     ! Neglect all the (ii|ii) integrals, as they move with the core.

     DO J=I+1,nbasis
        ! Set some variables to reduce access time for some of the more
        ! used quantities. (AGAIN)

        DENSEJI=DENSE(J,I)
        DENSEJJ=DENSE(J,J)
        IF (quick_method%unrst) THEN
           DENSEJI=DENSEJI+DENSEB(J,I)
           DENSEJJ=DENSEJJ+DENSEB(J,J)
        ENDIF

        ! Find  all the (ii|jj) integrals.

        constant = (DENSEII*DENSEJJ-.5d0*DENSEJI*DENSEJI)

        call hess2elec(I,I,J,J,constant)

        ! Find  all the (ij|jj) integrals.

        constant =  DENSEJJ*DENSEJI
        call hess2elec(I,J,J,J,constant)


        ! Find  all the (ii|ij) integrals.
        constant= DENSEJI*DENSEII
        call hess2elec(I,I,I,J,constant)

        ! Find all the (ij|ij) integrals
        constant =(1.5d0*DENSEJI*DENSEJI-0.50d0*DENSEJJ*DENSEII)
        call hess2elec(I,J,I,J,constant)

        DO K=J+1,nbasis
           ! Set some variables to reduce access time for some of the more
           ! used quantities. (AGAIN)

           DENSEKI=DENSE(K,I)
           DENSEKJ=DENSE(K,J)
           DENSEKK=DENSE(K,K)
           IF (quick_method%unrst) THEN
              DENSEKI=DENSEKI+DENSEB(K,I)
              DENSEKJ=DENSEKJ+DENSEB(K,J)
              DENSEKK=DENSEKK+DENSEB(K,K)
           ENDIF

           ! Find all the (ij|ik) integrals where j>i,k>j

           constant = (3.0d0*DENSEJI*DENSEKI-DENSEKJ*DENSEII)
           call hess2elec(I,J,I,K,constant)

           ! Find all the (ij|kk) integrals where j>i, k>j.

           constant=(2.d0*DENSEJI*DENSEKK-DENSEKI*DENSEKJ)
           call hess2elec(I,J,K,K,constant)

           ! Find all the (ik|jj) integrals where j>i, k>j.

           constant= (2.d0*DENSEKI*DENSEJJ-DENSEKJ*DENSEJI)
           call hess2elec(I,K,J,J,constant)

           ! Find all the (ii|jk) integrals where j>i, k>j.

           constant = (2.d0*DENSEKJ*DENSEII-DENSEJI*DENSEKI)
           call hess2elec(I,I,J,K,constant)
        ENDDO

        DO K=I+1,nbasis-1
           DENSEKI=DENSE(K,I)
           DENSEKJ=DENSE(K,J)
           DENSEKK=DENSE(K,K)
           IF (quick_method%unrst) THEN
              DENSEKI=DENSEKI+DENSEB(K,I)
              DENSEKJ=DENSEKJ+DENSEB(K,J)
              DENSEKK=DENSEKK+DENSEB(K,K)
           ENDIF

           DO L=K+1,nbasis
              DENSELJ=DENSE(L,J)
              DENSELI=DENSE(L,I)
              DENSELK=DENSE(L,K)
              IF (quick_method%unrst) THEN
                 DENSELJ=DENSELJ+DENSEB(L,J)
                 DENSELI=DENSELI+DENSEB(L,I)
                 DENSELK=DENSELK+DENSEB(L,K)
              ENDIF

              ! Find the (ij|kl) integrals where j>i,k>i,l>k. Note that k and j
              ! can be equal.

              constant = (4.d0*DENSEJI*DENSELK-DENSEKI*DENSELJ &
                   -DENSELI*DENSEKJ)
              call hess2elec(I,J,K,L,constant)

           ENDDO
        ENDDO
     ENDDO
  ENDDO

  !  DO I=1,natom*3
  !     DO J=I,natom*3
  !        print *,'INTEGRALS',HESSIAN(J,I)
  !     ENDDO
  !  ENDDO


  ! At this point we have all the second derivatives of the energy
  ! in terms of guassian basis integrals.
  ! Now we need to form the first derivative of the energy
  ! weighted density matrix and the bond order density matrix.  We will
  ! be following the procedure detailed in Pople et. al., Int. J. Quant.
  ! Chem. 13, 225-241, 1979.

  ! Basically we need to find some array that produces the derivatives of
  ! the MO coefficients.  This array is called B, and is calculated from
  ! the set of linear equations (I-A)B=B0.  Thus we are going to form A
  ! and all the B0, and then find B thus finding the the derivatives of
  ! the density matrix and the energy weighted density matrix.

  ! First we need to form the CPHFA array.  Before we define it, a note
  ! about notation.  a and b will be used to denote virtual orbitals,
  ! and i and j will denote occupied orbitals.  The subscript ai is one
  ! subscript referring to the pairing of a and i.  If there are 4 occ
  ! and 3 virtual orbitals, and we want to find location a=2 i=3, it would
  ! be the location 7.  (a=1 is 1-4, a=2 i=1 is 5, a=2 i=2 is 6, etc.)

  ! The CPHFA is:

  ! A(ai,bj) = (E(i)-E(a))^-1 (2 (ai|bj) ) - (aj|bi) - (ab|ij)

  ! The first step is building up the molecular orbitals from atomic orbitals
  ! and putting them in the right location.  We'll then go through and
  ! divide by the energies.

  ! Since we have to consider spin, set up COB and EB if this is RHF.

  IF ( .NOT. quick_method%unrst) THEN
     DO I=1,nbasis
        EB(I)=E(I)
        DO J=1,nbasis
           COB(J,I)=CO(J,I)
        ENDDO
     ENDDO
  ENDIF

  ! Blank out the CPHFA array.

  IF (quick_method%unrst) THEN
     idimA = (nbasis-nelec)*nelec + (nbasis-nelecB)*nelecB
  ELSE
     idimA = 2*(nbasis-(nelec/2))*(nelec/2)
  ENDIF
  
  allocate(CPHFA(idimA,idimA))
  allocate(CPHFB(idimA,natom*3))
  allocate(W(idimA,idimA))
  allocate(B0(idimA))
  allocate(BU(idimA))


  DO I=1,idimA
     DO J=1,idimA
        CPHFA(J,I)=0.d0
     ENDDO
  ENDDO

  ! We now pass all the nonredundant AO repulsion integrals to a subroutine.
  ! Note that this series of calls is not well commented, but is exactly the
  ! same order as is found in hfenergy or hfgrad.

  call cpu_time(t1)
  DO I=1,nbasis
     call formCPHFA(I,I,I,I)
     DO J=I+1,nbasis
        call formCPHFA(I,I,J,J)
        call formCPHFA(I,J,J,J)
        call formCPHFA(I,I,I,J)
        call formCPHFA(I,J,I,J)
        DO K=J+1,nbasis
           call formCPHFA(I,J,I,K)
           call formCPHFA(I,J,K,K)
           call formCPHFA(I,K,J,J)
           call formCPHFA(I,I,J,K)
        ENDDO
        DO K=I+1,nbasis-1
           DO L=K+1,nbasis
              call formCPHFA(I,J,K,L)
         ENDDO
        ENDDO
     ENDDO
  ENDDO
  ! At this point, CPHFA(ai,bj) = 2(ai|bj)-(aj|bi)-(ab|ij)
  ! We need to go through and divide by E(I)-E(A) (or EB(I)-EB(A))

  IF (quick_method%unrst) THEN
     lastAocc = nelec
     lastBocc = nelecb
  ELSE
     lastAocc = nelec/2
     lastBocc = lastAocc
  ENDIF
  iBetastart = lastAocc*(nbasis-lastAocc)

  DO iAvirt = lastAocc+1,nbasis
     DO iAocc = 1,lastAocc
        iaCPHFA = (iAvirt-lastAocc-1)*lastAocc + iAocc
        denom = E(iAocc)-E(iAvirt)
        DO jbCPHFA = 1,idimA
           CPHFA(iaCPHFA,jbCPHFA) = CPHFA(iaCPHFA,jbCPHFA)/denom
        ENDDO
     ENDDO
  ENDDO
  DO iBvirt = lastBocc+1,nbasis
     DO iBocc = 1,lastBocc
        iaCPHFA = (iBvirt-lastBocc-1)*lastBocc + iBocc+iBetastart
        denom = EB(iBocc)-EB(iBvirt)
        DO jbCPHFA = 1,idimA
           CPHFA(iaCPHFA,jbCPHFA) = CPHFA(iaCPHFA,jbCPHFA)/denom
        ENDDO
     ENDDO
  ENDDO
  call cpu_time(t2)
  !  print *,'FORM1',T2-T1

  ! APPEARS PERFECT TO HERE. DEBUGGED.
  ! Now we are going to form all of the B0 in an array.  There is one B0 for
  ! each possible nuclear perturbation.  Thus the array is organized with
  ! CPHFB(occ-virtual pair, dX)  This is done in a subroutine.

  DO I=1,natom*3
     DO J=1,idimA
        CPHFB(J,I)=0.d0
     ENDDO
  ENDDO

  call formCPHFB
  call cpu_time(t3)
!    print *,'FORM2',T3-T2

  ! NOTE:  THERE ARE BETTER CONVERGERS.  USE THEM LATER.  I WANT TO GRADUATE.


  ! Now we are going to solve the CPHF  equation.


  DO IdX=1,natom*3

     ! To solve the CPHF equations we are going to need the transpose of A many
     ! times.  Place it  in array W.

     DO I=1,idimA
        DO J=1,idimA
           W(J,I) = CPHFA(I,J)
        ENDDO
     ENDDO
     
     ! Place the correct column of CPHFB in B0.  Also, set BU to BO.

     DO J=1,idimA
        B0(J) = CPHFB(J,IDX)
        BU(J) = CPHFB(J,IDX)
     ENDDO

     ! Solve the equation. BU = (Sum over i=0,N) A^i BO.  This section assumes
     ! on first pass that CPHFA contains the actual A Array, W contain
     ! Transpose[CPHFA].
     change=1.D10
     iBUform=0
     DO WHILE (change.gt.1.D-10)
        iBUform=iBUform+1
        change=-1.d0
        DO I=1,idimA
           element = 0.0D0
           DO K=1,idimA
              element = element + W(K,I)*B0(K)
           ENDDO
           change=max(change,dabs(element))
           BU(I)=BU(I)+element
        ENDDO

        ! At this point the have formed the iBUformth BU array.
        ! If it hasn't converged,
        ! we need to form A^iBUform+1.  We do this by multiplying W
        ! (A^iBUform transposed)
        ! by A and storing it in the CPHFfile, then transposing the file onto W.


        IF (change > 1.D-10) THEN
           open(iCPHFfile,file=CPHFfilename,status='unknown')
           DO I=1,idimA
              DO J=1,idimA
                 W2IJ = 0.0D0
                 DO K=1,idimA
                    W2IJ = W2IJ + W(K,I)*CPHFA(K,J)
                 ENDDO
                 write (iCPHFfile,'(E30.20)') W2IJ
              ENDDO
           ENDDO
           close  (iCPHFfile)
           open(iCPHFfile,file=CPHFfilename,status='unknown')
           DO I=1,idimA
              DO J=1,idimA
                 read (iCPHFfile,'(E30.20)') W(J,I)
              ENDDO
           ENDDO
           close(iCPHFfile)
        ENDIF
        
     ENDDO

     ! BU is now filled with the u(ai) values need to for the first derivative
     ! of the density matrix and the first derivative of the energy weighted
     ! density matrix.  This is done in two subprograms.  The first of these
     ! forms the first derivative of the density matrix and adds the contribution
     ! to the Hessian.  This is fairly simple.
     
     call dmxderiv(IDX)

     ! At this point we now have the density matrix derivatives.  Now we need to
     ! use them with the first derivatives of the integrals to form another
     ! part of the hessian.

     call hfdmxderuse(IDX)

     call Ewtdmxder(IDX)
  ENDDO

  ! At this point some of the elements above the diagonal.  Sum those into
  ! below the diagonal and then set the whol thing to be symmetric.

  DO I=1,natom*3
     DO J=I+1,natom*3
        HESSIAN(I,J) = HESSIAN(J,I)
     ENDDO
  ENDDO

END subroutine hfhessian


! Ed Brothers. November 5, 2002.
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

subroutine hess2elec(Ibas,Jbas,IIbas,JJbas,coeff)
  use allmod
  implicit double precision(a-h,o-z)

  dimension itype2(3,4)
  logical :: same
  same = .false.

  ! The purpose of this subroutine is to calculate the second derivative of
  ! the 2-electron 4-center integrals and add them into the total hessian.
  ! This requires the use of a multiplicitave constant (coeff) which arises
  ! by the way the integral enters the total energy. For example, (ij|ij) is
  ! a repulsion integral and an exchange integral, thus it enters the energy
  ! with a coefficient of (1.5d0*DENSEJI*DENSEJI-0.50d0*DENSEJJ*DENSEII).


  ! First, find the centers the basis functions are located on.  If all the
  ! functions are on the same center, return as this is a zero result.

  iA = ncenter(Ibas)
  iB = ncenter(Jbas)
  iC = ncenter(IIbas)
  iD = ncenter(JJbas)

  same = iA.eq.iB
  same = same .and. iB.eq.iC
  same = same .and. iC.eq.iD

  IF (same) return

  iAstart = (iA-1)*3
  iBstart = (iB-1)*3
  iCstart = (iC-1)*3
  iDstart = (iD-1)*3

  ! The itype2 array was added because if Ibas=Jbas, the code raises two
  ! angular momentums instead of one.

  DO Imomentum=1,3
     itype2(Imomentum,1) = itype(Imomentum,Ibas)
     itype2(Imomentum,2) = itype(Imomentum,Jbas)
     itype2(Imomentum,3) = itype(Imomentum,IIbas)
     itype2(Imomentum,4) = itype(Imomentum,JJbas)
  ENDDO

  ! The first thing to calculate is the d2/dXA elements of the hessian.

  DO Imomentum=1,3
     D2I = 0.d0
     D2J = 0.d0
     D2II = 0.d0
     D2JJ = 0.d0

     ! Ibas' center first.

     itype2(Imomentum,1) = itype2(Imomentum,1)+2
     DO Icon=1,ncontract(Ibas)
        DO Jcon=1,ncontract(Jbas)
           DO IIcon=1,ncontract(IIbas)
              DO JJcon=1,ncontract(JJbas)
                 d2I = d2I + 4.d0*aexp(Icon,Ibas)*aexp(Icon,Ibas) &
                      *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                      *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                      *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                      aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                      itype2(1,1), itype2(2,1), itype2(3,1), &
                      itype2(1,2), itype2(2,2), itype2(3,2), &
                      itype2(1,3), itype2(2,3), itype2(3,3), &
                      itype2(1,4), itype2(2,4), itype2(3,4), &
                      xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                      xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                      xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                      xyz(1,iD), xyz(2,iD),xyz(3,iD))
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     itype2(Imomentum,1) = itype2(Imomentum,1)-2

     DO Icon=1,ncontract(Ibas)
        DO Jcon=1,ncontract(Jbas)
           DO IIcon=1,ncontract(IIbas)
              DO JJcon=1,ncontract(JJbas)
                 d2I = d2I - 2.d0*aexp(Icon,Ibas) &
                      *(1.d0+2.d0*dble(itype2(Imomentum,1))) &
                      *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                      *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                      *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                      aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                      itype2(1,1), itype2(2,1), itype2(3,1), &
                      itype2(1,2), itype2(2,2), itype2(3,2), &
                      itype2(1,3), itype2(2,3), itype2(3,3), &
                      itype2(1,4), itype2(2,4), itype2(3,4), &
                      xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                      xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                      xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                      xyz(1,iD), xyz(2,iD),xyz(3,iD))
              ENDDO
           ENDDO
        ENDDO
     ENDDO


     IF (itype2(Imomentum,1) >= 2) THEN
        const = dble(itype2(Imomentum,1)) &
             *dble(itype2(Imomentum,1)-1)
        itype2(Imomentum,1) = itype2(Imomentum,1)-2
        DO Icon=1,ncontract(Ibas)
           DO Jcon=1,ncontract(Jbas)
              DO IIcon=1,ncontract(IIbas)
                 DO JJcon=1,ncontract(JJbas)
                    d2I = d2I + const* &
                         dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                         *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                         aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                         itype2(1,1), itype2(2,1), itype2(3,1), &
                         itype2(1,2), itype2(2,2), itype2(3,2), &
                         itype2(1,3), itype2(2,3), itype2(3,3), &
                         itype2(1,4), itype2(2,4), itype2(3,4), &
                         xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                         xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                         xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                         xyz(1,iD), xyz(2,iD),xyz(3,iD))
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        itype2(Imomentum,1) = itype2(Imomentum,1)+2
     ENDIF
     Hessian(iASTART+Imomentum,iASTART+Imomentum) = &
          Hessian(iASTART+Imomentum,iASTART+Imomentum) &
          +d2I*coeff

     ! Jbas' center.

     itype2(Imomentum,2) = itype2(Imomentum,2)+2
     DO Icon=1,ncontract(Ibas)
        DO Jcon=1,ncontract(Jbas)
           DO IIcon=1,ncontract(IIbas)
              DO JJcon=1,ncontract(JJbas)
                 d2J = d2J + 4.d0*aexp(Jcon,Jbas)*aexp(Jcon,Jbas) &
                      *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                      *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                      *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                      aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                      itype2(1,1), itype2(2,1), itype2(3,1), &
                      itype2(1,2), itype2(2,2), itype2(3,2), &
                      itype2(1,3), itype2(2,3), itype2(3,3), &
                      itype2(1,4), itype2(2,4), itype2(3,4), &
                      xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                      xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                      xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                      xyz(1,iD), xyz(2,iD),xyz(3,iD))
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     itype2(Imomentum,2) = itype2(Imomentum,2)-2

     DO Icon=1,ncontract(Ibas)
        DO Jcon=1,ncontract(Jbas)
           DO IIcon=1,ncontract(IIbas)
              DO JJcon=1,ncontract(JJbas)
                 d2J = d2J - 2.d0*aexp(Jcon,Jbas) &
                      *(1.d0+2.d0*dble(itype2(Imomentum,2))) &
                      *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                      *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                      *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                      aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                      itype2(1,1), itype2(2,1), itype2(3,1), &
                      itype2(1,2), itype2(2,2), itype2(3,2), &
                      itype2(1,3), itype2(2,3), itype2(3,3), &
                      itype2(1,4), itype2(2,4), itype2(3,4), &
                      xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                      xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                      xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                      xyz(1,iD), xyz(2,iD),xyz(3,iD))
              ENDDO
           ENDDO
        ENDDO
     ENDDO


     IF (itype2(Imomentum,2) >= 2) THEN
        const = dble(itype2(Imomentum,2)) &
             *dble(itype2(Imomentum,2)-1)
        itype2(Imomentum,2) = itype2(Imomentum,2)-2
        DO Icon=1,ncontract(Ibas)
           DO Jcon=1,ncontract(Jbas)
              DO IIcon=1,ncontract(IIbas)
                 DO JJcon=1,ncontract(JJbas)
                    d2J = d2J + const* &
                         dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                         *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                         aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                         itype2(1,1), itype2(2,1), itype2(3,1), &
                         itype2(1,2), itype2(2,2), itype2(3,2), &
                         itype2(1,3), itype2(2,3), itype2(3,3), &
                         itype2(1,4), itype2(2,4), itype2(3,4), &
                         xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                         xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                         xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                         xyz(1,iD), xyz(2,iD),xyz(3,iD))
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        itype2(Imomentum,2) = itype2(Imomentum,2)+2
     ENDIF
     Hessian(iBSTART+Imomentum,iBSTART+Imomentum) = &
          Hessian(iBSTART+Imomentum,iBSTART+Imomentum) &
          +d2J*coeff
     ! IIbas' center.

     itype2(Imomentum,3) = itype2(Imomentum,3)+2
     DO Icon=1,ncontract(Ibas)
        DO Jcon=1,ncontract(Jbas)
           DO IIcon=1,ncontract(IIbas)
              DO JJcon=1,ncontract(JJbas)
                 d2II = d2II + 4.d0*aexp(IIcon,IIbas)*aexp(IIcon,IIbas) &
                      *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                      *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                      *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                      aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                      itype2(1,1), itype2(2,1), itype2(3,1), &
                      itype2(1,2), itype2(2,2), itype2(3,2), &
                      itype2(1,3), itype2(2,3), itype2(3,3), &
                      itype2(1,4), itype2(2,4), itype2(3,4), &
                      xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                      xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                      xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                      xyz(1,iD), xyz(2,iD),xyz(3,iD))
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     itype2(Imomentum,3) = itype2(Imomentum,3)-2

     DO Icon=1,ncontract(Ibas)
        DO Jcon=1,ncontract(Jbas)
           DO IIcon=1,ncontract(IIbas)
              DO JJcon=1,ncontract(JJbas)
                 d2II = d2II - 2.d0*aexp(IIcon,IIbas) &
                      *(1.d0+2.d0*dble(itype2(Imomentum,3))) &
                      *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                      *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                      *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                      aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                      itype2(1,1), itype2(2,1), itype2(3,1), &
                      itype2(1,2), itype2(2,2), itype2(3,2), &
                      itype2(1,3), itype2(2,3), itype2(3,3), &
                      itype2(1,4), itype2(2,4), itype2(3,4), &
                      xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                      xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                      xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                      xyz(1,iD), xyz(2,iD),xyz(3,iD))
              ENDDO
           ENDDO
        ENDDO
     ENDDO


     IF (itype2(Imomentum,3) >= 2) THEN
        const = dble(itype2(Imomentum,3)) &
             *dble(itype2(Imomentum,3)-1)
        itype2(Imomentum,3) = itype2(Imomentum,3)-2
        DO Icon=1,ncontract(Ibas)
           DO Jcon=1,ncontract(Jbas)
              DO IIcon=1,ncontract(IIbas)
                 DO JJcon=1,ncontract(JJbas)
                    d2II = d2II + const* &
                         dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                         *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                         aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                         itype2(1,1), itype2(2,1), itype2(3,1), &
                         itype2(1,2), itype2(2,2), itype2(3,2), &
                         itype2(1,3), itype2(2,3), itype2(3,3), &
                         itype2(1,4), itype2(2,4), itype2(3,4), &
                         xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                         xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                         xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                         xyz(1,iD), xyz(2,iD),xyz(3,iD))
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        itype2(Imomentum,3) = itype2(Imomentum,3)+2
     ENDIF
     Hessian(iCSTART+Imomentum,iCSTART+Imomentum) = &
          Hessian(iCSTART+Imomentum,iCSTART+Imomentum) &
          +d2II*coeff

     ! JJbas' center.

     itype2(Imomentum,4) = itype2(Imomentum,4)+2
     DO Icon=1,ncontract(Ibas)
        DO Jcon=1,ncontract(Jbas)
           DO IIcon=1,ncontract(IIbas)
              DO JJcon=1,ncontract(JJbas)
                 d2JJ = d2JJ + 4.d0*aexp(JJcon,JJbas)*aexp(JJcon,JJbas) &
                      *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                      *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                      *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                      aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                      itype2(1,1), itype2(2,1), itype2(3,1), &
                      itype2(1,2), itype2(2,2), itype2(3,2), &
                      itype2(1,3), itype2(2,3), itype2(3,3), &
                      itype2(1,4), itype2(2,4), itype2(3,4), &
                      xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                      xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                      xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                      xyz(1,iD), xyz(2,iD),xyz(3,iD))
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     itype2(Imomentum,4) = itype2(Imomentum,4)-2

     DO Icon=1,ncontract(Ibas)
        DO Jcon=1,ncontract(Jbas)
           DO IIcon=1,ncontract(IIbas)
              DO JJcon=1,ncontract(JJbas)
                 d2JJ = d2JJ - 2.d0*aexp(JJcon,JJbas) &
                      *(1.d0+2.d0*dble(itype2(Imomentum,4))) &
                      *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                      *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                      *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                      aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                      itype2(1,1), itype2(2,1), itype2(3,1), &
                      itype2(1,2), itype2(2,2), itype2(3,2), &
                      itype2(1,3), itype2(2,3), itype2(3,3), &
                      itype2(1,4), itype2(2,4), itype2(3,4), &
                      xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                      xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                      xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                      xyz(1,iD), xyz(2,iD),xyz(3,iD))
              ENDDO
           ENDDO
        ENDDO
     ENDDO


     IF (itype2(Imomentum,4) >= 2) THEN
        const = dble(itype2(Imomentum,4)) &
             *dble(itype2(Imomentum,4)-1)
        itype2(Imomentum,4) = itype2(Imomentum,4)-2
        DO Icon=1,ncontract(Ibas)
           DO Jcon=1,ncontract(Jbas)
              DO IIcon=1,ncontract(IIbas)
                 DO JJcon=1,ncontract(JJbas)
                    d2JJ = d2JJ + const* &
                         dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                         *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                         aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                         itype2(1,1), itype2(2,1), itype2(3,1), &
                         itype2(1,2), itype2(2,2), itype2(3,2), &
                         itype2(1,3), itype2(2,3), itype2(3,3), &
                         itype2(1,4), itype2(2,4), itype2(3,4), &
                         xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                         xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                         xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                         xyz(1,iD), xyz(2,iD),xyz(3,iD))
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        itype2(Imomentum,4) = itype2(Imomentum,4)+2
     ENDIF
     Hessian(iDSTART+Imomentum,iDSTART+Imomentum) = &
          Hessian(iDSTART+Imomentum,iDSTART+Imomentum) &
          +d2JJ*coeff


     ! Now do the d^2E/dXAdXB type integrals.

     DO Jmomentum=Imomentum+1,3
        d2I=0.d0
        d2J=0.d0
        d2II=0.d0
        d2JJ=0.d0

        ! Start  with Ibas

        itype2(Imomentum,1) = itype2(Imomentum,1)+1
        itype2(Jmomentum,1) = itype2(Jmomentum,1)+1
        DO Icon=1,ncontract(Ibas)
           DO Jcon=1,ncontract(Jbas)
              DO IIcon=1,ncontract(IIbas)
                 DO JJcon=1,ncontract(JJbas)
                    d2I = d2I + 4.d0*aexp(Icon,Ibas)*aexp(Icon,Ibas) &
                         *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                         *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                         aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                         itype2(1,1), itype2(2,1), itype2(3,1), &
                         itype2(1,2), itype2(2,2), itype2(3,2), &
                         itype2(1,3), itype2(2,3), itype2(3,3), &
                         itype2(1,4), itype2(2,4), itype2(3,4), &
                         xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                         xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                         xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                         xyz(1,iD), xyz(2,iD),xyz(3,iD))
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        itype2(Imomentum,1) = itype2(Imomentum,1)-1
        itype2(Jmomentum,1) = itype2(Jmomentum,1)-1

        IF (itype2(Imomentum,1) /= 0) THEN
           const=dble(itype2(Imomentum,1))
           itype2(Imomentum,1) = itype2(Imomentum,1)-1
           itype2(Jmomentum,1) = itype2(Jmomentum,1)+1
           DO Icon=1,ncontract(Ibas)
              DO Jcon=1,ncontract(Jbas)
                 DO IIcon=1,ncontract(IIbas)
                    DO JJcon=1,ncontract(JJbas)
                       d2I = d2I - 2.d0*aexp(Icon,Ibas)*const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                            *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype2(1,1), itype2(2,1), itype2(3,1), &
                            itype2(1,2), itype2(2,2), itype2(3,2), &
                            itype2(1,3), itype2(2,3), itype2(3,3), &
                            itype2(1,4), itype2(2,4), itype2(3,4), &
                            xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                            xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                            xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                            xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           itype2(Imomentum,1) = itype2(Imomentum,1)+1
           itype2(Jmomentum,1) = itype2(Jmomentum,1)-1
        ENDIF

        IF (itype2(Jmomentum,1) /= 0) THEN
           const=dble(itype2(Jmomentum,1))
           itype2(Imomentum,1) = itype2(Imomentum,1)+1
           itype2(Jmomentum,1) = itype2(Jmomentum,1)-1
           DO Icon=1,ncontract(Ibas)
              DO Jcon=1,ncontract(Jbas)
                 DO IIcon=1,ncontract(IIbas)
                    DO JJcon=1,ncontract(JJbas)
                       d2I = d2I - 2.d0*aexp(Icon,Ibas)*const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                            *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype2(1,1), itype2(2,1), itype2(3,1), &
                            itype2(1,2), itype2(2,2), itype2(3,2), &
                            itype2(1,3), itype2(2,3), itype2(3,3), &
                            itype2(1,4), itype2(2,4), itype2(3,4), &
                            xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                            xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                            xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                            xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           itype2(Imomentum,1) = itype2(Imomentum,1)-1
           itype2(Jmomentum,1) = itype2(Jmomentum,1)+1
        ENDIF

        IF (itype2(Jmomentum,1) /= 0 .AND. itype2(Imomentum,1) /= 0) THEN
           const=dble(itype2(Jmomentum,1)*itype2(Imomentum,1))
           itype2(Imomentum,1) = itype2(Imomentum,1)-1
           itype2(Jmomentum,1) = itype2(Jmomentum,1)-1
           DO Icon=1,ncontract(Ibas)
              DO Jcon=1,ncontract(Jbas)
                 DO IIcon=1,ncontract(IIbas)
                    DO JJcon=1,ncontract(JJbas)
                       d2I = d2I +const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                            *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype2(1,1), itype2(2,1), itype2(3,1), &
                            itype2(1,2), itype2(2,2), itype2(3,2), &
                            itype2(1,3), itype2(2,3), itype2(3,3), &
                            itype2(1,4), itype2(2,4), itype2(3,4), &
                            xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                            xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                            xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                            xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           itype2(Imomentum,1) = itype2(Imomentum,1)+1
           itype2(Jmomentum,1) = itype2(Jmomentum,1)+1
        ENDIF

        Hessian(iASTART+Jmomentum,iASTART+Imomentum) = &
             Hessian(iASTART+Jmomentum,iASTART+Imomentum) &
             +d2I*coeff

        ! Now do Jbas.

        itype2(Imomentum,2) = itype2(Imomentum,2)+1
        itype2(Jmomentum,2) = itype2(Jmomentum,2)+1
        DO Icon=1,ncontract(Ibas)
           DO Jcon=1,ncontract(Jbas)
              DO IIcon=1,ncontract(IIbas)
                 DO JJcon=1,ncontract(JJbas)
                    d2J = d2J + 4.d0*aexp(Jcon,Jbas)*aexp(Jcon,Jbas) &
                         *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                         *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                         aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                         itype2(1,1), itype2(2,1), itype2(3,1), &
                         itype2(1,2), itype2(2,2), itype2(3,2), &
                         itype2(1,3), itype2(2,3), itype2(3,3), &
                         itype2(1,4), itype2(2,4), itype2(3,4), &
                         xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                         xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                         xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                         xyz(1,iD), xyz(2,iD),xyz(3,iD))
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        itype2(Imomentum,2) = itype2(Imomentum,2)-1
        itype2(Jmomentum,2) = itype2(Jmomentum,2)-1

        IF (itype2(Imomentum,2) /= 0) THEN
           const=dble(itype2(Imomentum,2))
           itype2(Imomentum,2) = itype2(Imomentum,2)-1
           itype2(Jmomentum,2) = itype2(Jmomentum,2)+1
           DO Icon=1,ncontract(Ibas)
              DO Jcon=1,ncontract(Jbas)
                 DO IIcon=1,ncontract(IIbas)
                    DO JJcon=1,ncontract(JJbas)
                       d2J = d2J - 2.d0*aexp(Jcon,Jbas)*const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                            *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype2(1,1), itype2(2,1), itype2(3,1), &
                            itype2(1,2), itype2(2,2), itype2(3,2), &
                            itype2(1,3), itype2(2,3), itype2(3,3), &
                            itype2(1,4), itype2(2,4), itype2(3,4), &
                            xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                            xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                            xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                            xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           itype2(Imomentum,2) = itype2(Imomentum,2)+1
           itype2(Jmomentum,2) = itype2(Jmomentum,2)-1
        ENDIF

        IF (itype2(Jmomentum,2) /= 0) THEN
           const=dble(itype2(Jmomentum,2))
           itype2(Imomentum,2) = itype2(Imomentum,2)+1
           itype2(Jmomentum,2) = itype2(Jmomentum,2)-1
           DO Icon=1,ncontract(Ibas)
              DO Jcon=1,ncontract(Jbas)
                 DO IIcon=1,ncontract(IIbas)
                    DO JJcon=1,ncontract(JJbas)
                       d2J = d2J - 2.d0*aexp(Jcon,Jbas)*const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                            *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype2(1,1), itype2(2,1), itype2(3,1), &
                            itype2(1,2), itype2(2,2), itype2(3,2), &
                            itype2(1,3), itype2(2,3), itype2(3,3), &
                            itype2(1,4), itype2(2,4), itype2(3,4), &
                            xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                            xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                            xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                            xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           itype2(Imomentum,2) = itype2(Imomentum,2)-1
           itype2(Jmomentum,2) = itype2(Jmomentum,2)+1
        ENDIF

        IF (itype2(Jmomentum,2) /= 0 .AND. itype2(Imomentum,2) /= 0) THEN
           const=dble(itype2(Jmomentum,2)*itype2(Imomentum,2))
           itype2(Imomentum,2) = itype2(Imomentum,2)-1
           itype2(Jmomentum,2) = itype2(Jmomentum,2)-1
           DO Icon=1,ncontract(Ibas)
              DO Jcon=1,ncontract(Jbas)
                 DO IIcon=1,ncontract(IIbas)
                    DO JJcon=1,ncontract(JJbas)
                       d2J = d2J +const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                            *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype2(1,1), itype2(2,1), itype2(3,1), &
                            itype2(1,2), itype2(2,2), itype2(3,2), &
                            itype2(1,3), itype2(2,3), itype2(3,3), &
                            itype2(1,4), itype2(2,4), itype2(3,4), &
                            xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                            xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                            xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                            xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           itype2(Imomentum,2) = itype2(Imomentum,2)+1
           itype2(Jmomentum,2) = itype2(Jmomentum,2)+1
        ENDIF

        Hessian(iBSTART+Jmomentum,iBSTART+Imomentum) = &
             Hessian(iBSTART+Jmomentum,iBSTART+Imomentum) &
             +d2J*coeff


        ! Do IIbas

        itype2(Imomentum,3) = itype2(Imomentum,3)+1
        itype2(Jmomentum,3) = itype2(Jmomentum,3)+1
        DO Icon=1,ncontract(Ibas)
           DO Jcon=1,ncontract(Jbas)
              DO IIcon=1,ncontract(IIbas)
                 DO JJcon=1,ncontract(JJbas)
                    d2II = d2II + 4.d0*aexp(IIcon,IIbas)*aexp(IIcon,IIbas) &
                         *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                         *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                         aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                         itype2(1,1), itype2(2,1), itype2(3,1), &
                         itype2(1,2), itype2(2,2), itype2(3,2), &
                         itype2(1,3), itype2(2,3), itype2(3,3), &
                         itype2(1,4), itype2(2,4), itype2(3,4), &
                         xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                         xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                         xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                         xyz(1,iD), xyz(2,iD),xyz(3,iD))
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        itype2(Imomentum,3) = itype2(Imomentum,3)-1
        itype2(Jmomentum,3) = itype2(Jmomentum,3)-1

        IF (itype2(Imomentum,3) /= 0) THEN
           const=dble(itype2(Imomentum,3))
           itype2(Imomentum,3) = itype2(Imomentum,3)-1
           itype2(Jmomentum,3) = itype2(Jmomentum,3)+1
           DO Icon=1,ncontract(Ibas)
              DO Jcon=1,ncontract(Jbas)
                 DO IIcon=1,ncontract(IIbas)
                    DO JJcon=1,ncontract(JJbas)
                       d2II = d2II - 2.d0*aexp(IIcon,IIbas)*const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                            *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype2(1,1), itype2(2,1), itype2(3,1), &
                            itype2(1,2), itype2(2,2), itype2(3,2), &
                            itype2(1,3), itype2(2,3), itype2(3,3), &
                            itype2(1,4), itype2(2,4), itype2(3,4), &
                            xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                            xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                            xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                            xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           itype2(Imomentum,3) = itype2(Imomentum,3)+1
           itype2(Jmomentum,3) = itype2(Jmomentum,3)-1
        ENDIF

        IF (itype2(Jmomentum,3) /= 0) THEN
           const=dble(itype2(Jmomentum,3))
           itype2(Imomentum,3) = itype2(Imomentum,3)+1
           itype2(Jmomentum,3) = itype2(Jmomentum,3)-1
           DO Icon=1,ncontract(Ibas)
              DO Jcon=1,ncontract(Jbas)
                 DO IIcon=1,ncontract(IIbas)
                    DO JJcon=1,ncontract(JJbas)
                       d2II = d2II - 2.d0*aexp(IIcon,IIbas)*const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                            *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype2(1,1), itype2(2,1), itype2(3,1), &
                            itype2(1,2), itype2(2,2), itype2(3,2), &
                            itype2(1,3), itype2(2,3), itype2(3,3), &
                            itype2(1,4), itype2(2,4), itype2(3,4), &
                            xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                            xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                            xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                            xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           itype2(Imomentum,3) = itype2(Imomentum,3)-1
           itype2(Jmomentum,3) = itype2(Jmomentum,3)+1
        ENDIF

        IF (itype2(Jmomentum,3) /= 0 .AND. itype2(Imomentum,3) /= 0) THEN
           const=dble(itype2(Jmomentum,3)*itype2(Imomentum,3))
           itype2(Imomentum,3) = itype2(Imomentum,3)-1
           itype2(Jmomentum,3) = itype2(Jmomentum,3)-1
           DO Icon=1,ncontract(Ibas)
              DO Jcon=1,ncontract(Jbas)
                 DO IIcon=1,ncontract(IIbas)
                    DO JJcon=1,ncontract(JJbas)
                       d2II = d2II +const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                            *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype2(1,1), itype2(2,1), itype2(3,1), &
                            itype2(1,2), itype2(2,2), itype2(3,2), &
                            itype2(1,3), itype2(2,3), itype2(3,3), &
                            itype2(1,4), itype2(2,4), itype2(3,4), &
                            xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                            xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                            xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                            xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           itype2(Imomentum,3) = itype2(Imomentum,3)+1
           itype2(Jmomentum,3) = itype2(Jmomentum,3)+1
        ENDIF

        Hessian(iCSTART+Jmomentum,iCSTART+Imomentum) = &
             Hessian(iCSTART+Jmomentum,iCSTART+Imomentum) &
             +d2II*coeff

        ! Now do JJbas.

        itype2(Imomentum,4) = itype2(Imomentum,4)+1
        itype2(Jmomentum,4) = itype2(Jmomentum,4)+1
        DO Icon=1,ncontract(Ibas)
           DO Jcon=1,ncontract(Jbas)
              DO IIcon=1,ncontract(IIbas)
                 DO JJcon=1,ncontract(JJbas)
                    d2JJ = d2JJ + 4.d0*aexp(JJcon,JJbas)*aexp(JJcon,JJbas) &
                         *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                         *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                         aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                         itype2(1,1), itype2(2,1), itype2(3,1), &
                         itype2(1,2), itype2(2,2), itype2(3,2), &
                         itype2(1,3), itype2(2,3), itype2(3,3), &
                         itype2(1,4), itype2(2,4), itype2(3,4), &
                         xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                         xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                         xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                         xyz(1,iD), xyz(2,iD),xyz(3,iD))
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        itype2(Imomentum,4) = itype2(Imomentum,4)-1
        itype2(Jmomentum,4) = itype2(Jmomentum,4)-1

        IF (itype2(Imomentum,4) /= 0) THEN
           const=dble(itype2(Imomentum,4))
           itype2(Imomentum,4) = itype2(Imomentum,4)-1
           itype2(Jmomentum,4) = itype2(Jmomentum,4)+1
           DO Icon=1,ncontract(Ibas)
              DO Jcon=1,ncontract(Jbas)
                 DO IIcon=1,ncontract(IIbas)
                    DO JJcon=1,ncontract(JJbas)
                       d2JJ = d2JJ - 2.d0*aexp(JJcon,JJbas)*const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                            *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype2(1,1), itype2(2,1), itype2(3,1), &
                            itype2(1,2), itype2(2,2), itype2(3,2), &
                            itype2(1,3), itype2(2,3), itype2(3,3), &
                            itype2(1,4), itype2(2,4), itype2(3,4), &
                            xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                            xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                            xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                            xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           itype2(Imomentum,4) = itype2(Imomentum,4)+1
           itype2(Jmomentum,4) = itype2(Jmomentum,4)-1
        ENDIF

        IF (itype2(Jmomentum,4) /= 0) THEN
           const=dble(itype2(Jmomentum,4))
           itype2(Imomentum,4) = itype2(Imomentum,4)+1
           itype2(Jmomentum,4) = itype2(Jmomentum,4)-1
           DO Icon=1,ncontract(Ibas)
              DO Jcon=1,ncontract(Jbas)
                 DO IIcon=1,ncontract(IIbas)
                    DO JJcon=1,ncontract(JJbas)
                       d2JJ= d2JJ - 2.d0*aexp(JJcon,JJbas)*const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                            *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype2(1,1), itype2(2,1), itype2(3,1), &
                            itype2(1,2), itype2(2,2), itype2(3,2), &
                            itype2(1,3), itype2(2,3), itype2(3,3), &
                            itype2(1,4), itype2(2,4), itype2(3,4), &
                            xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                            xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                            xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                            xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           itype2(Imomentum,4) = itype2(Imomentum,4)-1
           itype2(Jmomentum,4) = itype2(Jmomentum,4)+1
        ENDIF

        IF (itype2(Jmomentum,4) /= 0 .AND. itype2(Imomentum,4) /= 0) THEN
           const=dble(itype2(Jmomentum,4)*itype2(Imomentum,4))
           itype2(Imomentum,4) = itype2(Imomentum,4)-1
           itype2(Jmomentum,4) = itype2(Jmomentum,4)-1
           DO Icon=1,ncontract(Ibas)
              DO Jcon=1,ncontract(Jbas)
                 DO IIcon=1,ncontract(IIbas)
                    DO JJcon=1,ncontract(JJbas)
                       d2JJ = d2JJ +const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                            *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype2(1,1), itype2(2,1), itype2(3,1), &
                            itype2(1,2), itype2(2,2), itype2(3,2), &
                            itype2(1,3), itype2(2,3), itype2(3,3), &
                            itype2(1,4), itype2(2,4), itype2(3,4), &
                            xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                            xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                            xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                            xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           itype2(Imomentum,4) = itype2(Imomentum,4)+1
           itype2(Jmomentum,4) = itype2(Jmomentum,4)+1
        ENDIF

        Hessian(iDSTART+Jmomentum,iDSTART+Imomentum) = &
             Hessian(iDSTART+Jmomentum,iDSTART+Imomentum) &
             +d2JJ*coeff

     ENDDO

     ! We have just closed the Jmomentum loop.
     ! Now we need to calculate the d2E/dX1dY2 type terms.  This requires 6
     ! different combinations.

     DO Jmomentum=1,3

        ! Do the Ibas with Jbas first.

        d2XY = 0.d0

        itype2(Imomentum,1) = itype2(Imomentum,1)+1
        itype2(Jmomentum,2) = itype2(Jmomentum,2)+1
        DO Icon=1,ncontract(Ibas)
           DO Jcon=1,ncontract(Jbas)
              DO IIcon=1,ncontract(IIbas)
                 DO JJcon=1,ncontract(JJbas)
                    d2XY = d2XY + 4.d0*aexp(Icon,Ibas)*aexp(Jcon,Jbas) &
                         *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                         *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                         aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                         itype2(1,1), itype2(2,1), itype2(3,1), &
                         itype2(1,2), itype2(2,2), itype2(3,2), &
                         itype2(1,3), itype2(2,3), itype2(3,3), &
                         itype2(1,4), itype2(2,4), itype2(3,4), &
                         xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                         xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                         xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                         xyz(1,iD), xyz(2,iD),xyz(3,iD))
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        itype2(Imomentum,1) = itype2(Imomentum,1)-1
        itype2(Jmomentum,2) = itype2(Jmomentum,2)-1

        IF (itype2(Imomentum,1) /= 0) THEN
           const=dble(itype2(Imomentum,1))
           itype2(Imomentum,1) = itype2(Imomentum,1)-1
           itype2(Jmomentum,2) = itype2(Jmomentum,2)+1
           DO Icon=1,ncontract(Ibas)
              DO Jcon=1,ncontract(Jbas)
                 DO IIcon=1,ncontract(IIbas)
                    DO JJcon=1,ncontract(JJbas)
                       d2XY = d2XY - 2.d0*aexp(Jcon,Jbas)*const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                            *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype2(1,1), itype2(2,1), itype2(3,1), &
                            itype2(1,2), itype2(2,2), itype2(3,2), &
                            itype2(1,3), itype2(2,3), itype2(3,3), &
                            itype2(1,4), itype2(2,4), itype2(3,4), &
                            xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                            xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                            xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                            xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           itype2(Imomentum,1) = itype2(Imomentum,1)+1
           itype2(Jmomentum,2) = itype2(Jmomentum,2)-1
        ENDIF

        IF (itype2(Jmomentum,2) /= 0) THEN
           const=dble(itype2(Jmomentum,2))
           itype2(Imomentum,1) = itype2(Imomentum,1)+1
           itype2(Jmomentum,2) = itype2(Jmomentum,2)-1
           DO Icon=1,ncontract(Ibas)
              DO Jcon=1,ncontract(Jbas)
                 DO IIcon=1,ncontract(IIbas)
                    DO JJcon=1,ncontract(JJbas)
                       d2XY = d2XY - 2.d0*aexp(Icon,Ibas)*const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                            *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype2(1,1), itype2(2,1), itype2(3,1), &
                            itype2(1,2), itype2(2,2), itype2(3,2), &
                            itype2(1,3), itype2(2,3), itype2(3,3), &
                            itype2(1,4), itype2(2,4), itype2(3,4), &
                            xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                            xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                            xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                            xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           itype2(Imomentum,1) = itype2(Imomentum,1)-1
           itype2(Jmomentum,2) = itype2(Jmomentum,2)+1
        ENDIF

        IF (itype2(Jmomentum,2) /= 0 .AND. itype2(Imomentum,1) /= 0) THEN
           const=dble(itype2(Jmomentum,2)*itype2(Imomentum,1))
           itype2(Imomentum,1) = itype2(Imomentum,1)-1
           itype2(Jmomentum,2) = itype2(Jmomentum,2)-1
           DO Icon=1,ncontract(Ibas)
              DO Jcon=1,ncontract(Jbas)
                 DO IIcon=1,ncontract(IIbas)
                    DO JJcon=1,ncontract(JJbas)
                       d2XY = d2XY +const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                            *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype2(1,1), itype2(2,1), itype2(3,1), &
                            itype2(1,2), itype2(2,2), itype2(3,2), &
                            itype2(1,3), itype2(2,3), itype2(3,3), &
                            itype2(1,4), itype2(2,4), itype2(3,4), &
                            xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                            xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                            xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                            xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           itype2(Imomentum,1) = itype2(Imomentum,1)+1
           itype2(Jmomentum,2) = itype2(Jmomentum,2)+1
        ENDIF
        IF (iA /= iB) THEN
           IF (iB > iA) THEN
              Hessian(iBSTART+Jmomentum,iASTART+Imomentum) = &
                   Hessian(iBSTART+Jmomentum,iASTART+Imomentum) &
                   +d2XY*coeff
           ELSE
              Hessian(iASTART+Imomentum,iBSTART+Jmomentum) = &
                   Hessian(iASTART+Imomentum,iBSTART+Jmomentum) &
                   +d2XY*coeff
           ENDIF
        ELSE
           IF (Imomentum == Jmomentum) THEN
              Hessian(iBSTART+Jmomentum,iASTART+Imomentum) = &
                   Hessian(iBSTART+Jmomentum,iASTART+Imomentum) &
                   +2.d0*d2XY*coeff
           ELSEIF (Jmomentum > Imomentum) THEN
              Hessian(iBSTART+Jmomentum,iASTART+Imomentum) = &
                   Hessian(iBSTART+Jmomentum,iASTART+Imomentum) &
                   +d2XY*coeff
           ELSE
              Hessian(iASTART+Imomentum,iBSTART+Jmomentum) = &
                   Hessian(iASTART+Imomentum,iBSTART+Jmomentum) &
                   +d2XY*coeff
           ENDIF
        ENDIF

        ! Ibas with IIbas.

        d2XY = 0.d0

        itype2(Imomentum,1) = itype2(Imomentum,1)+1
        itype2(Jmomentum,3) = itype2(Jmomentum,3)+1
        DO Icon=1,ncontract(Ibas)
           DO Jcon=1,ncontract(Jbas)
              DO IIcon=1,ncontract(IIbas)
                 DO JJcon=1,ncontract(JJbas)
                    d2XY = d2XY + 4.d0*aexp(Icon,Ibas)*aexp(IIcon,IIbas) &
                         *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                         *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                         aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                         itype2(1,1), itype2(2,1), itype2(3,1), &
                         itype2(1,2), itype2(2,2), itype2(3,2), &
                         itype2(1,3), itype2(2,3), itype2(3,3), &
                         itype2(1,4), itype2(2,4), itype2(3,4), &
                         xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                         xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                         xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                         xyz(1,iD), xyz(2,iD),xyz(3,iD))
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        itype2(Imomentum,1) = itype2(Imomentum,1)-1
        itype2(Jmomentum,3) = itype2(Jmomentum,3)-1

        IF (itype2(Imomentum,1) /= 0) THEN
           const=dble(itype2(Imomentum,1))
           itype2(Imomentum,1) = itype2(Imomentum,1)-1
           itype2(Jmomentum,3) = itype2(Jmomentum,3)+1
           DO Icon=1,ncontract(Ibas)
              DO Jcon=1,ncontract(Jbas)
                 DO IIcon=1,ncontract(IIbas)
                    DO JJcon=1,ncontract(JJbas)
                       d2XY = d2XY - 2.d0*aexp(IIcon,IIbas)*const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                            *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype2(1,1), itype2(2,1), itype2(3,1), &
                            itype2(1,2), itype2(2,2), itype2(3,2), &
                            itype2(1,3), itype2(2,3), itype2(3,3), &
                            itype2(1,4), itype2(2,4), itype2(3,4), &
                            xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                            xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                            xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                            xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           itype2(Imomentum,1) = itype2(Imomentum,1)+1
           itype2(Jmomentum,3) = itype2(Jmomentum,3)-1
        ENDIF

        IF (itype2(Jmomentum,3) /= 0) THEN
           const=dble(itype2(Jmomentum,3))
           itype2(Imomentum,1) = itype2(Imomentum,1)+1
           itype2(Jmomentum,3) = itype2(Jmomentum,3)-1
           DO Icon=1,ncontract(Ibas)
              DO Jcon=1,ncontract(Jbas)
                 DO IIcon=1,ncontract(IIbas)
                    DO JJcon=1,ncontract(JJbas)
                       d2XY = d2XY - 2.d0*aexp(Icon,Ibas)*const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                            *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype2(1,1), itype2(2,1), itype2(3,1), &
                            itype2(1,2), itype2(2,2), itype2(3,2), &
                            itype2(1,3), itype2(2,3), itype2(3,3), &
                            itype2(1,4), itype2(2,4), itype2(3,4), &
                            xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                            xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                            xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                            xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           itype2(Imomentum,1) = itype2(Imomentum,1)-1
           itype2(Jmomentum,3) = itype2(Jmomentum,3)+1
        ENDIF

        IF (itype2(Jmomentum,3) /= 0 .AND. itype2(Imomentum,1) /= 0) THEN
           const=dble(itype2(Jmomentum,3)*itype2(Imomentum,1))
           itype2(Imomentum,1) = itype2(Imomentum,1)-1
           itype2(Jmomentum,3) = itype2(Jmomentum,3)-1
           DO Icon=1,ncontract(Ibas)
              DO Jcon=1,ncontract(Jbas)
                 DO IIcon=1,ncontract(IIbas)
                    DO JJcon=1,ncontract(JJbas)
                       d2XY = d2XY +const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                            *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype2(1,1), itype2(2,1), itype2(3,1), &
                            itype2(1,2), itype2(2,2), itype2(3,2), &
                            itype2(1,3), itype2(2,3), itype2(3,3), &
                            itype2(1,4), itype2(2,4), itype2(3,4), &
                            xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                            xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                            xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                            xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           itype2(Imomentum,1) = itype2(Imomentum,1)+1
           itype2(Jmomentum,3) = itype2(Jmomentum,3)+1
        ENDIF
        IF (iA /= iC) THEN
           IF (iC > iA) THEN
              Hessian(iCSTART+Jmomentum,iASTART+Imomentum) = &
                   Hessian(iCSTART+Jmomentum,iASTART+Imomentum) &
                   +d2XY*coeff
           ELSE
              Hessian(iASTART+Imomentum,iCSTART+Jmomentum) = &
                   Hessian(iASTART+Imomentum,iCSTART+Jmomentum) &
                   +d2XY*coeff
           ENDIF
        ELSE
           IF (Imomentum == Jmomentum) THEN
              Hessian(iCSTART+Jmomentum,iASTART+Imomentum) = &
                   Hessian(iCSTART+Jmomentum,iASTART+Imomentum) &
                   +2.d0*d2XY*coeff
           ELSEIF (Jmomentum > Imomentum) THEN
              Hessian(iCSTART+Jmomentum,iASTART+Imomentum) = &
                   Hessian(iCSTART+Jmomentum,iASTART+Imomentum) &
                   +d2XY*coeff
           ELSE
              Hessian(iASTART+Imomentum,iCSTART+Jmomentum) = &
                   Hessian(iASTART+Imomentum,iCSTART+Jmomentum) &
                   +d2XY*coeff
           ENDIF
        ENDIF

        ! Ibas with JJbas.

        d2XY = 0.d0

        itype2(Imomentum,1) = itype2(Imomentum,1)+1
        itype2(Jmomentum,4) = itype2(Jmomentum,4)+1
        DO Icon=1,ncontract(Ibas)
           DO Jcon=1,ncontract(Jbas)
              DO IIcon=1,ncontract(IIbas)
                 DO JJcon=1,ncontract(JJbas)
                    d2XY = d2XY + 4.d0*aexp(Icon,Ibas)*aexp(JJcon,JJbas) &
                         *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                         *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                         aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                         itype2(1,1), itype2(2,1), itype2(3,1), &
                         itype2(1,2), itype2(2,2), itype2(3,2), &
                         itype2(1,3), itype2(2,3), itype2(3,3), &
                         itype2(1,4), itype2(2,4), itype2(3,4), &
                         xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                         xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                         xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                         xyz(1,iD), xyz(2,iD),xyz(3,iD))
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        itype2(Imomentum,1) = itype2(Imomentum,1)-1
        itype2(Jmomentum,4) = itype2(Jmomentum,4)-1

        IF (itype2(Imomentum,1) /= 0) THEN
           const=dble(itype2(Imomentum,1))
           itype2(Imomentum,1) = itype2(Imomentum,1)-1
           itype2(Jmomentum,4) = itype2(Jmomentum,4)+1
           DO Icon=1,ncontract(Ibas)
              DO Jcon=1,ncontract(Jbas)
                 DO IIcon=1,ncontract(IIbas)
                    DO JJcon=1,ncontract(JJbas)
                       d2XY = d2XY - 2.d0*aexp(JJcon,JJbas)*const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                            *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype2(1,1), itype2(2,1), itype2(3,1), &
                            itype2(1,2), itype2(2,2), itype2(3,2), &
                            itype2(1,3), itype2(2,3), itype2(3,3), &
                            itype2(1,4), itype2(2,4), itype2(3,4), &
                            xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                            xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                            xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                            xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           itype2(Imomentum,1) = itype2(Imomentum,1)+1
           itype2(Jmomentum,4) = itype2(Jmomentum,4)-1
        ENDIF

        IF (itype2(Jmomentum,4) /= 0) THEN
           const=dble(itype2(Jmomentum,4))
           itype2(Imomentum,1) = itype2(Imomentum,1)+1
           itype2(Jmomentum,4) = itype2(Jmomentum,4)-1
           DO Icon=1,ncontract(Ibas)
              DO Jcon=1,ncontract(Jbas)
                 DO IIcon=1,ncontract(IIbas)
                    DO JJcon=1,ncontract(JJbas)
                       d2XY = d2XY - 2.d0*aexp(Icon,Ibas)*const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                            *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype2(1,1), itype2(2,1), itype2(3,1), &
                            itype2(1,2), itype2(2,2), itype2(3,2), &
                            itype2(1,3), itype2(2,3), itype2(3,3), &
                            itype2(1,4), itype2(2,4), itype2(3,4), &
                            xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                            xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                            xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                            xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           itype2(Imomentum,1) = itype2(Imomentum,1)-1
           itype2(Jmomentum,4) = itype2(Jmomentum,4)+1
        ENDIF

        IF (itype2(Jmomentum,4) /= 0 .AND. itype2(Imomentum,1) /= 0) THEN
           const=dble(itype2(Jmomentum,4)*itype2(Imomentum,1))
           itype2(Imomentum,1) = itype2(Imomentum,1)-1
           itype2(Jmomentum,4) = itype2(Jmomentum,4)-1
           DO Icon=1,ncontract(Ibas)
              DO Jcon=1,ncontract(Jbas)
                 DO IIcon=1,ncontract(IIbas)
                    DO JJcon=1,ncontract(JJbas)
                       d2XY = d2XY +const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                            *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype2(1,1), itype2(2,1), itype2(3,1), &
                            itype2(1,2), itype2(2,2), itype2(3,2), &
                            itype2(1,3), itype2(2,3), itype2(3,3), &
                            itype2(1,4), itype2(2,4), itype2(3,4), &
                            xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                            xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                            xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                            xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           itype2(Imomentum,1) = itype2(Imomentum,1)+1
           itype2(Jmomentum,4) = itype2(Jmomentum,4)+1
        ENDIF
        IF (iA /= iD) THEN
           IF (iD > iA) THEN
              Hessian(iDSTART+Jmomentum,iASTART+Imomentum) = &
                   Hessian(iDSTART+Jmomentum,iASTART+Imomentum) &
                   +d2XY*coeff
           ELSE
              Hessian(iASTART+Imomentum,iDSTART+Jmomentum) = &
                   Hessian(iASTART+Imomentum,iDSTART+Jmomentum) &
                   +d2XY*coeff
           ENDIF
        ELSE
           IF (Imomentum == Jmomentum) THEN
              Hessian(iDSTART+Jmomentum,iASTART+Imomentum) = &
                   Hessian(iDSTART+Jmomentum,iASTART+Imomentum) &
                   +2.d0*d2XY*coeff
           ELSEIF (Jmomentum > Imomentum) THEN
              Hessian(iDSTART+Jmomentum,iASTART+Imomentum) = &
                   Hessian(iDSTART+Jmomentum,iASTART+Imomentum) &
                   +d2XY*coeff
           ELSE
              Hessian(iASTART+Imomentum,iDSTART+Jmomentum) = &
                   Hessian(iASTART+Imomentum,iDSTART+Jmomentum) &
                   +d2XY*coeff
           ENDIF
        ENDIF

        ! IIbas with Jbas.

        d2XY = 0.d0

        itype2(Imomentum,3) = itype2(Imomentum,3)+1
        itype2(Jmomentum,2) = itype2(Jmomentum,2)+1
        DO Icon=1,ncontract(Ibas)
           DO Jcon=1,ncontract(Jbas)
              DO IIcon=1,ncontract(IIbas)
                 DO JJcon=1,ncontract(JJbas)
                    d2XY = d2XY + 4.d0*aexp(IIcon,IIbas)*aexp(Jcon,Jbas) &
                         *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                         *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                         aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                         itype2(1,1), itype2(2,1), itype2(3,1), &
                         itype2(1,2), itype2(2,2), itype2(3,2), &
                         itype2(1,3), itype2(2,3), itype2(3,3), &
                         itype2(1,4), itype2(2,4), itype2(3,4), &
                         xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                         xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                         xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                         xyz(1,iD), xyz(2,iD),xyz(3,iD))
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        itype2(Imomentum,3) = itype2(Imomentum,3)-1
        itype2(Jmomentum,2) = itype2(Jmomentum,2)-1

        IF (itype2(Imomentum,3) /= 0) THEN
           const=dble(itype2(Imomentum,3))
           itype2(Imomentum,3) = itype2(Imomentum,3)-1
           itype2(Jmomentum,2) = itype2(Jmomentum,2)+1
           DO Icon=1,ncontract(Ibas)
              DO Jcon=1,ncontract(Jbas)
                 DO IIcon=1,ncontract(IIbas)
                    DO JJcon=1,ncontract(JJbas)
                       d2XY = d2XY - 2.d0*aexp(Jcon,Jbas)*const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                            *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype2(1,1), itype2(2,1), itype2(3,1), &
                            itype2(1,2), itype2(2,2), itype2(3,2), &
                            itype2(1,3), itype2(2,3), itype2(3,3), &
                            itype2(1,4), itype2(2,4), itype2(3,4), &
                            xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                            xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                            xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                            xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           itype2(Imomentum,3) = itype2(Imomentum,3)+1
           itype2(Jmomentum,2) = itype2(Jmomentum,2)-1
        ENDIF

        IF (itype2(Jmomentum,2) /= 0) THEN
           const=dble(itype2(Jmomentum,2))
           itype2(Imomentum,3) = itype2(Imomentum,3)+1
           itype2(Jmomentum,2) = itype2(Jmomentum,2)-1
           DO Icon=1,ncontract(Ibas)
              DO Jcon=1,ncontract(Jbas)
                 DO IIcon=1,ncontract(IIbas)
                    DO JJcon=1,ncontract(JJbas)
                       d2XY = d2XY - 2.d0*aexp(IIcon,IIbas)*const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                            *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype2(1,1), itype2(2,1), itype2(3,1), &
                            itype2(1,2), itype2(2,2), itype2(3,2), &
                            itype2(1,3), itype2(2,3), itype2(3,3), &
                            itype2(1,4), itype2(2,4), itype2(3,4), &
                            xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                            xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                            xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                            xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           itype2(Imomentum,3) = itype2(Imomentum,3)-1
           itype2(Jmomentum,2) = itype2(Jmomentum,2)+1
        ENDIF

        IF (itype2(Jmomentum,2) /= 0 .AND. itype2(Imomentum,3) /= 0) THEN
           const=dble(itype2(Jmomentum,2)*itype2(Imomentum,3))
           itype2(Imomentum,3) = itype2(Imomentum,3)-1
           itype2(Jmomentum,2) = itype2(Jmomentum,2)-1
           DO Icon=1,ncontract(Ibas)
              DO Jcon=1,ncontract(Jbas)
                 DO IIcon=1,ncontract(IIbas)
                    DO JJcon=1,ncontract(JJbas)
                       d2XY = d2XY +const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                            *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype2(1,1), itype2(2,1), itype2(3,1), &
                            itype2(1,2), itype2(2,2), itype2(3,2), &
                            itype2(1,3), itype2(2,3), itype2(3,3), &
                            itype2(1,4), itype2(2,4), itype2(3,4), &
                            xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                            xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                            xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                            xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           itype2(Imomentum,3) = itype2(Imomentum,3)+1
           itype2(Jmomentum,2) = itype2(Jmomentum,2)+1
        ENDIF
        IF (iC /= iB) THEN
           IF (iB > iC) THEN
              Hessian(iBSTART+Jmomentum,iCSTART+Imomentum) = &
                   Hessian(iBSTART+Jmomentum,iCSTART+Imomentum) &
                   +d2XY*coeff
           ELSE
              Hessian(iCSTART+Imomentum,iBSTART+Jmomentum) = &
                   Hessian(iCSTART+Imomentum,iBSTART+Jmomentum) &
                   +d2XY*coeff
           ENDIF
        ELSE
           IF (Imomentum == Jmomentum) THEN
              Hessian(iBSTART+Jmomentum,iCSTART+Imomentum) = &
                   Hessian(iBSTART+Jmomentum,iCSTART+Imomentum) &
                   +2.d0*d2XY*coeff
           ELSEIF (Jmomentum > Imomentum) THEN
              Hessian(iBSTART+Jmomentum,iCSTART+Imomentum) = &
                   Hessian(iBSTART+Jmomentum,iCSTART+Imomentum) &
                   +d2XY*coeff
           ELSE
              Hessian(iCSTART+Imomentum,iBSTART+Jmomentum) = &
                   Hessian(iCSTART+Imomentum,iBSTART+Jmomentum) &
                   +d2XY*coeff
           ENDIF
        ENDIF

        ! JJbas with Jbas.

        d2XY = 0.d0

        itype2(Imomentum,4) = itype2(Imomentum,4)+1
        itype2(Jmomentum,2) = itype2(Jmomentum,2)+1
        DO Icon=1,ncontract(Ibas)
           DO Jcon=1,ncontract(Jbas)
              DO IIcon=1,ncontract(IIbas)
                 DO JJcon=1,ncontract(JJbas)
                    d2XY = d2XY + 4.d0*aexp(JJcon,JJbas)*aexp(Jcon,Jbas) &
                         *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                         *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                         aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                         itype2(1,1), itype2(2,1), itype2(3,1), &
                         itype2(1,2), itype2(2,2), itype2(3,2), &
                         itype2(1,3), itype2(2,3), itype2(3,3), &
                         itype2(1,4), itype2(2,4), itype2(3,4), &
                         xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                         xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                         xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                         xyz(1,iD), xyz(2,iD),xyz(3,iD))
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        itype2(Imomentum,4) = itype2(Imomentum,4)-1
        itype2(Jmomentum,2) = itype2(Jmomentum,2)-1

        IF (itype2(Imomentum,4) /= 0) THEN
           const=dble(itype2(Imomentum,4))
           itype2(Imomentum,4) = itype2(Imomentum,4)-1
           itype2(Jmomentum,2) = itype2(Jmomentum,2)+1
           DO Icon=1,ncontract(Ibas)
              DO Jcon=1,ncontract(Jbas)
                 DO IIcon=1,ncontract(IIbas)
                    DO JJcon=1,ncontract(JJbas)
                       d2XY = d2XY - 2.d0*aexp(Jcon,Jbas)*const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                            *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype2(1,1), itype2(2,1), itype2(3,1), &
                            itype2(1,2), itype2(2,2), itype2(3,2), &
                            itype2(1,3), itype2(2,3), itype2(3,3), &
                            itype2(1,4), itype2(2,4), itype2(3,4), &
                            xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                            xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                            xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                            xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           itype2(Imomentum,4) = itype2(Imomentum,4)+1
           itype2(Jmomentum,2) = itype2(Jmomentum,2)-1
        ENDIF

        IF (itype2(Jmomentum,2) /= 0) THEN
           const=dble(itype2(Jmomentum,2))
           itype2(Imomentum,4) = itype2(Imomentum,4)+1
           itype2(Jmomentum,2) = itype2(Jmomentum,2)-1
           DO Icon=1,ncontract(Ibas)
              DO Jcon=1,ncontract(Jbas)
                 DO IIcon=1,ncontract(IIbas)
                    DO JJcon=1,ncontract(JJbas)
                       d2XY = d2XY - 2.d0*aexp(JJcon,JJbas)*const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                            *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype2(1,1), itype2(2,1), itype2(3,1), &
                            itype2(1,2), itype2(2,2), itype2(3,2), &
                            itype2(1,3), itype2(2,3), itype2(3,3), &
                            itype2(1,4), itype2(2,4), itype2(3,4), &
                            xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                            xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                            xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                            xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           itype2(Imomentum,4) = itype2(Imomentum,4)-1
           itype2(Jmomentum,2) = itype2(Jmomentum,2)+1
        ENDIF

        IF (itype2(Jmomentum,2) /= 0 .AND. itype2(Imomentum,4) /= 0) THEN
           const=dble(itype2(Jmomentum,2)*itype2(Imomentum,4))
           itype2(Imomentum,4) = itype2(Imomentum,4)-1
           itype2(Jmomentum,2) = itype2(Jmomentum,2)-1
           DO Icon=1,ncontract(Ibas)
              DO Jcon=1,ncontract(Jbas)
                 DO IIcon=1,ncontract(IIbas)
                    DO JJcon=1,ncontract(JJbas)
                       d2XY = d2XY +const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                            *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype2(1,1), itype2(2,1), itype2(3,1), &
                            itype2(1,2), itype2(2,2), itype2(3,2), &
                            itype2(1,3), itype2(2,3), itype2(3,3), &
                            itype2(1,4), itype2(2,4), itype2(3,4), &
                            xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                            xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                            xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                            xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           itype2(Imomentum,4) = itype2(Imomentum,4)+1
           itype2(Jmomentum,2) = itype2(Jmomentum,2)+1
        ENDIF
        IF (iD /= iB) THEN
           IF (iB > iD) THEN
              Hessian(iBSTART+Jmomentum,iDSTART+Imomentum) = &
                   Hessian(iBSTART+Jmomentum,iDSTART+Imomentum) &
                   +d2XY*coeff
           ELSE
              Hessian(iDSTART+Imomentum,iBSTART+Jmomentum) = &
                   Hessian(iDSTART+Imomentum,iBSTART+Jmomentum) &
                   +d2XY*coeff
           ENDIF
        ELSE
           IF (Imomentum == Jmomentum) THEN
              Hessian(iBSTART+Jmomentum,iDSTART+Imomentum) = &
                   Hessian(iBSTART+Jmomentum,iDSTART+Imomentum) &
                   +2.d0*d2XY*coeff
           ELSEIF (Jmomentum > Imomentum) THEN
              Hessian(iBSTART+Jmomentum,iDSTART+Imomentum) = &
                   Hessian(iBSTART+Jmomentum,iDSTART+Imomentum) &
                   +d2XY*coeff
           ELSE
              Hessian(iDSTART+Imomentum,iBSTART+Jmomentum) = &
                   Hessian(iDSTART+Imomentum,iBSTART+Jmomentum) &
                   +d2XY*coeff
           ENDIF
        ENDIF

        ! JJbas with IIbas.

        d2XY = 0.d0

        itype2(Imomentum,4) = itype2(Imomentum,4)+1
        itype2(Jmomentum,3) = itype2(Jmomentum,3)+1
        DO Icon=1,ncontract(Ibas)
           DO Jcon=1,ncontract(Jbas)
              DO IIcon=1,ncontract(IIbas)
                 DO JJcon=1,ncontract(JJbas)
                    d2XY = d2XY + 4.d0*aexp(JJcon,JJbas)*aexp(IIcon,IIbas) &
                         *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                         *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                         aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                         itype2(1,1), itype2(2,1), itype2(3,1), &
                         itype2(1,2), itype2(2,2), itype2(3,2), &
                         itype2(1,3), itype2(2,3), itype2(3,3), &
                         itype2(1,4), itype2(2,4), itype2(3,4), &
                         xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                         xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                         xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                         xyz(1,iD), xyz(2,iD),xyz(3,iD))
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        itype2(Imomentum,4) = itype2(Imomentum,4)-1
        itype2(Jmomentum,3) = itype2(Jmomentum,3)-1

        IF (itype2(Imomentum,4) /= 0) THEN
           const=dble(itype2(Imomentum,4))
           itype2(Imomentum,4) = itype2(Imomentum,4)-1
           itype2(Jmomentum,3) = itype2(Jmomentum,3)+1
           DO Icon=1,ncontract(Ibas)
              DO Jcon=1,ncontract(Jbas)
                 DO IIcon=1,ncontract(IIbas)
                    DO JJcon=1,ncontract(JJbas)
                       d2XY = d2XY - 2.d0*aexp(IIcon,IIbas)*const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                            *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype2(1,1), itype2(2,1), itype2(3,1), &
                            itype2(1,2), itype2(2,2), itype2(3,2), &
                            itype2(1,3), itype2(2,3), itype2(3,3), &
                            itype2(1,4), itype2(2,4), itype2(3,4), &
                            xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                            xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                            xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                            xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           itype2(Imomentum,4) = itype2(Imomentum,4)+1
           itype2(Jmomentum,3) = itype2(Jmomentum,3)-1
        ENDIF

        IF (itype2(Jmomentum,3) /= 0) THEN
           const=dble(itype2(Jmomentum,3))
           itype2(Imomentum,4) = itype2(Imomentum,4)+1
           itype2(Jmomentum,3) = itype2(Jmomentum,3)-1
           DO Icon=1,ncontract(Ibas)
              DO Jcon=1,ncontract(Jbas)
                 DO IIcon=1,ncontract(IIbas)
                    DO JJcon=1,ncontract(JJbas)
                       d2XY = d2XY - 2.d0*aexp(JJcon,JJbas)*const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                            *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype2(1,1), itype2(2,1), itype2(3,1), &
                            itype2(1,2), itype2(2,2), itype2(3,2), &
                            itype2(1,3), itype2(2,3), itype2(3,3), &
                            itype2(1,4), itype2(2,4), itype2(3,4), &
                            xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                            xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                            xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                            xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           itype2(Imomentum,4) = itype2(Imomentum,4)-1
           itype2(Jmomentum,3) = itype2(Jmomentum,3)+1
        ENDIF

        IF (itype2(Jmomentum,3) /= 0 .AND. itype2(Imomentum,4) /= 0) THEN
           const=dble(itype2(Jmomentum,3)*itype2(Imomentum,4))
           itype2(Imomentum,4) = itype2(Imomentum,4)-1
           itype2(Jmomentum,3) = itype2(Jmomentum,3)-1
           DO Icon=1,ncontract(Ibas)
              DO Jcon=1,ncontract(Jbas)
                 DO IIcon=1,ncontract(IIbas)
                    DO JJcon=1,ncontract(JJbas)
                       d2XY = d2XY +const &
                            *dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                            *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype2(1,1), itype2(2,1), itype2(3,1), &
                            itype2(1,2), itype2(2,2), itype2(3,2), &
                            itype2(1,3), itype2(2,3), itype2(3,3), &
                            itype2(1,4), itype2(2,4), itype2(3,4), &
                            xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                            xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                            xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                            xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           itype2(Imomentum,4) = itype2(Imomentum,4)+1
           itype2(Jmomentum,3) = itype2(Jmomentum,3)+1
        ENDIF
        IF (iD /= iC) THEN
           IF (iC > iD) THEN
              Hessian(iCSTART+Jmomentum,iDSTART+Imomentum) = &
                   Hessian(iCSTART+Jmomentum,iDSTART+Imomentum) &
                   +d2XY*coeff
           ELSE
              Hessian(iDSTART+Imomentum,iCSTART+Jmomentum) = &
                   Hessian(iDSTART+Imomentum,iCSTART+Jmomentum) &
                   +d2XY*coeff
           ENDIF
        ELSE
           IF (Imomentum == Jmomentum) THEN
              Hessian(iCSTART+Jmomentum,iDSTART+Imomentum) = &
                   Hessian(iCSTART+Jmomentum,iDSTART+Imomentum) &
                   +2.d0*d2XY*coeff
           ELSEIF (Jmomentum > Imomentum) THEN
              Hessian(iCSTART+Jmomentum,iDSTART+Imomentum) = &
                   Hessian(iCSTART+Jmomentum,iDSTART+Imomentum) &
                   +d2XY*coeff
           ELSE
              Hessian(iDSTART+Imomentum,iCSTART+Jmomentum) = &
                   Hessian(iDSTART+Imomentum,iCSTART+Jmomentum) &
                   +d2XY*coeff
           ENDIF
        ENDIF

     ENDDO
  ENDDO
  return
end subroutine hess2elec
! Ed Brothers. December 12, 2002.
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP



subroutine hfdmxderuse(IDX)
  use allmod
  implicit double precision(a-h,o-z)
  dimension GRADIENT2(natom*3)

  ! When calling this code, the derivative of the alpha and beta density
  ! matrices with respect to IDX are in HOLD and HOLD2.  We are going to
  ! form the contribution of the hessian that comes from the the first
  ! derivative of the density matrix combined with the first derivates
  ! of the integrals.  This is stored in Gradient2, and is then added
  ! to the hessian.

  ! Note that this is "uhfgradient.F" modified.  Thats why the temporary
  ! array is called gradient2.  The code is a little hard to follow
  ! for the 4center 2e- integrals, but if you work it out by hand and
  ! then compare to the integrals, it will become reasonably transparent.
  ! Note that for the two and three encer integrals, DENSEJI etc contains
  ! the density derivatives.

  DO Iatm=1,natom*3
     GRADIENT2(iatm)=0.d0
  ENDDO


  ! 1)  The derivative of the 1 electron kinetic energy term ij times
  ! the density matrix derivative element ij.

  DO Ibas=1,nbasis
     ISTART = (ncenter(Ibas)-1) *3
     DO Jbas=ilast(ncenter(IBAS))+1,nbasis
        JSTART = (ncenter(Jbas)-1) *3
        DENSEJI = HOLD(Jbas,Ibas)+HOLD2(Jbas,Ibas)

        ! We have selected our two basis functions, now loop over angular momentum.

        DO Imomentum=1,3
           dKEI = 0.d0
           dKEJ = 0.d0

           ! Do the Ibas derivatives first.

           itype(Imomentum,Ibas) = itype(Imomentum,Ibas)+1
           DO Icon=1,ncontract(Ibas)
              DO Jcon=1,ncontract(Jbas)
                 dKEI = dKEI + 2.d0*aexp(Icon,Ibas)* &
                      dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                      *ekinetic(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                      itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                      itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                      xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                      xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                      xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
              ENDDO
           ENDDO
           itype(Imomentum,Ibas) = itype(Imomentum,Ibas)-1
           IF (itype(Imomentum,Ibas) /= 0) THEN
              itype(Imomentum,Ibas) = itype(Imomentum,Ibas)-1
              DO Icon=1,ncontract(Ibas)
                 DO Jcon=1,ncontract(Jbas)
                    dKEI = dKEI - dble(itype(Imomentum,Ibas)+1)* &
                         dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *ekinetic(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                         itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                         itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                         xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                         xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                         xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                 ENDDO
              ENDDO
              itype(Imomentum,Ibas) = itype(Imomentum,Ibas)+1
           ENDIF
           GRADIENT2(ISTART+Imomentum) = GRADIENT2(ISTART+Imomentum) &
                +dKeI*DENSEJI*2.d0

           ! Now do the Jbas derivatives.

           itype(Imomentum,Jbas) = itype(Imomentum,Jbas)+1
           DO Icon=1,ncontract(Ibas)
              DO Jcon=1,ncontract(Jbas)
                 dKEJ = dKEJ + 2.d0*aexp(Jcon,Jbas)* &
                      dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                      *ekinetic(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                      itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                      itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                      xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                      xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                      xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
              ENDDO
           ENDDO
           itype(Imomentum,Jbas) = itype(Imomentum,Jbas)-1
           IF (itype(Imomentum,Jbas) /= 0) THEN
              itype(Imomentum,Jbas) = itype(Imomentum,Jbas)-1
              DO Icon=1,ncontract(Ibas)
                 DO Jcon=1,ncontract(Jbas)
                    dKEJ = dKEJ - dble(itype(Imomentum,Jbas)+1)* &
                         dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *ekinetic(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                         itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                         itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                         xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                         xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                         xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                 ENDDO
              ENDDO
              itype(Imomentum,Jbas) = itype(Imomentum,Jbas)+1
           ENDIF
           GRADIENT2(JSTART+Imomentum) = GRADIENT2(JSTART+Imomentum) &
                +dKEJ*DENSEJI*2.d0
        ENDDO
     ENDDO
  ENDDO


  ! 2)  The derivative of the 1 electron nuclear attraction term ij times
  ! the density matrix element ij.

  ! Please note that these are the three center terms.

  DO Ibas=1,nbasis
     iA=ncenter(Ibas)
     ISTART = (iA-1)*3

     DO Jbas=Ibas,nbasis
        iB = ncenter(Jbas)
        JSTART = (iB-1)*3

        DO iC = 1,natom
           iCSTART = (iC-1)*3

           ! As before, if all terms are on the same atom, they move with the
           ! atom and the dreivative is zero.

           IF (iA == iC .AND. iB == iC) THEN
              continue
           ELSE
              DENSEJI=HOLD(Jbas,Ibas)+HOLD2(Jbas,Ibas)


              ! If Ibas=Jbas, the term only shows up once in the energy, otherwise
              ! it shows up twice. This is not an issue with the 2-center terms above.

              IF (Ibas /= Jbas) DENSEJI=2.d0*DENSEJI


              dNAIX = 0.d0
              dNAIY = 0.d0
              dNAIZ = 0.d0
              dNAJX = 0.d0
              dNAJY = 0.d0
              dNAJZ = 0.d0
              dNACX = 0.d0
              dNACY = 0.d0
              dNACZ = 0.d0

              ! Do the Ibas derivatives.

              DO Icon=1,ncontract(Ibas)
                 DO Jcon=1,ncontract(Jbas)
                    dNAIX = dNAIX + 2.d0*aexp(Icon,Ibas)* &
                         dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                         itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                         itype(1,Ibas)+1,itype(2,Ibas),itype(3,Ibas), &
                         xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                         xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                         xyz(1,iC),xyz(2,iC),xyz(3,iC), chg(iC))
                    dNAIY = dNAIY + 2.d0*aexp(Icon,Ibas)* &
                         dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                         itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                         itype(1,Ibas),itype(2,Ibas)+1,itype(3,Ibas), &
                         xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                         xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                         xyz(1,iC),xyz(2,iC),xyz(3,iC), chg(iC))
                    dNAIZ = dNAIZ + 2.d0*aexp(Icon,Ibas)* &
                         dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                         itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                         itype(1,Ibas),itype(2,Ibas),itype(3,Ibas)+1, &
                         xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                         xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                         xyz(1,iC),xyz(2,iC),xyz(3,iC), chg(iC))

                 ENDDO
              ENDDO


              IF (itype(1,Ibas) /= 0) THEN
                 DO Icon=1,ncontract(Ibas)
                    DO Jcon=1,ncontract(Jbas)
                       dNAIX= dNAIX - dble(itype(1,Ibas))* &
                            dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                            itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                            itype(1,Ibas)-1,itype(2,Ibas),itype(3,Ibas), &
                            xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                            xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                            xyz(1,iC),xyz(2,iC),xyz(3,iC), chg(iC))
                    ENDDO
                 ENDDO
              ENDIF
              IF (itype(2,Ibas) /= 0) THEN
                 DO Icon=1,ncontract(Ibas)
                    DO Jcon=1,ncontract(Jbas)
                       dNAIY= dNAIY - dble(itype(2,Ibas))* &
                            dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                            itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                            itype(1,Ibas),itype(2,Ibas)-1,itype(3,Ibas), &
                            xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                            xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                            xyz(1,iC),xyz(2,iC),xyz(3,iC), chg(iC))
                    ENDDO
                 ENDDO
              ENDIF
              IF (itype(3,Ibas) /= 0) THEN
                 DO Icon=1,ncontract(Ibas)
                    DO Jcon=1,ncontract(Jbas)
                       dNAIZ= dNAIZ - dble(itype(3,Ibas))* &
                            dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                            itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                            itype(1,Ibas),itype(2,Ibas),itype(3,Ibas)-1, &
                            xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                            xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                            xyz(1,iC),xyz(2,iC),xyz(3,iC), chg(iC))
                    ENDDO
                 ENDDO
              ENDIF

              ! Do the Jbas derivatives.

              DO Icon=1,ncontract(Ibas)
                 DO Jcon=1,ncontract(Jbas)
                    dNAJX = dNAJX + 2.d0*aexp(Jcon,Jbas)* &
                         dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                         itype(1,Jbas)+1,itype(2,Jbas),itype(3,Jbas), &
                         itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                         xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                         xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                         xyz(1,iC),xyz(2,iC),xyz(3,iC), chg(iC))
                    dNAJY = dNAJY + 2.d0*aexp(Jcon,Jbas)* &
                         dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                         itype(1,Jbas),itype(2,Jbas)+1,itype(3,Jbas), &
                         itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                         xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                         xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                         xyz(1,iC),xyz(2,iC),xyz(3,iC), chg(iC))
                    dNAJZ = dNAJZ + 2.d0*aexp(Jcon,Jbas)* &
                         dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                         itype(1,Jbas),itype(2,Jbas),itype(3,Jbas)+1, &
                         itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                         xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                         xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                         xyz(1,iC),xyz(2,iC),xyz(3,iC), chg(iC))
                 ENDDO
              ENDDO

              IF (itype(1,Jbas) /= 0) THEN
                 DO Icon=1,ncontract(Ibas)
                    DO Jcon=1,ncontract(Jbas)
                       dNAJX= dNAJX - dble(itype(1,Jbas))* &
                            dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                            itype(1,Jbas)-1,itype(2,Jbas),itype(3,Jbas), &
                            itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                            xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                            xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                            xyz(1,iC),xyz(2,iC),xyz(3,iC), chg(iC))
                    ENDDO
                 ENDDO
              ENDIF
              IF (itype(2,Jbas) /= 0) THEN
                 DO Icon=1,ncontract(Ibas)
                    DO Jcon=1,ncontract(Jbas)
                       dNAJY= dNAJY - dble(itype(2,Jbas))* &
                            dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                            itype(1,Jbas),itype(2,Jbas)-1,itype(3,Jbas), &
                            itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                            xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                            xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                            xyz(1,iC),xyz(2,iC),xyz(3,iC), chg(iC))
                    ENDDO
                 ENDDO
              ENDIF
              IF (itype(3,Jbas) /= 0) THEN
                 DO Icon=1,ncontract(Ibas)
                    DO Jcon=1,ncontract(Jbas)
                       dNAJZ= dNAJZ - dble(itype(3,Jbas))* &
                            dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                            itype(1,Jbas),itype(2,Jbas),itype(3,Jbas)-1, &
                            itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                            xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                            xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                            xyz(1,iC),xyz(2,iC),xyz(3,iC), chg(iC))
                    ENDDO
                 ENDDO
              ENDIF

              ! Now do the derivative with respect to the atom the basis functions
              ! are attracted to.

              DO Icon=1,ncontract(Ibas)
                 DO Jcon=1,ncontract(Jbas)
                    dNACX= dNACX + dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *electricfld(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                         itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                         itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                         1,0,0, &
                         xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                         xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                         xyz(1,iC),xyz(2,iC),xyz(3,iC), chg(iC))
                    dNACY= dNACY + dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *electricfld(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                         itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                         itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                         0,1,0, &
                         xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                         xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                         xyz(1,iC),xyz(2,iC),xyz(3,iC), chg(iC))
                    dNACZ= dNACZ + dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *electricfld(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                         itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                         itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                         0,0,1, &
                         xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                         xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                         xyz(1,iC),xyz(2,iC),xyz(3,iC), chg(iC))
                 ENDDO
              ENDDO

              ! Now add these into the GRADIENT2.

              Gradient2(ISTART+1) = Gradient2(ISTART+1) &
                   +dNAIX*DENSEJI
              Gradient2(ISTART+2) = Gradient2(ISTART+2) &
                   +dNAIY*DENSEJI
              Gradient2(ISTART+3) = Gradient2(ISTART+3) &
                   +dNAIZ*DENSEJI
              Gradient2(JSTART+1) = Gradient2(JSTART+1) &
                   +dNAJX*DENSEJI
              Gradient2(JSTART+2) = Gradient2(JSTART+2) &
                   +dNAJY*DENSEJI
              Gradient2(JSTART+3) = Gradient2(JSTART+3) &
                   +dNAJZ*DENSEJI
              Gradient2(ICSTART+1) = Gradient2(ICSTART+1) &
                   +dNACX*DENSEJI
              Gradient2(ICSTART+2) = Gradient2(ICSTART+2) &
                   +dNACY*DENSEJI
              Gradient2(ICSTART+3) = Gradient2(ICSTART+3) &
                   +dNACZ*DENSEJI
           ENDIF
        ENDDO
     ENDDO
  ENDDO



  ! 3)  The derivative of the 4center 2e- terms with respect to X times
  ! the coefficient found in the energy. (i.e. the multiplicative
  ! constants from the density matrix that arise as these are both
  ! the exchange and correlation integrals.

  ! START HERE
  DO I=1,nbasis
     ! Neglect all the (ii|ii) integrals, as they move with the core.

     DO J=I+1,nbasis

        ! Find  all the (ii|jj) integrals.
        call graddmx2elec(I,I,J,J,GRADIENT2)

        ! Find  all the (ij|jj) integrals.
        call graddmx2elec(I,J,J,J,GRADIENT2)

        ! Find  all the (ii|ij) integrals.
        call graddmx2elec(I,I,I,J,GRADIENT2)

        ! Find all the (ij|ij) integrals
        call graddmx2elec(I,J,I,J,GRADIENT2)

        DO K=J+1,nbasis

           ! Find all the (ij|ik) integrals where j>i,k>j
           call graddmx2elec(I,J,I,K,GRADIENT2)

           ! Find all the (ij|kk) integrals where j>i, k>j.
           call graddmx2elec(I,J,K,K,GRADIENT2)

           ! Find all the (ik|jj) integrals where j>i, k>j.
           call graddmx2elec(I,K,J,J,GRADIENT2)

           ! Find all the (ii|jk) integrals where j>i, k>j.
           call graddmx2elec(I,I,J,K,GRADIENT2)

        ENDDO

        DO K=I+1,nbasis-1

           DO L=K+1,nbasis

              ! Find the (ij|kl) integrals where j>i,k>i,l>k. Note that k and j
              ! can be equal.

              call graddmx2elec(I,J,K,L,GRADIENT2)

           ENDDO
        ENDDO
     ENDDO
  ENDDO

  DO I=1,natom*3
     HESSIAN(I,idX) = HESSIAN(I,idX)+Gradient2(I)
  ENDDO

  return
end subroutine hfdmxderuse



! Ed Brothers. June 3, 2002.
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

subroutine grad2elec(Ibas,Jbas,IIbas,JJbas,coeff)
  use allmod
  implicit double precision(a-h,o-z)

  dimension itype2(3,4)
  logical :: same
  same = .false.

  ! The purpose of this subroutine is to calculate the gradient of
  ! the 2-electron 4-center integrals and add them into the total gradient.
  ! This requires the use of a multiplicitave constant (coeff) which arises
  ! by the way the integral enters the total energy. For example, (ij|ij) is
  ! a repulsion_prim integral and an exchange integral, thus it enters the energy
  ! with a coefficient of (1.5d0*DENSEJI*DENSEJI-0.50d0*DENSEJJ*DENSEII).


  ! First, find the centers the basis functions are located on.  If all the
  ! functions are on the same center, return as this is a zero result.

  iA = ncenter(Ibas)
  iB = ncenter(Jbas)
  iC = ncenter(IIbas)
  iD = ncenter(JJbas)

  same = iA.eq.iB
  same = same .and. iB.eq.iC
  same = same .and. iC.eq.iD

  IF (same) return

  iAstart = (iA-1)*3
  iBstart = (iB-1)*3
  iCstart = (iC-1)*3
  iDstart = (iD-1)*3

  ! The itype2 array was added because if Ibas=Jbas, the code raises two
  ! angular momentums instead of one.

  DO Imomentum=1,3
     itype2(Imomentum,1) = itype(Imomentum,Ibas)
     itype2(Imomentum,2) = itype(Imomentum,Jbas)
     itype2(Imomentum,3) = itype(Imomentum,IIbas)
     itype2(Imomentum,4) = itype(Imomentum,JJbas)
  ENDDO

  ! We have to calculate 12 quantities in this subprogram: the gradient in the
  ! X,Y, and Z directions for the 4 atom A,B,C,D.

  DO Imomentum=1,3
     Agrad=0.d0
     Bgrad=0.d0
     Cgrad=0.d0
     Dgrad=0.d0

     DO Icon = 1, ncontract(Ibas)
        DO Jcon = 1, ncontract(Jbas)
           DO IIcon = 1, ncontract(IIbas)
              DO JJcon = 1, ncontract(JJbas)
                 cntrctcoeff = dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                      *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)

                 itype2(Imomentum,1) = itype2(Imomentum,1)+1
                 Agrad = Agrad+2.d0*aexp(Icon,Ibas)*cntrctcoeff* &
                      repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                      aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                      itype2(1,1), itype2(2,1), itype2(3,1), &
                      itype2(1,2), itype2(2,2), itype2(3,2), &
                      itype2(1,3), itype2(2,3), itype2(3,3), &
                      itype2(1,4), itype2(2,4), itype2(3,4), &
                      xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                      xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                      xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                      xyz(1,iD), xyz(2,iD),xyz(3,iD))
                 itype2(Imomentum,1) = itype2(Imomentum,1)-1

                 IF (itype2(Imomentum,1) /= 0) THEN
                    itype2(Imomentum,1) = itype2(Imomentum,1)-1
                    temp = cntrctcoeff* &
                         repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                         aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                         itype2(1,1), itype2(2,1), itype2(3,1), &
                         itype2(1,2), itype2(2,2), itype2(3,2), &
                         itype2(1,3), itype2(2,3), itype2(3,3), &
                         itype2(1,4), itype2(2,4), itype2(3,4), &
                         xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                         xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                         xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                         xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    itype2(Imomentum,1) = itype2(Imomentum,1)+1
                    Agrad = Agrad-dble(itype2(Imomentum,1))*temp
                 ENDIF

                 itype2(Imomentum,2) = itype2(Imomentum,2)+1
                 Bgrad = Bgrad+2.d0*aexp(Jcon,Jbas)*cntrctcoeff* &
                      repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                      aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                      itype2(1,1), itype2(2,1), itype2(3,1), &
                      itype2(1,2), itype2(2,2), itype2(3,2), &
                      itype2(1,3), itype2(2,3), itype2(3,3), &
                      itype2(1,4), itype2(2,4), itype2(3,4), &
                      xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                      xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                      xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                      xyz(1,iD), xyz(2,iD),xyz(3,iD))
                 itype2(Imomentum,2) = itype2(Imomentum,2)-1
                 IF (itype2(Imomentum,2) /= 0) THEN
                    itype2(Imomentum,2) = itype2(Imomentum,2)-1
                    temp = cntrctcoeff* &
                         repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                         aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                         itype2(1,1), itype2(2,1), itype2(3,1), &
                         itype2(1,2), itype2(2,2), itype2(3,2), &
                         itype2(1,3), itype2(2,3), itype2(3,3), &
                         itype2(1,4), itype2(2,4), itype2(3,4), &
                         xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                         xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                         xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                         xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    itype2(Imomentum,2) = itype2(Imomentum,2)+1
                    Bgrad = Bgrad-dble(itype2(Imomentum,2))*temp
                 ENDIF


                 itype2(Imomentum,3) = itype2(Imomentum,3)+1
                 Cgrad = Cgrad+2.d0*aexp(IIcon,IIbas)*cntrctcoeff* &
                      repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                      aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                      itype2(1,1), itype2(2,1), itype2(3,1), &
                      itype2(1,2), itype2(2,2), itype2(3,2), &
                      itype2(1,3), itype2(2,3), itype2(3,3), &
                      itype2(1,4), itype2(2,4), itype2(3,4), &
                      xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                      xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                      xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                      xyz(1,iD), xyz(2,iD),xyz(3,iD))
                 itype2(Imomentum,3) = itype2(Imomentum,3)-1
                 IF (itype2(Imomentum,3) /= 0) THEN
                    itype2(Imomentum,3) = itype2(Imomentum,3)-1
                    temp = cntrctcoeff* &
                         repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                         aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                         itype2(1,1), itype2(2,1), itype2(3,1), &
                         itype2(1,2), itype2(2,2), itype2(3,2), &
                         itype2(1,3), itype2(2,3), itype2(3,3), &
                         itype2(1,4), itype2(2,4), itype2(3,4), &
                         xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                         xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                         xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                         xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    itype2(Imomentum,3) = itype2(Imomentum,3)+1
                    Cgrad = Cgrad-dble(itype2(Imomentum,3))*temp
                 ENDIF

                 !                        itype2(Imomentum,4) = itype2(Imomentum,4)+1
                 !                        Dgrad = Dgrad+2.d0*aexp(JJcon,JJbas)*cntrctcoeff* &
                 !                        repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                 !                        aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                 !                        itype2(1,1), itype2(2,1), itype2(3,1), &
                 !                        itype2(1,2), itype2(2,2), itype2(3,2), &
                 !                        itype2(1,3), itype2(2,3), itype2(3,3), &
                 !                        itype2(1,4), itype2(2,4), itype2(3,4), &
                 !                        xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                 !                        xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                 !                        xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                 !                        xyz(1,iD), xyz(2,iD),xyz(3,iD))
                 !                        itype2(Imomentum,4) = itype2(Imomentum,4)-1
                 !                        IF (itype2(Imomentum,4) /= 0) THEN
                 !                            itype2(Imomentum,4) = itype2(Imomentum,4)-1
                 !                            temp = cntrctcoeff* &
                 !                            repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                 !                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                 !                            itype2(1,1), itype2(2,1), itype2(3,1), &
                 !                            itype2(1,2), itype2(2,2), itype2(3,2), &
                 !                            itype2(1,3), itype2(2,3), itype2(3,3), &
                 !                            itype2(1,4), itype2(2,4), itype2(3,4), &
                 !                            xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                 !                            xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                 !                            xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                 !                            xyz(1,iD), xyz(2,iD),xyz(3,iD))
                 !                            itype2(Imomentum,4) = itype2(Imomentum,4)+1
                 !                            Dgrad = Dgrad-dble(itype2(Imomentum,4))*temp
                 !                        ENDIF
              ENDDO
           ENDDO
        ENDDO
     ENDDO

     ! Now we have the 4 gradients in a direction, e.g. the X gradient for
     ! atom A,B,C, and D.  Now add it into the gradient time the passed
     ! coefficient.

     Gradient(iASTART+imomentum) = Gradient(iASTART+imomentum)+ &
          AGrad*coeff
     Gradient(iBSTART+imomentum) = Gradient(iBSTART+imomentum)+ &
          BGrad*coeff
     Gradient(iCSTART+imomentum) = Gradient(iCSTART+imomentum)+ &
          CGrad*coeff
     Gradient(iDSTART+imomentum) = Gradient(iDSTART+imomentum) &
          !        DGrad*coeff
          -AGrad*coeff-BGrad*coeff-CGrad*coeff
  ENDDO


  return
end subroutine grad2elec



! Ed Brothers. December 13,2002.
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

subroutine graddmx2elec(Ibas,Jbas,IIbas,JJbas,GRADIENT2)
  use allmod
  implicit double precision(a-h,o-z)

  dimension GRADIENT2(natom*3)
  dimension deriv(4,3),icenter(4),isame(4,8)
  dimension itype2(3,4)
  logical :: same

  ! The purpose of this subroutine is to calculate the gradient of
  ! the 2-electron 4-center integrals and use those with the derivative
  ! of the alpha and beta density matrix derivatives to form the
  ! contribution to the Hessian.

  ! Note that this is basically CPHFB2elec, and could be used to replace
  ! it at a later time.

  iA = ncenter(Ibas)
  iB = ncenter(Jbas)
  iC = ncenter(IIbas)
  iD = ncenter(JJbas)
  icenter(1)=iA
  icenter(2)=iB
  icenter(3)=iC
  icenter(4)=iD

  ! The itype2 array was added because if Ibas=Jbas, the code raises two
  ! angular momentums instead of one.

  DO Imomentum=1,3
     itype2(Imomentum,1) = itype(Imomentum,Ibas)
     itype2(Imomentum,2) = itype(Imomentum,Jbas)
     itype2(Imomentum,3) = itype(Imomentum,IIbas)
     itype2(Imomentum,4) = itype(Imomentum,JJbas)
  ENDDO

  ! We have to calculate 12 quantities in this subprogram: the gradient in the
  ! X,Y, and Z directions for the 4 atom A,B,C,D.

  DO Imomentum=1,3
     Agrad=0.d0
     Bgrad=0.d0
     Cgrad=0.d0
     Dgrad=0.d0

     DO Icon = 1, ncontract(Ibas)
        DO Jcon = 1, ncontract(Jbas)
           DO IIcon = 1, ncontract(IIbas)
              DO JJcon = 1, ncontract(JJbas)
                 cntrctcoeff = dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                      *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)

                 itype2(Imomentum,1) = itype2(Imomentum,1)+1
                 Agrad = Agrad+2.d0*aexp(Icon,Ibas)*cntrctcoeff* &
                      repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                      aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                      itype2(1,1), itype2(2,1), itype2(3,1), &
                      itype2(1,2), itype2(2,2), itype2(3,2), &
                      itype2(1,3), itype2(2,3), itype2(3,3), &
                      itype2(1,4), itype2(2,4), itype2(3,4), &
                      xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                      xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                      xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                      xyz(1,iD), xyz(2,iD),xyz(3,iD))
                 itype2(Imomentum,1) = itype2(Imomentum,1)-1

                 IF (itype2(Imomentum,1) /= 0) THEN
                    itype2(Imomentum,1) = itype2(Imomentum,1)-1
                    temp = cntrctcoeff* &
                         repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                         aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                         itype2(1,1), itype2(2,1), itype2(3,1), &
                         itype2(1,2), itype2(2,2), itype2(3,2), &
                         itype2(1,3), itype2(2,3), itype2(3,3), &
                         itype2(1,4), itype2(2,4), itype2(3,4), &
                         xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                         xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                         xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                         xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    itype2(Imomentum,1) = itype2(Imomentum,1)+1
                    Agrad = Agrad-dble(itype2(Imomentum,1))*temp
                 ENDIF

                 itype2(Imomentum,2) = itype2(Imomentum,2)+1
                 Bgrad = Bgrad+2.d0*aexp(Jcon,Jbas)*cntrctcoeff* &
                      repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                      aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                      itype2(1,1), itype2(2,1), itype2(3,1), &
                      itype2(1,2), itype2(2,2), itype2(3,2), &
                      itype2(1,3), itype2(2,3), itype2(3,3), &
                      itype2(1,4), itype2(2,4), itype2(3,4), &
                      xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                      xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                      xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                      xyz(1,iD), xyz(2,iD),xyz(3,iD))
                 itype2(Imomentum,2) = itype2(Imomentum,2)-1
                 IF (itype2(Imomentum,2) /= 0) THEN
                    itype2(Imomentum,2) = itype2(Imomentum,2)-1
                    temp = cntrctcoeff* &
                         repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                         aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                         itype2(1,1), itype2(2,1), itype2(3,1), &
                         itype2(1,2), itype2(2,2), itype2(3,2), &
                         itype2(1,3), itype2(2,3), itype2(3,3), &
                         itype2(1,4), itype2(2,4), itype2(3,4), &
                         xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                         xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                         xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                         xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    itype2(Imomentum,2) = itype2(Imomentum,2)+1
                    Bgrad = Bgrad-dble(itype2(Imomentum,2))*temp
                 ENDIF


                 itype2(Imomentum,3) = itype2(Imomentum,3)+1
                 Cgrad = Cgrad+2.d0*aexp(IIcon,IIbas)*cntrctcoeff* &
                      repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                      aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                      itype2(1,1), itype2(2,1), itype2(3,1), &
                      itype2(1,2), itype2(2,2), itype2(3,2), &
                      itype2(1,3), itype2(2,3), itype2(3,3), &
                      itype2(1,4), itype2(2,4), itype2(3,4), &
                      xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                      xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                      xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                      xyz(1,iD), xyz(2,iD),xyz(3,iD))
                 itype2(Imomentum,3) = itype2(Imomentum,3)-1
                 IF (itype2(Imomentum,3) /= 0) THEN
                    itype2(Imomentum,3) = itype2(Imomentum,3)-1
                    temp = cntrctcoeff* &
                         repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                         aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                         itype2(1,1), itype2(2,1), itype2(3,1), &
                         itype2(1,2), itype2(2,2), itype2(3,2), &
                         itype2(1,3), itype2(2,3), itype2(3,3), &
                         itype2(1,4), itype2(2,4), itype2(3,4), &
                         xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                         xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                         xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                         xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    itype2(Imomentum,3) = itype2(Imomentum,3)+1
                    Cgrad = Cgrad-dble(itype2(Imomentum,3))*temp
                 ENDIF

                 itype2(Imomentum,4) = itype2(Imomentum,4)+1
                 Dgrad = Dgrad+2.d0*aexp(JJcon,JJbas)*cntrctcoeff* &
                      repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                      aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                      itype2(1,1), itype2(2,1), itype2(3,1), &
                      itype2(1,2), itype2(2,2), itype2(3,2), &
                      itype2(1,3), itype2(2,3), itype2(3,3), &
                      itype2(1,4), itype2(2,4), itype2(3,4), &
                      xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                      xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                      xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                      xyz(1,iD), xyz(2,iD),xyz(3,iD))
                 itype2(Imomentum,4) = itype2(Imomentum,4)-1
                 IF (itype2(Imomentum,4) /= 0) THEN
                    itype2(Imomentum,4) = itype2(Imomentum,4)-1
                    temp = cntrctcoeff* &
                         repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                         aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                         itype2(1,1), itype2(2,1), itype2(3,1), &
                         itype2(1,2), itype2(2,2), itype2(3,2), &
                         itype2(1,3), itype2(2,3), itype2(3,3), &
                         itype2(1,4), itype2(2,4), itype2(3,4), &
                         xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                         xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                         xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                         xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    itype2(Imomentum,4) = itype2(Imomentum,4)+1
                    Dgrad = Dgrad-dble(itype2(Imomentum,4))*temp
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     ! Now we have the 4 gradients in a direction, e.g. the X gradient for
     ! atom A,B,C, and D.  Now add it into the gradient time the passed
     ! coefficient.

     deriv(1,imomentum) = Agrad
     deriv(2,imomentum) = Bgrad
     deriv(3,imomentum) = Cgrad
     deriv(4,imomentum) = Dgrad
  ENDDO


  ! Now check to see if any of the centers are redundant, starting from
  ! the end of the array.

  DO J=4,2,-1
     DO I=J-1,1,-1
        IF (icenter(I) == icenter(J) .AND. icenter(J) > 0) THEN
           DO K=1,3
              deriv(I,K) = deriv(J,K) + deriv(I,K)
           ENDDO
           icenter(J) = 0
        ENDIF
     ENDDO
  ENDDO

  ! Now we need to examine the symmetry of the integral.
  ! For an integral (ij|kl):
  ! (ij|kl)=(ji|kl)=(ij|lk)=(ji|lk)=(kl|ij)=(lk|ij)=(kl|ji)=(lk|ji)

  ! set up (ij|kl)
  isame(1,1)=Ibas
  isame(2,1)=Jbas
  isame(3,1)=IIbas
  isame(4,1)=JJbas
  ! set up (ji|kl)
  isame(1,2)=Jbas
  isame(2,2)=Ibas
  isame(3,2)=IIbas
  isame(4,2)=JJbas
  ! set up (ij|lk)
  isame(1,3)=Ibas
  isame(2,3)=Jbas
  isame(3,3)=JJbas
  isame(4,3)=IIbas
  ! set up (ji|lk)
  isame(1,4)=Jbas
  isame(2,4)=Ibas
  isame(3,4)=JJbas
  isame(4,4)=IIbas
  ! set up (kl|ij)
  isame(1,5)=IIbas
  isame(2,5)=JJbas
  isame(3,5)=Ibas
  isame(4,5)=Jbas
  ! set up (lk|ij)
  isame(1,6)=JJbas
  isame(2,6)=IIbas
  isame(3,6)=Ibas
  isame(4,6)=Jbas
  ! set up (kl|ji)
  isame(1,7)=IIbas
  isame(2,7)=JJbas
  isame(3,7)=Jbas
  isame(4,7)=Ibas
  ! set up (lk|ji)
  isame(1,8)=JJbas
  isame(2,8)=IIbas
  isame(3,8)=Jbas
  isame(4,8)=Ibas

  ! Now we check for redundancy.

  DO Icheck=1,8
     IF (isame(1,Icheck) /= 0) THEN
        DO Jcheck=Icheck+1,8
           IF (isame(1,Jcheck) /= 0) THEN
              same = isame(1,Icheck).eq.isame(1,Jcheck)
              same = same.and.isame(2,Icheck).eq.isame(2,Jcheck)
              same = same.and.isame(3,Icheck).eq.isame(3,Jcheck)
              same = same.and.isame(4,Icheck).eq.isame(4,Jcheck)
              IF (same) THEN
                 DO Iblank=1,4
                    isame(Iblank,Jcheck)=0
                 ENDDO
              ENDIF
           ENDIF
        ENDDO
     ENDIF
  ENDDO

  ! Now we add it into gradient2.

  DO Iatm = 1,4
     IF (icenter(Iatm) /= 0) THEN
        DO Imomentum=1,3
           ISTART = (icenter(iatm)-1)*3+Imomentum
           currderiv=deriv(Iatm,Imomentum)

           ! Iatm and Imomentum define what atom is being moved in what cartesian
           ! direction, and thus where it goes in Gradient2.  Now we loop over
           ! integral symmetries to find the constant to multiply the integral
           ! derivative times.

           DO Isym = 1,8
              IF (isame(1,ISYM) /= 0) THEN
                 iA = isame(1,ISYM)
                 iB = isame(2,ISYM)
                 iC = isame(3,ISYM)
                 iD = isame(4,ISYM)
                 IF (quick_method%unrst) THEN
                    coeff = (DENSE(iC,iD) + DENSEB(iC,iD)) * &
                         (HOLD(iA,iB) + HOLD2(iA,iB))
                    coeff = coeff - (HOLD(iA,iC)*DENSE(iB,iD))
                    coeff = coeff - (HOLD2(iA,iC)*DENSEB(iB,iD))
                 ELSE
                    coeff = DENSE(iC,iD) *(HOLD(iA,iB) + HOLD2(iA,iB))
                    coeff = coeff - (HOLD(iA,iC)*.5d0*DENSE(iB,iD))
                    coeff = coeff - (HOLD2(iA,iC)*.5d0*DENSE(iB,iD))
                 ENDIF
                 GRADIENT2(ISTART)=GRADIENT2(ISTART) + currderiv*coeff
              ENDIF
           ENDDO
        ENDDO
     ENDIF
  ENDDO

  return
end subroutine graddmx2elec



! Ed Brothers. December 2, 2002.
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

subroutine dmxderiv(IDX)
  use allmod
  implicit double precision(a-h,o-z)


  ! The purpose of the subroutine is to calculate the first derivative of
  ! the alpha and beta density matrices given the u(ai) found in the
  ! CPHF procedure.  The results are stored in HOLD and HOLD2. Note that
  ! these arrays should be symmetric.

  ! DENSE(Jbas,Ibas) = (Sum over orbitals occ i) CO(jbas,I) CO(Ibas,i)
  ! dDENSE(Jbas,Ibas)/dy = (Sum over orbitals occ i) dCO(jbas,I)/dy CO(Ibas,i) +
  ! CO(Jbas,i) dCO(Ibas,i)/dy

  ! dCO(Jbas,i)/dy  = (Sum over orbitals r) CO(Jbas,r) u(i,r)
  ! THUS

  ! dDENSE(Jbas,Ibas)/dy = (Sum over orbitals occ i) (Sum over orbitals r)
  ! CO(Jbas,r) u(i,r) CO(Ibas,i) + CO(Jbas,i) CO(Ibas,r) u(i,r)

  ! dDENSE(Jbas,Ibas)/dy = (Sum over orbitals occ i) (Sum over orbitals r)
  ! u(i,r) (CO(Jbas,r) CO(Ibas,i) + CO(Jbas,i) CO(Ibas,r))


  ! First we blank out the arrays.

  DO Ibas=1,nbasis
     DO Jbas=1,nbasis
        HOLD(Jbas,Ibas)=0.d0
     ENDDO
  ENDDO
  DO Ibas=1,nbasis
     DO Jbas=1,nbasis
        HOLD2(Jbas,Ibas)=0.d0
     ENDDO
  ENDDO

  ! We also need to find out from IDX which atom we are moving and in
  ! which direction.

  Imomentum=mod(idx,3)
  IF (Imomentum == 0) THEN
     Iatom=IDX/3
     Imomentum = 3
  ELSE
     Iatom = (IDX-Imomentum)/3 + 1
  ENDIF


  ! Now we need to find out where the alpha and beta occupied/virtual
  ! lines are.

  IF (quick_method%unrst) THEN
     lastAocc = nelec
     lastBocc = nelecb
  ELSE
     lastAocc = nelec/2
     lastBocc = lastAocc
  ENDIF
  iBetastart = lastAocc*(nbasis-lastAocc)

  ! In the case of r and s both being occupied (call them k and l)

  ! u(kl) = -1/2  S'kl

  ! Thus we calculate this and add it to the density matrix derivative.


  ! Do the alpha first.

  DO iAocc = 1,lastAocc
     DO iAocc2 = 1,lastAocc

        ! K and L are selected.  IDX gave us atom and direction of perturbation.


        ! Now we loop over basis functions.  Note that Kbas functions are always
        ! on center Katom, and the Lbas functions are always not on that center.
        ! This actually calculates the Skl

        Skl  = 0.d0
        DO Kbas = ifirst(Iatom),ilast(Iatom)
           DO Lbas = 1,nbasis
              IF (Lbas < ifirst(Iatom) .OR. Lbas > ilast(Iatom)) THEN
                 dSK=0.d0
                 itype(Imomentum,Kbas) = itype(Imomentum,Kbas)+1
                 DO Kcon=1,ncontract(Kbas)
                    DO Lcon=1,ncontract(Lbas)
                       dSK = dSK + 2.d0*aexp(Kcon,Kbas)* &
                            dcoeff(Lcon,Lbas)*dcoeff(Kcon,Kbas) &
                            *overlap(aexp(Lcon,Lbas),aexp(Kcon,Kbas), &
                            itype(1,Lbas),itype(2,Lbas),itype(3,Lbas), &
                            itype(1,Kbas),itype(2,Kbas),itype(3,Kbas), &
                            xyz(1,ncenter(Lbas)),xyz(2,ncenter(Lbas)), &
                            xyz(3,ncenter(Lbas)),xyz(1,ncenter(Kbas)), &
                            xyz(2,ncenter(Kbas)),xyz(3,ncenter(Kbas)))
                    ENDDO
                 ENDDO
                 itype(Imomentum,Kbas) = itype(Imomentum,Kbas)-1
                 IF (itype(Imomentum,Kbas) /= 0) THEN
                    itype(Imomentum,Kbas) = itype(Imomentum,Kbas)-1
                    DO Kcon=1,ncontract(Kbas)
                       DO Lcon=1,ncontract(Lbas)
                          dSK = dSK - dble(itype(Imomentum,Kbas)+1)* &
                               dcoeff(Lcon,Lbas)*dcoeff(Kcon,Kbas) &
                               *overlap(aexp(Lcon,Lbas),aexp(Kcon,Kbas), &
                               itype(1,Lbas),itype(2,Lbas),itype(3,Lbas), &
                               itype(1,Kbas),itype(2,Kbas),itype(3,Kbas), &
                               xyz(1,ncenter(Lbas)),xyz(2,ncenter(Lbas)), &
                               xyz(3,ncenter(Lbas)),xyz(1,ncenter(Kbas)), &
                               xyz(2,ncenter(Kbas)),xyz(3,ncenter(Kbas)))
                       ENDDO
                    ENDDO
                    itype(Imomentum,Kbas) = itype(Imomentum,Kbas)+1
                 ENDIF
                 Skl=Skl+dSK*(CO(Kbas,iAocc)*CO(Lbas,iAocc2) + &
                      CO(Lbas,iAocc)*CO(Kbas,iAocc2))
              ENDIF
           ENDDO
        ENDDO

        ! Now that we have Skl, we need to add it into the HOLD array, i.e.
        ! the density matrix derivative.

        DO Ibas=1,nbasis
           DO Jbas=Ibas,nbasis
              HOLD(Jbas,Ibas)= HOLD(Jbas,Ibas)-.5d0*Skl* &
                   (CO(Jbas,iAocc)*CO(Ibas,iAocc2) + &
                   CO(Jbas,iAocc2)*CO(Ibas,iAocc))
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  ! Now we need to repeat the whole process for the beta kl pairs.


  DO iBocc = 1,lastBocc
     DO iBocc2 = 1,lastBocc
        ! Now we loop over basis functions.  Note that Ibas functions are always
        ! on center Iatom, and the Jbas functions are always not on that center.
        ! This actually calculates the Skl

        Skl  = 0.d0
        DO Kbas = ifirst(Iatom),ilast(Iatom)
           DO Lbas = 1,nbasis
              IF (Lbas < ifirst(Iatom) .OR. Lbas > ilast(Iatom)) THEN
                 dSK=0.d0
                 itype(Imomentum,Kbas) = itype(Imomentum,Kbas)+1
                 DO Kcon=1,ncontract(Kbas)
                    DO Lcon=1,ncontract(Lbas)
                       dSK = dSK + 2.d0*aexp(Kcon,Kbas)* &
                            dcoeff(Lcon,Lbas)*dcoeff(Kcon,Kbas) &
                            *overlap(aexp(Lcon,Lbas),aexp(Kcon,Kbas), &
                            itype(1,Lbas),itype(2,Lbas),itype(3,Lbas), &
                            itype(1,Kbas),itype(2,Kbas),itype(3,Kbas), &
                            xyz(1,ncenter(Lbas)),xyz(2,ncenter(Lbas)), &
                            xyz(3,ncenter(Lbas)),xyz(1,ncenter(Kbas)), &
                            xyz(2,ncenter(Kbas)),xyz(3,ncenter(Kbas)))
                    ENDDO
                 ENDDO
                 itype(Imomentum,Kbas) = itype(Imomentum,Kbas)-1
                 IF (itype(Imomentum,Kbas) /= 0) THEN
                    itype(Imomentum,Kbas) = itype(Imomentum,Kbas)-1
                    DO Kcon=1,ncontract(Kbas)
                       DO Lcon=1,ncontract(Lbas)
                          dSK = dSK - dble(itype(Imomentum,Kbas)+1)* &
                               dcoeff(Lcon,Lbas)*dcoeff(Kcon,Kbas) &
                               *overlap(aexp(Lcon,Lbas),aexp(Kcon,Kbas), &
                               itype(1,Lbas),itype(2,Lbas),itype(3,Lbas), &
                               itype(1,Kbas),itype(2,Kbas),itype(3,Kbas), &
                               xyz(1,ncenter(Lbas)),xyz(2,ncenter(Lbas)), &
                               xyz(3,ncenter(Lbas)),xyz(1,ncenter(Kbas)), &
                               xyz(2,ncenter(Kbas)),xyz(3,ncenter(Kbas)))
                       ENDDO
                    ENDDO
                    itype(Imomentum,Kbas) = itype(Imomentum,Kbas)+1
                 ENDIF
                 Skl=Skl+dSK*(COB(Kbas,iBocc)*COB(Lbas,iBocc2) + &
                      COB(Lbas,iBocc)*COB(Kbas,iBocc2))
              ENDIF
           ENDDO
        ENDDO

        ! Now that we have Skl, we need to add it into the HOLD2 array, i.e.
        ! the Beta density matrix derivative.

        DO Ibas=1,nbasis
           DO Jbas=Ibas,nbasis
              HOLD2(Jbas,Ibas)= HOLD2(Jbas,Ibas)-.5d0*Skl* &
                   (COB(Jbas,iBocc)*COB(Ibas,iBocc2) + &
                   COB(Jbas,iBocc2)*COB(Ibas,iBocc))
           ENDDO
        ENDDO

     ENDDO
  ENDDO


  ! At this point we have done the occupied-occupied section.  Now do the
  ! virtual occupied section.  Note that this gives us the completed
  ! derivative of the density matrix.


  ! Alpha first.

  DO iAvirt = lastAocc+1,nbasis
     DO iAocc = 1,lastAocc

        ! iAvirt and iAocc form an ai pair.  Find it's location.

        iaCPHF = (iAvirt-lastAocc-1)*lastAocc + iAocc

        ! Now add this into the dmx derivative.

        DO Ibas=1,nbasis
           DO Jbas=Ibas,nbasis
              HOLD(Jbas,Ibas)= HOLD(Jbas,Ibas)+BU(iaCPHF)* &
                   (CO(Ibas,iAocc)*CO(Jbas,iAvirt)+ &
                   CO(Ibas,iAvirt)*CO(Jbas,iAocc))
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  ! Now Beta.
  DO iBvirt = lastBocc+1,nbasis
     DO iBocc = 1,lastBocc

        ! iBvirt and iBocc form an ai pair.  Find it's location.

        iaCPHF = (iBvirt-lastBocc-1)*lastBocc + iBocc + ibetastart

        ! Now add this into the dmx derivative.

        DO Ibas=1,nbasis
           DO Jbas=Ibas,nbasis
              HOLD2(Jbas,Ibas)= HOLD2(Jbas,Ibas)+BU(iaCPHF)* &
                   (COB(Ibas,iBocc)*COB(Jbas,iBvirt)+ &
                   COB(Ibas,iBvirt)*COB(Jbas,iBocc))
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  ! Now copy the lower diagonal to the upper diagonal.

  DO I=1,nbasis
     DO J=I+1,nbasis
        HOLD(I,J) =  HOLD(J,I)
        HOLD2(I,J) =  HOLD2(J,I)
     ENDDO
  ENDDO

  return
end subroutine dmxderiv



! Ed Brothers. May 28, 2002
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

double precision function electricfld(a,b,i,j,k,ii,jj,kk, &
     idx,idy,idz,Ax,Ay,Az, &
     Bx,By,Bz,Cx,Cy,Cz,Z)
  use quick_constants_module
  implicit double precision(a-h,o-z)
  dimension aux(0:20)

  ! Variables needed later:
  !    pi=3.1415926535897932385

  g = a+b
  Px = (a*Ax + b*Bx)/g
  Py = (a*Ay + b*By)/g
  Pz = (a*Az + b*Bz)/g

  PCsquare = (Px-Cx)**2.d0 + (Py -Cy)**2.d0 + (Pz -Cz)**2.d0

  ! The purpose of this subroutine is to calculate the derivative of the
  ! nuclear attraction integral with respect to the nuclear displacement
  ! of the atom the electronic distribution is attracted to:

  ! d/dXC (Integral over all space) Phi(mu) Phi(nu) 1/rC

  ! The notation is the same used throughout: gtfs with orbital exponents a
  ! and b on A and B with angular momentums defined by i,j,k (a's x, y
  ! and z exponents, respectively) and ii,jj,k and kk on B with the core at
  ! (Cx,Cy,Cz) with charge Z. m is the "order" of the integral which
  ! arises from the recusion relationship. New to this are the idx, idy, and
  ! idz terms which denote derivatives in the x y and z direction for the
  ! C atom.

  ! The this is taken from the recursive relation found in Obara and Saika,
  ! J. Chem. Phys. 84 (7) 1986, 3963.

  ! The first step is generating all the necessary auxillary integrals.
  ! These are (0|1/rc|0)^(m) = 2 Sqrt (g/Pi) (0||0) Fm(g(Rpc)^2)
  ! The values of m range from 0 to i+j+k+ii+jj+kk+2. This is exactly the
  ! same as in the attraction code, and is necessary as we will be calling
  ! that code eventually.

  U = g* PCsquare
  Maxm = i+j+k+ii+jj+kk+2
  call FmT(Maxm,U,aux)
  constant = overlap(a,b,0,0,0,0,0,0,Ax,Ay,Az,Bx,By,Bz) &
       * 2.d0 * (g/Pi)**0.5d0
  DO L = 0,maxm
     aux(L) = aux(L)*constant
  ENDDO

  ! At this point all the auxillary integrals have been calculated.
  ! It is now time to decompase the attraction integral to it's
  ! auxillary integrals through the recursion scheme.  To do this we use
  ! a recursive function.

  electricfld = -1.d0*Z*elctfldrecurse(i,j,k,ii,jj,kk,idx,idy,idz, &
       0,aux,Ax,Ay,Az,Bx,By,Bz, &
       Cx,Cy,Cz,Px,Py,Pz,g)

  return
end function electricfld



! Ed Brothers. December 12, 2002
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

subroutine ewtdmxder(IDX)
  use allmod
  implicit double precision(a-h,o-z)
  dimension temp(nbasis,nbasis),ewtdmx(nbasis,nbasis)

  ! This program calculates and uses the first derivative of the energy
  ! weighted density matrix.

  ! The energy weighted density matrix is:

  ! W = ALPHA(DENSE.OPERATOR.DENSE)
  ! Wx = ALPHA(DENSEx.OPERATOR.DENSE + DENSE.OPERATORx.DENSE
  ! + DENSE.OPERATOR.DENSEx) + BETA(DENSEx.OPERATOR.DENSE
  ! + DENSE.OPERATORx.DENSE + DENSE.OPERATOR.DENSEx)

  ! First we are going to do something sloppy.  If this is a restricted run
  ! split DENSE into DENSE and DENSEB.

  IF ( .NOT. quick_method%unrst) THEN
     DO I=1,nbasis
        DO J=1,nbasis
           DENSEB(J,I) = .5d0*DENSE(J,I)
           DENSE(J,I) = .5d0*DENSE(J,I)
        ENDDO
     ENDDO
  ENDIF

  ! Now get the alpha operator matrix and start building the first derivative
  ! of the energy weighted density matrix.

  call uhfoperatorA

  DO I=1,nbasis
     DO J=1,nbasis
        tempIJ = 0.0D0
        DO K=1,nbasis
           tempIJ = tempIJ + O(K,I)*HOLD(K,J)
        ENDDO
        temp(I,J) = tempIJ
     ENDDO
  ENDDO

  DO I=1,nbasis
     DO J=1,nbasis
        ewtdmxIJ = 0.0D0
        DO K=1,nbasis
           ewtdmxIJ = ewtdmxIJ + DENSE(K,J)*temp(K,I)
        ENDDO
        ewtdmx(I,J) = ewtdmxIJ
     ENDDO
  ENDDO

  ! ewtdmx now contains ALPHA(DENSE.OPERATOR.DENSEx)

  DO I=1,nbasis
     DO J=1,nbasis
        tempIJ = 0.0D0
        DO K=1,nbasis
           tempIJ = tempIJ + O(K,I)*DENSE(K,J)
        ENDDO
        temp(I,J) = tempIJ
     ENDDO
  ENDDO

  DO I=1,nbasis
     DO J=1,nbasis
        ewtdmxIJ = 0.0D0
        DO K=1,nbasis
           ewtdmxIJ = ewtdmxIJ + HOLD(K,J)*temp(K,I)
        ENDDO
        ewtdmx(I,J) = ewtdmx(I,J)+ewtdmxIJ
     ENDDO
  ENDDO

  ! ewtdmx now contains ALPHA(DENSE.OPERATOR.DENSEx)
  ! +ALPHA(DENSEx.OPERATOR.DENSE)

  ! Now get the beta operator matrix.

  call uhfoperatorB

  DO I=1,nbasis
     DO J=1,nbasis
        tempIJ = 0.0D0
        DO K=1,nbasis
           tempIJ = tempIJ + O(K,I)*HOLD2(K,J)
        ENDDO
        temp(I,J) = tempIJ
     ENDDO
  ENDDO

  DO I=1,nbasis
     DO J=1,nbasis
        ewtdmxIJ = 0.0D0
        DO K=1,nbasis
           ewtdmxIJ = ewtdmxIJ + DENSEB(K,J)*temp(K,I)
        ENDDO
        ewtdmx(I,J) = ewtdmx(I,J)+ewtdmxIJ
     ENDDO
  ENDDO

  ! ewtdmx now contains ALPHA(DENSE.OPERATOR.DENSEx)
  ! +ALPHA(DENSEx.OPERATOR.DENSE)
  ! +BETA(DENSE.OPERATOR.DENSEx)

  DO I=1,nbasis
     DO J=1,nbasis
        tempIJ = 0.0D0
        DO K=1,nbasis
           tempIJ = tempIJ + O(K,I)*DENSEB(K,J)
        ENDDO
        temp(I,J) = tempIJ
     ENDDO
  ENDDO

  DO I=1,nbasis
     DO J=1,nbasis
        ewtdmxIJ = 0.0D0
        DO K=1,nbasis
           ewtdmxIJ = ewtdmxIJ + HOLD2(K,J)*temp(K,I)
        ENDDO
        ewtdmx(I,J) = ewtdmx(I,J)+ewtdmxIJ
     ENDDO
  ENDDO


  ! ewtdmx now contains ALPHA(DENSE.OPERATOR.DENSEx)
  ! +ALPHA(DENSEx.OPERATOR.DENSE)
  ! +BETA(DENSE.OPERATOR.DENSEx)
  ! +BETA(DENSEx.OPERATOR.DENSE)

  ! Now we need to calculate the 1st derivative of the operator matrices.

  call duhfoperatorA(IDX)

  DO I=1,nbasis
     DO J=1,nbasis
        tempIJ = 0.0D0
        DO K=1,nbasis
           tempIJ = tempIJ + O(K,I)*DENSE(K,J)
        ENDDO
        temp(I,J) = tempIJ
     ENDDO
  ENDDO

  DO I=1,nbasis
     DO J=1,nbasis
        ewtdmxIJ = 0.0D0
        DO K=1,nbasis
           ewtdmxIJ = ewtdmxIJ + DENSE(K,J)*temp(K,I)
        ENDDO
        ewtdmx(I,J) = ewtdmx(I,J)+ewtdmxIJ
     ENDDO
  ENDDO

  ! ewtdmx now contains ALPHA(DENSE.OPERATOR.DENSEx)
  ! +ALPHA(DENSEx.OPERATOR.DENSE)
  ! +BETA(DENSE.OPERATOR.DENSEx)
  ! +BETA(DENSEx.OPERATOR.DENSE)
  ! +ALPHA(DENSE.OPERATORx.DENSE)


  call duhfoperatorB(IDX)
  DO I=1,nbasis
     DO J=1,nbasis
        tempIJ = 0.0D0
        DO K=1,nbasis
           tempIJ = tempIJ + O(K,I)*DENSEB(K,J)
        ENDDO
        temp(I,J) = tempIJ
     ENDDO
  ENDDO

  DO I=1,nbasis
     DO J=1,nbasis
        ewtdmxIJ = 0.0D0
        DO K=1,nbasis
           ewtdmxIJ = ewtdmxIJ + DENSEB(K,J)*temp(K,I)
        ENDDO
        ewtdmx(I,J) = ewtdmx(I,J)+ewtdmxIJ
     ENDDO
  ENDDO

  ! ewtdmx now contains ALPHA(DENSE.OPERATOR.DENSEx)
  ! +ALPHA(DENSEx.OPERATOR.DENSE)
  ! +BETA(DENSE.OPERATOR.DENSEx)
  ! +BETA(DENSEx.OPERATOR.DENSE)
  ! +ALPHA(DENSE.OPERATORx.DENSE)
  ! +BETA(DENSE.OPERATORx.DENSE)

  ! This is the complete derivative of the energy weighted density matrix.
  ! Now we use it with the derivative of the overlap element to form the
  ! final hessian.

  DO Ibas=1,nbasis
     ISTART = (ncenter(Ibas)-1) *3
     DO Jbas=ilast(ncenter(IBAS))+1,nbasis
        JSTART = (ncenter(Jbas)-1) *3
        DENSEJI = DENSE(Jbas,Ibas)

        ! We have selected our two basis functions, now loop over angular momentum.

        DO Imomentum=1,3
           dSI = 0.d0
           dSJ =0.d0

           ! Do the Ibas derivatives first.

           itype(Imomentum,Ibas) = itype(Imomentum,Ibas)+1
           DO Icon=1,ncontract(Ibas)
              DO Jcon=1,ncontract(Jbas)
                 dSI = dSI + 2.d0*aexp(Icon,Ibas)* &
                      dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                      *overlap(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                      itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                      itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                      xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                      xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                      xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
              ENDDO
           ENDDO
           itype(Imomentum,Ibas) = itype(Imomentum,Ibas)-1
           IF (itype(Imomentum,Ibas) /= 0) THEN
              itype(Imomentum,Ibas) = itype(Imomentum,Ibas)-1
              DO Icon=1,ncontract(Ibas)
                 DO Jcon=1,ncontract(Jbas)
                    dSI = dSI - dble(itype(Imomentum,Ibas)+1)* &
                         dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *overlap(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                         itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                         itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                         xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                         xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                         xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                 ENDDO
              ENDDO
              itype(Imomentum,Ibas) = itype(Imomentum,Ibas)+1
           ENDIF
           HESSIAN(ISTART+Imomentum,IDX)=HESSIAN(ISTART+Imomentum,IDX) &
                -2.d0*dSI*ewtdmx(Jbas,Ibas)

           ! Now do the Jbas derivatives.

           itype(Imomentum,Jbas) = itype(Imomentum,Jbas)+1
           DO Icon=1,ncontract(Ibas)
              DO Jcon=1,ncontract(Jbas)
                 dSJ = dSJ + 2.d0*aexp(Jcon,Jbas)* &
                      dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                      *overlap(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                      itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                      itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                      xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                      xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                      xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
              ENDDO
           ENDDO
           itype(Imomentum,Jbas) = itype(Imomentum,Jbas)-1
           IF (itype(Imomentum,Jbas) /= 0) THEN
              itype(Imomentum,Jbas) = itype(Imomentum,Jbas)-1
              DO Icon=1,ncontract(Ibas)
                 DO Jcon=1,ncontract(Jbas)
                    dSJ = dSJ - dble(itype(Imomentum,Jbas)+1)* &
                         dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *overlap(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                         itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                         itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                         xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                         xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                         xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                 ENDDO
              ENDDO
              itype(Imomentum,Jbas) = itype(Imomentum,Jbas)+1
           ENDIF
           HESSIAN(JSTART+Imomentum,IDX)=HESSIAN(JSTART+Imomentum,IDX) &
                -2.d0*dSJ*ewtdmx(Jbas,Ibas)
        ENDDO
     ENDDO
  ENDDO


  ! Before returning, if this is a restricted run, reform the density matrix.

  IF ( .NOT. quick_method%unrst) THEN
     DO I=1,nbasis
        DO J=1,nbasis
           DENSE(J,I) = 2.d0*DENSE(J,I)
        ENDDO
     ENDDO
  ENDIF

END subroutine ewtdmxder





! Ed Brothers. December 16, 2002
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

subroutine duhfoperatora(IDX)
  use allmod
  implicit double precision(a-h,o-z)
  dimension igrad(3)
  logical :: IonMove, JonMove, ConMove

  ! The purpose of this subroutine is to form the first derivative
  ! of the alpha fock matrix for use in forming the first derivative
  ! of the energy weighted density matrix.

  ! Note that this matrix is symmetric.

  ! Blank out the array.

  DO I=1,nbasis
     DO J=1,nbasis
        O(J,I)=0.d0
     ENDDO
  ENDDO

  ! The first thing to do is find the atom and direction of the perturbation.

  Imomentum=mod(idx,3)
  IF (Imomentum == 0) THEN
     Iatom=IDX/3
     Imomentum = 3
  ELSE
     Iatom = (IDX-Imomentum)/3 + 1
  ENDIF

  ! We need to take the derivative of the one electron portion of the
  ! Hamiltonian.

  DO Ibas=1,nbasis
     IonMove = ncenter(Ibas).eq.Iatom
     DO Jbas = Ibas,nbasis
        JonMove = ncenter(Jbas).eq.Iatom

        ! If I is on Iatom (the moving atom) and J is not, the derivative  of
        ! the kinetic energy is non-zero.

        dJI=0.d0

        IF (IonMove .AND. .NOT. Jonmove) THEN
           DO Icon=1,ncontract(ibas)
              DO Jcon=1,ncontract(jbas)
                 itype(Imomentum,Ibas) = itype(Imomentum,Ibas)+1
                 dJI = dJI + 2.d0*aexp(Icon,Ibas)* &
                      dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                      *ekinetic(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                      itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                      itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                      xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                      xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                      xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                 itype(Imomentum,Ibas) = itype(Imomentum,Ibas)-1
                 IF (itype(Imomentum,Ibas) /= 0) THEN
                    itype(Imomentum,Ibas) = itype(Imomentum,Ibas)-1
                    dJI = dJI - dble(itype(Imomentum,Ibas)+1)* &
                         dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *ekinetic(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                         itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                         itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                         xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                         xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                         xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                    itype(Imomentum,Ibas) = itype(Imomentum,Ibas)+1
                 ENDIF
              ENDDO
           ENDDO
        ENDIF

        ! If J is on Iatom (the moving atom) and I is not, the derivative  of
        ! the kinetic energy is non-zero.

        IF (JonMove .AND. .NOT. Ionmove) THEN
           DO Icon=1,ncontract(ibas)
              DO Jcon=1,ncontract(jbas)
                 itype(Imomentum,Jbas) = itype(Imomentum,Jbas)+1
                 dJI = dJI + 2.d0*aexp(Jcon,Jbas)* &
                      dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                      *ekinetic(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                      itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                      itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                      xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                      xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                      xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                 itype(Imomentum,Jbas) = itype(Imomentum,Jbas)-1
                 IF (itype(Imomentum,Jbas) /= 0) THEN
                    itype(Imomentum,Jbas) = itype(Imomentum,Jbas)-1
                    dJI = dJI - dble(itype(Imomentum,Jbas)+1)* &
                         dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *ekinetic(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                         itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                         itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                         xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                         xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                         xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                    itype(Imomentum,Jbas) = itype(Imomentum,Jbas)+1
                 ENDIF
              ENDDO
           ENDDO
        ENDIF


        DO iCatom=1,natom
           ConMove = iCatom.eq.Iatom

           ! Check to see if they are all on Iatom.  If they are, skip this.
           ! Also check to see if none are.

           IF (IonMove .AND. JonMove .AND. ConMove) THEN
              continue
           ELSEIF ( .NOT. IonMove .AND. .NOT. JonMove .AND. .NOT. ConMove) THEN
              continue
           ELSE

              ! Note that for the moving atom we use igrad instead of itype.  This is
              ! because if Ibas=Jbas, we end up modifying both itype(Ibas) and itype(Jbas)
              ! instead of just one.  This only comes into play for moving atoms with
              ! basis functions on them.

              IF (IonMOVE) THEN
                 DO I=1,3
                    igrad(I) = itype(I,Ibas)
                 ENDDO
                 DO Icon=1,ncontract(ibas)
                    DO Jcon=1,ncontract(jbas)
                       igrad(Imomentum) = igrad(Imomentum)+1
                       dJI = dJI + 2.d0*aexp(Icon,Ibas)* &
                            dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas)* &
                            attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                            itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                            Igrad(1),Igrad(2),Igrad(3), &
                            xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                            xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                            xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)), &
                            xyz(1,iCatom),xyz(2,iCatom),xyz(3,iCatom), &
                            chg(iCatom))
                       igrad(Imomentum) = igrad(Imomentum)-1
                       IF (itype(Imomentum,Ibas) /= 0) THEN
                          igrad(Imomentum) = igrad(Imomentum)-1
                          dJI = dJI - dble(igrad(Imomentum)+1)* &
                               dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas)* &
                               attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                               itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                               Igrad(1),Igrad(2),Igrad(3), &
                               xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                               xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                               xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)), &
                               xyz(1,iCatom),xyz(2,iCatom),xyz(3,iCatom), &
                               chg(iCatom))
                          igrad(Imomentum) = igrad(Imomentum)+1
                       ENDIF
                    ENDDO
                 ENDDO
              ENDIF

              IF (JonMOVE) THEN
                 DO I=1,3
                    igrad(I) = itype(I,Jbas)
                 ENDDO
                 DO Icon=1,ncontract(ibas)
                    DO Jcon=1,ncontract(jbas)
                       igrad(Imomentum) = igrad(Imomentum)+1
                       dJI = dJI + 2.d0*aexp(Jcon,Jbas)* &
                            dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas)* &
                            attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                            Igrad(1),Igrad(2),Igrad(3), &
                            itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                            xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                            xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                            xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)), &
                            xyz(1,iCatom),xyz(2,iCatom),xyz(3,iCatom), &
                            chg(iCatom))
                       igrad(Imomentum) = igrad(Imomentum)-1
                       IF (itype(Imomentum,Jbas) /= 0) THEN
                          igrad(Imomentum) = igrad(Imomentum)-1
                          dJI = dJI - dble(igrad(Imomentum)+1)* &
                               dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas)* &
                               attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                               Igrad(1),Igrad(2),Igrad(3), &
                               itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                               xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                               xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                               xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)), &
                               xyz(1,iCatom),xyz(2,iCatom),xyz(3,iCatom), &
                               chg(iCatom))
                          igrad(Imomentum) = igrad(Imomentum)+1
                       ENDIF
                    ENDDO
                 ENDDO
              ENDIF

              ! We still have one more part to the derivative of the one center terms.
              ! If two basis functions are attracted to Iatom, the perturbation of
              ! Iatom is a non zero derivative.

              IF (ConMove) THEN
                 DO I=1,3
                    igrad(I)=0
                 ENDDO
                 igrad(imomentum)=1

                 iA = ncenter(Ibas)
                 iB = ncenter(Jbas)
                 DO Icon=1,ncontract(ibas)
                    DO Jcon=1,ncontract(jbas)
                       dJI =dJI +dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *electricfld(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                            itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                            itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                            igrad(1),igrad(2),igrad(3), &
                            xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                            xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                            xyz(1,iCatom),xyz(2,iCatom),xyz(3,iCatom),chg(iCAtom))
                    ENDDO
                 ENDDO
              ENDIF
           ENDIF
        ENDDO
        O(Jbas,Ibas)=O(Jbas,Ibas) + dJI
     ENDDO
  ENDDO

  ! The previous two terms are the one electron part of the Fock matrix.
  ! The next two terms define the two electron part.

  ! First we are going to loop over the repulsion and exchange integrals
  ! with the derivative of the density matrix.
  DO I=1,nbasis
     ! Set some variables to reduce access time for some of the more
     ! used quantities.

     xI = xyz(1,ncenter(I))
     yI = xyz(2,ncenter(I))
     zI = xyz(3,ncenter(I))
     itype1I=itype(1,I)
     itype2I=itype(2,I)
     itype3I=itype(3,I)
     DENSEII=HOLD(I,I) + HOLD2(I,I)
     DENSEIIX=HOLD(I,I)

     ! DO all the (ii|ii) integrals.
     Ibas=I
     Jbas=I
     IIbas=I
     JJbas=I
     repint=0.d0
     DO Icon=1,ncontract(ibas)
        DO Jcon=1,ncontract(jbas)
           DO IIcon=1,ncontract(iibas)
              DO JJcon=1,ncontract(jjbas)
                 repint = repint+ &
                      dcoeff(Icon,Ibas)*dcoeff(Jcon,Jbas) &
                      *dcoeff(IIcon,IIbas)*dcoeff(JJcon,JJbas)* &
                      (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                      aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                      itype1I,itype2I,itype3I,itype1I,itype2I,itype3I, &
                      itype1I,itype2I,itype3I,itype1I,itype2I,itype3I, &
                      xI,yI,zI,xI,yI,zI,xI,yI,zI,xI,yI,zI))
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     O(I,I) = O(I,I)+DENSEII*repint
     O(I,I) = O(I,I)-DENSEIIX*repint

     DO J=I+1,nbasis
        ! Set some variables to reduce access time for some of the more
        ! used quantities. (AGAIN)

        xJ = xyz(1,ncenter(J))
        yJ = xyz(2,ncenter(J))
        zJ = xyz(3,ncenter(J))
        itype1J=itype(1,J)
        itype2J=itype(2,J)
        itype3J=itype(3,J)
        DENSEJI=HOLD( J,I)+HOLD2( J,I)
        DENSEJJ=HOLD( J,J)+HOLD2( J,J)
        DENSEJIX=HOLD( J,I)
        DENSEJJX=HOLD( J,J)

        ! Find  all the (ii|jj) integrals.
        Ibas=I
        Jbas=I
        IIbas=J
        JJbas=J
        repint=0.d0
        DO Icon=1,ncontract(ibas)
           DO Jcon=1,ncontract(jbas)
              DO IIcon=1,ncontract(iibas)
                 DO JJcon=1,ncontract(jjbas)
                    repint = repint+ &
                         dcoeff(Icon,Ibas)*dcoeff(Jcon,Jbas) &
                         *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                         (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                         aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                         itype1I,itype2I,itype3I,itype1I,itype2I,itype3I, &
                         itype1J,itype2J,itype3J,itype1J,itype2J,itype3J, &
                         xI,yI,zI,xI,yI,zI,xJ,yJ,zJ,xJ,yJ,zJ))
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        O(I,I) = O(I,I)+DENSEJJ*repint
        O(J,J) = O(J,J)+DENSEII*repint
        O(J,I) = O(J,I)-DENSEJIX*repint

        ! Find  all the (ij|jj) integrals.
        Ibas=I
        Jbas=J
        IIbas=J
        JJbas=J
        repint=0.d0
        DO Icon=1,ncontract(ibas)
           DO Jcon=1,ncontract(jbas)
              DO IIcon=1,ncontract(iibas)
                 DO JJcon=1,ncontract(jjbas)
                    repint = repint+ &
                         dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                         (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                         aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                         itype1I,itype2I,itype3I,itype1J,itype2J,itype3J, &
                         itype1J,itype2J,itype3J,itype1J,itype2J,itype3J, &
                         xI,yI,zI,xJ,yJ,zJ,xJ,yJ,zJ,xJ,yJ,zJ))

                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        O(J,I) = O(J,I)+DENSEJJ*repint
        O(J,J) = O(J,J)+2.d0*DENSEJI*repint
        O(J,I) = O(J,I)-DENSEJJX*repint
        O(J,J) = O(J,J)-2.d0*DENSEJIX*repint

        ! Find  all the (ii|ij) integrals.
        Ibas=I
        Jbas=I
        iiBAS=i
        JJbas=J
        repint=0.d0
        DO Icon=1,ncontract(ibas)
           DO Jcon=1,ncontract(jbas)
              DO IIcon=1,ncontract(iibas)
                 DO JJcon=1,ncontract(jjbas)
                    repint = repint+ &
                         dcoeff(Icon,Ibas)*dcoeff(Jcon,Jbas) &
                         *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                         (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                         aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                         itype1I,itype2I,itype3I,itype1I,itype2I,itype3I, &
                         itype1I,itype2I,itype3I,itype1J,itype2J,itype3J, &
                         xI,yI,zI,xI,yI,zI,xI,yI,zI,xJ,yJ,zJ))

                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        O(J,I) = O(J,I)+DENSEII*repint
        O(I,I) = O(I,I)+2.d0*DENSEJI*repint
        O(J,I) = O(J,I)-DENSEIIX*repint
        O(I,I) = O(I,I)-2.d0*DENSEJIX*repint
        ! Find all the (ij|ij) integrals
        Ibas=I
        Jbas=J
        IIbas=I
        JJbas=J
        repint=0.d0
        DO Icon=1,ncontract(ibas)
           DO Jcon=1,ncontract(jbas)
              DO IIcon=1,ncontract(iibas)
                 DO JJcon=1,ncontract(jjbas)
                    repint = repint+ &
                         dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                         (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                         aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                         itype1I,itype2I,itype3I,itype1J,itype2J,itype3J, &
                         itype1I,itype2I,itype3I,itype1J,itype2J,itype3J, &
                         xI,yI,zI,xJ,yJ,zJ,xI,yI,zI,xJ,yJ,zJ))

                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        O(J,I) = O(J,I)+2.d0*DENSEJI*repint
        O(J,J) = O(J,J)-DENSEIIX*repint
        O(I,I) = O(I,I)-DENSEJJX*repint
        O(J,I) = O(J,I)-DENSEJIX*repint

        DO K=J+1,nbasis
           ! Set some variables to reduce access time for some of the more
           ! used quantities. (AGAIN)

           xK = xyz(1,ncenter(K))
           yK = xyz(2,ncenter(K))
           zK = xyz(3,ncenter(K))
           itype1K=itype(1,K)
           itype2K=itype(2,K)
           itype3K=itype(3,K)
           DENSEKI=HOLD(K,I)+HOLD2(K,I)
           DENSEKJ=HOLD(K,J)+HOLD2(K,J)
           DENSEKK=HOLD(K,K)+HOLD2(K,K)
           DENSEKIX=HOLD(K,I)
           DENSEKJX=HOLD(K,J)
           DENSEKKX=HOLD(K,K)

           ! Find all the (ij|ik) integrals where j>i,k>j
           Ibas=I
           Jbas=J
           IIbas=I
           JJbas=K
           repint=0.d0
           DO Icon=1,ncontract(ibas)
              DO Jcon=1,ncontract(jbas)
                 DO IIcon=1,ncontract(iibas)
                    DO JJcon=1,ncontract(jjbas)
                       repint = repint+ &
                            dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                            (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype1I,itype2I,itype3I,itype1J,itype2J,itype3J, &
                            itype1I,itype2I,itype3I,itype1K,itype2K,itype3K, &
                            xI,yI,zI,xJ,yJ,zJ,xI,yI,zI,xK,yK,zK))

                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           O(J,I) = O(J,I)+2.d0*DENSEKI*repint
           O(K,I) = O(K,I)+2.d0*DENSEJI*repint
           O(I,I) = O(I,I)-2.d0*DENSEKJX*repint
           O(J,I) = O(J,I)-DENSEKIX*repint
           O(K,I) = O(K,I)-DENSEJIX*repint
           O(K,J) = O(K,J)-DENSEIIX*repint

           ! Find all the (ij|kk) integrals where j>i, k>j.
           Ibas=I
           Jbas=J
           IIbas=K
           JJbas=K
           repint=0.d0
           DO Icon=1,ncontract(ibas)
              DO Jcon=1,ncontract(jbas)
                 DO IIcon=1,ncontract(iibas)
                    DO JJcon=1,ncontract(jjbas)
                       repint = repint+ &
                            dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                            (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype1I,itype2I,itype3I,itype1J,itype2J,itype3J, &
                            itype1K,itype2K,itype3K,itype1K,itype2K,itype3K, &
                            xI,yI,zI,xJ,yJ,zJ,xK,yK,zK,xK,yK,zK))

                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           O(J,I) = O(J,I)+DENSEKK*repint
           O(K,K) = O(K,K)+2.d0*DENSEJI*repint
           O(K,I) = O(K,I)-DENSEKJX*repint
           O(K,J) = O(K,J)-DENSEKIX*repint

           ! Find all the (ik|jj) integrals where j>i, k>j.
           Ibas=I
           Jbas=K
           IIbas=J
           JJbas=J
           repint=0.d0
           DO Icon=1,ncontract(ibas)
              DO Jcon=1,ncontract(jbas)
                 DO IIcon=1,ncontract(iibas)
                    DO JJcon=1,ncontract(jjbas)
                       repint = repint+ &
                            dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                            (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype1I,itype2I,itype3I,itype1K,itype2K,itype3K, &
                            itype1J,itype2J,itype3J,itype1J,itype2J,itype3J, &
                            xI,yI,zI,xK,yK,zK,xJ,yJ,zJ,xJ,yJ,zJ))

                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           O(K,I) = O(K,I)+DENSEJJ*repint
           O(J,J) = O(J,J)+2.d0*DENSEKI*repint
           O(K,J) = O(K,J)-DENSEJIX*repint
           O(J,I) = O(J,I)-DENSEKJX*repint

           ! Find all the (ii|jk) integrals where j>i, k>j.
           Ibas=I
           Jbas=I
           IIbas=J
           JJbas=K
           repint=0.d0
           DO Icon=1,ncontract(ibas)
              DO Jcon=1,ncontract(jbas)
                 DO IIcon=1,ncontract(iibas)
                    DO JJcon=1,ncontract(jjbas)
                       repint = repint+ &
                            dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                            (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype1I,itype2I,itype3I,itype1I,itype2I,itype3I, &
                            itype1J,itype2J,itype3J,itype1K,itype2K,itype3K, &
                            xI,yI,zI,xI,yI,zI,xJ,yJ,zJ,xK,yK,zK))

                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           O(K,J) = O(K,J)+DENSEII*repint
           O(I,I) = O(I,I)+2.d0*DENSEKJ*repint
           O(J,I) = O(J,I)-DENSEKIX*repint
           O(K,I) = O(K,I)-DENSEJIX*repint
        ENDDO

        DO K=I+1,nbasis-1
           xK = xyz(1,ncenter(K))
           yK = xyz(2,ncenter(K))
           zK = xyz(3,ncenter(K))
           itype1K=itype(1,K)
           itype2K=itype(2,K)
           itype3K=itype(3,K)
           DENSEKI=HOLD( K,I)+HOLD2( K,I)
           DENSEKJ=HOLD( K,J)+HOLD2( K,J)
           DENSEKK=HOLD( K,K)+HOLD2( K,K)
           DENSEKIX=HOLD( K,I)
           DENSEKJX=HOLD( K,J)
           DENSEKKX=HOLD( K,K)

           DO L=K+1,nbasis
              xL = xyz(1,ncenter(L))
              yL = xyz(2,ncenter(L))
              zL = xyz(3,ncenter(L))
              itype1L=itype(1,L)
              itype2L=itype(2,L)
              itype3L=itype(3,L)
              DENSELJ=HOLD( L,J)+HOLD2( L,J)
              DENSELI=HOLD( L,I)+HOLD2( L,I)
              DENSELK=HOLD( L,K)+HOLD2( L,K)
              DENSELJX=HOLD( L,J)
              DENSELIX=HOLD( L,I)
              DENSELKX=HOLD( L,K)

              ! Find the (ij|kl) integrals where j>i,k>i,l>k. Note that k and j
              ! can be equal.

              Ibas=I
              Jbas=J
              IIbas=K
              JJbas=L
              repint=0.d0
              DO Icon=1,ncontract(ibas)
                 DO Jcon=1,ncontract(jbas)
                    DO IIcon=1,ncontract(iibas)
                       DO JJcon=1,ncontract(jjbas)
                          repint = repint+ &
                               dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                               *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                               (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                               aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                               itype1I,itype2I,itype3I,itype1J,itype2J,itype3J, &
                               itype1K,itype2K,itype3K,itype1L,itype2L,itype3L, &
                               xI,yI,zI,xJ,yJ,zJ,xK,yK,zK,xL,yL,zL))


                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
              O(J,I) = O(J,I)+2.d0*DENSELK*repint
              O(L,K) = O(L,K)+2.d0*DENSEJI*repint
              O(K,I) = O(K,I)-DENSELJX*repint
              O(L,I) = O(L,I)-DENSEKJX*repint
              IF (J == K) THEN
                 O(J,K) = O(J,K)-2.d0*DENSELIX*repint
                 O(J,L) = O(J,L)-DENSEKIX*repint
                 O(L,J) = O(L,J)-DENSEKIX*repint
              ELSEIF (J == L) THEN
                 O(J,L) = O(J,L)-2.d0*DENSEKIX*repint
                 O(J,K) = O(J,K)-DENSELIX*repint
                 O(K,J) = O(K,J)-DENSELIX*repint
              ELSE
                 O(J,K) = O(J,K)-DENSELIX*repint
                 O(J,L) = O(J,L)-DENSEKIX*repint
                 O(K,J) = O(K,J)-DENSELIX*repint
                 O(L,J) = O(L,J)-DENSEKIX*repint
              ENDIF
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  ! Now we are going to loop over the derivatives of the repulsion and exchange
  ! integrals with the regular density matrix.

  DO I=1,nbasis
     ! Set some variables to reduce access time for some of the more
     ! used quantities.

     xI = xyz(1,ncenter(I))
     yI = xyz(2,ncenter(I))
     zI = xyz(3,ncenter(I))
     itype1I=itype(1,I)
     itype2I=itype(2,I)
     itype3I=itype(3,I)
     DENSEII=DENSE(I,I) + DENSEB(I,I)
     DENSEIIX=DENSE(I,I)

     ! Skip the (ii|ii) integralas, as they move with the core.

     DO J=I+1,nbasis
        ! Set some variables to reduce access time for some of the more
        ! used quantities. (AGAIN)

        xJ = xyz(1,ncenter(J))
        yJ = xyz(2,ncenter(J))
        zJ = xyz(3,ncenter(J))
        itype1J=itype(1,J)
        itype2J=itype(2,J)
        itype3J=itype(3,J)
        DENSEJI=DENSE(J,I)+DENSEB(J,I)
        DENSEJJ=DENSE(J,J)+DENSEB(J,J)
        DENSEJIX=DENSE(J,I)
        DENSEJJX=DENSE(J,J)

        ! Find  all the (ii|jj) integrals.
        Ibas=I
        Jbas=I
        IIbas=J
        JJbas=J
        repint=0.d0
        call move1twoe(I,I,J,J,Iatom,Imomentum,repint)
        O(I,I) = O(I,I)+DENSEJJ*repint
        O(J,J) = O(J,J)+DENSEII*repint
        O(J,I) = O(J,I)-DENSEJIX*repint

        ! Find  all the (ij|jj) integrals.
        Ibas=I
        Jbas=J
        IIbas=J
        JJbas=J
        repint=0.d0
        call move1twoe(I,J,J,J,Iatom,Imomentum,repint)
        O(J,I) = O(J,I)+DENSEJJ*repint
        O(J,J) = O(J,J)+2.d0*DENSEJI*repint
        O(J,I) = O(J,I)-DENSEJJX*repint
        O(J,J) = O(J,J)-2.d0*DENSEJIX*repint

        ! Find  all the (ii|ij) integrals.
        Ibas=I
        Jbas=I
        iiBAS=i
        JJbas=J
        repint=0.d0
        call move1twoe(I,I,I,J,Iatom,Imomentum,repint)
        O(J,I) = O(J,I)+DENSEII*repint
        O(I,I) = O(I,I)+2.d0*DENSEJI*repint
        O(J,I) = O(J,I)-DENSEIIX*repint
        O(I,I) = O(I,I)-2.d0*DENSEJIX*repint

        ! Find all the (ij|ij) integrals
        Ibas=I
        Jbas=J
        IIbas=I
        JJbas=J
        repint=0.d0
        call move1twoe(I,J,I,J,Iatom,Imomentum,repint)
        O(J,I) = O(J,I)+2.d0*DENSEJI*repint
        O(J,J) = O(J,J)-DENSEIIX*repint
        O(I,I) = O(I,I)-DENSEJJX*repint
        O(J,I) = O(J,I)-DENSEJIX*repint

        DO K=J+1,nbasis
           ! Set some variables to reduce access time for some of the more
           ! used quantities. (AGAIN)

           xK = xyz(1,ncenter(K))
           yK = xyz(2,ncenter(K))
           zK = xyz(3,ncenter(K))
           itype1K=itype(1,K)
           itype2K=itype(2,K)
           itype3K=itype(3,K)
           DENSEKI=DENSE(K,I)+DENSEB(K,I)
           DENSEKJ=DENSE(K,J)+DENSEB(K,J)
           DENSEKK=DENSE(K,K)+DENSEB(K,K)
           DENSEKIX=DENSE(K,I)
           DENSEKJX=DENSE(K,J)
           DENSEKKX=DENSE(K,K)

           ! Find all the (ij|ik) integrals where j>i,k>j
           Ibas=I
           Jbas=J
           IIbas=I
           JJbas=K
           repint=0.d0
           call move1twoe(I,J,I,K,Iatom,Imomentum,repint)
           O(J,I) = O(J,I)+2.d0*DENSEKI*repint
           O(K,I) = O(K,I)+2.d0*DENSEJI*repint
           O(I,I) = O(I,I)-2.d0*DENSEKJX*repint
           O(J,I) = O(J,I)-DENSEKIX*repint
           O(K,I) = O(K,I)-DENSEJIX*repint
           O(K,J) = O(K,J)-DENSEIIX*repint

           ! Find all the (ij|kk) integrals where j>i, k>j.
           Ibas=I
           Jbas=J
           IIbas=K
           JJbas=K
           repint=0.d0
           call move1twoe(I,J,K,K,Iatom,Imomentum,repint)
           O(J,I) = O(J,I)+DENSEKK*repint
           O(K,K) = O(K,K)+2.d0*DENSEJI*repint
           O(K,I) = O(K,I)-DENSEKJX*repint
           O(K,J) = O(K,J)-DENSEKIX*repint

           ! Find all the (ik|jj) integrals where j>i, k>j.
           Ibas=I
           Jbas=K
           IIbas=J
           JJbas=J
           repint=0.d0
           call move1twoe(I,K,J,J,Iatom,Imomentum,repint)
           O(K,I) = O(K,I)+DENSEJJ*repint
           O(J,J) = O(J,J)+2.d0*DENSEKI*repint
           O(K,J) = O(K,J)-DENSEJIX*repint
           O(J,I) = O(J,I)-DENSEKJX*repint

           ! Find all the (ii|jk) integrals where j>i, k>j.
           Ibas=I
           Jbas=I
           IIbas=J
           JJbas=K
           repint=0.d0
           call move1twoe(I,I,J,K,Iatom,Imomentum,repint)
           O(K,J) = O(K,J)+DENSEII*repint
           O(I,I) = O(I,I)+2.d0*DENSEKJ*repint
           O(J,I) = O(J,I)-DENSEKIX*repint
           O(K,I) = O(K,I)-DENSEJIX*repint
        ENDDO

        DO K=I+1,nbasis-1
           xK = xyz(1,ncenter(K))
           yK = xyz(2,ncenter(K))
           zK = xyz(3,ncenter(K))
           itype1K=itype(1,K)
           itype2K=itype(2,K)
           itype3K=itype(3,K)
           DENSEKI=DENSE(K,I)+DENSEB(K,I)
           DENSEKJ=DENSE(K,J)+DENSEB(K,J)
           DENSEKK=DENSE(K,K)+DENSEB(K,K)
           DENSEKIX=DENSE(K,I)
           DENSEKJX=DENSE(K,J)
           DENSEKKX=DENSE(K,K)

           DO L=K+1,nbasis
              xL = xyz(1,ncenter(L))
              yL = xyz(2,ncenter(L))
              zL = xyz(3,ncenter(L))
              itype1L=itype(1,L)
              itype2L=itype(2,L)
              itype3L=itype(3,L)
              DENSELJ=DENSE(L,J)+DENSEB(L,J)
              DENSELI=DENSE(L,I)+DENSEB(L,I)
              DENSELK=DENSE(L,K)+DENSEB(L,K)
              DENSELJX=DENSE(L,J)
              DENSELIX=DENSE(L,I)
              DENSELKX=DENSE(L,K)

              ! Find the (ij|kl) integrals where j>i,k>i,l>k. Note that k and j
              ! can be equal.

              Ibas=I
              Jbas=J
              IIbas=K
              JJbas=L
              repint=0.d0
              call move1twoe(I,J,K,L,Iatom,Imomentum,repint)
              O(J,I) = O(J,I)+2.d0*DENSELK*repint
              O(L,K) = O(L,K)+2.d0*DENSEJI*repint
              O(K,I) = O(K,I)-DENSELJX*repint
              O(L,I) = O(L,I)-DENSEKJX*repint
              IF (J == K) THEN
                 O(J,K) = O(J,K)-2.d0*DENSELIX*repint
                 O(J,L) = O(J,L)-DENSEKIX*repint
                 O(L,J) = O(L,J)-DENSEKIX*repint
              ELSEIF (J == L) THEN
                 O(J,L) = O(J,L)-2.d0*DENSEKIX*repint
                 O(J,K) = O(J,K)-DENSELIX*repint
                 O(K,J) = O(K,J)-DENSELIX*repint
              ELSE
                 O(J,K) = O(J,K)-DENSELIX*repint
                 O(J,L) = O(J,L)-DENSEKIX*repint
                 O(K,J) = O(K,J)-DENSELIX*repint
                 O(L,J) = O(L,J)-DENSEKIX*repint
              ENDIF
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  DO Ibas=1,nbasis
     DO Jbas=Ibas+1,nbasis
        O(Ibas,Jbas) = O(Jbas,Ibas)
     ENDDO
  ENDDO

END subroutine duhfoperatora


! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Ed Brothers. December 16, 2002
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

subroutine duhfoperatorb(IDX)
  use allmod
  implicit double precision(a-h,o-z)
  dimension igrad(3)
  logical :: IonMove, JonMove, ConMove

  ! The purpose of this subroutine is to form the first derivative
  ! of the beta fock matrix for use in forming the first derivative
  ! of the energy weighted density matrix.


  ! Note that this matrix is symmetric.

  ! Blank out the array.

  DO I=1,nbasis
     DO J=1,nbasis
        O(J,I)=0.d0
     ENDDO
  ENDDO

  ! The first thing to do is find the atom and direction of the perturbation.

  Imomentum=mod(idx,3)
  IF (Imomentum == 0) THEN
     Iatom=IDX/3
     Imomentum = 3
  ELSE
     Iatom = (IDX-Imomentum)/3 + 1
  ENDIF

  ! We need to take the derivative of the one electron portion of the
  ! Hamiltonian.

  DO Ibas=1,nbasis
     IonMove = ncenter(Ibas).eq.Iatom
     DO Jbas = Ibas,nbasis
        JonMove = ncenter(Jbas).eq.Iatom

        ! If I is on Iatom (the moving atom) and J is not, the derivative  of
        ! the kinetic energy is non-zero.

        dJI=0.d0
        IF (IonMove .AND. .NOT. Jonmove) THEN
           DO Icon=1,ncontract(ibas)
              DO Jcon=1,ncontract(jbas)
                 itype(Imomentum,Ibas) = itype(Imomentum,Ibas)+1
                 dJI = dJI + 2.d0*aexp(Icon,Ibas)* &
                      dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                      *ekinetic(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                      itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                      itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                      xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                      xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                      xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                 itype(Imomentum,Ibas) = itype(Imomentum,Ibas)-1
                 IF (itype(Imomentum,Ibas) /= 0) THEN
                    itype(Imomentum,Ibas) = itype(Imomentum,Ibas)-1
                    dJI = dJI - dble(itype(Imomentum,Ibas)+1)* &
                         dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *ekinetic(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                         itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                         itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                         xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                         xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                         xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                    itype(Imomentum,Ibas) = itype(Imomentum,Ibas)+1
                 ENDIF
              ENDDO
           ENDDO
        ENDIF

        ! If J is on Iatom (the moving atom) and I is not, the derivative  of
        ! the kinetic energy is non-zero.

        IF (JonMove .AND. .NOT. Ionmove) THEN
           DO Icon=1,ncontract(ibas)
              DO Jcon=1,ncontract(jbas)
                 itype(Imomentum,Jbas) = itype(Imomentum,Jbas)+1
                 dJI = dJI + 2.d0*aexp(Jcon,Jbas)* &
                      dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                      *ekinetic(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                      itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                      itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                      xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                      xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                      xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                 itype(Imomentum,Jbas) = itype(Imomentum,Jbas)-1
                 IF (itype(Imomentum,Jbas) /= 0) THEN
                    itype(Imomentum,Jbas) = itype(Imomentum,Jbas)-1
                    dJI = dJI - dble(itype(Imomentum,Jbas)+1)* &
                         dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *ekinetic(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                         itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                         itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                         xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                         xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                         xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
                    itype(Imomentum,Jbas) = itype(Imomentum,Jbas)+1
                 ENDIF
              ENDDO
           ENDDO
        ENDIF


        DO iCatom=1,natom
           ConMove = iCatom.eq.Iatom

           ! Check to see if they are all on Iatom.  If they are, skip this.
           ! Also check to see if none are.

           IF (IonMove .AND. JonMove .AND. ConMove) THEN
              continue
           ELSEIF ( .NOT. IonMove .AND. .NOT. JonMove .AND. .NOT. ConMove) THEN
              continue
           ELSE
              ! Note that for the moving atom we use igrad instead of itype.  This is
              ! because if Ibas=Jbas, we end up modifying both itype(Ibas) and itype(Jbas)
              ! instead of just one.  This only comes into play for moving atoms with
              ! basis functions on them.

              IF (IonMOVE) THEN
                 DO I=1,3
                    igrad(I) = itype(I,Ibas)
                 ENDDO
                 DO Icon=1,ncontract(ibas)
                    DO Jcon=1,ncontract(jbas)
                       igrad(Imomentum) = igrad(Imomentum)+1
                       dJI = dJI + 2.d0*aexp(Icon,Ibas)* &
                            dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas)* &
                            attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                            itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                            Igrad(1),Igrad(2),Igrad(3), &
                            xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                            xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                            xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)), &
                            xyz(1,iCatom),xyz(2,iCatom),xyz(3,iCatom), &
                            chg(iCatom))
                       igrad(Imomentum) = igrad(Imomentum)-1
                       IF (itype(Imomentum,Ibas) /= 0) THEN
                          igrad(Imomentum) = igrad(Imomentum)-1
                          dJI = dJI - dble(igrad(Imomentum)+1)* &
                               dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas)* &
                               attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                               itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                               Igrad(1),Igrad(2),Igrad(3), &
                               xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                               xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                               xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)), &
                               xyz(1,iCatom),xyz(2,iCatom),xyz(3,iCatom), &
                               chg(iCatom))
                          igrad(Imomentum) = igrad(Imomentum)+1
                       ENDIF
                    ENDDO
                 ENDDO
              ENDIF

              IF (JonMOVE) THEN
                 DO I=1,3
                    igrad(I) = itype(I,Jbas)
                 ENDDO
                 DO Icon=1,ncontract(ibas)
                    DO Jcon=1,ncontract(jbas)
                       igrad(Imomentum) = igrad(Imomentum)+1
                       dJI = dJI + 2.d0*aexp(Jcon,Jbas)* &
                            dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas)* &
                            attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                            Igrad(1),Igrad(2),Igrad(3), &
                            itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                            xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                            xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                            xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)), &
                            xyz(1,iCatom),xyz(2,iCatom),xyz(3,iCatom), &
                            chg(iCatom))
                       igrad(Imomentum) = igrad(Imomentum)-1
                       IF (itype(Imomentum,Jbas) /= 0) THEN
                          igrad(Imomentum) = igrad(Imomentum)-1
                          dJI = dJI - dble(igrad(Imomentum)+1)* &
                               dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas)* &
                               attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                               Igrad(1),Igrad(2),Igrad(3), &
                               itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                               xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                               xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                               xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)), &
                               xyz(1,iCatom),xyz(2,iCatom),xyz(3,iCatom), &
                               chg(iCatom))
                          igrad(Imomentum) = igrad(Imomentum)+1
                       ENDIF
                    ENDDO
                 ENDDO
              ENDIF

              ! We still have one more part to the derivative of the one center terms.
              ! If two basis functions are attracted to Iatom, the perturbation of
              ! Iatom is a non zero derivative.

              IF (ConMove) THEN
                 DO I=1,3
                    igrad(I)=0
                 ENDDO
                 igrad(imomentum)=1

                 iA = ncenter(Ibas)
                 iB = ncenter(Jbas)
                 DO Icon=1,ncontract(ibas)
                    DO Jcon=1,ncontract(jbas)
                       dJI =dJI +dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *electricfld(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                            itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                            itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                            igrad(1),igrad(2),igrad(3), &
                            xyz(1,iB),xyz(2,iB),xyz(3,iB), &
                            xyz(1,iA),xyz(2,iA),xyz(3,iA), &
                            xyz(1,iCatom),xyz(2,iCatom),xyz(3,iCatom),chg(iCAtom))
                    ENDDO
                 ENDDO
              ENDIF
           ENDIF
        ENDDO
        O(Jbas,Ibas)=O(Jbas,Ibas) + dJI
     ENDDO
  ENDDO

  ! The previous two terms are the one electron part of the Fock matrix.
  ! The next two terms define the two electron part.

  DO I=1,nbasis

     ! Set some variables to reduce access time for some of the more
     ! used quantities.

     xI = xyz(1,ncenter(I))
     yI = xyz(2,ncenter(I))
     zI = xyz(3,ncenter(I))
     itype1I=itype(1,I)
     itype2I=itype(2,I)
     itype3I=itype(3,I)
     DENSEII=HOLD( I,I) + HOLD2( I,I)
     DENSEIIX=HOLD2( I,I)

     ! DO all the (ii|ii) integrals.
     Ibas=I
     Jbas=I
     IIbas=I
     JJbas=I
     repint=0.d0
     DO Icon=1,ncontract(ibas)
        DO Jcon=1,ncontract(jbas)
           DO IIcon=1,ncontract(iibas)
              DO JJcon=1,ncontract(jjbas)
                 repint = repint+ &
                      dcoeff(Icon,Ibas)*dcoeff(Jcon,Jbas) &
                      *dcoeff(IIcon,IIbas)*dcoeff(JJcon,JJbas)* &
                      (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                      aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                      itype1I,itype2I,itype3I,itype1I,itype2I,itype3I, &
                      itype1I,itype2I,itype3I,itype1I,itype2I,itype3I, &
                      xI,yI,zI,xI,yI,zI,xI,yI,zI,xI,yI,zI))
              ENDDO
           ENDDO
        ENDDO
     ENDDO
     O(I,I) = O(I,I)+DENSEII*repint
     O(I,I) = O(I,I)-DENSEIIX*repint

     DO J=I+1,nbasis
        ! Set some variables to reduce access time for some of the more
        ! used quantities. (AGAIN)

        xJ = xyz(1,ncenter(J))
        yJ = xyz(2,ncenter(J))
        zJ = xyz(3,ncenter(J))
        itype1J=itype(1,J)
        itype2J=itype(2,J)
        itype3J=itype(3,J)
        DENSEJI=HOLD( J,I)+HOLD2( J,I)
        DENSEJJ=HOLD( J,J)+HOLD2( J,J)
        DENSEJIX=HOLD2( J,I)
        DENSEJJX=HOLD2( J,J)

        ! Find  all the (ii|jj) integrals.
        Ibas=I
        Jbas=I
        IIbas=J
        JJbas=J
        repint=0.d0
        DO Icon=1,ncontract(ibas)
           DO Jcon=1,ncontract(jbas)
              DO IIcon=1,ncontract(iibas)
                 DO JJcon=1,ncontract(jjbas)
                    repint = repint+ &
                         dcoeff(Icon,Ibas)*dcoeff(Jcon,Jbas) &
                         *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                         (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                         aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                         itype1I,itype2I,itype3I,itype1I,itype2I,itype3I, &
                         itype1J,itype2J,itype3J,itype1J,itype2J,itype3J, &
                         xI,yI,zI,xI,yI,zI,xJ,yJ,zJ,xJ,yJ,zJ))
                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        O(I,I) = O(I,I)+DENSEJJ*repint
        O(J,J) = O(J,J)+DENSEII*repint
        O(J,I) = O(J,I)-DENSEJIX*repint

        ! Find  all the (ij|jj) integrals.
        Ibas=I
        Jbas=J
        IIbas=J
        JJbas=J
        repint=0.d0
        DO Icon=1,ncontract(ibas)
           DO Jcon=1,ncontract(jbas)
              DO IIcon=1,ncontract(iibas)
                 DO JJcon=1,ncontract(jjbas)
                    repint = repint+ &
                         dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                         (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                         aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                         itype1I,itype2I,itype3I,itype1J,itype2J,itype3J, &
                         itype1J,itype2J,itype3J,itype1J,itype2J,itype3J, &
                         xI,yI,zI,xJ,yJ,zJ,xJ,yJ,zJ,xJ,yJ,zJ))

                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        O(J,I) = O(J,I)+DENSEJJ*repint
        O(J,J) = O(J,J)+2.d0*DENSEJI*repint
        O(J,I) = O(J,I)-DENSEJJX*repint
        O(J,J) = O(J,J)-2.d0*DENSEJIX*repint

        ! Find  all the (ii|ij) integrals.
        Ibas=I
        Jbas=I
        iiBAS=i
        JJbas=J
        repint=0.d0
        DO Icon=1,ncontract(ibas)
           DO Jcon=1,ncontract(jbas)
              DO IIcon=1,ncontract(iibas)
                 DO JJcon=1,ncontract(jjbas)
                    repint = repint+ &
                         dcoeff(Icon,Ibas)*dcoeff(Jcon,Jbas) &
                         *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                         (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                         aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                         itype1I,itype2I,itype3I,itype1I,itype2I,itype3I, &
                         itype1I,itype2I,itype3I,itype1J,itype2J,itype3J, &
                         xI,yI,zI,xI,yI,zI,xI,yI,zI,xJ,yJ,zJ))

                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        O(J,I) = O(J,I)+DENSEII*repint
        O(I,I) = O(I,I)+2.d0*DENSEJI*repint
        O(J,I) = O(J,I)-DENSEIIX*repint
        O(I,I) = O(I,I)-2.d0*DENSEJIX*repint
        ! Find all the (ij|ij) integrals
        Ibas=I
        Jbas=J
        IIbas=I
        JJbas=J
        repint=0.d0
        DO Icon=1,ncontract(ibas)
           DO Jcon=1,ncontract(jbas)
              DO IIcon=1,ncontract(iibas)
                 DO JJcon=1,ncontract(jjbas)
                    repint = repint+ &
                         dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                         (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                         aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                         itype1I,itype2I,itype3I,itype1J,itype2J,itype3J, &
                         itype1I,itype2I,itype3I,itype1J,itype2J,itype3J, &
                         xI,yI,zI,xJ,yJ,zJ,xI,yI,zI,xJ,yJ,zJ))

                 ENDDO
              ENDDO
           ENDDO
        ENDDO
        O(J,I) = O(J,I)+2.d0*DENSEJI*repint
        O(J,J) = O(J,J)-DENSEIIX*repint
        O(I,I) = O(I,I)-DENSEJJX*repint
        O(J,I) = O(J,I)-DENSEJIX*repint

        DO K=J+1,nbasis
           ! Set some variables to reduce access time for some of the more
           ! used quantities. (AGAIN)

           xK = xyz(1,ncenter(K))
           yK = xyz(2,ncenter(K))
           zK = xyz(3,ncenter(K))
           itype1K=itype(1,K)
           itype2K=itype(2,K)
           itype3K=itype(3,K)
           DENSEKI=HOLD( K,I)+HOLD2( K,I)
           DENSEKJ=HOLD( K,J)+HOLD2( K,J)
           DENSEKK=HOLD( K,K)+HOLD2( K,K)
           DENSEKIX=HOLD2( K,I)
           DENSEKJX=HOLD2( K,J)
           DENSEKKX=HOLD2( K,K)

           ! Find all the (ij|ik) integrals where j>i,k>j
           Ibas=I
           Jbas=J
           IIbas=I
           JJbas=K
           repint=0.d0
           DO Icon=1,ncontract(ibas)
              DO Jcon=1,ncontract(jbas)
                 DO IIcon=1,ncontract(iibas)
                    DO JJcon=1,ncontract(jjbas)
                       repint = repint+ &
                            dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                            (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype1I,itype2I,itype3I,itype1J,itype2J,itype3J, &
                            itype1I,itype2I,itype3I,itype1K,itype2K,itype3K, &
                            xI,yI,zI,xJ,yJ,zJ,xI,yI,zI,xK,yK,zK))

                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           O(J,I) = O(J,I)+2.d0*DENSEKI*repint
           O(K,I) = O(K,I)+2.d0*DENSEJI*repint
           O(I,I) = O(I,I)-2.d0*DENSEKJX*repint
           O(J,I) = O(J,I)-DENSEKIX*repint
           O(K,I) = O(K,I)-DENSEJIX*repint
           O(K,J) = O(K,J)-DENSEIIX*repint

           ! Find all the (ij|kk) integrals where j>i, k>j.
           Ibas=I
           Jbas=J
           IIbas=K
           JJbas=K
           repint=0.d0
           DO Icon=1,ncontract(ibas)
              DO Jcon=1,ncontract(jbas)
                 DO IIcon=1,ncontract(iibas)
                    DO JJcon=1,ncontract(jjbas)
                       repint = repint+ &
                            dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                            (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype1I,itype2I,itype3I,itype1J,itype2J,itype3J, &
                            itype1K,itype2K,itype3K,itype1K,itype2K,itype3K, &
                            xI,yI,zI,xJ,yJ,zJ,xK,yK,zK,xK,yK,zK))

                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           O(J,I) = O(J,I)+DENSEKK*repint
           O(K,K) = O(K,K)+2.d0*DENSEJI*repint
           O(K,I) = O(K,I)-DENSEKJX*repint
           O(K,J) = O(K,J)-DENSEKIX*repint

           ! Find all the (ik|jj) integrals where j>i, k>j.
           Ibas=I
           Jbas=K
           IIbas=J
           JJbas=J
           repint=0.d0
           DO Icon=1,ncontract(ibas)
              DO Jcon=1,ncontract(jbas)
                 DO IIcon=1,ncontract(iibas)
                    DO JJcon=1,ncontract(jjbas)
                       repint = repint+ &
                            dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                            (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype1I,itype2I,itype3I,itype1K,itype2K,itype3K, &
                            itype1J,itype2J,itype3J,itype1J,itype2J,itype3J, &
                            xI,yI,zI,xK,yK,zK,xJ,yJ,zJ,xJ,yJ,zJ))

                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           O(K,I) = O(K,I)+DENSEJJ*repint
           O(J,J) = O(J,J)+2.d0*DENSEKI*repint
           O(K,J) = O(K,J)-DENSEJIX*repint
           O(J,I) = O(J,I)-DENSEKJX*repint

           ! Find all the (ii|jk) integrals where j>i, k>j.
           Ibas=I
           Jbas=I
           IIbas=J
           JJbas=K
           repint=0.d0
           DO Icon=1,ncontract(ibas)
              DO Jcon=1,ncontract(jbas)
                 DO IIcon=1,ncontract(iibas)
                    DO JJcon=1,ncontract(jjbas)
                       repint = repint+ &
                            dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                            (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype1I,itype2I,itype3I,itype1I,itype2I,itype3I, &
                            itype1J,itype2J,itype3J,itype1K,itype2K,itype3K, &
                            xI,yI,zI,xI,yI,zI,xJ,yJ,zJ,xK,yK,zK))

                    ENDDO
                 ENDDO
              ENDDO
           ENDDO
           O(K,J) = O(K,J)+DENSEII*repint
           O(I,I) = O(I,I)+2.d0*DENSEKJ*repint
           O(J,I) = O(J,I)-DENSEKIX*repint
           O(K,I) = O(K,I)-DENSEJIX*repint
        ENDDO

        DO K=I+1,nbasis-1
           xK = xyz(1,ncenter(K))
           yK = xyz(2,ncenter(K))
           zK = xyz(3,ncenter(K))
           itype1K=itype(1,K)
           itype2K=itype(2,K)
           itype3K=itype(3,K)
           DENSEKI=HOLD( K,I)+HOLD2( K,I)
           DENSEKJ=HOLD( K,J)+HOLD2( K,J)
           DENSEKK=HOLD( K,K)+HOLD2( K,K)
           DENSEKIX=HOLD2( K,I)
           DENSEKJX=HOLD2( K,J)
           DENSEKKX=HOLD2( K,K)

           DO L=K+1,nbasis
              xL = xyz(1,ncenter(L))
              yL = xyz(2,ncenter(L))
              zL = xyz(3,ncenter(L))
              itype1L=itype(1,L)
              itype2L=itype(2,L)
              itype3L=itype(3,L)
              DENSELJ=HOLD( L,J)+HOLD2( L,J)
              DENSELI=HOLD( L,I)+HOLD2( L,I)
              DENSELK=HOLD( L,K)+HOLD2( L,K)
              DENSELJX=HOLD2( L,J)
              DENSELIX=HOLD2( L,I)
              DENSELKX=HOLD2( L,K)

              ! Find the (ij|kl) integrals where j>i,k>i,l>k. Note that k and j
              ! can be equal.

              Ibas=I
              Jbas=J
              IIbas=K
              JJbas=L
              repint=0.d0
              DO Icon=1,ncontract(ibas)
                 DO Jcon=1,ncontract(jbas)
                    DO IIcon=1,ncontract(iibas)
                       DO JJcon=1,ncontract(jjbas)
                          repint = repint+ &
                               dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                               *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                               (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                               aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                               itype1I,itype2I,itype3I,itype1J,itype2J,itype3J, &
                               itype1K,itype2K,itype3K,itype1L,itype2L,itype3L, &
                               xI,yI,zI,xJ,yJ,zJ,xK,yK,zK,xL,yL,zL))


                       ENDDO
                    ENDDO
                 ENDDO
              ENDDO
              O(J,I) = O(J,I)+2.d0*DENSELK*repint
              O(L,K) = O(L,K)+2.d0*DENSEJI*repint
              O(K,I) = O(K,I)-DENSELJX*repint
              O(L,I) = O(L,I)-DENSEKJX*repint
              IF (J == K) THEN
                 O(J,K) = O(J,K)-2.d0*DENSELIX*repint
                 O(J,L) = O(J,L)-DENSEKIX*repint
                 O(L,J) = O(L,J)-DENSEKIX*repint
              ELSEIF (J == L) THEN
                 O(J,L) = O(J,L)-2.d0*DENSEKIX*repint
                 O(J,K) = O(J,K)-DENSELIX*repint
                 O(K,J) = O(K,J)-DENSELIX*repint
              ELSE
                 O(J,K) = O(J,K)-DENSELIX*repint
                 O(J,L) = O(J,L)-DENSEKIX*repint
                 O(K,J) = O(K,J)-DENSELIX*repint
                 O(L,J) = O(L,J)-DENSEKIX*repint
              ENDIF
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  ! Now we are going to loop over the derivatives of the repulsion and exchange
  ! integrals with the regular density matrix.

  DO I=1,nbasis
     ! Set some variables to reduce access time for some of the more
     ! used quantities.

     xI = xyz(1,ncenter(I))
     yI = xyz(2,ncenter(I))
     zI = xyz(3,ncenter(I))
     itype1I=itype(1,I)
     itype2I=itype(2,I)
     itype3I=itype(3,I)
     DENSEII=DENSE(I,I) + DENSEB(I,I)
     DENSEIIX=DENSEB(I,I)

     ! Skip the (ii|ii) integrals, as they move with the core.

     DO J=I+1,nbasis
        ! Set some variables to reduce access time for some of the more
        ! used quantities. (AGAIN)

        xJ = xyz(1,ncenter(J))
        yJ = xyz(2,ncenter(J))
        zJ = xyz(3,ncenter(J))
        itype1J=itype(1,J)
        itype2J=itype(2,J)
        itype3J=itype(3,J)
        DENSEJI=DENSE(J,I)+DENSEB(J,I)
        DENSEJJ=DENSE(J,J)+DENSEB(J,J)
        DENSEJIX=DENSEB(J,I)
        DENSEJJX=DENSEB(J,J)

        ! Find  all the (ii|jj) integrals.
        Ibas=I
        Jbas=I
        IIbas=J
        JJbas=J
        repint=0.d0
        call move1twoe(I,I,J,J,Iatom,Imomentum,repint)
        O(I,I) = O(I,I)+DENSEJJ*repint
        O(J,J) = O(J,J)+DENSEII*repint
        O(J,I) = O(J,I)-DENSEJIX*repint

        ! Find  all the (ij|jj) integrals.
        Ibas=I
        Jbas=J
        IIbas=J
        JJbas=J
        repint=0.d0
        call move1twoe(I,J,J,J,Iatom,Imomentum,repint)
        O(J,I) = O(J,I)+DENSEJJ*repint
        O(J,J) = O(J,J)+2.d0*DENSEJI*repint
        O(J,I) = O(J,I)-DENSEJJX*repint
        O(J,J) = O(J,J)-2.d0*DENSEJIX*repint

        ! Find  all the (ii|ij) integrals.
        Ibas=I
        Jbas=I
        iiBAS=i
        JJbas=J
        repint=0.d0
        call move1twoe(I,I,I,J,Iatom,Imomentum,repint)
        O(J,I) = O(J,I)+DENSEII*repint
        O(I,I) = O(I,I)+2.d0*DENSEJI*repint
        O(J,I) = O(J,I)-DENSEIIX*repint
        O(I,I) = O(I,I)-2.d0*DENSEJIX*repint
        ! Find all the (ij|ij) integrals
        Ibas=I
        Jbas=J
        IIbas=I
        JJbas=J
        repint=0.d0
        call move1twoe(I,J,I,J,Iatom,Imomentum,repint)
        O(J,I) = O(J,I)+2.d0*DENSEJI*repint
        O(J,J) = O(J,J)-DENSEIIX*repint
        O(I,I) = O(I,I)-DENSEJJX*repint
        O(J,I) = O(J,I)-DENSEJIX*repint

        DO K=J+1,nbasis
           ! Set some variables to reduce access time for some of the more
           ! used quantities. (AGAIN)

           xK = xyz(1,ncenter(K))
           yK = xyz(2,ncenter(K))
           zK = xyz(3,ncenter(K))
           itype1K=itype(1,K)
           itype2K=itype(2,K)
           itype3K=itype(3,K)
           DENSEKI=DENSE(K,I)+DENSEB(K,I)
           DENSEKJ=DENSE(K,J)+DENSEB(K,J)
           DENSEKK=DENSE(K,K)+DENSEB(K,K)
           DENSEKIX=DENSEB(K,I)
           DENSEKJX=DENSEB(K,J)
           DENSEKKX=DENSEB(K,K)

           ! Find all the (ij|ik) integrals where j>i,k>j
           Ibas=I
           Jbas=J
           IIbas=I
           JJbas=K
           repint=0.d0
           call move1twoe(I,J,I,K,Iatom,Imomentum,repint)
           O(J,I) = O(J,I)+2.d0*DENSEKI*repint
           O(K,I) = O(K,I)+2.d0*DENSEJI*repint
           O(I,I) = O(I,I)-2.d0*DENSEKJX*repint
           O(J,I) = O(J,I)-DENSEKIX*repint
           O(K,I) = O(K,I)-DENSEJIX*repint
           O(K,J) = O(K,J)-DENSEIIX*repint

           ! Find all the (ij|kk) integrals where j>i, k>j.
           Ibas=I
           Jbas=J
           IIbas=K
           JJbas=K
           repint=0.d0
           call move1twoe(I,J,K,K,Iatom,Imomentum,repint)
           O(J,I) = O(J,I)+DENSEKK*repint
           O(K,K) = O(K,K)+2.d0*DENSEJI*repint
           O(K,I) = O(K,I)-DENSEKJX*repint
           O(K,J) = O(K,J)-DENSEKIX*repint

           ! Find all the (ik|jj) integrals where j>i, k>j.
           Ibas=I
           Jbas=K
           IIbas=J
           JJbas=J
           repint=0.d0
           call move1twoe(I,K,J,J,Iatom,Imomentum,repint)
           O(K,I) = O(K,I)+DENSEJJ*repint
           O(J,J) = O(J,J)+2.d0*DENSEKI*repint
           O(K,J) = O(K,J)-DENSEJIX*repint
           O(J,I) = O(J,I)-DENSEKJX*repint

           ! Find all the (ii|jk) integrals where j>i, k>j.
           Ibas=I
           Jbas=I
           IIbas=J
           JJbas=K
           repint=0.d0
           call move1twoe(I,I,J,K,Iatom,Imomentum,repint)
           O(K,J) = O(K,J)+DENSEII*repint
           O(I,I) = O(I,I)+2.d0*DENSEKJ*repint
           O(J,I) = O(J,I)-DENSEKIX*repint
           O(K,I) = O(K,I)-DENSEJIX*repint
        ENDDO

        DO K=I+1,nbasis-1
           xK = xyz(1,ncenter(K))
           yK = xyz(2,ncenter(K))
           zK = xyz(3,ncenter(K))
           itype1K=itype(1,K)
           itype2K=itype(2,K)
           itype3K=itype(3,K)
           DENSEKI=DENSE(K,I)+DENSEB(K,I)
           DENSEKJ=DENSE(K,J)+DENSEB(K,J)
           DENSEKK=DENSE(K,K)+DENSEB(K,K)
           DENSEKIX=DENSEB(K,I)
           DENSEKJX=DENSEB(K,J)
           DENSEKKX=DENSEB(K,K)

           DO L=K+1,nbasis
              xL = xyz(1,ncenter(L))
              yL = xyz(2,ncenter(L))
              zL = xyz(3,ncenter(L))
              itype1L=itype(1,L)
              itype2L=itype(2,L)
              itype3L=itype(3,L)
              DENSELJ=DENSE(L,J)+DENSEB(L,J)
              DENSELI=DENSE(L,I)+DENSEB(L,I)
              DENSELK=DENSE(L,K)+DENSEB(L,K)
              DENSELJX=DENSEB(L,J)
              DENSELIX=DENSEB(L,I)
              DENSELKX=DENSEB(L,K)

              ! Find the (ij|kl) integrals where j>i,k>i,l>k. Note that k and j
              ! can be equal.

              Ibas=I
              Jbas=J
              IIbas=K
              JJbas=L
              repint=0.d0
              call move1twoe(I,J,K,L,Iatom,Imomentum,repint)
              O(J,I) = O(J,I)+2.d0*DENSELK*repint
              O(L,K) = O(L,K)+2.d0*DENSEJI*repint
              O(K,I) = O(K,I)-DENSELJX*repint
              O(L,I) = O(L,I)-DENSEKJX*repint
              IF (J == K) THEN
                 O(J,K) = O(J,K)-2.d0*DENSELIX*repint
                 O(J,L) = O(J,L)-DENSEKIX*repint
                 O(L,J) = O(L,J)-DENSEKIX*repint
              ELSEIF (J == L) THEN
                 O(J,L) = O(J,L)-2.d0*DENSEKIX*repint
                 O(J,K) = O(J,K)-DENSELIX*repint
                 O(K,J) = O(K,J)-DENSELIX*repint
              ELSE
                 O(J,K) = O(J,K)-DENSELIX*repint
                 O(J,L) = O(J,L)-DENSEKIX*repint
                 O(K,J) = O(K,J)-DENSELIX*repint
                 O(L,J) = O(L,J)-DENSEKIX*repint
              ENDIF
           ENDDO
        ENDDO
     ENDDO
  ENDDO
  DO Ibas=1,nbasis
     DO Jbas=Ibas+1,nbasis
        O(Ibas,Jbas) = O(Jbas,Ibas)
     ENDDO
  ENDDO

END subroutine duhfoperatorb



! Ed Brothers. December 16, 2002.
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

subroutine move1twoe(Ibas,Jbas,IIbas,JJbas,Iatom,Imomentum,repint)
  use allmod
  implicit double precision(a-h,o-z)

  dimension itype2(3,4)
  logical :: same,doit
  same = .false.
  doit = .false.

  ! This calculates the derivative of 4 center integral (I J | II JJ) with
  ! respect to moving center Iatom, and stores the value in repint.


  ! Find the centers the basis functions are located on.  If all the
  ! functions are on the same center, return as this is a zero result.

  iA = ncenter(Ibas)
  iB = ncenter(Jbas)
  iC = ncenter(IIbas)
  iD = ncenter(JJbas)

  same = iA.eq.iB
  same = same .and. iB.eq.iC
  same = same .and. iC.eq.iD

  IF (same) return

  ! Now check to see if any of the functions are on center iatom.  If none
  ! are, return.

  doit = iA.eq.Iatom
  doit = doit .or. iB.eq.Iatom
  doit = doit .or. iC.eq.Iatom
  doit = doit .or. iD.eq.Iatom

  IF ( .NOT. doit) return

  Agrad=0.d0
  Bgrad=0.d0
  Cgrad=0.d0
  Dgrad=0.d0

  ! The itype2 array was added because if Ibas=Jbas, the code raises two
  ! angular momentums instead of one.

  DO Icopy=1,3
     itype2(Icopy,1) = itype(Icopy,Ibas)
     itype2(Icopy,2) = itype(Icopy,Jbas)
     itype2(Icopy,3) = itype(Icopy,IIbas)
     itype2(Icopy,4) = itype(Icopy,JJbas)
  ENDDO

  ! This is again a modification of grad2elec, but we are only looking at
  ! mone atom and one direction.  Check each atom to see if we need it.

  DO Icon = 1, ncontract(Ibas)
     DO Jcon = 1, ncontract(Jbas)
        DO IIcon = 1, ncontract(IIbas)
           DO JJcon = 1, ncontract(JJbas)
              cntrctcoeff = dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                   *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)

              doit = iA.eq.iatom

              IF (doit) THEN
                 itype2(Imomentum,1) = itype2(Imomentum,1)+1
                 Agrad = Agrad+2.d0*aexp(Icon,Ibas)*cntrctcoeff* &
                      repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                      aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                      itype2(1,1), itype2(2,1), itype2(3,1), &
                      itype2(1,2), itype2(2,2), itype2(3,2), &
                      itype2(1,3), itype2(2,3), itype2(3,3), &
                      itype2(1,4), itype2(2,4), itype2(3,4), &
                      xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                      xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                      xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                      xyz(1,iD), xyz(2,iD),xyz(3,iD))
                 itype2(Imomentum,1) = itype2(Imomentum,1)-1

                 IF (itype2(Imomentum,1) /= 0) THEN
                    itype2(Imomentum,1) = itype2(Imomentum,1)-1
                    temp = cntrctcoeff* &
                         repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                         aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                         itype2(1,1), itype2(2,1), itype2(3,1), &
                         itype2(1,2), itype2(2,2), itype2(3,2), &
                         itype2(1,3), itype2(2,3), itype2(3,3), &
                         itype2(1,4), itype2(2,4), itype2(3,4), &
                         xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                         xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                         xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                         xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    itype2(Imomentum,1) = itype2(Imomentum,1)+1
                    Agrad = Agrad-dble(itype2(Imomentum,1))*temp
                 ENDIF
              ENDIF

              doit = iB.eq.iatom

              IF (doit) THEN
                 itype2(Imomentum,2) = itype2(Imomentum,2)+1
                 Bgrad = Bgrad+2.d0*aexp(Jcon,Jbas)*cntrctcoeff* &
                      repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                      aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                      itype2(1,1), itype2(2,1), itype2(3,1), &
                      itype2(1,2), itype2(2,2), itype2(3,2), &
                      itype2(1,3), itype2(2,3), itype2(3,3), &
                      itype2(1,4), itype2(2,4), itype2(3,4), &
                      xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                      xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                      xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                      xyz(1,iD), xyz(2,iD),xyz(3,iD))
                 itype2(Imomentum,2) = itype2(Imomentum,2)-1
                 IF (itype2(Imomentum,2) /= 0) THEN
                    itype2(Imomentum,2) = itype2(Imomentum,2)-1
                    temp = cntrctcoeff* &
                         repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                         aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                         itype2(1,1), itype2(2,1), itype2(3,1), &
                         itype2(1,2), itype2(2,2), itype2(3,2), &
                         itype2(1,3), itype2(2,3), itype2(3,3), &
                         itype2(1,4), itype2(2,4), itype2(3,4), &
                         xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                         xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                         xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                         xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    itype2(Imomentum,2) = itype2(Imomentum,2)+1
                    Bgrad = Bgrad-dble(itype2(Imomentum,2))*temp
                 ENDIF
              ENDIF

              doit = iC.eq.iAtom

              IF (doit) THEN

                 itype2(Imomentum,3) = itype2(Imomentum,3)+1
                 Cgrad = Cgrad+2.d0*aexp(IIcon,IIbas)*cntrctcoeff* &
                      repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                      aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                      itype2(1,1), itype2(2,1), itype2(3,1), &
                      itype2(1,2), itype2(2,2), itype2(3,2), &
                      itype2(1,3), itype2(2,3), itype2(3,3), &
                      itype2(1,4), itype2(2,4), itype2(3,4), &
                      xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                      xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                      xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                      xyz(1,iD), xyz(2,iD),xyz(3,iD))
                 itype2(Imomentum,3) = itype2(Imomentum,3)-1
                 IF (itype2(Imomentum,3) /= 0) THEN
                    itype2(Imomentum,3) = itype2(Imomentum,3)-1
                    temp = cntrctcoeff* &
                         repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                         aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                         itype2(1,1), itype2(2,1), itype2(3,1), &
                         itype2(1,2), itype2(2,2), itype2(3,2), &
                         itype2(1,3), itype2(2,3), itype2(3,3), &
                         itype2(1,4), itype2(2,4), itype2(3,4), &
                         xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                         xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                         xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                         xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    itype2(Imomentum,3) = itype2(Imomentum,3)+1
                    Cgrad = Cgrad-dble(itype2(Imomentum,3))*temp
                 ENDIF
              ENDIF

              doit = iD.eq.iatom

              IF (doit) THEN

                 itype2(Imomentum,4) = itype2(Imomentum,4)+1
                 Dgrad = Dgrad+2.d0*aexp(JJcon,JJbas)*cntrctcoeff* &
                      repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                      aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                      itype2(1,1), itype2(2,1), itype2(3,1), &
                      itype2(1,2), itype2(2,2), itype2(3,2), &
                      itype2(1,3), itype2(2,3), itype2(3,3), &
                      itype2(1,4), itype2(2,4), itype2(3,4), &
                      xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                      xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                      xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                      xyz(1,iD), xyz(2,iD),xyz(3,iD))
                 itype2(Imomentum,4) = itype2(Imomentum,4)-1
                 IF (itype2(Imomentum,4) /= 0) THEN
                    itype2(Imomentum,4) = itype2(Imomentum,4)-1
                    temp = cntrctcoeff* &
                         repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                         aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                         itype2(1,1), itype2(2,1), itype2(3,1), &
                         itype2(1,2), itype2(2,2), itype2(3,2), &
                         itype2(1,3), itype2(2,3), itype2(3,3), &
                         itype2(1,4), itype2(2,4), itype2(3,4), &
                         xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                         xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                         xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                         xyz(1,iD), xyz(2,iD),xyz(3,iD))
                    itype2(Imomentum,4) = itype2(Imomentum,4)+1
                    Dgrad = Dgrad-dble(itype2(Imomentum,4))*temp
                 ENDIF
              ENDIF
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  repint = Agrad+Bgrad+Cgrad+Dgrad


end subroutine move1twoe



! Ed Brothers. May 29, 2002.
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

double precision recursive function elctfldrecurse(i,j,k,ii,jj,kk, &
     idx,idy,idz,m,aux,Ax,Ay,Az,Bx,By,Bz, &
     Cx,Cy,Cz,Px,Py,Pz,g) &
     result(elctfldrec)
  implicit double precision(a-h,o-z)
  dimension iexponents(6),center(12),aux(0:20)

  ! The this is taken from the recursive relation found in Obara and Saika,
  ! J. Chem. Phys. 84 (7) 1986, 3963.

  ! Check to see if the integral has become just the nuclear attraction
  ! integral.

  IF (idx+idy+idz == 0) THEN
     elctfldrec = attrecurse(i,j,k,ii,jj,kk,m,aux,Ax,Ay, &
          Az,Bx,By,Bz,Cx,Cy,Cz,Px,Py,Pz,g)

     ! If it hasn't, check to see if it has become a simple integral over s
     ! functions.

  ELSEIF (i+j+k+ii+jj+kk == 0) THEN
     IF (idx == 2) THEN
        elctfldrec = -2.d0*g*attrecurse(i,j,k,ii,jj,kk, &
             m+1,aux,Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Px,Py,Pz,g) &
             +4.d0*g*g*(Px-Cx)*(Px-Cx)*attrecurse(i,j,k,ii,jj,kk, &
             m+2,aux,Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Px,Py,Pz,g)
     ELSEIF (idy == 2) THEN
        elctfldrec = -2.d0*g*attrecurse(i,j,k,ii,jj,kk, &
             m+1,aux,Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Px,Py,Pz,g) &
             +4.d0*g*g*(Py-Cy)*(Py-Cy)*attrecurse(i,j,k,ii,jj,kk, &
             m+2,aux,Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Px,Py,Pz,g)
     ELSEIF (idz == 2) THEN
        elctfldrec = -2.d0*g*attrecurse(i,j,k,ii,jj,kk, &
             m+1,aux,Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Px,Py,Pz,g) &
             +4.d0*g*g*(Pz-Cz)*(Pz-Cz)*attrecurse(i,j,k,ii,jj,kk, &
             m+2,aux,Ax,Ay,Az,Bx,By,Bz,Cx,Cy,Cz,Px,Py,Pz,g)
     ELSEIF (idx == 1 .AND. idy == 1) THEN
        elctfldrec = 4.d0*g*g*(Px-Cx)*(Py-Cy)* &
             attrecurse(i,j,k,ii,jj,kk,m+2,aux,Ax,Ay,Az,Bx,By,Bz, &
             Cx,Cy,Cz,Px,Py,Pz,g)
     ELSEIF (idx == 1 .AND. idz == 1) THEN
        elctfldrec = 4.d0*g*g*(Px-Cx)*(Pz-Cz)* &
             attrecurse(i,j,k,ii,jj,kk,m+2,aux,Ax,Ay,Az,Bx,By,Bz, &
             Cx,Cy,Cz,Px,Py,Pz,g)
     ELSEIF (idy == 1 .AND. idz == 1) THEN
        elctfldrec = 4.d0*g*g*(Py-Cy)*(Pz-Cz)* &
             attrecurse(i,j,k,ii,jj,kk,m+2,aux,Ax,Ay,Az,Bx,By,Bz, &
             Cx,Cy,Cz,Px,Py,Pz,g)
     ELSEIF (idx == 1) THEN
        elctfldrec = 2.d0*g*(Px-Cx)* &
             attrecurse(i,j,k,ii,jj,kk,m+1,aux,Ax,Ay,Az,Bx,By,Bz, &
             Cx,Cy,Cz,Px,Py,Pz,g)
     ELSEIF (idy == 1) THEN
        elctfldrec = 2.d0*g*(Py-Cy)* &
             attrecurse(i,j,k,ii,jj,kk,m+1,aux,Ax,Ay,Az,Bx,By,Bz, &
             Cx,Cy,Cz,Px,Py,Pz,g)
     ELSEIF (idz == 1) THEN
        elctfldrec = 2.d0*g*(Pz-Cz)* &
             attrecurse(i,j,k,ii,jj,kk,m+1,aux,Ax,Ay,Az,Bx,By,Bz, &
             Cx,Cy,Cz,Px,Py,Pz,g)
     ENDIF

     ! Otherwise, use the recusion relation from Obara and Saika.  The first
     ! step is to find the lowest nonzero angular momentum exponent.  This is
     ! because the more exponents equal zero the fewer terms need to be
     ! calculated, and each recursive loop reduces the angular momentum
     ! exponents. This therefore reorders the atoms and sets the exponent
     ! to be reduced.

  ELSE
     iexponents(1) = i
     iexponents(2) = j
     iexponents(3) = k
     iexponents(4) = ii
     iexponents(5) = jj
     iexponents(6) = kk
     center(7) = Cx
     center(8) = Cy
     center(9) = Cz
     center(10)= Px
     center(11)= Py
     center(12)= Pz
     ilownum=300
     ilowex=300
     DO L=1,6
        IF (iexponents(L) < ilowex .AND. iexponents(L) /= 0) THEN
           ilowex=iexponents(L)
           ilownum=L
        ENDIF
     ENDDO
     IF (ilownum <= 3) THEN
        center(1)=Ax
        center(2)=Ay
        center(3)=Az
        center(4)=Bx
        center(5)=By
        center(6)=Bz
     ELSE
        center(4)=Ax
        center(5)=Ay
        center(6)=Az
        center(1)=Bx
        center(2)=By
        center(3)=Bz
        iexponents(4) = i
        iexponents(5) = j
        iexponents(6) = k
        iexponents(1) = ii
        iexponents(2) = jj
        iexponents(3) = kk
        ilownum = ilownum - 3
     ENDIF

     ! The first step is lowering the orbital exponent by one.

     iexponents(ilownum) = iexponents(ilownum)-1

     ! At this point, calculate the first two terms of the recusion
     ! relation.

     elctfldrec = 0.d0
     PA = center(9+ilownum)-center(ilownum)
     IF (PA /= 0) elctfldrec  = elctfldrec  + PA * &
          elctfldrecurse(iexponents(1),iexponents(2), &
          iexponents(3),iexponents(4), &
          iexponents(5),iexponents(6), &
          idx,idy,idz,m,aux, &
          center(1),center(2),center(3), &
          center(4),center(5),center(6), &
          center(7),center(8),center(9), &
          center(10),center(11),center(12),g)

     PC = center(9+ilownum)-center(6+ilownum)
     IF (PC /= 0) elctfldrec  = elctfldrec  - PC * &
          elctfldrecurse(iexponents(1),iexponents(2), &
          iexponents(3),iexponents(4), &
          iexponents(5),iexponents(6), &
          idx,idy,idz,m+1,aux, &
          center(1),center(2),center(3), &
          center(4),center(5),center(6), &
          center(7),center(8),center(9), &
          center(10),center(11),center(12),g)

     ! The next two terms only arise is the angual momentum of the dimension
     ! of A that has already been lowered is not zero.  In other words, if a
     ! (px|1/rc|px) was passed to this subroutine, we are now considering
     ! (s|1/rc|px), and the following term does not arise, as the x expoent
     ! on A is zero.

     IF (iexponents(ilownum) /= 0) THEN
        coeff = dble(iexponents(ilownum))/(2.d0*g)
        iexponents(ilownum) = iexponents(ilownum)-1
        elctfldrec  = elctfldrec  + coeff*( &
             elctfldrecurse(iexponents(1),iexponents(2), &
             iexponents(3),iexponents(4), &
             iexponents(5),iexponents(6), &
             idx,idy,idz,m,aux, &
             center(1),center(2),center(3), &
             center(4),center(5),center(6), &
             center(7),center(8),center(9), &
             center(10),center(11),center(12),g) &
             -elctfldrecurse(iexponents(1),iexponents(2), &
             iexponents(3),iexponents(4), &
             iexponents(5),iexponents(6), &
             idx,idy,idz,m+1,aux, &
             center(1),center(2),center(3), &
             center(4),center(5),center(6), &
             center(7),center(8),center(9), &
             center(10),center(11),center(12),g) &
             )
        iexponents(ilownum) = iexponents(ilownum)+1
     ENDIF

     ! The next two terms only arise is the angual momentum of the dimension
     ! of A that has already been lowered is not zero in B.  If a
     ! (px|1/rc|px) was passed to this subroutine, we are now considering
     ! (s|1/rc|px), and the following term does arise, as the x exponent on
     ! B is 1.

     IF (iexponents(ilownum+3) /= 0) THEN
        coeff = dble(iexponents(ilownum+3))/(2.d0*g)
        iexponents(ilownum+3) = iexponents(ilownum+3)-1
        elctfldrec = elctfldrec + coeff*( &
             elctfldrecurse(iexponents(1),iexponents(2), &
             iexponents(3),iexponents(4), &
             iexponents(5),iexponents(6), &
             idx,idy,idz,m,aux, &
             center(1),center(2),center(3), &
             center(4),center(5),center(6), &
             center(7),center(8),center(9), &
             center(10),center(11),center(12),g) &
             -elctfldrecurse(iexponents(1),iexponents(2), &
             iexponents(3),iexponents(4), &
             iexponents(5),iexponents(6), &
             idx,idy,idz,m+1,aux, &
             center(1),center(2),center(3), &
             center(4),center(5),center(6), &
             center(7),center(8),center(9), &
             center(10),center(11),center(12),g) &
             )
        iexponents(ilownum+3) = iexponents(ilownum+3)+1
     ENDIF

     ! Finally there is a lowering of the derivative term, which only occurs
     ! if the angular momentum being lowered in this step corresponds to the
     ! derivative direction.

     IF (ilownum == 1 .AND. idx > 0) THEN
        elctfldrec = elctfldrec + dble(idx)* &
             elctfldrecurse(iexponents(1),iexponents(2), &
             iexponents(3),iexponents(4), &
             iexponents(5),iexponents(6), &
             idx-1,idy,idz,m+1,aux, &
             center(1),center(2),center(3), &
             center(4),center(5),center(6), &
             center(7),center(8),center(9), &
             center(10),center(11),center(12),g)
     ELSEIF (ilownum == 2 .AND. idy > 0) THEN
        elctfldrec = elctfldrec + dble(idy)* &
             elctfldrecurse(iexponents(1),iexponents(2), &
             iexponents(3),iexponents(4), &
             iexponents(5),iexponents(6), &
             idx,idy-1,idz,m+1,aux, &
             center(1),center(2),center(3), &
             center(4),center(5),center(6), &
             center(7),center(8),center(9), &
             center(10),center(11),center(12),g)
     ELSEIF (ilownum == 3 .AND. idz > 0) THEN
        elctfldrec = elctfldrec + dble(idz)* &
             elctfldrecurse(iexponents(1),iexponents(2), &
             iexponents(3),iexponents(4), &
             iexponents(5),iexponents(6), &
             idx,idy,idz-1,m+1,aux, &
             center(1),center(2),center(3), &
             center(4),center(5),center(6), &
             center(7),center(8),center(9), &
             center(10),center(11),center(12),g)
     ENDIF
  ENDIF


  return
end function elctfldrecurse




! S. Dixon's diagonalization code from DivCon.
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

    SUBROUTINE hessDIAG(NDIM,A,NEVEC1,TOLERA,V,EVAL1,IDEGEN1,EVEC1, &
    IERROR)

! DRIVER ROUTINE FOR DIAGONALIZATION OF THE REAL, SYMMETRIC,
! MATRIX A.

! VARIABLES REQUIRED:
! ------------------

! NDIM = ORDER OF THE MATRIX A (I.E., A IS NDIM BY NDIM);

! A = REAL SYMMETRIC MATRIX TO BE DIAGONALIZED.  ONLY THE LOWER
! HALF OF A NEED BE FILLED WHEN CALLING.  A IS DESTROYED BY
! THIS ROUTINE.

! NEVEC1 = THE NUMBER OF EIGENVECTORS REQUIRED;

! TOLERA = TOLERANCE FACTOR USED IN THE QR ITERATION TO DETERMINE
! WHEN OFF-DIAGONAL ENTRIES ARE ESSENTIALLY ZERO.  (DEFAULT
! IS 1.0D-8).

! V = 3 BY NDIM WORKSPACE.


! VARIABLES RETURNED:
! ------------------

! EVAL1 = EIGENVALUES OF A (SORTED IN INCREASING ALGEBRAIC VALUE);

! IDEGEN1 = DEGENERACIES (NUMBER OF TIMES REPEATED) FOR EIGENVALUES;

! EVEC1 = EIGENVECTORS OF A (IN COLUMNS OF EVEC1);

! ERROR CODES:  IERROR=0 - SUCCESSFUL CALCULATION.
! IERROR=1 - NO CONVERGENCE IN QR ITERATION.


! PROGRAMMED BY S. L. DIXON, OCT., 1991.


    use allmod
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
! DIMENSION A(NDIM,*),V(3,*),EVAL1(*),IDEGEN1(*),EVEC1(NDIM,*)
    DIMENSION A(natom*3,natom*3),V(3,natom*3),EVAL1(natom*3), &
    IDEGEN1(natom*3),EVEC1(natom*3,natom*3)

! FLAG FOR WHETHER OR NOT TO COMPUTE EIGENVALUES:


    IF(NDIM == 1)THEN
        EVAL1(1) = A(1,1)
        IDEGEN1(1) = 1
        EVEC1(1,1) = 1.0D0
        RETURN
    ENDIF

! TRIDIAGONALIZE THE MATRIX A.  THIS WILL OVERWRITE THE DIAGONAL
! AND SUBDIAGONAL OF A WITH THE TRIDIAGONALIZED VERSION.  THE
! HOUSEHOLDER VECTORS ARE RETURNED IN THE ROWS ABOVE THE DIAGONAL,
! AND THE BETAHS ARE RETURNED BELOW THE SUBDIAGONAL.

    CALL hessTRIDI(NDIM,V,A)

! COMPUTE NORM OF TRIDIAGONAL MATRIX FROM THE "LARGEST" COLUMN.

    ANORM = ABS(A(1,1)) + ABS(A(2,1))
    IF(NDIM > 2)THEN
        DO 20 I=2,NDIM-1
            AICOL = ABS(A(I-1,I)) + ABS(A(I,I)) + ABS(A(I+1,I))
            ANORM = MAX(ANORM,AICOL)
        20 ENDDO
    ENDIF
    ANCOL = ABS(A(NDIM-1,NDIM)) + ABS(A(NDIM,NDIM))
    ANORM = MAX(ANORM,ANCOL)

! GET EIGENVALUES AND DEGENERACIES OF THE TRIDIAGONAL MATRIX A.
! IF THE CALLING ROUTINE HAS NOT SUPPLIED A TOLERANCE FACTOR FOR
! OFF-DIAGONAL ENTRIES IN THE QR ITERATION, A DEFAULT OF 1.0D-8
! WILL BE USED.

    TOLTMP = TOLERA
    IF(TOLTMP <= 0.0D0) TOLTMP = 1.0D-8
    CALL hessEIGVAL(NDIM,A,V,TOLTMP,ANORM,EVAL1,IERROR)
    IF(IERROR /= 0) RETURN

! DETERMINE DEGENERACIES OF EIGENVALUES.

    CALL hessDEGEN(NDIM,EVAL1,TOLTMP,ANORM,IDEGEN1)

! GET EIGENVECTORS OF TRIDIAGONALIZED VERSION OF A.

    IF(NEVEC1 <= 0) RETURN
    CALL hessEIGVEC(NDIM,NEVEC1,A,V,TOLTMP,ANORM,EVAL1,IDEGEN1,EVEC1)

! PREMULTIPLY EVEC1 BY THE HOUSEHOLDER MATRIX USED TO TRIDIAGONALIZE
! A.  THIS TRANSFORMS EIGENVECTORS OF THE TRIDIAGONAL MATRIX TO
! THOSE OF THE ORIGINAL MATRIX A.  SEE SUBROUTINE TRIDI FOR
! STORAGE OF HOUSEHOLDER TRANSFORMATION.

    IF(NDIM > 2)THEN
        DO 30 K=1,NDIM-2
            V(1,K) = A(K+2,K)
        30 ENDDO
    
    ! SWAP STORAGE SO THAT THE EXPENSIVE TRIPLE LOOP BELOW DOESN'T
    ! HAVE TO JUMP ACROSS COLUMNS TO GET ENTRIES OF A.
    
        DO 50 I=2,NDIM
            DO 40 J=1,I-1
                A(I,J) = A(J,I)
            40 ENDDO
        50 ENDDO
        DO 160 J=1,NEVEC1
            DO 140 M=2,NDIM-1
                K = NDIM - M
                BETAH = V(1,K)
                IF(ABS(BETAH) < 1.0D-50) CONTINUE !GO TO 140
                SUM = 0.0D0
                DO 100 I=K+1,NDIM
                    SUM = SUM + A(I,K)*EVEC1(I,J)
                100 ENDDO
                BSUM = BETAH*SUM
                DO 120 I=K+1,NDIM
                    EVEC1(I,J) = EVEC1(I,J) - A(I,K)*BSUM
                120 ENDDO
            140 ENDDO
        160 ENDDO
    ENDIF
    RETURN
    end SUBROUTINE hessDIAG



    SUBROUTINE hessTRIDI(NDIM,V,A)

! TRIDIAGONALIZES A REAL, SYMMETRIC MATRIX A BY THE METHOD OF
! HOUSEHOLDER (J. H. WILKINSON, THE COMPUTER JOURNAL, VOL. 3,
! P. 23 (1960)).  NDIM IS THE ORDER OF A.  THE DIAGONAL AND
! SUBDIAGONAL OF A ARE OVERWRITTEN WITH THE TRIDIAGONALIZED
! VERSION OF A.  THE VECTORS USED IN EACH HOUSEHOLDER
! TRANSFORMATION ARE STORED ABOVE THE DIAGONAL IN THE FIRST
! NDIM-2 ROWS OF A.  THE BETAHS ARE RETURNED BELOW THE SUBDIAGONAL
! OF A.  V IS A WORKSPACE ARRAY.

! PROGRAMMED BY S. L. DIXON, OCT., 1991.


    use allmod
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
! DIMENSION A(NDIM,*),V(3,*)
    DIMENSION A(natom*3,natom*3),V(3,natom*3)

! THRESH WILL BE USED AS A THRESHOLD TO DETERMINE IF A VALUE SHOULD
! BE CONSIDERED TO BE ZERO.  THIS CAN BE CHANGED BY THE USER.

    THRESH = 1.0D-50

! IF A IS 2 BY 2 OR SMALLER, THEN IT IS ALREADY TRIDIAGONAL -- NO
! NEED TO CONTINUE.

    IF(NDIM <= 2) GO TO 1000
    DO 500 K=1,NDIM-2
    
    ! DETERMINE THE VECTOR V USED IN THE HOUSEHOLDER TRANSFORMATION P.
    ! FOR EACH VALUE OF K THE HOUSEHOLDER MATRIX P IS DEFINED AS:
    
    ! P = I - BETAH*V*V'
    
    
    ! CONSTRUCT A HOUSEHOLDER TRANSFORMATION ONLY IF THERE IS A NONZERO
    ! OFF-DIAGONAL ELEMENT BELOW A(K,K).
    
        ALPHA2 = 0.0D0
        DO 60 I=K+1,NDIM
            V(1,I) = A(I,K)
            ALPHA2 = ALPHA2 + V(1,I)**2
        60 ENDDO
        APTEMP = ALPHA2 - V(1,K+1)**2
        ALPHA = DSQRT(ALPHA2)
        IF(ALPHA >= THRESH)THEN
            BETAH = 1.0D0/(ALPHA*(ALPHA + ABS(V(1,K+1))))
            SGN = SIGN(1.0D0,V(1,K+1))
            V(1,K+1) = V(1,K+1) + SGN*ALPHA
        
        ! NOW OVERWRITE A WITH P'*A*P.  THE ENTRIES BELOW THE SUBDIAGONAL
        ! IN THE KTH COLUMN ARE ZEROED BY THE PREMULTIPLICATION BY P'.
        ! THESE ENTRIES WILL BE LEFT ALONE TO SAVE TIME.
        
            AKV = APTEMP + A(K+1,K)*V(1,K+1)
            S = BETAH*AKV
            A(K+1,K) = A(K+1,K) - S*V(1,K+1)
        
        ! NOW THE SUBMATRIX CONSISTING OF ROWS K+1,NDIM AND COLUMNS K+1,NDIM
        ! MUST BE OVERWRITTEN WITH THE TRANSFORMATION.
        
            DOT12 = 0.0D0
            BHALF = BETAH*0.5D0
            DO 220 I=K+1,NDIM
                SUM = 0.0D0
                DO 100 J=K+1,I
                    SUM = SUM + A(I,J)*V(1,J)
                100 ENDDO
                IF(I < NDIM)THEN
                    DO 180 J=I+1,NDIM
                    
                    ! AN UPPER TRIANGULAR ENTRY OF A WILL BE REQUIRED.  MUST USE
                    ! THE SYMMETRIC ENTRY IN THE LOWER TRIANGULAR PART OF A.
                    
                        SUM = SUM + A(J,I)*V(1,J)
                    180 ENDDO
                ENDIF
                V(2,I) = BETAH*SUM
                DOT12 = DOT12 + V(1,I)*V(2,I)
            220 ENDDO
            BH12 = BHALF*DOT12
            DO 300 I=K+1,NDIM
                V(2,I) = V(2,I) - BH12*V(1,I)
            300 ENDDO
            DO 350 J=K+1,NDIM
                DO 310 I=J,NDIM
                    A(I,J) = A(I,J) - V(1,I)*V(2,J) - V(2,I)*V(1,J)
                310 ENDDO
            
            ! STORE V(1,J) ABOVE THE DIAGONAL IN ROW K OF A
            
                A(K,J) = V(1,J)
            350 ENDDO
        
        ! STORE BETAH BELOW THE SUBDIAGONAL OF A.
        
            A(K+2,K) = BETAH
        ELSE
        
        ! NO HOUSEHOLDER TRANSFORMATION IS NECESSARY BECAUSE THE OFF-
        ! DIAGONALS ARE ALL ESSENTIALLY ZERO.
        
            A(K+2,K) = 0.0D0
            DO 460 J=K+1,NDIM
                A(K,J) = 0.0D0
            460 ENDDO
        ENDIF
    500 ENDDO
    1000 RETURN
    end SUBROUTINE hessTRIDI



    SUBROUTINE hessEIGVAL(NDIM,A,BETAH,TOLERA,ANORM,EVAL1,IERROR)

! QR ROUTINE FOR THE DETERMINATION OF ALL THE EIGENVALUES
! OF THE NDIM BY NDIM SYMMETRIC, TRIDIAGONAL MATRIX A.

! INPUT:

! NDIM   = SIZE OF MATRIX A.
! A      = NDIM BY NDIM SYMMETRIC TRIDIAGONAL MATRIX.
! BETAH   = 3 BY NDIM WORKSPACE.
! TOLERA = SMALL NUMBER USED TO DETERMINE WHEN OFF-DIAGONAL
! ELEMENTS ARE ESSENTIALLY ZERO.
! ANORM  = ABSOLUTE COLUMN NORM OF TRIDIAGONAL MATRIX A.


! RETURNED:

! EVAL1   = EIGENVALUES OF A IN ASCENDING ORDER.
! IERROR = 1 IF QR ITERATION DID NOT CONVERGE; 0 OTHERWISE.

! PROGRAMMED BY S. L. DIXON.


    use allmod
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
! DIMENSION A(NDIM,*),BETAH(3,*),EVAL1(*)
    DIMENSION A(natom*3,natom*3),BETAH(3,natom*3),EVAL1(natom*3)
    IERROR = 0
! ITMAX = 20
    ITMAX = 200

! TOLERANCE FOR OFF-DIAGONAL ELEMENTS:

    EPSLON = TOLERA*ANORM

! COPY DIAGONAL ELEMENTS OF A TO EVAL1, AND SUBDIAGONAL ELEMENTS
! TO BETAH.

    EVAL1(1) = A(1,1)
    BETAH(1,1) = A(2,1)
    IF(NDIM > 2)THEN
        DO 50 I=2,NDIM-1
            EVAL1(I) = A(I,I)
            BETAH(1,I) = A(I+1,I)
        50 ENDDO
    ENDIF
    EVAL1(NDIM) = A(NDIM,NDIM)

! EACH QR ITERATION WILL OPERATE ON THE UNREDUCED TRIDIAGONAL
! SUBMATRIX WITH UPPER LEFT ELEMENT (L,L) AND LOWER RIGHT ELEMENT
! (N,N).

    L = 1
    N = NDIM
    ITER = 0

! FIND THE SMALLEST UNREDUCED SUBMATRIX WITH LOWER RIGHT CORNER AT
! (N,N).  I.E., SEARCH UPWARD FOR A BETAH THAT IS ZERO.

    80 KUPPER = N-L
    DO 100 K=1,KUPPER
        I = N-K
        IF(ABS(BETAH(1,I)) <= EPSLON)THEN
            L = I+1
            GO TO 150
        ENDIF
    100 ENDDO

! IF WE GET TO THE NEXT STATEMENT, THEN THERE ARE NO ZERO OFF-DIAGONALS
! FOR THE SUBMATRIX WITH UPPER LEFT A(L,L) AND LOWER RIGHT A(N,N).
! WE CAN STILL GET EIGENVALUES IF THE MATRIX IS 2 BY 2 OR 1 BY 1.
! OTHERWISE, DO ANOTHER QR ITERATION PROVIDED ITMAX CYCLES HAVE
! NOT OCCURRED.

    IF(L == N .OR. L == N-1)THEN
        GO TO 150
    ELSE
        IF(ITER == ITMAX)THEN
            IERROR = 1
            GO TO 1000
        ELSE
            GO TO 200
        ENDIF
    ENDIF

! IF WE GET TO 150 THEN A(L,L-1) IS ZERO AND THE UNREDUCED SUBMATRIX
! HAS UPPER LEFT AT A(L,L) AND LOWER RIGHT AT A(N,N).  WE CAN
! EXTRACT ONE EIGENVALUE IF THIS MATRIX IS 1 BY 1 AND 2 EIGENVALUES
! IF IT IS 2 BY 2.

    150 IF(L == N)THEN
    
    ! IT'S A 1 BY 1 AND EVAL1(N) IS AN EIGENVALUE.  IF L=2 OR 1 WE ARE
    ! DONE.  OTHERWISE, UPDATE N, RESET L AND ITER, AND REPEAT THE
    ! SEARCH.
    
        IF(L <= 2)THEN
            GO TO 500
        ELSE
            N = L-1
            L = 1
            ITER = 0
            GO TO 80
        ENDIF
    ELSEIF(L == N-1)THEN
    
    ! THE UNREDUCED SUBMATRIX IS A 2 BY 2.  OVERWRITE EVAL1(N-1)
    ! AND EVAL1(N) WITH THE EIGENVALUES OF THE LOWER RIGHT 2 BY 2.
    
        BTERM = EVAL1(N-1) + EVAL1(N)
        ROOT1 = BTERM*0.5D0
        ROOT2 = ROOT1
        DISCR = BTERM**2 - 4.0D0*(EVAL1(N-1)*EVAL1(N)-BETAH(1,N-1)**2)
        IF(DISCR > 0.0D0)THEN
            D = DSQRT(DISCR)*0.5D0
            ROOT1 = ROOT1 - D
            ROOT2 = ROOT2 + D
        ENDIF
        EVAL1(N-1) = ROOT1
        EVAL1(N) = ROOT2
    
    ! SEE IF WE ARE DONE.  IF NOT, RESET N, L, AND ITER AND LOOK
    ! FOR NEXT UNREDUCED SUBMATRIX.
    
        IF(L <= 2)THEN
            GO TO 500
        ELSE
            N = L-1
            L = 1
            ITER = 0
            GO TO 80
        ENDIF
    ELSE
    
    ! AN EIGENVALUE WAS FOUND AND THE NEW UNREDUCED MATRIX LIMITS
    ! N AND L ARE SET.  DO A QR ITERATION ON NEW MATRIX.
    
        ITER = 0
        GO TO 200
    ENDIF

! QR ITERATION BEGINS HERE.

    200 ITER = ITER + 1

! USE EIGENVALUES OF THE LOWER RIGHT 2 BY 2 TO COMPUTE SHIFT.  SHIFT
! BY THE EIGENVALUE CLOSEST TO EVAL1(N).

    D = (EVAL1(N-1) - EVAL1(N))*0.5D0
    SIGND = 1.0D0
    IF(D < 0.0D0) SIGND = -1.0D0
    SHIFT = EVAL1(N) + D - SIGND*DSQRT(D*D + BETAH(1,N-1)**2)
    P = EVAL1(L) - SHIFT
    R = BETAH(1,L)
    T = EVAL1(L)
    W = BETAH(1,L)

! OVERWRITE A WITH Q'*A*Q.

    DO 250 K=L,N-1
        D = DSQRT(P*P + R*R)
        C = P/D
        S = R/D
        IF(K /= L) BETAH(1,K-1) = D
        CC = C*C
        SS = 1.0D0 - CC
        CS = C*S
        CSW = 2.0D0*CS*W
        AK1 = EVAL1(K+1)
        EVAL1(K) = CC*T + CSW + SS*AK1
        P = (CC - SS)*W + CS*(AK1 - T)
        T = SS*T - CSW + CC*AK1
        R = S*BETAH(1,K+1)
        W = C*BETAH(1,K+1)
    250 ENDDO
    BETAH(1,N-1) = P
    EVAL1(N) = T

! GO BACK AND SEE IF L AND N NEED TO BE UPDATED.

    GO TO 80

! SORT EIGENVALUES IN ASCENDING ALGEBRAIC ORDER.

    500 DO 600 I=2,NDIM
        JMAX = NDIM-I+1
        ISORT = 0
        DO 550 J=1,JMAX
            IF(EVAL1(J) > EVAL1(J+1))THEN
                ETEMP = EVAL1(J)
                EVAL1(J) = EVAL1(J+1)
                EVAL1(J+1) = ETEMP
                ISORT = 1
            ENDIF
        550 ENDDO
        IF(ISORT == 0) GO TO 1000
    600 ENDDO
    1000 RETURN
    end SUBROUTINE hessEIGVAL



    SUBROUTINE hessDEGEN(NDIM,EVAL1,TOLERA,ANORM,IDEGEN1)

! DETERMINES DEGENERACIES OF THE EIGENVALUES.

! INPUT:

! NDIM   = SIZE OF MATRIX BEING DIAGONALIZED.
! EVAL1   = SORTED EIGENVALUES (INCREASING VALUE).
! TOLERA = SAME TOLERANCE USED TO DETERMINE EIGENVALUES.
! ANORM  = ABSOLUTE COLUMN NORM OF TRIDIAGONAL MATRIX.


! RETURNED:

! IDEGEN1 = DEGENERACIES OF EIGENVALUES.

    use allmod
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
! DIMENSION EVAL1(*),IDEGEN1(*)
    DIMENSION EVAL1(natom*3),IDEGEN1(natom*3)

! DETERMINE DEGENERACIES OF EIGENVALUES.  ADJACENT EIGENVALUES
! WILL BE CONSIDERED TO BE DEGENERATE WHEN THEY DIFFER BY LESS
! THAN DTOLER.

    DTOLER = MAX(ANORM*DSQRT(TOLERA),1.0D-8)
    NSAME = 1
    DO 200 I=2,NDIM
        DIFF = ABS(EVAL1(I-1) - EVAL1(I))
        IF(DIFF <= DTOLER)THEN
        
        ! EIGENVALUES I-1 AND I ARE DEGENERATE.
        
            NSAME = NSAME + 1
            IF(I == NDIM)THEN
            
            ! WE'VE COME TO THE LAST REQUESTED EIGENVALUE, AND IT'S TIME
            ! TO ASSIGN DEGENERACIES FOR THE BLOCK ENDING WITH THE ITH
            ! EIGENVALUE.
            
                DO 100 J=I-NSAME+1,I
                    IDEGEN1(J) = NSAME
                100 ENDDO
            ENDIF
        
        ! GO TO THE NEXT EIGENVALUE (IF THERE ARE ANY LEFT) AND SEE IF
        ! IT'S DEGENERATE WITH THE NSAME EIGENVALUES WE'VE ALREADY
        ! FOUND.
        
            GO TO 200
        ELSE
        
        ! EITHER EIGENVALUE I-1 IS NONDEGENERATE OR IT'S THE LAST
        ! EIGENVALUE IN A DEGENERATE BLOCK.  CORRESPONDINGLY, ASSIGN THE
        ! PROPER DEGENERACY TO I-1 OR TO EACH EIGENVALUE IN THE BLOCK.
        
            DO 150 J=I-NSAME,I-1
                IDEGEN1(J) = NSAME
            150 ENDDO
            NSAME = 1
        
        ! IF I=NDIM THEN IT MUST BE THE CASE THAT THIS LAST EIGENVALUE
        ! IS NONDEGENERATE.
        
            IF(I == NDIM) IDEGEN1(I) = 1
        ENDIF
    200 ENDDO
    RETURN
    end SUBROUTINE hessDEGEN



    SUBROUTINE hessEIGVEC(NDIM,NEVEC1,A,AWORK,TOLERA,ANORM,EVAL1,IDEGEN1, &
    EVEC1)

! INVERSE ITERATION ROUTINE FOR EIGENVECTOR DETERMINATION.
! CALCULATES THE EIGENVECTORS OF AN NDIM BY NDIM SYMMETRIC,
! TRIDIAGONAL MATRIX A.

! INPUT:

! NDIM   = SIZE OF MATRIX A.
! NEVEC1  = NUMBER OF EIGENVECTORS REQUIRED.
! A      = NDIM BY NDIM TRIDIAGONAL MATRIX.
! TOLERA = SAME TOLERANCE USED TO DETERMINE EIGENVALUES.
! ANORM  = ABSOLUTE COLUMN NORM OF TRIDIAGONAL MATRIX A.
! EVAL1   = SORTED EIGENVALUES OF A.
! IDEGEN1 = DEGENERACIES OF EIGENVALUES.


! RETURNED:

! EVEC1   = EIGENVECTORS OF TRIDIAGONAL MATRIX (IN COLUMNS).

! PROGRAMMED BY S. L. DIXON, OCT., 1991.


    use allmod
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
! DIMENSION A(NDIM,*),AWORK(3,*),EVAL1(*),IDEGEN1(*),EVEC1(NDIM,*)
    DIMENSION A(natom*3,natom*3),AWORK(3,natom*3),EVAL1(natom*3), &
    IDEGEN1(natom*3),EVEC1(natom*3,natom*3)
    LOGICAL :: ORTH
    IRAND = 13876532

! COMPUTE THRESHOLD EPSLON WHICH WILL BE USED IF THE INVERSE ITERATION
! MATRIX IS SINGULAR.

    EPSLON = ANORM*TOLERA

! WHEN DEGENERACIES OCCUR, THERE ARE RARE INSTANCES WHEN THE
! DEGENERATE BLOCK OF EIGENVECTORS ARE NOT LINEARLY INDEPENDENT.
! IN THESE CASES, AN ADDITIONAL PASS THROUGH THE INVERSE ITERATION
! (WITH A NEW SET OF RANDOM NUMBERS) IS CARRIED OUT, I.E., CONTROL
! PASSES TO STATEMENT 40.  IVECT WILL KEEP TRACK OF THE CURRENT
! STARTING EIGENVECTOR WHEN ADDITIONAL PASSES ARE NECESSARY.

    IVECT = 1
    NFAIL = 0
    NPRTRB = 0

! DO ONE ITERATION FOR EACH EIGENVECTOR.

    40 ORTH = .TRUE.
    NDEGEN = 0
    MSTART = IVECT
    DO 380 M=MSTART,NEVEC1
        IF(IDEGEN1(M) > 1) NDEGEN = NDEGEN + 1
        Z = EVAL1(M)
    
    ! IF THE INVERSE ITERATION HAS FAILED TWICE DUE TO NON-ORTHOGONALITY
    ! OF DEGENERATE EIGENVECTORS, PERTURB THE EIGENVALUE BY A SMALL
    ! AMOUNT.
    
        IF(NFAIL >= 2)THEN
            NPRTRB = NPRTRB + 1
            Z = Z + 0.001D0*DBLE(NPRTRB)*TOLERA
        ENDIF
    
    ! STORE THE TRIDIAGONAL ENTRIES OF THE INVERSE ITERATION MATRIX IN
    ! THE 3 BY NDIM WORKSPACE AWORK.
    
        AWORK(1,1) = 0.0D0
        AWORK(2,1) = A(1,1) - Z
        AWORK(3,1) = A(2,1)
        IF(NDIM > 2)THEN
            DO 80 I=2,NDIM-1
                AWORK(1,I) = A(I,I-1)
                AWORK(2,I) = A(I,I) - Z
                AWORK(3,I) = A(I+1,I)
            80 ENDDO
        ENDIF
        AWORK(1,NDIM) = A(NDIM,NDIM-1)
        AWORK(2,NDIM) = A(NDIM,NDIM) - Z
        AWORK(3,NDIM) = 0.0D0
    
    ! ASSIGN INVERSE ITERATION VECTOR FROM RANDOM NUMBERS.
    
        DO 120 I=1,NDIM
            CALL RANDOM(IRAND,RNDOM)
            RNDOM = 2.0D0*(RNDOM - 0.5D0)
            EVEC1(I,M) = RNDOM*TOLERA
        120 ENDDO
    
    ! CARRY OUT FORWARD GAUSSIAN ELIMINATION WITH ROW PIVOTING
    ! ON THE INVERSE ITERATION MATRIX.
    
        DO 160 K=1,NDIM-1
            ADIAG = ABS(AWORK(2,K))
            ASUB = ABS(AWORK(1,K+1))
            IF(ADIAG >= ASUB)THEN
            
            ! USE PIVOTAL ELEMENT FROM ROW K.
            
                IF(AWORK(2,K) == 0.0D0) AWORK(2,K) = EPSLON
                T = AWORK(1,K+1)/AWORK(2,K)
                AWORK(1,K+1) = 0.0D0
                AWORK(2,K+1) = AWORK(2,K+1) - T*AWORK(3,K)
            
            ! LEFT-JUSTIFY EQUATION K SO THAT DIAGONAL ENTRY IS STORED
            ! IN AWORK(1,K).
            
                AWORK(1,K) = AWORK(2,K)
                AWORK(2,K) = AWORK(3,K)
                AWORK(3,K) = 0.0D0
            
            ! OPERATE ON VECTOR AS WELL.
            
                EVEC1(K+1,M) = EVEC1(K+1,M) - T*EVEC1(K,M)
            ELSE
            
            ! USE PIVOTAL ELEMENT FROM ROW K+1 AND SWAP ROWS K AND K+1.
            
                IF(AWORK(1,K+1) == 0.0D0) AWORK(1,K+1) = EPSLON
                T = AWORK(2,K)/AWORK(1,K+1)
                ATEMP = AWORK(3,K) - T*AWORK(2,K+1)
                AWORK(1,K) = AWORK(1,K+1)
                AWORK(2,K) = AWORK(2,K+1)
                AWORK(3,K) = AWORK(3,K+1)
                AWORK(1,K+1) = 0.0D0
                AWORK(2,K+1) = ATEMP
                AWORK(3,K+1) = -T*AWORK(3,K+1)
            
            ! OPERATE ON VECTOR AND SWAP ENTRIES.
            
                ETEMP = EVEC1(K+1,M)
                EVEC1(K+1,M) = EVEC1(K,M) - ETEMP*T
                EVEC1(K,M) = ETEMP
            ENDIF
        160 ENDDO
    
    ! FORWARD ELIMINATION COMPLETE.  BACK SUBSTITUTE TO GET SOLUTION.
    ! OVERWRITE COLUMN M OF EVEC1 WITH SOLUTION.
    
        IF(AWORK(2,NDIM) == 0.0D0) AWORK(2,NDIM) = EPSLON
        EVEC1(NDIM,M) = EVEC1(NDIM,M)/AWORK(2,NDIM)
        ETEMP = EVEC1(NDIM-1,M) - AWORK(2,NDIM-1)*EVEC1(NDIM,M)
        EVEC1(NDIM-1,M) = ETEMP/AWORK(1,NDIM-1)
        ENORM = EVEC1(NDIM,M)**2 + EVEC1(NDIM-1,M)**2
        IF(NDIM > 2)THEN
        
        ! CAUTION: PROBLEM LOOP FOR SOME IBM RS/6000 COMPILERS.  VALUE
        ! OF K CAN GET LOST WHEN OPTIMIZE FLAG IS USED.
        
            DO 200 L=1,NDIM-2
                K = NDIM-L-1
                ETEMP = EVEC1(K,M) - AWORK(2,K)*EVEC1(K+1,M) &
                - AWORK(3,K)*EVEC1(K+2,M)
                EVEC1(K,M) = ETEMP/AWORK(1,K)
                ENORM = ENORM + EVEC1(K,M)**2
            200 ENDDO
        ENDIF
        EINV = 1.0D0/DSQRT(ENORM)
    
    ! NORMALIZE EIGENVECTOR.
    
        DO 240 I=1,NDIM
            EVEC1(I,M) = EVEC1(I,M)*EINV
        240 ENDDO
    
    ! IF WE HAVE COME TO THE END OF A DEGENERATE BLOCK OF EIGENVECTORS,
    ! ORTHOGONALIZE THE BLOCK.
    
        IF(NDEGEN > 1)THEN
            IF(NDEGEN == IDEGEN1(M) .OR. M == NEVEC1)THEN
                JSTART = M-NDEGEN+1
                CALL hessORTHOG(NDIM,NDEGEN,JSTART,EVEC1,ORTH)
                IF(ORTH)THEN
                    NFAIL = 0
                    NPRTRB = 0
                
                ! THE DEGENERATE VECTORS WERE LINEARLY INDEPENDENT AND WERE
                ! SUCCESSFULLY ORTHOGONALIZED.
                
                    IVECT = IVECT + NDEGEN
                    NDEGEN = 0
                ELSE
                
                ! THE BLOCK IS APPARENTLY NOT LINEARLY INDEPENDENT.  GO BACK
                ! AND REPEAT THE INVERSE ITERATION FOR THESE VECTORS.  AFTER
                ! AN INDEPENDENT SET HAS BEEN FOUND, ANY ADDITIONAL EIGENVECTORS
                ! WILL BE DETERMINED.
                
                    NFAIL = NFAIL + 1
                    GO TO 40
                ENDIF
            ENDIF
        ENDIF
    
    ! THE CURRENT EIGENVECTOR SHOULD BE OKAY IF IT IS NONDEGENERATE.
    
        IF(IDEGEN1(M) == 1) IVECT = IVECT + 1
    380 ENDDO
    RETURN
    end SUBROUTINE hessEIGVEC



    SUBROUTINE hessORTHOG(NDIM,NVECT,JSTART,VECT,ORTH)

! CONSTRUCTS A SET OF ORTHONORMAL VECTORS FROM THE NVECT LINEARLY
! INDEPENDENT, NORMALIZED VECTORS IN THE ARRAY VECT.  THE VECTORS
! SHOULD BE STORED COLUMNWISE, STARTING IN COLUMN JSTART.  VECT IS
! OVERWRITTEN WITH THE ORTHONORMAL SET.  ALL VECTORS ARE NDIM BY 1.
! ORTH IS RETURNED WITH A VALUE OF .TRUE. IF THE SET WAS LINEARLY
! INDEPENDENT AND .FALSE. OTHERWISE.

! PROGRAMMED BY S. L. DIXON.


    use allmod
    IMPLICIT DOUBLE PRECISION (A-H,O-Z)
! DIMENSION VECT(NDIM,*)
    DIMENSION VECT(natom*3,natom*3)
    LOGICAL :: ORTH

    ORTH = .TRUE.
    ORTEST = 1.0D-8

! BEGIN ORTHOGONALIZATION.

    JSTOP = JSTART + NVECT - 1
    DO 120 J=JSTART,JSTOP
        IF(J > JSTART)THEN
        
        ! SUBTRACT OFF COMPONENTS OF PREVIOUSLY DETERMINED ORTHOGONAL
        ! VECTORS FROM THE VECTOR IN COLUMN J.
        
            DO 60 JPREV=JSTART,J-1
                DOT = 0.0D0
                DO 20 I=1,NDIM
                    DOT = DOT + VECT(I,JPREV)*VECT(I,J)
                20 ENDDO
                DO 40 I=1,NDIM
                    VECT(I,J) = VECT(I,J) - DOT*VECT(I,JPREV)
                40 ENDDO
            60 ENDDO
        ENDIF
    
    ! NORMALIZE COLUMN J.
    
        VJNORM = 0.0D0
        DO 80 I=1,NDIM
            VJNORM = VJNORM + VECT(I,J)**2
        80 ENDDO
        VJNORM = DSQRT(VJNORM)
    
    ! IF THE NORM OF THIS VECTOR IS TOO SMALL THEN THE VECTORS ARE
    ! NOT LINEARLY INDEPENDENT.
    
        IF(VJNORM < ORTEST)THEN
            ORTH = .FALSE.
            GO TO 1000
        ENDIF
        DO 100 I=1,NDIM
            VECT(I,J) = VECT(I,J)/VJNORM
        100 ENDDO
    120 ENDDO
    1000 RETURN
    end SUBROUTINE hessORTHOG


subroutine PriHessian(io,n,mat,fm) ! format: f(x.y) x>7 sugg 12.5,12.7,14.9
  implicit none
  integer j,jj,n,io,n5,nf,x,y,ini,ifi,k,i,iatom,imom
  real*8 mat(n,n)
  character fm*(*),ch,fm2*10
  character*40 fmt1,fmt2,fmt3,fmt4
  character(len=1) cartsym(3)
  character(len=5) name(n)
  cartsym(1) = 'X'
  cartsym(2) = 'Y'
  cartsym(3) = 'Z'

  iatom=1
  imom=0
  do i=1,n
    imom=imom+1
    if (imom.eq.4) then
        imom=1
        iatom=iatom+1
    endif
    write(name(i),'(i4,a1)') iatom,cartsym(imom)
  enddo
    
  
  n5=n/5
  nf=mod(n,5)
  fm2=fm
  ch=fm2(1:1)
  k=index(fm2,'.')
  read(fm2(2:k-1),*) x
  read(fm2(k+1:10),*) y

  write(fmt1,101) ch,x,y
  write(fmt2,102) nf,ch,x,y
101 format('(a5,5',a1,i2,'.',i2,')')
102 format('(a5,',i2,a1,i2,'.',i2,')')
  write(fmt3,103) x-7
  write(fmt4,104) nf
103 format('(3x,5(',i2,'x,i7))')
104 format('(1x,',i2,'(7x,a5))')

  do jj=1,n5
     ini=1+(jj-1)*5
     write(io,'(8x,5(a5,7x))') (name(j),j=ini,jj*5)
     do k=1+(jj-1)*5,n
        ifi=min(jj*5,k)
        write(io,fmt1) name(k),(mat(k,j),j=ini,ifi)
     enddo
     !         if (jj.ne.n5.or.nf.ne.0) write(io,*)
  enddo

  if (nf.ne.0) then
     ini=n-nf+1
     write(io,fmt4)(name(j),j=ini,n)
     do k=ini,n
        write(io,fmt2) name(k),(mat(k,j),j=ini,k)
     enddo
  endif
  call flush(io)

end subroutine PriHessian
