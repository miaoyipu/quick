
! Ed Brothers. November 14, 2002.
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

subroutine formCPHFA(Ibas,Jbas,IIbas,JJbas)
  use allmod
  implicit double precision(a-h,o-z)

  dimension isame(4,8)
  logical :: same

  ! The purpose of the subroutine is to calculate an AO repulsion
  ! integral, determine it's contribution to an MO repulsion integral,
  ! and put it in the correct location in the CPHF A matrix.

  ! First, calculate the integral.

  iA = ncenter(Ibas)
  iB = ncenter(Jbas)
  iC = ncenter(IIbas)
  iD = ncenter(JJbas)
  AOint=0.d0
  DO Icon=1,ncontract(Ibas)
     DO Jcon=1,ncontract(Jbas)
        DO IIcon=1,ncontract(IIbas)
           DO JJcon=1,ncontract(JJbas)
              AOint = AOint + &
                   dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                   *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                   *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                   aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                   itype(1,Ibas), itype(2,Ibas), itype(3,Ibas), &
                   itype(1,Jbas), itype(2,Jbas), itype(3,Jbas), &
                   itype(1,IIbas),itype(2,IIbas),itype(3,IIbas), &
                   itype(1,JJbas),itype(2,JJbas),itype(3,JJbas), &
                   xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                   xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                   xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                   xyz(1,iD), xyz(2,iD),xyz(3,iD))

           ENDDO
        ENDDO
     ENDDO
  ENDDO
  
  ! Now we need to find how many times the AO integral appears by examining
  ! it's symmetry.  For an integral (ij|kl):
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

  ! Now we can start filling up the CPHFA array.


  DO Iunique=1,8
     IF (isame(1,Iunique) /= 0) THEN

        ! Set up some dummy variables.

        Ival = isame(1,Iunique)
        Jval = isame(2,Iunique)
        IIval = isame(3,Iunique)
        JJval = isame(4,Iunique)

        ! Loop over alpha pairs.

        DO iAvirt = lastAocc+1,nbasis
           DO iAocc = 1,lastAocc

              ! iAvirt and iAocc form an ai pair.  Find it's location.

              iaCPHFA = (iAvirt-lastAocc-1)*lastAocc + iAocc

              ! Loop over alpha again.  This means we are filling the alpha-alpha
              ! portion of the matrix.

              DO iAvirt2 = lastAocc+1,nbasis
                 DO iAocc2 = 1,lastAocc

                    ! iAvirt2 and iAocc2 form an bj pair.  Find it's location.

                    jbCPHFA = (iAvirt2-lastAocc-1)*lastAocc + iAocc2

                    ! CPHFA(ai,jb) = 2(ai|jb) - (aj|bi) - (ab|ji)
                    ! Since all these are alpha, no elements drop out.
                    ! Calculate the value of the AO repulsions contribution to the MO repulsion.

                    CPHFA(iaCPHFA,jbCPHFA) = CPHFA(iaCPHFA,jbCPHFA) + &
                         2.d0*CO(Ival,IAvirt)*CO(Jval,IAocc)* &
                         CO(IIval,IAvirt2)*CO(JJval,iAocc2)*AOint &
                         -CO(Ival,IAvirt)*CO(Jval,iAocc2)* &
                         CO(IIval,IAvirt2)*CO(JJval,iAocc)*AOint &
                         -CO(Ival,IAvirt)*CO(Jval,iAvirt2)* &
                         CO(IIval,IAocc2)*CO(JJval,iAocc)*AOint

                 ENDDO
              ENDDO

              ! Now loop over beta for the bj pairs.

              DO iBvirt2 = lastBocc+1,nbasis
                 DO iBocc2 = 1,lastBocc

                    ! iBvirt2 and iBocc2 form an bj pair.  Find it's location.

                    jbCPHFA = (iBvirt2-lastBocc-1)*lastBocc +iBocc2 +iBetastart

                    ! CPHFA(ai,jb) = 2(ai|jb) - (aj|bi) - (ab|ji)
                    ! j and b are beta, thus it becomes:
                    ! CPHFA(ai,jb) = 2(ai|jb)
                    ! Calculate the value of the AO repulsions contribution to the MO repulsion.

                    CPHFA(iaCPHFA,jbCPHFA) = CPHFA(iaCPHFA,jbCPHFA) + &
                         2.d0*CO(Ival,IAvirt)*CO(Jval,IAocc)* &
                         COB(IIval,IBvirt2)*COB(JJval,iBocc2)*AOint

                 ENDDO
              ENDDO
           ENDDO
        ENDDO

        ! Loop over beta pairs.

        DO iBvirt = lastBocc+1,nbasis
           DO iBocc = 1,lastBocc

              ! iBvirt and iBocc form an ai pair.  Find it's location.

              iaCPHFA = (iBvirt-lastBocc-1)*lastBocc +iBocc +iBetastart

              ! Loop over beta again.  This means we are filling the beta-beta
              ! portion of the matrix.

              DO iBvirt2 = lastBocc+1,nbasis
                 DO iBocc2 = 1,lastBocc

                    ! iBvirt2 and iBocc2 form an bj pair.  Find it's location.

                    jbCPHFA = (iBvirt2-lastBocc-1)*lastBocc +iBocc2 +iBetastart

                    ! CPHFA(ai,jb) = 2(ai|jb) - (aj|bi) - (ab|ji)
                    ! Since all these are beta, no elements drop out.
                    ! Calculate the value of the AO repulsions contribution to the MO repulsion.

                    CPHFA(iaCPHFA,jbCPHFA) = CPHFA(iaCPHFA,jbCPHFA) + &
                         2.d0*COB(Ival,IBvirt)*COB(Jval,IBocc)* &
                         COB(IIval,IBvirt2)*COB(JJval,iBocc2)*AOint &
                         -COB(Ival,IBvirt)*COB(Jval,iBocc2)* &
                         COB(IIval,IBvirt2)*COB(JJval,iBocc)*AOint &
                         -COB(Ival,IBvirt)*COB(Jval,iBvirt2)* &
                         COB(IIval,IBocc2)*COB(JJval,iBocc)*AOint

                 ENDDO
              ENDDO

              ! Now loop over alpha for the bj pairs.

              DO iAvirt2 = lastAocc+1,nbasis
                 DO iAocc2 = 1,lastAocc

                    ! iAvirt2 and iAocc2 form an bj pair.  Find it's location.

                    jbCPHFA = (iAvirt2-lastAocc-1)*lastAocc + iAocc2

                    ! CPHFA(ai,jb) = 2(ai|jb) - (aj|bi) - (ab|ji)
                    ! j and b are alpha, thus it becomes:
                    ! CPHFA(ai,jb) = 2(ai|jb)
                    ! Calculate the value of the AO repulsions contribution to the MO repulsion.

                    CPHFA(iaCPHFA,jbCPHFA) = CPHFA(iaCPHFA,jbCPHFA) + &
                         2.d0*COB(Ival,IBvirt)*COB(Jval,IBocc)* &
                         CO(IIval,IAvirt2)*CO(JJval,iAocc2)*AOint

                 ENDDO
              ENDDO
           ENDDO
        ENDDO

     ENDIF
  ENDDO



  return
end subroutine formCPHFA



! Ed Brothers. November 18, 2002.
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

subroutine formCPHFB
  use allmod
  implicit double precision(a-h,o-z)

  ! First some quick setup.

  IF (quick_method%unrst) THEN
     lastAocc = nelec
     lastBocc = nelecb
  ELSE
     lastAocc = nelec/2
     lastBocc = lastAocc
  ENDIF
  iBetastart = lastAocc*(nbasis-lastAocc)

  ! The purpose of this subroutine is to form B0 for the CPHF. There is
  ! one B0 for each nuclear perturbation, and this subroutine forms all
  ! of them.  The elements of B0(ai) = Qai(1)/(E(i)-E(a))

  ! Now some math for the B0 where the prreturbation is the y component
  ! of atom1.  (Y1)  This is similar for any perturbation.

  ! Qai(1) = Hai(1) - Sai(1) E(i) - (Sum over k and l) Skl(1) [(ai|lk)-(ak|li)]
  ! + (Sum over mu,nu,lambda,sigma) C(mu,a) C(nu,i) DENSE(lamba,sigma)
  ! * d/dY1 [(mu nu|lamba sigma) - (mu sigma| lamba nu)]

  ! Now Hai(1) = (Sum over mu,nu) C(mu,a) C(nu,i) Hmunu(1)
  ! = (Sum over mu,nu) C(mu,a) C(nu,i) d/dY1 H(mu,nu)

  ! And Sai(1) = (Sum over mu,nu) C(mu,a) C(nu,i) Smunu(1)
  ! = (Sum over mu,nu) C(mu,a) C(nu,i) d/dY1 S(mu,nu)

  ! We are going to calculate the first two terms: Hai(1) - Sai(1) E(i)

  DO Ibas=1,nbasis
     ISTART = (ncenter(Ibas)-1) *3
     DO Jbas=ilast(ncenter(IBAS))+1,nbasis
        JSTART = (ncenter(Jbas)-1) *3

        ! We have selected our two basis functions, now loop over angular momentum.

        DO Imomentum=1,3
           dSI = 0.d0
           dSJ =0.d0
           dKEI = 0.d0
           dKEJ = 0.d0

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
                    dSI = dSI - dble(itype(Imomentum,Ibas)+1)* &
                         dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *overlap(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                         itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                         itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                         xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                         xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                         xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
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
                    dSJ = dSJ - dble(itype(Imomentum,Jbas)+1)* &
                         dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                         *overlap(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                         itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                         itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                         xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                         xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                         xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))
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

           ! At this point we have the derivatives.  Now we need to put them in
           ! the CPHFB array.
           ! ALPHA first.

           DO iAvirt = lastAocc+1,nbasis
              DO iAocc = 1,lastAocc
                 iaCPHF = (iAvirt-lastAocc-1)*lastAocc + iAocc
                 CPHFB(iaCPHF,ISTART+Imomentum) = &
                      CPHFB(iaCPHF,ISTART+Imomentum) &
                      +dKEI*(CO(Ibas,iAvirt)*CO(Jbas,iAocc) &
                      +CO(Ibas,iAocc)*CO(Jbas,iAvirt)) &
                      -dSI* (CO(Ibas,iAvirt)*CO(Jbas,iAocc) &
                      +CO(Ibas,iAocc)*CO(Jbas,iAvirt))*E(iAocc)
                 CPHFB(iaCPHF,JSTART+Imomentum) = &
                      CPHFB(iaCPHF,JSTART+Imomentum) &
                      +dKEJ*(CO(Jbas,iAvirt)*CO(Ibas,iAocc) &
                      +CO(Jbas,iAocc)*CO(Ibas,iAvirt)) &
                      -dSJ* (CO(Jbas,iAvirt)*CO(Ibas,iAocc) &
                      +CO(Jbas,iAocc)*CO(Ibas,iAvirt))*E(iAocc)
              ENDDO
           ENDDO
           ! BETA
           DO iBvirt = lastBocc+1,nbasis
              DO iBocc = 1,lastBocc
                 iaCPHF = (iBvirt-lastBocc-1)*lastBocc + iBocc +iBetaStart
                 CPHFB(iaCPHF,ISTART+Imomentum) = &
                      CPHFB(iaCPHF,ISTART+Imomentum) &
                      +dKEI*(COB(Ibas,iBvirt)*COB(Jbas,iBocc) &
                      +COB(Ibas,iBocc)*COB(Jbas,iBvirt)) &
                      -dSI* (COB(Ibas,iBvirt)*COB(Jbas,iBocc) &
                      +COB(Ibas,iBocc)*COB(Jbas,iBvirt))*EB(iBocc)
                 CPHFB(iaCPHF,JSTART+Imomentum) = &
                      CPHFB(iaCPHF,JSTART+Imomentum) &
                      +dKEJ*(COB(Jbas,iBvirt)*COB(Ibas,iBocc) &
                      +COB(Jbas,iBocc)*COB(Ibas,iBvirt)) &
                      -dSJ* (COB(Jbas,iBvirt)*COB(Ibas,iBocc) &
                      +COB(Jbas,iBocc)*COB(Ibas,iBvirt))*EB(iBocc)
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO


  ! We still have to calculate the three center term that arises in the
  ! core Hamiltonian.

  DO Ibas=1,nbasis
     iA=ncenter(Ibas)
     ISTART = (iA-1)*3

     DO Jbas=Ibas,nbasis
        iB = ncenter(Jbas)
        JSTART = (iB-1)*3

        DO iC = 1,natom
           iCSTART = (iC-1)*3

           ! As before, if all terms are on the same atom, they move with the
           ! atom and the derivative is zero.

           IF (iA == iC .AND. iB == iC) THEN
              continue
           ELSE


              ! If Ibas=Jbas, the term only shows up once in the energy, otherwise
              ! it shows up twice. This is not an issue with the 2-center terms above.


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

              ! Now add these into the CPHFB array.

              DO iAvirt = lastAocc+1,nbasis
                 DO iAocc = 1,lastAocc
                    iaCPHF = (iAvirt-lastAocc-1)*lastAocc + iAocc
                    CPHFB(iaCPHF,ISTART+1) = CPHFB(iaCPHF,ISTART+1) &
                         +dNAIX*(CO(Ibas,iAvirt)*CO(Jbas,iAocc))
                    CPHFB(iaCPHF,ISTART+2) = CPHFB(iaCPHF,ISTART+2) &
                         +dNAIY*(CO(Ibas,iAvirt)*CO(Jbas,iAocc))
                    CPHFB(iaCPHF,ISTART+3) = CPHFB(iaCPHF,ISTART+3) &
                         +dNAIZ*(CO(Ibas,iAvirt)*CO(Jbas,iAocc))
                    CPHFB(iaCPHF,JSTART+1) = CPHFB(iaCPHF,JSTART+1) &
                         +dNAJX*(CO(Ibas,iAvirt)*CO(Jbas,iAocc))
                    CPHFB(iaCPHF,JSTART+2) = CPHFB(iaCPHF,JSTART+2) &
                         +dNAJY*(CO(Ibas,iAvirt)*CO(Jbas,iAocc))
                    CPHFB(iaCPHF,JSTART+3) = CPHFB(iaCPHF,JSTART+3) &
                         +dNAJZ*(CO(Ibas,iAvirt)*CO(Jbas,iAocc))
                    CPHFB(iaCPHF,ICSTART+1) = CPHFB(iaCPHF,ICSTART+1) &
                         +dNACX*(CO(Ibas,iAvirt)*CO(Jbas,iAocc))
                    CPHFB(iaCPHF,ICSTART+2) = CPHFB(iaCPHF,ICSTART+2) &
                         +dNACY*(CO(Ibas,iAvirt)*CO(Jbas,iAocc))
                    CPHFB(iaCPHF,ICSTART+3) = CPHFB(iaCPHF,ICSTART+3) &
                         +dNACZ*(CO(Ibas,iAvirt)*CO(Jbas,iAocc))
                    IF (Ibas /= Jbas) THEN
                       CPHFB(iaCPHF,ISTART+1) = CPHFB(iaCPHF,ISTART+1) &
                            +dNAIX*(CO(Jbas,iAvirt)*CO(Ibas,iAocc))
                       CPHFB(iaCPHF,ISTART+2) = CPHFB(iaCPHF,ISTART+2) &
                            +dNAIY*(CO(Jbas,iAvirt)*CO(Ibas,iAocc))
                       CPHFB(iaCPHF,ISTART+3) = CPHFB(iaCPHF,ISTART+3) &
                            +dNAIZ*(CO(Jbas,iAvirt)*CO(Ibas,iAocc))
                       CPHFB(iaCPHF,JSTART+1) = CPHFB(iaCPHF,JSTART+1) &
                            +dNAJX*(CO(Jbas,iAvirt)*CO(Ibas,iAocc))
                       CPHFB(iaCPHF,JSTART+2) = CPHFB(iaCPHF,JSTART+2) &
                            +dNAJY*(CO(Jbas,iAvirt)*CO(Ibas,iAocc))
                       CPHFB(iaCPHF,JSTART+3) = CPHFB(iaCPHF,JSTART+3) &
                            +dNAJZ*(CO(Jbas,iAvirt)*CO(Ibas,iAocc))
                       CPHFB(iaCPHF,ICSTART+1) = CPHFB(iaCPHF,ICSTART+1) &
                            +dNACX*(CO(Jbas,iAvirt)*CO(Ibas,iAocc))
                       CPHFB(iaCPHF,ICSTART+2) = CPHFB(iaCPHF,ICSTART+2) &
                            +dNACY*(CO(Jbas,iAvirt)*CO(Ibas,iAocc))
                       CPHFB(iaCPHF,ICSTART+3) = CPHFB(iaCPHF,ICSTART+3) &
                            +dNACZ*(CO(Jbas,iAvirt)*CO(Ibas,iAocc))
                    ENDIF
                 ENDDO
              ENDDO
              ! BETA
              DO iBvirt = lastBocc+1,nbasis
                 DO iBocc = 1,lastBocc
                    iaCPHF = (iBvirt-lastBocc-1)*lastBocc + iBocc + IbetaSTART
                    CPHFB(iaCPHF,ISTART+1) = CPHFB(iaCPHF,ISTART+1) &
                         +dNAIX*(COB(Ibas,iBvirt)*COB(Jbas,iBocc))
                    CPHFB(iaCPHF,ISTART+2) = CPHFB(iaCPHF,ISTART+2) &
                         +dNAIY*(COB(Ibas,iBvirt)*COB(Jbas,iBocc))
                    CPHFB(iaCPHF,ISTART+3) = CPHFB(iaCPHF,ISTART+3) &
                         +dNAIZ*(COB(Ibas,iBvirt)*COB(Jbas,iBocc))
                    CPHFB(iaCPHF,JSTART+1) = CPHFB(iaCPHF,JSTART+1) &
                         +dNAJX*(COB(Ibas,iBvirt)*COB(Jbas,iBocc))
                    CPHFB(iaCPHF,JSTART+2) = CPHFB(iaCPHF,JSTART+2) &
                         +dNAJY*(COB(Ibas,iBvirt)*COB(Jbas,iBocc))
                    CPHFB(iaCPHF,JSTART+3) = CPHFB(iaCPHF,JSTART+3) &
                         +dNAJZ*(COB(Ibas,iBvirt)*COB(Jbas,iBocc))
                    CPHFB(iaCPHF,ICSTART+1) = CPHFB(iaCPHF,ICSTART+1) &
                         +dNACX*(COB(Ibas,iBvirt)*COB(Jbas,iBocc))
                    CPHFB(iaCPHF,ICSTART+2) = CPHFB(iaCPHF,ICSTART+2) &
                         +dNACY*(COB(Ibas,iBvirt)*COB(Jbas,iBocc))
                    CPHFB(iaCPHF,ICSTART+3) = CPHFB(iaCPHF,ICSTART+3) &
                         +dNACZ*(COB(Ibas,iBvirt)*COB(Jbas,iBocc))
                    IF (Ibas /= Jbas) THEN
                       CPHFB(iaCPHF,ISTART+1) = CPHFB(iaCPHF,ISTART+1) &
                            +dNAIX*(COB(Jbas,iBvirt)*COB(Ibas,iBocc))
                       CPHFB(iaCPHF,ISTART+2) = CPHFB(iaCPHF,ISTART+2) &
                            +dNAIY*(COB(Jbas,iBvirt)*COB(Ibas,iBocc))
                       CPHFB(iaCPHF,ISTART+3) = CPHFB(iaCPHF,ISTART+3) &
                            +dNAIZ*(COB(Jbas,iBvirt)*COB(Ibas,iBocc))
                       CPHFB(iaCPHF,JSTART+1) = CPHFB(iaCPHF,JSTART+1) &
                            +dNAJX*(COB(Jbas,iBvirt)*COB(Ibas,iBocc))
                       CPHFB(iaCPHF,JSTART+2) = CPHFB(iaCPHF,JSTART+2) &
                            +dNAJY*(COB(Jbas,iBvirt)*COB(Ibas,iBocc))
                       CPHFB(iaCPHF,JSTART+3) = CPHFB(iaCPHF,JSTART+3) &
                            +dNAJZ*(COB(Jbas,iBvirt)*COB(Ibas,iBocc))
                       CPHFB(iaCPHF,ICSTART+1) = CPHFB(iaCPHF,ICSTART+1) &
                            +dNACX*(COB(Jbas,iBvirt)*COB(Ibas,iBocc))
                       CPHFB(iaCPHF,ICSTART+2) = CPHFB(iaCPHF,ICSTART+2) &
                            +dNACY*(COB(Jbas,iBvirt)*COB(Ibas,iBocc))
                       CPHFB(iaCPHF,ICSTART+3) = CPHFB(iaCPHF,ICSTART+3) &
                            +dNACZ*(COB(Jbas,iBvirt)*COB(Ibas,iBocc))
                    ENDIF
                 ENDDO
              ENDDO

           ENDIF
        ENDDO
     ENDDO
  ENDDO


  ! We've now done all the two and three center terms.  Now we need to
  ! add the 4 center terms into the B0 array.  These are:
  ! - (Sum over k and l) Skl(1) [(ai|lk)-(ak|li)]
  ! + (Sum over mu,nu,lambda,sigma) C(mu,a) C(nu,i) DENSE(lamba,sigma)
  ! * d/dY1 [(mu nu|lamba sigma) - (mu sigma| lamba nu)]
  ! We'll do this in a subprogram which is called for each unique integral.


  DO I=1,nbasis
     call CPHFB4cnt(I,I,I,I)
     DO J=I+1,nbasis
        call CPHFB4cnt(I,I,J,J)
        call CPHFB4cnt(I,J,J,J)
        call CPHFB4cnt(I,I,I,J)
        call CPHFB4cnt(I,J,I,J)
        DO K=J+1,nbasis
           call CPHFB4cnt(I,J,I,K)
           call CPHFB4cnt(I,J,K,K)
           call CPHFB4cnt(I,K,J,J)
           call CPHFB4cnt(I,I,J,K)
        ENDDO
        DO K=I+1,nbasis-1
           DO L=K+1,nbasis
              call CPHFB4cnt(I,J,K,L)
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  ! Now we need to go through and divide by (Ei-Ea).


  DO iAvirt = lastAocc+1,nbasis
     DO iAocc = 1,lastAocc
        iaCPHF = (iAvirt-lastAocc-1)*lastAocc + iAocc
        denom = E(iAocc)-E(iAvirt)
        DO IDX = 1,natom*3
           CPHFB(iaCPHF,IDX) = CPHFB(iaCPHF,idX)/denom
        ENDDO
     ENDDO
  ENDDO

  DO iBvirt = lastBocc+1,nbasis
     DO iBocc = 1,lastBocc
        iaCPHF = (iBvirt-lastBocc-1)*lastBocc + iBocc+iBetastart
        denom = EB(iBocc)-EB(iBvirt)
        DO IDX = 1,natom*3
           CPHFB(iaCPHF,IDX) = CPHFB(iaCPHF,IDX)/denom
        ENDDO
     ENDDO
  ENDDO

  ! Now the B0 array is complete.

END subroutine formcphfb

! Ed Brothers. November 18, 2002.
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

subroutine CPHFB4cnt(Ibas,Jbas,IIbas,JJbas)
  use allmod
  implicit double precision(a-h,o-z)

  dimension isame(4,8),deriv(4,3),icenter(4)
  logical :: same

  ! The purpose of the subroutine is to calculate an AO repulsion
  ! integral, determine it's contribution to an MO repulsion integral,
  ! and put it in the correct location in the CPHFB array, and then do the
  ! same with the integrals derivatives.

  ! First, calculate the integral.

  iA = ncenter(Ibas)
  iB = ncenter(Jbas)
  iC = ncenter(IIbas)
  iD = ncenter(JJbas)
  icenter(1)=iA
  icenter(2)=iB
  icenter(3)=iC
  icenter(4)=iD

  AOint=0.d0
  DO Icon=1,ncontract(Ibas)
     DO Jcon=1,ncontract(Jbas)
        DO IIcon=1,ncontract(IIbas)
           DO JJcon=1,ncontract(JJbas)
              AOint = AOint + &
                   dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                   *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas) &
                   *repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                   aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                   itype(1,Ibas), itype(2,Ibas), itype(3,Ibas), &
                   itype(1,Jbas), itype(2,Jbas), itype(3,Jbas), &
                   itype(1,IIbas),itype(2,IIbas),itype(3,IIbas), &
                   itype(1,JJbas),itype(2,JJbas),itype(3,JJbas), &
                   xyz(1,iA), xyz(2,iA),xyz(3,iA), &
                   xyz(1,iB), xyz(2,iB),xyz(3,iB), &
                   xyz(1,iC), xyz(2,iC),xyz(3,iC), &
                   xyz(1,iD), xyz(2,iD),xyz(3,iD))
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  ! Now we need to find how many times the AO integral appears by examining
  ! it's symmetry.  For an integral (ij|kl):
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

  ! Now we can start filling up the CPHFB array.
  ! Note we are first doing the term:
  ! - (Sum over k and l) Skl(1) [(ai|lk)-(ak|li)]

  ! Note k and l are both occupied occupied molecular orbitals.

  ! Do the alpha first.

  DO iAocc = 1,lastAocc
     DO iAocc2 = 1,lastAocc

        ! K and L are selected.  Find atom and direction of perturbation.


        DO Iatom = 1,natom
           DO Imomentum=1,3

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

              ! At this point we have the SKl value.  Now we need to loop over
              ! unique AO repulsions and and add the values into the array.
              ! we also need to find the location in the array defined by the Iatom
              ! Imomentum pair.

              ISTART=(Iatom-1)*3+Imomentum
              DO Iunique=1,8
                 IF (isame(1,Iunique) /= 0) THEN
                    ! Set up some dummy variables.

                    Ival = isame(1,Iunique)
                    Jval = isame(2,Iunique)
                    IIval = isame(3,Iunique)
                    JJval = isame(4,Iunique)

                    ! Loop over alpha pairs.

                    DO iAvirt = lastAocc+1,nbasis
                       DO iAocc3 = 1,lastAocc
                          iaCPHF = (iAvirt-lastAocc-1)*lastAocc + iAocc3
                          CPHFB(iaCPHF,Istart) = CPHFB(iaCPHF,Istart) &
                               -Skl*AOint*(CO(Ival,iAvirt)*CO(Jval,iAocc3)* &
                               CO(IIval,iAocc2)*CO(JJval,iAocc)-CO(Ival,iAvirt)* &
                               CO(Jval,iAocc)*CO(IIval,iAocc2)*CO(JJval,iAocc3))

                       ENDDO
                    ENDDO

                    ! Loop over beta pairs.

                    DO iBvirt = lastBocc+1,nbasis
                       DO iBocc3 = 1,lastBocc
                          iaCPHF=(iBvirt-lastBocc-1)*lastBocc + iBocc3+iBetastart
                          CPHFB(iaCPHF,Istart) = CPHFB(iaCPHF,Istart) &
                               -Skl*AOint*(COB(Ival,iBvirt)*COB(Jval,iBocc3)* &
                               CO(IIval,iAocc2)*CO(JJval,iAocc))
                       ENDDO
                    ENDDO
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  ! Thats a lot of loop closures.  From top to bottow we are closing
  ! iBocc3,iBvirt,The iunique if, iunique, Imomentum,Iatom,IAocc2,IAocc.

  ! Now we need to repeat the whole process for the beta kl pairs.

  DO iBocc = 1,lastBocc
     DO iBocc2 = 1,lastBocc

        ! K and L are selected.  Find atom and direction of perturbation.


        DO Iatom = 1,natom
           DO Imomentum=1,3

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

              ! At this point we have the SKl value.  Now we need to loop over
              ! unique AO repulsions and and add the values into the array.
              ! we also need to find the location in the array defined by the Iatom
              ! Imomentum pair.

              ISTART=(Iatom-1)*3+Imomentum
              DO Iunique=1,8
                 IF (isame(1,Iunique) /= 0) THEN
                    ! Set up some dummy variables.

                    Ival = isame(1,Iunique)
                    Jval = isame(2,Iunique)
                    IIval = isame(3,Iunique)
                    JJval = isame(4,Iunique)

                    ! Loop over beta pairs.

                    DO iBvirt = lastBocc+1,nbasis
                       DO iBocc3 = 1,lastBocc
                          iaCPHF =(iBvirt-lastBocc-1)*lastBocc + iBocc3+ibetastart
                          CPHFB(iaCPHF,Istart) = CPHFB(iaCPHF,Istart) &
                               -Skl*AOint*(COB(Ival,iBvirt)*COB(Jval,iBocc3)* &
                               COB(IIval,iBocc2)*COB(JJval,iBocc)-COB(Ival,iBvirt)* &
                               COB(Jval,iBocc)*COB(IIval,iBocc2)*COB(Jval,iBocc3))
                       ENDDO
                    ENDDO

                    ! Loop over alpha pairs.

                    DO iAvirt = lastAocc+1,nbasis
                       DO iAocc3 = 1,lastAocc
                          iaCPHF = (iAvirt-lastAocc-1)*lastAocc + iAocc3
                          CPHFB(iaCPHF,Istart) = CPHFB(iaCPHF,Istart) &
                               -Skl*AOint*(CO(Ival,iAvirt)*CO(Jval,iAocc3)* &
                               COB(IIval,iBocc2)*COB(JJval,iBocc))
                       ENDDO
                    ENDDO
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
     ENDDO
  ENDDO

  ! Thats a lot of loop closures.  From top to bottow we are closing
  ! iAocc3,iAvirt,The iunique if, iunique, Imomentum,Iatom,IBocc2,IBocc.

  ! Now we calculate the final term:  (For the case where a and i are alpha)

  ! CO(Ibas,A)*CO(Jbas,I)*DENSE(TOTAL)(IIbas,JJbas)*d/dy(IbasJbas|IIbasJJbas)
  ! - CO(Ibas,A)*CO(JJbas,I)*DENSE(A)(IIbas,Jbas)*d/dy(IbasJbas|IIbasJJbas)

  ! Now, if all the gaussians are on the same center, ther derivative is zero.
  ! Check this.

  same=icenter(1).eq.icenter(2)
  same=same.and.icenter(1).eq.icenter(3)
  same=same.and.icenter(1).eq.icenter(4)
  IF (same) return

  ! Otherwise, calculate the derivative.  This returns two arrays, one filled
  ! with derivatives and one filled with the center identities.  Note that
  ! this removes redundant centers, i.e. if this is a two center repulsion,
  ! only two slots are filled in the icenter and 6 slots in the deriv array.

  call CPHFB2Egrad(Ibas,Jbas,IIbas,JJbas,deriv,icenter)

  ! Now loop over atoms in icenter and momentums.  This will give us our
  ! Istart for the array.

  DO Iatom = 1,4
     IF (Icenter(Iatom) /= 0) THEN
        DO Imomentum=1,3
           currderiv=deriv(Iatom,Imomentum)
           Istart = (Icenter(Iatom)-1)*3 + Imomentum
           DO Iunique=1,8
              IF (isame(1,Iunique) /= 0) THEN

                 ! Set up some dummy variables.

                 Ival = isame(1,Iunique)
                 Jval = isame(2,Iunique)
                 IIval = isame(3,Iunique)
                 JJval = isame(4,Iunique)

                 ! A quick note about the density used below.  If the integral was (ij|kl),
                 ! we need the total density matrix element k,l and the alpha and beta
                 ! density elements k,j.

                 ! Why is this?  The proof is left to the reader.

                 ! Seriously, (ij|kl) is the exchange integral that occurrs with the
                 ! coulombic integral (il|kj), and the indices on the exchange elements
                 ! density matrix are always the last two of the coloumbic matrix.

                 IF (quick_method%unrst) THEN
                    DENSEKL = DENSE(IIval,JJval)+DENSEB(IIval,JJval)
                    DENSEAKJ = DENSE(IIval,Jval)
                    DENSEBKJ = DENSEB(IIval,Jval)
                 ELSE
                    DENSEKL = DENSE(IIval,JJval)
                    DENSEAKJ = DENSE(IIval,Jval)*.5d0
                    DENSEBKJ = DENSE(IIval,Jval)*.5d0
                 ENDIF

                 ! Loop over alpha pairs.

                 DO iAvirt = lastAocc+1,nbasis
                    DO iAocc = 1,lastAocc
                       iaCPHF =(iAvirt-lastAocc-1)*lastAocc + iAocc
                       CPHFB(iaCPHF,Istart) = CPHFB(iaCPHF,Istart) &
                            +currderiv*CO(Ival,iAvirt)* &
                            (CO(Jval,iAocc)*DENSEKL-CO(JJval,iAocc)*DENSEAKJ)
                    ENDDO
                 ENDDO

                 ! Loop over beta pairs.

                 DO iBvirt = lastBocc+1,nbasis
                    DO iBocc = 1,lastBocc
                       iaCPHF =(iBvirt-lastBocc-1)*lastBocc + iBocc+Ibetastart
                       CPHFB(iaCPHF,Istart) = CPHFB(iaCPHF,Istart) &
                            +currderiv*COB(Ival,iBvirt)* &
                            (COB(Jval,iBocc)*DENSEKL-COB(JJval,iBocc)*DENSEBKJ)
                    ENDDO
                 ENDDO

              ENDIF
           ENDDO
        ENDDO
     ENDIF
  ENDDO

  return
end subroutine CPHFB4cnt



! Ed Brothers. November 27,2002.
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

subroutine CPHFB2Egrad(Ibas,Jbas,IIbas,JJbas,deriv,icenter)
  use allmod
  implicit double precision(a-h,o-z)

  dimension deriv(4,3),icenter(4)
  dimension itype2(3,4)

  ! The purpose of this subroutine is to calculate the gradient of
  ! the 2-electron 4-center integrals used in forming the B array of the
  ! CPHF.

  ! Note that this is basically grad2elec, and could be used to replace
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

END subroutine CPHFB2Egrad

