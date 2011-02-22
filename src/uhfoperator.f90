! Ed Brothers. December 20, 2001
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

    subroutine uhfoperatora
    use allmod
    implicit double precision(a-h,o-z)

! The purpose of this subroutine is to form the operator matrix
! for a full Hartree-Fock calculation, i.e. the Fock matrix.  The
! Fock matrix is as follows:

! O(I,J) =  F(I,J) = KE(I,J) + IJ attraction to each atom + repulsion_prim
! with each possible basis  - 1/2 exchange with each
! possible basis.

! Note that the Fock matrix is symmetric.

    do Ibas=1,nbasis
        do Jbas=Ibas,nbasis
            O(Jbas,Ibas) = 0.d0
            do Icon=1,ncontract(ibas)
                do Jcon=1,ncontract(jbas)

                ! The first part is the kinetic energy.

                    O(Jbas,Ibas)=O(Jbas,Ibas)+ &
                    dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas)* &
                    ekinetic(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                    itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                    itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                    xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                    xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                    xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))

                ! Next is a loop over atoms to contruct the attraction terms.

                    do iatom = 1,natom
                        O(Jbas,Ibas)=O(Jbas,Ibas)+ &
                        dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas)* &
                        attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                        itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                        itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                        xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                        xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                        xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)), &
                        xyz(1,iatom),xyz(2,iatom),xyz(3,iatom), &
                        chg(iatom))
                    enddo
                enddo
            enddo
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
! The next two terms define the two electron part.
    do I=1,nbasis
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

    ! do all the (ii|ii) integrals.
        Ibas=I
        Jbas=I
        IIbas=I
        JJbas=I
        repint=0.d0
        do Icon=1,ncontract(ibas)
            do Jcon=1,ncontract(jbas)
                do IIcon=1,ncontract(iibas)
                    do JJcon=1,ncontract(jjbas)
                        repint = repint+ &
                        dcoeff(Icon,Ibas)*dcoeff(Jcon,Jbas) &
                        *dcoeff(IIcon,IIbas)*dcoeff(JJcon,JJbas)* &
                        (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                        aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                        itype1I,itype2I,itype3I,itype1I,itype2I,itype3I, &
                        itype1I,itype2I,itype3I,itype1I,itype2I,itype3I, &
                        xI,yI,zI,xI,yI,zI,xI,yI,zI,xI,yI,zI))
                    enddo
                enddo
            enddo
        enddo
        O(I,I) = O(I,I)+DENSEII*repint
        O(I,I) = O(I,I)-DENSEIIX*repint

        do J=I+1,nbasis
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
            do Icon=1,ncontract(ibas)
                do Jcon=1,ncontract(jbas)
                    do IIcon=1,ncontract(iibas)
                        do JJcon=1,ncontract(jjbas)
                            repint = repint+ &
                            dcoeff(Icon,Ibas)*dcoeff(Jcon,Jbas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                            (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype1I,itype2I,itype3I,itype1I,itype2I,itype3I, &
                            itype1J,itype2J,itype3J,itype1J,itype2J,itype3J, &
                            xI,yI,zI,xI,yI,zI,xJ,yJ,zJ,xJ,yJ,zJ))
                        enddo
                    enddo
                enddo
            enddo
            O(I,I) = O(I,I)+DENSEJJ*repint
            O(J,J) = O(J,J)+DENSEII*repint
            O(J,I) = O(J,I)-DENSEJIX*repint

        ! Find  all the (ij|jj) integrals.
            Ibas=I
            Jbas=J
            IIbas=J
            JJbas=J
            repint=0.d0
            do Icon=1,ncontract(ibas)
                do Jcon=1,ncontract(jbas)
                    do IIcon=1,ncontract(iibas)
                        do JJcon=1,ncontract(jjbas)
                            repint = repint+ &
                            dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                            (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype1I,itype2I,itype3I,itype1J,itype2J,itype3J, &
                            itype1J,itype2J,itype3J,itype1J,itype2J,itype3J, &
                            xI,yI,zI,xJ,yJ,zJ,xJ,yJ,zJ,xJ,yJ,zJ))

                        enddo
                    enddo
                enddo
            enddo
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
            do Icon=1,ncontract(ibas)
                do Jcon=1,ncontract(jbas)
                    do IIcon=1,ncontract(iibas)
                        do JJcon=1,ncontract(jjbas)
                            repint = repint+ &
                            dcoeff(Icon,Ibas)*dcoeff(Jcon,Jbas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                            (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype1I,itype2I,itype3I,itype1I,itype2I,itype3I, &
                            itype1I,itype2I,itype3I,itype1J,itype2J,itype3J, &
                            xI,yI,zI,xI,yI,zI,xI,yI,zI,xJ,yJ,zJ))

                        enddo
                    enddo
                enddo
            enddo
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
            do Icon=1,ncontract(ibas)
                do Jcon=1,ncontract(jbas)
                    do IIcon=1,ncontract(iibas)
                        do JJcon=1,ncontract(jjbas)
                            repint = repint+ &
                            dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                            (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype1I,itype2I,itype3I,itype1J,itype2J,itype3J, &
                            itype1I,itype2I,itype3I,itype1J,itype2J,itype3J, &
                            xI,yI,zI,xJ,yJ,zJ,xI,yI,zI,xJ,yJ,zJ))

                        enddo
                    enddo
                enddo
            enddo
            O(J,I) = O(J,I)+2.d0*DENSEJI*repint
            O(J,J) = O(J,J)-DENSEIIX*repint
            O(I,I) = O(I,I)-DENSEJJX*repint
            O(J,I) = O(J,I)-DENSEJIX*repint

            do K=J+1,nbasis
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
                do Icon=1,ncontract(ibas)
                    do Jcon=1,ncontract(jbas)
                        do IIcon=1,ncontract(iibas)
                            do JJcon=1,ncontract(jjbas)
                                repint = repint+ &
                                dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                                *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                                (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                                aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                                itype1I,itype2I,itype3I,itype1J,itype2J,itype3J, &
                                itype1I,itype2I,itype3I,itype1K,itype2K,itype3K, &
                                xI,yI,zI,xJ,yJ,zJ,xI,yI,zI,xK,yK,zK))

                            enddo
                        enddo
                    enddo
                enddo
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
                do Icon=1,ncontract(ibas)
                    do Jcon=1,ncontract(jbas)
                        do IIcon=1,ncontract(iibas)
                            do JJcon=1,ncontract(jjbas)
                                repint = repint+ &
                                dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                                *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                                (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                                aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                                itype1I,itype2I,itype3I,itype1J,itype2J,itype3J, &
                                itype1K,itype2K,itype3K,itype1K,itype2K,itype3K, &
                                xI,yI,zI,xJ,yJ,zJ,xK,yK,zK,xK,yK,zK))

                            enddo
                        enddo
                    enddo
                enddo
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
                do Icon=1,ncontract(ibas)
                    do Jcon=1,ncontract(jbas)
                        do IIcon=1,ncontract(iibas)
                            do JJcon=1,ncontract(jjbas)
                                repint = repint+ &
                                dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                                *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                                (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                                aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                                itype1I,itype2I,itype3I,itype1K,itype2K,itype3K, &
                                itype1J,itype2J,itype3J,itype1J,itype2J,itype3J, &
                                xI,yI,zI,xK,yK,zK,xJ,yJ,zJ,xJ,yJ,zJ))

                            enddo
                        enddo
                    enddo
                enddo
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
                do Icon=1,ncontract(ibas)
                    do Jcon=1,ncontract(jbas)
                        do IIcon=1,ncontract(iibas)
                            do JJcon=1,ncontract(jjbas)
                                repint = repint+ &
                                dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                                *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                                (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                                aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                                itype1I,itype2I,itype3I,itype1I,itype2I,itype3I, &
                                itype1J,itype2J,itype3J,itype1K,itype2K,itype3K, &
                                xI,yI,zI,xI,yI,zI,xJ,yJ,zJ,xK,yK,zK))

                            enddo
                        enddo
                    enddo
                enddo
                O(K,J) = O(K,J)+DENSEII*repint
                O(I,I) = O(I,I)+2.d0*DENSEKJ*repint
                O(J,I) = O(J,I)-DENSEKIX*repint
                O(K,I) = O(K,I)-DENSEJIX*repint
            enddo

            do K=I+1,nbasis-1
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

                do L=K+1,nbasis
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
                    do Icon=1,ncontract(ibas)
                        do Jcon=1,ncontract(jbas)
                            do IIcon=1,ncontract(iibas)
                                do JJcon=1,ncontract(jjbas)
                                    repint = repint+ &
                                    dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                                    *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                                    (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                                    aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                                    itype1I,itype2I,itype3I,itype1J,itype2J,itype3J, &
                                    itype1K,itype2K,itype3K,itype1L,itype2L,itype3L, &
                                    xI,yI,zI,xJ,yJ,zJ,xK,yK,zK,xL,yL,zL))


                                enddo
                            enddo
                        enddo
                    enddo
                    O(J,I) = O(J,I)+2.d0*DENSELK*repint
                    O(L,K) = O(L,K)+2.d0*DENSEJI*repint
                    O(K,I) = O(K,I)-DENSELJX*repint
                    O(L,I) = O(L,I)-DENSEKJX*repint
                    if (J == K) then
                        O(J,K) = O(J,K)-2.d0*DENSELIX*repint
                        O(J,L) = O(J,L)-DENSEKIX*repint
                        O(L,J) = O(L,J)-DENSEKIX*repint
                    elseif (J == L) then
                        O(J,L) = O(J,L)-2.d0*DENSEKIX*repint
                        O(J,K) = O(J,K)-DENSELIX*repint
                        O(K,J) = O(K,J)-DENSELIX*repint
                    else
                        O(J,K) = O(J,K)-DENSELIX*repint
                        O(J,L) = O(J,L)-DENSEKIX*repint
                        O(K,J) = O(K,J)-DENSELIX*repint
                        O(L,J) = O(L,J)-DENSEKIX*repint
                    endif
                enddo
            enddo
        enddo
    enddo

    do Ibas=1,nbasis
        do Jbas=Ibas+1,nbasis
            O(Ibas,Jbas) = O(Jbas,Ibas)
        enddo
    enddo

    END subroutine uhfoperatora

! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Ed Brothers. December 20, 2001
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

    subroutine uhfoperatorb
    use allmod
    implicit double precision(a-h,o-z)

! The purpose of this subroutine is to form the operator matrix
! for the beta portion of the unrestricted HF erquation.
! Fock matrix is as follows:

! O(I,J) =  F(I,J) = KE(I,J) + IJ attraction to each atom + repulsion_prim
! with each possible basis  - 1/2 exchange with each
! possible basis.

! Note that the Fock matrix is symmetric.

    do Ibas=1,nbasis
        do Jbas=Ibas,nbasis
            O(Jbas,Ibas) = 0.d0
            do Icon=1,ncontract(ibas)
                do Jcon=1,ncontract(jbas)

                ! The first part is the kinetic energy.

                    O(Jbas,Ibas)=O(Jbas,Ibas)+ &
                    dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas)* &
                    ekinetic(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                    itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                    itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                    xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                    xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                    xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)))

                ! Next is a loop over atoms to contruct the attraction terms.

                    do iatom = 1,natom
                        O(Jbas,Ibas)=O(Jbas,Ibas)+ &
                        dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas)* &
                        attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                        itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                        itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                        xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                        xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                        xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)), &
                        xyz(1,iatom),xyz(2,iatom),xyz(3,iatom), &
                        chg(iatom))
                    enddo
                enddo
            enddo
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
! The next two terms define the two electron part.

    do I=1,nbasis

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

    ! do all the (ii|ii) integrals.
        Ibas=I
        Jbas=I
        IIbas=I
        JJbas=I
        repint=0.d0
        do Icon=1,ncontract(ibas)
            do Jcon=1,ncontract(jbas)
                do IIcon=1,ncontract(iibas)
                    do JJcon=1,ncontract(jjbas)
                        repint = repint+ &
                        dcoeff(Icon,Ibas)*dcoeff(Jcon,Jbas) &
                        *dcoeff(IIcon,IIbas)*dcoeff(JJcon,JJbas)* &
                        (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                        aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                        itype1I,itype2I,itype3I,itype1I,itype2I,itype3I, &
                        itype1I,itype2I,itype3I,itype1I,itype2I,itype3I, &
                        xI,yI,zI,xI,yI,zI,xI,yI,zI,xI,yI,zI))
                    enddo
                enddo
            enddo
        enddo
        O(I,I) = O(I,I)+DENSEII*repint
        O(I,I) = O(I,I)-DENSEIIX*repint

        do J=I+1,nbasis
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
            do Icon=1,ncontract(ibas)
                do Jcon=1,ncontract(jbas)
                    do IIcon=1,ncontract(iibas)
                        do JJcon=1,ncontract(jjbas)
                            repint = repint+ &
                            dcoeff(Icon,Ibas)*dcoeff(Jcon,Jbas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                            (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype1I,itype2I,itype3I,itype1I,itype2I,itype3I, &
                            itype1J,itype2J,itype3J,itype1J,itype2J,itype3J, &
                            xI,yI,zI,xI,yI,zI,xJ,yJ,zJ,xJ,yJ,zJ))
                        enddo
                    enddo
                enddo
            enddo
            O(I,I) = O(I,I)+DENSEJJ*repint
            O(J,J) = O(J,J)+DENSEII*repint
            O(J,I) = O(J,I)-DENSEJIX*repint

        ! Find  all the (ij|jj) integrals.
            Ibas=I
            Jbas=J
            IIbas=J
            JJbas=J
            repint=0.d0
            do Icon=1,ncontract(ibas)
                do Jcon=1,ncontract(jbas)
                    do IIcon=1,ncontract(iibas)
                        do JJcon=1,ncontract(jjbas)
                            repint = repint+ &
                            dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                            (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype1I,itype2I,itype3I,itype1J,itype2J,itype3J, &
                            itype1J,itype2J,itype3J,itype1J,itype2J,itype3J, &
                            xI,yI,zI,xJ,yJ,zJ,xJ,yJ,zJ,xJ,yJ,zJ))

                        enddo
                    enddo
                enddo
            enddo
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
            do Icon=1,ncontract(ibas)
                do Jcon=1,ncontract(jbas)
                    do IIcon=1,ncontract(iibas)
                        do JJcon=1,ncontract(jjbas)
                            repint = repint+ &
                            dcoeff(Icon,Ibas)*dcoeff(Jcon,Jbas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                            (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype1I,itype2I,itype3I,itype1I,itype2I,itype3I, &
                            itype1I,itype2I,itype3I,itype1J,itype2J,itype3J, &
                            xI,yI,zI,xI,yI,zI,xI,yI,zI,xJ,yJ,zJ))

                        enddo
                    enddo
                enddo
            enddo
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
            do Icon=1,ncontract(ibas)
                do Jcon=1,ncontract(jbas)
                    do IIcon=1,ncontract(iibas)
                        do JJcon=1,ncontract(jjbas)
                            repint = repint+ &
                            dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                            (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                            itype1I,itype2I,itype3I,itype1J,itype2J,itype3J, &
                            itype1I,itype2I,itype3I,itype1J,itype2J,itype3J, &
                            xI,yI,zI,xJ,yJ,zJ,xI,yI,zI,xJ,yJ,zJ))

                        enddo
                    enddo
                enddo
            enddo
            O(J,I) = O(J,I)+2.d0*DENSEJI*repint
            O(J,J) = O(J,J)-DENSEIIX*repint
            O(I,I) = O(I,I)-DENSEJJX*repint
            O(J,I) = O(J,I)-DENSEJIX*repint

            do K=J+1,nbasis
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
                do Icon=1,ncontract(ibas)
                    do Jcon=1,ncontract(jbas)
                        do IIcon=1,ncontract(iibas)
                            do JJcon=1,ncontract(jjbas)
                                repint = repint+ &
                                dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                                *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                                (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                                aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                                itype1I,itype2I,itype3I,itype1J,itype2J,itype3J, &
                                itype1I,itype2I,itype3I,itype1K,itype2K,itype3K, &
                                xI,yI,zI,xJ,yJ,zJ,xI,yI,zI,xK,yK,zK))

                            enddo
                        enddo
                    enddo
                enddo
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
                do Icon=1,ncontract(ibas)
                    do Jcon=1,ncontract(jbas)
                        do IIcon=1,ncontract(iibas)
                            do JJcon=1,ncontract(jjbas)
                                repint = repint+ &
                                dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                                *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                                (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                                aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                                itype1I,itype2I,itype3I,itype1J,itype2J,itype3J, &
                                itype1K,itype2K,itype3K,itype1K,itype2K,itype3K, &
                                xI,yI,zI,xJ,yJ,zJ,xK,yK,zK,xK,yK,zK))

                            enddo
                        enddo
                    enddo
                enddo
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
                do Icon=1,ncontract(ibas)
                    do Jcon=1,ncontract(jbas)
                        do IIcon=1,ncontract(iibas)
                            do JJcon=1,ncontract(jjbas)
                                repint = repint+ &
                                dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                                *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                                (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                                aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                                itype1I,itype2I,itype3I,itype1K,itype2K,itype3K, &
                                itype1J,itype2J,itype3J,itype1J,itype2J,itype3J, &
                                xI,yI,zI,xK,yK,zK,xJ,yJ,zJ,xJ,yJ,zJ))

                            enddo
                        enddo
                    enddo
                enddo
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
                do Icon=1,ncontract(ibas)
                    do Jcon=1,ncontract(jbas)
                        do IIcon=1,ncontract(iibas)
                            do JJcon=1,ncontract(jjbas)
                                repint = repint+ &
                                dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                                *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                                (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                                aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                                itype1I,itype2I,itype3I,itype1I,itype2I,itype3I, &
                                itype1J,itype2J,itype3J,itype1K,itype2K,itype3K, &
                                xI,yI,zI,xI,yI,zI,xJ,yJ,zJ,xK,yK,zK))

                            enddo
                        enddo
                    enddo
                enddo
                O(K,J) = O(K,J)+DENSEII*repint
                O(I,I) = O(I,I)+2.d0*DENSEKJ*repint
                O(J,I) = O(J,I)-DENSEKIX*repint
                O(K,I) = O(K,I)-DENSEJIX*repint
            enddo

            do K=I+1,nbasis-1
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

                do L=K+1,nbasis
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
                    do Icon=1,ncontract(ibas)
                        do Jcon=1,ncontract(jbas)
                            do IIcon=1,ncontract(iibas)
                                do JJcon=1,ncontract(jjbas)
                                    repint = repint+ &
                                    dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                                    *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                                    (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                                    aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                                    itype1I,itype2I,itype3I,itype1J,itype2J,itype3J, &
                                    itype1K,itype2K,itype3K,itype1L,itype2L,itype3L, &
                                    xI,yI,zI,xJ,yJ,zJ,xK,yK,zK,xL,yL,zL))


                                enddo
                            enddo
                        enddo
                    enddo
                    O(J,I) = O(J,I)+2.d0*DENSELK*repint
                    O(L,K) = O(L,K)+2.d0*DENSEJI*repint
                    O(K,I) = O(K,I)-DENSELJX*repint
                    O(L,I) = O(L,I)-DENSEKJX*repint
                    if (J == K) then
                        O(J,K) = O(J,K)-2.d0*DENSELIX*repint
                        O(J,L) = O(J,L)-DENSEKIX*repint
                        O(L,J) = O(L,J)-DENSEKIX*repint
                    elseif (J == L) then
                        O(J,L) = O(J,L)-2.d0*DENSEKIX*repint
                        O(J,K) = O(J,K)-DENSELIX*repint
                        O(K,J) = O(K,J)-DENSELIX*repint
                    else
                        O(J,K) = O(J,K)-DENSELIX*repint
                        O(J,L) = O(J,L)-DENSEKIX*repint
                        O(K,J) = O(K,J)-DENSELIX*repint
                        O(L,J) = O(L,J)-DENSEKIX*repint
                    endif
                enddo
            enddo
        enddo
    enddo

    do Ibas=1,nbasis
        do Jbas=Ibas+1,nbasis
            O(Ibas,Jbas) = O(Jbas,Ibas)
        enddo
    enddo

    END subroutine uhfoperatorb

