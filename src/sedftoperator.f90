! Ed Brothers. November 27, 2001
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

    subroutine sedftoperator
    use allmod
    implicit double precision(a-h,o-z)

! The purpose of this subroutine is to form the operator matrix
! for an sedft calculation.  Note this operator matrix is symmetric.

! Blank out the array.

    do Ibas=1,nbasis
        do Jbas=Ibas,nbasis
            O(Jbas,Ibas)=0.d0
        enddo
    enddo


! do the one center terms first:

! There are two 1e- 1center terms.  They are the kinetic energy and
! a subset of the Nuclear attractions.

    do Iatm=1,natom
        Iatmfirst = Ifirst(Iatm)
        Iatmlast = Ilast(Iatm)
        xIatm = xyz(1,Iatm)
        yIatm = xyz(2,Iatm)
        zIatm = xyz(3,Iatm)
        chgIatm = chg(Iatm )
        ITyp = iattype(Iatm)
        do Ibas=Iatmfirst,Iatmlast
            do Jbas=Ibas,Iatmlast
                param=(EK1prm(itype(1,Jbas),itype(2,Jbas),itype(3,Jbas),ITyp)+ &
                EK1prm(itype(1,Ibas),itype(2,Ibas),itype(3,Ibas),ITyp)) &
                /2.d0
                param2=(At1prm(itype(1,Jbas),itype(2,Jbas),itype(3,Jbas),ITyp)+ &
                At1prm(itype(1,Ibas),itype(2,Ibas),itype(3,Ibas),ITyp)) &
                /2.d0
                do Icon=1,ncontract(ibas)
                    do Jcon=1,ncontract(jbas)

                        O(Jbas,Ibas)=O(Jbas,Ibas)+dcoeff(Jcon,Jbas)* &
                        dcoeff(Icon,Ibas)* &
                        (param* &
                        ekinetic(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                        itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                        itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                        xIatm,yIatm,zIatm,xIatm,yIatm,zIatm) &
                        +param2* &
                        attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                        itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                        itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                        xIatm,yIatm,zIatm,xIatm,yIatm,zIatm, &
                        xIatm,yIatm,zIatm,chgIatm))
                    enddo
                enddo
            enddo
        enddo
    enddo

! The next term is the two electron 1 center repulsion_prims.

    do Iatm=1,natom
        Iatmfirst = Ifirst(Iatm)
        Iatmlast = Ilast(Iatm)
        do I=Iatmfirst,iatmlast
        ! Set some variables to reduce access time for some of the more
        ! used quantities.

            xI = xyz(1,ncenter(I))
            yI = xyz(2,ncenter(I))
            zI = xyz(3,ncenter(I))
            itype1I=itype(1,I)
            itype2I=itype(2,I)
            itype3I=itype(3,I)
            DENSEII=DENSE(I,I)

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

            do J=I+1,iatmlast
            ! Set some variables to reduce access time for some of the more
            ! used quantities. (AGAIN)

                xJ = xyz(1,ncenter(J))
                yJ = xyz(2,ncenter(J))
                zJ = xyz(3,ncenter(J))
                itype1J=itype(1,J)
                itype2J=itype(2,J)
                itype3J=itype(3,J)
                DENSEJI=DENSE(J,I)
                DENSEJJ=DENSE(J,J)

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

            ! Find  all the (ii|ij) integrals.
                Ibas=I
                Jbas=I
                IIbas=I
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

                do K=J+1,iatmlast
                ! Set some variables to reduce access time for some of the more
                ! used quantities. (AGAIN)

                    xK = xyz(1,ncenter(K))
                    yK = xyz(2,ncenter(K))
                    zK = xyz(3,ncenter(K))
                    itype1K=itype(1,K)
                    itype2K=itype(2,K)
                    itype3K=itype(3,K)
                    DENSEKI=DENSE(K,I)
                    DENSEKJ=DENSE(K,J)
                    DENSEKK=DENSE(K,K)

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
                enddo

                do K=I+1,iatmlast-1
                    xK = xyz(1,ncenter(K))
                    yK = xyz(2,ncenter(K))
                    zK = xyz(3,ncenter(K))
                    itype1K=itype(1,K)
                    itype2K=itype(2,K)
                    itype3K=itype(3,K)

                    do L=K+1,iatmlast
                        xL = xyz(1,ncenter(L))
                        yL = xyz(2,ncenter(L))
                        zL = xyz(3,ncenter(L))
                        itype1L=itype(1,L)
                        itype2L=itype(2,L)
                        itype3L=itype(3,L)
                        DENSELK=DENSE(L,K)

                    ! Find the (ij|kl) integrals where j>i,k>i,l>k.
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
                    enddo
                enddo
            enddo
        enddo

    enddo


! The next 2 terms are the two center nuclear attractions.
! They are Ibas,Jbas,Jatm, where Ibas and Jbas are on Iatm, and
! Ibas,Jbas,Iatm, where Ibas is on Iatm.

    do Ibas=1,nbasis
        Iatm = ncenter(Ibas)
        Iatmlast = Ilast(Iatm)
        ITyp = iattype(Iatm)
        do Jbas = Ibas,Iatmlast
            param=(At2prm(itype(1,Jbas),itype(2,Jbas),itype(3,Jbas),ITyp)+ &
            At2prm(itype(1,Ibas),itype(2,Ibas),itype(3,Ibas),ITyp)) &
            /2.d0
            do Jatm=1,natom
                if (Jatm /= Iatm) then
                    do Icon=1,ncontract(ibas)
                        do Jcon=1,ncontract(jbas)

                            O(Jbas,Ibas)=O(Jbas,Ibas)+param* &
                            dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas)* &
                            attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                            itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                            itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                            xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                            xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                            xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)), &
                            xyz(1,Jatm),xyz(2,Jatm),xyz(3,Jatm), &
                            chg(Jatm))
                        enddo
                    enddo
                endif
            enddo
        enddo
        do Jbas = Iatmlast+1,nbasis
            Jatm=ncenter(Jbas)
            JTyp = iattype(Jatm)
            param=(Bndprm(itype(1,Jbas),itype(2,Jbas),itype(3,Jbas),JTyp)+ &
            Bndprm(itype(1,Ibas),itype(2,Ibas),itype(3,Ibas),ITyp)) &
            /2.d0
            do Icon=1,ncontract(ibas)
                do Jcon=1,ncontract(jbas)
                    O(Jbas,Ibas)=O(Jbas,Ibas)+param* &
                    dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas)* &
                    attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                    itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                    itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                    xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                    xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                    xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)), &
                    xyz(1,Jatm),xyz(2,Jatm),xyz(3,Jatm), &
                    chg(Jatm))
                    O(Jbas,Ibas)=O(Jbas,Ibas)+param* &
                    dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas)* &
                    attraction(aexp(Jcon,Jbas),aexp(Icon,Ibas), &
                    itype(1,Jbas),itype(2,Jbas),itype(3,Jbas), &
                    itype(1,Ibas),itype(2,Ibas),itype(3,Ibas), &
                    xyz(1,ncenter(Jbas)),xyz(2,ncenter(Jbas)), &
                    xyz(3,ncenter(Jbas)),xyz(1,ncenter(Ibas)), &
                    xyz(2,ncenter(Ibas)),xyz(3,ncenter(Ibas)), &
                    xyz(1,Iatm),xyz(2,Iatm),xyz(3,Iatm), &
                    chg(Iatm))
                enddo
            enddo
        enddo
    enddo



! Now we do the two center 2e- repulsion_prim terms.

    do Iatm=1,natom
        Iatmfirst = Ifirst(Iatm)
        Iatmlast = Ilast(Iatm)
        do Ibas = Iatmfirst,Iatmlast
            do Jbas = Ibas,Iatmlast
                do Jatm = Iatm+1,natom
                    Jatmfirst = Ifirst(Jatm)
                    Jatmlast = Ilast(Jatm)
                    do IIbas = Jatmfirst,Jatmlast
                        do JJbas = IIbas,Jatmlast
                            repint = 0.d0
                            do Icon=1,ncontract(ibas)
                                do Jcon=1,ncontract(jbas)
                                    do IIcon=1,ncontract(iibas)
                                        do JJcon=1,ncontract(jjbas)
                                            repint = repint+ &
                                            dcoeff(Jcon,Jbas)*dcoeff(Icon,Ibas) &
                                            *dcoeff(JJcon,JJbas)*dcoeff(IIcon,IIbas)* &
                                            (repulsion_prim(aexp(Icon,Ibas),aexp(Jcon,Jbas), &
                                            aexp(IIcon,IIbas),aexp(JJcon,JJbas), &
                                            itype(1,Ibas),itype(2,Ibas),Itype(3,Ibas), &
                                            itype(1,Jbas),itype(2,Jbas),Itype(3,Jbas), &
                                            itype(1,IIbas),itype(2,IIbas),Itype(3,IIbas), &
                                            itype(1,JJbas),itype(2,JJbas),Itype(3,JJbas), &
                                            xyz(1,iatm),xyz(2,iatm),xyz(3,iatm), &
                                            xyz(1,iatm),xyz(2,iatm),xyz(3,iatm), &
                                            xyz(1,jatm),xyz(2,jatm),xyz(3,jatm), &
                                            xyz(1,jatm),xyz(2,jatm),xyz(3,jatm)))
                                        enddo
                                    enddo
                                enddo
                            enddo
                            if (JJbas == IIbas) then
                                O(Jbas,Ibas)=O(Jbas,Ibas)+DENSE(JJbas,IIbas)*repint
                            else
                                O(Jbas,Ibas)=O(Jbas,Ibas)+2.d0*DENSE(JJbas,IIbas)*repint
                            endif
                            if (Jbas == Ibas) then
                                O(JJbas,IIbas)=O(JJbas,IIbas)+DENSE(Jbas,Ibas)*repint
                            else
                                O(JJbas,IIbas)=O(JJbas,IIbas)+2.d0*DENSE(Jbas,Ibas)*repint
                            endif
                        enddo
                    enddo
                enddo
            enddo
        enddo
    enddo

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

    do Ireg=1,iregion
        call gridform(iangular(Ireg))
        do Irad=iradial(ireg-1)+1,iradial(ireg)
            do Iatm=1,natom
                rad = radii(iattype(iatm))
                rad3 = rad*rad*rad
                do Iang=1,iangular(Ireg)
                    gridx=xyz(1,Iatm)+rad*RGRID(Irad)*XANG(Iang)
                    gridy=xyz(2,Iatm)+rad*RGRID(Irad)*YANG(Iang)
                    gridz=xyz(3,Iatm)+rad*RGRID(Irad)*ZANG(Iang)

                ! Next, calculate the weight of the grid point in the SSW scheme.  If
                ! the grid point has a zero weight, we can skip it.

                    weight=SSW(gridx,gridy,gridz,Iatm) &
                    *WTANG(Iang)*RWT(Irad)*rad3

                    if (weight < tol ) then
                        continue
                    else

                    ! Next, evaluate the densities at the grid point and the gradient
                    ! at that grid point.

                        call denspt(gridx,gridy,gridz,density,densityb,gax,gay,gaz, &
                        gbx,gby,gbz)

                        if (density < tol ) then
                            continue
                        else
                            aelec = weight*density+aelec
                            belec = weight*densityb+belec

                        ! This allows the calculation of the derivative of the functional
                        ! with regard to the density (dfdr), with regard to the alpha-alpha
                        ! density invariant (df/dgaa), and the alpha-beta density invariant.

                            call mpw91(density,gax,gay,gaz,gbx,gby,gbz, &
                            dfdr,dfdgaa,dfdgab)
                            call lyp(density,densityb,gax,gay,gaz,gbx,gby,gbz, &
                            dfdr2,dfdgaa2,dfdgab2)
                            dfdr = param7*dfdr+param8*dfdr2
                            dfdgaa = param7*dfdgaa + param8*dfdgaa2
                            dfdgab = param7*dfdgab + param8*dfdgab2

                        ! Calculate the first term in the dot product shown above,i.e.:
                        ! (2 df/dgaa Grad(rho a) + df/dgab Grad(rho b)) doT Grad(Phimu Phinu))

                            xdot=2.d0*dfdgaa*gax+dfdgab*gbx
                            ydot=2.d0*dfdgaa*gay+dfdgab*gby
                            zdot=2.d0*dfdgaa*gaz+dfdgab*gbz

                        ! Now loop over basis functions and compute the addition to the matrix
                        ! element.

                            do Ibas=1,nbasis
                                call pteval(gridx,gridy,gridz,phi,dphidx,dphidy, &
                                dphidz,Ibas)
                                quicktest = DABS(dphidx)+DABS(dphidy)+DABS(dphidz) &
                                +DABS(phi)
                                if (quicktest < tol ) then
                                    continue
                                else
                                    do Jbas=Ibas,nbasis
                                        call pteval(gridx,gridy,gridz,phi2,dphi2dx,dphi2dy, &
                                        dphi2dz,Jbas)
                                        quicktest = DABS(dphi2dx)+DABS(dphi2dy)+DABS(dphi2dz) &
                                        +DABS(phi2)
                                        if (quicktest < tol ) then
                                            continue
                                        else
                                            temp = phi*phi2
                                            tempgx = phi*dphi2dx + phi2*dphidx
                                            tempgy = phi*dphi2dy + phi2*dphidy
                                            tempgz = phi*dphi2dz + phi2*dphidz
                                            O(Jbas,Ibas)=O(Jbas,Ibas)+(temp*dfdr+ &
                                            xdot*tempgx+ydot*tempgy+zdot*tempgz)*weight
                                        endif
                                    enddo
                                endif
                            enddo
                        endif
                    endif
                enddo
            enddo
        enddo
    enddo

! Finally, copy lower diagonal to upper diagonal.

    do Ibas=1,nbasis
        do Jbas=Ibas,nbasis
            O(Ibas,Jbas) = O(Jbas,Ibas)
        enddo
    enddo

    END subroutine sedftoperator

