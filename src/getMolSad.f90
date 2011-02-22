! Ed Brothers. 11/26/01
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

    subroutine getmolsad()
    use allmod
    implicit double precision(a-h,o-z)

    logical :: present,HFsaved,UNRSTsaved,MPIsaved,ZMATsaved,divconsaved,coresaved
    integer imultsaved
    character(len=80) :: keywd
    character(len=20) :: tempstring
    

! Constants
!    pi=3.1415926535897932385d0
!    pito3half=5.568327996831707845284817982118835702014d0
    pi32 = pi**(1.5)

    istart = 1
    ifinal = 80
    ibasisstart = 1
    ibasisend = 80

    ! Save varibles first
    MPIsaved=bMPI
    HFsaved=quick_method%HF
    UNRSTsaved=quick_method%UNRST
    divconsaved=quick_method%divcon
    imultsaved=imult
    ZMATsaved=quick_method%ZMAT
    coresaved=quick_method%core
    
    ! Then give them new value
    bMPI=.false.    
    imult=0
    quick_method%HF=.true.
    quick_method%UNRST=.true.
    quick_method%ZMAT=.false.
    

    if (master) then
!        open(infile,file=inFileName,status='old')
!        read (infile,'(A80)') keywd
!        call upcase(keywd,80)
!        close(infile)
    
        call PrtAct(ioutfile,"Begin SAD initial guess")
        call rdword(basisdir,ibasisstart,ibasisend)
!AG 03/05/2007
        if (quick_method%ecp) call rdword(ecpdir,iecpstart,iecpend)
    
        !-------------------------------------------
        ! First, find atom type and initialize
        !-------------------------------------------

        natom=1
        idivcon=1
        do iitemp=1,iatomtype
            write(ioutfile,'(" For Atom Kind = ",i4)') iitemp
            do i=1,90 
                if(symbol(i).eq.atomxiao(iitemp))then
                    if(mod(i,2).eq.0)then
                        imult=1
                    else
                        imult=2
                    endif
                    if(symbol(i).eq.'N ')imult=4
                    if(symbol(i).eq.'O ')imult=3
                    if(symbol(i).eq.'C ')imult=3
                    if(symbol(i).eq.'S ')imult=3
                    if(symbol(i).eq.'P ')imult=4
                    chg(1)=i
                    iiatom=i
                    iattype(1)=i
                    write(ioutfile,'(" ELEMENT = ",a)') symbol(i)
                endif
            enddo
    
            if (ibasis == 2 .AND. .NOT. quick_method%core) then
                write (ioutfile,'("VO BASIS FUNCTIONS ARE USED WITH THE CORE APPROXIMATION.")')
                quick_method%core=.true.
            endif
    
            tol=1.d-12
            if (imult /= 1) quick_method%UNRST= .TRUE. 
            nelec = iiatom
            do I=1,3
                xyz(I,1) = 0.0d0
            enddo


            if (quick_method%DFT .OR. quick_method%SEDFT) then
                ! Xiao HE 01/09/07 SG1 SG0 grid
                if (quick_method%iSG.eq.1) then
                    itemp=50
                    do I=1,itemp
                        RGRID(I)=(I**2.d0)/dble((itemp+1-I)*(itemp+1-I))
                        RWT(I)=2.d0*dble(itemp+1)*(dble(I)**5.d0) &
                                *dble(itemp+1-I)**(-7.d0)
                    enddo
                else
                    continue
                endif
            endif
    
            if (quick_method%core) then
                do Iatm=1,natom
                    if (iattype(Iatm) >= 3) then
                        chg(Iatm)=chg(Iatm)-2.d0
                        nelec=nelec-2
                    endif
                    if (iattype(Iatm) >= 11) then
                        chg(Iatm)=chg(Iatm)-8.d0
                        nelec=nelec-8
                    endif
                    if (iattype(Iatm) >= 18) write (ioutfile,'("ATOM OUT OF RANGE FOR CORE")')
                enddo
            endif
        
            !-------------------------------------------
            ! Alessandro GENONI 03/05/2007
            ! Only for ECP calculations:
            ! * Allocate arrays whose dimensions depend on NATOM (allocateatoms_ecp)
            ! * Read the Effective Core Potentials (ECPs), modify the atomic charges 
            !   and the total number of electrons (readecp)
            !-------------------------------------------
        
            if (quick_method%ecp) then
                call allocateatoms_ecp
                call readecp
            END if

            ! Modify the electron numbers based on charge and multiplicity.
            nelec = nelec - molchg
            xmulttest = mod(dble(nelec),2.d0)

            if (xmulttest /= 0 .AND. .NOT. quick_method%unrst) then
                write (ioutfile,'("WARNING: UNPAIRED ELECTRONS REQUIRE UNRESTRICTED CALCULATIONS.")')
                quick_method%UNRST=.true.
            endif
            
            if (imult /= 1 .AND. .NOT. quick_method%unrst) then 
                write (ioutfile,'("WARNING: HIGHER MULTIPLICITIES REQUIRE  UNRESTRICTED CALCULATIONS.")')
            endif

            ! If unrestricted obtain alpha and beta electrons
            if (quick_method%unrst) then
                nelecb = nelec
                nelec = 0

                do WHILE (nelec.lt.nelecb)
                    nelecb = nelecb-1
                    nelec = nelec +1
                enddo

                if (imult == 2 .AND. nelec-1 /= nelecb) write (ioutfile,'("WARNING: INCORRECT NUMBER OF ELECTRONS FOR A doUBLET.")')
                if (imult == 3) then
                    nelec = nelec+1
                    nelecb = nelecb - 1
                    if (nelec-2 /= nelecb) write (ioutfile,'("WARNING: INCORRECT NUMBER OF ELECTRONS FOR A TRIPLET.")')
                endif
                if (imult == 4) then
                    nelec = nelec+1
                    nelecb = nelecb - 1
                    if (nelec-3 /= nelecb) write (ioutfile,'("WARNING: INCORRECT NUMBER OF ELECTRONS FOR A QUADRUPLET.")')
                endif
                write (ioutfile,'("NUMBER OF ALPHA ELECTRONS = ",I4)')nelec
                write (ioutfile,'("NUMBER OF BETA ELECTRONS  = ",I4)')nelecb
            else
                write (ioutfile,'(/"NUMBER OF ELECTRONS = ",I4)')nelec
            endif
    
            !-------------------------------------------
            ! At this point we have the positions and identities of the atoms. We also
            ! have the number of electrons. Now we must assign basis functions. This
            ! is done in a subroutine.
            !-------------------------------------------
            nsenhai=1
            call readbasis(nsenhai,0,0,0,0)
            atombasis(iitemp)=nbasis
            write (ioutfile,'("BASIS FUNCTIONS = ",I4)') nbasis

            ! Include the normalization constant in the coefficient.
            do Jbas=1,nbasis
                do Jcon=1,ncontract(jbas)
                    dcoeff(Jcon,Jbas)=dcoeff(Jcon,Jbas) *xnorm(aexp(Jcon,Jbas), &
                    itype(1,Jbas),itype(2,Jbas),itype(3,Jbas))
                enddo
            enddo
            do Ibas=1,nbasis
                dconew = 0.d0
                nxponent=-itype(1,Ibas)-itype(2,Ibas)-itype(3,Ibas)
                xponent=-1.5d0+dble(nxponent)
                do Icon1 = 1,ncontract(Ibas)
                    do Icon2 = 1,ncontract(Ibas)
                        dconew = dconew + dcoeff(Icon1,Ibas)*dcoeff(Icon2,Ibas) &
                        *(aexp(Icon1,Ibas)+aexp(Icon2,Ibas))**xponent
                    enddo
                enddo
                gamma=1.d0
                do L=1,itype(1,Ibas)
                    gamma = gamma * (dble(itype(1,Ibas) - L) + .5d0)
                enddo
                do L=1,itype(2,Ibas)
                    gamma = gamma * (dble(itype(2,Ibas) - L) + .5d0)
                enddo
                do L=1,itype(3,Ibas)
                    gamma = gamma * (dble(itype(3,Ibas) - L) + .5d0)
                enddo
                dconew = (dconew*gamma*pi32)**(-.5d0)
                do Icon1 = 1,ncontract(Ibas)
                    dcoeff(Icon1,Ibas) = dconew*dcoeff(Icon1,Ibas)
                enddo
            enddo

            !-------------------------------------------
            ! Alessandro GENONI 03/21/2007
            ! Store the normalized primitives coefficients for ECP calculations
            !-------------------------------------------
            if (quick_method%ecp) then
                iicont=0
                icontb=1
                do i=1,nshell
                    do j=1,kprim(i)
                        iicont=iicont+1
                        eta(iicont)=dcoeff(j,icontb)
                    end do
                    icontb=icontb+ktype(i)
                end do
            end if

            !-------------------------------------------
            ! Now that the basis functions have been normalized, calculate the radius
            ! of the sphere of basis function signifigance. (See Stratmann,Scuseria,
            ! and Frisch, Chem. Phys. Lett., 257, 1996, page 213-223 Section 5.)
            ! Also, the radius of the sphere comes from the spherical average of
            ! the basis function, from Perez-Jorda and Yang, Chem. Phys. Lett., 241,
            ! 1995, pg 469-76.
            ! The spherical average of a gaussian function is:

            ! (1 + 2 L)/4  (3 + 2 L)/4  L
            ! 2            a            r
            ! ave  = ---------------------------------
            ! 2
            ! a r                      3
            ! E     Sqrt[Pi] Sqrt[Gamma[- + L]]
            ! 2
            ! where a is the most diffuse (smallest) orbital exponent and L is the
            ! sum of the angular momentum exponents.  This code finds the r value where
            ! the average is the signifigance threshold (signif) and this r value is
            ! called the target below. Rearranging gives us:

            ! -(1 + 2 L)/4   -(3 + 2 L)/4                           3
            !r^L E^-ar^2= 2               a           Sqrt[Pi] signif Sqrt[Gamma[- + L]]
            ! 2
            ! Which is our function to work with.
            
            if (quick_method%DFT .OR. quick_method%SEDFT) then
                write (ioutfile, &
                    & '(/"RADII OF SIGNifIGANCE FOR THE BASIS FUNCTIONS")')

                do Ibas=1,nbasis
                
                    ! Find the minimum gaussian exponent.
                    amin=10.D10
                    do Icon=1,ncontract(Ibas)
                        amin=min(amin,aexp(Icon,Ibas))
                    enddo

                    ! calculate L.
                    L = itype(1,Ibas)+ itype(2,Ibas)+ itype(3,Ibas)

                    ! calculate 2 Pi Gamma[L+3/2]
                    ! Remember that Gamma[i+1/2]=(i-1+1/2) Gamma[i-1+1/2] until you get to
                    ! Gamma[1/2] = Sqrt[Pi]
                    gamma=1.d0
                    do i=1,L+1
                        gamma = gamma * (dble(L+1-i) + .5d0)
                    enddo
                    gamma2pi=gamma*11.13665599366341569

                    ! Now put it all together to get the target value.
                    target = signif* &
                        (((2.d0*amin)**(dble(L)+1.5))/gamma2pi)**(-.5d0)

                    ! Now search to find the correct radial value.
                    stepsize=1.d0
                    radial=0.d0
                    do WHILE (stepsize.gt.1.d-4)
                        radial=radial+stepsize
                        current=Dexp(-amin*radial*radial)*radial**(dble(L))
                        if (current < target) then
                            radial=radial-stepsize
                            stepsize=stepsize/10.d0
                        endif
                    enddo

                    ! Store the square of the radii of signifigance as this is what the
                    ! denisty calculator works in.
                    sigrad2(Ibas)=radial*radial
                    write (ioutfile,'(I4,7x,F12.6)') Ibas,radial
                enddo
            endif

            ! Initialize Density arrays. Create initial density matrix guess.

            do Ibas=1,nbasis
                do Jbas=1,nbasis
                    DENSE(Jbas,Ibas)=0.d0
                    DENSEB(Jbas,Ibas)=0.d0
                enddo
            enddo
            present = .false.
            if (quick_method%readdmx) inquire (file=dmxfilename,exist=present)
            if (present) then
                open(idmxfile,file=dmxfilename,status='old')
                read (idmxfile,'(A80)') keywd
                do WHILE (index(keywd,'END').eq.0.and.index(keywd,'BETA').eq.0.)
                    read (keywd,'(I4,I4,3x,E15.8)') iread,jread,temp
                    DENSE(iread,jread)=temp
                    DENSE(jread,iread)=temp
                    read (idmxfile,'(A80)') keywd
                enddo
                if (index(keywd,'BETA') /= 0) read (idmxfile,'(A80)') keywd
                do WHILE (index(keywd,'END').eq.0)
                    read (keywd,'(I4,I4,3x,E15.8)') iread,jread,temp
                    DENSEB(iread,jread)=temp
                    DENSEB(jread,iread)=temp
                read (idmxfile,'(A80)') keywd
            enddo
            close(idmxfile)
            
            if (DENSEB(1,1) == 0.d0 .AND. quick_method%unrst) then
                write (ioutfile, &
                & '(/,"CONVERTING RESTRICTED DENSITY TO UNRESTRICTED")')
                do I=1,nbasis
                    do J =1,nbasis
                        DENSE(J,I) = DENSE(J,I)/2.d0
                        DENSEB(J,I) = DENSE(J,I)
                    enddo
                enddo
            endif
        else
            ! Initial Guess
            diagelement=dble(nelec)/dble(nbasis)
            diagelementb=dble(nelecb)/dble(nbasis)+1.d-8
            do I=1,nbasis
                DENSE(I,I)=diagelement
                DENSEB(I,I)=diagelementb
            enddo
        endif

        ! From SCF calculation to get initial density guess
        if(atomxiao(iitemp).ne.'ZN')then ! If not ZN
            call getenergy(failed)
            do i=1,nbasis
                do j=1,nbasis
                    atomdens(iitemp,i,j)=DENSE(i,j)+DENSEB(i,j)
                enddo
            enddo
        else
            open(213,file='znsad.txt')  !Read Zn
            do i=1,39
                do j=1,39
                    read(213,*)ixiao1,ixiao2,xiaotemp
                    atomdens(iitemp,ixiao1,ixiao2)=xiaotemp
                enddo
            enddo
            close(213)
        endif
        call deallocateall
    enddo
    
    call PrtAct(ioutfile,"Finish SAD initial guess")
    endif
    
    quick_method%HF=HFsaved
    quick_method%UNRST=UNRSTsaved
    imult=imultsaved
    quick_method%divcon=divconsaved
    bMPI=MPIsaved
    quick_method%ZMAT=ZMATsaved
    quick_method%core=coresaved
    return
    end subroutine getmolsad

