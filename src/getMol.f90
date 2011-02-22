!*******************************************************
! getAtoms(iAtomType)
!-------------------------------------------------------
! Ed Brothers. 11/26/01
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP
!-------------------------------------------------------
! This subroutine is to get molecule information 
! and assign basis function.
!
    subroutine getMol(nAtomSave)
    use allmod
    implicit double precision(a-h,o-z)

    include 'mpif.h'
    
    logical :: present
    character(len=80) :: keyWD
    character(len=20) :: tempString
!----------------------------------------
! Xiao HE Superposition of Atomic density
!
!   1 1   H  1S          0.13943
!   2        2S          0.26785   0.51455
!----------------------------------------
    double precision H321(2,2),C321(9,9),N321(9,9),O321(9,9),tmpxiao(20,20),cl321(13,13)
! Constants
!    pi=3.1415926535897932385d0
    pi32 = pi**(1.5)
    
    
!-----------MPI/MASTER------------------------    
    masterwork:if (master) then
!-----------END MPI/MASTER------------------------

!----------------
! NORMAL QUICK
!----------------

    iattype=0
!   natom = 0
    iiatom = 0
    nelec = 0

    natom=natomsave   

! Read info from AMBER
    if (amber_interface_logic) then
        call read_AMBER_crd
    else

    open(infile,file=inFileName,status='old')
    call PrtAct(iOutfile,"Begin to Read Mol Info")

    ! The first line is Keyword
    ! the second is bland line
    read (infile,'(A80)') keywd
    read (infile,'(A80)') keywd

    ! Now it is time to read in the coordinates, and assign basis functions.
    ! Note this is blank line terminated.
    istart = 1
    ifinal = 80
    read (infile,'(A80)') keywd

    call rdword(keywd,ISTART,ifINAL)

    do WHILE ((istart /=  0) .and. (ifinal /= 0))
!++++++++++++++++++++++++++++++++++++++++++
! First, find the atom type.
!------------------------------------------
        iiatom = iiatom + 1 
        iattype(iiatom) = -1
        I=0
        do WHILE (iattype(iiatom) .eq. -1)
          if (keywd(ISTART:ifINAL) == symbol(I)) iattype(iiatom)=I
            I=I+1
        enddo
        chg(iiatom)=iattype(iiatom)
!
!++++++++++++++++++++++++++++++++++++++++++
! Next, find the xyz coordinates of the atom and convert to bohr.
!------------------------------------------

        do I=1,3
            istart=ifinal+1
            ifinal=80
            call rdword(keywd,ISTART,ifINAL)
            call rdnum(keywd,istart,temp,ierror)
            xyz(I,iiatom) = temp/bohr
        enddo

    ! Keep a running tally of electrons. This is a total not counting the charge.
    ! We'll continue this after everything else is read in.

!        nelec = nelec+iattype(natom)
        nelec = nelec+iattype(iiatom)

        read (infile,'(A80)') keywd
        call upcase(keywd,80)
        istart=1
        ifinal = 80
        call rdword(keywd,ISTART,ifINAL)
    enddo
    
    close(inFile)
    endif
! Finish reading molecular cooridate
!++++++++++++++++++++++++++++++++++++++++++


!++++++++++++++++++++++++++++++++++++++++++
! Check for errors.
!------------------------------------------
    if (natom == 1 .AND. quick_method%opt) then
        write (iOutFile,'(//" ONE ATOM = NO OPTIMIZATION!!!!"//)')
        quick_method%opt=.false.
    endif
!++++++++++++++++++++++++++++++++++++++++++
! Define the number of optimization cycles if not set. default is 9*natom
!------------------------------------------
    if (iopt == 0 .AND. quick_method%opt) then
        iopt=3*3*natom
        write (iOutFile,'("MAX OPTIMIZATION CYCLES = ", &
        & I5,"  (DEFAULT)")') iopt
    endif
!
!++++++++++++++++++++++++++++++++++++++++++
! At this point a blank line has been read in.  (That is what
! teminated the geometry read)  Now read in another line
! to see the specific grid has been requested,assuming this is a
! DFT job.
!-------------------------------------------
    if (quick_method%DFT .OR. quick_method%SEDFT) then
! Xiao HE 01/09/07 SG1 SG0 grid
       if (quick_method%ISG.eq.1) then
       itemp=50
       do I=1,itemp
        RGRID(I)=(I**2.d0)/dble((itemp+1-I)*(itemp+1-I))
        RWT(I)=2.d0*dble(itemp+1)*(dble(I)**5.d0) &
                           *dble(itemp+1-I)**(-7.d0)
       enddo
       else
       Continue
       endif

        if (natom > 1) then
            do Iatom=1,natom
                distnbor(Iatom)=1.D30
                do Jatom=1,Iatom-1
                    DIST=(xyz(1,Iatom)-xyz(1,Jatom))**2.d0
                    DIST=DIST+(xyz(2,Iatom)-xyz(2,Jatom))**2.d0
                    DIST=DIST+(xyz(3,Iatom)-xyz(3,Jatom))**2.d0
                    DIST=DIST**.5d0
                    distnbor(Iatom)=Min(distnbor(Iatom),DIST)
                enddo
                do Jatom=Iatom+1,natom
                    DIST=(xyz(1,Iatom)-xyz(1,Jatom))**2.d0
                    DIST=DIST+(xyz(2,Iatom)-xyz(2,Jatom))**2.d0
                    DIST=DIST+(xyz(3,Iatom)-xyz(3,Jatom))**2.d0
                    DIST=DIST**.5d0
                    distnbor(Iatom)=Min(distnbor(Iatom),DIST)
                enddo
            enddo
        endif
    endif
!
! That's all the DFT stuff.  Now back to molecule set up.
!
!!++++++++++++++++++++++++++++++++++++++++++
! If this is a core approximatation calculation, modify the atomic
! charges and the number of electrons.
!
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
            if (iattype(Iatm) >= 18) &
            write (iOutFile,'("ATOM OUT OF RANGE FOR CORE")')
        enddo
    endif
    
    write(iOutFile,'("TOTOAL ATOM NUMBER = ",i4)') natom
    j=0
    do i=1,natom
      if (iattype(I).eq.1) j=j+1
    enddo
    nHAtom=j
    nNonHAtom=natom-j
    write(iOutFile,'("NUMBER OF HYDROGEN ATOM = ",i4)') nhatom
    write(iOutFile,'("NUMBER OF NON-HYDROGEN ATOM = ",i4)') nNonHAtom
    write(iOutFile,'("DEGREE OF FREEdoM = ",i5)') natom*3-6

!++++++++++++++++++++++++++++++++++++++++++
! Modify the electron numbers based on charge and multiplicity.
!
    nelec = nelec - molchg
    xmulttest = mod(dble(nelec),2.d0)

    if (xmulttest /= 0 .AND. .NOT. quick_method%unrst) then
        write (iOutFile,'("WARNING: UNPAIRED ELECTRONS REQUIRE", &
        & " UNRESTRICTED CALCULATIONS.")')
       quick_method%UNRST=.true.
    endif

    if (imult /= 1 .AND. .NOT. quick_method%unrst) &
    write (iOutFile,'("WARNING: HIGHER MULTIPLICITIES REQUIRE", &
    & " UNRESTRICTED CALCULATIONS.")')

    if (quick_method%unrst) then
        nelecb = nelec
        nelec = 0

        do WHILE (nelec.lt.nelecb)
            nelecb = nelecb-1
            nelec = nelec +1
        enddo

        if (imult == 2 .AND. nelec-1 /= nelecb) &
        write (iOutFile,'("WARNING: INCORRECT NUMBER OF ELECTRONS", &
        & " FOR A doUBLET.")')

        if (imult == 3) then
            nelec = nelec+1
            nelecb = nelecb - 1
            if (nelec-2 /= nelecb) &
            write (iOutFile,'("WARNING: INCORRECT NUMBER OF ELECTRONS", &
            & " FOR A TRIPLET.")')

        endif

        if (imult == 4) then
            nelec = nelec+1
            nelecb = nelecb - 1
            if (nelec-3 /= nelecb) &
            write (iOutFile,'("WARNING: INCORRECT NUMBER OF ELECTRONS", &
            & " FOR A QUADRUPLET.")')

        endif
        write (iOutFile,'("NUMBER OF ALPHA ELECTRONS = ",I4)')nelec
        write (iOutFile,'("NUMBER OF BETA ELECTRONS  = ",I4)')nelecb
    else
        write (iOutFile,'(/"NUMBER OF ELECTRONS = ",I4)')nelec
    endif
    
!-----------MPI/MASTER------------------------
    endif masterwork
!-----------END MPI/MASTER------------------------  


!-----------MPI/ALL NODES------------------------
    if (bMPI) then
      call mpi_setup_mol1()
      call MPI_BARRIER(MPI_COMM_WORLD,mpierror)
    endif
    
!-----------END MPI/ALL NODES------------------------


!++++++++++++++++++++++++++++++++++++++++++
! Alessandro GENONI 03/05/2007
! Only for ECP calculations:
! * Allocate arrays whose dimensions depend on NATOM (allocateatoms_ecp)
! * Read the Effective Core Potentials (ECPs), modify the atomic charges 
!   and the total number of electrons (readecp)
!
    
    if (quick_method%ecp) then
      call allocateatoms_ecp
      call readecp
    END if
!++++++++++++++++++++++++++++++++++++++++++
! At this point we have the positions and identities of the atoms. We also
! have the number of electrons. Now we must assign basis functions. This
! is done in a subroutine.
!++++++++++++++++++++++++++++++++++++++++++

    call readbasis(natom,0,0,0,0)
    
    allocate(Apri(jbasis,jbasis))
    allocate(Kpri(jbasis,jbasis))
    allocate(cutprim(jbasis,jbasis))
    allocate(Ppri(3,jbasis,jbasis))
    allocate(Xcoeff(jbasis,jbasis,0:3,0:3))
    If(quick_method%DFT)then
       allocate(phiXiao(nbasis))
       allocate(dPhidXXiao(nbasis))
       allocate(dPhidYXiao(nbasis))
       allocate(dPhidZXiao(nbasis))
    endif

!-----------MPI/MASTER------------------------
    if (master) then
!-----------END MPI/MASTER------------------------


!++++++++++++++++++++++++++++++++++++++++++
    if (quick_method%debug) then
        do I=1,nbasis
            write(iOutFile,'(/"BASIS FUNCTON ",I4," ON ATOM ",I4)') &
            I,ncenter(I)
            write(iOutFile,'("THIS IS AN ",I1,I1,I1," FUNCTION")') &
            itype(1,I),itype(2,I),itype(3,I)
            write(iOutFile,'("THERE ARE ",I4," CONTRACTED GAUSSIANS")') &
            ncontract(I)
            do J=1,ncontract(I)
                write(iOutFile,'(F10.6,6x,F10.6)') aexp(J,I),dcoeff(J,I)
            enddo
        enddo
    endif
!++++++++++++++++++++++++++++++++++++++++++
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

    write (iOutFile,'("BASIS FUNCTIONS = ",I4)') nbasis
!    stop


!
! Alessandro GENONI 03/21/2007
! Store the normalized primitives coefficients for ECP calculations
!
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
!++++++++++++++++++++++++++++++++++++++++++
!    do I=1,nbasis
!        gauss(i)%dcoeff = dcoeff(:,i)
!    enddo

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

!++++++++++++++++++++++++++++++++++++++++++
! Which is our function to work with.
    if (quick_method%DFT .OR. quick_method%SEDFT) then
        write (iOutFile, &
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
            write (iOutFile,'(I4,7x,F12.6)') Ibas,radial
        enddo
    endif
!+++++++++++++++++++++++++++++++++++++++++
! Read params for semi-emipeircal DFT
    if (quick_method%SEDFT) then
        open(17,file='PARAMS',status='old')
        read (17,'(F12.7)') At2prm(0,0,0,1)
        read (17,'(F12.7)') At2prm(0,0,0,6)
        read (17,'(F12.7)') At2prm(1,0,0,6)
        At2prm(0,1,0,6) = At2prm(1,0,0,6)
        At2prm(0,0,1,6) = At2prm(1,0,0,6)
        read (17,'(F12.7)') Bndprm(0,0,0,1)
        read (17,'(F12.7)') Bndprm(0,0,0,6)
        read (17,'(F12.7)') Bndprm(1,0,0,6)
        Bndprm(0,1,0,6) = Bndprm(1,0,0,6)
        Bndprm(0,0,1,6) = Bndprm(1,0,0,6)
        close (17)
        write (iOutFile,'(/"HYDROGEN PARAMS")')
        write (iOutFile,'(F12.7,F12.7)') At2prm(0,0,0,1), &
        Bndprm(0,0,0,1)
        write (iOutFile,'(/"CARBON PARAMS")')
        write (iOutFile,'(F12.7,F12.7,F12.7,F12.7)') At2prm(0,0,0,6), &
        At2prm(1,0,0,6),Bndprm(0,0,0,6),Bndprm(1,0,0,6)
     endif
   
!++++++++++++++++++++++++++++++++++++++++++
!  initialize density matrix

    call initialGuess
    if (quick_method%debug) then
        write(iOutFile,*) "DENSITY MATRIX AFTER INITIAL GUESS"
        call PriSym(iOutFile,nbasis,dense,'f14.8')
    endif
!++++++++++++++++++++++++++++++++++++++++++
! Then read external charges
call flush(6)
    if(quick_method%extCharges) then
        if (amber_interface_logic) then
            call read_Amber_Charge
        else
            call readExtCharges   
        endif
    endif
call flush(6)
!++++++++++++++++++++++++++++++++++++++++++
! Now it's time to output mol infor
    Write (iOutFile,'(//,"INPUT GEOMETRY:")')
    do I=1,natom
        Write (iOutFile,'(A2,6x,F10.4,3x,F10.4,3x,F10.4)') &
        symbol(iattype(I)),xyz(1,I)*bohr, &
        xyz(2,I)*bohr,xyz(3,I)*bohr
    enddo
    if(quick_method%extcharges)then
        write(iOutFile,'(/"EXTERNAL POINT CHARGES: (Q,X,Y,Z)")')
        do i=1,nextatom
            write(iOutFile,'(F7.4,3(F10.4,1x))') extchg(i),extxyz(:,i)*bohr
        enddo
    endif
!+++++++++++++++++++++++++++++++++++++++++++
! Output Distance Matrix    
    do i=1,natom
      do j=1,natom
        atomdistance(i,j)=0d0
        do k=1,3
          atomdistance(i,j)=atomdistance(i,j)+(xyz(k,i)-xyz(k,j))**2
        enddo
        atomdistance(i,j)=sqrt(atomdistance(i,j))*bohr
      enddo
    enddo
    
    ! if no. of atom is less than 30, then output them
    if (natom.le.30) then
        write(iOutFile,*)
        write(iOutFile,'("DISTANCE MATRIX:")')
        call PriSym(iOutFile,natom,atomdistance,'f10.4')
    endif
    call PrtAct(iOutfile,"End Reading Mol Info ")
!++++++++++++++++++++++++++++++++++++++++++

!-----------MPI/MASTER------------------------
    endif
!-----------END MPI/MASTER------------------------  
    if (bMPI) then
      call mpi_setup_mol2()
      call MPI_BARRIER(MPI_COMM_WORLD,mpierror)
    endif
    
    
    return
    end subroutine getmol


!********************************************
! Read External Charges
!--------------------------------------------
!JF 9/2009
!
!Reads external MM Point charges to be included in the one electron part of
!the Hamiltonian.
!
!Invoke with EXTCHARGES keyword.
!Put the charges in the input file after QM atoms and a blank line in this format:
!CHARGE X Y Z
!
!First I will run through the inputfile and count the charges, then I'll allocate
!the arrays and then fill them. 
!

subroutine readextcharges
    use allmod
    implicit double precision(a-h,o-z)

    character(len=80) :: keywd
     
    open(infile,file=infilename,status='old')
!Need to skip keyword lines
    istart = 1
    ifinal = 80   
    do 
        read(infile,'(A80)') keywd  
        call rdword(keywd,istart,ifinal)
        if(istart==0 .or. ifinal==0) exit
    enddo
!Need to skip QM atom list
    istart = 1
    ifinal = 80
    do 
        read(infile,'(A80)') keywd  
        call rdword(keywd,istart,ifinal)
        if(istart==0 .or. ifinal==0) exit
    enddo
!Ready to count point charges
    istart = 1
    ifinal = 80
    nextatom=0
    do 
        read(infile,'(A80)') keywd
        call rdword(keywd,istart,ifinal)
        if(istart==0 .or. ifinal==0) exit
        nextatom=nextatom+1 
    enddo
    close(infile)

    allocate(extchg(nextatom),extxyz(3,nextatom))

    open(infile,file=infilename,status='old')
!Need to skip keyword lines
    istart = 1
    ifinal = 80   
    do 
        read(infile,'(A80)') keywd  
        call rdword(keywd,istart,ifinal)
        if(istart==0 .or. ifinal==0) exit
    enddo
!Need to skip QM atom list
    istart = 1
    ifinal = 80
    do 
        read(infile,'(A80)') keywd  
        call rdword(keywd,istart,ifinal)
        if(istart==0 .or. ifinal==0) exit
    enddo
!Ready to read point charges into arrays
    iextatom=1
    do 
        istart = 1
        ifinal = 80
        read(infile,'(A80)') keywd
        call rdword(keywd,istart,ifinal)
        if(istart==0 .or. ifinal==0) exit
        call rdnum(keywd,istart,extchg(iextatom),ierr)
        if(ierr==1)then
            print *, 'error reading extchg ',iextatom
            stop
        endif
        do i=1,3
            istart=ifinal+1
            ifinal=80
            call rdword(keywd,istart,ifinal)
            call rdnum(keywd,istart,temp,ierr)
            if(ierr==1)then
                print *, 'error reading extxyz',i,iextatom 
                stop
            endif
            extxyz(i,iextatom) = temp/bohr
        enddo 
        iextatom=iextatom+1
    enddo
    close(infile)
    
end subroutine readextcharges


!********************************************
! Initial Densitry Matrix
!--------------------------------------------
subroutine initialGuess
    use allmod
    implicit double precision(a-h,o-z)
    logical :: present
    character(len=80) :: keyWD


!++++++++++++++++++++++++++++++++++++++++++
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
            write (iOutFile, &
            & '(/,"CONVERTING RESTRICTED DENSITY TO UNRESTRICTED")')
            do I=1,nbasis
                do J =1,nbasis
                    DENSE(J,I) = DENSE(J,I)/2.d0
                    DENSEB(J,I) = DENSE(J,I)
                enddo
            enddo
        endif
    else
!++++++++++++++++++++++++++++++++++++++++++
! Xiao HE Initial Guess
!      If(ISAD.eq.1)then
!        diagelement=dble(nelec)/dble(nbasis)
!        diagelementb=dble(nelecb)/dble(nbasis)+1.d-8
!        do I=1,nbasis
!            DENSE(I,I)=diagelement
!            DENSEB(I,I)=diagelementb
!        enddo
!
!      else

       
        do ixiao=1,npmfcc
          print*,mfccbases(ixiao),mfccbasef(ixiao),mfccbasescap(ixiao),mfccbasefcap(ixiao)
          print*,matombases(ixiao),matombasescap(ixiao)
        enddo        

        do ixiao=1,kxiaoconnect
          print*,'check connection',ixiao
          print*,mfccbasescon(ixiao),mfccbasefcon(ixiao),matombasescon(ixiao)
          print*,mfccbasescon2(ixiao),mfccbasefcon2(ixiao),matombasescon2(ixiao)
        enddo
          

!         print*,'IMFCC=',IMFCC
!++++++++++++++++++++++++++++++++++++++++++
! MFCC Initial Guess 
         if(quick_method%MFCC)then

           do i=1,nbasis
             do j=1,nbasis
               DENSE(i,j)=0.0d0
             enddo
           enddo
        
           do ixiao=1,npmfcc
            do i=mfccbases(ixiao),mfccbasef(ixiao)
             do j=mfccbases(ixiao),mfccbasef(ixiao)
               DENSE(matombases(ixiao)+i-mfccbases(ixiao),matombases(ixiao)+j-mfccbases(ixiao)) &
                    =DENSE(matombases(ixiao)+i-mfccbases(ixiao),matombases(ixiao)+j-mfccbases(ixiao))+ &
                     mfccdens(ixiao,i-mfccbases(ixiao)+1,j-mfccbases(ixiao)+1)
            if(mfccdens(ixiao,i-mfccbases(ixiao)+1,j-mfccbases(ixiao)+1).gt.0.3d0)then
               print*,'fragment',ixiao,matombases(ixiao)+i-mfccbases(ixiao), &
                 matombases(ixiao)+j-mfccbases(ixiao),mfccdens(ixiao,i-mfccbases(ixiao)+1, &
                                    j-mfccbases(ixiao)+1)
            endif
             enddo
            enddo
           enddo

           do ixiao=1,npmfcc-1
            do i=mfccbasescap(ixiao),mfccbasefcap(ixiao)
             do j=mfccbasescap(ixiao),mfccbasefcap(ixiao)
               DENSE(matombasescap(ixiao)+i-mfccbasescap(ixiao),matombasescap(ixiao)+j-mfccbasescap(ixiao))= &
               DENSE(matombasescap(ixiao)+i-mfccbasescap(ixiao),matombasescap(ixiao)+j-mfccbasescap(ixiao)) &
                     -mfccdenscap(ixiao,i-mfccbasescap(ixiao)+1,j-mfccbasescap(ixiao)+1)
             if(mfccdenscap(ixiao,i-mfccbasescap(ixiao)+1,j-mfccbasescap(ixiao)+1).gt.0.3d0)then
               print*,'cap',ixiao,matombasescap(ixiao)+i-mfccbasescap(ixiao), &
                    matombasescap(ixiao)+j-mfccbasescap(ixiao),mfccdenscap(ixiao,i-mfccbasescap(ixiao)+1, &
                          j-mfccbasescap(ixiao)+1)
             endif
             enddo
            enddo
           enddo

!    do i=1,nbasis
!      do j=1,nbasis
!!        if(DENSE(i,j).gt.0.001d0)then
!          DENSE(i,j)=0.0d0
!!        endif
!      enddo
!    enddo

           do ixiao=1,kxiaoconnect
            do i=mfccbasesconi(ixiao),mfccbasefconi(ixiao)
             do j=mfccbasesconi(ixiao),mfccbasefconi(ixiao)
          DENSE(matombasesconi(ixiao)+i-mfccbasesconi(ixiao),matombasesconi(ixiao)+j-mfccbasesconi(ixiao))= &
          DENSE(matombasesconi(ixiao)+i-mfccbasesconi(ixiao),matombasesconi(ixiao)+j-mfccbasesconi(ixiao)) &
                     -mfccdensconi(ixiao,i-mfccbasesconi(ixiao)+1,j-mfccbasesconi(ixiao)+1)
             if(mfccdensconi(ixiao,i-mfccbasesconi(ixiao)+1,j-mfccbasesconi(ixiao)+1).gt.0.3d0)then
               print*,'connect-I',ixiao,matombasesconi(ixiao)+i-mfccbasesconi(ixiao), &
               matombasesconi(ixiao)+j-mfccbasesconi(ixiao),mfccdensconi(ixiao,i-mfccbasesconi(ixiao)+1, &
                          j-mfccbasesconi(ixiao)+1)
             endif
             enddo
            enddo
           enddo

           do ixiao=1,kxiaoconnect
            do i=mfccbasesconj(ixiao),mfccbasefconj(ixiao)
             do j=mfccbasesconj(ixiao),mfccbasefconj(ixiao)
          DENSE(matombasesconj(ixiao)+i-mfccbasesconj(ixiao),matombasesconj(ixiao)+j-mfccbasesconj(ixiao))= &
          DENSE(matombasesconj(ixiao)+i-mfccbasesconj(ixiao),matombasesconj(ixiao)+j-mfccbasesconj(ixiao)) &
                     -mfccdensconj(ixiao,i-mfccbasesconj(ixiao)+1,j-mfccbasesconj(ixiao)+1)
             if(mfccdensconj(ixiao,i-mfccbasesconj(ixiao)+1,j-mfccbasesconj(ixiao)+1).gt.0.3d0)then
               print*,'connect-J',ixiao,matombasesconj(ixiao)+i-mfccbasesconj(ixiao), &
               matombasesconj(ixiao)+j-mfccbasesconj(ixiao),mfccdensconj(ixiao,i-mfccbasesconj(ixiao)+1, &
                          j-mfccbasesconj(ixiao)+1)
             endif
             enddo
            enddo
           enddo

           do ixiao=1,kxiaoconnect
            do i=mfccbasesconi(ixiao),mfccbasefconi(ixiao)
             do j=mfccbasesconi(ixiao),mfccbasefconi(ixiao)
          DENSE(matombasesconi(ixiao)+i-mfccbasesconi(ixiao),matombasesconi(ixiao)+j-mfccbasesconi(ixiao))= &
          DENSE(matombasesconi(ixiao)+i-mfccbasesconi(ixiao),matombasesconi(ixiao)+j-mfccbasesconi(ixiao)) &
                     +mfccdenscon(ixiao,i-mfccbasesconi(ixiao)+1,j-mfccbasesconi(ixiao)+1)
             if(mfccdenscon(ixiao,i-mfccbasesconi(ixiao)+1,j-mfccbasesconi(ixiao)+1).gt.0.3d0)then
               print*,'connect-IJ',ixiao,matombasesconi(ixiao)+i-mfccbasesconi(ixiao), &
               matombasesconi(ixiao)+j-mfccbasesconi(ixiao),mfccdenscon(ixiao,i-mfccbasesconi(ixiao)+1, &
                          j-mfccbasesconi(ixiao)+1)
             endif
             enddo
            enddo
           enddo


           do ixiao=1,kxiaoconnect
            do i=mfccbasesconj(ixiao),mfccbasefconj(ixiao)
             do j=mfccbasesconj(ixiao),mfccbasefconj(ixiao)
                
          iixiaotemp=mfccbasefconi(ixiao)-mfccbasesconi(ixiao)+1

          DENSE(matombasesconj(ixiao)+i-mfccbasesconj(ixiao),matombasesconj(ixiao)+j-mfccbasesconj(ixiao))= &
          DENSE(matombasesconj(ixiao)+i-mfccbasesconj(ixiao),matombasesconj(ixiao)+j-mfccbasesconj(ixiao)) &
                     +mfccdenscon(ixiao,iixiaotemp+i-mfccbasesconj(ixiao)+1, &
                      iixiaotemp+j-mfccbasesconj(ixiao)+1)
             if(mfccdenscon(ixiao,iixiaotemp+i-mfccbasesconj(ixiao)+1, &
                iixiaotemp+j-mfccbasesconj(ixiao)+1).gt.0.3d0)then
               print*,'connect-IJ',ixiao,matombasesconj(ixiao)+i-mfccbasesconj(ixiao), &
!                     iixiaotemp+i-mfccbasesconj(ixiao)+1,iixiaotemp+j-mfccbasesconj(ixiao)+1, &
               matombasesconj(ixiao)+j-mfccbasesconj(ixiao),mfccdenscon(ixiao,iixiaotemp+ &
                          i-mfccbasesconj(ixiao)+1, &
                          iixiaotemp+j-mfccbasesconj(ixiao)+1)
             endif
             enddo
            enddo
           enddo


         else
!++++++++++++++++++++++++++++++++++++++++++
!  SAD inital guess
         nincrease=0

         do Iatm=1,natom
            do ixiaosad=1,10
             If(symbol(iattype(Iatm)).eq.atomxiao(ixiaosad))then
!                 print*,'hexiao',iatm,ixiaosad,atomxiao(ixiaosad),atombasis(ixiaosad)
                 do i=1,atombasis(ixiaosad)
                   do j=1,atombasis(ixiaosad)
                     DENSE(i+nincrease,j+nincrease)=atomdens(ixiaosad,i,j)
!                     print*,i+nincrease,j+nincrease,DENSE(i+nincrease,j+nincrease)
                   enddo
                 enddo
                 nincrease=nincrease+atombasis(ixiaosad)
             endif
            enddo
               
        enddo
       endif
    endif
end subroutine initialGuess
