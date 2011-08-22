! Ed Brothers. 01/22/03
! 3456789012345678901234567890123456789012345678901234567890123456789012<<STOP

    subroutine frequency
    use allmod
    implicit double precision(a-h,o-z)

!    double precision :: kB
! dimension V(3,maxatm*3),HARMONIC(maxatm*3),EV(maxatm*3,maxatm*3)
    dimension V(3,natom*3),HARMONIC(natom*3),EV(natom*3,natom*3)

    call prtAct(ioutfile,"Begin Frequency calculation")
!    pi=3.1415926535897932385d0

! Note:  This can be changed.

    TempK = 298.15d0

! Conversion factor to get from atomic units:
! from atomic charge to ESU = 4.803242D-10
! from particle to mole = 6.0221367D23
! from Bohrs to cm = .529177249D-8
! And since the eigenvalues are 4 Pi^2 (freq)^2 you also need the speed
! of light in cm = 2.99792458D10

    convfact = (4.803242D-10)**2.d0
    convfact = convfact* 6.0221367D23
    convfact = convfact/ (4 *pi*pi)
    convfact = convfact/ (2.99792458D10)**2.d0
    convfact = convfact/ (.529177249D-8)**3.d0

! This procedure calculates the frequencies of the molecule given the
! hessian.

! First, mass weight the hessian.

    write (ioutfile,'(/" MASS WEIGHT THE HESSIAN. ")')
    DO I=1,natom
        DO J=1,natom
            denom = (emass(quick_molspec%iattype(I))*emass(quick_molspec%iattype(J)))**(-.5d0)
            ISTART = (I-1)*3
            JSTART = (J-1)*3
            DO K=1,3
                DO L=1,3
                    quick_qm_struct%hessian(JSTART+L,ISTART+K)=quick_qm_struct%hessian(JSTART+L,ISTART+K)*denom
                ENDDO
            ENDDO
        ENDDO
    ENDDO

    call PriHessian(ioutfile,3*natom,quick_qm_struct%hessian,'f12.6')

! Diagonalize the Hessian. Also write out the frequencies.
    non=0
    call hessDIAG(natom*3,quick_qm_struct%hessian,non,quick_method%DMCutoff,V,HARMONIC,quick_qm_struct%idegen,EV,IERROR)
    write (ioutfile,'(/" THE HARMONIC FREQUENCIES (1/cm): ")')
    write (ioutfile,*)
    DO I=1,natom*3
        HARMONIC(I) = convfact*HARMONIC(I)
        IF (HARMONIC(I) < 0.d0) THEN
            HARMONIC(I) = -1.d0* DABS(HARMONIC(I))**.5d0
        ELSE
            HARMONIC(I) = HARMONIC(I)**.5d0
        ENDIF
        write (ioutfile,'(6x,F15.5)') HARMONIC(I)
    ENDDO
    write (ioutfile,*)
    
! Now we have the frequencies.  Before continuing, we need to see if the
! molecule is linear.  Note that this assumes the molecule is polyatomic.

    IF (natom == 2) THEN
        ignore=5
    ELSE
        ignore=5
        DO I=3,natom
            CALL BNDANG(1,2,I,ANGLE)
            ANGLE = ANGLE*180.d0/pi
            IF (ANGLE > 5.d0 .AND. ANGLE < 175.d0) ignore=6
        ENDDO
    ENDIF
    IF (ignore == 5) write (ioutfile,'(/" MOLECULE IS LINEAR!! ")')

! Now we can calculate the zero point energy.
! ZPE = 1/2 (Sum over i) Harmonic(i)
! This is in wavenumbers, so:
! 1/cm * h * c
! 1/cm * J/sec * cm/sec = Joules/particle, and convert to Hartrees.
! This is from Radom and Scott, J. Chem. Phys. 1996, 100, 16502-13.

    write(ioutfile,'(2x,"------------------------")')
    write(ioutfile,'(5x,"Thermodynamic Data")')
    write(ioutfile,'(2x,"------------------------")')
    write (ioutfile,'("TEMPERATURE          = ",F15.7,"K")') tempK
    Ezp = 0.d0
    DO I=ignore+1,natom*3
        Ezp = Ezp + HARMONIC(I)
    ENDDO
    Ezp = Ezp*.5d0
    Ezp = Ezp*6.6260755D-34*2.99792458D10/4.3597482D-18
    write (ioutfile,'("ZERO POINT ENERGY    = ",F15.7)') Ezp

! We can also calculate the vibrational contribution to internal energy.

    Evib = 0.d0
    DO I=ignore+1,natom*3
        vibtemp = Harmonic(I)*(6.6260755D-34)/1.380658D-23*2.99792458D10
        ratio = vibtemp/TempK
        Evib = Evib +vibtemp*(.5d0 + 1.d0/(DEXP(ratio)-1.d0))
    ENDDO
    Evib = Evib*1.380658D-23/4.3597482D-18
    write (ioutfile,'("VIBRATIONAL ENERGY   = ",F15.7)') Evib

! Now we need the translation and rotation.

    Etrans=1.5d0*TempK*1.380658D-23/4.3597482D-18
    write (ioutfile,'("TRANSLATIONAL ENERGY = ",F15.7)') Etrans

    IF (ignore == 5) THEN
        Erot = TempK*1.380658D-23/4.3597482D-18
    ELSE
        Erot = 1.5d0*TempK*1.380658D-23/4.3597482D-18
    ENDIF
    write (ioutfile,'("ROTATIONAL ENERGY    = ",F15.7)') Etrans


    write (ioutfile,'("INTERNAL ENERGY      = ",F15.7)') &
    quick_qm_struct%Etot+Ezp+Etrans+Erot+Evib
    write (ioutfile,'("INTERNAL ENTHALPY    = ",F15.7)') &
    quick_qm_struct%Etot+Ezp+Etrans+Erot+Evib+TempK*1.380658D-23/4.3597482D-18

    call prtAct(ioutfile,"Finish Frequency calculation")
    end subroutine frequency

