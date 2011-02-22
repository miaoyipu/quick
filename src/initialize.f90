! --Yipu Miao 05/09/2010
!********************************************************
! This file list initial step of Quick
! Subroutine List:
! initialize1(ierror)
! outputCopyright(ierror)
! allocateAtoms
! allocateAtoms_ECP
! allocatebasis


    subroutine initialize1(ierr)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Initalize
!-----------------------------------------------------------------------
!     
!  initalize variables
!

    use allMod
    implicit none
    integer ierr

    !--------------------MPI/ALL NODES--------------------------------
    ! MPI Initializer
    call MPI_initialize()
    !------------------- End MPI  -----------------------------------
    
    call set_quick_method(quick_method)

    pi =  (4.d0*atan(1.0d0))
    pito3half = pi**1.5
    X0 = 2.0d0*(PI)**(2.5d0)
    X00 = 1.0d0
    integralCutoff=1.0d0/(10.0d0**6.0d0)    
    Primlimit=1.0d0/(10.0d0**6.0d0)
    gradcutoff=1.0d0/(10.0d0**7.0d0)
    
    molchg=100
    imult=0
    iscf=0
    pmaxrms=0.d0
    MAXDIISSCF=10
    NCYC=2
    signif=0.d0
    tol=0.d0

    allocate(atomDens(40,100,100))
    allocate(atomBasis(100))   
    allocate(MFCCDens(40,600,600))
    allocate(MFCCDensCap(40,400,400))
    allocate(MFCCDensCon(40,200,200))
    allocate(MFCCDensCon2(40,200,200))
    allocate(MFCCDensConI(40,200,200))
    allocate(MFCCDensConJ(40,200,200))
            
    ierr=0

    call cpu_time(timer_begin%TTotal) ! Trigger time counter    

    return

    end
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    subroutine outputCopyright(ierr)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ output_copyright
!-----------------------------------------------------------------------
!     
!  Output Copyright information
!

    use quick_files_module
    implicit none
    integer ierror,ierr
 
    write(iOutFile,*) " **************************************************************************"
    write(iOutFile,*) " **                      QUICK 0.01-2010-Jul                             **"
    write(iOutFile,*) " **                                                                      **"
    write(iOutFile,*) " **                        Copyright (c) 2010                            **"
    write(iOutFile,*) " **                Regents of the University of Florida                  **"
    write(iOutFile,*) " **                       All Rights Reserved.                           **"
    write(iOutFile,*) " **                                                                      **"
    write(iOutFile,*) " **  This software provided pursuant to a license agreement containing   **"
    write(iOutFile,*) " **  restrictions on its disclosure, duplication, and use. This software **"
    write(iOutFile,*) " **  contains confidential and proprietary information, and may not be   **"
    write(iOutFile,*) " **  extracted or distributed, in whole or in part, for any purpose      **"
    write(iOutFile,*) " **  whatsoever, without the express written permission of the authors.  **"
    write(iOutFile,*) " **  This notice, and the associated author list, must be attached to    **"
    write(iOutFile,*) " **  all copies, or extracts, of this software. Any additional           **"
    write(iOutFile,*) " **  restrictions set forth in the license agreement also apply to this  **"
    write(iOutFile,*) " **  software.                                                           **"
    write(iOutFile,*) " **************************************************************************"
    write(iOutFile,*)
    write(iOutFile,*) " Cite this work as:"
    write(iOutFile,*) " Miao, Y.: He, X.: Ayers,K: Brothers, E.: Merz,K. M. QUICK;"
    write(iOutFile,*) " University of Florida, Gainesville, FL, 2010"
    write(iOutFile,*)
    
    ierr=0

    return

    end
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


    subroutine allocateAtoms
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! allocateAtoms
! Ken Ayers 05/26/04
!-----------------------------------------------------------------
! This subroutine is to allocate matricies that will be used.
! Subroutines for allocation of the various matricies in quick
! These routines are not the ideal way to deal with allocation of
! variables.  Large sized arrays should only be allocated when
! they are needed.  Eventually someone will deal with this.
!-----------------------------------------------------------------
!
     use quick_molspec_module
     use quick_calculated_module
     use quick_basis_module
     implicit none

     allocate(xyz(3,natom))
     allocate(distnbor(natom))
     allocate(iattype(natom))
     allocate(chg(natom))
     allocate(gradient(3*natom))
     allocate(hessian(3*natom,3*natom))
     allocate(ifirst(natom))
     allocate(ilast(natom))
     allocate(ishellfirst(natom))
     allocate(ishelllast(natom))
     allocate(AtomDistance(natom,natom))
     
     allocate(Mulliken(natom))
     allocate(Lowdin(natom))

     end subroutine allocateAtoms
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    Subroutine allocateAtoms_ECP
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! allocateAtoms_ECP
! Alessandro GENONI 03/06/2007
!-----------------------------------------------------------------
! ECP Allocations
!----------------------------------------------------------------
!
     use quick_molspec_module
     use quick_ecp_module
     implicit none
!
     allocate(nelecp(natom))
     allocate(lmaxecp(natom))
!
    end subroutine allocateatoms_ecp
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    subroutine allocatebasis
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! allocateBasis
!-----------------------------------------------------------------
! Basis Set Allocations
!----------------------------------------------------------------
!    
    use quick_basis_module
    use quick_molspec_module
    use quick_gridpoints_module
    use quick_calculated_module
    use quick_scratch_module
    implicit none

    integer :: i

    allocate(aexp(maxcontract,nbasis))
    allocate(dcoeff(maxcontract,nbasis))
    allocate(gauss(nbasis))
    do i=1,nbasis
        allocate(gauss(i)%aexp(maxcontract))
        allocate(gauss(i)%dcoeff(maxcontract))
    enddo
    allocate(itype(3,nbasis))
    allocate(ncenter(nbasis))
    allocate(ncontract(nbasis))

    allocate(sigrad2(nbasis))

    allocate(Smatrix(nbasis,nbasis))
    allocate(X(nbasis,nbasis))
    allocate(O(nbasis,nbasis))
    allocate(CO(nbasis,nbasis))
    allocate(COB(nbasis,nbasis))
    allocate(VEC(nbasis,nbasis))
    allocate(DENSE(nbasis,nbasis))
    allocate(DENSEB(nbasis,nbasis))
    allocate(V2(3,nbasis))
    allocate(E(nbasis))
    allocate(EB(nbasis))
    allocate(idegen(nbasis))
! allocate(CPHFA(2*(maxbasis/2)**2,2*(maxbasis/2)**2))
! allocate(CPHFB(2*(maxbasis/2)**2,maxatm*3))
! allocate( B0(2*(maxbasis/2)**2))
! allocate(BU(2*(maxbasis/2)**2))

    allocate(hold(nbasis,nbasis))
    allocate(hold2(nbasis,nbasis))

    end subroutine allocatebasis
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++