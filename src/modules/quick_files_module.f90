!
!	quick_files_module.f90
!	new_quick
!
!	Created by Yipu Miao on 2/18/11.
!	Copyright 2011 University of Florida. All rights reserved.
!
!********************************************************
! File module.
!--------------------------------------------------------
!
    module quick_files_module
    implicit none

    character(len=80) :: inFileName,outFileName,dmxFileName,rstFileName, &
    CPHFFileName,basisDir,basisFileName,ECPDir,ECPFileName,basisCustName, &
    PDBFileName
    integer :: inFile = 15            ! input file
    integer :: iOutFile = 16          ! output file
    integer :: iDmxFile = 17          ! density matrix file
    integer :: iRstFile = 18          ! Restricted file
    integer :: iCPHFFile = 19         ! CPHF file
    integer :: iBasisFile = 20        ! basis set file
    integer :: iECPFile = 21          ! ECP file
    integer :: iBasisCustFile = 22
    integer :: iPDBFile=23            ! PDB file for D&C

    end module quick_files_module
!
!********************************************************

