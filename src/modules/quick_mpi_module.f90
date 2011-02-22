!
!	quick_mpi_module.f90
!	new_quick
!
!	Created by Yipu Miao on 2/18/11.
!	Copyright 2011 University of Florida. All rights reserved.
!
    module quick_mpi_module
    
    integer :: mpierror
    integer :: mpirank
    integer :: myid
    integer :: namelen
    integer :: mpisize
    character(len=80) pname
    logical :: master = .true.  ! flag to show if the node is master node
    logical :: bMPI = .true.    ! flag to show if MPI is turn on

    integer, allocatable :: MPI_STATUS(:)
    integer, parameter :: MIN_1E_MPI_BASIS=200

    end module quick_mpi_module