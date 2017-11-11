!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

! mpi_morton module



      module mpi_morton


      private

! 
! variables for storing the morton environment
      public :: pe_remote

      integer, save,dimension(:),allocatable :: pe_remote

      end module mpi_morton
!-----------------------------------------------------------------------
