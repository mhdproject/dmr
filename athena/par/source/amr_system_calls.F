!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"
!
! Some useful system calls which may or may not be available on some systems.
!
!
! $RCSfile: amr_system_calls.F,v $
! $Revision: 1.7 $
! $Date: 2003/07/17 12:06:22 $
!
!--------------------------------------------------------------
!
! abort
! This is an SGI IRIX or UNICOS command. use amr_abort instead.

!
! flush
! use amr_flush instead

!--------------------------------------------------------------




      subroutine amr_abort()
      implicit none
#ifdef MPI_USED
      include 'mpif.h'

      integer :: ierrorcode,ierr
#endif

#ifdef SHMEM_USED
#ifdef SGI_SHMEM
      call abort()
#endif
#ifdef SHMEM_DUMMYLIB
      stop
#endif
#endif

#ifdef MPI_USED
      call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
#endif

      return
      end subroutine amr_abort





      subroutine amr_flush(iunit)
      implicit none
      integer,intent(in) :: iunit

#ifdef SHMEM_USED
      call flush(iunit)
#endif

      return
      end subroutine amr_flush

