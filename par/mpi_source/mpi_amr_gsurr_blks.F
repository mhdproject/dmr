!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"

!#define DEBUG

      subroutine mpi_amr_gsurr_blks(mype,nprocs)



!------------------------------------------------------------------------
!
! This routine calculates the addresses of surrounding blocks of 
! the block lb stored on processor pe. It is necessary that the
! routine ?? be called before this routine to make sure that
! sufficient tree data is available to compute the surrounding blocks.
!
!
! Written :     Peter MacNeice          June 2000
!------------------------------------------------------------------------
!
! Arguments:
!      mype           rank of local processor
!
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use timings
      use mpi_morton

      use paramesh_mpi_interfaces, only : mpi_amr_local_surr_blks

      implicit none

      include 'mpif.h'

      integer, intent(in)    ::  mype,nprocs

      integer ::  surrblks(3,3,3,3)
      integer ::  psurrblks(3,3,3,3)

!------------------------------------------------------------------------
! local variables

      integer :: lb,local_lb
      integer :: ierrorcode,ierr
      integer ::  max_no_of_blocks, ierror

#ifdef TIMING_MPI
      double precision :: time1
#endif /* TIMING_MPI */
!------------------------------------------------------------------------
!
!
! This routine assumes that the grid blocks are ordered by morton
! number and that any blocks with different refinement levels but
! the same morton number are ordered from coarse to fine.

!------------------------------------------------------------------------

#ifdef TIMING_MPI
      time1 = mpi_wtime()
#endif /* TIMING_MPI */


      if(.not.morton_limits_set) then
        write(*,*) 'Error : mpi_amr_local_surr_blks : morton info '
     &            ,'is out of date or was never set up. Make sure '
     &            ,'there is a call to mpi_morton_bnd before this '
     &            ,'routine is called. '
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif

      call MPI_ALLREDUCE(lnblocks,
     .                   max_no_of_blocks,
     .                   1,
     .                   MPI_INTEGER,
     .                   MPI_MAX,
     .                   MPI_COMM_WORLD,
     .                   ierror)

!----------------------------
      do lb = 1,lnblocks

      local_lb = lb

      if(local_lb.lt.1) then
        write(*,*) 'Error mpi_amr_surr_blks : block ',lb,mype,
     &             ' not found in buffer space on this proc '
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      else
        call mpi_amr_local_surr_blks(mype,local_lb,nprocs,
     &               max_no_of_blocks,
     &               surrblks,.false.,psurrblks)
        surr_blks(:,:,1:1+2*k2d,1:1+2*k3d,lb) =
     &               surrblks(:,:,2-k2d:2+k2d,2-k3d:2+k3d)
#ifdef DEBUG
        write(*,*) 'gsurr: pe ',mype,' block ',lb
        do k = 1,1+2*k3d
        do j = 1,1+2*k2d
        write(*,*) '(',lb,'/',mype,')  ',j,k,' : ',
     .          surr_blks(:,:,j,k,lb)
        enddo
        enddo
#endif /* DEBUG */
      endif

      enddo
!----------------------------

      gsurrblks_set = +1


#ifdef TIMING_MPI
              timer_amr_gsurr_blks =  timer_amr_gsurr_blks
     .                          + mpi_wtime() - time1
#endif /* TIMING_MPI */

      return
      end subroutine mpi_amr_gsurr_blks
