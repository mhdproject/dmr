!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"
!#define DEBUG_AMR
!#define DEBUG
      subroutine amr_reorder_grid


      use paramesh_dimensions
      use physicaldata
      use tree


      use paramesh_interfaces, only :
     &                                amr_morton_order

      implicit none


#ifdef MPI_USED
      include 'mpif.h'
#endif


#ifdef TIMINGS
#include "timer.fh"
      integer :: itimer1,itimer2
#endif
      integer :: nprocs,mype
      integer :: lnblocks_old
      integer :: ierr, i
      integer :: lnblocks2, tot_blocks, tot_blocksa
      integer :: max_blocks, min_blocks
      logical :: l_move_solution

#ifdef SHMEM_USED
      integer,external :: shmem_n_pes,shmem_my_pe
#endif


!pmn temporary fix to be tested for Rick
!    same fix is necessary in amr_refine_derefine
      lnblocks_old = lnblocks
!pmn end temporary fix to be tested for Rick



#ifdef SHMEM_USED
      nprocs = shmem_n_pes()
      mype   = shmem_my_pe()

      work_block(:) = 0.
      do i = 1,lnblocks
         if (nodetype(i).eq.1) work_block(i) = 2.        !<<< USER EDIT
         if (nodetype(i).ge.2) work_block(i) = 1.        !<<< USER EDIT
      end do

      l_move_solution = .false.
      call amr_morton_order (nprocs,lnblocks_old,
     .                       l_move_solution)

#endif
#ifdef MPI_USED

      call MPI_COMM_SIZE (MPI_COMM_WORLD,nprocs,ierr)
      call MPI_COMM_RANK (MPI_COMM_WORLD,mype,ierr)



      newchild(:) = .FALSE.


      call MPI_ALLREDUCE (lnblocks,tot_blocks,1,MPI_INTEGER,
     &                    MPI_SUM,MPI_COMM_WORLD,ierr)

#ifdef DEBUG_AMR
      write(*,*) 'refderef: tot_blocks ',tot_blocks,mype
      if (mype.eq.0) then
         print *,' tot_blocks before ',tot_blocks
      end if

! I copy lnblocks to lnblocks2 since lnblocks2 can be put in a save statement.
      lnblocks2 = lnblocks 
      call MPI_ALLREDUCE (lnblocks2,max_blocks,1,MPI_INTEGER,
     &                    MPI_MAX,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE (lnblocks2,min_blocks,1,MPI_INTEGER,
     &                    MPI_MIN,MPI_COMM_WORLD,ierr)

      if (mype.eq.0) then
         print *, ' max_blocks 1',max_blocks
         print *, ' min_blocks 1',min_blocks
      end if
#endif

! set work values

      if (mype.eq.0) then
         print *, ' starting MORTON ORDERING '
      end if
      
      work_block(:) = 0.
      do i = 1,lnblocks
         if (nodetype(i).eq.1) work_block(i) = 2.        !<<< USER EDIT
         if (nodetype(i).ge.2) work_block(i) = 1.        !<<< USER EDIT
      end do

      l_move_solution = .false.
      call amr_morton_order (lnblocks_old,nprocs,mype,
     .                       l_move_solution)


      
! I copy lnblocks to lnblocks2 since lnblocks2 can be put in a save statement.
      lnblocks2 = lnblocks 
      call MPI_ALLREDUCE (lnblocks2,tot_blocksa,1,MPI_INTEGER,
     &                    MPI_SUM,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE (lnblocks2,max_blocks,1,MPI_INTEGER,
     &                    MPI_MAX,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE (lnblocks2,min_blocks,1,MPI_INTEGER,
     &                    MPI_MIN,MPI_COMM_WORLD,ierr)


      if (mype.eq.0) then
         print *, ' tot_blocks after ',tot_blocksa
         print *, ' max_blocks 2',max_blocks
         print *, ' min_blocks 2',min_blocks
      end if


      call amr_morton_process()

#endif

!
! set grid modification flag
      grid_changed = 1


      return
      end subroutine amr_reorder_grid
