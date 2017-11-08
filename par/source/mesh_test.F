!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"


      subroutine mesh_test(mype)





!------------------------------------------------------------------------
!
! This routine tests the mesh to see if there is a jump of more
! than one refinement level anywhere in the grid. It does it
! by looping over all parents of leaf nodes, getting the list of blocks
! surrounding them, and checking that the appropriate
! parts of this list exist.
!
! Arguments:
!      mype             local processor
!
!
! Written :     Peter MacNeice          August 1998
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace


#ifdef MPI_USED
      use paramesh_mpi_interfaces, only : mpi_amr_local_surr_blks,
     &                                    mpi_morton_bnd
#endif /* MPI_USED */

      implicit none

#include "amr_shmem.fh"
      include 'mpif.h'

#ifdef TIMINGS
#include "timer.fh"
#endif

      integer, intent(in) :: mype

!-------------------------


! local arrays

      real :: xtest,ytest,ztest

      logical ldiag,l_parent
      integer shmem_n_pes
      integer psurround(3,3,3,3)
      integer psurrblks(3,3,3,3)
      integer :: tag_offset
      integer :: nprocs,lb,i,j,k,ierror,max_no_of_blocks

!-------------------------

      nprocs = shmem_n_pes()

#ifdef MPI_USED
      tag_offset = 100
      call mpi_morton_bnd(mype,nprocs,tag_offset)
      call MPI_ALLREDUCE(lnblocks,
     .                   max_no_of_blocks,
     .                   1,
     .                   MPI_INTEGER,
     .                   MPI_MAX,
     .                   MPI_COMM_WORLD,
     .                   ierror)
#endif /* MPI_USED */

      call shmem_barrier_all()

! Loop over blocks.
      if(lnblocks.gt.0) then
      do lb = 1,lnblocks
      if(nodetype(lb).eq.2) then


        ldiag = .true.
#ifdef MPI_USED
        l_parent = .false.
        call mpi_amr_local_surr_blks(mype,lb,nprocs,
     &                          max_no_of_blocks,
     &                          psurround,l_parent,psurrblks)
        psurround(:,:,2-k2d:2+k2d,2-k3d:2+k3d) =
     &       surr_blks(:,:,1:1+2*k2d,1:1+2*k3d,lb)
#else
        psurround(:,:,2-k2d:2+k2d,2-k3d:2+k3d) =
     &       surr_blks(:,:,1:1+2*k2d,1:1+2*k3d,lb)
#endif /* MPI_USED */

 
       if (ndim == 1) then
       do i = 1,3,2
        if(psurround(1,i,2,2).gt.-20.and.psurround(1,i,2,2).lt.0) then
         write(*,*) 'Bad grid : 2 level jump near block  (',
     .        mype,lb,') ','  surrounding block (',i,',2,2)'
        endif
       enddo
       end if

       if (ndim == 2) then
       do j = 2-k2d,2+k2d
       do i = 1,3
        if(i.ne.2.or.j.ne.2) then
        if(psurround(1,i,j,2).gt.-20.and.psurround(1,i,j,2).lt.0) then
         write(*,*) 'Bad grid : 2 level jump near block  (',
     .        mype,lb,') ','  surrounding block (',i,j,',2)'
        endif
        endif
       enddo
       enddo
       end if

       if (ndim == 3) then
       do k = 2-k3d,2+k3d
       do j = 2-k2d,2+k2d
       do i = 1,3
        if(i.ne.2.or.j.ne.2.or.k.ne.2) then
        if(psurround(1,i,j,k).gt.-20.and.psurround(1,i,j,k).lt.0) then
         write(*,*) 'Bad grid : 2 level jump near block  (',
     .        mype,lb,') ','  surrounding block (',i,j,k,')'
        endif
        endif
       enddo
       enddo
       enddo 
       end if


      endif
      enddo
      endif


      call shmem_barrier_all()

! Check consistency of bounding box and coordinate info.
      if(lnblocks.gt.0) then
      do lb = 1,lnblocks
        xtest = (coord(1,lb)-bnd_box(1,1,lb))*
     &          (coord(1,lb)-bnd_box(2,1,lb))
        if(xtest.ge.0.) write(*,*) 'coord and bnd_box are ',
     &'inconsistent in x direction for block ',lb,' proc ',mype
        if (ndim >= 2) then
        ytest = (coord(2,lb)-bnd_box(1,2,lb))*
     &          (coord(2,lb)-bnd_box(2,2,lb))
        if(ytest.ge.0.) write(*,*) 'coord and bnd_box are ',
     &'inconsistent in y direction for block ',lb,' proc ',mype
        end if
        if (ndim == 3) then
        ztest = (coord(3,lb)-bnd_box(1,3,lb))*
     &          (coord(3,lb)-bnd_box(2,3,lb))
        if(ztest.ge.0.) write(*,*) 'coord and bnd_box are ',
     &'inconsistent in z direction for block ',lb,' proc ',mype
        end if
      enddo
      endif

      call shmem_barrier_all()

      return
      end subroutine mesh_test
