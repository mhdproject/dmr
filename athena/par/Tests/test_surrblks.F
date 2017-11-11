!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

      program test_surrblks






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! include file to define physical qualities of the model and mesh
      use paramesh_dimensions
      use physicaldata

! include file defining the tree
      use tree
      use workspace

! include file required for shmem library.
      include 'mpp/shmem.fh'
      include 'shmem_reduce.fh'

#include "test_defs.fh"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local amr variables
      integer nprocs,mype
      integer shmem_my_pe,shmem_n_pes,num_procs

      save mype

      integer surrblks(2,3,3,3)
!
! application specific variables

      integer iopt,nlayers
      integer ierror_sum,ierror_tot
      logical lrefine_again
      logical rflags(maxblocks_tr)

      save ierror_sum,ierror_tot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef SGI_SHMEM
      call comm_start(maxprocs,nprocs,mype)
#endif
      call amr_initialize


      mype = shmem_my_pe()
      nprocs = shmem_n_pes()

      print *,' nprocs = ',nprocs,mype


! set default value of dz and z0 to cater for 2D case.
      z0 = 0.
      dz = 0.

      rflags(:) = .true.

      ierror_sum = 0
      ierror_tot = 0

      iopt = 1
      nlayers = nguard
      if(mype.eq.0) write(*,*) 'nlayers = ',nlayers

!

! set a limit on the refinement level
      lrefine_max = 5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! set up initial grid state.

! set up step function on a single block covering the whole cubic domain
      lnblocks = 0
      if(mype.eq.0.) then
                lnblocks = 1
                bsize(:,1)=1.
                coord(:,1) = .5
                bnd_box(1,:,1) = .0
                bnd_box(2,:,1) = 1.0
                nodetype(1) = 1
                lrefine(1) = 1

                neigh(:,:,1) = -21

                refine(1)=.true.
      endif

      call shmem_barrier_all()


      loop_count=0
! Now cycle over blocks adjusting refinement of initial setup as required
      do while(loop_count.lt.2)

! refine grid and apply morton reordering to grid blocks if necessary
      call amr_refine_derefine

      refine(:) = .true.

      loop_count=loop_count+1

      enddo
      call shmem_barrier_all()


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!start test

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      do ii=0,nprocs-1
      if(mype.eq.ii) then
      do l=1,lnblocks
         write(*,*) 'proc ',ii,' block ',l,' coord= ',
     .                  (coord(icoord,l),icoord=1,ndim),
     .       ' size = ',bsize(1,l)
         write(*,*) 'proc ',ii,' block ',l,' parent= ',
     .       parent(1,l),parent(2,l)
         write(*,*) 'proc ',ii,' block ',l,' nodety= ',
     .       nodetype(l)
      enddo
      endif
      call shmem_barrier_all()
      enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call shmem_barrier_all()

      if(mype.eq.0) then
      neigh(:,1,35) = -1
      neigh(:,2,26) = -1
      neigh(:,3,36) = -1
      neigh(:,4,18) = -1
      neigh(:,5,66) = -1
      neigh(:,6,30) = -1
      endif

      call shmem_barrier_all()

      if(mype.eq.0) then
      lb = 69
      call amr_surrounding_blks(mype,surrblks,lb)
      endif

      call amr_close()

         
      end
