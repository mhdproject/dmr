!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"


      subroutine amr_1blk_copy_soln(level)



!------------------------------------------------------------------------
!
! This routine copies a global solution update from the time 
! synchronized global solution arrays, into the arrays used
! during the solution update, as is required when 
! using NO_PERMANENT_GUARDCELLS and the amr_1blk_guardcell routines.
!
! Arguments:
!      level        integer           if -1 then blocks at all refinement
!                                     levels are copied, otherwise only blocks
!                                     at level are copied.
!
! Written :     Peter MacNeice          May 1999
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use timings

#ifdef SHMEM_USED
      use paramesh_interfaces, only :
     .                       amr_gsurrounding_blks
#endif
#ifdef MPI_USED
      use paramesh_mpi_interfaces, only :
     .                       mpi_amr_global_domain_limits,
     .                       mpi_amr_morton_limits,
     .                       mpi_morton_bnd,
     .                       mpi_amr_gsurr_blks

#endif

      implicit none

      integer, intent(in) :: level

#include "amr_shmem.fh"


      integer,external :: shmem_n_pes,shmem_my_pe
      integer :: mype,nprocs,lb
      integer :: tag_offset, ivar

#ifdef TIMING_MPI
      include 'mpif.h'
      double precision :: time1
#endif /* TIMING_MPI */

!-------------------------

#ifdef TIMING_MPI
      time1 = mpi_wtime()
#endif /* TIMING_MPI */

#ifdef NO_PERMANENT_GUARDCELLS

        nprocs = shmem_n_pes()



        if(level.eq.-1) then

        if(nvar.gt.0) then
          do ivar=1,nvar
          if(int_gcell_on_cc(ivar)) then
          gt_unk(ivar,:,:,:,1:lnblocks) = unk(ivar,:,:,:,1:lnblocks)
          endif
          enddo
        endif

        if(nfacevar.gt.0) then
          do ivar=1,nfacevar
          if(int_gcell_on_fc(1,ivar)) then
          gt_facevarx(ivar,:,:,:,1:lnblocks) = 
     .       facevarx(ivar,:,:,:,1:lnblocks)
          endif
          enddo
!!!#if N_DIM >= 2
          if (ndim >= 2) then
          do ivar=1,nfacevar
          if(int_gcell_on_fc(2,ivar)) then
          gt_facevary(ivar,:,:,:,1:lnblocks) = 
     .       facevary(ivar,:,:,:,1:lnblocks)
          endif
          enddo
!!!#endif
          end if
!!!#if N_DIM == 3
          if (ndim == 3) then
          do ivar=1,nfacevar
          if(int_gcell_on_fc(3,ivar)) then
          gt_facevarz(ivar,:,:,:,1:lnblocks) = 
     .       facevarz(ivar,:,:,:,1:lnblocks)
          endif
          enddo
!!!#endif
          end if
        endif

        if(nvaredge.gt.0) then
!!!#if N_DIM > 1
          if (ndim > 1) then
          do ivar=1,nvaredge
          if(int_gcell_on_ec(1,ivar)) then
          gt_unk_e_x(ivar,:,:,:,1:lnblocks) = 
     .       unk_e_x(ivar,:,:,:,1:lnblocks)
          endif
          enddo
          do ivar=1,nvaredge
          if(int_gcell_on_ec(2,ivar)) then
          gt_unk_e_y(ivar,:,:,:,1:lnblocks) = 
     .       unk_e_y(ivar,:,:,:,1:lnblocks)
          endif
          enddo
!!!#if N_DIM == 3
          if (ndim == 3) then
          do ivar=1,nvaredge
          if(int_gcell_on_ec(3,ivar)) then
          gt_unk_e_z(ivar,:,:,:,1:lnblocks) = 
     .       unk_e_z(ivar,:,:,:,1:lnblocks)
          endif
          enddo
!!!#endif
          end if
!!!#endif
          end if
        endif

        if(nvarcorn.gt.0) then
          do ivar=1,nvarcorn
          if(int_gcell_on_nc(ivar)) then
          gt_unk_n(ivar,:,:,:,1:lnblocks) = unk_n(ivar,:,:,:,1:lnblocks)
          endif
          enddo
        endif

      else

        do lb=1,lnblocks
        if(lrefine(lb).eq.level) then

          if(nvar.gt.0) then
          do ivar=1,nvar
          if(int_gcell_on_cc(ivar)) then
            gt_unk(ivar,:,:,:,lb) = unk(ivar,:,:,:,lb)
          endif
          enddo
          endif
          if(nfacevar.gt.0) then
          do ivar=1,nfacevar
          if(int_gcell_on_fc(1,ivar)) then
            gt_facevarx(ivar,:,:,:,lb) = facevarx(ivar,:,:,:,lb)
          endif
          enddo
!!!#if N_DIM >= 2
            if (ndim >= 2) then
          do ivar=1,nfacevar
          if(int_gcell_on_fc(2,ivar)) then
            gt_facevary(ivar,:,:,:,lb) = facevary(ivar,:,:,:,lb)
          endif
          enddo
!!!#endif
            end if
!!!#if N_DIM == 3
            if (ndim == 3) then
          do ivar=1,nfacevar
          if(int_gcell_on_fc(3,ivar)) then
            gt_facevarz(ivar,:,:,:,lb) = facevarz(ivar,:,:,:,lb)
          endif
          enddo
!!!#endif
            end if
          endif
          if(nvaredge.gt.0) then
!!!#if N_DIM > 1
            if (ndim > 1) then
          do ivar=1,nvaredge
          if(int_gcell_on_ec(1,ivar)) then
            gt_unk_e_x(ivar,:,:,:,lb) = unk_e_x(ivar,:,:,:,lb)
          endif
          enddo
          do ivar=1,nvaredge
          if(int_gcell_on_ec(2,ivar)) then
            gt_unk_e_y(ivar,:,:,:,lb) = unk_e_y(ivar,:,:,:,lb)
          endif
          enddo
!!!#if N_DIM == 3
            if (ndim == 3) then
          do ivar=1,nvaredge
          if(int_gcell_on_ec(3,ivar)) then
            gt_unk_e_z(ivar,:,:,:,lb) = unk_e_z(ivar,:,:,:,lb)
          endif
          enddo
!!!#endif
            end if
!!!#endif
            end if
          endif
          if(nvarcorn.gt.0) then
          do ivar=1,nvarcorn
          if(int_gcell_on_nc(ivar)) then
            gt_unk_n(ivar,:,:,:,lb) = unk_n(ivar,:,:,:,lb)
          endif
          enddo
          endif

        endif
        enddo

      endif                           ! end of level iftest

      call shmem_barrier_all()

! Since this routine is a prelude to calling amr_1blk_guardcell
! we check here to make sure that surrblks has been computed and
! stored for each block.
      if(gsurrblks_set.ne.1) then
        mype = shmem_my_pe()
#ifdef SHMEM_USED
        call amr_gsurrounding_blks(mype,.true.)
#endif
#ifdef MPI_USED

!--------
! This is a temporary section.
! call mpi_amr_gsurr_blks immediately after call to mpi_morton_bnd
! because it uses pe_source and r_mortonbnd which are reset in the
! other morton_bnd_?? routines. Note this is temporary. Will
! redesign so mpi_amr_gsurr_blks is less context sensitive.

! Find the coordinate ranges
         call mpi_amr_global_domain_limits
!
! Compute and save morton number range for each processor
      call mpi_amr_morton_limits(mype)
!
! Set up surrounding blocks of all local blocks (must not precede
! setting of grid_xmin,... etc)
      tag_offset = 100
      if(nprocs.gt.1)
     .    call mpi_morton_bnd(mype,nprocs,tag_offset)
!--------
        nprocs = shmem_n_pes()
        call mpi_amr_gsurr_blks(mype,nprocs)
#endif

      endif

#else  /* NO_PERMANENT_GUARDCELLS */

#ifdef FORCE_CONSISTENCY_AT_SRL_INTERFACES
      if(nfacevar.gt.0) then
        do lb = 1,lnblocks
          do ivar=1,nfacevar
          if(int_gcell_on_fc(1,ivar)) then
          gt_facevarx(ivar,1,:,:,lb) = facevarx(ivar,1+nguard,:,:,lb)
          gt_facevarx(ivar,2,:,:,lb) = 
     .                             facevarx(ivar,nxb+1+nguard,:,:,lb)
          endif
          enddo

          if(ndim.ge.2) then
          do ivar=1,nfacevar
          if(int_gcell_on_fc(2,ivar)) then
          gt_facevary(ivar,:,1,:,lb) = 
     .                             facevary(ivar,:,1+nguard*k2d,:,lb)
          gt_facevary(ivar,:,1+k2d,:,lb) =
     .                       facevary(ivar,:,nyb+(1+nguard)*k2d,:,lb)
          endif
          enddo
          endif

          if(ndim.eq.3) then
          do ivar=1,nfacevar
          if(int_gcell_on_fc(3,ivar)) then
          gt_facevarz(ivar,:,:,1,lb) = 
     .                           facevarz(ivar,:,:,1+nguard*k3d,lb)
          gt_facevarz(ivar,:,:,1+k3d,lb) =
     .                     facevarz(ivar,:,:,nzb+(1+nguard)*k3d,lb)
          endif
          enddo
          endif

        enddo
      endif
#endif  /* FORCE_CONSISTENCY_AT_SRL_INTERFACES */

#endif /* NO_PERMANENT_GUARDCELLS */

      call shmem_barrier_all()

#ifdef TIMING_MPI
              timer_amr_1blk_copy_soln =  timer_amr_1blk_copy_soln
     .                          + mpi_wtime() - time1
#endif /* TIMING_MPI */
      return
      end subroutine amr_1blk_copy_soln
