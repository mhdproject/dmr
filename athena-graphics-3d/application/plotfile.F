#include "paramesh_preprocessor.fh"
#include "athena_preprocessor.fh"
#undef CRAY_T3E
      subroutine plotfile (time, iunit1)


! $RCSfile: plotfile.F,v $
! $Revision: 1.1.1.1 $
! $Date: 2003/11/07 16:04:24 $


! Subroutine to checkpoint runs using AMR package.
! Writes out tree data structure and data stored in blocks
! Currently writes are done serially by processor 0 collecting data from
! other processors and then writing it out.
! USES UNFORMATTED DIRECT I/O

      use physcons
      use paramesh_dimensions
      use physicaldata
      use tree
      use strongfields
      use interior_gbc
      use grid

      use paramesh_interfaces
#ifdef MPI_USED
      use paramesh_mpi_interfaces
#endif /* MPI_USED */

      implicit none

#include "amr_shmem.fh"
#ifdef MPI_USED
      include "mpif.h"

      integer :: tag_offset,max_blks_sent
#endif /* MPI_USED */


      integer loop_start_index
      integer nguard0
      parameter(nguard0 = nguard*npgs)


      integer block_no,shmem_n_pes,shmem_my_pe
      integer jproc,i,j,ivar,iunit1,ix,iy,iz,nprocs,mype,iunit2

      integer lnblockst
      integer ierr

! TEMPORARIES WHICH ARE WRITTEN !
      
      integer  tot_blocks
      integer  lrefinet(maxblocks)
      real  coordt(mdim,maxblocks),sizet(mdim,maxblocks)
      real  unkt(nvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
      real  unkt2(nvar,nxb,nyb,nzb)

      save lrefinet,lnblockst,coordt,sizet

      real xleft,deltax,xpos
      real simtime,dt
      save unkt
      real :: dtold,time
      common/old_timestep/dtold
      real :: rl,vxl,vyl,vzl,pl,bxl,byl,bzl
      common/left_state/rl,vxl,vyl,vzl,pl,bxl,byl,bzl

      character*4 fnum_string
      character*80 compNames(nvar)
      character*80 filename
      integer minlevel, maxlevel, numLevels
      integer min_coord(3),max_coord(3)
      integer no_at_level(100), no_at_level_tot(100)
      integer first, level
      real dx(100)
      real x, y, z
      real tmp_real
      real min_domain(mdim)
      integer tmp_int
      integer k, ii, jj, kk
      integer errcode
      integer block, block_tot, block_loc

      save tot_blocks, tmp_real, min_domain, tmp_int, minlevel, maxlevel
      save dx, no_at_level_tot

      call shmem_barrier_all()

      nprocs = shmem_n_pes()
      mype  = shmem_my_pe()

c set the variables we are going to store

!      do i = 1,nvar
!         write (fnum_string, '(i4.4)') i
!         compNames(i) = 'var#' // fnum_string
!      end do

      compNames(1) = 'Density'
      compNames(2) = 'X-momentum'
      compNames(3) = 'Y-momentum'
      compNames(4) = 'Z-momentum'
      compNames(5) = 'Energy-Density'
      compNames(6) = 'Bx'
      compNames(7) = 'By'
      compNames(8) = 'Bz'

      call shmem_barrier_all()

C COMPUTE TOTAL NO. OF BLOCKS

      call comm_int_sum_to_all(tot_blocks, lnblocks)
      call shmem_barrier_all()

C Find min domain coordinates

      min_domain(:) = 1.e20
      do block = 1, lnblocks
         do i = 1,ndim
            x = coord(i,block) - bsize(i,block)/2.
            min_domain(i) = min(x,min_domain(i))
         end do
      end do
      call shmem_barrier_all()
      do i = 1,ndim
         call comm_real_min_to_all(tmp_real, min_domain(i))
         call shmem_barrier_all()
         min_domain(i) = tmp_real
      end do

C COMPUTE minlevel, maxlevel, and min and max coordinates that chombo wants

      min_coord(:) = 10000000
      max_coord(:) = -10000000
      minlevel = 1000000
      maxlevel = -1
      no_at_level_tot(:) = 0
      do block = 1, lnblocks

         if (lrefine(block) < minlevel) minlevel = lrefine(block)
         if (lrefine(block) > maxlevel) maxlevel = lrefine(block)

         no_at_level_tot(lrefine(block)) =
     &        no_at_level_tot(lrefine(block)) + 1

         if (lrefine(block) == 1) then
 
            x = coord(1,block) - bsize(1,block)/2.
            x = x - min_domain(1)
            x = x + bsize(1,block)/(2*nxb)
            ix = int(x*nxb/bsize(1,block))
 
            y = coord(2,block) - bsize(2,block)/2.
            y = y - min_domain(2)
            y = y + bsize(2,block)/(2*nyb)
            iy = int(y*nyb/bsize(2,block))
 
#if N_DIM == 3
            z = coord(3,block) - bsize(3,block)/2.
            z = z - min_domain(3)
            z = z + bsize(3,block)/(2*nzb)
            iz = int(z*nzb/bsize(3,block))
#endif
            if (ix <= min_coord(1)) min_coord(1) = ix
            if (ix >= max_coord(1)) max_coord(1) = ix
            if (iy <= min_coord(2)) min_coord(2) = iy
            if (iy >= max_coord(2)) max_coord(2) = iy
#if N_DIM == 3
            if (iz <= min_coord(3)) min_coord(3) = iz
            if (iz >= max_coord(3)) max_coord(3) = iz
#endif
         end if            
      end do

      call shmem_barrier_all()
      call comm_int_min_to_all(tmp_int, minlevel)
      minlevel = tmp_int

      call shmem_barrier_all()
      call comm_int_max_to_all(tmp_int, maxlevel)
      maxlevel = tmp_int

      numLevels = maxlevel - minlevel + 1

      do i = 1,ndim

         call shmem_barrier_all()
         call comm_int_min_to_all(tmp_int, min_coord(i))
         min_coord(i) = tmp_int

         call shmem_barrier_all()
         call comm_int_max_to_all(tmp_int, max_coord(i))
         max_coord(i) = tmp_int

       end do

       do level = 1,numLevels
          call shmem_barrier_all()
          call comm_int_sum_to_all(tmp_int, no_at_level_tot(level))
          no_at_level_tot(level) = tmp_int
       end do

       dx(1) = 1./nxb
       do level = 2,numLevels
          dx(level) = dx(level-1)/2.
       end do

       call shmem_barrier_all()

! Add in strong fields
       do block_no = 1,lnblocks
          call grid_variables(block_no,dt)
          call strongfields_wb(mype,block_no)
          kk = nguard*k3d+1
          do k = nguard0*k3d+1,nguard0*k3d+nzb
             jj = nguard*k2d+1
             do j = nguard0*k2d+1,nguard0*k2d+nyb
                ii = nguard+1
                do i = nguard0+1,nguard0+nxb
                   do ivar = 6,8
                      unk(ivar,i,j,k,block_no) = 
     &                     unk(ivar,i,j,k,block_no) +
     &                     b_sf(ivar-5,ii,jj,kk)
                   end do
                  ii = ii + 1
                end do
                jj = jj + 1
             end do
             kk = kk + 1
          end do
       end do
       call shmem_barrier_all()

! NOW WRITE OUT THE DATA FROM PROC 0

       if (mype .eq. 0) then
          
          write (fnum_string, '(i4.4)') iunit1
          filename = 'chombo_' // fnum_string
          
          block_tot = 0
          no_at_level(:) = 0
          first = 1
          
          do jproc = 0,nprocs-1

! fetch lnblocks from other processors

            call SHMEM_INTEGER_GET (lnblockst,lnblocks,1,jproc)
            call SHMEM_INTEGER_GET (lrefinet,lrefine,lnblockst,jproc)
            call SHMEM_REAL_GET (coordt,coord,mdim*lnblockst,jproc)
            call SHMEM_REAL_GET (sizet,bsize,mdim*lnblockst,jproc)
            
            do block_no = 1,lnblockst

! fetch data for this block
               call SHMEM_REAL_GET (unkt,unk(1,1,1,1,block_no),
     $              len_block,jproc)

               kk = nguard0*k3d+1
               do k = 1,nzb
                  jj = nguard0*k2d+1
                  do j = 1,nyb
                     ii = nguard0+1
                     do i = 1,nxb
                        do ivar = 1,nvar
                           unkt2(ivar,i,j,k) = unkt(ivar,ii,jj,kk)
                        end do
                        ii = ii + 1
                     end do
                     jj = jj + 1
                  end do
                  kk = kk + 1
               end do

               no_at_level(lrefinet(block_no)) =
     &              no_at_level(lrefinet(block_no)) + 1
               block_loc = block_no - 1                                            
               call writeblocktochombofile (block_loc, block_tot,
     &              no_at_level(lrefinet(block_no)), filename,
     &              tot_blocks,
     &              numLevels, ndim, nvar, nxb, nyb, nzb, lrefinet,
     &              min_coord, max_coord, coordt, sizet, unkt2, time,
     &              dt, first, no_at_level_tot, dx, compNames,
     &              min_domain)
               block_tot = block_tot + 1

            end do              ! end loop over fetched blocks
         end do                 ! end loop over jproc
      end if                    ! if mype == 0
      
      call shmem_barrier_all()

! Subtract off strong fields
      do block_no = 1,lnblocks
         call grid_variables(block_no,dt)
         call strongfields_wb(mype,block_no)
         kk = nguard*k3d+1
         do k = nguard0*k3d+1,nguard0*k3d+nzb
            jj = nguard*k2d+1
            do j = nguard0*k2d+1,nguard0*k2d+nyb
               ii = nguard+1
               do i = nguard0+1,nguard0+nxb
                  do ivar = 6,8
                     unk(ivar,i,j,k,block_no) = 
     &                    unk(ivar,i,j,k,block_no) -
     &                    b_sf(ivar-5,ii,jj,kk)
                  end do
                  ii = ii + 1
               end do
               jj = jj + 1
            end do
            kk = kk + 1
         end do
      end do
      call shmem_barrier_all()

      return
      end subroutine plotfile

