#include "paramesh_preprocessor.fh"
#include "athena_preprocessor.fh"
C******************************************************************************
      
C  Routine:     plotfile()
      
C  Description: Creates a plotfile file for display using the chombVis package
C               written @ LBL.  The chomboVis package can be downloaded from
C               http://SEESAR.LBL.GOV/anag/chombo/ .

      
C  Author:     K. Olson, 4/30/2001
      
      subroutine plotfile_chombo (NumPes, mype, filenum, simtime)

C Use of this software should credit the chombo project as well as the ASCI/FLASH */
C center                                                                          

C**   _______              __
C**  / ___/ /  ___  __ _  / /  ___
C** / /__/ _ \/ _ \/  ' \/ _ \/ _ \
C** \___/_//_/\___/_/_/_/_.__/\___/ 
C**
C**
C** This software is copyright (C) by the Lawrence Berkeley
C** National Laboratory.  Permission is granted to reproduce
C** this software for non-commercial purposes provided that
C** this notice is left intact.
C** 
C** It is acknowledged that the U.S. Government has rights to
C** this software under Contract DE-AC03-765F00098 between
C** the U.S. Department of Energy and the University of
C** California.
C**
C** This software is provided as a professional and academic
C** contribution for joint exchange.  Thus it is experimental,
C** is provided ``as is'', with no warranties of any kind
C** whatsoever, no support, no promise of updates, or printed
C** documentation.  By using this software, you acknowledge
C** that the Lawrence Berkeley National Laboratory and
C** Regents of the University of California shall have no
C** liability with respect to the infringement of other
C** copyrights by any part of this software.
C**

C==============================================================================

      use physcons
      use physicaldata
      use tree
      use paramesh_dimensions
      use strongfields
      use interior_gbc
      use grid

      use paramesh_interfaces
#ifdef MPI_USED
      use paramesh_mpi_interfaces
#endif /* MPI_USED */


      include 'mpif.h'
      
      integer           filenum
      real              simtime
      
      integer           block_no, jproc, i, j, k
      integer           ivar, iunit1, ix, iy, iz
      integer           iunit2, lnblockst, ngid, ierr
      integer           n_to_left(0:16383) ! must extend from 0 to NumPEs-1
!      integer           gid(mfaces+1+mchild,maxblocks_tr)
      integer, allocatable :: gid(:,:)
      character*4       fnum_string
      character*80      filename
      
      integer           tot_blocks
      integer           status(MPI_STATUS_SIZE)
      integer loop_start_index
      integer nguard0
      parameter(nguard0 = nguard*npgs)

!      integer           lrefinet(maxblocks_tr), nodetypet(maxblocks_tr)
!      integer           gidt(mfaces+1+mchild,maxblocks_tr)
!      real              coordt(mdim,maxblocks_tr)
!      real              sizet(mdim,maxblocks_tr)
!      real              bnd_boxt(2,mdim,maxblocks_tr)

      integer, allocatable :: lrefinet(:), nodetypet(:)
      integer, allocatable :: gidt(:,:)
      real   , allocatable :: coordt(:,:)
      real   , allocatable :: sizet(:,:)
      real   , allocatable :: bnd_boxt(:,:,:)

c set the number of variables we are going to write out
      integer num_out, size_out,bsize_out
      parameter (num_out  = 8,
     +           size_out = (num_out+2)*iu_bnd*ju_bnd*ku_bnd,
     +          bsize_out = (        3)*iu_bnd*ju_bnd*ku_bnd)

      real bnkt(3,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
      real bnkt2(3,nxb,nyb,nzb)

      real unkt(num_out+2,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
      real unkt2(num_out+2,nxb,nyb,nzb)
      real          facevarxt(nbndvar,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,
     &     kl_bnd:ku_bnd)
      real          facevaryt(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,
     &     kl_bnd:ku_bnd)
      real          facevarzt(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,
     &     kl_bnd:ku_bnd+k3d)
      real              xleft, deltax, xpos

      character*80 compNames(num_out+2)

c hold pointers to the location in unk of the variables we are writing out
      integer iout(num_out)

      integer ii, jj, kk
      integer minlevel, maxlevel, numLevels
      integer min_coord(3),max_coord(3)
      integer no_at_level(100), no_at_level_tot(100)
      integer block_tot, block, block_loc
      integer first, level
      real dx(100)
      real x, y, z
      real tmp_real
      integer tmp_int
      real min_domain(mdim)
      real :: rh,ri,vx,vy,vz,pt,p,ee,ek,bx,by,bz,b2,v2,bv

C=============================================================================

      allocate(gid(mfaces+1+mchild,maxblocks_tr))
      allocate(lrefinet(maxblocks_tr))
      allocate(nodetypet(maxblocks_tr))
      allocate(gidt(mfaces+1+mchild,maxblocks_tr))
      allocate(coordt(mdim,maxblocks_tr))
      allocate(sizet(mdim,maxblocks_tr))
      allocate(bnd_boxt(2,mdim,maxblocks_tr))




c set the variables we are going to store

      iout(1) = 1
      iout(2) = 2
      iout(3) = 3
      iout(4) = 4
      iout(5) = 5
      iout(6) = 6
      iout(7) = 7
      iout(8) = 8

      compNames(1) = 'Density'
      compNames(2) = 'X-velocity'
      compNames(3) = 'Y-velocity'
      compNames(4) = 'Z-velocity'
      compNames(5) = 'Pressure'
      compNames(6) = 'Bx'
      compNames(7) = 'By'
      compNames(8) = 'Bz'

      compNames(num_out+1) = 'Block'
      compNames(num_out+2) = 'Proc'

C COMPUTE TOTAL NO. OF BLOCKS

      lnblockst = 0
      do block_no = 1,lnblocks
	 if (.true..or.nodetype(block_no) == 1) lnblockst = lnblockst + 1
      end do
      call MPI_ALLREDUCE (lnblockst, tot_blocks, 1,
     &     MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)

C Find min domain coordinates

      min_domain(:) = 1.e20
      do block = 1, lnblocks
         do i = 1,ndim
            x = coord(i,block) - bsize(i,block)/2.
            min_domain(i) = min(x,min_domain(i))
         end do
      end do
      do i = 1,ndim
          call MPI_ALLREDUCE (min_domain(i), tmp_real, 1,
     &        MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierr)
         min_domain(i) = tmp_real
      end do

C COMPUTE minlevel, maxlevel, and min and max coordinates that chombo wants

      min_coord(:) =  10000000
      max_coord(:) = -10000000
      minlevel     =  1000000
      maxlevel     = -1
      no_at_level_tot(:) = 0
      do block = 1, lnblocks

         if (.true..or.nodetype(block) == 1) then

         if (lrefine(block) < minlevel) minlevel = lrefine(block)
         if (lrefine(block) > maxlevel) maxlevel = lrefine(block)


         no_at_level_tot(lrefine(block)) = 
     &        no_at_level_tot(lrefine(block)) + 1

	 end if

         if (lrefine(block) == 1) then

            x  = coord(1,block) - bsize(1,block)/2.
            x  = x - min_domain(1)
            x  = x + bsize(1,block)/(2*nxb)
            ix = int(x*nxb/bsize(1,block))

            y  = coord(2,block) - bsize(2,block)/2.
            y  = y - min_domain(2)
            y  = y + bsize(2,block)/(2*nyb)
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

      call MPI_ALLREDUCE (minlevel, tmp_int, 1,
     &     MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
      minlevel = tmp_int

      call MPI_ALLREDUCE (maxlevel, tmp_int, 1,
     &     MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
      maxlevel = tmp_int

      numLevels = maxlevel - minlevel + 1

      do i = 1,ndim

         call MPI_ALLREDUCE (min_coord(i), tmp_int, 1,
     &        MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierr)
         min_coord(i) = tmp_int

         call MPI_ALLREDUCE (max_coord(i), tmp_int, 1,
     &        MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
         max_coord(i) = tmp_int

      end do

      do level = 1,numLevels
         call MPI_ALLREDUCE (no_at_level_tot(level), tmp_int, 1,
     &        MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
         no_at_level_tot(level) = tmp_int
      end do
      
      dx(1) = 1./real(nxb)
      do level = 2,numLevels
         dx(level) = dx(level-1)/2.
      end do

! Add in strong fields
       do block_no = 1,lnblocks
          if (.true..or.nodetype(block_no) == 1) then
          call grid_variables(block_no,dt)
          call strongfields_wb(mype,block_no)
CYYC
c                  b_sf=0.

          kk = nguard*k3d+1
          do k = nguard0*k3d+1,nguard0*k3d+nzb
             jj = nguard*k2d+1
             do j = nguard0*k2d+1,nguard0*k2d+nyb
                ii = nguard+1
                do i = nguard0+1,nguard0+nxb
                   do ivar = 6,8
                      unk(ivar,i,j,k,block_no) = unk(ivar,i,j,k,block_no) !cyyc+b_sf(ivar-5,ii,jj,kk)

                      
                   end do
                  ii = ii + 1
                end do
                jj = jj + 1
             end do
             kk = kk + 1
          end do
          end if
       end do
       call shmem_barrier_all()
      
C NOW WRITE OUT THE DATA FROM PROC 0

      if (MyPe .eq. 0) then

C       Determine the filename. 
         write(*,*)'within plotfile_chombo: filenum=',filenum
         write (fnum_string, '(i4.4)') filenum
         filename = './chombo_plt_cnt_' // fnum_string // '.hdf5'

      end if

! loop through the blocks and write them out !!!
      block_tot = 0
      no_at_level(:) = 0
      first = 1
      do jproc = 0,NumPEs-1

         if (MyPe.eq.0) then

C fetch lnblocks from other processors

            if (jproc.ne.0) then

            call MPI_RECV (lnblockst,1,MPI_INTEGER,jproc,
     &              1,MPI_COMM_WORLD,status,ierr)
            if (lnblockst.gt.0) then
               call MPI_RECV (lrefinet(1),lnblockst,MPI_INTEGER,jproc,
     &              2,MPI_COMM_WORLD,status,ierr)
               call MPI_RECV (nodetypet(1),lnblockst,MPI_INTEGER,jproc,
     &              3,MPI_COMM_WORLD,status,ierr)
               call MPI_RECV (coordt(1,1),mdim*lnblockst,
     &              MPI_DOUBLE_PRECISION,
     &              jproc,4,MPI_COMM_WORLD,status,ierr)
               call MPI_RECV (sizet(1,1),mdim*lnblockst,
     &              MPI_DOUBLE_PRECISION,
     &              jproc,5,MPI_COMM_WORLD,status,ierr)
               call MPI_RECV (bnd_boxt(1,1,1),2*mdim*lnblockst,
     &              MPI_DOUBLE_PRECISION,
     &              jproc,6,MPI_COMM_WORLD,status,ierr)
               call MPI_RECV (gidt(1,1),lnblockst*(mfaces+1+mchild),
     &              MPI_INTEGER,
     &              jproc,7,MPI_COMM_WORLD,status,ierr)
            end if

            else
               
               lnblockst = lnblocks
               lrefinet(1:lnblockst)      = lrefine(1:lnblockst)
               nodetypet(1:lnblockst)     = nodetype(1:lnblockst)
               coordt(:,1:lnblockst)      = coord(:,1:lnblockst)
               sizet(:,1:lnblockst)       = bsize(:,1:lnblockst)
               bnd_boxt(:,:,1:lnblockst)  = bnd_box(:,:,1:lnblockst)
               gidt(:,1:lnblockst)        = gid(:,1:lnblockst)

            end if

            do block_no = 1,lnblockst

C fetch data for this block

               if (jproc.ne.0) then
                  call MPI_RECV (unkt(1,1,1,1),
     &                 size_out,
     &                 MPI_DOUBLE_PRECISION,
     &                 jproc,7+block_no,
     &                 MPI_COMM_WORLD,status,ierr)

                  call MPI_RECV (bnkt(1,1,1,1),
     &                 bsize_out,
     &                 MPI_DOUBLE_PRECISION,
     &                 jproc,100*lnblockst+block_no,
     &                 MPI_COMM_WORLD,status,ierr)

!                  if (nfacevar.gt.0) then
!                  call MPI_RECV (facevarxt(1,1,1,1),
!     &                 nbndvar*len_blockfx,
!     &                 MPI_DOUBLE_PRECISION,
!     &                 jproc,7+block_no+lnblockst,
!     &                 MPI_COMM_WORLD,status,ierr)
!                  call MPI_RECV (facevaryt(1,1,1,1),
!     &                 nbndvar*len_blockfy,
!     &                 MPI_DOUBLE_PRECISION,
!     &                 jproc,7+block_no+2*lnblockst,
!     &                 MPI_COMM_WORLD,status,ierr)
!                  call MPI_RECV (facevarzt(1,1,1,1),
!     &                 nbndvar*len_blockfz,
!     &                 MPI_DOUBLE_PRECISION,
!     &                 jproc,7+block_no+3*lnblockst,
!     &                 MPI_COMM_WORLD,status,ierr)
!                  end if
               else



cxxc                  do i = 1, num_out+1
cxxc                     if (i <= num_out) then
cxxc                        unkt(i,:,:,:) = unk(iout(i),:,:,:,block_no)
cxxc                     else
cxxc                        unkt(i,:,:,:) = block_no
cxxc                        unkt(i+1,:,:,:) = jproc
cxxc                     endif
cxxc                  enddo

                   do i = 1, num_out
                     unkt(i,:,:,:) = unk(iout(i),:,:,:,block_no)
                  end do

                  call grid_variables(block_no,dt)
                  call strongfields_wb(mype,block_no)
                  kk = nguard*k3d+1
                  do k = nguard0*k3d+1,nguard0*k3d+nzb
                     jj = nguard*k2d+1
                     do j = nguard0*k2d+1,nguard0*k2d+nyb
                        ii = nguard+1
                        do i = nguard0+1,nguard0+nxb
                           do ivar = 1,3
                              bnkt(ivar,i,j,k) = b_sf(ivar,ii,jj,kk)

                           end do
                           ii = ii + 1
                        end do
                        jj = jj + 1
                     end do
                     kk = kk + 1
                  end do
                
                  if (nfacevar.gt.0) then
                     facevarxt(:,:,:,:) = facevarx(:,:,:,:,block_no)
                     facevaryt(:,:,:,:) = facevary(:,:,:,:,block_no)
                     facevarzt(:,:,:,:) = facevarz(:,:,:,:,block_no)
                  end if
               end if
               
               kk = nguard*k3d*npgs+1
               do k = 1,nzb
                  jj = nguard*k2d*npgs+1
                  do j = 1,nyb
                     ii = nguard*npgs+1
                     do i = 1,nxb
cxxc                        do ivar = 1,num_out+2
                        do ivar = 1,num_out
                           unkt2(ivar,i,j,k) = unkt(ivar,ii,jj,kk)
                        end do
                        do ivar = 1,3
                           bnkt2(ivar,i,j,k) = bnkt(ivar,ii,jj,kk)
                        end do
                        ii = ii + 1
                     end do
                     jj = jj + 1
                  end do
                  kk = kk + 1
               end do
	       

               if (.true..and.nodetypet(block_no) == 1) then
c
c------convert to primitives
c
c               write(*,*)'entering conversion to primitives'
c               write(*,*)'k1d,k2,k3d=',k1d,k2d,k3d

c                  call grid_variables(block_no,dt)
c                  call strongfields_wb(mype,block_no)

CYYC
c                  b_sf=0.

	       kk = nguard*k3d+1
                do k = 1,nzb
                 jj = nguard*k2d+1
                  do j = 1,nyb
                   ii = nguard+1
                    do i = 1,nxb
                       rh = unkt2(1,i,j,k)
                       px = unkt2(2,i,j,k) 
                       py = unkt2(3,i,j,k) 
                       pz = unkt2(4,i,j,k) 
                       ee = unkt2(5,i,j,k) 
c                       write(*,*) 'plotfile, b_sf', mype,block_no, i,j,bnkt2(1,i,j,k)
                       bx = unkt2(6,i,j,k) + bnkt2(1,i,j,k)
                       by = unkt2(7,i,j,k) + bnkt2(2,i,j,k)
                       bz = unkt2(8,i,j,k) + bnkt2(3,i,j,k)
                       bxd= unkt2(6,i,j,k)
                       byd= unkt2(7,i,j,k)
                       bzd= unkt2(8,i,j,k)
                       b2 = bx*bx+by*by+bz*bz 
                       b2d= bxd*bxd+byd*byd+bzd*bzd 
                       ri = 1.0/(rh + afact*b2)
                       bv = (px*bx+py*by+pz*bz)*ri
                       vx = (px + afact*bv*bx)*ri
                       vy = (py + afact*bv*by)*ri
                       vz = (pz + afact*bv*bz)*ri
                       v2 = vx*vx+vy*vy+vz*vz 
                       ek = 0.5*rh*v2
                       b2d= pie8i*b2d
                       pt = ee - ek - b2d 
                       p  = max(gammam1*pt,pmin)
                       unkt2(2,i,j,k) = vx
                       unkt2(3,i,j,k) = vy
                       unkt2(4,i,j,k) = vz
                       unkt2(5,i,j,k) = p
                       unkt2(6,i,j,k) = bx
                       unkt2(7,i,j,k) = by
                       unkt2(8,i,j,k) = bz
CYYC
c                       unkt2(6,i,j,k) = bxd
c                       unkt2(7,i,j,k) = byd
c                       unkt2(8,i,j,k) = bzd

cyyc
c                       if(jproc.eq.4.and.block_no.eq.2)then
c                          write(*,*) 'plotfile_chombo', p
c                       endif


                       ii = ii + 1
                    end do
                    jj = jj + 1
                 end do
                 kk = kk + 1
              end do

c              write(*,*)'leaving conversion to primitives'
              end if

	      if (.true..or.nodetypet(block_no) == 1) then

               no_at_level(lrefinet(block_no)) =
     &              no_at_level(lrefinet(block_no)) + 1
               block_loc = block_no - 1
               
c               write(*,*) nxb,nyb,nzb,
c     &              lrefinet,min_coord,max_coord
c               stop


               if (nodetypet(block_no) .ne. 1) then
c                  unkt2=4.77
               end if

               call writeblocktochombofile (block_loc,
     &              block_tot,
     &              no_at_level(lrefinet(block_no)),
     &              filename,tot_blocks,numLevels,
     &              ndim,num_out+2,nxb,nyb,nzb,
     &              lrefinet,min_coord,max_coord,
     &              coordt, sizet,unkt2,simtime,
     &              dt,first,
     &              no_at_level_tot,
     &              dx,
     &              compNames,
     &              min_domain)
               block_tot = block_tot + 1

	       end if
               
            end do              ! end loop over lnblockst
            
         else                   ! if (MyPe.eq.0

            if (jproc.eq.MyPe) then

               call MPI_SEND (lnblocks,1,MPI_INTEGER,0,
     &              1,MPI_COMM_WORLD,ierr)
               if (lnblocks.gt.0) then
                  call MPI_SEND (lrefine(1),lnblocks,MPI_INTEGER,0,
     &                 2,MPI_COMM_WORLD,ierr)
                  call MPI_SEND (nodetype(1),lnblocks,MPI_INTEGER,0,
     &                 3,MPI_COMM_WORLD,ierr)
                  call MPI_SEND (coord(1,1),mdim*lnblocks,
     &                 MPI_DOUBLE_PRECISION,0,
     &                 4,MPI_COMM_WORLD,ierr)
                  call MPI_SEND (bsize(1,1),mdim*lnblocks,
     &                 MPI_DOUBLE_PRECISION,0,
     &                 5,MPI_COMM_WORLD,ierr)
                  call MPI_SEND (bnd_box(1,1,1),2*mdim*lnblocks,
     &                 MPI_DOUBLE_PRECISION,
     &                 0,6,MPI_COMM_WORLD,ierr)
                  call MPI_SEND (gid(1,1),lnblocks*(mfaces+1+mchild),
     &                 MPI_INTEGER,
     &                 0,7,MPI_COMM_WORLD,ierr)
                  do block_no = 1,lnblocks
                     do i = 1, num_out+1
                        if (i <= num_out) then
                           unkt(i,:,:,:) = unk(iout(i),:,:,:,block_no)
                        else
                           unkt(i,:,:,:) = block_no
                           unkt(i+1,:,:,:) = jproc
                        endif
                     enddo
                     call MPI_SEND (unkt(1,1,1,1),
     &                    size_out,
     &                    MPI_DOUBLE_PRECISION,
     &                    0,7+block_no,
     &                    MPI_COMM_WORLD,ierr)


                     call grid_variables(block_no,dt)
                     call strongfields_wb(mype,block_no)
                     kk = nguard*k3d+1
                     do k = nguard0*k3d+1,nguard0*k3d+nzb
                        jj = nguard*k2d+1
                        do j = nguard0*k2d+1,nguard0*k2d+nyb
                           ii = nguard+1
                           do i = nguard0+1,nguard0+nxb
                              do ivar = 1,3
                                 bnkt(ivar,i,j,k) = b_sf(ivar,ii,jj,kk)
                              end do
                              ii = ii + 1
                           end do
                           jj = jj + 1
                        end do
                        kk = kk + 1
                     end do
                     call MPI_SEND (bnkt(1,1,1,1),
     &                    bsize_out,
     &                    MPI_DOUBLE_PRECISION,
     &                    0,100*lnblockst+block_no,
     &                    MPI_COMM_WORLD,ierr)

!                     if (nfacevar.gt.0) then
!                        call MPI_SEND (facevarx(1,1,1,1,block_no),
!     &                       nbndvar*len_blockfx,
!     &                       MPI_DOUBLE_PRECISION,
!     &                       0,7+block_no+lnblocks,
!     &                       MPI_COMM_WORLD,status,ierr)
!                        call MPI_SEND (facevary(1,1,1,1,block_no),
!     &                       nbndvar*len_blockfy,
!     &                       MPI_DOUBLE_PRECISION,
!     &                       0,7+block_no+2*lnblocks,
!     &                       MPI_COMM_WORLD,status,ierr)
!                        call MPI_SEND (facevarz(1,1,1,1,block_no),
!     &                       nbndvar*len_blockfz,
!     &                       MPI_DOUBLE_PRECISION,
!     &                       0,7+block_no+3*lnblocks,
!     &                       MPI_COMM_WORLD,status,ierr)
!                     end if
                  end do
               end if
               
            end if

         end if                 ! if MyPe == 0
      
         call MPI_BARRIER(MPI_COMM_WORLD,ierr)

      end do  ! end loop over processors

!      if (MyPe.eq.0) then
!         close(iunit1)
!      end if
! substract in strong fields
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
                      unk(ivar,i,j,k,block_no) = unk(ivar,i,j,k,block_no) !cyyc-b_sf(ivar-5,ii,jj,kk)
                   end do
                  ii = ii + 1
                end do
                jj = jj + 1
             end do
             kk = kk + 1
          end do
       end do
      if (mype.eq.0) then
         print *,'**** wrote file ',filename
      end if

      deallocate(gid)
      deallocate(lrefinet)
      deallocate(nodetypet)
      deallocate(gidt)
      deallocate(coordt)
      deallocate(sizet)
      deallocate(bnd_boxt)

      call MPI_BARRIER(MPI_COMM_WORLD,ierr)


      return

      end



