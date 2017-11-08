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

      subroutine amr_checkpoint_wr(iunit1,l_with_guardcells)


! Subroutine to checkpoint runs using AMR package.
! Writes out tree data structure and data stored in blocks
! Currently writes are done serially by processor 0 collecting data from
! other processors and then writing it out.
! USES UNFORMATTED DIRECT I/O

! Written: K. Olson and C. Mobarry 7/97
! Modified: P. MacNeice 4/99


      use paramesh_dimensions
      use physicaldata
      use tree
      use timings

      implicit none

#include "amr_shmem.fh"

#ifdef MPI_USED
      include 'mpif.h'
#endif /* MPI_USED */

      integer, intent(in) :: iunit1
      logical, optional, intent(in)  :: l_with_guardcells

      integer nguard0
      parameter(nguard0 = nguard*npgs)


      integer block_no,shmem_n_pes,shmem_my_pe
      integer jproc,i,j,ivar,ix,iy,iz,nprocs,iproc

      integer lnblockst
      integer ngid
      integer,dimension (:),  allocatable :: n_to_left

      integer gid(nfaces+1+nchild,maxblocks_tr)
      integer lrefinet(maxblocks_tr),nodetypet(maxblocks_tr)
      integer which_childt(maxblocks_tr)
      integer gidt(nfaces+1+nchild,maxblocks_tr)
      integer bflagst(mflags,maxblocks_tr)
      real  coordt(mdim,maxblocks_tr)
      real  work_blockt(maxblocks_tr)
      real  bnd_boxt(2,mdim,maxblocks_tr)
      integer ierr

! TEMPORARIES WHICH ARE WRITTEN !
      
      integer  tot_blocks

      real  unkt(nvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
      real  facevarxt(nbndvar,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,
     &                kl_bnd:ku_bnd)
      real  facevaryt(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,
     &                kl_bnd:ku_bnd)
      real  facevarzt(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,
     &                kl_bnd:ku_bnd+k3d)
      real  :: unk_nt(nbndvarc,
     &               il_bnd:iu_bnd+1,
     &               jl_bnd:ju_bnd+k2d,
     &               kl_bnd:ku_bnd+k3d)
      real  :: unk_e_xt(nbndvare,
     &                  il_bnd:iu_bnd,
     &                  jl_bnd:ju_bnd+k2d,
     &                  kl_bnd:ku_bnd+k3d)
      real  :: unk_e_yt(nbndvare,
     &                  il_bnd:iu_bnd+1,
     &                  jl_bnd:ju_bnd,
     &                  kl_bnd:ku_bnd+k3d)
      real  :: unk_e_zt(nbndvare,
     &                  il_bnd:iu_bnd+1,
     &                  jl_bnd:ju_bnd+k2d,
     &                  kl_bnd:ku_bnd)
      integer :: il0,iu0,jl0,ju0,kl0,ku0
      integer :: ion_c,ion_f,ion_e,ion_n,iv_c,iv_f,iv_e,iv_n


#ifdef MPI_USED
      integer,dimension (:),  allocatable :: glnblocks
      integer :: isrc,idest,itag,isize,ierror,position

      integer status(MPI_STATUS_SIZE)

      integer :: buf_dim_int
      integer :: buf_dim_real
      integer :: buf_dim1, buf_dim2 
      integer :: no_of_bytes_per_real,no_of_bytes_per_integer
      integer :: buf_dim_bytes1,buf_dim_bytes2 

      real,allocatable :: CS_buffer1(:), CR_buffer1(:)
      real,allocatable :: CS_buffer2(:), CR_buffer2(:)

      integer :: nvar_chk_cc,nvar_chk_fc,nvar_chk_ec,nvar_chk_nc


      buf_dim_int = 
     .              maxblocks*( 3+mflags+(nfaces+1+nchild))
      buf_dim_real = 
     .              maxblocks*( 3*mdim + 1)

      buf_dim1 = buf_dim_real + buf_dim_int

      allocate(CS_buffer1(buf_dim1))
      allocate(CR_buffer1(buf_dim1))
#endif /* MPI_USED */

#ifdef CHECKP_ON
      nvar_chk_cc =  0
      do i=1,nvar
        if(checkp_on_cc(i)) nvar_chk_cc = nvar_chk_cc + 1
      enddo
      nvar_chk_fc =  0
      do i=1,nfacevar
        if(checkp_on_fc(1,i)) nvar_chk_fc = nvar_chk_fc + 1
      enddo
      nvar_chk_ec =  0
      do i=1,nvaredge
        if(checkp_on_ec(1,i)) nvar_chk_ec = nvar_chk_ec + 1
      enddo
      nvar_chk_nc =  0
      do i=1,nvarcorn
        if(checkp_on_nc(i)) nvar_chk_nc = nvar_chk_nc + 1
      enddo
#endif /* CHECKP_ON */

#ifdef MPI_USED
      buf_dim2 = len_block 
     .       + nbndvar*(len_blockfx + len_blockfy + len_blockfz)
     .       + nbndvarc*len_blockn 
     .       + nbndvare*(len_blockex + len_blockey + len_blockez)

      allocate(CS_buffer2(buf_dim2))
      allocate(CR_buffer2(buf_dim2))

#endif /* MPI_USED */

#ifdef DEBUG
      write(*,*) 'entered checkpoint_wr '
#endif /* DEBUG */

      nprocs = shmem_n_pes()
      iproc  = shmem_my_pe()


      if(ndim.lt.3) then
        bnd_box(1,3,:) = 0.
        bnd_box(2,3,:) = 1.
        coord(3,:) = .5*(bnd_box(2,3,:)+bnd_box(1,3,:))
      endif
      if(ndim.lt.2) then
        bnd_box(1,2,:) = 0.
        bnd_box(2,2,:) = 1.
        coord(2,:) = .5*(bnd_box(2,2,:)+bnd_box(1,2,:))
      endif



! COMPUTE TOTAL NO. OF BLOCKS STORED TO THE 'LEFT' OF THIS PROCESSOR

      if(allocated(n_to_left)) deallocate( n_to_left )
      allocate ( n_to_left(0:nprocs-1) )

#ifdef MPI_USED
      if(allocated(glnblocks)) deallocate( glnblocks )
      allocate ( glnblocks(0:nprocs-1) )

      call mpi_type_size(MPI_INTEGER,no_of_bytes_per_integer,ierr)
      call mpi_type_size(MPI_REAL
     .                  ,no_of_bytes_per_real   ,ierr)
      buf_dim_bytes1 = buf_dim_real*no_of_bytes_per_real +
     .                 buf_dim_int *no_of_bytes_per_integer
      buf_dim_bytes2 = buf_dim2*no_of_bytes_per_real

      glnblocks(iproc) = lnblocks
      call MPI_Allgather(glnblocks(iproc), 1,MPI_INTEGER,
     .                   glnblocks,1,MPI_INTEGER,
     .                   MPI_COMM_WORLD,ierror)
      n_to_left = glnblocks
#ifdef DEBUG
      write(*,*) 'n_to_left ',n_to_left
#endif /* DEBUG */
#else /* MPI_USED */
      do i = 0,nprocs-1
         call SHMEM_INTEGER_GET(n_to_left(i),lnblocks,1,i)
      end do
#endif /* MPI_USED */

      tot_blocks = 0
      do i = 0,nprocs-1
         tot_blocks = tot_blocks + n_to_left(i)
      end do
#ifdef DEBUG
      write(*,*) 'pe ',iproc,' tot_blocks ',tot_blocks
#endif /* DEBUG */
            
      do i = nprocs-1,1,-1
         n_to_left(i) = n_to_left(i-1)
      end do

      n_to_left(0) = 0
      do i = 2,nprocs-1
         n_to_left(i) = n_to_left(i) + n_to_left(i-1)
      end do

#ifdef DEBUG
      write(*,*) 'pe ',iproc,' n_to_left ',n_to_left
#endif /* DEBUG */

! COMPUTE GLOBAL INDIRECT ADDRESSES FOR TREE DATA (gid)

      do block_no = 1,lnblocks

         ngid = 0
         do j = 1,nfaces
            ngid = ngid + 1
            if (neigh(1,j,block_no).gt.0) then
               gid(ngid,block_no) = neigh(1,j,block_no) + 
     $              n_to_left(neigh(2,j,block_no))
            else
               gid(ngid,block_no) = neigh(1,j,block_no)
            end if
         end do
         
         ngid = ngid + 1
         if (parent(1,block_no).gt.0) then
            gid(ngid,block_no) = parent(1,block_no) + 
     $           n_to_left(parent(2,block_no))
         else
            gid(ngid,block_no) = parent(1,block_no)
         end if
         
         do j = 1,nchild
            ngid = ngid + 1
            if (child(1,j,block_no).gt.0) then
               gid(ngid,block_no) = child(1,j,block_no) + 
     $              n_to_left(child(2,j,block_no))
            else
               gid(ngid,block_no) = child(1,j,block_no)
            end if
         end do

      end do

! NOW WRITE OUT THE DATA FROM PROC 0

! set limits on data arrays
       il0 = nguard0
       iu0 = nxb-1+nguard0
       jl0 = nguard0*k2d
       ju0 = (nyb-1+nguard0)*k2d
       kl0 = nguard0*k3d
       ku0 = (nzb-1+nguard0)*k3d
#ifndef NO_PERMANENT_GUARDCELLS
       if(l_with_guardcells) then
         il0 = 0
         iu0 = nxb+2*nguard0-1
         jl0 = 0
         ju0 = (nyb+2*nguard0-1)*k2d
         kl0 = 0
         ku0 = (nzb+2*nguard0-1)*k3d
       endif
#endif /* NO_PERMANENT_GUARDCELLS */

! cell centered data
      iv_c = max(1,nvar)
      ion_c = min(nvar,1)
! cell face-centered data
      iv_f = max(1,nfacevar)
      ion_f = min(nfacevar,1)
! cell face-centered data
      iv_e = max(1,nvaredge)
      ion_e = min(nvaredge,1)
! cell corner data
      iv_n = max(1,nvarcorn)
      ion_n = min(nvarcorn,1)

#ifdef CHECKP_ON
#ifdef MPI_USED
      iv_c = max(1,nvar_chk_cc)
      ion_c = min(nvar_chk_cc,1)
      iv_f = max(1,nvar_chk_fc)
      ion_f = min(nvar_chk_fc,1)
      iv_e = max(1,nvar_chk_ec)
      ion_e = min(nvar_chk_ec,1)
      iv_n = max(1,nvar_chk_nc)
      ion_n = min(nvar_chk_nc,1)
#endif /* MPI_USED  */
#endif /* CHECKP_ON */


#ifndef MPI_USED

      if (iproc .eq. 0) then

         open(unit=iunit1,form='unformatted',status='unknown'
     $        )

         write (iunit1) tot_blocks

         do jproc = 0,nprocs-1
            
! fetch lnblocks from other processors

            call SHMEM_INTEGER_GET (lnblockst,lnblocks,1,jproc)
            call SHMEM_INTEGER_GET (lrefinet,lrefine,lnblockst,jproc)
            call SHMEM_INTEGER_GET (nodetypet,nodetype,lnblockst,jproc)
            call SHMEM_INTEGER_GET (which_childt,which_child,
     $                              lnblockst,jproc)
            call SHMEM_INTEGER_GET (bflagst,bflags,lnblockst*mflags,
     $                              jproc)
            call SHMEM_REAL_GET (coordt,coord,mdim*lnblockst,jproc)
            call SHMEM_REAL_GET (bnd_boxt,bnd_box,2*mdim*lnblockst,
     $                           jproc)
            call SHMEM_REAL_GET (work_blockt,work_block,
     $                           lnblockst,jproc)
            call SHMEM_INTEGER_GET (gidt,gid,
     $           lnblockst*(nfaces+1+nchild),jproc)

            do block_no = 1,lnblockst

! fetch data for this block
               if(nvar.gt.0)
     $         call SHMEM_REAL_GET (unkt,unk(1,1,1,1,block_no),
     $              len_block,jproc)
               if(nfacevar.gt.0) then
               call shmem_real_get 
     $              (facevarxt,facevarx(1,1,1,1,block_no),
     $              nbndvar*len_blockfx,jproc)
               call shmem_real_get 
     $              (facevaryt,facevary(1,1,1,1,block_no),
     $              nbndvar*len_blockfy,jproc)
               call shmem_real_get 
     $              (facevarzt,facevarz(1,1,1,1,block_no),
     $              nbndvar*len_blockfz,jproc)
               endif
               if(nvaredge.gt.0) then
               call shmem_real_get 
     $              (unk_e_xt,unk_e_x(1,1,1,1,block_no),
     $              nbndvare*len_blockex,jproc)
               call shmem_real_get 
     $              (unk_e_yt,unk_e_y(1,1,1,1,block_no),
     $              nbndvare*len_blockey,jproc)
               call shmem_real_get 
     $              (unk_e_zt,unk_e_z(1,1,1,1,block_no),
     $              nbndvare*len_blockez,jproc)
               endif
               if(nvarcorn.gt.0) 
     $          call shmem_real_get 
     $              (unk_nt,unk_n(1,1,1,1,block_no),
     $              nvarcorn*len_blockn,jproc)

               write (iunit1) 
     &              lrefinet(block_no),
     &              nodetypet(block_no),
     &              which_childt(block_no),
     &              (gidt(j,block_no),j=1,nfaces+1+nchild),
     &              (bflagst(j,block_no),j=1,mflags),
     &              (coordt(j,block_no),j=1,ndim),
     &              (bnd_boxt(1,j,block_no),j=1,ndim),
     &              (bnd_boxt(2,j,block_no),j=1,ndim),
     &              work_blockt(block_no)
#ifndef CHECKP_ON
               if(nvar.gt.0) then
               write (iunit1) 
     &              ((((unkt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_c),
     &                 ix = 1+il0*ion_c,1+iu0*ion_c),
     &                 iy = 1+jl0*ion_c,1+ju0*ion_c),
     &                 iz = 1+kl0*ion_c,1+ku0*ion_c)
               endif
               if(nfacevar.gt.0) then
               write (iunit1) 
     &              ((((facevarxt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_f),
     &                 ix = 1+il0*ion_f,1+(iu0+1)*ion_f),
     &                 iy = 1+jl0*ion_f,1+ju0*ion_f),
     &                 iz = 1+kl0*ion_f,1+ku0*ion_f)
               write (iunit1) 
     &              ((((facevaryt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_f),
     &                 ix = 1+il0*ion_f,1+iu0*ion_f),
     &                 iy = 1+jl0*ion_f,1+(ju0+k2d)*ion_f),
     &                 iz = 1+kl0*ion_f,1+ku0*ion_f)
               write (iunit1) 
     &              ((((facevarzt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_f),
     &                 ix = 1+il0*ion_f,1+iu0*ion_f),
     &                 iy = 1+jl0*ion_f,1+ju0*ion_f),
     &                 iz = 1+kl0*ion_f,1+(ku0+k3d)*ion_f)
               endif
               if(nvaredge.gt.0) then
               write (iunit1) 
     &              ((((unk_e_xt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_e),
     &                 ix = 1+il0*ion_e,1+iu0*ion_e),
     &                 iy = 1+jl0*ion_e,1+(ju0+k2d)*ion_e),
     &                 iz = 1+kl0*ion_e,1+(ku0+k3d)*ion_e)
               write (iunit1) 
     &              ((((unk_e_yt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_e),
     &                 ix = 1+il0*ion_e,1+(iu0+1)*ion_e),
     &                 iy = 1+jl0*ion_e,1+ju0*ion_e),
     &                 iz = 1+kl0*ion_e,1+(ku0+k3d)*ion_e)
               write (iunit1) 
     &              ((((unk_e_zt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_e),
     &                 ix = 1+il0*ion_e,1+(iu0+1)*ion_e),
     &                 iy = 1+jl0*ion_e,1+(ju0+k2d)*ion_e),
     &                 iz = 1+kl0*ion_e,1+ku0*ion_e)
               endif
               if(nvarcorn.gt.0) then
               write (iunit1) 
     &              ((((unk_nt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_n),
     &                 ix = 1+il0*ion_n,1+(iu0+1)*ion_n),
     &                 iy = 1+jl0*ion_n,1+(ju0+k2d)*ion_n),
     &                 iz = 1+kl0*ion_n,1+(ku0+k3d)*ion_n)
               endif
#else /*  CHECKP_ON */
               if(nvar_chk_cc.gt.0) then
                 do ivar=1,nvar
                 if(checkp_on_cc(ivar)) then
                 do iz=1+ik0*ion_c,1+ku0*ion_c
                 do iy=1+ij0*ion_c,1+ju0*ion_c
                 do ix=1+il0*ion_c,1+iu0*ion_c
                   write (iunit1) unkt(ivar,ix,iy,iz)
                 enddo
                 enddo
                 enddo
                 endif
                 enddo
               endif

               if(nvar_chk_fc.gt.0) then
                 do ivar=1,nfacevar
                 if(checkp_on_fc(1,ivar)) then
                 do iz = 1+kl0*ion_f,1+ku0*ion_f
                 do iy = 1+jl0*ion_f,1+ju0*ion_f
                 do ix = 1+il0*ion_f,1+(iu0+1)*ion_f
                   write (iunit1) facevarxt(ivar,ix,iy,iz)
                 enddo
                 enddo
                 enddo
                 endif
                 enddo

                 do ivar=1,nfacevar
                 if(checkp_on_fc(2,ivar)) then
                 do iz = 1+kl0*ion_f,1+ku0*ion_f
                 do iy = 1+jl0*ion_f,1+(ju0+k2d)*ion_f
                 do ix = 1+il0*ion_f,1+iu0*ion_f
                   write (iunit1) facevaryt(ivar,ix,iy,iz)
                 enddo
                 enddo
                 enddo
                 endif
                 enddo

                 do ivar=1,nfacevar
                 if(checkp_on_fc(3,ivar)) then
                 do iz = 1+kl0*ion_f,1+(ku0+k3d)*ion_f
                 do iy = 1+jl0*ion_f,1+ju0*ion_f
                 do ix = 1+il0*ion_f,1+iu0*ion_f
                   write (iunit1) facevarzt(ivar,ix,iy,iz)
                 enddo
                 enddo
                 enddo
                 endif
                 enddo

               endif

               if(nvar_chk_ec.gt.0) then
                 do ivar=1,nvaredge
                 if(checkp_on_ec(1,ivar)) then
                 do iz = 1+kl0*ion_e,1+(ku0+k3d)*ion_e
                 do iy = 1+jl0*ion_e,1+(ju0+k2d)*ion_e
                 do ix = 1+il0*ion_e,1+iu0*ion_e
                   write (iunit1) unk_e_xt(ivar,ix,iy,iz)
                 enddo
                 enddo
                 enddo
                 endif
                 enddo

                 do ivar=1,nvaredge
                 if(checkp_on_ec(2,ivar)) then
                 do iz = 1+kl0*ion_e,1+(ku0+k3d)*ion_e
                 do iy = 1+jl0*ion_e,1+ju0*ion_e
                 do ix = 1+il0*ion_e,1+(iu0+1)*ion_e
                   write (iunit1) unk_e_yt(ivar,ix,iy,iz)
                 enddo
                 enddo
                 enddo
                 endif
                 enddo

                 do ivar=1,nvaredge
                 if(checkp_on_ec(3,ivar)) then
                 do iz = 1+kl0*ion_e,1+ku0*ion_e
                 do iy = 1+jl0*ion_e,1+(ju0+k2d)*ion_e
                 do ix = 1+il0*ion_e,1+(iu0+1)*ion_e
                   write (iunit1) unk_e_zt(ivar,ix,iy,iz)
                 enddo
                 enddo
                 enddo
                 endif
                 enddo

               endif

               if(nvar_chk_nc.gt.0) then
                 do ivar=1,nvarcorn
                 if(checkp_on_nc(ivar)) then
                 do iz = 1+kl0*ion_n,1+(ku0+k3d)*ion_n
                 do iy = 1+jl0*ion_n,1+(ju0+k2d)*ion_n
                 do ix = 1+il0*ion_n,1+(iu0+1)*ion_n
                   write (iunit1) unk_nt(ivar,ix,iy,iz)
                 enddo
                 enddo
                 enddo
                 endif
                 enddo
               endif

#endif /* CHECKP_ON */
            end do
         enddo

         close(iunit1)

      end if ! if iproc == 0
      


#else /* MPI_USED */



#ifdef DEBUG
      write(*,*) 'pe ',iproc,' starting packing '
#endif /* DEBUG */

      CR_buffer1 = 0
      CS_buffer1 = 0
      CR_buffer2 = 0
      CS_buffer2 = 0

#ifdef DEBUG
      write(*,*) 'pe ',iproc,' initialized buffers '
#endif /* DEBUG */

      if (iproc .eq. 0) then

         open(unit=iunit1,form='unformatted',status='unknown'
     $        )

         write (iunit1) tot_blocks

#ifdef DEBUG
      write(*,*) 'pe ',iproc,' writing data from proc 0 '
#endif /* DEBUG */


! Write data from processor 0
         do block_no = 1,lnblocks

! fetch data for this block
               write (iunit1) 
     &              lrefine(block_no),
     &              nodetype(block_no),
     &              which_child(block_no),
     &              (gid(j,block_no),j=1,nfaces+1+nchild),
     &              (bflags(j,block_no),j=1,mflags),
     &              (coord(j,block_no),j=1,ndim),
     &              (bnd_box(1,j,block_no),j=1,ndim),
     &              (bnd_box(2,j,block_no),j=1,ndim),
     &              work_block(block_no)
#ifndef CHECKP_ON
               if(nvar.gt.0) then
               write (iunit1)
     &              ((((unk(ivar,ix,iy,iz,block_no),
     &                 ivar = 1,iv_c),
     &                 ix = 1+il0*ion_c,1+iu0*ion_c),
     &                 iy = 1+jl0*ion_c,1+ju0*ion_c),
     &                 iz = 1+kl0*ion_c,1+ku0*ion_c)
               endif
               if(nfacevar.gt.0) then
               write (iunit1)
     &              ((((facevarx(ivar,ix,iy,iz,block_no),
     &                 ivar = 1,iv_f),
     &                 ix = 1+il0*ion_f,1+(iu0+1)*ion_f),
     &                 iy = 1+jl0*ion_f,1+ju0*ion_f),
     &                 iz = 1+kl0*ion_f,1+ku0*ion_f)
               write (iunit1)
     &              ((((facevary(ivar,ix,iy,iz,block_no),
     &                 ivar = 1,iv_f),
     &                 ix = 1+il0*ion_f,1+iu0*ion_f),
     &                 iy = 1+jl0*ion_f,1+(ju0+k2d)*ion_f),
     &                 iz = 1+kl0*ion_f,1+ku0*ion_f)
               write (iunit1)
     &              ((((facevarz(ivar,ix,iy,iz,block_no),
     &                 ivar = 1,iv_f),
     &                 ix = 1+il0*ion_f,1+iu0*ion_f),
     &                 iy = 1+jl0*ion_f,1+ju0*ion_f),
     &                 iz = 1+kl0*ion_f,1+(ku0+k3d)*ion_f)
               endif
               if(nvaredge.gt.0) then
               write (iunit1)
     &              ((((unk_e_x(ivar,ix,iy,iz,block_no),
     &                 ivar = 1,iv_e),
     &                 ix = 1+il0*ion_e,1+iu0*ion_e),
     &                 iy = 1+jl0*ion_e,1+(ju0+k2d)*ion_e),
     &                 iz = 1+kl0*ion_e,1+(ku0+k3d)*ion_e)
               write (iunit1)
     &              ((((unk_e_y(ivar,ix,iy,iz,block_no),
     &                 ivar = 1,iv_e),
     &                 ix = 1+il0*ion_e,1+(iu0+1)*ion_e),
     &                 iy = 1+jl0*ion_e,1+ju0*ion_e),
     &                 iz = 1+kl0*ion_e,1+(ku0+k3d)*ion_e)
               write (iunit1)
     &              ((((unk_e_z(ivar,ix,iy,iz,block_no),
     &                 ivar = 1,iv_e),
     &                 ix = 1+il0*ion_e,1+(iu0+1)*ion_e),
     &                 iy = 1+jl0*ion_e,1+(ju0+k2d)*ion_e),
     &                 iz = 1+kl0*ion_e,1+ku0*ion_e)
               endif
               if(nvarcorn.gt.0) then
               write (iunit1)
     &              ((((unk_n(ivar,ix,iy,iz,block_no),
     &                 ivar = 1,iv_n),
     &                 ix = 1+il0*ion_n,1+(iu0+1)*ion_n),
     &                 iy = 1+jl0*ion_n,1+(ju0+k2d)*ion_n),
     &                 iz = 1+kl0*ion_n,1+(ku0+k3d)*ion_n)
               endif
#else /*  CHECKP_ON */
               if(nvar_chk_cc.gt.0) then
                 do ivar=1,nvar
                 if(checkp_on_cc(ivar)) then
                 do iz=1+ik0*ion_c,1+ku0*ion_c
                 do iy=1+ij0*ion_c,1+ju0*ion_c
                 do ix=1+il0*ion_c,1+iu0*ion_c
                   write (iunit1) unk(ivar,ix,iy,iz,block_no)
                 enddo
                 enddo
                 enddo
                 endif
                 enddo
               endif

               if(nvar_chk_fc.gt.0) then
                 do ivar=1,nfacevar
                 if(checkp_on_fc(1,ivar)) then
                 do iz = 1+kl0*ion_f,1+ku0*ion_f
                 do iy = 1+jl0*ion_f,1+ju0*ion_f
                 do ix = 1+il0*ion_f,1+(iu0+1)*ion_f
                   write (iunit1) facevarx(ivar,ix,iy,iz,block_no)
                 enddo
                 enddo
                 enddo
                 endif
                 enddo

                 do ivar=1,nfacevar
                 if(checkp_on_fc(2,ivar)) then
                 do iz = 1+kl0*ion_f,1+ku0*ion_f
                 do iy = 1+jl0*ion_f,1+(ju0+k2d)*ion_f
                 do ix = 1+il0*ion_f,1+iu0*ion_f
                   write (iunit1) facevary(ivar,ix,iy,iz,block_no)
                 enddo
                 enddo
                 enddo
                 endif
                 enddo

                 do ivar=1,nfacevar
                 if(checkp_on_fc(3,ivar)) then
                 do iz = 1+kl0*ion_f,1+(ku0+k3d)*ion_f
                 do iy = 1+jl0*ion_f,1+ju0*ion_f
                 do ix = 1+il0*ion_f,1+iu0*ion_f
                   write (iunit1) facevarz(ivar,ix,iy,iz,block_no)
                 enddo
                 enddo
                 enddo
                 endif
                 enddo

               endif

               if(nvar_chk_ec.gt.0) then
                 do ivar=1,nvaredge
                 if(checkp_on_ec(1,ivar)) then
                 do iz = 1+kl0*ion_e,1+(ku0+k3d)*ion_e
                 do iy = 1+jl0*ion_e,1+(ju0+k2d)*ion_e
                 do ix = 1+il0*ion_e,1+iu0*ion_e
                   write (iunit1) unk_e_x(ivar,ix,iy,iz,block_no)
                 enddo
                 enddo
                 enddo
                 endif
                 enddo

                 do ivar=1,nvaredge
                 if(checkp_on_ec(2,ivar)) then
                 do iz = 1+kl0*ion_e,1+(ku0+k3d)*ion_e
                 do iy = 1+jl0*ion_e,1+ju0*ion_e
                 do ix = 1+il0*ion_e,1+(iu0+1)*ion_e
                   write (iunit1) unk_e_y(ivar,ix,iy,iz,block_no)
                 enddo
                 enddo
                 enddo
                 endif
                 enddo

                 do ivar=1,nvaredge
                 if(checkp_on_ec(3,ivar)) then
                 do iz = 1+kl0*ion_e,1+ku0*ion_e
                 do iy = 1+jl0*ion_e,1+(ju0+k2d)*ion_e
                 do ix = 1+il0*ion_e,1+(iu0+1)*ion_e
                   write (iunit1) unk_e_z(ivar,ix,iy,iz,block_no)
                 enddo
                 enddo
                 enddo
                 endif
                 enddo

               endif

               if(nvar_chk_nc.gt.0) then
                 do ivar=1,nvarcorn
                 if(checkp_on_nc(ivar)) then
                 do iz = 1+kl0*ion_n,1+(ku0+k3d)*ion_n
                 do iy = 1+jl0*ion_n,1+(ju0+k2d)*ion_n
                 do ix = 1+il0*ion_n,1+(iu0+1)*ion_n
                   write (iunit1) unk_n(ivar,ix,iy,iz,block_no)
                 enddo
                 enddo
                 enddo
                 endif
                 enddo
               endif

#endif /* CHECKP_ON */
         enddo


         do jproc = 1,nprocs-1

! Post receives on pe 0 for messages from all other procs

! fetch lnblocks from other processors

            lnblockst = glnblocks(jproc)
            isrc = jproc
            idest= 0
            itag = (isrc+1)*(maxblocks+1)
            isize = lnblockst*( 3+mflags+3*mdim+(nfaces+1+nchild)
     .                + 1) 

#ifdef DEBUG
      write(*,*) 'pe ',iproc,' posting recv from proc ',jproc,
     .      ' tag ',itag
#endif /* DEBUG */

            call Mpi_recv(CR_buffer1,isize,MPI_REAL,
     .       isrc,itag,MPI_COMM_WORLD,status,ierr)

            position = 0
            call Mpi_int_unpack(CR_buffer1,buf_dim_bytes1,position,
     &          lrefinet,lnblockst,MPI_INTEGER,MPI_COMM_WORLD,ierr)
            call Mpi_int_unpack(CR_buffer1,buf_dim_bytes1,position,
     &          nodetypet,lnblockst,MPI_INTEGER,MPI_COMM_WORLD,ierr)
            call Mpi_int_unpack(CR_buffer1,buf_dim_bytes1,position,
     &          which_childt,lnblockst,MPI_INTEGER,MPI_COMM_WORLD,ierr)
            call Mpi_int_unpack(CR_buffer1,buf_dim_bytes1,position,
     &          bflagst,lnblockst*mflags,MPI_INTEGER,
     &          MPI_COMM_WORLD,ierr)
            call Mpi_real_unpack(CR_buffer1,buf_dim_bytes1,position,
     &          coordt(1,1),mdim*lnblockst,MPI_REAL,
     &          MPI_COMM_WORLD,ierr)
            call Mpi_real_unpack(CR_buffer1,buf_dim_bytes1,position,
     &          bnd_boxt(1,1,1),2*mdim*lnblockst,
     &          MPI_REAL,
     &          MPI_COMM_WORLD,ierr)
            call Mpi_real_unpack(CR_buffer1,buf_dim_bytes1,position,
     &          work_blockt(1),lnblockst,
     &          MPI_REAL,
     &          MPI_COMM_WORLD,ierr)
            call Mpi_int_unpack(CR_buffer1,buf_dim_bytes1,position,
     &          gidt(1,1),lnblockst*(nfaces+1+nchild),MPI_INTEGER,
     &          MPI_COMM_WORLD,ierr)

            do block_no = 1,lnblockst

            position = 0
            isize = len_block 
     .        + nbndvar*(len_blockfx + len_blockfy + len_blockfz) 
     .        + nbndvare*(len_blockex + len_blockey + len_blockez) 
     .        + nbndvarc*len_blockn 
            itag = block_no

            call Mpi_recv(CR_buffer2,isize,MPI_REAL,
     .       isrc,itag,MPI_COMM_WORLD,status,ierr)

! fetch data for this block
               if(nvar.gt.0) then
                 call Mpi_real_unpack(CR_buffer2,
     &             buf_dim_bytes2,position,
     &             unkt(1,1,1,1),len_block,MPI_REAL,
     &             MPI_COMM_WORLD,ierr)
               endif
               if(nfacevar.gt.0) then
                 call Mpi_real_unpack(CR_buffer2,
     &             buf_dim_bytes2,position,
     &             facevarxt(1,1,1,1),nbndvar*len_blockfx,
     &             MPI_REAL,
     &             MPI_COMM_WORLD,ierr)
                 call Mpi_real_unpack(CR_buffer2,
     &             buf_dim_bytes2,position,
     &             facevaryt(1,1,1,1),nbndvar*len_blockfy,
     &             MPI_REAL,
     &             MPI_COMM_WORLD,ierr)
                 call Mpi_real_unpack(CR_buffer2,
     &             buf_dim_bytes2,position,
     &             facevarzt(1,1,1,1),nbndvar*len_blockfz,
     &             MPI_REAL,
     &             MPI_COMM_WORLD,ierr)
               endif
               if(nvaredge.gt.0) then
                 call Mpi_real_unpack(CR_buffer2,
     &             buf_dim_bytes2,position,
     &             unk_e_xt(1,1,1,1),nbndvare*len_blockex,
     &             MPI_REAL,
     &             MPI_COMM_WORLD,ierr)
                 call Mpi_real_unpack(CR_buffer2,
     &             buf_dim_bytes2,position,
     &             unk_e_yt(1,1,1,1),nbndvare*len_blockey,
     &             MPI_REAL,
     &             MPI_COMM_WORLD,ierr)
                 call Mpi_real_unpack(CR_buffer2,
     &             buf_dim_bytes2,position,
     &             unk_e_zt(1,1,1,1),nbndvare*len_blockez,
     &             MPI_REAL,
     &             MPI_COMM_WORLD,ierr)
               endif
               if(nvarcorn.gt.0) then
                 call Mpi_real_unpack(CR_buffer2,
     &             buf_dim_bytes2,position,
     &             unk_nt(1,1,1,1),nbndvarc*len_blockn,
     &             MPI_REAL,
     &             MPI_COMM_WORLD,ierr)
               endif

               write (iunit1) 
     &              lrefinet(block_no),
     &              nodetypet(block_no),
     &              which_childt(block_no),
     &              (gidt(j,block_no),j=1,nfaces+1+nchild),
     &              (bflagst(j,block_no),j=1,mflags),
     &              (coordt(j,block_no),j=1,ndim),
     &              (bnd_boxt(1,j,block_no),j=1,ndim),
     &              (bnd_boxt(2,j,block_no),j=1,ndim),
     &              work_blockt(block_no)
#ifndef  CHECKP_ON 
               if(nvar.gt.0) then
               write (iunit1) 
     &              ((((unkt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_c),
     &                 ix = 1+il0*ion_c,1+iu0*ion_c),
     &                 iy = 1+jl0*ion_c,1+ju0*ion_c),
     &                 iz = 1+kl0*ion_c,1+ku0*ion_c)
               endif
               if(nfacevar.gt.0) then
               write (iunit1) 
     &              ((((facevarxt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_f),
     &                 ix = 1+il0*ion_f,1+(iu0+1)*ion_f),
     &                 iy = 1+jl0*ion_f,1+ju0*ion_f),
     &                 iz = 1+kl0*ion_f,1+ku0*ion_f)
               write (iunit1) 
     &              ((((facevaryt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_f),
     &                 ix = 1+il0*ion_f,1+iu0*ion_f),
     &                 iy = 1+jl0*ion_f,1+(ju0+k2d)*ion_f),
     &                 iz = 1+kl0*ion_f,1+ku0*ion_f)
               write (iunit1) 
     &              ((((facevarzt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_f),
     &                 ix = 1+il0*ion_f,1+iu0*ion_f),
     &                 iy = 1+jl0*ion_f,1+ju0*ion_f),
     &                 iz = 1+kl0*ion_f,1+(ku0+k3d)*ion_f)
               endif
               if(nvaredge.gt.0) then
               write (iunit1)
     &              ((((unk_e_xt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_e),
     &                 ix = 1+il0*ion_e,1+iu0*ion_e),
     &                 iy = 1+jl0*ion_e,1+(ju0+k2d)*ion_e),
     &                 iz = 1+kl0*ion_e,1+(ku0+k3d)*ion_e)
               write (iunit1) 
     &              ((((unk_e_yt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_e),
     &                 ix = 1+il0*ion_e,1+(iu0+1)*ion_e),
     &                 iy = 1+jl0*ion_e,1+ju0*ion_e),
     &                 iz = 1+kl0*ion_e,1+(ku0+k3d)*ion_e)
               write (iunit1) 
     &              ((((unk_e_zt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_e),
     &                 ix = 1+il0*ion_e,1+(iu0+1)*ion_e),
     &                 iy = 1+jl0*ion_e,1+(ju0+k2d)*ion_e),
     &                 iz = 1+kl0*ion_e,1+ku0*ion_e)
               endif
               if(nvarcorn.gt.0) then
               write (iunit1) 
     &              ((((unk_nt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_n),
     &                 ix = 1+il0*ion_n,1+(iu0+1)*ion_n),
     &                 iy = 1+jl0*ion_n,1+(ju0+k2d)*ion_n),
     &                 iz = 1+kl0*ion_n,1+(ku0+k3d)*ion_n)
               endif
#else /*  CHECKP_ON */
               if(nvar_chk_cc.gt.0) then
                 do ivar=1,nvar
                 if(checkp_on_cc(ivar)) then
                 do iz=1+ik0*ion_c,1+ku0*ion_c
                 do iy=1+ij0*ion_c,1+ju0*ion_c
                 do ix=1+il0*ion_c,1+iu0*ion_c
                   write (iunit1) unkt(ivar,ix,iy,iz)
                 enddo
                 enddo
                 enddo
                 endif
                 enddo
               endif

               if(nvar_chk_fc.gt.0) then
                 do ivar=1,nfacevar
                 if(checkp_on_fc(1,ivar)) then
                 do iz = 1+kl0*ion_f,1+ku0*ion_f
                 do iy = 1+jl0*ion_f,1+ju0*ion_f
                 do ix = 1+il0*ion_f,1+(iu0+1)*ion_f
                   write (iunit1) facevarxt(ivar,ix,iy,iz)
                 enddo
                 enddo
                 enddo
                 endif
                 enddo

                 do ivar=1,nfacevar
                 if(checkp_on_fc(2,ivar)) then
                 do iz = 1+kl0*ion_f,1+ku0*ion_f
                 do iy = 1+jl0*ion_f,1+(ju0+k2d)*ion_f
                 do ix = 1+il0*ion_f,1+iu0*ion_f
                   write (iunit1) facevaryt(ivar,ix,iy,iz)
                 enddo
                 enddo
                 enddo
                 endif
                 enddo

                 do ivar=1,nfacevar
                 if(checkp_on_fc(3,ivar)) then
                 do iz = 1+kl0*ion_f,1+(ku0+k3d)*ion_f
                 do iy = 1+jl0*ion_f,1+ju0*ion_f
                 do ix = 1+il0*ion_f,1+iu0*ion_f
                   write (iunit1) facevarzt(ivar,ix,iy,iz)
                 enddo
                 enddo
                 enddo
                 endif
                 enddo

               endif

               if(nvar_chk_ec.gt.0) then
                 do ivar=1,nvaredge
                 if(checkp_on_ec(1,ivar)) then
                 do iz = 1+kl0*ion_e,1+(ku0+k3d)*ion_e
                 do iy = 1+jl0*ion_e,1+(ju0+k2d)*ion_e
                 do ix = 1+il0*ion_e,1+iu0*ion_e
                   write (iunit1) unk_e_xt(ivar,ix,iy,iz)
                 enddo
                 enddo
                 enddo
                 endif
                 enddo

                 do ivar=1,nvaredge
                 if(checkp_on_ec(2,ivar)) then
                 do iz = 1+kl0*ion_e,1+(ku0+k3d)*ion_e
                 do iy = 1+jl0*ion_e,1+ju0*ion_e
                 do ix = 1+il0*ion_e,1+(iu0+1)*ion_e
                   write (iunit1) unk_e_yt(ivar,ix,iy,iz)
                 enddo
                 enddo
                 enddo
                 endif
                 enddo

                 do ivar=1,nvaredge
                 if(checkp_on_ec(3,ivar)) then
                 do iz = 1+kl0*ion_e,1+ku0*ion_e
                 do iy = 1+jl0*ion_e,1+(ju0+k2d)*ion_e
                 do ix = 1+il0*ion_e,1+(iu0+1)*ion_e
                   write (iunit1) unk_e_zt(ivar,ix,iy,iz)
                 enddo
                 enddo
                 enddo
                 endif
                 enddo

               endif

               if(nvar_chk_nc.gt.0) then
                 do ivar=1,nvarcorn
                 if(checkp_on_nc(ivar)) then
                 do iz = 1+kl0*ion_n,1+(ku0+k3d)*ion_n
                 do iy = 1+jl0*ion_n,1+(ju0+k2d)*ion_n
                 do ix = 1+il0*ion_n,1+(iu0+1)*ion_n
                   write (iunit1) unk_nt(ivar,ix,iy,iz)
                 enddo
                 enddo
                 enddo
                 endif
                 enddo
               endif

#endif /* CHECKP_ON */

            enddo
         enddo
         close(iunit1)

       end if ! if iproc == 0

       if(iproc.gt.0) then

!
! Post sends to pe 0 for messages from all other procs

! fetch lnblocks from other processors

            lnblockst = lnblocks
            isrc = iproc
            idest= 0
            itag = (isrc+1)*(maxblocks+1)
            isize = lnblockst*( 3+mflags+3*mdim+(nfaces+1+nchild)
     .                          + 1 )
#ifdef DEBUG
      write(*,*) 'pe ',iproc,' packing data for send '
#endif /* DEBUG */
            position = 0
            call Mpi_int_pack(lrefine(1),lnblocks,MPI_INTEGER,
     &          CS_buffer1,
     &          buf_dim_bytes1,position,MPI_COMM_WORLD,ierr)
            call Mpi_int_pack(nodetype(1),lnblocks,MPI_INTEGER,
     &          CS_buffer1,
     &          buf_dim_bytes1,position,MPI_COMM_WORLD,ierr)
            call Mpi_int_pack(which_child(1),lnblocks,MPI_INTEGER,
     &          CS_buffer1,
     &          buf_dim_bytes1,position,MPI_COMM_WORLD,ierr)
            call Mpi_int_pack(bflags(1,1),lnblocks*mflags,MPI_INTEGER,
     &          CS_buffer1,buf_dim_bytes1,position,MPI_COMM_WORLD,ierr)
            call Mpi_real_pack(coord(1,1),mdim*lnblocks,
     &          MPI_REAL,
     &          CS_buffer1,buf_dim_bytes1,position,MPI_COMM_WORLD,ierr)
            call Mpi_real_pack(bnd_box(1,1,1),2*mdim*lnblocks,
     &          MPI_REAL,
     &          CS_buffer1,buf_dim_bytes1,position,MPI_COMM_WORLD,ierr)
            call Mpi_real_pack(work_block(1),lnblocks,
     &          MPI_REAL,
     &          CS_buffer1,buf_dim_bytes1,position,MPI_COMM_WORLD,ierr)
            call Mpi_int_pack(gid(1,1),lnblocks*(nfaces+1+nchild),
     &          MPI_INTEGER,CS_buffer1,buf_dim_bytes1,position,
     &          MPI_COMM_WORLD,ierr)

            call Mpi_send(CS_buffer1,isize,MPI_REAL,
     .       idest,itag,MPI_COMM_WORLD,ierr)

#ifdef DEBUG
      write(*,*) 'pe ',iproc,' tree data packed for send '
#endif /* DEBUG */

            do block_no= 1,lnblocks

              position = 0

              if(nvar.gt.0)
     &          call Mpi_real_pack(unk(1,1,1,1,block_no),
     &            len_block,MPI_REAL,CS_buffer2,
     &            buf_dim_bytes2,position,MPI_COMM_WORLD,ierr)

              if(nfacevar.gt.0) then
                call Mpi_real_pack(facevarx(1,1,1,1,block_no),
     &            len_blockfx*nbndvar,
     &            MPI_REAL,CS_buffer2,
     &            buf_dim_bytes2,position,MPI_COMM_WORLD,ierr)
                call Mpi_real_pack(facevary(1,1,1,1,block_no),
     &            len_blockfy*nbndvar,
     &            MPI_REAL,CS_buffer2,
     &            buf_dim_bytes2,position,MPI_COMM_WORLD,ierr)
                call Mpi_real_pack(facevarz(1,1,1,1,block_no),
     &            len_blockfz*nbndvar,
     &            MPI_REAL,CS_buffer2,
     &            buf_dim_bytes2,position,MPI_COMM_WORLD,ierr)
               endif
              if(nvaredge.gt.0) then
                call Mpi_real_pack(unk_e_x(1,1,1,1,block_no),
     &            len_blockex*nbndvare,
     &            MPI_REAL,CS_buffer2,
     &            buf_dim_bytes2,position,MPI_COMM_WORLD,ierr)
                call Mpi_real_pack(unk_e_y(1,1,1,1,block_no),
     &            len_blockey*nbndvare,
     &            MPI_REAL,CS_buffer2,
     &            buf_dim_bytes2,position,MPI_COMM_WORLD,ierr)
                call Mpi_real_pack(unk_e_z(1,1,1,1,block_no),
     &            len_blockez*nbndvare,
     &            MPI_REAL,CS_buffer2,
     &            buf_dim_bytes2,position,MPI_COMM_WORLD,ierr)
               endif
              if(nvarcorn.gt.0) then
                call Mpi_real_pack(unk_n(1,1,1,1,block_no),
     &            len_blockn*nbndvarc,
     &            MPI_REAL,CS_buffer2,
     &            buf_dim_bytes2,position,MPI_COMM_WORLD,ierr)
               endif

               itag = block_no
               isize = len_block 
     .           + nbndvar*(len_blockfx + len_blockfy + len_blockfz) 
     .           + nbndvare*(len_blockex + len_blockey + len_blockez) 
     .           + nbndvarc*len_blockn 

               call Mpi_send(CS_buffer2,isize,MPI_REAL,
     .                       idest,itag,MPI_COMM_WORLD,ierr)

            enddo
 
#ifdef DEBUG
      write(*,*) 'pe ',iproc,' sending data : tag ',itag
#endif /* DEBUG */

#ifdef DEBUG
      write(*,*) 'pe ',iproc,' data sent '
#endif /* DEBUG */
       end if ! if iproc > 0

#endif /* MPI_USED */
      
      if(allocated(n_to_left)) deallocate( n_to_left )
      if(allocated(glnblocks)) deallocate( glnblocks )
      if(allocated(CS_buffer1)) deallocate( CS_buffer1 )
      if(allocated(CR_buffer1)) deallocate( CR_buffer1 )
      if(allocated(CS_buffer2)) deallocate( CS_buffer2 )
      if(allocated(CR_buffer2)) deallocate( CR_buffer2 )

      return
      end subroutine amr_checkpoint_wr

!-----------------------------------------------------------------------
      

      subroutine amr_checkpoint_re (iunit1,l_with_guardcells)

! Subroutine to read checkpoint file using AMR package.
! Writes out tree data structure and data stored in blocks
! Currently reads are done serially by processor 0 and data sent to
! other processors.
! Uses UNFORMATTED, DIRECT ACCESS

! Written: K. Olson 7/97


      use paramesh_dimensions
      use physicaldata
      use tree
      use timings
#ifdef SAVE_MORTS
      use mpi_morton
#endif /*  SAVE_MORTS */

      use paramesh_interfaces, only : amr_morton_order,
     &                                amr_guardcell
      use paramesh_mpi_interfaces, only : mpi_amr_global_domain_limits,
     &                                    mpi_amr_boundary_block_info


      implicit none

#include "amr_shmem.fh"

#ifdef MPI_USED
      include 'mpif.h'
#endif /* MPI_USED */


      integer, intent(in) :: iunit1
      logical, optional, intent(in)  :: l_with_guardcells

      integer nguard0
      parameter(nguard0 = nguard*npgs)

      integer block_no,shmem_n_pes,shmem_my_pe
      integer jproc,i,j,ivar,ix,iy,iz,nprocs,iproc

      integer lnblockst,ngid
      integer alnblocks,alnblockst
      integer ierr_read,ierr_readt
      integer ierr,mype

! TEMPORARIES WHICH ARE READ

      integer  tot_blocks

      integer gid(nfaces+1+nchild,maxblocks_tr)
      integer  lrefinet(maxblocks_tr),nodetypet(maxblocks_tr)
      integer  which_childt(maxblocks_tr)
      integer  gidt(nfaces+1+nchild,maxblocks_tr)
      integer  bflagst(mflags,maxblocks_tr)
      real  coordt(mdim,maxblocks_tr)
      real  work_blockt(maxblocks_tr)
      real  bnd_boxt(2,mdim,maxblocks_tr)

      real  unkt(nvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
      real  facevarxt(nbndvar,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,
     &                kl_bnd:ku_bnd)
      real  facevaryt(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,
     &                kl_bnd:ku_bnd)
      real  facevarzt(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,
     &                kl_bnd:ku_bnd+k3d)
      real  :: unk_nt(nbndvarc,
     &               il_bnd:iu_bnd+1,
     &               jl_bnd:ju_bnd+k2d,
     &               kl_bnd:ku_bnd+k3d)
      real  :: unk_e_xt(nbndvare,
     &                  il_bnd:iu_bnd,
     &                  jl_bnd:ju_bnd+k2d,
     &                  kl_bnd:ku_bnd+k3d)
      real  :: unk_e_yt(nbndvare,
     &                  il_bnd:iu_bnd+1,
     &                  jl_bnd:ju_bnd,
     &                  kl_bnd:ku_bnd+k3d)
      real  :: unk_e_zt(nbndvare,
     &                  il_bnd:iu_bnd+1,
     &                  jl_bnd:ju_bnd+k2d,
     &                  kl_bnd:ku_bnd)

      integer :: il0,iu0,jl0,ju0,kl0,ku0
      integer :: ion_c,ion_f,ion_e,ion_n,iv_c,iv_f,iv_e,iv_n

#ifdef MPI_USED
      integer,dimension (:),  allocatable :: glnblocks
      integer :: isrc,idest,itag,isize,position,position2
      integer :: ierrorcode
      integer :: lnblocks_old
      integer status(MPI_STATUS_SIZE)

      integer :: buf_dim_int 
      integer :: buf_dim_real
      integer :: buf_dim1,buf_dim2
      integer :: no_of_bytes_per_real,no_of_bytes_per_integer
      integer :: buf_dim_bytes1,buf_dim_bytes2

      integer :: nvar_chk_cc,nvar_chk_fc,nvar_chk_ec,nvar_chk_nc

      logical :: l_move_solution

      real,allocatable :: CS_buffer1(:), CR_buffer1(:)
      real,allocatable :: CS_buffer2(:), CR_buffer2(:)
      integer,dimension (1) :: ibcast_data


#ifdef SAVE_MORTS
      integer :: lb
      integer :: mort_neigh(6,3,3,3)
      real    :: xmin,ymin,zmin,xmax,ymax,zmax
#endif /* SAVE_MORTS */
#ifdef TIMING_MPI
      double precision :: time1
#endif /* TIMING_MPI */

      buf_dim_int = 
     .         maxblocks*( 3+mflags+(nfaces+1+nchild))
      buf_dim_real = 
     .         maxblocks*( 3*mdim + 1 )
      buf_dim1 = buf_dim_real + buf_dim_int

      allocate(CS_buffer1(buf_dim1))
      allocate(CR_buffer1(buf_dim1))
#endif /* MPI_USED */

#ifdef CHECKP_ON
      nvar_chk_cc =  0
      do i=1,nvar
        if(checkp_on_cc(i)) nvar_chk_cc = nvar_chk_cc + 1
      enddo
      nvar_chk_fc =  0
      do i=1,nfacevar
        if(checkp_on_fc(1,i)) nvar_chk_fc = nvar_chk_fc + 1
      enddo
      nvar_chk_ec =  0
      do i=1,nvaredge
        if(checkp_on_ec(1,i)) nvar_chk_ec = nvar_chk_ec + 1
      enddo
      nvar_chk_nc =  0
      do i=1,nvarcorn
        if(checkp_on_nc(i)) nvar_chk_nc = nvar_chk_nc + 1
      enddo
#endif /* CHECKP_ON */

#ifdef MPI_USED
      buf_dim2 = len_block
     .       + nbndvar*(len_blockfx + len_blockfy + len_blockfz)
     .       + nbndvarc*len_blockn
     .       + nbndvare*(len_blockex + len_blockey + len_blockez)

      allocate(CS_buffer2(buf_dim2))
      allocate(CR_buffer2(buf_dim2))

#endif /* MPI_USED */

      nprocs = shmem_n_pes()
      iproc  = shmem_my_pe()
      mype = iproc
      ierr_read = 0
      ierr_readt = 0


#ifdef MPI_USED
      if(allocated(glnblocks)) deallocate( glnblocks )
      allocate ( glnblocks(0:nprocs-1) )

      call mpi_type_size(MPI_INTEGER,no_of_bytes_per_integer,ierr)
      call mpi_type_size(MPI_REAL
     .                  ,no_of_bytes_per_real   ,ierr)
      buf_dim_bytes1 = buf_dim_real*no_of_bytes_per_real +
     .                 buf_dim_int *no_of_bytes_per_integer
      buf_dim_bytes2 = buf_dim2*no_of_bytes_per_real 

#endif /* MPI_USED */

! set limits on data arrays
       il0 = nguard0
       iu0 = nxb-1+nguard0
       jl0 = nguard0*k2d
       ju0 = (nyb-1+nguard0)*k2d
       kl0 = nguard0*k3d
       ku0 = (nzb-1+nguard0)*k3d
#ifndef NO_PERMANENT_GUARDCELLS
       if(l_with_guardcells) then
         il0 = 0
         iu0 = nxb+2*nguard0-1
         jl0 = 0
         ju0 = (nyb+2*nguard0-1)*k2d
         kl0 = 0
         ku0 = (nzb+2*nguard0-1)*k3d
       endif
#endif /* NO_PERMANENT_GUARDCELLS */


! cell centered data
      iv_c = max(1,nvar)
      ion_c = min(nvar,1)
! cell face-centered data
      iv_f = max(1,nfacevar)
      ion_f = min(nfacevar,1)
! cell face-centered data
      iv_e = max(1,nvaredge)
      ion_e = min(nvaredge,1)
! cell corner data
      iv_n = max(1,nvarcorn)
      ion_n = min(nvarcorn,1)

#ifdef CHECKP_ON
#ifdef MPI_USED
      iv_c = max(1,nvar_chk_cc)
      ion_c = min(nvar_chk_cc,1)
      iv_f = max(1,nvar_chk_fc)
      ion_f = min(nvar_chk_fc,1)
      iv_e = max(1,nvar_chk_ec)
      ion_e = min(nvar_chk_ec,1)
      iv_n = max(1,nvar_chk_nc)
      ion_n = min(nvar_chk_nc,1)
#endif /* MPI_USED  */
#endif /* CHECKP_ON */

#ifndef MPI_USED

      if (iproc .eq. 0) then

         open(unit=iunit1,form='unformatted',status='unknown'
     $        )

         read (iunit1) tot_blocks
         write(*,*) 'blocks to be input ',tot_blocks

! compute approximate lnblocks (this will be the number of blocks stored on
! processors 0 -> nprocs-2, nprocs-1 gets tot_blocks - the total number on the
! rest of the blocks)

         alnblocks = int(tot_blocks/nprocs)

! check for error
         if (tot_blocks-(alnblocks*(nprocs-1)).gt.maxblocks) then

          print *,' ******** ERROR in checkpoint_re: ********'
          print *,' No. of blocks per processor exceeds maxblocks.'
          print *,' Suggest you reset maxblocks to a larger number or '
          print *,' run on a larger no. of processors. '

            ierr_read = 1

            go to 2

         end if

         do jproc = 0,nprocs-1

            if (jproc.lt.nprocs-1) then
               lnblockst = alnblocks
            else
               lnblockst = tot_blocks - (alnblocks*(nprocs-1))
            end if

            do block_no = 1,lnblockst

! Read in data for this block

               read (iunit1) 
     &              lrefinet(block_no),
     &              nodetypet(block_no),
     &              which_childt(block_no),
     &              (gidt(j,block_no),j=1,nfaces+1+nchild),
     &              (bflagst(j,block_no),j=1,mflags),
     &              (coordt(j,block_no),j=1,ndim),
     &              (bnd_boxt(1,j,block_no),j=1,ndim),
     &              (bnd_boxt(2,j,block_no),j=1,ndim),
     &              work_blockt(block_no)
#ifndef CHECKP_ON
               if(nvar.gt.0) then
               read (iunit1) 
     &              ((((unkt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_c),
     &                 ix = 1+il0*ion_c,1+iu0*ion_c),
     &                 iy = 1+jl0*ion_c,1+ju0*ion_c),
     &                 iz = 1+kl0*ion_c,1+ku0*ion_c)
               endif
               if(nfacevar.gt.0) then
               read (iunit1) 
     &              ((((facevarxt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_f),
     &                 ix = 1+il0*ion_f,1+(iu0+1)*ion_f),
     &                 iy = 1+jl0*ion_f,1+ju0*ion_f),
     &                 iz = 1+kl0*ion_f,1+ku0*ion_f) 
               read (iunit1) 
     &              ((((facevaryt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_f),
     &                 ix = 1+il0*ion_f,1+iu0*ion_f),
     &                 iy = 1+jl0*ion_f,1+(ju0+k2d)*ion_f),
     &                 iz = 1+kl0*ion_f,1+ku0*ion_f)   
               read (iunit1) 
     &              ((((facevarzt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_f),
     &                 ix = 1+il0*ion_f,1+iu0*ion_f),
     &                 iy = 1+jl0*ion_f,1+ju0*ion_f),
     &                 iz = 1+kl0*ion_f,1+(ku0+k3d)*ion_f),
               endif
               if(nvaredge.gt.0) then
               read (iunit1) 
     &              ((((unk_e_xt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_e),
     &                 ix = 1+il0*ion_e,1+iu0*ion_e),
     &                 iy = 1+jl0*ion_e,1+(ju0+k2d)*ion_e),
     &                 iz = 1+kl0*ion_e,1+(ku0+k3d)*ion_e)
               read (iunit1) 
     &              ((((unk_e_yt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_e),
     &                 ix = 1+il0*ion_e,1+(iu0+1)*ion_e),
     &                 iy = 1+jl0*ion_e,1+ju0*ion_e),
     &                 iz = 1+kl0*ion_e,1+(ku0+k3d)*ion_e)
               read (iunit1) 
     &              ((((unk_e_zt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_e),
     &                 ix = 1+il0*ion_e,1+(iu0+1)*ion_e),
     &                 iy = 1+jl0*ion_e,1+(ju0+k2d)*ion_e),
     &                 iz = 1+kl0*ion_e,1+ku0*ion_e)
               endif
               if(nvarcorn.gt.0) then
               read (iunit1) 
     &              ((((unk_nt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_n),
     &                 ix = 1+il0*ion_n,1+(iu0+1)*ion_n),
     &                 iy = 1+jl0*ion_n,1+(ju0+k2d)*ion_n),
     &                 iz = 1+kl0*ion_n,1+(ku0+k3d)*ion_n) 
               endif
#else /*  CHECKP_ON */
               if(nvar_chk_cc.gt.0) then
                 do ivar=1,nvar
                 if(checkp_on_cc(ivar)) then
                 do iz=1+ik0*ion_c,1+ku0*ion_c
                 do iy=1+ij0*ion_c,1+ju0*ion_c
                 do ix=1+il0*ion_c,1+iu0*ion_c
                   read (iunit1) unkt(ivar,ix,iy,iz)
                 enddo
                 enddo
                 enddo
                 else
                   unkt(ivar,:,:,:) = 0.
                 endif
                 enddo
               endif

               if(nvar_chk_fc.gt.0) then
                 do ivar=1,nfacevar
                 if(checkp_on_fc(1,ivar)) then
                 do iz = 1+kl0*ion_f,1+ku0*ion_f
                 do iy = 1+jl0*ion_f,1+ju0*ion_f
                 do ix = 1+il0*ion_f,1+(iu0+1)*ion_f
                   read (iunit1) facevarxt(ivar,ix,iy,iz)
                 enddo
                 enddo
                 enddo
                 else
                   facevarxt(ivar,:,:,:) = 0.
                 endif
                 enddo

                 do ivar=1,nfacevar
                 if(checkp_on_fc(2,ivar)) then
                 do iz = 1+kl0*ion_f,1+ku0*ion_f
                 do iy = 1+jl0*ion_f,1+(ju0+k2d)*ion_f
                 do ix = 1+il0*ion_f,1+iu0*ion_f
                   read (iunit1) facevaryt(ivar,ix,iy,iz)
                 enddo
                 enddo
                 enddo
                 else
                   facevaryt(ivar,:,:,:) = 0.
                 endif
                 enddo

                 do ivar=1,nfacevar
                 if(checkp_on_fc(3,ivar)) then
                 do iz = 1+kl0*ion_f,1+(ku0+k3d)*ion_f
                 do iy = 1+jl0*ion_f,1+ju0*ion_f
                 do ix = 1+il0*ion_f,1+iu0*ion_f
                   read (iunit1) facevarzt(ivar,ix,iy,iz)
                 enddo
                 enddo
                 enddo
                 else
                   facevarzt(ivar,:,:,:) = 0.
                 endif
                 enddo

               endif

               if(nvar_chk_ec.gt.0) then
                 do ivar=1,nvaredge
                 if(checkp_on_ec(1,ivar)) then
                 do iz = 1+kl0*ion_e,1+(ku0+k3d)*ion_e
                 do iy = 1+jl0*ion_e,1+(ju0+k2d)*ion_e
                 do ix = 1+il0*ion_e,1+iu0*ion_e
                   read (iunit1) unk_e_xt(ivar,ix,iy,iz)
                 enddo
                 enddo
                 enddo
                 else
                   unk_e_xt(ivar,:,:,:) = 0.
                 endif
                 enddo

                 do ivar=1,nvaredge
                 if(checkp_on_ec(2,ivar)) then
                 do iz = 1+kl0*ion_e,1+(ku0+k3d)*ion_e
                 do iy = 1+jl0*ion_e,1+ju0*ion_e
                 do ix = 1+il0*ion_e,1+(iu0+1)*ion_e
                   read (iunit1) unk_e_yt(ivar,ix,iy,iz)
                 enddo
                 enddo
                 enddo
                 else
                   unk_e_yt(ivar,:,:,:) = 0.
                 endif
                 enddo

                 do ivar=1,nvaredge
                 if(checkp_on_ec(3,ivar)) then
                 do iz = 1+kl0*ion_e,1+ku0*ion_e
                 do iy = 1+jl0*ion_e,1+(ju0+k2d)*ion_e
                 do ix = 1+il0*ion_e,1+(iu0+1)*ion_e
                   read (iunit1) unk_e_zt(ivar,ix,iy,iz)
                 enddo
                 enddo
                 enddo
                 else
                   unk_e_zt(ivar,:,:,:) = 0.
                 endif
                 enddo

               endif

               if(nvar_chk_nc.gt.0) then
                 do ivar=1,nvarcorn
                 if(checkp_on_nc(ivar)) then
                 do iz = 1+kl0*ion_n,1+(ku0+k3d)*ion_n
                 do iy = 1+jl0*ion_n,1+(ju0+k2d)*ion_n
                 do ix = 1+il0*ion_n,1+(iu0+1)*ion_n
                   read (iunit1) unk_nt(ivar,ix,iy,iz)
                 enddo
                 enddo
                 enddo
                 else
                   unk_nt(ivar,:,:,:) = 0.
                 endif
                 enddo
               endif

#endif /* CHECKP_ON */

      if(ndim.lt.3) then
        bnd_box(1,3,:) = 0.
        bnd_box(2,3,:) = 1.
        coord(3,:) = .5*(bnd_box(2,3,:)+bnd_box(1,3,:))
        bsize(3,:) =     bnd_box(2,3,:)-bnd_box(1,3,:)
      endif
      if(ndim.lt.2) then
        bnd_box(1,2,:) = 0.
        bnd_box(2,2,:) = 1.
        coord(2,:) = .5*(bnd_box(2,2,:)+bnd_box(1,2,:))
        bsize(2,:) =     bnd_box(2,2,:)-bnd_box(1,2,:)
      endif

! put data on proc jproc for this block
               if(nvar.gt.0)
     $         call SHMEM_REAL_PUT (unk(1,1,1,1,block_no),unkt,
     $              len_block,jproc)
               if(nfacevar.gt.0) then
               call shmem_real_put 
     $              (facevarx(1,1,1,1,block_no),facevarxt,
     $              nbndvar*len_blockfx,jproc)
               call shmem_real_put 
     $              (facevary(1,1,1,1,block_no),facevaryt,
     $              nbndvar*len_blockfy,jproc)
               call shmem_real_put 
     $              (facevarz(1,1,1,1,block_no),facevarzt,
     $              nbndvar*len_blockfz,jproc)
               endif
               if(nvaredge.gt.0) then
               call shmem_real_put 
     $              (unk_e_x(1,1,1,1,block_no),unk_e_xt,
     $              nbndvare*len_blockex,jproc)
               call shmem_real_put 
     $              (unk_e_y(1,1,1,1,block_no),unk_e_yt,
     $              nbndvare*len_blockey,jproc)
               call shmem_real_put 
     $              (unk_e_z(1,1,1,1,block_no),unk_e_zt,
     $              nbndvare*len_blockez,jproc)
               endif
               if(nvarcorn.gt.0) 
     $          call shmem_real_put 
     $              (unk_n(1,1,1,1,block_no),unk_nt,
     $              nbndvarc*len_blockn,jproc)

            end do

! put no. of blocks on processor jproc

            call SHMEM_INTEGER_PUT (lnblocks,lnblockst,1,jproc)

! put tree data on processor jproc

            call SHMEM_REAL_PUT (coord,coordt,mdim*lnblockst,jproc)
            call SHMEM_REAL_PUT (bnd_box,bnd_boxt,2*mdim*lnblockst,
     $                           jproc)
            call SHMEM_REAL_PUT (work_block,work_blockt,lnblockst,
     $                           jproc)
            call SHMEM_INTEGER_PUT (lrefine,lrefinet,lnblockst,jproc)
            call SHMEM_INTEGER_PUT (nodetype,nodetypet,lnblockst,jproc)
            call SHMEM_INTEGER_PUT (which_child,which_childt,
     $                              lnblockst,jproc)
            call SHMEM_INTEGER_PUT (bflags,bflagst,lnblockst*mflags,
     $                              jproc)
            call SHMEM_INTEGER_PUT (gid,gidt,
     $           lnblockst*(nfaces+1+nchild),jproc)

         end do

         close(iunit1)

      end if ! if iproc == 0

 2    continue


#else /* MPI_USED */



      if (iproc .eq. 0) then

         open(unit=iunit1,form='unformatted',status='unknown'
     $        )

         read (iunit1) tot_blocks
         write(*,*) 'blocks to be input ',tot_blocks

! compute approximate lnblocks (this will be the number of blocks stored on
! processors 0 -> nprocs-2, nprocs-1 gets tot_blocks - the total number on the
! rest of the blocks)

         alnblocks = int(tot_blocks/nprocs)

! check for error
         if (tot_blocks-(alnblocks*(nprocs-1)).gt.maxblocks) then

          print *,' ******** ERROR in checkpoint_re: ********'
          print *,' No. of blocks per processor exceeds maxblocks.'
          print *,' Suggest you reset maxblocks to a larger number or '
          print *,' run on a larger no. of processors. '

            ierr_read = 1

            go to 2

         end if

         do jproc = 0,nprocs-1

            if (jproc.lt.nprocs-1) then
               lnblockst = alnblocks
            else
               lnblockst = tot_blocks - (alnblocks*(nprocs-1))
            end if
            glnblocks(jproc) = lnblockst

         enddo

      endif
      call mpi_bcast(glnblocks,nprocs,MPI_INTEGER,0,
     &                 MPI_COMM_WORLD,ierr)
      lnblocks = glnblocks(iproc)

#ifdef DEBUG
      write(*,*) 'pe ',iproc,' glnblocks ',glnblocks
      write(*,*) 'putting ',lnblocks,' on pe ',iproc
#endif /* DEBUG */

      if (iproc .eq. 0) then
            do block_no = 1,lnblocks

! Read in data for this block

               read (iunit1) 
     &              lrefine(block_no),
     &              nodetype(block_no),
     &              which_child(block_no),
     &              (gid(j,block_no),j=1,nfaces+1+nchild),
     &              (bflags(j,block_no),j=1,mflags),
     &              (coord(j,block_no),j=1,ndim),
     &              (bnd_box(1,j,block_no),j=1,ndim),
     &              (bnd_box(2,j,block_no),j=1,ndim),
     &              work_block(block_no)
#ifndef CHECKP_ON
               if(nvar.gt.0) then
               read (iunit1) 
     &              ((((unk(ivar,ix,iy,iz,block_no),
     &                 ivar = 1,iv_c),
     &                 ix = 1+il0*ion_c,1+iu0*ion_c),
     &                 iy = 1+jl0*ion_c,1+ju0*ion_c),
     &                 iz = 1+kl0*ion_c,1+ku0*ion_c)
               endif
               if(nfacevar.gt.0) then
               read (iunit1) 
     &              ((((facevarx(ivar,ix,iy,iz,block_no),
     &                 ivar = 1,iv_f),
     &                 ix = 1+il0*ion_f,1+(iu0+1)*ion_f),
     &                 iy = 1+jl0*ion_f,1+ju0*ion_f),
     &                 iz = 1+kl0*ion_f,1+ku0*ion_f)    
               read (iunit1) 
     &              ((((facevary(ivar,ix,iy,iz,block_no),
     &                 ivar = 1,iv_f),
     &                 ix = 1+il0*ion_f,1+iu0*ion_f),
     &                 iy = 1+jl0*ion_f,1+(ju0+k2d)*ion_f),
     &                 iz = 1+kl0*ion_f,1+ku0*ion_f)    
               read (iunit1) 
     &              ((((facevarz(ivar,ix,iy,iz,block_no),
     &                 ivar = 1,iv_f),
     &                 ix = 1+il0*ion_f,1+iu0*ion_f),
     &                 iy = 1+jl0*ion_f,1+ju0*ion_f),
     &                 iz = 1+kl0*ion_f,1+(ku0+k3d)*ion_f)
               endif
               if(nvaredge.gt.0) then
               read (iunit1) 
     &              ((((unk_e_x(ivar,ix,iy,iz,block_no),
     &                 ivar = 1,iv_e),
     &                 ix = 1+il0*ion_e,1+iu0*ion_e),
     &                 iy = 1+jl0*ion_e,1+(ju0+k2d)*ion_e),
     &                 iz = 1+kl0*ion_e,1+(ku0+k3d)*ion_e)
               read (iunit1) 
     &              ((((unk_e_y(ivar,ix,iy,iz,block_no),
     &                 ivar = 1,iv_e),
     &                 ix = 1+il0*ion_e,1+(iu0+1)*ion_e),
     &                 iy = 1+jl0*ion_e,1+ju0*ion_e),
     &                 iz = 1+kl0*ion_e,1+(ku0+k3d)*ion_e)
               read (iunit1) 
     &              ((((unk_e_z(ivar,ix,iy,iz,block_no),
     &                 ivar = 1,iv_e),
     &                 ix = 1+il0*ion_e,1+(iu0+1)*ion_e),
     &                 iy = 1+jl0*ion_e,1+(ju0+k2d)*ion_e),
     &                 iz = 1+kl0*ion_e,1+ku0*ion_e)
               endif
               if(nvarcorn.gt.0) then
               read (iunit1) 
     &              ((((unk_n(ivar,ix,iy,iz,block_no),
     &                 ivar = 1,iv_n),
     &                 ix = 1+il0*ion_n,1+(iu0+1)*ion_n),
     &                 iy = 1+jl0*ion_n,1+(ju0+k2d)*ion_n),
     &                 iz = 1+kl0*ion_n,1+(ku0+k3d)*ion_n) 
               endif
#else /*  CHECKP_ON */
               if(nvar_chk_cc.gt.0) then
                 do ivar=1,nvar
                 if(checkp_on_cc(ivar)) then
                 do iz=1+ik0*ion_c,1+ku0*ion_c
                 do iy=1+ij0*ion_c,1+ju0*ion_c
                 do ix=1+il0*ion_c,1+iu0*ion_c
                   read (iunit1) unk(ivar,ix,iy,iz,block_no)
                 enddo
                 enddo
                 enddo
                 else
                   unk(ivar,:,:,:,block_no) = 0.
                 endif
                 enddo
               endif

               if(nvar_chk_fc.gt.0) then
                 do ivar=1,nfacevar
                 if(checkp_on_fc(1,ivar)) then
                 do iz = 1+kl0*ion_f,1+ku0*ion_f
                 do iy = 1+jl0*ion_f,1+ju0*ion_f
                 do ix = 1+il0*ion_f,1+(iu0+1)*ion_f
                   read (iunit1) facevarx(ivar,ix,iy,iz,block_no)
                 enddo
                 enddo
                 enddo
                 else
                   facevarx(ivar,:,:,:,block_no) = 0.
                 endif
                 enddo

                 do ivar=1,nfacevar
                 if(checkp_on_fc(2,ivar)) then
                 do iz = 1+kl0*ion_f,1+ku0*ion_f
                 do iy = 1+jl0*ion_f,1+(ju0+k2d)*ion_f
                 do ix = 1+il0*ion_f,1+iu0*ion_f
                   read (iunit1) facevary(ivar,ix,iy,iz,block_no)
                 enddo
                 enddo
                 enddo
                 else
                   facevary(ivar,:,:,:,block_no) = 0.
                 endif
                 enddo

                 do ivar=1,nfacevar
                 if(checkp_on_fc(3,ivar)) then
                 do iz = 1+kl0*ion_f,1+(ku0+k3d)*ion_f
                 do iy = 1+jl0*ion_f,1+ju0*ion_f
                 do ix = 1+il0*ion_f,1+iu0*ion_f
                   read (iunit1) facevarz(ivar,ix,iy,iz,block_no)
                 enddo
                 enddo
                 enddo
                 else
                   facevarz(ivar,:,:,:,block_no) = 0.
                 endif
                 enddo

               endif

               if(nvar_chk_ec.gt.0) then
                 do ivar=1,nvaredge
                 if(checkp_on_ec(1,ivar)) then
                 do iz = 1+kl0*ion_e,1+(ku0+k3d)*ion_e
                 do iy = 1+jl0*ion_e,1+(ju0+k2d)*ion_e
                 do ix = 1+il0*ion_e,1+iu0*ion_e
                   read (iunit1) unk_e_x(ivar,ix,iy,iz,block_no)
                 enddo
                 enddo
                 enddo
                 else
                   unk_e_x(ivar,:,:,:,block_no) = 0.
                 endif
                 enddo

                 do ivar=1,nvaredge
                 if(checkp_on_ec(2,ivar)) then
                 do iz = 1+kl0*ion_e,1+(ku0+k3d)*ion_e
                 do iy = 1+jl0*ion_e,1+ju0*ion_e
                 do ix = 1+il0*ion_e,1+(iu0+1)*ion_e
                   read (iunit1) unk_e_y(ivar,ix,iy,iz,block_no)
                 enddo
                 enddo
                 enddo
                 else
                   unk_e_y(ivar,:,:,:,block_no) = 0.
                 endif
                 enddo

                 do ivar=1,nvaredge
                 if(checkp_on_ec(3,ivar)) then
                 do iz = 1+kl0*ion_e,1+ku0*ion_e
                 do iy = 1+jl0*ion_e,1+(ju0+k2d)*ion_e
                 do ix = 1+il0*ion_e,1+(iu0+1)*ion_e
                   read (iunit1) unk_e_z(ivar,ix,iy,iz,block_no)
                 enddo
                 enddo
                 enddo
                 else
                   unk_e_z(ivar,:,:,:,block_no) = 0.
                 endif
                 enddo

               endif

               if(nvar_chk_nc.gt.0) then
                 do ivar=1,nvarcorn
                 if(checkp_on_nc(ivar)) then
                 do iz = 1+kl0*ion_n,1+(ku0+k3d)*ion_n
                 do iy = 1+jl0*ion_n,1+(ju0+k2d)*ion_n
                 do ix = 1+il0*ion_n,1+(iu0+1)*ion_n
                   read (iunit1) unk_n(ivar,ix,iy,iz,block_no)
                 enddo
                 enddo
                 enddo
                 else
                   unk_n(ivar,:,:,:,block_no) = 0.
                 endif
                 enddo
               endif

#endif /* CHECKP_ON */

            enddo
      if(ndim.lt.3) then
        bnd_box(1,3,:) = 0.
        bnd_box(2,3,:) = 1.
        coord(3,:) = .5*(bnd_box(2,3,:)+bnd_box(1,3,:))
      endif
      if(ndim.lt.2) then
        bnd_box(1,2,:) = 0.
        bnd_box(2,2,:) = 1.
        coord(2,:) = .5*(bnd_box(2,2,:)+bnd_box(1,2,:))
      endif

      endif

      if (iproc .gt. 0) then
!
! Post receives on pe iproc for messages from proc 0

         lnblockst = glnblocks(iproc)

         do block_no = 1, lnblockst

            isrc = 0
            itag = block_no
            isize = len_block 
     .        + nbndvar*(len_blockfx + len_blockfy + len_blockfz) 
     .        + nbndvare*(len_blockex + len_blockey + len_blockez) 
     .        + nbndvarc*len_blockn 

            call Mpi_recv(CR_buffer2,isize,MPI_REAL,
     .       isrc,itag,MPI_COMM_WORLD,status,ierr)

            position2 = 0

               if(nvar.gt.0)
     &           call Mpi_real_unpack(CR_buffer2,
     &             buf_dim_bytes2,position2,
     &             unk(1,1,1,1,block_no),len_block,
     &             MPI_REAL,
     &             MPI_COMM_WORLD,ierr)
               if(nfacevar.gt.0) then
                 call Mpi_real_unpack(CR_buffer2,
     &             buf_dim_bytes2,position2,
     &             facevarx(1,1,1,1,block_no),nbndvar*len_blockfx,
     &             MPI_REAL,MPI_COMM_WORLD,ierr)
                 call Mpi_real_unpack(CR_buffer2,
     &             buf_dim_bytes2,position2,
     &             facevary(1,1,1,1,block_no),nbndvar*len_blockfy,
     &             MPI_REAL,MPI_COMM_WORLD,ierr)
                 call Mpi_real_unpack(CR_buffer2,
     &             buf_dim_bytes2,position2,
     &             facevarz(1,1,1,1,block_no),nbndvar*len_blockfz,
     &             MPI_REAL,MPI_COMM_WORLD,ierr)
               endif
               if(nvaredge.gt.0) then
                 call Mpi_real_unpack(CR_buffer2,
     &             buf_dim_bytes2,position2,
     &             unk_e_x(1,1,1,1,block_no),nbndvare*len_blockex,
     &             MPI_REAL,MPI_COMM_WORLD,ierr)
                 call Mpi_real_unpack(CR_buffer2,
     &             buf_dim_bytes2,position2,
     &             unk_e_y(1,1,1,1,block_no),nbndvare*len_blockey,
     &             MPI_REAL,MPI_COMM_WORLD,ierr)
                 call Mpi_real_unpack(CR_buffer2,
     &             buf_dim_bytes2,position2,
     &             unk_e_z(1,1,1,1,block_no),nbndvare*len_blockez,
     &             MPI_REAL,MPI_COMM_WORLD,ierr)
               endif
               if(nvarcorn.gt.0) then
                 call Mpi_real_unpack(CR_buffer2,
     &             buf_dim_bytes2,position2,
     &             unk_n(1,1,1,1,block_no),nbndvarc*len_blockn,
     &             MPI_REAL,MPI_COMM_WORLD,ierr)
               endif
#ifdef DEBUG
               write(*,*) 'unk for blk ',block_no,' put into '
     .         ,'pe ',iproc
#endif /* DEBUG */
            enddo

            isrc = 0
            idest= iproc
            itag = (idest+1)*(maxblocks+1)
            isize = lnblockst*( 3+mflags+3*mdim+(nfaces+1+nchild)
     .                          + 1 )
#ifdef DEBUG
       write(*,*) 'pe ',iproc,' about to post receive from pe 0 ',
     .                ' tag ',itag
#endif /* DEBUG */
            call Mpi_recv(CR_buffer1,isize,MPI_REAL,
     .       isrc,itag,MPI_COMM_WORLD,status,ierr)
#ifdef DEBUG
            write(*,*) 'pe ',iproc,' posted receive from pe 0 ',
     .                ' tag ',itag
#endif /* DEBUG */

            position = 0
            do block_no = 1,lnblockst
#ifdef DEBUG
            write(*,*) 'pe ',iproc,' begin unpack blk ',block_no
#endif /* DEBUG */

! fetch data for this block
            call Mpi_int_unpack(CR_buffer1,buf_dim_bytes1,position,
     &        lrefine(block_no),1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
            call Mpi_int_unpack(CR_buffer1,buf_dim_bytes1,position,
     &        nodetype(block_no),1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
            call Mpi_int_unpack(CR_buffer1,buf_dim_bytes1,position,
     &        which_child(block_no),1,MPI_INTEGER,MPI_COMM_WORLD,
     &             ierr)
            call Mpi_int_unpack(CR_buffer1,buf_dim_bytes1,position,
     &        bflags(1,block_no),mflags,MPI_INTEGER,
     &        MPI_COMM_WORLD,ierr)
            call Mpi_real_unpack(CR_buffer1,buf_dim_bytes1,position,
     &        coord(1,block_no),mdim,MPI_REAL,
     &        MPI_COMM_WORLD,ierr)
            call Mpi_real_unpack(CR_buffer1,buf_dim_bytes1,position,
     &        bnd_box(1,1,block_no),2*mdim,
     &        MPI_REAL,
     &        MPI_COMM_WORLD,ierr)
            call Mpi_real_unpack(CR_buffer1,buf_dim_bytes1,position,
     &        work_block(block_no),1,MPI_REAL,
     &        MPI_COMM_WORLD,ierr)
            call Mpi_int_unpack(CR_buffer1,buf_dim_bytes1,position,
     &        gid(1,block_no),(nfaces+1+nchild),MPI_INTEGER,
     &        MPI_COMM_WORLD,ierr)

         end do

      endif


      if(iproc.eq.0) then

! Post sends from pe 0 for messages to all other procs

        do jproc = 1,nprocs-1
            lnblockst = glnblocks(jproc)

            position = 0

#ifdef DEBUG
      write(*,*) 'pe ',iproc,' computed isize '
#endif /* DEBUG */
            do block_no = 1,lnblockst

! Read in data for this block

               read (iunit1) 
     &              lrefinet(block_no),
     &              nodetypet(block_no),
     &              which_childt(block_no),
     &              (gidt(j,block_no),j=1,nfaces+1+nchild),
     &              (bflagst(j,block_no),j=1,mflags),
     &              (coordt(j,block_no),j=1,ndim),
     &              (bnd_boxt(1,j,block_no),j=1,ndim),
     &              (bnd_boxt(2,j,block_no),j=1,ndim),
     &              work_blockt(block_no)
#ifndef CHECKP_ON
               if(nvar.gt.0) then
               read (iunit1) 
     &              ((((unkt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_c),
     &                 ix = 1+il0*ion_c,1+iu0*ion_c),
     &                 iy = 1+jl0*ion_c,1+ju0*ion_c),
     &                 iz = 1+kl0*ion_c,1+ku0*ion_c)
               endif
               if(nfacevar.gt.0) then
               read (iunit1) 
     &              ((((facevarxt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_f),
     &                 ix = 1+il0*ion_f,1+(iu0+1)*ion_f),
     &                 iy = 1+jl0*ion_f,1+ju0*ion_f),
     &                 iz = 1+kl0*ion_f,1+ku0*ion_f)     
               read (iunit1) 
     &              ((((facevaryt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_f),
     &                 ix = 1+il0*ion_f,1+iu0*ion_f),
     &                 iy = 1+jl0*ion_f,1+(ju0+k2d)*ion_f),
     &                 iz = 1+kl0*ion_f,1+ku0*ion_f)    
               read (iunit1) 
     &              ((((facevarzt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_f),
     &                 ix = 1+il0*ion_f,1+iu0*ion_f),
     &                 iy = 1+jl0*ion_f,1+ju0*ion_f),
     &                 iz = 1+kl0*ion_f,1+(ku0+k3d)*ion_f)
               endif
               if(nvaredge.gt.0) then
               read (iunit1) 
     &              ((((unk_e_xt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_e),
     &                 ix = 1+il0*ion_e,1+iu0*ion_e),
     &                 iy = 1+jl0*ion_e,1+(ju0+k2d)*ion_e),
     &                 iz = 1+kl0*ion_e,1+(ku0+k3d)*ion_e)
               read (iunit1) 
     &              ((((unk_e_yt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_e),
     &                 ix = 1+il0*ion_e,1+(iu0+1)*ion_e),
     &                 iy = 1+jl0*ion_e,1+ju0*ion_e),
     &                 iz = 1+kl0*ion_e,1+(ku0+k3d)*ion_e)
               read (iunit1) 
     &              ((((unk_e_zt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_e),
     &                 ix = 1+il0*ion_e,1+(iu0+1)*ion_e),
     &                 iy = 1+jl0*ion_e,1+(ju0+k2d)*ion_e),
     &                 iz = 1+kl0*ion_e,1+ku0*ion_e)
               endif
               if(nvarcorn.gt.0) then
               read (iunit1) 
     &              ((((unk_nt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_n),
     &                 ix = 1+il0*ion_n,1+(iu0+1)*ion_n),
     &                 iy = 1+jl0*ion_n,1+(ju0+k2d)*ion_n),
     &                 iz = 1+kl0*ion_n,1+(ku0+k3d)*ion_n) 
               endif
#else /*  CHECKP_ON */
               if(nvar_chk_cc.gt.0) then
                 do ivar=1,nvar
                 if(checkp_on_cc(ivar)) then
                 do iz=1+ik0*ion_c,1+ku0*ion_c
                 do iy=1+ij0*ion_c,1+ju0*ion_c
                 do ix=1+il0*ion_c,1+iu0*ion_c
                   read (iunit1) unkt(ivar,ix,iy,iz)
                 enddo
                 enddo
                 enddo
                 else
                   unkt(ivar,:,:,:) = 0.
                 endif
                 enddo
               endif

               if(nvar_chk_fc.gt.0) then
                 do ivar=1,nfacevar
                 if(checkp_on_fc(1,ivar)) then
                 do iz = 1+kl0*ion_f,1+ku0*ion_f
                 do iy = 1+jl0*ion_f,1+ju0*ion_f
                 do ix = 1+il0*ion_f,1+(iu0+1)*ion_f
                   read (iunit1) facevarxt(ivar,ix,iy,iz)
                 enddo
                 enddo
                 enddo
                 else
                   facevarxt(ivar,:,:,:) = 0.
                 endif
                 enddo

                 do ivar=1,nfacevar
                 if(checkp_on_fc(2,ivar)) then
                 do iz = 1+kl0*ion_f,1+ku0*ion_f
                 do iy = 1+jl0*ion_f,1+(ju0+k2d)*ion_f
                 do ix = 1+il0*ion_f,1+iu0*ion_f
                   read (iunit1) facevaryt(ivar,ix,iy,iz)
                 enddo
                 enddo
                 enddo
                 else
                   facevaryt(ivar,:,:,:) = 0.
                 endif
                 enddo

                 do ivar=1,nfacevar
                 if(checkp_on_fc(3,ivar)) then
                 do iz = 1+kl0*ion_f,1+(ku0+k3d)*ion_f
                 do iy = 1+jl0*ion_f,1+ju0*ion_f
                 do ix = 1+il0*ion_f,1+iu0*ion_f
                   read (iunit1) facevarzt(ivar,ix,iy,iz)
                 enddo
                 enddo
                 enddo
                 else
                   facevarzt(ivar,:,:,:) = 0.
                 endif
                 enddo

               endif

               if(nvar_chk_ec.gt.0) then
                 do ivar=1,nvaredge
                 if(checkp_on_ec(1,ivar)) then
                 do iz = 1+kl0*ion_e,1+(ku0+k3d)*ion_e
                 do iy = 1+jl0*ion_e,1+(ju0+k2d)*ion_e
                 do ix = 1+il0*ion_e,1+iu0*ion_e
                   read (iunit1) unk_e_xt(ivar,ix,iy,iz)
                 enddo
                 enddo
                 enddo
                 else
                   unk_e_xt(ivar,:,:,:) = 0.
                 endif
                 enddo

                 do ivar=1,nvaredge
                 if(checkp_on_ec(2,ivar)) then
                 do iz = 1+kl0*ion_e,1+(ku0+k3d)*ion_e
                 do iy = 1+jl0*ion_e,1+ju0*ion_e
                 do ix = 1+il0*ion_e,1+(iu0+1)*ion_e
                   read (iunit1) unk_e_yt(ivar,ix,iy,iz)
                 enddo
                 enddo
                 enddo
                 else
                   unk_e_yt(ivar,:,:,:) = 0.
                 endif
                 enddo

                 do ivar=1,nvaredge
                 if(checkp_on_ec(3,ivar)) then
                 do iz = 1+kl0*ion_e,1+ku0*ion_e
                 do iy = 1+jl0*ion_e,1+(ju0+k2d)*ion_e
                 do ix = 1+il0*ion_e,1+(iu0+1)*ion_e
                   read (iunit1) unk_e_zt(ivar,ix,iy,iz)
                 enddo
                 enddo
                 enddo
                 else
                   unk_e_zt(ivar,:,:,:) = 0.
                 endif
                 enddo

               endif

               if(nvar_chk_nc.gt.0) then
                 do ivar=1,nvarcorn
                 if(checkp_on_nc(ivar)) then
                 do iz = 1+kl0*ion_n,1+(ku0+k3d)*ion_n
                 do iy = 1+jl0*ion_n,1+(ju0+k2d)*ion_n
                 do ix = 1+il0*ion_n,1+(iu0+1)*ion_n
                   read (iunit1) unk_nt(ivar,ix,iy,iz)
                 enddo
                 enddo
                 enddo
                 else
                   unk_nt(ivar,:,:,:) = 0.
                 endif
                 enddo
               endif

#endif /* CHECKP_ON */


      if(ndim.lt.3) then
        bnd_boxt(1,3,:) = 0.
        bnd_boxt(2,3,:) = 1.
        coordt(3,:) = .5*(bnd_boxt(2,3,:)+bnd_boxt(1,3,:))
      endif
      if(ndim.lt.2) then
        bnd_boxt(1,2,:) = 0.
        bnd_boxt(2,2,:) = 1.
        coordt(2,:) = .5*(bnd_boxt(2,2,:)+bnd_boxt(1,2,:))
      endif

               call Mpi_int_pack(lrefinet(block_no),1,MPI_INTEGER,
     &           CS_buffer1,buf_dim_bytes1,position,
     &           MPI_COMM_WORLD,ierr)
               call Mpi_int_pack(nodetypet(block_no),1,MPI_INTEGER,
     &           CS_buffer1,buf_dim_bytes1,position,
     &           MPI_COMM_WORLD,ierr)
               call Mpi_int_pack(which_childt(block_no),1,MPI_INTEGER,
     &           CS_buffer1,buf_dim_bytes1,position,
     &           MPI_COMM_WORLD,ierr)
               call Mpi_int_pack(bflagst(1,block_no),mflags,MPI_INTEGER,
     &           CS_buffer1,buf_dim_bytes1,position,
     &           MPI_COMM_WORLD,ierr)
               call Mpi_real_pack(coordt(1,block_no),mdim,
     &           MPI_REAL,
     &           CS_buffer1,buf_dim_bytes1,position,
     &           MPI_COMM_WORLD,ierr)
               call Mpi_real_pack(bnd_boxt(1,1,block_no),2*mdim,
     &           MPI_REAL,
     &           CS_buffer1,buf_dim_bytes1,position,
     &           MPI_COMM_WORLD,ierr)
               call Mpi_real_pack(work_blockt(block_no),1,
     &           MPI_REAL,
     &           CS_buffer1,buf_dim_bytes1,position,
     &           MPI_COMM_WORLD,ierr)
               call Mpi_int_pack(gidt(1,block_no),(nfaces+1+nchild),
     &           MPI_INTEGER,CS_buffer1,buf_dim_bytes1,position,
     &           MPI_COMM_WORLD,ierr)

              position2 = 0

              if(nvar.gt.0)
     &          call Mpi_real_pack(unkt(1,1,1,1),
     &            len_block,MPI_REAL,CS_buffer2,
     &            buf_dim_bytes2,position2,MPI_COMM_WORLD,ierr)
              if(nfacevar.gt.0) then
                call Mpi_real_pack(facevarxt(1,1,1,1),
     &            len_blockfx*nbndvar,
     &            MPI_REAL,CS_buffer2,
     &            buf_dim_bytes2,position2,MPI_COMM_WORLD,ierr)
                call Mpi_real_pack(facevaryt(1,1,1,1),
     &            len_blockfy*nbndvar,
     &            MPI_REAL,CS_buffer2,
     &            buf_dim_bytes2,position2,MPI_COMM_WORLD,ierr)
                call Mpi_real_pack(facevarzt(1,1,1,1),
     &            len_blockfz*nbndvar,
     &            MPI_REAL,CS_buffer2,
     &            buf_dim_bytes2,position2,MPI_COMM_WORLD,ierr)
               endif
              if(nvaredge.gt.0) then
                call Mpi_real_pack(unk_e_xt(1,1,1,1),
     &            len_blockex*nbndvare,
     &            MPI_REAL,CS_buffer2,
     &            buf_dim_bytes2,position2,MPI_COMM_WORLD,ierr)
                call Mpi_real_pack(unk_e_yt(1,1,1,1),
     &            len_blockey*nbndvare,
     &            MPI_REAL,CS_buffer2,
     &            buf_dim_bytes2,position2,MPI_COMM_WORLD,ierr)
                call Mpi_real_pack(unk_e_zt(1,1,1,1),
     &            len_blockez*nbndvare,
     &            MPI_REAL,CS_buffer2,
     &            buf_dim_bytes2,position2,MPI_COMM_WORLD,ierr)
               endif
              if(nvarcorn.gt.0) then
                call Mpi_real_pack(unk_nt(1,1,1,1),
     &            len_blockn*nbndvarc,
     &            MPI_REAL,CS_buffer2,
     &            buf_dim_bytes2,position2,MPI_COMM_WORLD,ierr)
               endif

               idest= jproc
               itag = block_no
               isize = len_block 
     .           + nbndvar*(len_blockfx + len_blockfy + len_blockfz) 
     .           + nbndvare*(len_blockex + len_blockey + len_blockez) 
     .           + nbndvarc*len_blockn 

               call Mpi_send(CS_buffer2,isize,MPI_REAL,
     .                       idest,itag,MPI_COMM_WORLD,ierr)

            enddo   ! end loop over blocks

            idest= jproc
            itag = (idest+1)*(maxblocks+1)
            isize = lnblockst*( 3+mflags+3*mdim+(nfaces+1+nchild)
     .                          + 1 )

            call Mpi_send(CS_buffer1,isize,MPI_REAL,
     .                    idest,itag,MPI_COMM_WORLD,ierr)

        enddo

        close(iunit1)

      endif

 2    continue

#ifdef DEBUG
      write(*,*) 'pe ',iproc,' starting bsize calc ',
     .                ' lnblocks ',lnblocks
#endif /* DEBUG */
 
      do block_no = 1,lnblocks
        bsize(:,block_no) = bnd_box(2,:,block_no)-
     $                     bnd_box(1,:,block_no)
      enddo

#ifdef DEBUG
      write(*,*) 'pe ',iproc,' completed bsize calc '
#endif /* DEBUG */
      

#endif /* MPI_USED */
      call shmem_barrier_all()

#ifdef DEBUG
      write(*,*) 'pe ',iproc,' starting tree rebuild '
#endif /* DEBUG */
 
! all processors fetch error code from proc 0

#ifdef MPI_USED
      ibcast_data(1) = ierr_read
      call mpi_bcast(ibcast_data,1,MPI_INTEGER,0,
     &                 MPI_COMM_WORLD,ierr)
      ierr_read = ibcast_data(1)
      ierr_readt = ierr_read
      if (ierr_readt.eq.1) 
     &   call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
#else
      call SHMEM_INTEGER_GET (ierr_readt,ierr_read,1,0)
      if (ierr_readt.eq.1) call amr_abort
#endif /* MPI_USED */


! COMPUTE TREE DATA FROM gid

#ifdef MPI_USED
      ibcast_data(1) = alnblocks
      call mpi_bcast(ibcast_data,1,MPI_INTEGER,0,
     &                 MPI_COMM_WORLD,ierr)
      alnblocks = ibcast_data(1)
      alnblockst = alnblocks
#else
      call SHMEM_INTEGER_GET (alnblockst,alnblocks,1,0)
      alnblocks = alnblockst
#endif /* MPI_USED */

      do block_no = 1,lnblocks

! neighbor data
         ngid = 0
         do j = 1,nfaces
            ngid = ngid + 1
            if (gid(ngid,block_no).gt.0) then
               neigh(2,j,block_no) = 
     $              int((gid(ngid,block_no)-1)/alnblocks)
               if (neigh(2,j,block_no).gt.nprocs-1) 
     $              neigh(2,j,block_no) = nprocs - 1
               neigh(1,j,block_no) = gid(ngid,block_no) - 
     $              (alnblocks*neigh(2,j,block_no))
            else
               neigh(1,j,block_no) = gid(ngid,block_no)
               neigh(2,j,block_no) = gid(ngid,block_no)
            end if
         end do
         
! parent data
         ngid = ngid + 1
         if (gid(ngid,block_no).gt.0) then
            parent(2,block_no) = 
     $           int((gid(ngid,block_no)-1)/alnblocks)
            if (parent(2,block_no).gt.nprocs-1) 
     $           parent(2,block_no) = nprocs - 1
            parent(1,block_no) = gid(ngid,block_no) - 
     $           (alnblocks*parent(2,block_no))
         else
            parent(1,block_no) = gid(ngid,block_no)
            parent(2,block_no) = gid(ngid,block_no)
         end if

! children data
         do j = 1,nchild
            ngid = ngid + 1
            if (gid(ngid,block_no).gt.0) then
               child(2,j,block_no) = 
     $              int((gid(ngid,block_no)-1)/alnblocks)
               if (child(2,j,block_no).gt.nprocs-1) 
     $              child(2,j,block_no) = nprocs - 1
               child(1,j,block_no) = gid(ngid,block_no) - 
     $              (alnblocks*child(2,j,block_no))
            else
               child(1,j,block_no) = gid(ngid,block_no)
               child(2,j,block_no) = gid(ngid,block_no)
            end if
         end do
         
      end do



      if(allocated(glnblocks)) deallocate( glnblocks )
      if(allocated(CS_buffer1)) deallocate( CS_buffer1 )
      if(allocated(CR_buffer1)) deallocate( CR_buffer1 )
      if(allocated(CS_buffer2)) deallocate( CS_buffer2 )
      if(allocated(CR_buffer2)) deallocate( CR_buffer2 )



! Now reorder blocks such that they are better balanced
! NOTE: this assumes that the total number of blocks is > nprocs

#ifdef MPI_USED
      lnblocks_old = lnblocks
      l_move_solution = .true.
      call amr_morton_order (lnblocks_old,nprocs,iproc,
     .                       l_move_solution)

      write(*,*) 'after amr_morton_order : pe lnb ',mype,lnblocks
#else /* MPI_USED */
      call amr_morton_order (nprocs,lnblocks)
#endif /* MPI_USED */


!---------------------------------------------
! compute grid_xmax, etc
#ifdef TIMING_MPI
      time1 = mpi_wtime()
#endif /* TIMING_MPI */
      call mpi_amr_global_domain_limits
#ifdef TIMING_MPI
      timer_amr_global_domain_limits =
     .     timer_amr_global_domain_limits + mpi_wtime() - time1
#endif /* TIMING_MPI */

#ifdef SAVE_MORTS
! Compute xmin,ymin,zmin,xmax,ymax,zmax or get them from storage
      xmin = grid_xmin
      ymin = grid_ymin
      zmin = grid_zmin
      xmax = grid_xmax
      ymax = grid_ymax
      zmax = grid_zmax

      write(*,*) 'checkre xmin etc ',xmin,ymin,zmin,xmax,ymax,zmax
      write(*,*) 'checkre lperiodicxtc ',
     .                        lperiodicx,lperiodicy,lperiodicz
      do lb = 1,lnblocks
        call morton_neighbors(xmin,ymin,zmin,xmax,ymax,zmax,
     .                        lperiodicx,lperiodicy,lperiodicz,
     .                        coord(:,lb),bsize(:,lb),ndim,
     .                        lrefine(lb),lrefine_max,mort_neigh)
        surr_morts(:,:,:,:,lb) =
     .                    mort_neigh(:,:,2-k2d:2+k2d,2-k3d:2+k3d)
      write(*,*) 'check_re set up surr_morts for initial blk ',lb,mype,
     .            ' surr_morts(6,:,:,1,lb) ',surr_morts(6,:,:,1,lb)
      enddo
#endif /* SAVE_MORTS */

      call amr_morton_process()

!
! St up an array of cell sizes for each grid refinement level.
! These can be used to minimize variation due to roundoff, but
! should ONLY be used with a uniformly spaced grid.
      level_cell_sizes = 0.
      level_cell_sizes(1,1) = (grid_xmax-grid_xmin)/real(nxb)
      if(ndim.gt.1)
     .  level_cell_sizes(2,1) = (grid_ymax-grid_ymin)/real(nyb)
      if(ndim.eq.3)
     .  level_cell_sizes(3,1) = (grid_zmax-grid_zmin)/real(nzb)
      do i=2,lrefine_max
        level_cell_sizes(1:ndim,i) = .5*level_cell_sizes(1:ndim,i-1)
      enddo
!---------------------------------------------


!
! mark grid as changed
      grid_changed = 1
      grid_analysed_mpi = 1


! Now make sure guardcell information is up to date

      call amr_guardcell(iproc,1,nguard)

      call mpi_amr_boundary_block_info(mype,nprocs)

      return
      end subroutine amr_checkpoint_re
