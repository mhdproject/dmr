!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"


! Development note : sizet is not required to be communicated
! This has already been removed from the mpi version.




      subroutine amr_checkpoint_wr (iunit1)




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

      use paramesh_interfaces, only : comm_real_min_to_all,
     &                                comm_real_max_to_all

      implicit none

#include "amr_shmem.fh"

      integer, intent(in) :: iunit1

      integer nguard0
      parameter(nguard0 = nguard*npgs)


      integer block_no,shmem_n_pes,shmem_my_pe
      integer jproc,i,j,ivar,ix,iy,iz,nprocs,iproc,iunit2

      integer lnblockst
      integer nrec,nrecl,ngid
      integer*4 inrecl
      integer n_to_left(0:2047) ! this array must extend from 0 to nprocs-1
      integer gid(nfaces+1+nchild,maxblocks)
      integer ierr

! TEMPORARIES WHICH ARE WRITTEN !
      
      integer  tot_blocks
      integer  lrefinet(maxblocks),nodetypet(maxblocks)
      integer  which_childt(maxblocks)
      integer  gidt(nfaces+1+nchild,maxblocks)
      integer  bflagst(mflags,maxblocks)
      real  coordt(mdim,maxblocks),sizet(mdim,maxblocks)
      real  work_blockt(maxblocks)
      real  bnd_boxt(2,mdim,maxblocks)
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

      integer,save :: idatain(1),idataout(1)

      save gid,gidt,lrefinet,nodetypet,lnblockst,coordt,sizet
      save which_childt,work_blockt
      save n_to_left,bnd_boxt

      real xleft,deltax,xpos
      real simtime,dt
      save unkt,facevarxt,facevaryt,facevarzt
      save unk_nt,unk_e_xt,unk_e_yt,unk_e_zt

      call shmem_barrier_all()

      nprocs = shmem_n_pes()
      iproc  = shmem_my_pe()

! COMPUTE TOTAL NO. OF BLOCKS STORED TO THE 'LEFT' OF THIS PROCESSOR

      call shmem_barrier_all()

      do i = 0,nprocs-1
         idatain(1) = lnblocks
         call shmem_barrier_all()
         call SHMEM_INTEGER_GET(n_to_left(i),idatain,1,i)
      end do

      tot_blocks = 0
      do i = 0,nprocs-1
         tot_blocks = tot_blocks + n_to_left(i)
      end do
      
      do i = nprocs-1,1,-1
         n_to_left(i) = n_to_left(i-1)
      end do

      n_to_left(0) = 0
      do i = 2,nprocs-1
         n_to_left(i) = n_to_left(i) + n_to_left(i-1)
      end do

      call shmem_barrier_all()

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

! store lnblocks for later getting
      idatain(1) = lnblocks


      call shmem_barrier_all()

! NOW WRITE OUT THE DATA FROM PROC 0

      nrec = 0

! set limits on data arrays
       il0 = nguard0
       iu0 = nxb-1+nguard0
       jl0 = nguard0*k2d
       ju0 = (nyb-1+nguard0)*k2d
       kl0 = nguard0*k3d
       ku0 = (nzb-1+nguard0)*k3d

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

      idatain(1) = lnblocks
      call shmem_barrier_all()

      if (iproc .eq. 0) then

         nrecl = 2 + nfaces + nchild + 1 + ndim + ndim +
     &        (nvar*nxb*nyb*nzb)
     &        + (nbndvar*((nxb+1)*nyb*nzb))
     &        + (nbndvar*(nxb*(nyb+k2d)*nzb))
     &        + (nbndvar*(nxb*nyb*(nzb+k3d)))
     &        + (nbndvare*(nxb*(nyb+k2d)*(nzb+k3d)))
     &        + (nbndvare*((nxb+1)*nyb*(nzb+k3d)))
     &        + (nbndvare*((nxb+1)*(nyb+k2d)*nzb))
     &        + nbndvarc*((nxb+1)*(nyb+k2d)*(nzb+1))
  
         inrecl = 8*nrecl
         open(unit=iunit1,form='unformatted',status='unknown'
     $        )
!    $        ,access='direct',recl=inrecl)

         nrec = nrec + 1
!        write (iunit1,rec=nrec) tot_blocks
         write (iunit1) tot_blocks

         do jproc = 0,nprocs-1
            
! fetch lnblocks from other processors

       call SHMEM_INTEGER_GET(idataout(1),idatain(1),1,jproc)
       lnblockst = idataout(1)
       call SHMEM_INTEGER_GET(lrefinet(1),lrefine(1),lnblockst,jproc)
       call SHMEM_INTEGER_GET(nodetypet(1),nodetype(1),lnblockst,jproc)
       call SHMEM_INTEGER_GET(which_childt(1),which_child(1),
     $                             lnblockst,jproc)
       call SHMEM_INTEGER_GET(bflagst(1,1),bflags(1,1),
     $                         lnblockst*mflags,jproc)
           call SHMEM_REAL_GET(coordt,coord,mdim*lnblockst,jproc)
           call SHMEM_REAL_GET(bnd_boxt,bnd_box,2*mdim*lnblockst,
     $                          jproc)
           do block_no = 1,lnblockst
             sizet(:,block_no) = bnd_boxt(2,:,block_no) - 
     $                           bnd_boxt(1,:,block_no)
           enddo
           call SHMEM_REAL_GET (work_blockt,work_block,
     $                          lnblockst,jproc)
           call SHMEM_INTEGER_GET (gidt(1,1),gid(1,1),
     $          lnblockst*(nfaces+1+nchild),jproc)

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

!               write (iunit1,rec=nrec) 
               write (iunit1) 
     &              lrefinet(block_no),
     &              nodetypet(block_no),
     &              which_childt(block_no),
     &              (gidt(j,block_no),j=1,nfaces+1+nchild),
     &              (bflagst(j,block_no),j=1,mflags),
     &              (coordt(j,block_no),j=1,ndim),
     &              (sizet(j,block_no),j=1,ndim),
     &              (bnd_boxt(1,j,block_no),j=1,ndim),
     &              (bnd_boxt(2,j,block_no),j=1,ndim),
     &              work_blockt(block_no),
     &              ((((unkt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_c),
     &                 ix = 1+il0*ion_c,1+iu0*ion_c),
     &                 iy = 1+jl0*ion_c,1+ju0*ion_c),
     &                 iz = 1+kl0*ion_c,1+ku0*ion_c),
     &              ((((facevarxt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_f),
     &                 ix = 1+il0*ion_f,1+(iu0+1)*ion_f),
     &                 iy = 1+jl0*ion_f,1+ju0*ion_f),
     &                 iz = 1+kl0*ion_f,1+ku0*ion_f),
     &              ((((facevaryt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_f),
     &                 ix = 1+il0*ion_f,1+iu0*ion_f),
     &                 iy = 1+jl0*ion_f,1+(ju0+k2d)*ion_f),
     &                 iz = 1+kl0*ion_f,1+ku0*ion_f),
     &              ((((facevarzt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_f),
     &                 ix = 1+il0*ion_f,1+iu0*ion_f),
     &                 iy = 1+jl0*ion_f,1+ju0*ion_f),
     &                 iz = 1+kl0*ion_f,1+(ku0+k3d)*ion_f),
     &              ((((unk_e_xt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_e),
     &                 ix = 1+il0*ion_e,1+iu0*ion_e),
     &                 iy = 1+jl0*ion_e,1+(ju0+k2d)*ion_e),
     &                 iz = 1+kl0*ion_e,1+(ku0+k3d)*ion_e),
     &              ((((unk_e_yt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_e),
     &                 ix = 1+il0*ion_e,1+(iu0+1)*ion_e),
     &                 iy = 1+jl0*ion_e,1+ju0*ion_e),
     &                 iz = 1+kl0*ion_e,1+(ku0+k3d)*ion_e),
     &              ((((unk_e_zt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_e),
     &                 ix = 1+il0*ion_e,1+(iu0+1)*ion_e),
     &                 iy = 1+jl0*ion_e,1+(ju0+k2d)*ion_e),
     &                 iz = 1+kl0*ion_e,1+ku0*ion_e),
     &              ((((unk_nt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_n),
     &                 ix = 1+il0*ion_n,1+(iu0+1)*ion_n),
     &                 iy = 1+jl0*ion_n,1+(ju0+k2d)*ion_n),
     &                 iz = 1+kl0*ion_n,1+(ku0+k3d)*ion_n)
               nrec = nrec + 1
            end do
         enddo

         close(iunit1)

      end if ! if iproc == 0
      
      call shmem_barrier_all()

      return
      end subroutine amr_checkpoint_wr

!-----------------------------------------------------------------------
      

      subroutine amr_checkpoint_re (iunit1)

! Subroutine to read checkpoint file using AMR package.
! Writes out tree data structure and data stored in blocks
! Currently reads are done serially by processor 0 and data sent to
! other processors.
! Uses UNFORMATTED, DIRECT ACCESS

! Written: K. Olson 7/97


      use paramesh_dimensions
      use physicaldata
      use tree

      use paramesh_interfaces, only : amr_morton_order,
     &                                amr_gsurrounding_blks,
     &                                amr_guardcell


      implicit none

#include "amr_shmem.fh"

      integer, intent(in) :: iunit1

      integer nguard0
      parameter(nguard0 = nguard*npgs)


      integer l,ii,icoord

      integer block_no,shmem_n_pes,shmem_my_pe
      integer jproc,i,j,ivar,ix,iy,iz,nprocs,iproc,iunit2

      integer lnblockst,nrec,nrecl,ngid,lb
      integer*4 inrecl
      integer gid(nfaces+1+nchild,maxblocks)
      integer alnblocks,alnblockst
      integer ierr,ierrt

! TEMPORARIES WHICH ARE READ

      integer  tot_blocks
      integer  lrefinet(maxblocks),nodetypet(maxblocks)
      integer  which_childt(maxblocks)
      integer  gidt(nfaces+1+nchild,maxblocks)
      integer  bflagst(mflags,maxblocks)
      real  coordt(mdim,maxblocks),sizet(mdim,maxblocks)
      real  work_blockt(maxblocks)
      real  bnd_boxt(2,mdim,maxblocks)
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

      save gid,gidt,coordt,sizet,lrefinet,nodetypet,lnblockst
      save which_childt,work_blockt
      save alnblocks,alnblockst,bnd_boxt
      save unkt,facevarxt,facevaryt,facevarzt
      save unk_nt,unk_e_xt,unk_e_yt,unk_e_zt
      save ierr,ierrt

      integer,save :: idatain(1),idataout(1)
     
      integer parent_pe,parent_blk
      real    par_coord(3)
      save    par_coord,parent_pe,parent_blk


      real :: xmin,xmax,ymin,ymax,zmin,zmax
      real :: xmin1,xmax1,ymin1,ymax1,zmin1,zmax1

      logical :: l_move_solution

 
      call shmem_barrier_all()

      nprocs = shmem_n_pes()
      iproc  = shmem_my_pe()
      ierr = 0
      ierrt = 0

      call shmem_barrier_all()


! set limits on data arrays
       il0 = nguard0
       iu0 = nxb-1+nguard0
       jl0 = nguard0*k2d
       ju0 = (nyb-1+nguard0)*k2d
       kl0 = nguard0*k3d
       ku0 = (nzb-1+nguard0)*k3d

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


      if (iproc .eq. 0) then

         nrec = 0

         nrecl = 2 + nfaces + nchild + 1 + ndim + ndim +
     &        (nvar*nxb*nyb*nzb)
     &        + (nbndvar*((nxb+1)*nyb*nzb))
     &        + (nbndvar*(nxb*(nyb+k2d)*nzb))
     &        + (nbndvar*(nxb*nyb*(nzb+k3d)))
     &        + (nbndvare*(nxb*(nyb+k2d)*(nzb+k3d)))
     &        + (nbndvare*((nxb+1)*nyb*(nzb+k3d)))
     &        + (nbndvare*((nxb+1)*(nyb+k2d)*nzb))
     &        + nbndvarc*((nxb+1)*(nyb+k2d)*(nzb+1))
  

         inrecl=8*nrecl
         open(unit=iunit1,form='unformatted',status='unknown'
     $        )
!     $        ,access='direct',recl=inrecl)

         nrec = nrec + 1
!         read (iunit1,rec=nrec) tot_blocks
         read (iunit1) tot_blocks
         write(*,*) 'blocks to be input ',tot_blocks

! compute approximate lnblocks (this will be the number of blocks stored on
! processors 0 -> nprocs-2, nprocs-1 gets tot_blocks - the total number on the
! rest of the blocks)

         alnblocks = int(tot_blocks/nprocs)

! check for error
         if (tot_blocks-(alnblocks*(nprocs-1)).gt.maxblocks) then

          print *,' ************* ERROR in inhale_u: ***************'
          print *,' No. of blocks per processor exceeds maxblocks.'
          print *,' Suggest you reset maxblocks to a larger number or '
          print *,' run on a larger no. of processors. '

            ierr = 1

            go to 2

         end if

         do jproc = 0,nprocs-1

            if (jproc.lt.nprocs-1) then
               lnblockst = alnblocks
            else
               lnblockst = tot_blocks - (alnblocks*(nprocs-1))
            end if

            do block_no = 1,lnblockst

!            write(*,*) 'reading block no ',block_no
!     .              ,' jproc ',jproc,' nprocs ',nprocs

! Read in data for this block

               nrec = nrec + 1

!               read (iunit1,rec=nrec) 
               read (iunit1) 
     &              lrefinet(block_no),
     &              nodetypet(block_no),
     &              which_childt(block_no),
     &              (gidt(j,block_no),j=1,nfaces+1+nchild),
     &              (bflagst(j,block_no),j=1,mflags),
     &              (coordt(j,block_no),j=1,ndim),
     &              (sizet(j,block_no),j=1,ndim),
     &              (bnd_boxt(1,j,block_no),j=1,ndim),
     &              (bnd_boxt(2,j,block_no),j=1,ndim),
     &              work_blockt(block_no),
     &              ((((unkt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_c),
     &                 ix = 1+il0*ion_c,1+iu0*ion_c),
     &                 iy = 1+jl0*ion_c,1+ju0*ion_c),
     &                 iz = 1+kl0*ion_c,1+ku0*ion_c),
     &              ((((facevarxt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_f),
     &                 ix = 1+il0*ion_f,1+(iu0+1)*ion_f),
     &                 iy = 1+jl0*ion_f,1+ju0*ion_f),
     &                 iz = 1+kl0*ion_f,1+ku0*ion_f),     
     &              ((((facevaryt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_f),
     &                 ix = 1+il0*ion_f,1+iu0*ion_f),
     &                 iy = 1+jl0*ion_f,1+(ju0+k2d)*ion_f),
     &                 iz = 1+kl0*ion_f,1+ku0*ion_f),     
     &              ((((facevarzt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_f),
     &                 ix = 1+il0*ion_f,1+iu0*ion_f),
     &                 iy = 1+jl0*ion_f,1+ju0*ion_f),
     &                 iz = 1+kl0*ion_f,1+(ku0+k3d)*ion_f),
     &              ((((unk_e_xt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_e),
     &                 ix = 1+il0*ion_e,1+iu0*ion_e),
     &                 iy = 1+jl0*ion_e,1+(ju0+k2d)*ion_e),
     &                 iz = 1+kl0*ion_e,1+(ku0+k3d)*ion_e),
     &              ((((unk_e_yt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_e),
     &                 ix = 1+il0*ion_e,1+(iu0+1)*ion_e),
     &                 iy = 1+jl0*ion_e,1+ju0*ion_e),
     &                 iz = 1+kl0*ion_e,1+(ku0+k3d)*ion_e),
     &              ((((unk_e_zt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_e),
     &                 ix = 1+il0*ion_e,1+(iu0+1)*ion_e),
     &                 iy = 1+jl0*ion_e,1+(ju0+k2d)*ion_e),
     &                 iz = 1+kl0*ion_e,1+ku0*ion_e) ,
     &              ((((unk_nt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_n),
     &                 ix = 1+il0*ion_n,1+(iu0+1)*ion_n),
     &                 iy = 1+jl0*ion_n,1+(ju0+k2d)*ion_n),
     &                 iz = 1+kl0*ion_n,1+(ku0+k3d)*ion_n) 

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
     $              nvarcorn*len_blockn,jproc)

            end do

! put no. of blocks on processor jproc

            idatain(1) = lnblockst
            call SHMEM_INTEGER_PUT (idataout(1),idatain(1),1,jproc)

! put tree data on processor jproc

            call SHMEM_REAL_PUT (coord,coordt,mdim*lnblockst,jproc)
!           call SHMEM_REAL_PUT (size,sizet,mdim*lnblockst,jproc)
            call SHMEM_REAL_PUT (bnd_box,bnd_boxt,2*mdim*lnblockst,
     $                           jproc)
            call SHMEM_REAL_PUT (work_block,work_blockt,lnblockst,
     $                           jproc)
            call SHMEM_INTEGER_PUT (lrefine(1),lrefinet(1),
     $                              lnblockst,jproc)
            call SHMEM_INTEGER_PUT (nodetype(1),nodetypet(1),
     $                              lnblockst,jproc)
            call SHMEM_INTEGER_PUT (which_child(1),which_childt(1),
     $                              lnblockst,jproc)
            call SHMEM_INTEGER_PUT (bflags(1,1),bflagst(1,1),
     $                              lnblockst*mflags,jproc)
            call SHMEM_INTEGER_PUT (gid(1,1),gidt(1,1),
     $           lnblockst*(nfaces+1+nchild),jproc)

         end do

         close(iunit1)

      end if ! if iproc == 0

 2    continue
      call shmem_barrier_all()

      lnblocks = idataout(1)

      do block_no = 1,lnblocks
        bsize(:,block_no) = bnd_box(2,:,block_no)-
     $                     bnd_box(1,:,block_no)
      enddo
      

! COMPUTE TREE DATA FROM gid


      idatain(1) = alnblocks
      call shmem_barrier_all()
      call SHMEM_INTEGER_GET (idataout(1),idatain(1),1,0)
      alnblocks = idataout(1)


! all processors fetch error code from proc 0

      idatain(1) = ierr
      call shmem_barrier_all()
      call SHMEM_INTEGER_GET (idataout(1),idatain(1),1,0)
      ierrt = idataout(1)

      call shmem_barrier_all()
      if (ierrt.eq.1) then
        write(*,*) iproc,'checkpoint error : ierrt ',ierrt
        call amr_abort
      endif

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

      call shmem_barrier_all()


! Now reorder blocks such that they are better balanced
! NOTE: this assumes that the total number of blocks is > nprocs

      l_move_solution = .true.
      call amr_morton_order (nprocs,lnblocks,l_move_solution)

      call shmem_barrier_all()

#ifdef MAYBEOBSOLETE
! Reconstruct which_child array
      if(lnblocks.gt.0) then
      do lb = 1,lnblocks
        parent_blk = parent(1,lb)
        parent_pe  = parent(2,lb)
        if(parent_blk.gt.0) then
          call shmem_real_get( par_coord(1),coord(1,parent_blk),
     .                       3,parent_pe )
          ix = 1
          iy = 1
          iz = 1
          if(coord(1,lb).gt.par_coord(1)) ix = 2
          if(coord(2,lb).gt.par_coord(2).and.ndim.ge.2) iy = 2
          if(coord(3,lb).gt.par_coord(3).and.ndim.ge.3) iz = 2
          which_child(lb) = ix + 2*(iy-1) + 4*(iz-1)
        else
          which_child(lb) = -1
        endif
      enddo
      endif
#endif /* MAYBEOBSOLETE */


!
! Fill the surr_blks arrays
      call amr_gsurrounding_blks(iproc,.true.)

! Find the coordinate ranges
      xmin1 = minval(bnd_box(1,1,1:lnblocks))
      ymin1 = minval(bnd_box(1,2,1:lnblocks))
      zmin1 = minval(bnd_box(1,3,1:lnblocks))
      xmax1 = maxval(bnd_box(2,1,1:lnblocks))
      ymax1 = maxval(bnd_box(2,2,1:lnblocks))
      zmax1 = maxval(bnd_box(2,3,1:lnblocks))
      xmin = min(1.e30, xmin1)
      ymin = min(1.e30, ymin1)
      zmin = min(1.e30, zmin1)
      xmax = max(-1.e30, xmax1)
      ymax = max(-1.e30, ymax1)
      zmax = max(-1.e30, zmax1)
      call comm_real_min_to_all(grid_xmin,xmin)
      call comm_real_min_to_all(grid_ymin,ymin)
      call comm_real_min_to_all(grid_zmin,zmin)
      call comm_real_max_to_all(grid_xmax,xmax)
      call comm_real_max_to_all(grid_ymax,ymax)
      call comm_real_max_to_all(grid_zmax,zmax)
      if(iproc.eq.0) then
      write(*,*) ' grid_xmin ',grid_xmin
      write(*,*) ' grid_xmax ',grid_xmax
      write(*,*) ' grid_ymin ',grid_ymin
      write(*,*) ' grid_ymax ',grid_ymax
      write(*,*) ' grid_zmin ',grid_zmin
      write(*,*) ' grid_zmax ',grid_zmax
      endif


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


! mark grid as changed
      grid_changed = 1

      call shmem_barrier_all()

! Now make sure guardcell information is up to date

#ifndef NO_PERMANENT_GUARDCELL
      call amr_guardcell(iproc,1,nguard)
      call shmem_barrier_all()
#endif

      return
      end subroutine amr_checkpoint_re
