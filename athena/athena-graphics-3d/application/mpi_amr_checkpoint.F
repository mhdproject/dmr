#include "paramesh_preprocessor.fh"
      subroutine amr_checkpoint_wr_a (dt_in,time_in,loop_index_start,nchombo,
     &                                                       iunit1)

!#define DEBUG


! Subroutine to checkpoint runs using AMR package.
! Writes out tree data structure and data stored in blocks
! Currently writes are done serially by processor 0 collecting data from
! other processors and then writing it out.
! USES UNFORMATTED DIRECT I/O

! Written: K. Olson and C. Mobarry 7/97
! Modified: P. MacNeice 4/99


      use physcons
      use paramesh_dimensions
      use physicaldata
      use tree
      use interior_gbc

      use paramesh_interfaces
#ifdef MPI_USED
      use paramesh_mpi_interfaces
#endif /* MPI_USED */


      implicit none

#include "amr_shmem.fh"

#ifdef MPI_USED
      include 'mpif.h'
#endif /* MPI_USED */

      integer, intent(in) :: iunit1
      integer, intent(in) :: loop_index_start,nchombo
      real,    intent(in) :: dt_in,time_in

      integer nguard0
      parameter(nguard0 = nguard*npgs)


      integer block_no,shmem_n_pes,shmem_my_pe
      integer jproc,i,j,ivar,ix,iy,iz,nprocs,iproc,iunit2

      integer lnblockst
      integer nrec,nrecl,ngid
      integer*4 inrecl
      integer,dimension (:),  allocatable :: n_to_left
!      integer n_to_left(0:2047) ! this array must extend from 0 to nprocs-1
      integer gid(nfaces+1+nchild,maxblocks)
      integer ierr

! TEMPORARIES WHICH ARE WRITTEN !
      
      integer  tot_blocks
      integer  lrefinet(maxblocks),nodetypet(maxblocks)
      integer  which_childt(maxblocks)
      integer  gidt(nfaces+1+nchild,maxblocks)
      integer  bflagst(mflags,maxblocks)
      real  coordt(mdim,maxblocks)
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

      save gid,gidt,lrefinet,nodetypet,lnblockst,coordt
      save which_childt,work_blockt
      save n_to_left,bnd_boxt

      real xleft,deltax,xpos
      real simtime,dt
      save unkt,facevarxt,facevaryt,facevarzt
      save unk_nt,unk_e_xt,unk_e_yt,unk_e_zt

      real :: dtold
      common/old_timestep/dtold
      real :: rl,vxl,vyl,vzl,pl,bxl,byl,bzl
      common/left_state/rl,vxl,vyl,vzl,pl,bxl,byl,bzl


#ifdef MPI_USED
      integer,dimension (:),  allocatable :: glnblocks
      integer :: isrc,idest,itag,isize,ierror,position

      integer status(MPI_STATUS_SIZE)

      integer,dimension (:),  allocatable :: recvrequest
      integer,dimension (:),  allocatable :: sendrequest
      integer,dimension (:,:),allocatable :: recvstatus
      integer,dimension (:,:),allocatable :: sendstatus
      integer, parameter :: buf_dim_int = 
     .              maxblocks*( 3+mflags+(nfaces+1+nchild))
      integer, parameter :: buf_dim_real = 
     .              maxblocks*( 3*mdim + 1 + len_block 
     .       + nbndvar*(len_blockfx + len_blockfy + len_blockfz)
     .       + nbndvarc*len_blockn 
     .       + nbndvare*(len_blockex + len_blockey + len_blockez) )
      integer, parameter :: buf_dim = buf_dim_real + buf_dim_int
      integer :: no_of_bytes_per_real,no_of_bytes_per_integer
      integer :: buf_dim_bytes

      real :: CS_buffer(buf_dim), CR_buffer(buf_dim)
#endif /* MPI_USED */

      call shmem_barrier_all()

      nprocs = shmem_n_pes()
      iproc  = shmem_my_pe()


! COMPUTE TOTAL NO. OF BLOCKS STORED TO THE 'LEFT' OF THIS PROCESSOR

      call shmem_barrier_all()


      if(allocated(n_to_left)) deallocate( n_to_left )
      allocate ( n_to_left(0:nprocs-1) )

#ifdef MPI_USED
      if(allocated(glnblocks)) deallocate( glnblocks )
      allocate ( glnblocks(0:nprocs-1) )

      call mpi_type_size(MPI_INTEGER,no_of_bytes_per_integer,ierr)
      call mpi_type_size(MPI_DOUBLE_PRECISION
     .                  ,no_of_bytes_per_real   ,ierr)
      buf_dim_bytes = buf_dim_real*no_of_bytes_per_real +
     .                buf_dim_int *no_of_bytes_per_integer

      call MPI_Barrier(MPI_COMM_WORLD,ierror)

      glnblocks(iproc) = lnblocks
      call MPI_Allgather(glnblocks(iproc), 1,MPI_INTEGER,
     .                    glnblocks,1,MPI_INTEGER,
     .                    MPI_COMM_WORLD,ierror)
      n_to_left = glnblocks
#ifdef DEBUG
      write(*,*) 'n_to_left ',n_to_left
#endif /* DEBUG */
      call MPI_Barrier   (MPI_COMM_WORLD,ierror)
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

      call shmem_barrier_all()

! NOW WRITE OUT THE DATA FROM PROC 0

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

#ifndef MPI_USED

      nrec = 0

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
         write (iunit1) time_in 
         write (iunit1) loop_index_start,nchombo 

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

!               write (iunit1,rec=nrec) 
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
               write (iunit1) 
     &              ((((unkt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_c),
     &                 ix = 1+il0*ion_c,1+iu0*ion_c),
     &                 iy = 1+jl0*ion_c,1+ju0*ion_c),
     &                 iz = 1+kl0*ion_c,1+ku0*ion_c)
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
               if(nvarcorn.gt.0) 
               write (iunit1) 
     &              ((((unk_nt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_n),
     &                 ix = 1+il0*ion_n,1+(iu0+1)*ion_n),
     &                 iy = 1+jl0*ion_n,1+(ju0+k2d)*ion_n),
     &                 iz = 1+kl0*ion_n,1+(ku0+k3d)*ion_n)
               endif
               nrec = nrec + 1
            end do
         enddo

         close(iunit1)

      end if ! if iproc == 0
      
      call shmem_barrier_all()



#else /* MPI_USED */



#ifdef DEBUG
      write(*,*) 'pe ',iproc,' starting packing '
#endif /* DEBUG */
      call MPI_Barrier  (MPI_COMM_WORLD,ierror)

      if(allocated(recvrequest)) deallocate( recvrequest )
      allocate ( recvrequest(nprocs) )
      if(allocated(sendrequest)) deallocate( sendrequest )
      allocate ( sendrequest(nprocs) )
      if(allocated(recvstatus)) deallocate( recvstatus )
      allocate ( recvstatus(MPI_STATUS_SIZE,nprocs) )
      if(allocated(sendstatus)) deallocate( sendstatus )
      allocate ( sendstatus(MPI_STATUS_SIZE,nprocs) )

      CR_buffer = 0
      CS_buffer = 0

#ifdef DEBUG
      write(*,*) 'pe ',iproc,' initialized buffers '
#endif /* DEBUG */

      nrec = 0

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
         write (iunit1) time_in 
         write (iunit1) loop_index_start,nchombo 


! Write data from processor 0
         do block_no = 1,lnblocks

! fetch data for this block
!               write (iunit1,rec=nrec) 
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
               write (iunit1)
     &              ((((unk(ivar,ix,iy,iz,block_no),
     &                 ivar = 1,iv_c),
     &                 ix = 1+il0*ion_c,1+iu0*ion_c),
     &                 iy = 1+jl0*ion_c,1+ju0*ion_c),
     &                 iz = 1+kl0*ion_c,1+ku0*ion_c)
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
               nrec = nrec + 1
         enddo


         do jproc = 1,nprocs-1

!
! Post receives on pe 0 for messages from all other procs

! fetch lnblocks from other processors

            lnblockst = glnblocks(jproc)
            isrc = jproc
            idest= 0
            itag = isrc*nprocs + idest+1
            isize = lnblockst*( 3+mflags+3*mdim+(nfaces+1+nchild)
     .                          + 1 + len_block 
     .        + nbndvar*(len_blockfx + len_blockfy + len_blockfz) 
     .        + nbndvare*(len_blockex + len_blockey + len_blockez) 
     .        + nbndvarc*len_blockn )
#ifdef DEBUG
      write(*,*) 'pe ',iproc,' posting recv from proc ',jproc,
     .      ' tag ',itag
#endif /* DEBUG */
            call Mpi_recv(CR_buffer,isize,MPI_DOUBLE_PRECISION,
     .       isrc,itag,MPI_COMM_WORLD,status,ierr)

            position = 0
            call Mpi_unpack(CR_buffer,buf_dim_bytes,position,
     &          lrefinet,lnblockst,MPI_INTEGER,MPI_COMM_WORLD,ierr)
            call Mpi_unpack(CR_buffer,buf_dim_bytes,position,
     &          nodetypet,lnblockst,MPI_INTEGER,MPI_COMM_WORLD,ierr)
            call Mpi_unpack(CR_buffer,buf_dim_bytes,position,
     &          which_childt,lnblockst,MPI_INTEGER,MPI_COMM_WORLD,ierr)
            call Mpi_unpack(CR_buffer,buf_dim_bytes,position,
     &          bflagst,lnblockst*mflags,MPI_INTEGER,
     &          MPI_COMM_WORLD,ierr)
            call Mpi_unpack(CR_buffer,buf_dim_bytes,position,
     &          coordt,mdim*lnblockst,MPI_DOUBLE_PRECISION,
     &          MPI_COMM_WORLD,ierr)
            call Mpi_unpack(CR_buffer,buf_dim_bytes,position,
     &          bnd_boxt,2*mdim*lnblockst,
     &          MPI_DOUBLE_PRECISION,
     &          MPI_COMM_WORLD,ierr)
            call Mpi_unpack(CR_buffer,buf_dim_bytes,position,
     &          work_blockt,lnblockst,
     &          MPI_DOUBLE_PRECISION,
     &          MPI_COMM_WORLD,ierr)
            call Mpi_unpack(CR_buffer,buf_dim_bytes,position,
     &          gidt,lnblockst*(nfaces+1+nchild),MPI_INTEGER,
     &          MPI_COMM_WORLD,ierr)


            do block_no = 1,lnblockst

! fetch data for this block
               if(nvar.gt.0)
     &           call Mpi_unpack(CR_buffer,buf_dim_bytes,position,
     &             unkt,len_block,MPI_DOUBLE_PRECISION,
     &             MPI_COMM_WORLD,ierr)
               if(nfacevar.gt.0) then
                 call Mpi_unpack(CR_buffer,buf_dim_bytes,position,
     &             facevarxt,nbndvar*len_blockfx,
     &             MPI_DOUBLE_PRECISION,
     &             MPI_COMM_WORLD,ierr)
                 call Mpi_unpack(CR_buffer,buf_dim_bytes,position,
     &             facevaryt,nbndvar*len_blockfy,
     &             MPI_DOUBLE_PRECISION,
     &             MPI_COMM_WORLD,ierr)
                 call Mpi_unpack(CR_buffer,buf_dim_bytes,position,
     &             facevarzt,nbndvar*len_blockfz,
     &             MPI_DOUBLE_PRECISION,
     &             MPI_COMM_WORLD,ierr)
               endif
               if(nvaredge.gt.0) then
                 call Mpi_unpack(CR_buffer,buf_dim_bytes,position,
     &             unk_e_xt,nbndvare*len_blockex,
     &             MPI_DOUBLE_PRECISION,
     &             MPI_COMM_WORLD,ierr)
                 call Mpi_unpack(CR_buffer,buf_dim_bytes,position,
     &             unk_e_yt,nbndvare*len_blockey,
     &             MPI_DOUBLE_PRECISION,
     &             MPI_COMM_WORLD,ierr)
                 call Mpi_unpack(CR_buffer,buf_dim_bytes,position,
     &             unk_e_zt,nbndvare*len_blockez,
     &             MPI_DOUBLE_PRECISION,
     &             MPI_COMM_WORLD,ierr)
               endif
               if(nvarcorn.gt.0) then
                 call Mpi_unpack(CR_buffer,buf_dim_bytes,position,
     &             unk_nt,nbndvarc*len_blockn,
     &             MPI_DOUBLE_PRECISION,
     &             MPI_COMM_WORLD,ierr)
               endif

!               write (iunit1,rec=nrec) 
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
               write (iunit1) 
     &              ((((unkt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_c),
     &                 ix = 1+il0*ion_c,1+iu0*ion_c),
     &                 iy = 1+jl0*ion_c,1+ju0*ion_c),
     &                 iz = 1+kl0*ion_c,1+ku0*ion_c)
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

               nrec = nrec + 1

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
            itag = isrc*nprocs + idest+1
            isize = lnblockst*( 3+mflags+3*mdim+(nfaces+1+nchild)
     .                          + 1 + len_block 
     .        + nbndvar*(len_blockfx + len_blockfy + len_blockfz) 
     .        + nbndvare*(len_blockex + len_blockey + len_blockez) 
     .        + nbndvarc*len_blockn )

#ifdef DEBUG
      write(*,*) 'pe ',iproc,' packing data for send '
#endif /* DEBUG */
            position = 0
            call Mpi_pack(lrefine(1),lnblocks,MPI_INTEGER,CS_buffer,
     &          buf_dim_bytes,position,MPI_COMM_WORLD,ierr)
            call Mpi_pack(nodetype(1),lnblocks,MPI_INTEGER,CS_buffer,
     &          buf_dim_bytes,position,MPI_COMM_WORLD,ierr)
            call Mpi_pack(which_child(1),lnblocks,MPI_INTEGER,CS_buffer,
     &          buf_dim_bytes,position,MPI_COMM_WORLD,ierr)
            call Mpi_pack(bflags(1,1),lnblocks*mflags,MPI_INTEGER,
     &          CS_buffer,buf_dim_bytes,position,MPI_COMM_WORLD,ierr)
            call Mpi_pack(coord(1,1),mdim*lnblocks,
     &          MPI_DOUBLE_PRECISION,
     &          CS_buffer,buf_dim_bytes,position,MPI_COMM_WORLD,ierr)
!            call Mpi_pack(bsize(1,1),mdim*lnblocks,
!     &          MPI_DOUBLE_PRECISION,
!     &          CS_buffer,buf_dim_bytes,position,MPI_COMM_WORLD,ierr)
            call Mpi_pack(bnd_box(1,1,1),2*mdim*lnblocks,
     &          MPI_DOUBLE_PRECISION,
     &          CS_buffer,buf_dim_bytes,position,MPI_COMM_WORLD,ierr)
            call Mpi_pack(work_block(1),lnblocks,
     &          MPI_DOUBLE_PRECISION,
     &          CS_buffer,buf_dim_bytes,position,MPI_COMM_WORLD,ierr)
            call Mpi_pack(gid(1,1),lnblocks*(nfaces+1+nchild),
     &          MPI_INTEGER,CS_buffer,buf_dim_bytes,position,
     &          MPI_COMM_WORLD,ierr)
#ifdef DEBUG
      write(*,*) 'pe ',iproc,' tree data packed for send '
#endif /* DEBUG */

            do block_no= 1,lnblocks
              if(nvar.gt.0)
     &          call Mpi_pack(unk(1,1,1,1,block_no),
     &            len_block,MPI_DOUBLE_PRECISION,CS_buffer,
     &            buf_dim_bytes,position,MPI_COMM_WORLD,ierr)

              if(nfacevar.gt.0) then
                call Mpi_pack(facevarx(1,1,1,1,block_no),
     &            len_blockfx*nbndvar,
     &            MPI_DOUBLE_PRECISION,CS_buffer,
     &            buf_dim_bytes,position,MPI_COMM_WORLD,ierr)
                call Mpi_pack(facevary(1,1,1,1,block_no),
     &            len_blockfy*nbndvar,
     &            MPI_DOUBLE_PRECISION,CS_buffer,
     &            buf_dim_bytes,position,MPI_COMM_WORLD,ierr)
                call Mpi_pack(facevarz(1,1,1,1,block_no),
     &            len_blockfz*nbndvar,
     &            MPI_DOUBLE_PRECISION,CS_buffer,
     &            buf_dim_bytes,position,MPI_COMM_WORLD,ierr)
               endif
              if(nvaredge.gt.0) then
                call Mpi_pack(unk_e_x(1,1,1,1,block_no),
     &            len_blockex*nbndvare,
     &            MPI_DOUBLE_PRECISION,CS_buffer,
     &            buf_dim_bytes,position,MPI_COMM_WORLD,ierr)
                call Mpi_pack(unk_e_y(1,1,1,1,block_no),
     &            len_blockey*nbndvare,
     &            MPI_DOUBLE_PRECISION,CS_buffer,
     &            buf_dim_bytes,position,MPI_COMM_WORLD,ierr)
                call Mpi_pack(unk_e_z(1,1,1,1,block_no),
     &            len_blockez*nbndvare,
     &            MPI_DOUBLE_PRECISION,CS_buffer,
     &            buf_dim_bytes,position,MPI_COMM_WORLD,ierr)
               endif
              if(nvarcorn.gt.0) then
                call Mpi_pack(unk_n(1,1,1,1,block_no),
     &            len_blockn*nbndvarc,
     &            MPI_DOUBLE_PRECISION,CS_buffer,
     &            buf_dim_bytes,position,MPI_COMM_WORLD,ierr)
               endif


            enddo
 
#ifdef DEBUG
      write(*,*) 'pe ',iproc,' sending data : tag ',itag
#endif /* DEBUG */
            call Mpi_send(CS_buffer,isize,MPI_DOUBLE_PRECISION,
     .       idest,itag,MPI_COMM_WORLD,ierr)

#ifdef DEBUG
      write(*,*) 'pe ',iproc,' data sent '
#endif /* DEBUG */
       end if ! if iproc > 0

      call MPI_Barrier  (MPI_COMM_WORLD,ierror)
#endif /* MPI_USED */

!     write(*,701) ebox_xmin, ebox_xmax, iproc
!     write(*,702) ebox_ymin, ebox_ymax, iproc
!     write(*,703) ebox_zmin, ebox_zmax, iproc
 701  format(2x,"before: write: ebox-x ",f10.5,3x,f10.5, 5x,i2)
 702  format(2x,"before: write: ebox-y ",f10.5,3x,f10.5, 5x,i2)
 703  format(2x,"before: write: ebox-z ",f10.5,3x,f10.5, 5x,i2)
      if(iproc.eq.0) then
         open(unit=65,form='unformatted',status='unknown')
         write(65)ebox_xmin,ebox_xmax,ebox_ymin,ebox_ymax,
     &            ebox_zmin,ebox_zmax,dt_in,dtold,
     &            rl,vxl,vyl,vzl,pl,bxl,byl,bzl
         close(65)
      endif
!     write(*,601) ebox_xmin, ebox_xmax, iproc
!     write(*,602) ebox_ymin, ebox_ymax, iproc
!     write(*,603) ebox_zmin, ebox_zmax, iproc
 601  format(2x,"write: ebox-x ",f10.5,3x,f10.5, 5x,i2)
 602  format(2x,"write: ebox-y ",f10.5,3x,f10.5, 5x,i2)
 603  format(2x,"write: ebox-z ",f10.5,3x,f10.5, 5x,i2)

      
      if(allocated(n_to_left)) deallocate( n_to_left )
      if(allocated(glnblocks)) deallocate( glnblocks )
      if(allocated(recvrequest)) deallocate( recvrequest )
      if(allocated(sendrequest)) deallocate( sendrequest )
      if(allocated(recvstatus)) deallocate( recvstatus )
      if(allocated(sendstatus)) deallocate( sendstatus )

      return
      end subroutine amr_checkpoint_wr_a

!-----------------------------------------------------------------------
      

      subroutine amr_checkpoint_re_a (dt_in,time_in,loop_index_start,nchombo,
     &                                                    iunit1)

! Subroutine to read checkpoint file using AMR package.
! Writes out tree data structure and data stored in blocks
! Currently reads are done serially by processor 0 and data sent to
! other processors.
! Uses UNFORMATTED, DIRECT ACCESS

! Written: K. Olson 7/97


      use physcons
      use paramesh_dimensions
      use physicaldata
      use tree
      use interior_gbc


      use paramesh_interfaces
#ifdef MPI_USED
      use paramesh_mpi_interfaces
#endif /* MPI_USED */


      implicit none

#include "amr_shmem.fh"

#ifdef MPI_USED
      include 'mpif.h'
#endif /* MPI_USED */


      integer, intent(in)  :: iunit1
      integer, intent(out) :: loop_index_start,nchombo
      real,    intent(out) :: dt_in,time_in

      integer nguard0
      parameter(nguard0 = nguard*npgs)


      integer l,ii,icoord

      integer block_no,shmem_n_pes,shmem_my_pe
      integer jproc,i,j,ivar,ix,iy,iz,nprocs,iproc,iunit2

      integer lnblockst,nrec,nrecl,ngid,lb
      integer*4 inrecl
      integer gid(nfaces+1+nchild,maxblocks)
      integer alnblocks,alnblockst
      integer ierr_read,ierr_readt
      integer ierr,mype

! TEMPORARIES WHICH ARE READ

      integer  tot_blocks
      integer  lrefinet(maxblocks),nodetypet(maxblocks)
      integer  which_childt(maxblocks)
      integer  gidt(nfaces+1+nchild,maxblocks)
      integer  bflagst(mflags,maxblocks)
      real  coordt(mdim,maxblocks)
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

      save gid,gidt,coordt,lrefinet,nodetypet,lnblockst
      save which_childt,work_blockt
      save alnblocks,alnblockst,bnd_boxt
      save unkt,facevarxt,facevaryt,facevarzt
      save unk_nt,unk_e_xt,unk_e_yt,unk_e_zt
      save ierr_read,ierr_readt
     
      integer parent_pe,parent_blk
      real    par_coord(3)
      save    par_coord,parent_pe,parent_blk

!     real     :: dtold,dt,time
      real     :: dtold
      common/old_timestep/dtold
      real :: rl,vxl,vyl,vzl,pl,bxl,byl,bzl
      common/left_state/rl,vxl,vyl,vzl,pl,bxl,byl,bzl
      common/tot/tot_blocks


      logical :: l_move_solution
      integer :: iopt
      logical :: lcc,lfc,lec,lnc,lfulltree

#ifdef MPI_USED
      integer,dimension (:),  allocatable :: glnblocks
      integer :: isrc,idest,itag,isize,ierror,position
      integer :: ierrorcode,tag_offset
      integer :: new_loc(2,maxblocks_tr),lnblocks_old
      integer status(MPI_STATUS_SIZE)
      real    :: xmin,xmax,ymin,ymax,zmin,zmax
      real    :: xmin1,xmax1,ymin1,ymax1,zmin1,zmax1

      integer,dimension (:),  allocatable :: recvrequest
      integer,dimension (:),  allocatable :: sendrequest
      integer,dimension (:,:),allocatable :: recvstatus
      integer,dimension (:,:),allocatable :: sendstatus
      integer, parameter :: buf_dim_int = 
     .              maxblocks*( 3+mflags+(nfaces+1+nchild))
      integer, parameter :: buf_dim_real = 
     .              maxblocks*( 3*mdim + 1 + len_block 
     .       + nbndvar*(len_blockfx + len_blockfy + len_blockfz)
     .       + nbndvarc*len_blockn
     .       + nbndvare*(len_blockex + len_blockey + len_blockez) )
      integer, parameter :: buf_dim = buf_dim_real + buf_dim_int
      integer :: no_of_bytes_per_real,no_of_bytes_per_integer
      integer :: buf_dim_bytes


      logical :: lguard,lprolong,lflux,ledge,lrestrict

      real :: CS_buffer(buf_dim), CR_buffer(buf_dim)
#endif /* MPI_USED */


      call shmem_barrier_all()

      nprocs = shmem_n_pes()
      iproc  = shmem_my_pe()
      mype = iproc
      ierr_read = 0
      ierr_readt = 0

      call shmem_barrier_all()


#ifdef MPI_USED
      if(allocated(glnblocks)) deallocate( glnblocks )
      allocate ( glnblocks(0:nprocs-1) )

      call mpi_type_size(MPI_INTEGER,no_of_bytes_per_integer,ierr)
      call mpi_type_size(MPI_DOUBLE_PRECISION
     .                  ,no_of_bytes_per_real   ,ierr)
      buf_dim_bytes = buf_dim_real*no_of_bytes_per_real +
     .                buf_dim_int *no_of_bytes_per_integer

#endif /* MPI_USED */

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


#ifndef MPI_USED

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
         read (iunit1)  time_in 
         write(*,*) 're-start time ',time_in
         read (iunit1)  loop_index_start,nchombo 
         write(*,*) 're-start index,nchombo: ',loop_index_start,nchombo

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

               nrec = nrec + 1

!               read (iunit1,rec=nrec) 
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
               read (iunit1) 
     &              ((((unkt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_c),
     &                 ix = 1+il0*ion_c,1+iu0*ion_c),
     &                 iy = 1+jl0*ion_c,1+ju0*ion_c),
     &                 iz = 1+kl0*ion_c,1+ku0*ion_c)
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

2    call shmem_barrier_all()



#else /* MPI_USED */



      call MPI_Barrier  (MPI_COMM_WORLD,ierror)

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
         read (iunit1)  time_in 
         write(*,*) 're-start time ',time_in
         read (iunit1)  loop_index_start,nchombo 
         write(*,*) 're-start index,nchombo: ',loop_index_start,nchombo

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

               nrec = nrec + 1

!               read (iunit1,rec=nrec) 
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
               read (iunit1) 
     &              ((((unk(ivar,ix,iy,iz,block_no),
     &                 ivar = 1,iv_c),
     &                 ix = 1+il0*ion_c,1+iu0*ion_c),
     &                 iy = 1+jl0*ion_c,1+ju0*ion_c),
     &                 iz = 1+kl0*ion_c,1+ku0*ion_c)
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
            enddo
      endif

      if (iproc .gt. 0) then
!
! Post receive on pe iproc for message from proc 0

            lnblockst = glnblocks(iproc)
            isrc = 0
            idest= iproc
            itag = isrc*nprocs + idest+1
            isize = lnblockst*( 3+mflags+3*mdim+(nfaces+1+nchild)
     .                          + 1 + len_block 
     .        + nbndvar*(len_blockfx + len_blockfy + len_blockfz) 
     .        + nbndvare*(len_blockex + len_blockey + len_blockez) 
     .        + nbndvarc*len_blockn )
#ifdef DEBUG
       write(*,*) 'pe ',iproc,' about to post receive from pe 0 ',
     .                ' tag ',itag
#endif /* DEBUG */
            call Mpi_recv(CR_buffer,isize,MPI_DOUBLE_PRECISION,
     .       isrc,itag,MPI_COMM_WORLD,status,ierr)
#ifdef DEBUG
            write(*,*) 'pe ',iproc,' posted receive from pe 0 ',
     .                ' tag ',itag
#endif /* DEBUG */

            position = 0
            do block_no = 1,lnblockst

! fetch data for this block
            call Mpi_unpack(CR_buffer,buf_dim_bytes,position,
     &        lrefine(block_no),1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
            call Mpi_unpack(CR_buffer,buf_dim_bytes,position,
     &        nodetype(block_no),1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
            call Mpi_unpack(CR_buffer,buf_dim_bytes,position,
     &        which_child(block_no),1,MPI_INTEGER,MPI_COMM_WORLD,
     &             ierr)
            call Mpi_unpack(CR_buffer,buf_dim_bytes,position,
     &        bflags(1,block_no),mflags,MPI_INTEGER,
     &        MPI_COMM_WORLD,ierr)
            call Mpi_unpack(CR_buffer,buf_dim_bytes,position,
     &        coord(1,block_no),mdim,MPI_DOUBLE_PRECISION,
     &        MPI_COMM_WORLD,ierr)
            call Mpi_unpack(CR_buffer,buf_dim_bytes,position,
     &        bnd_box(1,1,block_no),2*mdim,
     &        MPI_DOUBLE_PRECISION,
     &        MPI_COMM_WORLD,ierr)
            call Mpi_unpack(CR_buffer,buf_dim_bytes,position,
     &        work_block(block_no),1,MPI_DOUBLE_PRECISION,
     &        MPI_COMM_WORLD,ierr)
            call Mpi_unpack(CR_buffer,buf_dim_bytes,position,
     &        gid(1,block_no),(nfaces+1+nchild),MPI_INTEGER,
     &        MPI_COMM_WORLD,ierr)

               if(nvar.gt.0)
     &           call Mpi_unpack(CR_buffer,buf_dim_bytes,position,
     &             unk(1,1,1,1,block_no),len_block,
     &             MPI_DOUBLE_PRECISION,
     &             MPI_COMM_WORLD,ierr)
               if(nfacevar.gt.0) then
                 call Mpi_unpack(CR_buffer,buf_dim_bytes,position,
     &             facevarx(1,1,1,1,block_no),nbndvar*len_blockfx,
     &             MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
                 call Mpi_unpack(CR_buffer,buf_dim_bytes,position,
     &             facevary(1,1,1,1,block_no),nbndvar*len_blockfy,
     &             MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
                 call Mpi_unpack(CR_buffer,buf_dim_bytes,position,
     &             facevarz(1,1,1,1,block_no),nbndvar*len_blockfz,
     &             MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
               endif
               if(nvaredge.gt.0) then
                 call Mpi_unpack(CR_buffer,buf_dim_bytes,position,
     &             unk_e_x(1,1,1,1,block_no),nbndvare*len_blockex,
     &             MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
                 call Mpi_unpack(CR_buffer,buf_dim_bytes,position,
     &             unk_e_y(1,1,1,1,block_no),nbndvare*len_blockey,
     &             MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
                 call Mpi_unpack(CR_buffer,buf_dim_bytes,position,
     &             unk_e_z(1,1,1,1,block_no),nbndvare*len_blockez,
     &             MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
               endif
               if(nvarcorn.gt.0) then
                 call Mpi_unpack(CR_buffer,buf_dim_bytes,position,
     &             unk_n(1,1,1,1,block_no),nbndvarc*len_blockn,
     &             MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
               endif
#ifdef DEBUG
               write(*,*) 'unk for blk ',block_no,' put into '
     .         ,'pe ',iproc,unk(1,:,:,1,block_no)
#endif /* DEBUG */
            enddo
      endif

      if(iproc.eq.0) then

! Post sends from pe 0 for messages to all other procs

        do jproc = 1,nprocs-1
            lnblockst = glnblocks(jproc)

            isrc = 0
            idest= jproc
            itag = isrc*nprocs + idest+1
            isize = lnblockst*( 3+mflags+3*mdim+(nfaces+1+nchild)
     .                          + 1 + len_block 
     .        + nbndvar*(len_blockfx + len_blockfy + len_blockfz) 
     .        + nbndvare*(len_blockex + len_blockey + len_blockez) 
     .        + nbndvarc*len_blockn )

            position = 0

            do block_no = 1,lnblockst

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
     &              (bnd_boxt(1,j,block_no),j=1,ndim),
     &              (bnd_boxt(2,j,block_no),j=1,ndim),
     &              work_blockt(block_no)
               read (iunit1) 
     &              ((((unkt(ivar,ix,iy,iz),
     &                 ivar = 1,iv_c),
     &                 ix = 1+il0*ion_c,1+iu0*ion_c),
     &                 iy = 1+jl0*ion_c,1+ju0*ion_c),
     &                 iz = 1+kl0*ion_c,1+ku0*ion_c)
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


               call Mpi_pack(lrefinet(block_no),1,MPI_INTEGER,
     &           CS_buffer,buf_dim_bytes,position,
     &           MPI_COMM_WORLD,ierr)
               call Mpi_pack(nodetypet(block_no),1,MPI_INTEGER,
     &           CS_buffer,buf_dim_bytes,position,
     &           MPI_COMM_WORLD,ierr)
               call Mpi_pack(which_childt(block_no),1,MPI_INTEGER,
     &           CS_buffer,buf_dim_bytes,position,
     &           MPI_COMM_WORLD,ierr)
               call Mpi_pack(bflagst(1,block_no),mflags,MPI_INTEGER,
     &           CS_buffer,buf_dim_bytes,position,
     &           MPI_COMM_WORLD,ierr)
               call Mpi_pack(coordt(1,block_no),mdim,
     &           MPI_DOUBLE_PRECISION,
     &           CS_buffer,buf_dim_bytes,position,
     &           MPI_COMM_WORLD,ierr)
               call Mpi_pack(bnd_boxt(1,1,block_no),2*mdim,
     &           MPI_DOUBLE_PRECISION,
     &           CS_buffer,buf_dim_bytes,position,
     &           MPI_COMM_WORLD,ierr)
               call Mpi_pack(work_blockt(block_no),1,
     &           MPI_DOUBLE_PRECISION,
     &           CS_buffer,buf_dim_bytes,position,
     &           MPI_COMM_WORLD,ierr)
               call Mpi_pack(gidt(1,block_no),(nfaces+1+nchild),
     &           MPI_INTEGER,CS_buffer,buf_dim_bytes,position,
     &           MPI_COMM_WORLD,ierr)


              if(nvar.gt.0)
     &          call Mpi_pack(unkt(1,1,1,1),
     &            len_block,MPI_DOUBLE_PRECISION,CS_buffer,
     &            buf_dim_bytes,position,MPI_COMM_WORLD,ierr)

              if(nfacevar.gt.0) then
                call Mpi_pack(facevarxt(1,1,1,1),
     &            len_blockfx*nbndvar,
     &            MPI_DOUBLE_PRECISION,CS_buffer,
     &            buf_dim_bytes,position,MPI_COMM_WORLD,ierr)
                call Mpi_pack(facevaryt(1,1,1,1),
     &            len_blockfy*nbndvar,
     &            MPI_DOUBLE_PRECISION,CS_buffer,
     &            buf_dim_bytes,position,MPI_COMM_WORLD,ierr)
                call Mpi_pack(facevarzt(1,1,1,1),
     &            len_blockfz*nbndvar,
     &            MPI_DOUBLE_PRECISION,CS_buffer,
     &            buf_dim_bytes,position,MPI_COMM_WORLD,ierr)
               endif
              if(nvaredge.gt.0) then
                call Mpi_pack(unk_e_xt(1,1,1,1),
     &            len_blockex*nbndvare,
     &            MPI_DOUBLE_PRECISION,CS_buffer,
     &            buf_dim_bytes,position,MPI_COMM_WORLD,ierr)
                call Mpi_pack(unk_e_yt(1,1,1,1),
     &            len_blockey*nbndvare,
     &            MPI_DOUBLE_PRECISION,CS_buffer,
     &            buf_dim_bytes,position,MPI_COMM_WORLD,ierr)
                call Mpi_pack(unk_e_zt(1,1,1,1),
     &            len_blockez*nbndvare,
     &            MPI_DOUBLE_PRECISION,CS_buffer,
     &            buf_dim_bytes,position,MPI_COMM_WORLD,ierr)
               endif
              if(nvarcorn.gt.0) then
                call Mpi_pack(unk_nt(1,1,1,1),
     &            len_blockn*nbndvarc,
     &            MPI_DOUBLE_PRECISION,CS_buffer,
     &            buf_dim_bytes,position,MPI_COMM_WORLD,ierr)
               endif


            enddo

            call Mpi_send(CS_buffer,isize,MPI_DOUBLE_PRECISION,
     .       idest,itag,MPI_COMM_WORLD,ierr)

        enddo

      endif

2     call MPI_Barrier  (MPI_COMM_WORLD,ierror)



      do block_no = 1,lnblocks
        bsize(:,block_no) = bnd_box(2,:,block_no)-
     $                     bnd_box(1,:,block_no)
      enddo

      

#endif /* MPI_USED */

#ifdef DEBUG
      write(*,*) 'pe ',iproc,' starting tree rebuild '
#endif /* DEBUG */
 
! all processors fetch error code from proc 0

#ifdef MPI_USED
      call mpi_bcast(ierr_read,1,MPI_INTEGER,0,
     &                 MPI_COMM_WORLD,ierr)
      ierr_readt = ierr_read
      if (ierr_readt.eq.1) 
     &   call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
#else
      call SHMEM_INTEGER_GET (ierr_readt,ierr_read,1,0)
      if (ierr_readt.eq.1) call abort
#endif /* MPI_USED */


! COMPUTE TREE DATA FROM gid

#ifdef MPI_USED
      call mpi_bcast(alnblocks,1,MPI_INTEGER,0,
     &                 MPI_COMM_WORLD,ierr)
      alnblockst = alnblocks
#else
      call shmem_barrier_all()
      call SHMEM_INTEGER_GET (alnblockst,alnblocks,1,0)
      alnblocks = alnblockst
#endif /* MPI_USED */

      call shmem_barrier_all()

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


      call shmem_barrier_all()


! Now reorder blocks such that they are better balanced
! NOTE: this assumes that the total number of blocks is > nprocs

#ifdef MPI_USED
      lnblocks_old = lnblocks
      l_move_solution=.true.
!     call amr_morton_order (lnblocks_old,new_loc,nprocs,iproc)
      call amr_morton_order (lnblocks_old,nprocs,iproc,l_move_solution)
#else /* MPI_USED */
      call amr_morton_order (nprocs,lnblocks)
#endif /* MPI_USED */

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

      call shmem_barrier_all()


#ifdef MPI_USED
! Find the coordinate ranges
         call mpi_amr_global_domain_limits

      if(mype.eq.0) then
      write(*,*) ' grid_xmin ',grid_xmin
      write(*,*) ' grid_xmax ',grid_xmax
      write(*,*) ' grid_ymin ',grid_ymin
      write(*,*) ' grid_ymax ',grid_ymax
      write(*,*) ' grid_zmin ',grid_zmin
      write(*,*) ' grid_zmax ',grid_zmax
      endif


!
! Compute and save morton number range for each processor
      call mpi_amr_morton_limits(mype)
      write(*,*) 'pe ',mype,' after mpi_amr_morton_limits'

!
! Set up surrounding blocks of all local blocks (must not precede
! setting of grid_xmin,... etc)
      tag_offset = 100
      call mpi_morton_bnd(mype,nprocs,tag_offset)
      write(*,*) 'pe ',mype,' after mpi_morton_bnd'
      call mpi_amr_gsurr_blks(mype,nprocs)
      write(*,*) 'pe ',mype,' after mpi_amr_gsurr_blks'
      tag_offset = 100
      call mpi_morton_bnd_prolong1(mype,nprocs,tag_offset)
      write(*,*) 'pe ',mype,' after mpi_morton_bnd_prol'
      tag_offset = 100
      call mpi_morton_bnd_fluxcon(mype,nprocs,tag_offset)
      write(*,*) 'pe ',mype,' after mpi_morton_bnd_flux'
      tag_offset = 100
      lec = .false.
      lnc = .false.
      if(nvaredge.gt.0) lec = .true.
      if(nvarcorn.gt.0) lnc = .true.
      lfulltree=.false.
      call mpi_morton_bnd_restrict(mype,nprocs,lfulltree,lec,lnc,tag_offset)
      write(*,*) 'pe ',mype,' after mpi_morton_bnd_restrict'


      call mpi_setup(mype,nprocs)
#endif /* MPI_USED */




! Now make sure guardcell information is up to date

#ifndef NO_PERMANENT_GUARDCELLS
#ifdef MPI_USED
      lcc = .true.
      lfc = .true.
      lec = .true.
      lnc = .true.
      iopt = 1
      tag_offset = 100
      lguard    = .true.
      lprolong  = .false.
      lflux     = .false.
      ledge     = .false.
      lrestrict = .false.
      lfulltree=.false.
      call mpi_amr_comm_setup(mype,nprocs,
     .                        lguard,lprolong,lflux,ledge,lrestrict,
     .                        lfulltree,
     .                        iopt,lcc,lfc,lec,lnc,tag_offset)

#endif /* MPI_USED */

      call amr_guardcell(iproc,1,nguard)

      call shmem_barrier_all()
#endif

      if(allocated(glnblocks)) deallocate( glnblocks )

      call shmem_barrier_all()

      if(iproc.eq.0) then
         open(unit=65,form='unformatted',status='unknown')
         read(65)ebox_xmin,ebox_xmax,ebox_ymin,ebox_ymax,
     &            ebox_zmin,ebox_zmax,dt_in,dtold,
     &            rl,vxl,vyl,vzl,pl,bxl,byl,bzl
         close(65)
      endif
      
      write(*,601) ebox_xmin, ebox_xmax, mype
      write(*,602) ebox_ymin, ebox_ymax, mype
      write(*,603) ebox_zmin, ebox_zmax, mype
 601  format(2x,"read: ebox-x ",f10.5,3x,f10.5, 5x,i2)
 602  format(2x,"read: ebox-y ",f10.5,3x,f10.5, 5x,i2)
 603  format(2x,"read: ebox-z ",f10.5,3x,f10.5, 5x,i2)

      call shmem_barrier_all()

#ifdef MPI_USED
      call mpi_bcast(loop_index_start,1,MPI_INTEGER,
     .                                    0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(nchombo,1,MPI_INTEGER,
     .                                    0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(ebox_xmin,1,MPI_DOUBLE_PRECISION,
     .                                    0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(ebox_xmax,1,MPI_DOUBLE_PRECISION,
     .                                    0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(ebox_ymin,1,MPI_DOUBLE_PRECISION,
     .                                    0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(ebox_ymax,1,MPI_DOUBLE_PRECISION,
     .                                    0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(ebox_zmin,1,MPI_DOUBLE_PRECISION,
     .                                    0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(ebox_zmax,1,MPI_DOUBLE_PRECISION,
     .                                    0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(time_in     ,1,MPI_DOUBLE_PRECISION,
     .                                    0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(dt_in       ,1,MPI_DOUBLE_PRECISION,
     .                                    0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(dtold    ,1,MPI_DOUBLE_PRECISION,
     .                                    0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(rl       ,1,MPI_DOUBLE_PRECISION,
     .                                    0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(vxl      ,1,MPI_DOUBLE_PRECISION,
     .                                    0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(vyl      ,1,MPI_DOUBLE_PRECISION,
     .                                    0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(vzl      ,1,MPI_DOUBLE_PRECISION,
     .                                    0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(pl       ,1,MPI_DOUBLE_PRECISION,
     .                                    0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(bxl      ,1,MPI_DOUBLE_PRECISION,
     .                                    0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(byl      ,1,MPI_DOUBLE_PRECISION,
     .                                    0,MPI_COMM_WORLD,ierr)
      call mpi_bcast(bzl      ,1,MPI_DOUBLE_PRECISION,
     .                                    0,MPI_COMM_WORLD,ierr)
      call MPI_Barrier  (MPI_COMM_WORLD,ierror)
      write(*,*) 'leaving MPI_USED:3'
#else
      call shmem_integer_get(loop_index_start,loop_index_start,1,0)
      call shmem_integer_get(nchombo,nchombo,1,0)
      call shmem_real_get(ebox_xmin,ebox_xmin,1,0)
      call shmem_real_get(ebox_xmax,ebox_xmax,1,0)
      call shmem_real_get(ebox_ymin,ebox_ymin,1,0)
      call shmem_real_get(ebox_ymax,ebox_ymax,1,0)
      call shmem_real_get(ebox_zmin,ebox_zmin,1,0)
      call shmem_real_get(ebox_zmax,ebox_zmax,1,0)
      call shmem_real_get(time_in,time_in,1,0)
      call shmem_real_get(dt_in,dt_in,1,0)
      call shmem_real_get(dtold,dtold,1,0)
      call shmem_real_get(rl,rl,1,0)
      call shmem_real_get(vxl,vxl,1,0)
      call shmem_real_get(vyl,vyl,1,0)
      call shmem_real_get(vzl,vzl,1,0)
      call shmem_real_get(pl,pl,1,0)
      call shmem_real_get(bxl,bxl,1,0)
      call shmem_real_get(byl,byl,1,0)
      call shmem_real_get(bzl,bzl,1,0)
      call shmem_barrier_all()
#endif /* MPI_USED */
      write(*,*) 'leaving subroutine amr_checkpoint_re_a'

      return
      end subroutine amr_checkpoint_re_a
