#include "paramesh_preprocessor.fh"
      subroutine amr_symmetry_check(mype)
      use physcons
      use paramesh_dimensions
      use physicaldata
      use workspace
      use athena_specific_parameters
! include file defining the tree
      use tree


      use paramesh_interfaces
#ifdef MPI_USED
      use paramesh_mpi_interfaces
#endif /* MPI_USED */


      implicit real(a-h,o-z)



#ifdef MPI_USED

      include "mpif.h"
      integer :: tag_offset,max_blks_sent
      integer status(MPI_STATUS_SIZE)
#endif /* MPI_USED */


#include "amr_shmem.fh"

! header file declaring arrays used to test symmetry properties.
	integer nvarsym
	parameter(nvarsym=ne)
        integer lengthsym
        parameter(lengthsym=nvarsym*iu_bnd*ju_bnd*ku_bnd)

        real cutoff
        parameter(cutoff = 1.e-20)

        integer, parameter ::  nguard0 = nguard*npgs
        integer :: errcode
!-------------------------

      integer lb1,lb2,pe1,pe2,iv

! local arrays
      real symm1(nvarsym,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
      real symm2(nvarsym,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
      real fv(nvar)

        real            recvx(nbndvar,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,
     .                                          kl_bnd:ku_bnd)
        real            recvy(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd+1,
     .                                          kl_bnd:ku_bnd)
        real            recvz(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,
     .                                          kl_bnd:ku_bnd+k3d)
        real, parameter ::  error = 1.0e-7
        save recvx,recvy,recvz

      save symm1,symm2,lb2,pe2
      integer idest, itag, num, ierr

!-------------------------



#ifdef MPI_USED
        call mpi_barrier (MPI_COMM_WORLD, errcode)
#else
        call shmem_barrier_all()
#endif /*MPI_USED*/

! Set symmetries for each variable
      fv(:) =  1.
      fv(3) = -1.        ! asymmetry in y mom.
      fv(6) = -1.        ! asymmetry in Bx.

! first test block
      pe1 = 0
      lb1 = 10
      lb1 = 9
      lb1 = 6
      lb1 = 229
c      lb1 = 99  ! parent

! second test block
      pe2 = 0
      lb2 = 37
      lb2 = 27
      lb2 = 445
c      lb2 = 177 ! parent

#ifdef MPI_USED
!send data
      if(mype.eq.pe1) then
          idest=0
          itag=5000+10*(pe1-1)+1; num=lengthsym;
          call MPI_SEND(unk(1,1,1,1,lb1),num,
     &    MPI_DOUBLE_PRECISION,idest,
     &    itag,MPI_COMM_WORLD,ierr)
      endif
      if(mype.eq.pe2) then
          idest=0
          itag=5000+10*(pe2-1)+1; num=lengthsym;
          call MPI_SEND(unk(1,1,1,1,lb2),num,
     &    MPI_DOUBLE_PRECISION,idest,
     &    itag,MPI_COMM_WORLD,ierr)

          itag=5000+10*(pe2-1)+2; num=len_blockfx*nbndvar;
          call MPI_SEND(facevarx(:,:,:,:,lb2),num,
     &    MPI_DOUBLE_PRECISION,idest,
     &    itag,MPI_COMM_WORLD,ierr)
          
          itag=5000+10*(pe2-1)+3; num=len_blockfy*nbndvar;
          call MPI_SEND(facevary(:,:,:,:,lb2),num,
     &    MPI_DOUBLE_PRECISION,idest,
     &    itag,MPI_COMM_WORLD,ierr)

          itag=5000+10*(pe2-1)+4; num=len_blockfz*nbndvar;
          call MPI_SEND(facevarz(:,:,:,:,lb2),num,
     &    MPI_DOUBLE_PRECISION,idest,
     &    itag,MPI_COMM_WORLD,ierr)
      endif
#endif /*MPI_USED*/

      if(mype.eq.0) then

#ifdef SHMEM_USED
! Copy complete remote blocks into a buffer block called recv.
      call shmem_real_get(symm1(1,1,1,1),unk(1,1,1,1,lb1),
     &       lengthsym,pe1) 
      call shmem_real_get(symm2(1,1,1,1),unk(1,1,1,1,lb2),
     &       lengthsym,pe2) 

! Copy complete remote block into a buffer block.
        call shmem_real_get(recvx,facevarx(:,:,:,:,lb2),
     &                          len_blockfx*nbndvar,pe2) 
        call shmem_real_get(recvy,facevary(:,:,:,:,lb2),
     &                          len_blockfy*nbndvar,pe2) 
        call shmem_real_get(recvz,facevarz(:,:,:,:,lb2),
     &                          len_blockfz*nbndvar,pe2) 
#endif /*SHMEM_USED*/
#ifdef MPI_USED
!receive data
          icomm = MPI_COMM_WORLD
          isource=pe1
          itag=5000+10*(pe1-1)+1; num=lengthsym;
          call MPI_RECV(symm1(1,1,1,1),num,
     &    MPI_DOUBLE_PRECISION,isource,
     &                itag,icomm,status,ierr)

          isource=pe2
          itag=5000+10*(pe2-1)+1; num=lengthsym;
          call MPI_RECV(symm2(1,1,1,1),num,
     &    MPI_DOUBLE_PRECISION,isource,
     &                itag,icomm,status,ierr)

          itag=5000+10*(pe2-1)+2; num=len_blockfx*nbndvar;
          call MPI_RECV(recvx,num,
     &    MPI_DOUBLE_PRECISION,isource,
     &                itag,icomm,status,ierr)

          itag=5000+10*(pe2-1)+3; num=len_blockfy*nbndvar;
          call MPI_RECV(recvy,num,
     &    MPI_DOUBLE_PRECISION,isource,
     &                itag,icomm,status,ierr)

          itag=5000+10*(pe2-1)+4; num=len_blockfz*nbndvar;
          call MPI_RECV(recvz,num,
     &    MPI_DOUBLE_PRECISION,isource,
     &                itag,icomm,status,ierr)
#endif /*MPI_USED*/

      do k=1,nzb+2*nguard0*k3d
      do j=1,nyb+2*nguard0*k2d
      do i=1,nxb+2*nguard0
!
! symmetry through origin
      kk = nzb+(2*nguard0+1 - k)*k3d
      jj = nyb+(2*nguard0+1 - j)*k2d
      ii = nxb+2*nguard0+1 - i
!
! symmetry about y mid point
      kk = k
      jj = nyb+(2*nguard0+1 - j)*k2d
      ii = i

      aa = symm1(1,i ,j ,k )
      bb = symm2(1,ii,jj,kk)
      cc = abs(aa-bb)/(abs(aa+bb)+1.e-30)
      if(cc.gt.error) write(*,*) 
     .      'asymm : symm1(1)',i,j,k,symm1(1,i,j,k),ii,jj,kk,
     .       symm2(1,ii,jj,kk)

!mom
      do iv = 2,4
      aa = symm1(iv,i,j,k)
      bb = fv(iv)*symm2(iv,ii,jj,kk)
      cc = abs(aa-bb)/(abs(aa)+abs(bb)+1.e-30)
      if((cc.gt.error).and.abs(aa).gt.cutoff) write(*,*) 
     .      'asymm : symm1(',iv,')',i,j,k,symm1(iv,i,j,k),
     .     ii,jj,kk, symm2(iv,ii,jj,kk)
      enddo
!B
      do iv = 6,8
      aa = symm1(iv,i,j,k)
      bb = fv(iv)*symm2(iv,ii,jj,kk)
      cc = abs(aa-bb)/(abs(aa)+abs(bb)+1.e-30)
      if((cc.gt.error).and.abs(aa).gt.cutoff) write(*,*) 
     .      'asymm : symm1(',iv,')',i,j,k,symm1(iv,i,j,k),
     .      ii,jj,kk,symm2(iv,ii,jj,kk)
      enddo
      enddo
      enddo
      enddo


      do k=1,nzb+2*nguard0*k3d
      do j=1,nyb+2*nguard0*k2d
      do i=1,nxb+2*nguard0+1

      kk = k
      jj = nyb+2*nguard0+1 - j
      ii = i
      aa = facevarx(1,i,j,k,lb1)
      bb = fv(6)*recvx(1,i,jj,kk)
      cc = abs(aa-bb)/(abs(aa)+abs(bb)+1.e-30)
      if((cc.gt.error).and.abs(aa).gt.cutoff) write(*,*) 
     .       'asymm : facevarx(1)',i,j,k,facevarx(1,i,j,k,lb1),
     .       ii,jj,kk,recvx(1,ii,jj,kk)
      enddo
      enddo
      enddo


      do k=1,nzb+2*nguard0*k3d
      do j=1,nyb+(2*nguard0+1)*k2d
      do i=1,nxb+2*nguard0
      kk = k
      jj = nyb+2*nguard0+1 - j + 1
      ii = i
      aa = facevary(1,i,j,k,lb1)
      bb = fv(7)*recvy(1,i,jj,kk)
      cc = abs(aa-bb)/(abs(aa)+abs(bb)+1.e-30)
      if((cc.gt.error).and.abs(aa).gt.cutoff) write(*,*) 
     .       'asymm : facevary(1)',i,j,k,facevary(1,i,j,k,lb1),
     .       i,jj,kk,recvy(1,ii,jj,kk)
      enddo
      enddo
      enddo

#if N_DIM == 3
!      do k=1+nguard*k3d,nzb+(nguard+1)*k3d
!      do j=1+nguard*k2d,nyb+nguard*k2d
!      do i=1+nguard,nxb+nguard+1

      do k=1,nzb+(2*nguard0+1)*k3d
      do j=1,nyb+2*nguard0*k2d
      do i=1,nxb+2*nguard0
      kk = k
      jj = nyb+2*nguard0+1 - j
      ii = i
      aa = facevarz(1,i,j,k,lb1)
      bb = fv(8)*recvz(1,i,jj,kk)
      cc = abs(aa-bb)/(abs(aa)+abs(bb)+1.e-30)
      if((cc.gt.error).and.abs(aa).gt.cutoff) write(*,*) 
     .       'asymm : facevarz(1)',i,j,k,facevarz(1,i,j,k,lb1),
     .       ii,jj,kk,recvz(1,ii,jj,kk)
      enddo
      enddo
      enddo
#endif

      endif
#ifdef MPI_USED
        call mpi_barrier (MPI_COMM_WORLD, errcode)
#else
        call shmem_barrier_all()
#endif /*MPI_USED*/

      return
      end
