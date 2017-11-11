!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"


      subroutine guardcell_test(mype)





!------------------------------------------------------------------------
!
! This routine tests the guardcell filling routine for omissions
! (but not necessarily inaccuracies). It does it by initializing
! unk(1,...) with a constant, gunk(...) with a different constant,
! filling the block interiors of both with the block address, 
! then calling amr_guardcell and comparing guardcell results.
! Any inconsistencies flag an error.
!
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

      use paramesh_interfaces, only : amr_1blk_guardcell_reset,
     &                                amr_guardcell

      implicit real(a-h,o-z)

#include "amr_shmem.fh"

      integer, intent(in)    ::  mype

      integer :: nprocs


      integer nguard0
#ifdef NO_PERMANENT_GUARDCELLS
      parameter(nguard0 = nguard*npgs)
#endif
#ifndef NO_PERMANENT_GUARDCELLS
      parameter(nguard0 = nguard)
#endif

#ifdef TIMINGS
#include "timer.fh"
#endif

!-------------------------
! common block storing the solution for cell-face-centered quantities.
#ifdef NO_PERMANENT_GUARDCELLS
        common/solution_gcc/
     .           gunk(il_bnd:iu_bnd,jl_bnd:ju_bnd,
     .                          kl_bnd:ku_bnd)
        real gunk

        common/solution_gfc/
     .           gfacevarx(il_bnd:iu_bnd+1,jl_bnd:ju_bnd,
     .                          kl_bnd:ku_bnd)
     .          ,gfacevary(il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,
     .                          kl_bnd:ku_bnd)
     .          ,gfacevarz(il_bnd:iu_bnd,jl_bnd:ju_bnd,
     .                          kl_bnd:ku_bnd+k3d)
        real gfacevarx,gfacevary,gfacevarz

        common/solution_gec/
     .           gunk_e_x(il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,
     .                          kl_bnd:ku_bnd+k3d)
     .          ,gunk_e_y(il_bnd:iu_bnd+1,jl_bnd:ju_bnd,
     .                          kl_bnd:ku_bnd+k3d)
     .          ,gunk_e_z(il_bnd:iu_bnd+1,jl_bnd:ju_bnd+k2d,
     .                          kl_bnd:ku_bnd)
        real :: gunk_e_x,gunk_e_y,gunk_e_z

        common/solution_gnc/
     .           gunk_n(il_bnd:iu_bnd+1,jl_bnd:ju_bnd+k2d,
     .                          kl_bnd:ku_bnd+k3d)
        real :: gunk_n
#else
        common/solution_gcc/
     .           gunk(il_bnd:iu_bnd,jl_bnd:ju_bnd,
     .                          kl_bnd:ku_bnd,maxblocks)
        real gunk

        common/solution_gfc/
     .           gfacevarx(il_bnd:iu_bnd+1,jl_bnd:ju_bnd,
     .                          kl_bnd:ku_bnd,maxblocksf)
     .          ,gfacevary(il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,
     .                          kl_bnd:ku_bnd,maxblocksf)
     .          ,gfacevarz(il_bnd:iu_bnd,jl_bnd:ju_bnd,
     .                          kl_bnd:ku_bnd+k3d,maxblocksf)
        real gfacevarx,gfacevary,gfacevarz

        common/solution_gec/
     .           gunk_e_x(il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,
     .                          kl_bnd:ku_bnd+k3d,maxblocksue)
     .          ,gunk_e_y(il_bnd:iu_bnd+1,jl_bnd:ju_bnd,
     .                          kl_bnd:ku_bnd+k3d,maxblocksue)
     .          ,gunk_e_z(il_bnd:iu_bnd+1,jl_bnd:ju_bnd+k2d,
     .                          kl_bnd:ku_bnd,maxblocksue)
        real :: gunk_e_x,gunk_e_y,gunk_e_z

        common/solution_gnc/
     .           gunk_n(il_bnd:iu_bnd+1,jl_bnd:ju_bnd+k2d,
     .                          kl_bnd:ku_bnd+k3d,maxblocksn)
        real :: gunk_n
#endif

      integer remote_pe,remote_block
      logical ldiag

! local arrays

      integer,external ::  shmem_n_pes
      integer :: surrblks(3,3,3,3)

!-------------------------

      nprocs = shmem_n_pes()


      if(nvar.lt.1.and.nfacevar.lt.1) then
        write(*,*) 'guardcell_test : ignored. nvar or nfacevar is ',
     .             'too small ! '
        return
      endif

!-------------------------

      call shmem_barrier_all()

      call amr_1blk_guardcell_reset

      call shmem_barrier_all()

#ifndef NO_PERMANENT_GUARDCELLS

! Loop over blocks, setting block interior data to the local block number,
! and the block guardcells to selected constants, differing for each data
! word.
      if(lnblocks.gt.0) then
      do lb = 1,lnblocks

        if(nvar.ge.1) then
        unk(1,:,:,:,lb) = -100.
        unk(1,1+nguard0:nxb+nguard0,1+nguard0*k2d:nyb+nguard0*k2d,
     .                1+nguard0*k3d:nzb+nguard0*k3d,lb) = real(lb)
        endif

        if(nfacevar.ge.1) then
        facevarx(1,:,:,:,lb) = -100.
        facevarx(1,1+nguard0:nxb+nguard0+1,
     .               1+nguard0*k2d:nyb+nguard0*k2d,
     .               1+nguard0*k3d:nzb+nguard0*k3d,lb) = real(lb)
        if(ndim.ge.2) then
        facevary(1,:,:,:,lb) = -100.
        facevary(1,1+nguard0:nxb+nguard0,
     .               1+nguard0*k2d:nyb+(nguard0+1)*k2d,
     .               1+nguard0*k3d:nzb+nguard0*k3d,lb) = real(lb)
        endif
        if(ndim.eq.3) then
        facevarz(1,:,:,:,lb) = -100.
        facevarz(1,1+nguard0:nxb+nguard0,
     .               1+nguard0*k2d:nyb+nguard0*k2d,
     .               1+nguard0*k3d:nzb+(nguard0+1)*k3d,lb) = real(lb)
        endif
        endif

        if(nvaredge.ge.1) then
        unk_e_x(1,:,:,:,lb) = -100.
        unk_e_y(1,:,:,:,lb) = -100.
        unk_e_z(1,:,:,:,lb) = -100.
        unk_e_x(1,1+nguard0:nxb+nguard0,
     .            1+nguard0*k2d:nyb+(nguard0+k2d)*k2d,
     .            1+nguard0*k3d:nzb+(nguard0+k3d)*k3d,lb) = real(lb)
        unk_e_y(1,1+nguard0:nxb+nguard0+1,
     .            1+nguard0*k2d:nyb+nguard0*k2d,
     .            1+nguard0*k3d:nzb+(nguard0+k3d)*k3d,lb) = real(lb)
        unk_e_z(1,1+nguard0:nxb+nguard0+1,
     .            1+nguard0*k2d:nyb+(nguard0+k2d)*k2d,
     .            1+nguard0*k3d:nzb+nguard0*k3d,lb) = real(lb)
        endif

        if(nvarcorn.ge.1) then
        unk_n(1,:,:,:,lb) = -100.
        unk_n(1,1+nguard0:nxb+nguard0+1,
     .          1+nguard0*k2d:nyb+(nguard0+k2d)*k2d,
     .          1+nguard0*k3d:nzb+(nguard0+k3d)*k3d,lb) = real(lb)
        endif

      enddo                             ! end of loop over blocks
      endif

      call shmem_barrier_all()
!-------------------------

      iopt = 1
      nlayers = nguard0
      call amr_guardcell(mype,iopt,nlayers)

      call shmem_barrier_all()
!-------------------------

      if(nvar.ge.1) then
        gunk(:,:,:,:) = unk(1,:,:,:,:)
      endif

      if(nfacevar.ge.1) then
        gfacevarx(:,:,:,:) = facevarx(1,:,:,:,:)
        if(ndim.ge.2) gfacevary(:,:,:,:) = facevary(1,:,:,:,:)
        if(ndim.eq.3) gfacevarz(:,:,:,:) = facevarz(1,:,:,:,:)
      endif

      if(nvaredge.ge.1) then
        gunk_e_x(:,:,:,:) = unk_e_x(1,:,:,:,:)
        gunk_e_y(:,:,:,:) = unk_e_y(1,:,:,:,:)
        gunk_e_z(:,:,:,:) = unk_e_z(1,:,:,:,:)
      endif

      if(nvarcorn.ge.1) then
        gunk_n(:,:,:,:) = unk_n(1,:,:,:,:)
      endif


      if(lnblocks.gt.0) then
      do lb = 1,lnblocks
      if(nodetype(lb).eq.1) then

        if(nvar.ge.1) then
        unk(1,:,:,:,lb) = -50.
        unk(1,1+nguard0:nxb+nguard0,1+nguard0*k2d:nyb+nguard0*k2d,
     .                1+nguard0*k3d:nzb+nguard0*k3d,lb) = real(lb)
        endif

        if(nfacevar.ge.1) then
        facevarx(1,:,:,:,lb) = -50.
        facevarx(1,1+nguard0:nxb+nguard0+1,
     .               1+nguard0*k2d:nyb+nguard0*k2d,
     .               1+nguard0*k3d:nzb+nguard0*k3d,lb) = real(lb)
        if(ndim.ge.2) then
        facevary(1,:,:,:,lb) = -50.
        facevary(1,1+nguard0:nxb+nguard0,
     .               1+nguard0*k2d:nyb+(nguard0+1)*k2d,
     .               1+nguard0*k3d:nzb+nguard0*k3d,lb) = real(lb)
        endif
        if(ndim.eq.3) then
        facevarz(1,:,:,:,lb) = -50.
        facevarz(1,1+nguard0:nxb+nguard0,
     .               1+nguard0*k2d:nyb+nguard0*k2d,
     .               1+nguard0*k3d:nzb+(nguard0+1)*k3d,lb) = real(lb)
        endif
        endif

        if(nvaredge.ge.1) then
        unk_e_x(1,:,:,:,lb) = -50.
        unk_e_y(1,:,:,:,lb) = -50.
        unk_e_z(1,:,:,:,lb) = -50.
        unk_e_x(1,1+nguard0:nxb+nguard0,
     .            1+nguard0*k2d:nyb+(nguard0+k2d)*k2d,
     .            1+nguard0*k3d:nzb+(nguard0+k3d)*k3d,lb) = real(lb)
        unk_e_y(1,1+nguard0:nxb+nguard0+1,
     .            1+nguard0*k2d:nyb+nguard0*k2d,
     .            1+nguard0*k3d:nzb+(nguard0+k3d)*k3d,lb) = real(lb)
        unk_e_z(1,1+nguard0:nxb+nguard0+1,
     .            1+nguard0*k2d:nyb+(nguard0+k2d)*k2d,
     .            1+nguard0*k3d:nzb+nguard0*k3d,lb) = real(lb)
        endif

        if(nvarcorn.ge.1) then
        unk_n(1,:,:,:,lb) = -50.
        unk_n(1,1+nguard0:nxb+nguard0+1,
     .          1+nguard0*k2d:nyb+(nguard0+k2d)*k2d,
     .          1+nguard0*k3d:nzb+(nguard0+k3d)*k3d,lb) = real(lb)
        endif

      endif
      enddo                             ! end of loop over blocks
      endif

!-------------------------

      call amr_guardcell(mype,iopt,nlayers)

      call shmem_barrier_all()
!-------------------------

      if(lnblocks.gt.0) then
      do lb = 1,lnblocks
      if(nodetype(lb).eq.1) then


        minneigh = 1
        do if = 1,nfaces
        minneigh = min(minneigh,neigh(1,if,lb))
        enddo

        iprint = 0
        if(minneigh.gt.-20) then


!
! Test cell centered data
        if(nvar.ge.1) then

        do k=1,nzb+2*nguard0*k3d
        do j=1,nyb+2*nguard0*k2d
        do i=1,nxb+2*nguard0
        if(unk(1,i,j,k,lb).ne.gunk(i,j,k,lb)) then
        write(*,*) 'guard_test - unk : pe/blk ',mype,lb,i,j,k,
     .       unk(1,i,j,k,lb),gunk(i,j,k,lb)
        iprint = 1
        endif
        enddo
        enddo
        enddo

        endif

        if(nfacevar.ge.1) then
!
! Test facevarx data
        do k=1,nzb+2*nguard0*k3d
        do j=1,nyb+2*nguard0*k2d
        do i=1,nxb+2*nguard0+1
        if(facevarx(1,i,j,k,lb).ne.gfacevarx(i,j,k,lb)) then
        write(*,*) 'guard_test - fx : pe/blk ',mype,lb,i,j,k,
     .      facevarx(1,i,j,k,lb),gfacevarx(i,j,k,lb)
        iprint = 1
        endif
        enddo
        enddo
        enddo

!
! Test facevary data
        if(ndim.ge.2) then
        do k=1,nzb+2*nguard0*k3d
        do j=1,nyb+2*nguard0*k2d+k2d
        do i=1,nxb+2*nguard0
        if(facevary(1,i,j,k,lb).ne.gfacevary(i,j,k,lb)) then
        write(*,*) 'guard_test - fy : pe/blk ',mype,lb,i,j,k,
     .      facevary(1,i,j,k,lb),gfacevary(i,j,k,lb)
        iprint = 1
        endif
        enddo
        enddo
        enddo
        endif

!
! Test facevarx data
        if(ndim.eq.3) then
        do k=1,nzb+2*nguard0*k3d+k3d
        do j=1,nyb+2*nguard0*k2d
        do i=1,nxb+2*nguard0
        if(facevarz(1,i,j,k,lb).ne.gfacevarz(i,j,k,lb)) then
        write(*,*) 'guard_test - fz : pe/blk ',mype,lb,i,j,k,
     .      facevarz(1,i,j,k,lb),gfacevarz(i,j,k,lb)
        iprint = 1
        endif
        enddo
        enddo
        enddo
        endif

        endif                           ! end of nfacevar if test

!
! Test cell edge centered data
        if(nvaredge.ge.1) then

        do k=1,nzb+(2*nguard0+1)*k3d
        do j=1,nyb+(2*nguard0+1)*k2d
        do i=1,nxb+2*nguard0
        if(unk_e_x(1,i,j,k,lb).ne.gunk_e_x(i,j,k,lb)) then
        write(*,*) 'guard_test - unk_e_x : pe/blk ',mype,lb,i,j,k,
     .       unk_e_x(1,i,j,k,lb),gunk_e_x(i,j,k,lb)
        iprint = 1
        endif
        enddo
        enddo
        enddo

        do k=1,nzb+(2*nguard0+1)*k3d
        do j=1,nyb+(2*nguard0)*k2d
        do i=1,nxb+2*nguard0+1
        if(unk_e_y(1,i,j,k,lb).ne.gunk_e_y(i,j,k,lb)) then
        write(*,*) 'guard_test - unk_e_y : pe/blk ',mype,lb,i,j,k,
     .       unk_e_y(1,i,j,k,lb),gunk_e_y(i,j,k,lb)
        iprint = 1
        endif
        enddo
        enddo
        enddo

        do k=1,nzb+(2*nguard0)*k3d
        do j=1,nyb+(2*nguard0+1)*k2d
        do i=1,nxb+2*nguard0+1
        if(unk_e_z(1,i,j,k,lb).ne.gunk_e_z(i,j,k,lb)) then
        write(*,*) 'guard_test - unk_e_z : pe/blk ',mype,lb,i,j,k,
     .       unk_e_z(1,i,j,k,lb),gunk_e_z(i,j,k,lb)
        iprint = 1
        endif
        enddo
        enddo
        enddo

        endif

!
! Test cell corner data
        if(nvarcorn.ge.1) then

        do k=1,nzb+(2*nguard0+1)*k3d
        do j=1,nyb+(2*nguard0+1)*k2d
        do i=1,nxb+2*nguard0+1
        if(unk_n(1,i,j,k,lb).ne.gunk_n(i,j,k,lb)) then
        write(*,*) 'guard_test - unk_n : pe/blk ',mype,lb,i,j,k,
     .       unk_n(1,i,j,k,lb),gunk_n(i,j,k,lb)
        iprint = 1
        endif
        enddo
        enddo
        enddo

        endif

        endif                           ! end of minneigh if test


      endif
      enddo
      endif

      call shmem_barrier_all()
#endif


#ifdef NO_PERMANENT_GUARDCELLS

#endif



      return
      end subroutine guardcell_test
