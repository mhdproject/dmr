!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

      subroutine amr_prolong_face_fun_init




!------------------------------------------------------------------------
!
! This routine computes the values of dx,dy and dz used during the
! interpolation process. These are used inside the prolongation routines
! saving needless repetitive computation at the cost of minimal storage
! space.
!
! This particular prolongation is simple linear interpolation. It can
! only be used for blocks with an even number of grid cells.
!
! Written :     Peter MacNeice          July 1997
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace
      use prolong_arrays

      implicit none

      real :: xc0,yc0,zc0,xc,yc,zc
      integer :: i,j,k,ii,jj,kk,loop,jchild,ioff,joff,koff
      integer :: i1,j1,k1,i1p,j1p,k1p

!------------------------------------


      xc0 = .5
      yc0 = .5
      zc0 = .5


! Initialize the values of dx,dy,dz needed in prolong_face_fun.
      do k=kl_bnd1,ku_bnd1+k3d
       kk = k-nguard+2*nzb
       zc = zc0+real(kk)*.5
       prol_f_dz(k) = mod(zc,1.)
      enddo
      do j=jl_bnd1,ju_bnd1+k2d
       jj = j-nguard+2*nyb
       yc = yc0+real(jj)*.5
       prol_f_dy(j) = mod(yc,1.)
      enddo
      do i=il_bnd1,iu_bnd1+1
       ii = i-nguard+2*nxb
       xc = xc0+real(ii)*.5
       prol_f_dx(i) = mod(xc,1.)
      enddo


! Compute the indeces used in the interpolation
! This includes the special conservative treatment near the block boundary.
! The outer loop selects which offset value is being used for indexing (ie
! which side of the parent block the child will be on.)

      do loop = 1,2

! compute the offset in the parent block appropriate for the different children
      if(loop.eq.1) jchild=1
      if(loop.eq.2) jchild=nchild
      ioff = mod(jchild-1,2)*nxb/2
      joff = mod((jchild-1)/2,2)*nyb/2
      koff = mod((jchild-1)/4,2)*nzb/2

                                                 ! note the 2*nxb and
      do i=il_bnd1,iu_bnd1+1                       ! nxb components in the
       ii = i-nguard                             ! expression for i1 are
       i1 = (ii+nxb*2)/2-nxb+ioff+nguard         ! included so i1 will
       i1p= i1+1                                 ! be correct for -ve
       prol_f_indexx(1,i,loop) = i1              ! values of i also.
       prol_f_indexx(2,i,loop) = i1p             ! (true also for j1,k1)
       enddo

      prol_f_indexy(:,:,loop)=jl_bnd1
      if(ndim.ge.2) then
      do j=jl_bnd1,ju_bnd1+1
       jj = j-nguard
       j1 = (jj+nyb*2)/2-nyb+joff+nguard
       j1p= j1+1
       prol_f_indexy(1,j,loop) = j1
       prol_f_indexy(2,j,loop) = j1p
      enddo
      endif

      prol_f_indexz(:,:,loop)=kl_bnd1
      if(ndim.eq.3) then
      do k=kl_bnd1,ku_bnd1+1
       kk = k-nguard
       k1 = (kk+nzb*2)/2-nzb+koff+nguard
       k1p= k1+1
       prol_f_indexz(1,k,loop) = k1
       prol_f_indexz(2,k,loop) = k1p
      enddo
      endif

      enddo

! set flag to pass error check at the start of prolong_face_fun
      prol_f_init = 100



      call shmem_barrier_all()
      return
      end subroutine amr_prolong_face_fun_init
