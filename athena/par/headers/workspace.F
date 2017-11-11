!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

! workspace module


      module workspace

#include "paramesh_preprocessor.fh"

      use paramesh_dimensions

      private


! workspace arrays
      public :: work, recvw, sendw, tempw, interp_mask_work
#ifndef LIBRARY
      real, save :: work(ilw:iuw,jlw:juw,klw:kuw,maxblocksw,             &
     &                   nvar_work)
      real, save :: recvw(ilw:iuw,jlw:juw,klw:kuw)
      real, save :: sendw(ilw:iuw,jlw:juw,klw:kuw)
      real, save :: tempw(ilw:iuw,jlw:juw,klw:kuw)
      integer, save :: interp_mask_work(nvar_work)
#else
      real, allocatable, save :: work(:,:,:,:,:)
      real, allocatable, save :: recvw(:,:,:)
      real, allocatable, save :: sendw(:,:,:)
      real, allocatable, save :: tempw(:,:,:)
      integer, allocatable, save :: interp_mask_work(:)
#endif

! common block storing the solution for cell-centered quantities.
      public :: work1, recvw1, tempw1, tempw2
#ifndef LIBRARY
      real, save :: work1(ilw1:iuw1,jlw1:juw1,klw1:kuw1,npblks)
      real, save :: recvw1(ilw1:iuw1,jlw1:juw1,klw1:kuw1,npblks)
      real, save :: tempw1(ilw1:iuw1,jlw1:juw1,klw1:kuw1)
#else
      real, allocatable, save :: work1(:,:,:,:)
      real, allocatable, save :: recvw1(:,:,:,:)
      real, allocatable, save :: tempw1(:,:,:)
#endif

#ifdef LIBRARY
      real, allocatable, save,public :: work1_fl(:,:,:)
#else
      real, save,public :: work1_fl(ilw1:iuw1+nxb+2*nguard_work,         &
     .                              jlw1:juw1+(nyb+2*nguard_work)*k2d,   &
     .                              klw1:kuw1+(nzb+2*nguard_work)*k3d)
#endif


#ifdef CURVILINEAR
! arrays used to store geometry information for the working block
      public :: cell_vol_w
#ifndef LIBRARY
      real :: cell_vol_w(ilw1:iuw1,jlw1:juw1,klw1:kuw1)
#else
      real, allocatable :: cell_vol_w(:,:,:)
#endif
#endif /* CURVILINEAR */


! Index arrays used to record destination data values for fine layer
! neighbor guardcells
      integer,public :: f2c_ind_work(2,3,27)


      end module workspace


