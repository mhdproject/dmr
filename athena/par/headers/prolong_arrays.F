!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"
!-----------------------------------------------------------------
! prolong_arrays module



      module prolong_arrays

      use paramesh_dimensions

      private


      public :: prol_dx,prol_dy,prol_dz
      public :: prol_indexx,prol_indexy,prol_indexz,prol_init
      public :: prol_f_dx,prol_f_dy,prol_f_dz
      public :: prol_f_indexx,prol_f_indexy,prol_f_indexz
      public :: prol_f_init
      public :: prolw_dx,prolw_dy,prolw_dz
      public :: prolw_indexx,prolw_indexy,prolw_indexz,prolw_init
#ifndef LIBRARY
      real, save :: prol_dx(il_bnd1:iu_bnd1)
      real, save :: prol_dy(jl_bnd1:ju_bnd1) 
      real, save :: prol_dz(kl_bnd1:ku_bnd1)
      integer, save :: prol_indexx(2,il_bnd1:iu_bnd1,2)
      integer, save :: prol_indexy(2,jl_bnd1:ju_bnd1,2)
      integer, save :: prol_indexz(2,kl_bnd1:ku_bnd1,2)
      real, save :: prol_f_dx(il_bnd1:iu_bnd1+1)
      real, save :: prol_f_dy(jl_bnd1:ju_bnd1+k2d)
      real, save :: prol_f_dz(kl_bnd1:ku_bnd1+k3d)
      integer, save :: prol_f_indexx(2,il_bnd1:iu_bnd1+1,2)
      integer, save :: prol_f_indexy(2,jl_bnd1:ju_bnd1+k2d,2)
      integer, save :: prol_f_indexz(2,kl_bnd1:ku_bnd1+k3d,2)
      real, save :: prolw_dx(ilw1:iuw1)
      real, save :: prolw_dy(jlw1:juw1)
      real, save :: prolw_dz(klw1:kuw1)
      integer, save :: prolw_indexx(2,ilw1:iuw1,2)
      integer, save :: prolw_indexy(2,jlw1:juw1,2)
      integer, save :: prolw_indexz(2,klw1:kuw1,2)
#else
      real, save, allocatable :: prol_dx(:)
      real, save, allocatable :: prol_dy(:) 
      real, save, allocatable :: prol_dz(:)
      integer, save, allocatable :: prol_indexx(:,:,:)
      integer, save, allocatable :: prol_indexy(:,:,:)
      integer, save, allocatable :: prol_indexz(:,:,:)
      real, save, allocatable :: prol_f_dx(:)
      real, save, allocatable :: prol_f_dy(:)
      real, save, allocatable :: prol_f_dz(:)
      integer, save, allocatable :: prol_f_indexx(:,:,:)
      integer, save, allocatable :: prol_f_indexy(:,:,:)
      integer, save, allocatable :: prol_f_indexz(:,:,:)
      real, save, allocatable :: prolw_dx(:)
      real, save, allocatable :: prolw_dy(:)
      real, save, allocatable :: prolw_dz(:)
      integer, save, allocatable :: prolw_indexx(:,:,:)
      integer, save, allocatable :: prolw_indexy(:,:,:)
      integer, save, allocatable :: prolw_indexz(:,:,:)
#endif
      integer, save :: prol_init
      integer, save :: prol_f_init
      integer, save :: prolw_init

      public :: prol_fc_dbz, prol_fc_dbz_ivar, prol_fc_dbz_n
      logical, save :: prol_fc_dbz = .false.
      integer, allocatable, save :: prol_fc_dbz_ivar(:,:)
      integer, save :: prol_fc_dbz_n = 0

      public :: prol_fc_clean_divb, prol_fc_clean_divb_ivar, 
     &     prol_fc_clean_divb_n
      logical, save :: prol_fc_clean_divb = .false.
      integer, allocatable, save :: prol_fc_clean_divb_ivar(:,:)
      integer, save :: prol_fc_clean_divb_n = 0

      end module prolong_arrays
