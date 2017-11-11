!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"


!------------------------------------------------------------------------------
! physicaldata module
!------------------------------------------------------------------------------




      module physicaldata

      use paramesh_dimensions

      private


!----------------------
! Solution Variables
!----------------------

!---------------------------------------
! Allocate memory for solution variables
!---------------------------------------

! the solution for cell-centered quantities.
      public :: unk, interp_mask_unk
      public :: gcell_on_cc,int_gcell_on_cc
      public :: ngcell_on_cc
      public :: checkp_on_cc
      public :: gcell_on_cc_pointer
#ifndef LIBRARY
      real,save ::  unk(nvar,                                           &
     &                  il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,      &
     &                  maxblocks)
      integer, save :: interp_mask_unk(nvar)
      integer, save :: gcell_on_cc_pointer(nvar)
      logical, save :: gcell_on_cc(nvar)
      logical, save :: int_gcell_on_cc(nvar)
      logical, save :: checkp_on_cc(nvar)
#else
      real,allocatable,save ::  unk(:,:,:,:,:)
      integer,allocatable,save :: interp_mask_unk(:)
      integer,allocatable,save :: gcell_on_cc_pointer(:)
      logical,allocatable,save :: gcell_on_cc(:)
      logical,allocatable,save :: int_gcell_on_cc(:)
      logical,allocatable,save :: checkp_on_cc(:)
#endif
      integer, save :: ngcell_on_cc
      target :: unk

! the solution for cell-face-centered quantities.
      public :: facevarx,facevary,facevarz
      public :: interp_mask_facex,interp_mask_facey,interp_mask_facez
      public :: gcell_on_fc,int_gcell_on_fc
      public :: ngcell_on_fc
      public :: gcell_on_fc_pointer
      public :: checkp_on_fc
#ifndef LIBRARY
      real,save ::  facevarx(nbndvar,                                   &
     &                       il_bnd:iu_bnd+1,jl_bnd:ju_bnd,             &
     &                       kl_bnd:ku_bnd,                             &
     &                       maxblocksf)
      real,save ::  facevary(nbndvar,                                   &
     &                       il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,           &
     &                       kl_bnd:ku_bnd,                             &
     &                       maxblocksf)
      real,save ::  facevarz(nbndvar,                                   &
     &                       il_bnd:iu_bnd,jl_bnd:ju_bnd,               &
     &                       kl_bnd:ku_bnd+k3d,                         &
     &                       maxblocksf)
      integer,save :: interp_mask_facex(nbndvar)
      integer,save :: interp_mask_facey(nbndvar)
      integer,save :: interp_mask_facez(nbndvar)
      integer, save :: gcell_on_fc_pointer(3,nbndvar)
      logical, save :: gcell_on_fc(3,nbndvar)
      logical, save :: int_gcell_on_fc(3,nbndvar)
      logical, save :: checkp_on_fc(3,nbndvar)
#else
      real,allocatable,save ::  facevarx(:,:,:,:,:)
      real,allocatable,save ::  facevary(:,:,:,:,:)
      real,allocatable,save ::  facevarz(:,:,:,:,:)
      integer,allocatable,save :: interp_mask_facex(:)
      integer,allocatable,save :: interp_mask_facey(:)
      integer,allocatable,save :: interp_mask_facez(:)
      integer,allocatable,save :: gcell_on_fc_pointer(:,:)
      logical,allocatable,save :: gcell_on_fc(:,:)
      logical,allocatable,save :: int_gcell_on_fc(:,:)
      logical,allocatable,save :: checkp_on_fc(:,:)
#endif

      integer, save :: ngcell_on_fc(3)
      target :: facevarx,facevary,facevarz

! the solution for cell-edge-centered quantities.
      public :: unk_e_x,unk_e_y,unk_e_z,interp_mask_ec
      public :: gcell_on_ec,int_gcell_on_ec
      public :: ngcell_on_ec
      public :: gcell_on_ec_pointer
      public :: checkp_on_ec
#ifndef LIBRARY
      real,save ::  unk_e_x(nbndvare,                                   &
     &                       il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,           &
     &                       kl_bnd:ku_bnd+k3d,                         &
     &                       maxblocksue)
      real,save ::  unk_e_y(nbndvare,                                   &
     &                       il_bnd:iu_bnd+1,jl_bnd:ju_bnd,             &
     &                       kl_bnd:ku_bnd+k3d,                         &
     &                       maxblocksue)
      real,save ::  unk_e_z(nbndvare,                                   &
     &                       il_bnd:iu_bnd+1,jl_bnd:ju_bnd+k2d,         &
     &                       kl_bnd:ku_bnd,                             &
     &                       maxblocksue)
      integer,save :: interp_mask_ec(nbndvare)
      integer, save :: gcell_on_ec_pointer(3,nbndvare)
      logical, save :: gcell_on_ec(3,nbndvare)
      logical, save :: int_gcell_on_ec(3,nbndvare)
      logical, save :: checkp_on_ec(3,nbndvare)
#else
      real,allocatable,save ::  unk_e_x(:,:,:,:,:)
      real,allocatable,save ::  unk_e_y(:,:,:,:,:)
      real,allocatable,save ::  unk_e_z(:,:,:,:,:)
      integer,allocatable,save :: interp_mask_ec(:)
      integer,allocatable,save :: gcell_on_ec_pointer(:,:)
      logical,allocatable,save :: gcell_on_ec(:,:)
      logical,allocatable,save :: int_gcell_on_ec(:,:)
      logical,allocatable,save :: checkp_on_ec(:,:)
#endif
      integer, save :: ngcell_on_ec(3)
      target :: unk_e_x,unk_e_y,unk_e_z

! the solution for cell-corner based quantities.
      public :: unk_n, interp_mask_nc
      public :: gcell_on_nc,int_gcell_on_nc
      public :: ngcell_on_nc
      public :: gcell_on_nc_pointer
      public :: checkp_on_nc
#ifndef LIBRARY
      real,save ::  unk_n(nbndvarc,                                     &
     &                    il_bnd:iu_bnd+1,jl_bnd:ju_bnd+k2d,            &
     &                    kl_bnd:ku_bnd+k3d,                            &
     &                    maxblocksn)
      integer,save :: interp_mask_nc(nbndvarc)
      integer, save :: gcell_on_nc_pointer(nbndvarc)
      logical, save :: gcell_on_nc(nbndvarc)
      logical, save :: int_gcell_on_nc(nbndvarc)
      logical, save :: checkp_on_nc(nbndvarc)
#else
      real,allocatable,save ::  unk_n(:,:,:,:,:)
      integer,allocatable,save :: interp_mask_nc(:)
      integer,allocatable,save :: gcell_on_nc_pointer(:)
      logical,allocatable,save :: gcell_on_nc(:)
      logical,allocatable,save :: int_gcell_on_nc(:)
      logical,allocatable,save :: checkp_on_nc(:)
#endif

      integer, save :: ngcell_on_nc
      target :: unk_n

!-----------------
! Timestep control
!-----------------

! arrays used for timestep control
      public :: time_loc,dtlevel,phase_dt,loc_cycle,ncyc_local
      public :: ldtcomplete
      real, save    :: dtlevel(maxlevels)
      integer, save :: phase_dt(maxlevels),loc_cycle(maxlevels)
      integer, save :: ncyc_local(maxlevels)
#ifndef LIBRARY
      logical, save :: ldtcomplete(maxblocks_alloc)
      real, save    :: time_loc(maxblocks_alloc)
#else
      logical, allocatable,save :: ldtcomplete(:)
      real, allocatable,save    :: time_loc(:)
#endif


#if defined(VAR_DT) || defined(PRED_CORR)
      public :: t_unk,tfacevarx,tfacevary,tfacevarz
      public :: t_unk_e_x,t_unk_e_y,t_unk_e_z,t_unk_n
#ifndef LIBRARY
      real, save    :: t_unk(nvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,          &
     &                            kl_bnd:ku_bnd,                        &
     &                            maxblocks)
      real, save  :: tfacevarx(nbndvar,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,   &
     &                          kl_bnd:ku_bnd,maxblocksf) 
      real, save  :: tfacevary(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d, &
     &                          kl_bnd:ku_bnd,maxblocksf)
      real, save  :: tfacevarz(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,     &
     &                          kl_bnd:ku_bnd+k3d,maxblocksf)

      real, save  :: t_unk_e_x(nbndvare,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,&
     &                          kl_bnd:ku_bnd+k3d,maxblocksue) 
      real, save  :: t_unk_e_y(nbndvare,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,  &
     &                          kl_bnd:ku_bnd+k3d,maxblocksue) 
      real, save  :: t_unk_e_z(nbndvare,il_bnd:iu_bnd+1,                &
     &                          jl_bnd:ju_bnd+k2d,                      &
     &                          kl_bnd:ku_bnd,maxblocksue) 
      real, save  :: t_unk_n(nbndvarc,il_bnd:iu_bnd+1,                  &
     &                          jl_bnd:ju_bnd+k2d,                      &
     &                          kl_bnd:ku_bnd+k3d,maxblocksn) 
#else
      real, allocatable,save  :: t_unk(:,:,:,:,:)
      real, allocatable,save  :: tfacevarx(:,:,:,:,:)
      real, allocatable,save  :: tfacevary(:,:,:,:,:)
      real, allocatable,save  :: tfacevarz(:,:,:,:,:)
      real, allocatable,save  :: t_unk_e_x(:,:,:,:,:)
      real, allocatable,save  :: t_unk_e_y(:,:,:,:,:)
      real, allocatable,save  :: t_unk_e_z(:,:,:,:,:)
      real, allocatable,save  :: t_unk_n(:,:,:,:,:)
#endif /* SHMEM_USED */
#endif


!-----------------------------------------------------------------
! include header file defining 1blk data structure

!------------------------------------------------------------------------------
! data_1blk
!------------------------------------------------------------------------------
!
! This file declares the storage space used to handle the `current
! working block' when the user decides not to reserve permanent
! storage space for guardcells for all blocks, but instead to 
! fill guardcells as needed. This strategy requires 2 working blocks,
! one for the leaf node and one for its parent.

      public :: unk1,facevarx1,facevary1,facevarz1
      public :: unk_e_x1,unk_e_y1,unk_e_z1
      public :: unk_n1

#ifndef LIBRARY
! the solution for cell-centered quantities.
      real, save :: unk1(nvar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,          &
     &                        kl_bnd1:ku_bnd1,                          &
     &                        npblks)
! the solution for cell-face-centered quantities.
      real, save :: facevarx1(nbndvar,il_bnd1:iu_bnd1+1,                &
     &                                jl_bnd1:ju_bnd1,                  &
     &                                kl_bnd1:ku_bnd1,npblks)
      real, save :: facevary1(nbndvar,il_bnd1:iu_bnd1,                  &
     &                                jl_bnd1:ju_bnd1+k2d,              &
     &                                kl_bnd1:ku_bnd1,npblks)
      real, save :: facevarz1(nbndvar,il_bnd1:iu_bnd1,                  &
     &                                jl_bnd1:ju_bnd1,                  &
     &                                kl_bnd1:ku_bnd1+k3d,npblks)  
! the solution for cell-edge-centered quantities.
      real,save ::  unk_e_x1(nbndvare,                                  &
     &                       il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,       &
     &                       kl_bnd1:ku_bnd1+k3d,                       &
     &                       npblks)
      real,save ::  unk_e_y1(nbndvare,                                  &
     &                       il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,         &
     &                       kl_bnd1:ku_bnd1+k3d,                        &
     &                       npblks)
      real,save ::  unk_e_z1(nbndvare,                                  &
     &                       il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,     &
     &                       kl_bnd1:ku_bnd1,                           &
     &                       npblks)
! the solution for cell-corner based quantities.
      real,save ::  unk_n1(nbndvarc,                                    &
     &                       il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,     &
     &                       kl_bnd1:ku_bnd1+k3d,                       &
     &                       npblks)
#else
! the solution for cell-centered quantities.
      real, save, allocatable :: unk1(:,:,:,:,:)
! the solution for cell-face-centered quantities.
      real, save, allocatable :: facevarx1(:,:,:,:,:)
      real, save, allocatable :: facevary1(:,:,:,:,:)
      real, save, allocatable :: facevarz1(:,:,:,:,:)
! the solution for cell-edge-centered quantities.
      real, save, allocatable ::  unk_e_x1(:,:,:,:,:)
      real, save, allocatable ::  unk_e_y1(:,:,:,:,:)
      real, save, allocatable ::  unk_e_z1(:,:,:,:,:)
! the solution for cell-corner based quantities.
      real, save, allocatable ::  unk_n1(:,:,:,:,:)
#endif

      target :: unk1
      target :: facevarx1, facevary1, facevarz1
      target :: unk_e_x1, unk_e_y1, unk_e_z1
      target :: unk_n1

!------------------------------------------------------------------------------
! workspace arrays for a finer block data layer, to be used in 
! amr_1blk_guardcell_f_to_c
!------------------------------------------------------------------------------

#ifdef LIBRARY
! for cell-centered quantities.
      real, save, allocatable,public ::  unk1_fl(:,:,:,:)
! for cell-face quantities.
      real, save, allocatable,public ::  facevarx1_fl(:,:,:,:)
      real, save, allocatable,public ::  facevary1_fl(:,:,:,:)
      real, save, allocatable,public ::  facevarz1_fl(:,:,:,:)
! for cell-edge quantities.
      real, save, allocatable,public ::  unk_e_x1_fl(:,:,:,:)
      real, save, allocatable,public ::  unk_e_y1_fl(:,:,:,:)
      real, save, allocatable,public ::  unk_e_z1_fl(:,:,:,:)
! for cell-corner quantities.
      real, save, allocatable,public ::  unk_n1_fl(:,:,:,:)
#else
      real, save, public :: unk1_fl(nvar,                               &
     .                              il_bnd1:iu_bnd1+nxb+2*nguard,       &
     .                              jl_bnd1:ju_bnd1+(nyb+2*nguard)*k2d, &
     .                              kl_bnd1:ku_bnd1+(nzb+2*nguard)*k3d )
      real, save, public :: facevarx1_fl(nbndvar,                       &
     .                            il_bnd1:iu_bnd1+nxb+2*nguard+1,       &
     .                            jl_bnd1:ju_bnd1+(nyb+2*nguard)*k2d,   &
     .                            kl_bnd1:ku_bnd1+(nzb+2*nguard)*k3d )
      real, save, public :: facevary1_fl(nbndvar,                       &
     .                            il_bnd1:iu_bnd1+nxb+2*nguard,         &
     .                            jl_bnd1:ju_bnd1+(nyb+2*nguard+1)*k2d, &
     .                            kl_bnd1:ku_bnd1+(nzb+2*nguard)*k3d )
      real, save, public :: facevarz1_fl(nbndvar,                       &
     .                            il_bnd1:iu_bnd1+nxb+2*nguard,         &
     .                            jl_bnd1:ju_bnd1+(nyb+2*nguard)*k2d,   &
     .                            kl_bnd1:ku_bnd1+(nzb+2*nguard+1)*k3d )
      real, save, public :: unk_e_x1_fl(nbndvare,                       &
     .                            il_bnd1:iu_bnd1+nxb+2*nguard,         &
     .                            jl_bnd1:ju_bnd1+(nyb+2*nguard+1)*k2d, &
     .                            kl_bnd1:ku_bnd1+(nzb+2*nguard+1)*k3d )
      real, save, public :: unk_e_y1_fl(nbndvare,                       &
     .                            il_bnd1:iu_bnd1+nxb+2*nguard+1,       &
     .                            jl_bnd1:ju_bnd1+(nyb+2*nguard)*k2d,   &
     .                            kl_bnd1:ku_bnd1+(nzb+2*nguard+1)*k3d )
      real, save, public :: unk_e_z1_fl(nbndvare,                       &
     .                            il_bnd1:iu_bnd1+nxb+2*nguard+1,       &
     .                            jl_bnd1:ju_bnd1+(nyb+2*nguard+1)*k2d, &
     .                            kl_bnd1:ku_bnd1+(nzb+2*nguard)*k3d )
      real, save, public :: unk_n1_fl(nbndvarc,
     .                            il_bnd1:iu_bnd1+nxb+2*nguard+1,       &
     .                            jl_bnd1:ju_bnd1+(nyb+2*nguard+1)*k2d, &
     .                            kl_bnd1:ku_bnd1+(nzb+2*nguard+1)*k3d )
#endif

!------------------------------------------------------------------------------
! temporary copy of solution to be used when storing solution prior
! to use of amr_1blk_guardcell
!------------------------------------------------------------------------------

      public :: gt_unk,gt_facevarx,gt_facevary,gt_facevarz
      public :: gt_unk_e_x,gt_unk_e_y,gt_unk_e_z,gt_unk_n
#ifndef LIBRARY
      real, save :: gt_unk(nvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,            &
     &                          kl_bnd:ku_bnd,maxblocks_gt)
#else
      real, allocatable, save :: gt_unk(:,:,:,:,:)
#endif
#ifdef NO_PERMANENT_GUARDCELLS
#ifndef LIBRARY
      real, save :: gt_facevarx(nbndvar,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,  &
     &                         kl_bnd:ku_bnd,maxblocksf_gt) 
      real, save :: gt_facevary(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,&
     &                         kl_bnd:ku_bnd,maxblocksf_gt)
      real, save :: gt_facevarz(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,    &
     &                         kl_bnd:ku_bnd+k3d,maxblocksf_gt)
#else
      real, allocatable, save :: gt_facevarx(:,:,:,:,:)
      real, allocatable, save :: gt_facevary(:,:,:,:,:)
      real, allocatable, save :: gt_facevarz(:,:,:,:,:)
#endif
#else /*  NO_PERMANENT_GUARDCELLS */
#ifndef LIBRARY
      real, save :: gt_facevarx(nbndvar,1:2,jl_bnd:ju_bnd,              &
     &                         kl_bnd:ku_bnd,maxblocksf) 
      real, save :: gt_facevary(nbndvar,il_bnd:iu_bnd,1:1+k2d,          &
     &                         kl_bnd:ku_bnd,maxblocksf)
      real, save :: gt_facevarz(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,    &
     &                         1:1+k3d,maxblocksf)
#else
      real, allocatable, save :: gt_facevarx(:,:,:,:,:)
      real, allocatable, save :: gt_facevary(:,:,:,:,:)
      real, allocatable, save :: gt_facevarz(:,:,:,:,:)
#endif
#endif

#ifndef LIBRARY
      real, save :: gt_unk_e_x(nbndvare,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,&
     &                         kl_bnd:ku_bnd+k3d,maxblocksue_gt) 
      real, save :: gt_unk_e_y(nbndvare,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,  &
     &                         kl_bnd:ku_bnd+k3d,maxblocksue_gt)
      real, save :: gt_unk_e_z(nbndvare,il_bnd:iu_bnd+1,                &
     &                         jl_bnd:ju_bnd+k2d,                       &
     &                         kl_bnd:ku_bnd,maxblocksue_gt)
      real, save :: gt_unk_n(nbndvarc,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+k2d,&
     &                         kl_bnd:ku_bnd+k3d,maxblocksn_gt)
#else
      real, allocatable, save :: gt_unk_e_x(:,:,:,:,:)
      real, allocatable, save :: gt_unk_e_y(:,:,:,:,:)
      real, allocatable, save :: gt_unk_e_z(:,:,:,:,:)
      real, allocatable, save :: gt_unk_n(:,:,:,:,:)
#endif

!
! Variables used to control data caching. This helps to avoid unnecessary
! repetition of some guardcell and surrounding-block mapping operations
      public ::  pcache_blk_u,pcache_pe_u,pcache_blk_w,pcache_pe_w
      public ::  lnew_parent
      integer, save :: pcache_blk_u,pcache_pe_u
      integer, save :: pcache_blk_w,pcache_pe_w
      logical, save :: lnew_parent

!-----------------------------------------------------------------
! include header file defining data structure on cell faces

!------------------------------------------------------------------------
! block_boundary_data


!
! This file defines a data structure to be used for quantities
! which may need to be defined at grid block interfaces, eg fluxes,
! pressures.
!

! The convention for relating varaibles associated with cell faces to the
! variables defined at cell centers is as follows:

! If iface_off=0 :
!         the array facevarx(:,i,j,k,:) for example defines data
!         on the x(i-1/2) face of the (i,j,k)-th mesh cell. 
! If iface_off=-1 :
!         the array facevarx(:,i,j,k,:) for example defines data
!         on the x(i+1/2) face of the (i,j,k)-th mesh cell. 


! storage used for fluxes at block boundaries. This is used when conservation
! constraints need to be imposed.

      public :: nfluxvar,nfluxes,maxblocksfl
      public :: flux_x,flux_y,flux_z
      public :: tflux_x,tflux_y,tflux_z
      integer, parameter :: nfluxvar =  N_FLUX_VAR

      integer, parameter :: nfluxes=max(1,nfluxvar)

#ifndef LIBRARY
      integer, parameter :: maxblocksfl=                                &
     &                       1+(maxblocks_alloc-1)*min(1,nfluxvar)
      real, save ::  flux_x(nfluxes,1:2,                                &
     &                      jl_bnd:ju_bnd,kl_bnd:ku_bnd,maxblocksfl)
      real, save ::  flux_y(nfluxes,il_bnd:iu_bnd,                      &
     &                              1:2,kl_bnd:ku_bnd,maxblocksfl)
      real, save ::  flux_z(nfluxes,il_bnd:iu_bnd,                      &
     &                              jl_bnd:ju_bnd,1:2,maxblocksfl)
      real, save :: tflux_x(nfluxes,1:2,                                &
     &                              jl_bnd:ju_bnd,                      &
     &                              kl_bnd:ku_bnd,maxblocksfl)
      real, save :: tflux_y(nfluxes,il_bnd:iu_bnd,                      &
     &                              1:2,kl_bnd:ku_bnd,maxblocksfl)
      real, save :: tflux_z(nfluxes,il_bnd:iu_bnd,                      &
     &                              jl_bnd:ju_bnd,1:2,maxblocksfl)
#else
      integer :: maxblocksfl
      real, allocatable, save ::  flux_x(:,:,:,:,:)
      real, allocatable, save ::  flux_y(:,:,:,:,:)
      real, allocatable, save ::  flux_z(:,:,:,:,:)
      real, allocatable, save ::  tflux_x(:,:,:,:,:)
      real, allocatable, save ::  tflux_y(:,:,:,:,:)
      real, allocatable, save ::  tflux_z(:,:,:,:,:)
#endif
      target :: flux_x, flux_y, flux_z 


#ifdef VAR_DT
! temporary flux storage needed inside amr_flux_conserve when using
! variable timestep
      public :: ttflux_x,ttflux_y,ttflux_z
#ifndef LIBRARY
      real, save :: ttflux_x(nfluxes,1:2,jl_bnd:ju_bnd,                 &
     &                               kl_bnd:ku_bnd,maxblocksfl)
      real, save :: ttflux_y(nfluxes,il_bnd:iu_bnd,                     &
     &                       1:2,kl_bnd:ku_bnd,maxblocksfl)
      real, save :: ttflux_z(nfluxes,il_bnd:iu_bnd,                     &
     &                               jl_bnd:ju_bnd,1:2,maxblocksfl)
#else
      real, allocatable, save :: ttflux_x(:,:,:,:,:)
      real, allocatable, save :: ttflux_y(:,:,:,:,:)
      real, allocatable, save :: ttflux_z(:,:,:,:,:)
#endif
#endif

! storage used for cell edges at block boundaries. 
! This is used when quantities located at cell edge centers need to
! be used consistently at the boundaries between blocks at different
! refinement levels.

      public ::  nedgevar,nedges,maxblockse
      public :: bedge_facex_y,bedge_facex_z,bedge_facey_x
      public :: bedge_facey_z,bedge_facez_x,bedge_facez_y
      public :: recvarx1e,recvary1e,recvarz1e
      public :: recvarx2e,recvary2e,recvarz2e

#ifndef LIBRARY
      integer, parameter :: nedgevar1=N_EDGE_VAR
      integer, parameter :: nedgevar=max(nedgevar1,nvaredge)
      integer, parameter :: nedges=max(1,nedgevar)

      integer, parameter :: maxblockse=                                 &
     &                           1+(maxblocks_alloc-1)*min(1,nedgevar)
      real, save :: bedge_facex_y(nedges,1:2,jl_bnd:ju_bnd+1,           &
     &                                   kl_bnd:ku_bnd+1,maxblockse)
      real, save ::  bedge_facex_z(nedges,1:2,jl_bnd:ju_bnd+1,          &
     &                                    kl_bnd:ku_bnd+1,maxblockse)
      real, save ::  bedge_facey_x(nedges,il_bnd:iu_bnd+1,1:2,          &
     &                                    kl_bnd:ku_bnd+1,maxblockse)
      real, save ::  bedge_facey_z(nedges,il_bnd:iu_bnd+1,1:2,          &
     &                                    kl_bnd:ku_bnd+1,maxblockse)
      real, save ::  bedge_facez_x(nedges,il_bnd:iu_bnd+1,              &
     &                             jl_bnd:ju_bnd+1,1:2,maxblockse)
      real, save ::  bedge_facez_y(nedges,il_bnd:iu_bnd+1,              &
     &                                    jl_bnd:ju_bnd+1,              &
     &                                    1:2,maxblockse)
#else
      integer :: nedgevar1
      integer :: nedgevar
      integer :: nedges

      integer :: maxblockse
      real, allocatable, save ::  bedge_facex_y(:,:,:,:,:)
      real, allocatable, save ::  bedge_facex_z(:,:,:,:,:)
      real, allocatable, save ::  bedge_facey_x(:,:,:,:,:)
      real, allocatable, save ::  bedge_facey_z(:,:,:,:,:)
      real, allocatable, save ::  bedge_facez_x(:,:,:,:,:)
      real, allocatable, save ::  bedge_facez_y(:,:,:,:,:)
#endif
#ifndef LIBRARY
      real, save ::  recvarx1e(nedges,1:2,jl_bnd:ju_bnd+1,              &
     &                                kl_bnd:ku_bnd+1)
      real, save ::  recvary1e(nedges,il_bnd:iu_bnd+1,1:2,              &
     &                                kl_bnd:ku_bnd+1)
      real, save ::  recvarz1e(nedges,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1,  &
     &                                1:2)
      real, save ::  recvarx2e(nedges,1:2,jl_bnd:ju_bnd+1,              &
     &                                kl_bnd:ku_bnd+1)
      real, save ::  recvary2e(nedges,il_bnd:iu_bnd+1,1:2,              &
     &                                kl_bnd:ku_bnd+1)
      real, save ::  recvarz2e(nedges,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+1,  &
     &                                1:2)
#else
      real, allocatable, save ::  recvarx1e(:,:,:,:)
      real, allocatable, save ::  recvary1e(:,:,:,:)
      real, allocatable, save ::  recvarz1e(:,:,:,:)
      real, allocatable, save ::  recvarx2e(:,:,:,:)
      real, allocatable, save ::  recvary2e(:,:,:,:)
      real, allocatable, save ::  recvarz2e(:,:,:,:)
#endif

      public :: tbedge_facex_y,tbedge_facex_z,tbedge_facey_x
      public :: tbedge_facey_z,tbedge_facez_x,tbedge_facez_y      
#ifndef LIBRARY
      real, save :: tbedge_facex_y(nedges,1:2,jl_bnd:ju_bnd+1,          &
     &                                    kl_bnd:ku_bnd+1,maxblockse)
      real, save :: tbedge_facex_z(nedges,1:2,jl_bnd:ju_bnd+1,          &
     &                                    kl_bnd:ku_bnd+1,maxblockse)
      real, save :: tbedge_facey_x(nedges,il_bnd:iu_bnd+1,1:2,          &
     &                                    kl_bnd:ku_bnd+1,maxblockse)
      real, save :: tbedge_facey_z(nedges,il_bnd:iu_bnd+1,1:2,          &
     &                                    kl_bnd:ku_bnd+1,maxblockse)
      real, save :: tbedge_facez_x(nedges,il_bnd:iu_bnd+1,              &
     &                                    jl_bnd:ju_bnd+1,              &
     &                                    1:2,maxblockse)
      real, save :: tbedge_facez_y(nedges,il_bnd:iu_bnd+1,              &
     &                                    jl_bnd:ju_bnd+1,              &
     &                                    1:2,maxblockse)
#else
      real, allocatable, save ::  tbedge_facex_y(:,:,:,:,:)
      real, allocatable, save ::  tbedge_facex_z(:,:,:,:,:)
      real, allocatable, save ::  tbedge_facey_x(:,:,:,:,:)
      real, allocatable, save ::  tbedge_facey_z(:,:,:,:,:)
      real, allocatable, save ::  tbedge_facez_x(:,:,:,:,:)
      real, allocatable, save ::  tbedge_facez_y(:,:,:,:,:)
#endif
#ifdef VAR_DT
      public :: ttbedge_facex_y,ttbedge_facex_z,ttbedge_facey_x
      public :: ttbedge_facey_z,ttbedge_facez_x,ttbedge_facez_y      
#ifndef LIBRARY
      real, save :: ttbedge_facex_y(nedges,1:2,jl_bnd:ju_bnd+1,         &
     &                                     kl_bnd:ku_bnd+1,maxblockse)
      real, save :: ttbedge_facex_z(nedges,1:2,jl_bnd:ju_bnd+1,         &
     &                                     kl_bnd:ku_bnd+1,maxblockse)
      real, save :: ttbedge_facey_x(nedges,il_bnd:iu_bnd+1,1:2,         &
     &                                     kl_bnd:ku_bnd+1,maxblockse)
      real, save :: ttbedge_facey_z(nedges,il_bnd:iu_bnd+1,1:2,         &
     &                                     kl_bnd:ku_bnd+1,maxblockse)
      real, save :: ttbedge_facez_x(nedges,il_bnd:iu_bnd+1,             &
     &                                     jl_bnd:ju_bnd+1,             &
     &                                     1:2,maxblockse)
      real, save :: ttbedge_facez_y(nedges,il_bnd:iu_bnd+1,             &
     &                                     jl_bnd:ju_bnd+1,             &   
     &                                     1:2,maxblockse)
#else
      real, allocatable, save ::  tbedge_facex_y(:,:,:,:,:)
      real, allocatable, save ::  tbedge_facex_z(:,:,:,:,:)
      real, allocatable, save ::  tbedge_facey_x(:,:,:,:,:)
      real, allocatable, save ::  tbedge_facey_z(:,:,:,:,:)
      real, allocatable, save ::  tbedge_facez_x(:,:,:,:,:)
      real, allocatable, save ::  tbedge_facez_y(:,:,:,:,:)
#endif
#endif




#ifdef CURVILINEAR
! arrays used to store geometry information for the working block
      public :: cell_vol
      public :: cell_area1,cell_area2,cell_area3
      public :: cell_leng1,cell_leng2,cell_leng3
      public :: cell_face_coord1,cell_face_coord2,cell_face_coord3
#ifndef LIBRARY
      real :: cell_vol(il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1)
      real :: cell_area1(il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,             &
     &                                        kl_bnd1:ku_bnd1)
      real :: cell_area2(il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,           &
     &                                        kl_bnd1:ku_bnd1)
      real :: cell_area3(il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,               &
     &                                      kl_bnd1:ku_bnd1+k3d)
      real :: cell_leng1(il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,           &
     &                                        kl_bnd1:ku_bnd1+k3d)
      real :: cell_leng2(il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,             &
     &                                        kl_bnd1:ku_bnd1+k3d)
      real :: cell_leng3(il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,         &
     &                                        kl_bnd1:ku_bnd1)
      real :: cell_face_coord1(il_bnd1:iu_bnd1+1)
      real :: cell_face_coord2(jl_bnd1:ju_bnd1+k2d)
      real :: cell_face_coord3(kl_bnd1:ku_bnd1+k3d)
#else
      real, allocatable :: cell_vol(:,:,:)
      real, allocatable  :: cell_area1(:,:,:)
      real, allocatable  :: cell_area2(:,:,:)
      real, allocatable  :: cell_area3(:,:,:)
      real, allocatable  :: cell_leng1(:,:,:)
      real, allocatable  :: cell_leng2(:,:,:)
      real, allocatable  :: cell_leng3(:,:,:)
      real, allocatable  :: cell_face_coord1(:)
      real, allocatable  :: cell_face_coord2(:)
      real, allocatable  :: cell_face_coord3(:)
#endif
#endif /* CURVILINEAR */

! workspace arrays used for inter-block communications
      public :: nbndmax
      public :: recvarx1,recvary1,recvarz1
      public :: recvarxf,recvaryf,recvarzf
      public :: bndtempx1,bndtempy1,bndtempz1
#ifndef LIBRARY
      integer, parameter :: nbndmax=max(nbndvar,nfluxes)
      real, save :: recvarx1(nbndmax,1:2,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
      real, save :: recvary1(nbndmax,il_bnd:iu_bnd,1:2,kl_bnd:ku_bnd)
      real, save :: recvarz1(nbndmax,il_bnd:iu_bnd,jl_bnd:ju_bnd,1:2)
      real, save :: recvarxf(nfluxes,1:2,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
      real, save :: recvaryf(nfluxes,il_bnd:iu_bnd,1:2,kl_bnd:ku_bnd)
      real, save :: recvarzf(nfluxes,il_bnd:iu_bnd,jl_bnd:ju_bnd,1:2)
      real, save :: bndtempx1(nfluxes,1:2,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
      real, save :: bndtempy1(nfluxes,il_bnd:iu_bnd,1:2,kl_bnd:ku_bnd)
      real, save :: bndtempz1(nfluxes,il_bnd:iu_bnd,jl_bnd:ju_bnd,1:2)
#else
      integer :: nbndmax
      real, allocatable, save :: recvarx1(:,:,:,:)
      real, allocatable, save :: recvary1(:,:,:,:)
      real, allocatable, save :: recvarz1(:,:,:,:)
      real, allocatable, save :: recvarxf(:,:,:,:)
      real, allocatable, save :: recvaryf(:,:,:,:)
      real, allocatable, save :: recvarzf(:,:,:,:)
      real, allocatable, save :: bndtempx1(:,:,:,:)
      real, allocatable, save :: bndtempy1(:,:,:,:)
      real, allocatable, save :: bndtempz1(:,:,:,:)
#endif



! parameters used in communication calls
      public :: len_block_bndx,len_block_bndy,len_block_bndz
      public :: len_block_ex,len_block_ey,len_block_ez
#ifndef LIBRARY
      integer, parameter :: len_block_bndx=2*ju_bnd*ku_bnd
      integer, parameter :: len_block_bndy=2*iu_bnd*ku_bnd
      integer, parameter :: len_block_bndz=2*iu_bnd*ju_bnd
      integer, parameter :: len_block_ex=2*(ju_bnd+k2d)*(ku_bnd+k3d)
      integer, parameter :: len_block_ey=2*(iu_bnd+1  )*(ku_bnd+k3d)
      integer, parameter :: len_block_ez=2*(iu_bnd+1  )*(ju_bnd+k2d)
#else
      integer :: len_block_bndx
      integer :: len_block_bndy
      integer :: len_block_bndz
      integer :: len_block_ex
      integer :: len_block_ey
      integer :: len_block_ez
#endif

!-----------------------------------------------------------------
! Array used to store variables which make up any divergence
! free fields
      integer,allocatable,save,public :: i_divf_fc_vars(:,:)

!-----------------------------------------------------------------
! Index arrays used in boundary condition routines.
      integer,public :: bc_index_i(2,3,5)
      integer,public :: bc_index_j(2,3,5)
      integer,public :: bc_index_k(2,3,5)


!-----------------------------------------------------------------
! Logical flags required to signal algorithmic states
      logical, public :: lrestrict_in_progress
      logical, public :: lprolong_in_progress
      logical, public :: lguard_in_progress

!-----------------------------------------------------------------
! Logical flag to indicate if restrictionless guardcell filling
! has been selected
      logical,public :: l_f_to_c

! Index arrays used to record destination data values for fine layer
! neighbor guardcells
      integer,public :: f2c_ind_unk(2,3,27)
      integer,public :: f2c_ind_facex(2,3,27)
      integer,public :: f2c_ind_facey(2,3,27)
      integer,public :: f2c_ind_facez(2,3,27)
      integer,public :: f2c_ind_unkex(2,3,27)
      integer,public :: f2c_ind_unkey(2,3,27)
      integer,public :: f2c_ind_unkez(2,3,27)
      integer,public :: f2c_ind_unkn(2,3,27)

!-----------------------------------------------------------------
! Mpi communication pattern identifier
      integer, public :: mpi_pattern_id

!-----------------------------------------------------------------
! Error trapping and management

! To record whether amr_gsurrounding_blks has been called.
      public :: gsurrblks_set
      integer :: gsurrblks_set

! a counter which can be used to keep track of calls to routines
      public :: instance
      integer :: instance

!-----------------------------------------------------------------

      end module physicaldata
