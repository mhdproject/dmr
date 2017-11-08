!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

! paramesh_dimensions module
!------------------------------------------------------------------------------


      module paramesh_dimensions

#include "paramesh_preprocessor.fh"

!------------------
! Model dimensions
!------------------

! set physical dimension of model and number of edges on each grid block

#ifndef LIBRARY
      integer, parameter ::  ndim = N_DIM
#else
      integer :: ndim
#endif

!pmn      integer, parameter ::  l2p5d = L_2P5DIM
#ifndef LIBRARY
      integer, parameter ::  l2p5d = L_2P5DIM
      integer, parameter ::  nbedges = ndim*2**(ndim-1)
#else
      integer ::  l2p5d
      integer :: nbedges
#endif



! an increment variable for the z dimension to enable the same code to
! work for 2D or 3D models.
#ifndef LIBRARY
      integer, parameter :: k3d=(ndim-1)/2
      integer, parameter :: k2d=ndim/2
      integer, parameter :: k1d=1
#else
      integer :: k3d, k2d, k1d
#endif





!----------------------
! Grid Block Properties
!----------------------

! set size of grid blocks
! The spatial relationship between grid points in parents
! and their children is subtley different depending on whether
! the block dimensions (ie nxb,nyb,nzb) are even or odd. This has 
! significant consequences when defining interpolation within
! restriction or prolongation operations.
#ifndef LIBRARY
      integer, parameter :: nxb = NX_B
      integer, parameter :: nyb = NY_B
      integer, parameter :: nzb = NZ_B
      integer, parameter :: maxdim=max(nxb,nyb,nzb)
#else
      integer :: nxb
      integer :: nyb
      integer :: nzb
      integer :: maxdim
#endif


! these guard cell offsets are required to accomodate differences
! in cases when block dimensions are odd or even
#ifndef LIBRARY
      integer, parameter :: gc_off_x=mod(nxb,2)
      integer, parameter :: gc_off_y=mod(nyb,2)
      integer, parameter :: gc_off_z=mod(nzb,2)
#else
      integer :: gc_off_x
      integer :: gc_off_y
      integer :: gc_off_z
#endif
! set the maximum number of blocks per processor
#ifndef LIBRARY
      integer, parameter :: maxblocks = MAX_BLOCKS
!      integer, parameter :: mxblks_buf = max(1,maxblocks/10)
!      integer, parameter :: maxblocks_alloc = maxblocks + mxblks_buf
      integer, parameter :: maxblocks_alloc = maxblocks*4
#else
      integer :: maxblocks
      integer :: mxblks_buf
      integer :: maxblocks_alloc
#endif

! set the number of guard cell layers at each boundary
      integer, parameter :: nguard = N_GUARD_CELLS


! common block storing the solution for cell-centered quantities.
#ifdef NO_PERMANENT_GUARDCELLS
      integer, parameter :: npgs = 0
#else
      integer, parameter :: npgs = 1
#endif





!----------------------
! Solution Variables
!----------------------

! Set number of unknowns associated with each grid cell.
! If you are using a multi-step timestep integration algorithm 
! (eg predictor-corrector) then the recommended value for nvar is 
!             nvarp*(nphase + 1)
! where nvarp denotes the number of physical variables (ie. 1D hydro would
! have nvarp=3, for mass, momentum and energy), nphase is the number of 
! stages in an integration timestep (ie predictor-corrector would have 
! nphase=2). Similar considerations apply for nfacevar.
#ifndef LIBRARY
      integer, parameter :: nvar = N_VAR

! The number of data words needed on a cell face is set by nfacevar.
      integer, parameter :: nfacevar = N_FACEVAR

! The number of data words needed on cell edges is set by nvaredge.
      integer, parameter :: nvaredge = N_VAR_EDGE

! The number of data words needed at cell corners is set by nvarcorn.
      integer, parameter :: nvarcorn = N_VAR_CORN

#else
      integer :: nvar

! The number of data words needed on a cell face is set by nfacevar.
      integer :: nfacevar

! The number of data words needed on cell edges is set by nvaredge.
      integer :: nvaredge

! The number of data words needed at cell corners is set by nvarcorn.
      integer :: nvarcorn

#endif

! The convention for relating variables associated with cell faces to the
! variables defined at cell centers is as follows:

! If iface_off=0 :
!         the array facevarx(:,i,j,k,:) for example defines data
!         on the x(i-1/2) face of the (i,j,k)-th mesh cell.
! If iface_off=-1 :
!         the array facevarx(:,i,j,k,:) for example defines data
!         on the x(i+1/2) face of the (i,j,k)-th mesh cell.

      integer, parameter :: iface_off= FACE_INDEX_OFFSET





!------------------------------------------------
! ! Declare dimensions for the solution variables
!------------------------------------------------

! Cell centered data bounds
#ifndef LIBRARY
      integer, parameter :: il_bnd=1, iu_bnd=nxb+2*nguard*npgs
      integer, parameter :: jl_bnd=1, ju_bnd=nyb+2*nguard*npgs*k2d
      integer, parameter :: kl_bnd=1, ku_bnd=nzb+2*nguard*npgs*k3d
#else
      integer :: il_bnd, iu_bnd
      integer :: jl_bnd, ju_bnd
      integer :: kl_bnd, ku_bnd
#endif

#ifndef LIBRARY
      integer, parameter :: nbndvar=max(1,nfacevar)
      integer, parameter :: nbndvare=max(1,nvaredge)
      integer, parameter :: nbndvarc=max(1,nvarcorn)
      integer, parameter :: maxblocksf= 1+(maxblocks-1)                 &
     &                                             *min(1,nfacevar)
      integer, parameter :: maxblocksue= 1+(maxblocks-1)                &
     &                                             *min(1,nvaredge)
      integer, parameter :: maxblocksn= 1+(maxblocks-1)                 &
     &                                             *min(1,nvarcorn)
#else
      integer :: nbndvar
      integer :: nbndvare
      integer :: nbndvarc
      integer :: maxblocksf
      integer :: maxblocksue
      integer :: maxblocksn
#endif


! set data length of grid blocks
#ifndef LIBRARY
! cell center
      integer, parameter :: len_block=iu_bnd*ju_bnd*ku_bnd*nvar
! cell face centers
      integer, parameter :: len_blockfx=(iu_bnd+1)*ju_bnd*ku_bnd
      integer, parameter :: len_blockfy=iu_bnd*(ju_bnd+k2d)*ku_bnd
      integer, parameter :: len_blockfz=iu_bnd*ju_bnd*(ku_bnd+k3d)
! cell edge centers
      integer, parameter :: len_blockex=iu_bnd*(ju_bnd+k2d)*            &
     &                                         (ku_bnd+k3d)
      integer, parameter :: len_blockey=(iu_bnd+1)*ju_bnd*              &
     &                                         (ku_bnd+k3d)
      integer, parameter :: len_blockez=(iu_bnd+1)*(ju_bnd+k2d)*        &
     &                                         ku_bnd
! cell corner
      integer, parameter :: len_blockn=(iu_bnd+1)*(ju_bnd+k2d)*         &
     &                                         (ku_bnd+k3d)
! cell face centers for exchange of block faces only
      integer, parameter :: len_blockfxf=2*ju_bnd*ku_bnd
      integer, parameter :: len_blockfyf=iu_bnd*2*ku_bnd
      integer, parameter :: len_blockfzf=iu_bnd*ju_bnd*2
#else
! cell center
      integer :: len_block
! cell face centers
      integer :: len_blockfx 
      integer :: len_blockfy
      integer :: len_blockfz
! cell edge centers
      integer :: len_blockex
      integer :: len_blockey
      integer :: len_blockez
! cell corner
      integer :: len_blockn
! cell face centers for exchange of block faces only
      integer :: len_blockfxf
      integer :: len_blockfyf
      integer :: len_blockfzf
#endif


!
! Set the number of padded blocks required for the 1blk guardcell routines.
! This should be 2 in almost all circumstances, one block for the current
! working leaf node, and one for its parent.
      integer, parameter :: npblks=2


!
! Set index bounds with guardcells included.
#ifndef LIBRARY
      integer, parameter :: il_bnd1=1,iu_bnd1=nxb+2*nguard
      integer, parameter :: jl_bnd1=1,ju_bnd1=nyb+2*nguard*k2d
      integer, parameter :: kl_bnd1=1,ku_bnd1=nzb+2*nguard*k3d
#else
      integer :: il_bnd1,iu_bnd1
      integer :: jl_bnd1,ju_bnd1
      integer :: kl_bnd1,ku_bnd1
#endif

!
! Set length of messages required when working blocks are to be passed.
#ifndef LIBRARY
      integer, parameter :: len_block1 = iu_bnd1*ju_bnd1*ku_bnd1*nvar
      integer, parameter :: len_blockfx1 = (iu_bnd1+1)*ju_bnd1*ku_bnd1
      integer, parameter :: len_blockfy1 = iu_bnd1*(ju_bnd1+k2d)*ku_bnd1
      integer, parameter :: len_blockfz1 = iu_bnd1*ju_bnd1*(ku_bnd1+k3d)
      integer, parameter :: len_blockex1 = iu_bnd1*(ju_bnd1+k2d)        &
     &                                            *(ku_bnd1+k3d)
      integer, parameter :: len_blockey1 = (iu_bnd1+1)*ju_bnd1          &
     &                                            *(ku_bnd1+k3d)
      integer, parameter :: len_blockez1 = (iu_bnd1+1)*(ju_bnd1+1)      &
     &                                            *ku_bnd1
      integer, parameter :: len_blockn1 = (iu_bnd1+1)*(ju_bnd1+1)       &
     &                                            *(ku_bnd1+k3d)
      integer, parameter :: maxblocks_gt=(maxblocks-1)*(1-npgs)+1
      integer, parameter :: maxblocksf_gt=(maxblocksf-1)*(1-npgs)+1
      integer, parameter :: maxblocksue_gt=(maxblocksue-1)*(1-npgs)+1
      integer, parameter :: maxblocksn_gt=(maxblocksn-1)*(1-npgs)+1
#else
      integer :: len_block1
      integer :: len_blockfx1 
      integer :: len_blockfy1
      integer :: len_blockfz1 
      integer :: len_blockex1
      integer :: len_blockey1
      integer :: len_blockez1 
      integer :: len_blockn1
      integer :: maxblocks_gt
      integer :: maxblocksf_gt
      integer :: maxblocksue_gt
      integer :: maxblocksn_gt
#endif


!---------------------------------------
! Declare dimensions for the WORK arrays
!---------------------------------------

! Set number of guard cells associated with the workspace array.
      integer, parameter :: nguard_work = N_GUARD_CELLS_WORK


! Set number of variables which the workspace array can handle.
#ifndef LIBRARY
      integer, parameter :: nvar_work = N_VAR_WORK

      integer, parameter :: maxblocksw= 1+(maxblocks-1)                 &
     &                                             *min(1,nvar_work)
#else
      integer :: nvar_work
      integer :: maxblocksw
#endif

      integer, parameter :: ngw2=2*nguard_work

#ifndef LIBRARY
      integer, parameter :: ilw=1,iuw=nxb+ngw2*npgs
      integer, parameter :: jlw=1,juw=nyb+ngw2*npgs*k2d
      integer, parameter :: klw=1,kuw=nzb+ngw2*npgs*k3d
      integer, parameter :: len_wblock=iuw*juw*kuw
! Set index bounds with guardcells included.
      integer, parameter :: ilw1=1,iuw1=nxb+ngw2
      integer, parameter :: jlw1=1,juw1=nyb+ngw2*k2d
      integer, parameter :: klw1=1,kuw1=nzb+ngw2*k3d
      integer, parameter :: len_wblock1=iuw1*juw1*kuw1
#else
      integer :: ilw,iuw
      integer :: jlw,juw
      integer :: klw,kuw
      integer :: len_wblock
! Set index bounds with guardcells included.
      integer :: ilw1,iuw1
      integer :: jlw1,juw1
      integer :: klw1,kuw1
      integer :: len_wblock1
#endif



!-----------------
! Timestep control
!-----------------

! common block for timestep control
      integer, parameter :: maxlevels=50


#if defined(VAR_DT) || defined(PRED_CORR)
!      integer, parameter :: maxblocks_t=(maxblocks_alloc-1)*ivar_dt+1
!      integer, parameter :: nvar_t=(nvar-1)*ivar_dt+1

#endif





!---------------------------------
! Conservation at refinement jumps
!---------------------------------
#ifdef CONSV_FLUX_DENSITIES

#ifndef LIBRARY
      real, parameter :: red_f = 0.25
#else
      real :: red_f
#endif

#endif
#ifdef CONSV_FLUXES
!    if N_DIM = 3 set red_f = 1.0
!    if N_DIM = 2 set red_f = 0.5
!    if N_DIM = 1 set red_f = 0.25

#ifndef LIBRARY
#if N_DIM == 1
      real, parameter :: red_f = 0.25
#endif
#if N_DIM == 2
      real, parameter :: red_f = 0.5
#endif
#if N_DIM == 3
      real, parameter :: red_f = 1.0
#endif
#else
      real :: red_f
#endif

#endif
#ifndef CONSV_FLUX_DENSITIES
#ifndef CONSV_FLUXES
! set default value of red_f

#ifndef LIBRARY
      real, parameter :: red_f = 0.25
#else
      real :: red_f
#endif

#endif
#endif


!---------------------------------
! mpi message type control
!---------------------------------
#ifndef LIBRARY
#if N_DIM == 1
      integer, parameter :: nmax_lays = nxb/2
#endif
#if N_DIM == 2
      integer, parameter :: nmax_lays = min(nxb/2,nyb/2)
#endif
#if N_DIM == 3
      integer, parameter :: nmax_lays = min(nxb/2,nyb/2,nzb/2)
#endif
#else
      integer :: nmax_lays
#endif

!---------------------------------
! No. of fields for which divergence free prolongation
! will be required
      integer, parameter :: nfield_divf = NFIELD_DIVF

!---------------------------------


      end module paramesh_dimensions
!-----------------------------------------------------------------
