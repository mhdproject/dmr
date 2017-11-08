!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"


      subroutine amr_1blk_to_2copy( lcc,lfc,lec,lnc,lb,idest)




!------------------------------------------------------------------------
!
! This routine copies data from the 1-block working arrays with guardcells
! to the permanent data arrays, which may or may not have permanent
! guardcells, depending on whether NO_PERMANENT_GUARDCELLS is defined 
! in physicaldata.fh.
!
!
! Arguments :
!      lcc          logical       copies cell centered data if true
!      lfc          logical       copies cell face-centered data if true
!      lec          logical       copies cell edge-centered data if true
!      lnc          logical       copies cell corner data if true
!      lb           integer       block into which data is to be copied
!      idest        integer       sets value for last dimension index
!                                  in the 1-blk data arrays
!
!
! Written :     Peter MacNeice          February 1999
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace


#include "amr_shmem.fh"


!------------------------------------

      integer, intent(in) :: lb,idest
      logical, intent(in) :: lcc,lfc,lec,lnc

!------------------------------------


! cell-centered data
       if(lcc) then

#ifdef NO_PERMANENT_GUARDCELLS
           gt_unk(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,lb)
     .       = unk1(:,il_bnd+nguard:iu_bnd+nguard,
     .               jl_bnd+nguard*k2d:ju_bnd+nguard*k2d,
     .               kl_bnd+nguard*k3d:ku_bnd+nguard*k3d,idest)
#else
           gt_unk(:,:,:,:,lb) = unk1(:,:,:,:,idest)
#endif

       endif                             ! end of lcc if test



! cell face-centered data
       if(lfc) then

! x-face
#ifdef NO_PERMANENT_GUARDCELLS
         gt_facevarx(1:nfacevar,il_bnd:iu_bnd+1,
     .                       jl_bnd:ju_bnd,kl_bnd:ku_bnd,lb)
     .       = facevarx1(1:nfacevar,il_bnd+nguard:iu_bnd+nguard+1,
     .               jl_bnd+nguard*k2d:ju_bnd+nguard*k2d,
     .               kl_bnd+nguard*k3d:ku_bnd+nguard*k3d,idest)
#else
         gt_facevarx(1:nfacevar,1,:,:,lb) =
     .        facevarx1(1:nfacevar,1+nguard,:,:,idest)
         gt_facevarx(1:nfacevar,2,:,:,lb) =
     .        facevarx1(1:nfacevar,nxb+1+nguard,:,:,idest)
#endif

!!!#if N_DIM > 1
         if (ndim > 1) then
! y-face
#ifdef NO_PERMANENT_GUARDCELLS
         gt_facevary(1:nfacevar,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,
     .                         kl_bnd:ku_bnd,lb)
     .       = facevary1(1:nfacevar,il_bnd+nguard:iu_bnd+nguard,
     .               jl_bnd+nguard*k2d:ju_bnd+nguard*k2d+k2d,
     .               kl_bnd+nguard*k3d:ku_bnd+nguard*k3d,idest)
#else
         gt_facevary(1:nfacevar,:,1,:,lb) =
     .        facevary1(1:nfacevar,:,1+nguard*k2d,:,idest)
         gt_facevary(1:nfacevar,:,2,:,lb) =
     .        facevary1(1:nfacevar,:,nyb+(1+nguard)*k2d,:,idest)
#endif
!!!#endif
         end if

!!!#if N_DIM == 3
         if (ndim == 3) then
! z-face
#ifdef NO_PERMANENT_GUARDCELLS
         gt_facevarz(1:nfacevar,il_bnd:iu_bnd,jl_bnd:ju_bnd,
     .                         kl_bnd:ku_bnd+k3d,lb)
     .       = facevarz1(1:nfacevar,il_bnd+nguard:iu_bnd+nguard,
     .               jl_bnd+nguard*k2d:ju_bnd+nguard*k2d,
     .               kl_bnd+nguard*k3d:ku_bnd+nguard*k3d+k3d,idest)
#else
         gt_facevarz(1:nfacevar,:,:,1,lb) =
     .        facevarz1(1:nfacevar,:,:,1+nguard*k3d,idest)
         gt_facevarz(1:nfacevar,:,:,1+k3d,lb) =
     .        facevarz1(1:nfacevar,:,:,nzb+(1+nguard)*k3d,idest)
#endif
!!!#endif
         end if

        endif                     ! end of lfc if test


! cell edge-centered data
       if(lec) then
!!!#if N_DIM > 1
         if (ndim > 1) then
! x-edge
#ifdef NO_PERMANENT_GUARDCELLS
         gt_unk_e_x(1:nvaredge,il_bnd:iu_bnd,
     .               jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd+k3d,lb)
     .       = unk_e_x1(1:nvaredge,il_bnd+nguard:iu_bnd+nguard,
     .               jl_bnd+nguard*k2d:ju_bnd+(nguard+1)*k2d,
     .               kl_bnd+nguard*k3d:ku_bnd+(nguard+1)*k3d,idest)
#else
         gt_unk_e_x(1:nvaredge,:,:,:,lb) =
     .                        unk_e_x1(1:nvaredge,:,:,:,idest)
#endif
! y-edge
#ifdef NO_PERMANENT_GUARDCELLS
         gt_unk_e_y(1:nvaredge,il_bnd:iu_bnd+1,
     .               jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d,lb)
     .       = unk_e_y1(1:nvaredge,il_bnd+nguard:iu_bnd+1+nguard,
     .               jl_bnd+nguard*k2d:ju_bnd+nguard*k2d,
     .               kl_bnd+nguard*k3d:ku_bnd+(nguard+1)*k3d,idest)
#else
         gt_unk_e_y(1:nvaredge,:,:,:,lb) =
     .                        unk_e_y1(1:nvaredge,:,:,:,idest)
#endif
!!!#if N_DIM == 3
         if (ndim == 3) then
! z-edge
#ifdef NO_PERMANENT_GUARDCELLS
         gt_unk_e_z(1:nvaredge,il_bnd:iu_bnd+1,
     .               jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd,lb)
     .       = unk_e_z1(1:nvaredge,il_bnd+nguard:iu_bnd+1+nguard,
     .               jl_bnd+nguard*k2d:ju_bnd+(nguard+1)*k2d,
     .               kl_bnd+nguard*k3d:ku_bnd+nguard*k3d,idest)
#else
         gt_unk_e_z(1:nvaredge,:,:,:,lb) =
     .                        unk_e_z1(1:nvaredge,:,:,:,idest)
#endif
!!!#endif /* N_DIM == 3 */
         end if

!!!#endif /* N_DIM > 1 */
         end if

        endif                     ! end of lec if test

! cell corner data
       if(lnc) then
#ifdef NO_PERMANENT_GUARDCELLS
         gt_unk_n(1:nvarcorn,il_bnd:iu_bnd+1,
     .               jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd+k3d,lb)
     .       = unk_n1(1:nvarcorn,il_bnd+nguard:iu_bnd+1+nguard,
     .               jl_bnd+nguard*k2d:ju_bnd+(nguard+1)*k2d,
     .               kl_bnd+nguard*k3d:ku_bnd+(nguard+1)*k3d,idest)
#else
         gt_unk_n(1:nvarcorn,:,:,:,lb) =
     .                        unk_n1(1:nvarcorn,:,:,:,idest)
#endif
        endif                     ! end of lnc if test


      return
      end subroutine amr_1blk_to_2copy
