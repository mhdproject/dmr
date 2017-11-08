!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"

      subroutine amr_1blk_save_soln





!------------------------------------------------------------------------
!
! This routine saves a global solution update into the time 
! synchronized global solution arrays, as is required when 
! using NO_PERMANENT_GUARDCELLS and the amr_1blk_guardcell routines.
!
! Written :     Peter MacNeice          May 1999
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree

      implicit none

#include "amr_shmem.fh"


!-------------------------

#ifdef NO_PERMANENT_GUARDCELLS

        if(nvar.gt.0) then
          unk(:,:,:,:,:) = gt_unk(:,:,:,:,:)
        endif
        if(nbndvar.gt.0) then
          facevarx(:,:,:,:,:) = gt_facevarx(:,:,:,:,:)
          if (ndim >= 2) then
          facevary(:,:,:,:,:) = gt_facevary(:,:,:,:,:)
          end if
          if (ndim == 3) then
          facevarz(:,:,:,:,:) = gt_facevarz(:,:,:,:,:)
          end if
        endif

        if (ndim > 1) then
        if(nvaredge.gt.0) then
          unk_e_x(:,:,:,:,:) = gt_unk_e_x(:,:,:,:,:)
          unk_e_y(:,:,:,:,:) = gt_unk_e_y(:,:,:,:,:)
          if (ndim == 3) then
          unk_e_z(:,:,:,:,:) = gt_unk_e_z(:,:,:,:,:)
          end if
        endif
        end if

        if(nvarcorn.gt.0) then
          unk_n(:,:,:,:,:) = gt_unk_n(:,:,:,:,:)
        endif

#endif

      return
      end subroutine amr_1blk_save_soln
