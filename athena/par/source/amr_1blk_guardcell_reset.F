!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

      subroutine amr_1blk_guardcell_reset




!------------------------------------------------------------------------
!
! This routine resets some variables which manage the guardcell
! filling operation when operating in 1-block mode. It should be called
! from amr_initialize, and at a synchronization point which separates
! guardcell filling of data with a given time stamp, and guardcell
! filling at the next time stamp.
!
! Written :     Peter MacNeice          February 1999
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata

      implicit none

!------------------------------------

!
! reset id recording parent block which is currently in cache
          lnew_parent = .true.
          pcache_blk_u = -1
          pcache_pe_u  =  -1
          pcache_blk_w = -1
          pcache_pe_w  =  -1


      return
      end subroutine amr_1blk_guardcell_reset
