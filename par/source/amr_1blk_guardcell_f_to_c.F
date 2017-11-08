!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

      subroutine amr_1blk_guardcell_f_to_c(
     .                                  mype,pe,lb,iblock,iopt,nlayers,
     .                                  surrblks,lcc,lfc,lec,lnc,
     .                                  icoord,ldiag,
     .                                  nlayersx,nlayersy,nlayersz)


      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace

      use paramesh_interfaces, only : amr_1blk_guardcell_f_to_c_set,
     .                                amr_1blk_guardcell_f_to_c_fil

!------------------------------------------------------------------------
!
! This routine manages the exchange of guard cell information between
! blocks required to fill guard cells on block (pe,lb), assuming that 
! exchange is only required to the current block, from neighbor blocks
! at a higher refinement level.
! The actual exchanges are performed with calls to the routines 
! amr_1blk_cc_cp_remote and amr_1blk_fc_cp_remote, called from
! amr_1blk_guardcell_f_to_c_set.
! Once the 3 data layers are completely setup, the call to
! amr_1blk_guardcell_f_to_c_fil interpolates values to fill the
! guardcell data covered by any more refined neighbors of leaf blocks.
!
!
! Written :     Peter MacNeice          January 2003
!------------------------------------------------------------------------
!
! Arguments:
!      mype           local processor number
!      pe             processor address of the selected block
!      lb             local address on proc. pe of the selected block
!      iblock         selects the storage space in data_1blk.fh which is to
!                      be used in this call. If the leaf node is having its
!                      guardcells filled then set this to 1, if its parent
!                      is being filled set it to 2.
!      iopt           a switch to control which data source is to be used
!                      iopt=1 will use 'unk'
!                      iopt=2 will use 'work'
!      nlayers        the number of guard cell layers at each boundary
!      surrblks       the list of addresses of blocks surrounding block lb
!      lcc            a logical switch controlling whether unk or work data
!                      is filled
!      lfc            a logical switch controlling whether facevar data
!                      is filled
!      lec            a logical switch controlling whether unk_e_x(y)(z) data
!                      is filled
!      lnc            a logical switch controlling whether unk_n data
!                      is filled
!      icoord         an integer switch used to select which faces of
!                      the block are to be considered. If icoord=0 all
!                      faces are considered. If icoord=1 only faces perp.
!                      to the y-axis are considered, if icoord=2 only faces
!                      perp. to the x-axis are considered, and if icoord=3
!                      only faces perp. to the z-axis are considered.
!      ldiag          a logical switch which controls whether guardcells
!                      corresponding to neighbor blocks diagonally opposite
!                      block edges and corners are filled.
!
!------------------------------------

      implicit none

#include "amr_shmem.fh"


      integer, intent(in) :: mype,iopt,nlayers,lb,pe,iblock,icoord
      integer, intent(in) :: surrblks(:,:,:,:)
      logical, intent(in) :: lcc,lfc,lec,lnc,ldiag
      integer, intent(in) :: nlayersx,nlayersy,nlayersz

!------------------------------------
! local arrays

!------------------------------------
      if(ndim.gt.1) then                  ! ??? Is this restriction necessary??


! do not execute this routine if the current block is not a leaf block
      if(surrblks(3,2,2,2).ne.1) return


      l_f_to_c = .true.


      call amr_1blk_guardcell_f_to_c_set(
     .                                  mype,pe,lb,iblock,iopt,nlayers,
     .                                  surrblks,lcc,lfc,lec,lnc,
     .                                  icoord,ldiag,
     .                                  nlayersx,nlayersy,nlayersz)

      call amr_1blk_guardcell_f_to_c_fil(
     .                                  mype,pe,lb,iblock,iopt,nlayers,
     .                                  surrblks,lcc,lfc,lec,lnc,
     .                                  icoord,ldiag)



      endif                                 ! end of ndim > 1 if test
      l_f_to_c = .false.

!------------------------------------


      return
      end subroutine amr_1blk_guardcell_f_to_c
