!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------


#include "paramesh_preprocessor.fh"
!#define DEBUG
      subroutine mpi_unpack_blocks(mype,iopt,
     .                             lcc,lfc,lec,lnc,
     .                             buf_dim,R_buffer,
     .                             nlayersx,nlayersy,nlayersz)

!------------------------------------------------------------------------
!
! This subroutine unpacks all blocks which are to be received on mype.
! It further stores the local (receiving) block id, the neighboring remote 
! (sending) block id, and the local guard block id into the array laddress 
! which is to be used in the subroutine mpi_1blk_guardcell.
!
!
! Written :     Maharaj Bhat & Michael Gehmeyr          March 2000
!------------------------------------------------------------------------
!
! Arguments:
!      mype           current processor id
!      iopt           option setting for work array
!      lcc            if true include unk data in buffer
!      lfc            if true include facevar data in buffer
!      lec            if true include unk_e_? data in buffer
!      lnc            if true include unk_n data in buffer
!      buf_dim        dimension of buffer
!      R_buffer       receive buffer 
!      ir_buf         starting & ending indices of buffer
!
!------------------------------------------------------------------------
      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton

      use paramesh_mpi_interfaces, only : mpi_put_buffer

      implicit none

      include 'mpif.h'

      integer, intent(in) :: mype,buf_dim,iopt
      logical, intent(in) :: lcc,lfc,lec,lnc
      real,    intent(inout) ::  R_buffer(buf_dim)
      integer, intent(in), optional :: nlayersx,nlayersy,nlayersz


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local variables

      integer :: lblk, lnumb, lb
      integer :: index
      integer :: ierrorcode,ierr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      lnumb = sum(commatrix_recv(:))
      if(lnumb.gt.maxblocks_alloc) then
            call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif
      index = 1

 
      do lblk=1,lnumb
        lb = lblk + strt_buffer - 1

                                  ! unpack all arrays from buffer into lb
#ifdef DEBUG
        write(*,*) 'pe ',mype,' lblk ',lblk,' unpacking starting ',
     .        ' at index ',index,' buf_dim ',buf_dim
#endif /* DEBUG */
        call mpi_put_buffer(
     .         lb,iopt,index,lcc,lfc,lec,lnc,buf_dim,R_buffer,
     .         nlayersx,nlayersy,nlayersz)
#ifdef DEBUG
        write(*,*) 'pe ',mype,' lblk ',lblk,' unpacked into ',lb
#endif /* DEBUG */

      enddo

#undef DEBUG
      return
      end subroutine mpi_unpack_blocks
