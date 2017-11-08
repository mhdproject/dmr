!#define DEBUG
#define NEWNEW

!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"


      subroutine amr_mpi_find_blk_in_buffer(
     .        mype,remote_block,remote_pe,idest,dtype,index,lfound)



!------------------------------------------------------------------------
!
! This routine finds where data for a remote block with address
! (remote_block,remote_pe) is in the recv buffer. It returns
! the message type, dtype, and the address in recv of the data word
! preceeding the first word of real data.
!
! Arguments:
!      mype             local processor
!      remote_pe0       remote processor
!      remote_block0    remote block 
!      dtype            message type - an integer between 1 and 27
!                       indicating the section of a blck contained in
!                       the message segment from (remote_block,remote_pe).
!      index            the address index+1 is where you should start
!                       reading the data from this message.
!
!
!
! Written :     Peter MacNeice          May 2001
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
#ifdef MPI_USED
      use mpi_morton
#endif /* MPI_USED */

      implicit none

#include "amr_shmem.fh"
      include 'mpif.h'


!-------------------------

      integer, intent(in)  :: mype,remote_pe,remote_block,idest
      integer, intent(out) :: dtype,index
      logical, intent(out) :: lfound

!-------------------------
! local arrays

      integer :: jseg,seg_no,iaddress,no_of_comms,jpe,jpe0
      integer :: ierrorcode,ierr,no_of_segments
      integer :: rem_pe,rem_blk,seg_offset
#ifdef NEWNEW
      integer :: iseg_no,jj
      logical :: llfound
#endif /* NEWNEW */


!-------------------------


#ifdef MPI_USED
      if(remote_pe.ne.mype) then

        rem_blk = remote_block
        rem_pe  = remote_pe

#ifdef NEWNEW
        llfound = .false.
        jj = ladd_strt(rem_pe)
        do while(.not.llfound.and.jj.le.ladd_end(rem_pe))
          if(rem_blk.eq.laddress(1,jj).and.
     .       rem_pe.eq.laddress(2,jj)) then
            llfound = .true.
            iseg_no = jj - strt_buffer + 1
          else
            jj = jj+1
          endif
        enddo
#endif /* NEWNEW */

      elseif(remote_pe.eq.mype.and.remote_block.gt.lnblocks) then

        rem_blk = laddress(1,remote_block)
        rem_pe  = laddress(2,remote_block)

#ifdef NEWNEW
        iseg_no = remote_block - strt_buffer + 1
#endif /* NEWNEW */
      endif

      if(rem_pe.ne.mype) then

#ifdef DEBUG
! locate rem_pe in the list of sending processors
        no_of_comms = size(pe_source)
        lfound = .false.
        jpe  = 0
        jpe0 = 0
        do while((.not.lfound).and.(jpe0.lt.no_of_comms))
          jpe0 = jpe0+1
          if(rem_pe.eq.pe_source(jpe0)-1) then
            lfound = .true.
            jpe = jpe0
          endif
        enddo
! If rem_pe is not located stop with error message
        if(jpe.eq.0) then
          if(idest.eq.2) return
          write(*,*) 'Paramesh error : pe ',mype,
     .     ' pe address of required data is not in the list of ',
     .     'communicating pes. ',
     .     ' remote_block ',remote_block,
     .     ' remote_pe ',remote_pe,
     .     ' rem_pe ',rem_pe,
     .     ' laddress ',laddress
          call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
        endif


        no_of_segments = size(to_be_received,2)
        lfound = .false.
        jseg = 0
        seg_no = jseg
        do while((.not.lfound).and.(jseg.lt.no_of_segments))
          jseg = jseg+1
          if(to_be_received(1,jseg,jpe).eq.rem_blk .and.
     .       to_be_received(2,jseg,jpe).eq.rem_pe+1 ) then
            lfound = .true.
            seg_no = jseg
          endif
        enddo

! Now compute where the list of segments from proc rem_pe actually begins
! in the complete list of message segments received on this processor
        seg_offset = 0
        if(rem_pe.gt.0) seg_offset = sum(commatrix_recv(1:rem_pe))
        seg_no = seg_no + seg_offset
#ifdef DEBUGX
        if(seg_no.ne.iseg_no) then
          write(*,*) 'seg_no and iseg_no are different ',
     .                seg_no,iseg_no
          call amr_abort()
        endif
#endif /* DEBUGX */
#endif /* DEBUG */

#ifdef NEWNEW
        seg_no = iseg_no
        lfound = llfound
#endif /* NEWNEW */

! If the requested segment is not located stop with error message
        if(seg_no.eq.0) then
          if(idest.eq.2) return
          write(*,*) 'Paramesh error : ',
     .     'message segment required is not in the list of ',
     .     'segments received.: proc ',mype,' to_be_received ',
     .     to_be_received(:,:,jpe),' mpi_pattern_id ',mpi_pattern_id
          call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
        endif
! set start address for this segment in R_buffer
        iaddress = mess_segment_loc(seg_no)

! Read out message into appropriate part of unk1 or work1
        dtype = anint(temprecv_buf(iaddress+2))

! We must write the message into recv first in case it is larger
! than the range we actually need.

        index = iaddress+2

! start address of data is now index+1

      endif
#endif /* MPI_USED */

      return
      end subroutine amr_mpi_find_blk_in_buffer
