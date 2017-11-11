!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"

!------------------
! Guardcell filling
!#define DEBUG

      subroutine mpi_amr_write_guard_comm(nprocs)

      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton

      implicit none
      include 'mpif.h'

      integer, intent(in) :: nprocs

      integer             :: i,istat,ierrorcode,ierr
       
      if(.not.allocated(pe_source_guard)) 
     .   allocate(pe_source_guard(1:nprocs),stat = istat) 

      if(allocated(to_be_received_guard)) 
     .                 deallocate(to_be_received_guard)
      if(allocated(to_be_received)) then
        allocate(to_be_received_guard(3,size(to_be_received,2),
     .         size(to_be_received,3)),stat = istat) 
      else
        write(*,*)'Paramesh error: write_guard ',
     .            'to_be_received is not allocated.'
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif
      to_be_received_guard = to_be_received

      if(allocated(to_be_sent_guard)) 
     .                 deallocate(to_be_sent_guard)
      if(allocated(to_be_sent)) then
        allocate(to_be_sent_guard(3,size(to_be_sent,2),
     .         size(to_be_sent,3)),stat = istat) 
      else
        write(*,*)'Paramesh error: write_guard ',
     .            'to_be_sent is not allocated.'
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif
      to_be_sent_guard = to_be_sent

      if(allocated(commatrix_guard)) deallocate(commatrix_guard)
      allocate(commatrix_guard(1:nprocs,2),
     .           stat = istat) 
       if(istat.ne.0) then
         write(*,*) 'mpi_amr_write_guard_comm: allocation error'
         call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
       endif

      pe_source_guard = pe_source
      commatrix_guard = 0

      do i = 1,nprocs
         commatrix_guard(i,1) = commatrix_recv(i)
         commatrix_guard(i,2) = commatrix_send(i)
      enddo

      strt_guard = strt_buffer
      largest_no_of_blocks_guard = largest_no_of_blocks
      max_no_to_send_guard = max_no_to_send

      laddress_guard = laddress

      return
      end subroutine mpi_amr_write_guard_comm





      subroutine mpi_amr_read_guard_comm(nprocs)

      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton

      implicit none
      include 'mpif.h'

      integer, intent(in) :: nprocs

      integer             :: i,istat,ierrorcode,ierr

!-
      if(grid_analysed_mpi.ne.1) then
        write(*,*) 'PARAMESH ERROR: communication control info is ',
     .       'being read, but it was never set up. You are ',
     .       'probably missing a call to amr_checkpoint_re or ',
     .       'amr_refine_derefine.'
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif
!-

      pe_source = pe_source_guard
      commatrix_recv(:) = 0
      commatrix_send(:) = 0

      do i = 1,nprocs
        commatrix_recv(i) = commatrix_guard(i,1)
        commatrix_send(i) = commatrix_guard(i,2)
      enddo

      strt_buffer = strt_guard
      largest_no_of_blocks = largest_no_of_blocks_guard
      max_no_to_send = max_no_to_send_guard

      laddress = laddress_guard


      if(largest_no_of_blocks.gt.0.and.max_no_to_send.gt.0) then

      if(allocated(to_be_received)) 
     .                 deallocate(to_be_received)
      allocate(to_be_received(3,size(to_be_received_guard,2),
     .         size(to_be_received_guard,3)),stat = istat) 
      to_be_received = to_be_received_guard

      if(allocated(to_be_sent)) deallocate(to_be_sent)
      allocate(
     .  to_be_sent(3,largest_no_of_blocks,max_no_to_send),
     .           stat = istat) 
       if(istat.ne.0) then
         write(*,*) 'mpi_amr_read_guard_comm: allocation error'
         call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
       endif

      to_be_sent(1:3,1:largest_no_of_blocks,1:max_no_to_send)=
     .    to_be_sent_guard(1:3,1:largest_no_of_blocks,
     .                     1:max_no_to_send)

      endif


      mpi_pattern_id = 10

      call amr_1blk_guardcell_reset

      return
      end subroutine mpi_amr_read_guard_comm


!-----------------
! Prolongation
  
      subroutine mpi_amr_write_prol_comm(nprocs)

      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton

      implicit none
      include 'mpif.h'

      integer, intent(in) :: nprocs

      integer             :: i,istat,ierrorcode,ierr


      if(.not.allocated(pe_source_prol))
     .      allocate(pe_source_prol(1:nprocs),stat = istat)

      if(allocated(to_be_received_prol))
     .                 deallocate(to_be_received_prol)
      if(allocated(to_be_received)) then
       allocate(to_be_received_prol(3,size(to_be_received,2),
     .         size(to_be_received,3)),stat = istat)
      else
        write(*,*)'Paramesh error: write_prol ',
     .            'to_be_received is not allocated.'
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif
      to_be_received_prol = to_be_received

      if(allocated(to_be_sent_prol))
     .                 deallocate(to_be_sent_prol)
      if(allocated(to_be_sent)) then
       allocate(to_be_sent_prol(3,size(to_be_sent,2),
     .         size(to_be_sent,3)),stat = istat)
      else
        write(*,*)'Paramesh error: write_prol ',
     .            'to_be_sent is not allocated.'
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif
      to_be_sent_prol = to_be_sent

      if(allocated(commatrix_prol)) deallocate(commatrix_prol)
      allocate(commatrix_prol(1:nprocs,2),
     .           stat = istat)
       if(istat.ne.0) then
         write(*,*) 'mpi_amr_write_prol_comm: allocation error'
         call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
       endif

      pe_source_prol = pe_source

      commatrix_prol = 0

      do i = 1 ,nprocs
         commatrix_prol(i,1) = commatrix_recv(i)
         commatrix_prol(i,2) = commatrix_send(i)
      enddo

      strt_prol = strt_buffer
      largest_no_of_blocks_prol = largest_no_of_blocks
      max_no_to_send_prol = max_no_to_send

      laddress_prol = laddress

      return
      end subroutine mpi_amr_write_prol_comm





      subroutine mpi_amr_read_prol_comm(nprocs)

      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton

      implicit none
      include 'mpif.h'

      integer, intent(in) :: nprocs

      integer             :: i,istat,ierrorcode,ierr

!-
      if(grid_analysed_mpi.ne.1) then
        write(*,*) 'PARAMESH ERROR: communication control info is ',
     .       'being read, but it was never set up. You are ',
     .       'probably missing a call to amr_checkpoint_re or ',
     .       'amr_refine_derefine.'
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif
!-

      pe_source = pe_source_prol

      commatrix_recv(:) = 0
      commatrix_send(:) = 0

      do i = 1,nprocs
        commatrix_recv(i) = commatrix_prol(i,1)
        commatrix_send(i) = commatrix_prol(i,2)
      enddo

      strt_buffer = strt_prol
      largest_no_of_blocks = largest_no_of_blocks_prol
      max_no_to_send = max_no_to_send_prol

      laddress = laddress_prol

      if(largest_no_of_blocks.gt.0.and.max_no_to_send.gt.0) then

      if(allocated(to_be_received))
     .                 deallocate(to_be_received)
      if(allocated(to_be_received_prol)) then
        allocate(to_be_received(3,size(to_be_received_prol,2),
     .         size(to_be_received_prol,3)),stat = istat)
      else
        write(*,*)'Paramesh error: to_be_received_prol ',
     .            'is not allocated.'
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif
      to_be_received = to_be_received_prol

      if(allocated(to_be_sent)) deallocate(to_be_sent)
      allocate(
     .  to_be_sent(3,largest_no_of_blocks,max_no_to_send),
     .           stat = istat) 
       if(istat.ne.0) then
         write(*,*) 'mpi_amr_read_prol_comm: allocation error'
         call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
       endif

      to_be_sent(1:3,1:largest_no_of_blocks,1:max_no_to_send)=
     .    to_be_sent_prol(1:3,1:largest_no_of_blocks,
     .                     1:max_no_to_send)

      endif

      mpi_pattern_id = 20

      return
      end subroutine mpi_amr_read_prol_comm
 

 
!------------------
! Flux Conservation
  
      subroutine mpi_amr_write_flux_comm(nprocs)

      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton

      implicit none
      include 'mpif.h'

      integer, intent(in) :: nprocs

      integer             :: i,istat,ierrorcode,ierr


      if(.not.allocated(pe_source_flux))
     .      allocate(pe_source_flux(1:nprocs),stat = istat)

      if(allocated(to_be_received_flux))
     .                 deallocate(to_be_received_flux)
      if(allocated(to_be_received)) then
       allocate(to_be_received_flux(3,size(to_be_received,2),
     .         size(to_be_received,3)),stat = istat)
      else
        write(*,*)'Paramesh error: write_flux ',
     .            'to_be_received is not allocated.'
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif
      to_be_received_flux = to_be_received

      if(allocated(to_be_sent_flux))
     .                 deallocate(to_be_sent_flux)
      if(allocated(to_be_sent)) then
       allocate(to_be_sent_flux(3,size(to_be_sent,2),
     .         size(to_be_sent,3)),stat = istat)
      else
        write(*,*)'Paramesh error: write_flux ',
     .            'to_be_sent is not allocated.'
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif
      to_be_sent_flux = to_be_sent

      if(allocated(commatrix_flux)) deallocate(commatrix_flux)
      allocate(commatrix_flux(1:nprocs,2),
     .           stat = istat) 
       if(istat.ne.0) then
         write(*,*) 'mpi_amr_write_flux_comm: allocation error'
         call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
       endif

      pe_source_flux = pe_source

      commatrix_flux = 0

      do i = 1 ,nprocs
         commatrix_flux(i,1) = commatrix_recv(i)
         commatrix_flux(i,2) = commatrix_send(i)
      enddo

      strt_flux = strt_buffer
      largest_no_of_blocks_flux = largest_no_of_blocks
      max_no_to_send_flux = max_no_to_send

      laddress_flux = laddress


      return
      end subroutine mpi_amr_write_flux_comm





      subroutine mpi_amr_read_flux_comm(nprocs)

      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton

      implicit none
      include 'mpif.h'

      integer, intent(in) :: nprocs

      integer             :: i,istat,ierrorcode,ierr

!-
      if(grid_analysed_mpi.ne.1) then
        write(*,*) 'PARAMESH ERROR: communication control info is ',
     .       'being read, but it was never set up. You are ',
     .       'probably missing a call to amr_checkpoint_re or ',
     .       'amr_refine_derefine.'
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif
!-

      pe_source = pe_source_flux

      commatrix_recv(:) = 0
      commatrix_send(:) = 0

      do i = 1,nprocs
        commatrix_recv(i) = commatrix_flux(i,1)
        commatrix_send(i) = commatrix_flux(i,2)
      enddo

      strt_buffer = strt_flux
      largest_no_of_blocks = largest_no_of_blocks_flux
      max_no_to_send = max_no_to_send_flux

      laddress = laddress_flux

      if(largest_no_of_blocks.gt.0.and.max_no_to_send.gt.0) then

      if(allocated(to_be_received))
     .                 deallocate(to_be_received)
      if(allocated(to_be_received_flux)) then
        allocate(to_be_received(3,size(to_be_received_flux,2),
     .         size(to_be_received_flux,3)),stat = istat)
      else
        write(*,*)'Paramesh error: to_be_received_flux ',
     .            'is not allocated.'
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif
      to_be_received = to_be_received_flux

      if(allocated(to_be_sent)) deallocate(to_be_sent)
      allocate(
     .  to_be_sent(3,largest_no_of_blocks,max_no_to_send),
     .           stat = istat) 
       if(istat.ne.0) then
         write(*,*) 'mpi_amr_read_flux_comm: allocation error'
         call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
       endif

      to_be_sent(1:3,1:largest_no_of_blocks,1:max_no_to_send)=
     .    to_be_sent_flux(1:3,1:largest_no_of_blocks,
     .                     1:max_no_to_send)

      endif

      mpi_pattern_id = 30

      return
      end subroutine mpi_amr_read_flux_comm
 


!------------------
! Restriction

 
      subroutine mpi_amr_write_restrict_comm(nprocs)

      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton

      implicit none
      include 'mpif.h'

      integer, intent(in) :: nprocs

      integer             :: i,istat,ierrorcode,ierr

      if(.not.allocated(pe_source_restrict))
     .      allocate(pe_source_restrict(1:nprocs),
     .               stat = istat)
      pe_source_restrict(:) = pe_source(:)


      if(allocated(to_be_received)) then
         if (allocated(to_be_received_restrict)) then
            deallocate(to_be_received_restrict)
         end if
         allocate(to_be_received_restrict(3,
     .                                    size(to_be_received,2),
     .                                    size(to_be_received,3)),
     .                                    stat = istat)
         to_be_received_restrict(:,:,:) = to_be_received(:,:,:)
      else
        write(*,*)'Paramesh error: write_restrict ',
     .            'to_be_received is not allocated.'
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif

      if(allocated(to_be_sent)) then
         if (allocated(to_be_sent_restrict)) then
             deallocate(to_be_sent_restrict)
         end if
         allocate(to_be_sent_restrict(3,
     .                                size(to_be_sent,2),
     .                                size(to_be_sent,3)),
     .                                stat = istat)
         to_be_sent_restrict(:,:,:) = to_be_sent(:,:,:)
      else
        write(*,*)'Paramesh error: write_restrict ',
     .            'to_be_sent is not allocated.'
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif


      if(.not.allocated(commatrix_restrict)) then
         allocate(commatrix_restrict(1:nprocs,2),
     .            stat = istat)
         if(istat.ne.0) then
            write(*,*) 'mpi_amr_write_restrict_comm: allocation error'
            call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
         endif
      end if
      do i = 1 ,nprocs
         commatrix_restrict(i,:) = 0
         commatrix_restrict(i,1) = commatrix_recv(i)
         commatrix_restrict(i,2) = commatrix_send(i)
      enddo


      strt_restrict = strt_buffer
      largest_no_of_blocks_restrict = largest_no_of_blocks
      max_no_to_send_restrict = max_no_to_send

      laddress_restrict(:,:) = laddress(:,:)

      return
      end subroutine mpi_amr_write_restrict_comm





      subroutine mpi_amr_read_restrict_comm(nprocs)

      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton

      implicit none
      include 'mpif.h'

      integer, intent(in) :: nprocs

      integer             :: i,istat,ierrorcode,ierr

!-
      if(grid_analysed_mpi.ne.1) then
        write(*,*) 'PARAMESH ERROR: communication control info is ',
     .       'being read, but it was never set up. You are ',
     .       'probably missing a call to amr_checkpoint_re or ',
     .       'amr_refine_derefine.'
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif
!-

      pe_source(:) = pe_source_restrict(:)

      commatrix_recv(:) = 0
      commatrix_send(:) = 0

      do i = 1,nprocs
        commatrix_recv(i) = commatrix_restrict(i,1)
        commatrix_send(i) = commatrix_restrict(i,2)
      enddo

      strt_buffer = strt_restrict
      largest_no_of_blocks = largest_no_of_blocks_restrict
      max_no_to_send = max_no_to_send_restrict


      laddress(:,:) = laddress_restrict(:,:)

      if(largest_no_of_blocks.gt.0.and.max_no_to_send.gt.0) then


      if(allocated(to_be_received))
     .                 deallocate(to_be_received)         
      if(allocated(to_be_received_restrict)) then
         allocate(to_be_received(3,
     .            size(to_be_received_restrict,2),
     .            size(to_be_received_restrict,3)),
     .            stat = istat)
         to_be_received(:,:,:) = to_be_received_restrict(:,:,:)
      end if

      if(allocated(to_be_sent))
     .                 deallocate(to_be_sent)
      if(allocated(to_be_sent_restrict)) then
         allocate(to_be_sent(3,
     .                       largest_no_of_blocks,
     .                       max_no_to_send),
     .                       stat = istat)
         to_be_sent(:,1:largest_no_of_blocks,1:max_no_to_send) = 
     .    to_be_sent_restrict(:,1:largest_no_of_blocks,1:max_no_to_send)
      end if

      call amr_1blk_guardcell_reset

      endif

      mpi_pattern_id = 40

      return
      end subroutine mpi_amr_read_restrict_comm
