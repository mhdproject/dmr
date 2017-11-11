#include "paramesh_preprocessor.fh"

      subroutine mpi_amr_boundary_block_info(mype,nprocs)


!-----------------------------------------------------------------
!
! This routine constructs a list of block faces which are external
! boundaries.
!
! Written : Peter MacNeice        November 2002
!
!-----------------------------------------------------------------
      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton

      implicit none

      include 'mpif.h'

      integer,intent(in) :: mype,nprocs

      integer,allocatable :: ib_global(:)
      integer,allocatable :: ib_count_send(:),ib_count_recv(:)
      integer,allocatable :: ib_global_recv(:)
      integer,allocatable :: ib_recv_index(:)
      integer,allocatable :: ib_send_index(:)

      integer :: ib,ib0,ib1,lb,ib_sum
      integer, allocatable :: ib_list(:,:)
      integer :: i,j,k,i1,i2,iproc,ierror

!-----------------------------------------------------------------
#ifdef MPI_USED

! Step 0.
! compute size for ib_list and then allocate it.
      ib0 = 0
      do lb = 1,lnblocks
        ib = minval(surr_blks(1,:,:,:,lb))
        if(ib.le.-20) then
          do k = 1,1+2*k3d
          do j = 1,1+2*k2d
          do i = 1,3
            ib1 = surr_blks(1,i,j,k,lb)
            if(ib1.le.-20) then
               ib0 = ib0 + 1
            end if
          enddo
          enddo
          enddo
        endif
      enddo
      if(allocated(ib_list)) deallocate(ib_list)
      allocate(ib_list(6,ib0))


! Step 1.
! construct on pe list of blocks next to boundary
      ib0 = 0
      do lb = 1,lnblocks
        ib = minval(surr_blks(1,:,:,:,lb))
        if(ib.le.-20) then
          do k = 1,1+2*k3d
          do j = 1,1+2*k2d
          do i = 1,3
            ib1 = surr_blks(1,i,j,k,lb)
            if(ib1.le.-20) then
            ib0 = ib0 + 1
            ib_list(1,ib0) = lb
            ib_list(2,ib0) = mype
            ib_list(3,ib0) = i
            ib_list(4,ib0) = j
            ib_list(5,ib0) = k
            ib_list(6,ib0) = ib
            endif
          enddo
          enddo
          enddo
        endif
      enddo

      if(.not.allocated(ib_global)) 
     .                         allocate(ib_global(0:nprocs-1))
      if(.not.allocated(ib_global_recv)) 
     .                         allocate(ib_global_recv(0:nprocs-1))
      if(.not.allocated(ib_count_send)) 
     .                         allocate(ib_count_send(0:nprocs-1))
      if(.not.allocated(ib_count_recv)) 
     .                         allocate(ib_count_recv(0:nprocs-1))
      if(.not.allocated(ib_recv_index)) 
     .                  allocate(ib_recv_index(0:nprocs-1))
      if(.not.allocated(ib_send_index)) 
     .                  allocate(ib_send_index(0:nprocs-1))
      ib_global_recv = 0
      ib_global = ib0

      call MPI_AlltoAll (ib_global     ,1,MPI_INTEGER,
     .                   ib_global_recv,1,MPI_INTEGER,
     .                   MPI_COMM_WORLD,ierror)
      ib_count_recv = ib_global_recv*6
      do iproc = 0,nprocs-1
        ib_count_send(iproc) = ib_count_recv(mype) 
      enddo

!
! Compute displacements in all to all message buffer, in bytes.
      ib_recv_index(0) = 0
      if(nprocs.gt.1) then
      do iproc = 1,nprocs-1
        ib_recv_index(iproc) = ib_recv_index(iproc-1) 
     .                        + ib_global_recv(iproc-1)*6
      enddo
      endif
      do iproc = 0,nprocs-1
        ib_send_index(iproc) = ib_recv_index(mype) 
      enddo



! Step 2.
! exchange no of boundary blocks on each processor
      call comm_int_sum_to_all(ib_sum,ib0)
      bc_block_neighs_length = ib_sum

! Step 3.
! allocate storage for global list
!!!      call mpi_barrier(MPI_COMM_WORLD,ierror)
      if(allocated(bc_block_neighs)) deallocate(bc_block_neighs)
      allocate(bc_block_neighs(6,ib_sum))
      if(allocated(bc_block_neighs_send)) 
     .                               deallocate(bc_block_neighs_send)
      allocate(bc_block_neighs_send(6,ib_sum))
!!!      call mpi_barrier(MPI_COMM_WORLD,ierror)


! Step 4.
! exchange info between procs

! Put local data into is correct place on the list
      bc_block_neighs_send = -1
      bc_block_neighs = -5
      i1 = ib_recv_index(mype)/6 + 1
      i2 = i1 + ib0 - 1
      
      if(ib0.gt.0) bc_block_neighs_send(:,i1:i2) = ib_list(:,1:ib0)


!!!      call mpi_barrier(MPI_COMM_WORLD,ierror)

! Exchange lists between all procs
      call MPI_AlltoAllv (bc_block_neighs_send,ib_count_send,
     .                                    ib_send_index,MPI_INTEGER,
     .                    bc_block_neighs,ib_count_recv,
     .                                    ib_recv_index,MPI_INTEGER,
     .                    MPI_COMM_WORLD,ierror)



      deallocate(ib_list)
      deallocate(ib_global)
      deallocate(ib_global_recv)
      deallocate(ib_count_send)
      deallocate(ib_count_recv)
      deallocate(ib_recv_index)
      deallocate(ib_send_index)

!
! Set status flag
      bc_block_neighs_status = 100

!!!      call mpi_barrier(MPI_COMM_WORLD,ierror)

#endif /* MPI_USED */
!-----------------------------------------------------------------

      return
      end subroutine mpi_amr_boundary_block_info
