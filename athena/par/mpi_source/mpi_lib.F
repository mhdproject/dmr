!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"

!#define DEBUGX
!#define DEBUGY
!#define DEBUG



        subroutine comm_start(MaxProcs,nprocs,mype)

        use paramesh_mpi_interfaces, only : mpi_array_allocate

        implicit none

        integer, intent(out) :: nprocs,mype
        integer, intent(in)  :: MaxProcs

        include 'mpif.h'
        integer :: ierr
        integer, external :: shmem_n_pes,shmem_my_pe
        logical :: initialized

        call mpi_initialized(initialized, ierr)
        if (.not.initialized) then
          call mpi_init(ierr)
        end if

        nprocs = shmem_n_pes()
        mype = shmem_my_pe()

        call mpi_array_allocate(nprocs)

        return
        end subroutine comm_start

        subroutine comm_finish()
        use paramesh_mpi_interfaces, only : mpi_array_deallocate

        implicit none

        integer :: ierr

        call mpi_array_deallocate
        call mpi_finalize(ierr)

        return
        end subroutine comm_finish



        integer function shmem_n_pes()
        implicit none
        include 'mpif.h'
        integer ::  numb_processors, errcode

        call MPI_COMM_SIZE(MPI_COMM_WORLD, numb_processors, errcode)
        shmem_n_pes = numb_processors

        return
        end



        integer function shmem_my_pe()
        implicit none
        include 'mpif.h'
        integer ::  id, errcode

        call MPI_COMM_RANK(MPI_COMM_WORLD, id, errcode)
        shmem_my_pe = id
        return
        end


        subroutine shmem_udcflush()
        return
        end subroutine shmem_udcflush


        subroutine shmem_barrier_all()
        implicit none
        include 'mpif.h'
        integer ::  errcode

        call mpi_barrier (MPI_COMM_WORLD, errcode)

        return
        end subroutine shmem_barrier_all






        subroutine comm_logical_or_to_all(target,source)
        implicit none
        logical, intent(in)  :: source
        logical, intent(out) :: target
        integer :: nred, ierr
        common/comm_lcommon/ reduce_datain,reduce_dataout
        logical :: reduce_datain(1),reduce_dataout(1)

        include 'mpif.h'

        nred = 1
        reduce_datain(1) = source
        call mpi_logical_allreduce(reduce_datain(1),reduce_dataout(1),
     &                     nred,
     &                     MPI_LOGICAL,
     -                     MPI_LOR,MPI_COMM_WORLD,ierr)

        target = reduce_dataout(1)
        return
        end subroutine comm_logical_or_to_all

        subroutine comm_real_sum_to_all(target,source)
        implicit none
        real, intent(in)  :: source
        real, intent(out) :: target
        integer :: nred, ierr
        common/comm_rcommon/ reduce_datain,reduce_dataout
        real :: reduce_datain(1),reduce_dataout(1)

        include 'mpif.h'

        nred = 1
        reduce_datain(1) = source
        call mpi_real_allreduce(reduce_datain(1),reduce_dataout(1),
     &                     nred,
     &                     MPI_REAL,
     -                     MPI_SUM,MPI_COMM_WORLD,ierr)

        target = reduce_dataout(1)
        return
        end subroutine comm_real_sum_to_all


        subroutine comm_dble_sum_to_all(target,source)
        implicit none
        double precision, intent(in)  :: source
        double precision, intent(out) :: target
        integer :: nred, ierr
        common/comm_dcommon/ reduce_datain,reduce_dataout
        double precision :: reduce_datain(1),reduce_dataout(1)

        include 'mpif.h'

        nred = 1
        reduce_datain(1) = source
        call mpi_real_allreduce(reduce_datain(1),reduce_dataout(1),
     &                     nred,
     &                     MPI_DOUBLE_PRECISION,
     -                     MPI_SUM,MPI_COMM_WORLD,ierr)

        target = reduce_dataout(1)
        return
        end subroutine comm_dble_sum_to_all





        subroutine comm_real_min_to_all(target,source)
        implicit none
        real, intent(in)  :: source
        real, intent(out) :: target
        integer :: nred, ierr
        common/comm_rcommon/ reduce_datain,reduce_dataout
        real :: reduce_datain(1),reduce_dataout(1)

        include 'mpif.h'

        nred = 1
        reduce_datain(1) = source
        call mpi_real_allreduce(reduce_datain(1),reduce_dataout(1),
     &                     nred,
     -                     MPI_REAL,
     -                     MPI_MIN,MPI_COMM_WORLD,ierr)

        target = reduce_dataout(1)

        return
        end subroutine comm_real_min_to_all





        subroutine comm_real_max_to_all(target,source)
        implicit none
        real, intent(in)  :: source
        real, intent(out) :: target
        integer :: nred, ierr
        common/comm_rcommon/ reduce_datain,reduce_dataout
        real :: reduce_datain(1),reduce_dataout(1)

        include 'mpif.h'

        nred = 1
        reduce_datain(1) = source
        call mpi_real_allreduce(reduce_datain(1),reduce_dataout(1),
     &                     nred,
     -                     MPI_REAL,
     -                     MPI_MAX,MPI_COMM_WORLD,ierr)


        target = reduce_dataout(1)

        return
        end subroutine comm_real_max_to_all






        subroutine comm_int_sum_to_all(target,source)
        implicit none
        integer, intent(in)  :: source
        integer, intent(out) :: target
        integer :: nred, ierr
        common/comm_icommon/ ireduce_datain,ireduce_dataout
        integer :: ireduce_datain(1),ireduce_dataout(1)

        include 'mpif.h'

        nred = 1
        ireduce_datain(1) = source
        call mpi_int_allreduce(ireduce_datain(1),ireduce_dataout(1),
     &                  nred,
     -			MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
        target = ireduce_dataout(1)

        return
        end subroutine comm_int_sum_to_all



        subroutine comm_int_min_to_all(target,source)
        implicit none
        integer, intent(in)  :: source
        integer, intent(out) :: target
        integer :: nred, ierr
        common/comm_icommon/ ireduce_datain,ireduce_dataout
        integer :: ireduce_datain(1),ireduce_dataout(1)

        include 'mpif.h'

        nred = 1
        ireduce_datain(1) = source
        call mpi_int_allreduce(ireduce_datain(1),ireduce_dataout(1),
     &                  nred,
     -                  MPI_INTEGER,MPI_MIN,MPI_COMM_WORLD,ierr)

        target = ireduce_dataout(1)

        return
        end subroutine comm_int_min_to_all






        subroutine comm_int_max_to_all(target,source)
        implicit none
        integer, intent(in)  :: source
        integer, intent(out) :: target
        integer :: nred, ierr
        common/comm_icommon/ ireduce_datain,ireduce_dataout
        integer :: ireduce_datain(1),ireduce_dataout(1)

        include 'mpif.h'

        nred = 1
        ireduce_datain(1) = source
        call mpi_int_allreduce(ireduce_datain(1),ireduce_dataout(1),
     &                  nred,
     -                  MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
        target = ireduce_dataout(1)


        return
        end subroutine comm_int_max_to_all



      subroutine mpi_setup(mype,nprocs)


!------------------------------------------------------------------------
!
! This routine sets up some variables which are needed to control
! the MPI message passing. 
! This routine should be called every time the grid is rebuilt.
!
!
! Written :     Peter MacNeice and Michael Gehmeyr          February 2000
!------------------------------------------------------------------------
!
! Arguments:
!      mype           rank of local processor
!      nprocs         number of processors
!
!------------------------------------------------------------------------
      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton

      implicit none

      integer, intent(in) :: mype,nprocs

      include 'mpif.h'


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local variables

      integer :: lb
      logical :: l_periodx,l_periody,l_periodz
      integer :: ierr, ierror

      logical :: lreduce_datain(1),lreduce_dataout(1)

!
! Set up error handler
       call mpi_errhandler_set(MPI_COMM_WORLD,MPI_ERRORS_RETURN,
     .                         ierr)

! Detect if grid is periodic in any direction
       l_periodx = .false.
       l_periody = .false.
       l_periodz = .false.
       do lb = 1,lnblocks
         if( (bnd_box(1,1,lb)-.1*bsize(1,lb).lt.grid_xmin)
     .                           .and.
     .       (neigh(1,1,lb).gt.0) )
     .        l_periodx = .true.
         if( (bnd_box(1,2,lb)-.1*bsize(2,lb).lt.grid_ymin)
     .       .and. (neigh(1,3,lb).gt.0) .and. (ndim.ge.2) )
     .        l_periody = .true.
         if( (bnd_box(1,3,lb)-.1*bsize(3,lb).lt.grid_zmin)
     .       .and. (neigh(1,5,lb).gt.0) .and. (ndim.eq.3) )
     .        l_periodz = .true.
       enddo


       lreduce_datain(1) = l_periodx
       call mpi_logical_allreduce(lreduce_datain(1),lreduce_dataout(1),
     .                 1, MPI_LOGICAL,
     .                 MPI_LOR,MPI_COMM_WORLD,ierror)
       l_periodx = lreduce_dataout(1)
       lreduce_datain(1) = l_periody
       call mpi_logical_allreduce(lreduce_datain(1),lreduce_dataout(1),
     .                 1, MPI_LOGICAL,
     .                 MPI_LOR,MPI_COMM_WORLD,ierror)
       l_periody = lreduce_dataout(1)
       lreduce_datain(1) = l_periodz
       call mpi_logical_allreduce(lreduce_datain(1),lreduce_dataout(1),
     .                 1, MPI_LOGICAL,
     .                 MPI_LOR,MPI_COMM_WORLD,ierror)
       l_periodz = lreduce_dataout(1)


        lperiodicx = l_periodx
        lperiodicy = l_periody
        lperiodicz = l_periodz

!       write(*,*) 'pe ',mype,' lperiodicx ',lperiodicx
!       write(*,*) 'pe ',mype,' lperiodicy ',lperiodicy
!       write(*,*) 'pe ',mype,' lperiodicz ',lperiodicz


! Initialize logical array which records the type of data in the last
! message which was packed.
        l_datapacked = .false.

      return
      end subroutine mpi_setup


#ifdef MOVED
      subroutine send_real(rv,n,idest,itag,icom,ierr)

!------------------------------------------------------------------------
!
! This routine encapsulats the MPI_Send for real arrays.
!
!
! Written :     Maharaj Bhat          February 2000
!------------------------------------------------------------------------
!
! Arguments:
!      rv             real valued array
!      n              dimension of rv
!      idest          MPI destination of msg
!      itag           MPI tag of msg
!      icom           MPI communicator
!      ierr           MPI error flag
!
!------------------------------------------------------------------------

      implicit none
      include 'mpif.h'
    
      real rv(*)
      integer n,idest,itag,icom,ierr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call MPI_SEND(rv,n,MPI_REAL,
     .                   idest,itag,icom,ierr)

      return
      end subroutine send_real

      subroutine send_int(iv,n,idest,itag,icom,ierr)

!------------------------------------------------------------------------
!
! This routine encapsulats the MPI_Send for integer arrays.
!
!
! Written :     Maharaj Bhat          February 2000
!------------------------------------------------------------------------
!
! Arguments:
!      iv             integer valued array
!      n              dimension of iv
!      idest          MPI destination of msg
!      itag           MPI tag of msg
!      icom           MPI communicator
!      ierr           MPI error flag
!
!------------------------------------------------------------------------

      implicit none
      include 'mpif.h'
     
      integer iv(*)
      integer n,idest,itag,icom,ierr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call MPI_SEND(iv,n,MPI_INTEGER,
     .                   idest,itag,icom,ierr)

      return
      end subroutine send_int

      subroutine recv_real(rv,n,isource,itag,icom,istatus,ierr)

!------------------------------------------------------------------------
!
! This routine encapsulats the MPI_Send for integer arrays.
!
!
! Written :     Maharaj Bhat          February 2000
!------------------------------------------------------------------------
!
! Arguments:
!      rv             real valued array
!      n              dimension of rv
!      isource        MPI source of msg
!      itag           MPI tag of msg
!      icom           MPI communicator
!      istatus        MPI status flag
!      ierr           MPI error flag
!
!------------------------------------------------------------------------

      implicit none
      include 'mpif.h'
     
      real rv(*)
      integer n,isource,itag,icom,ierr
      integer istatus(MPI_STATUS_SIZE)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call MPI_RECV(rv,n,MPI_REAL,
     .                   isource,itag,icom,istatus,ierr)

      return
      end subroutine recv_real

      subroutine recv_int(iv,n,isource,itag,icom,istatus,ierr)

!------------------------------------------------------------------------
!
! This routine encapsulats the MPI_Send for integer arrays.
!
!
! Written :     Maharaj Bhat          February 2000
!------------------------------------------------------------------------
!
! Arguments:
!      iv             integer valued array
!      n              dimension of iv
!      isource        MPI source of msg
!      itag           MPI tag of msg
!      icom           MPI communicator
!      istatus        MPI status flag
!      ierr           MPI error flag
!
!------------------------------------------------------------------------

      implicit none
      include 'mpif.h'
     
      integer iv(*)
      integer n,isource,itag,icom,ierr
      integer istatus(MPI_STATUS_SIZE)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      call MPI_RECV(iv,n,MPI_INTEGER,
     .                   isource,itag,icom,istatus,ierr)

      return
      end subroutine recv_int


#endif /*  MOVED */

      subroutine mpi_array_allocate(nprocs)


!------------------------------------------------------------------------
!
! Some of the arrays which are used by the routines which manage MPI
! communications need to be dimensioned with sizes determined at
! run time. These arrays are declared as allocatable in  the
! module mpi_morton, and the actual allocation is performed here.
!
! Written :     Peter MacNeice          June 2000
!------------------------------------------------------------------------
!
! Arguments:
!      nprocs         number of processors
!
!------------------------------------------------------------------------
      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton


      integer, intent(in) :: nprocs


!------------------------------------------------------------------------

      allocate( pe_remote(1:nprocs)             )
      allocate( pe_source(1:nprocs)             )
      allocate( pe_destination(1:nprocs)        )
      allocate( morton_limits(6,1:2,1:2,1:nprocs) )
      allocate( commatrix_send(1:nprocs)        )
      allocate( commatrix_recv(1:nprocs)        )
      allocate( ir_buf(2,1:nprocs)              )
      allocate( is_buf(2,1:nprocs)              )


      return
      end subroutine mpi_array_allocate



      subroutine mpi_array_deallocate


!------------------------------------------------------------------------
!
! Some of the arrays which are used by the routines which manage MPI
! communications need to be dimensioned with sizes determined at
! run time. These arrays are declared as allocatable in the
! module mpi_morton. This routine deallocates them for a clean
! exit from the program.
!
! Written :     Peter MacNeice          June 2000
!------------------------------------------------------------------------
!
! Arguments:
!
!------------------------------------------------------------------------
      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton

!------------------------------------------------------------------------

      if(allocated(pe_remote)) deallocate( pe_remote     )
      if(allocated(pe_source)) deallocate( pe_source     )
      if(allocated(pe_destination)) deallocate( pe_destination)
      if(allocated(morton_limits)) deallocate( morton_limits )
      if(allocated(commatrix_send)) deallocate( commatrix_send)
      if(allocated(commatrix_recv)) deallocate( commatrix_recv)
      if(allocated(ir_buf)) deallocate( ir_buf        )
      if(allocated(is_buf)) deallocate( is_buf        )


      return
      end subroutine mpi_array_deallocate


      subroutine morton_number(
     .          x0,y0,z0,bbsize,ndim,lrefine_max,lrefine,
     .          mort)

!---------------------------------------------------------------------
!
! This subroutine computes a morton number for a given set of
! coordinates and refinement level.
!
! Written :     December 2002     Peter MacNeice
! Modified for new morton numbering scheme: March 2003  Kevin Olson
!---------------------------------------------------------------------

      integer,intent(in ) :: lrefine,lrefine_max,ndim
      real,intent(in)     :: bbsize(3),x0,y0,z0
      integer,intent(out) :: mort(6)
      logical, save :: first = .true.
      integer, save :: nbits,nbits2
!------------------------------------------------------------

      integer, parameter :: ikind = selected_int_kind(18)
      integer (kind=selected_int_kind(18)) :: ix, iy, iz

      if (first) then
       first = .false.
       nbits = bit_size(mort(1))
       if (nbits == 32) then
           nbits = nbits-2
       elseif (nbits == 64) then
           nbits = nbits-4
       end if
       nbits2 = bit_size(ix)-4
      end if

! compute morton number
      ix = int(x0/bbsize(1),kind=ikind)
      if (ndim >= 2) then
         iy = int(y0/bbsize(2),kind=ikind)
      else
         iy = 0
      end if
      if (ndim == 3) then
         iz = int(z0/bbsize(3),kind=ikind)
      else
         iz = 0
      end if
     
      mort(:) = 0

      ipos2 = 0
      j = 6
      do while (ipos2 < nbits2-1)

      ipos = 0
      do while (ipos < nbits-1)

         if (btest(ix,ipos2)) then
            mort(j) = ibset(mort(j),ipos)
         else
            mort(j) = ibclr(mort(j),ipos)
         end if
         ipos = ipos + 1
         if (ndim >= 2) then
            if (btest(iy,ipos2)) then
               mort(j) = ibset(mort(j),ipos)
            else
               mort(j) = ibclr(mort(j),ipos)
            end if
            ipos = ipos + 1
         end if
         if (ndim == 3) then
            if (btest(iz,ipos2)) then
               mort(j) = ibset(mort(j),ipos)
            else
               mort(j) = ibclr(mort(j),ipos)
            end if
            ipos = ipos + 1
         end if

         ipos2 = ipos2 + 1

      end do

      j = j - 1

      end do

c now shift bits to the left by max_levels - level

      nbitshft = ndim*(lrefine_max-lrefine)
      if(nbitshft.gt.0) then

        nshifts = 0
        do while (nbitshft > nbits)
           nbitshft = nbitshft - nbits
           nshifts = nshifts + 1
        end do

        do k = 1,nshifts
        mort(1) = ishft(mort(1),nbits)
        do j = 1,5
        ipos2 = 0
        do ipos = 0,nbits-1
           call mvbits (mort(j+1),ipos,1,mort(j),ipos2)
           mort(j+1) = ibclr(mort(j+1),ipos)
           ipos2 = ipos2 + 1
        end do
        mort(j+1) = ishft(mort(j+1),nbits)
        end do
        end do

        mort(1) = ishft(mort(1),nbitshft)
        do j = 1,5
        ipos2 = 0
        do ipos = nbits-1-(nbitshft-1),nbits-1
           call mvbits (mort(j+1),ipos,1,mort(j),ipos2)
           mort(j+1) = ibclr(mort(j+1),ipos)
           ipos2 = ipos2 + 1
        end do
        mort(j+1) = ishft(mort(j+1),nbitshft)
        end do

      end if

      return
      end subroutine morton_number


      subroutine morton_parent(mort,lrefine,lrefine_max,ndim,
     .                         morton_p)

!---------------------------------------------------------------------
!
! This function returns the morton number for the parent of the
! tree node MORT.
!
! Written :     April 1999     Peter MacNeice
!---------------------------------------------------------------------

      implicit none

      integer, intent(in) :: mort(6), lrefine, lrefine_max, ndim
      integer, intent(out) :: morton_p(6)

      integer :: lbit, nbits_zero, j, nbits_to_zero
      integer, save :: nbits
      logical, save :: first = .true.

      nbits_zero = (lrefine_max+1-lrefine)*ndim
      if (first) then
      nbits = bit_size(mort(1))
      if (nbits == 32) then
         nbits = nbits-2
      else
         nbits = nbits-4
      end if
      first = .false.
      end if

      j = 6
      morton_p(:) = mort(:)
      do while (nbits_zero > 0 .or. j == 6)
         nbits_to_zero = min(nbits,nbits_zero)
         do lbit = 0,nbits_to_zero-1
           morton_p(j) = ibclr(morton_p(j),lbit)
         end do
         nbits_zero = nbits_zero-nbits
         j = j - 1
      end do

      return
      end




      subroutine morton_children(mort,lrefine,lrefine_max,ndim,
     .                           mort_child,nchild)
!---------------------------------------------------------------------
! Computes the morton numbers for children of a node with given
! morton number, refinement level, and where the maximum refinement
! level in use is known.
!
! Arguments:
! mort          integer        morton number of current node
! lrefine       integer        refinement level of current node
! lrefine_max   integer        highest refinement level currently in use
! ndim          integer        physical dimension
! mort_child    integer array  array of morton numbers of this nodes 
!                                children
! nchild        integer        number of children of any node
! 
! Written :     April 1999     Peter MacNeice
!
!---------------------------------------------------------------------

      implicit none

      integer, intent(in)  :: lrefine,lrefine_max,ndim
      integer, intent(in)  :: mort(6), nchild
      integer, intent(out) :: mort_child(6,nchild)

      integer :: lbit,ichild
      integer :: i, j, jj
      integer, save :: nbits
      logical, save :: first = .true.

      if (first) then
      nbits = bit_size(mort(1))
      if (nbits == 32) then
         nbits = nbits - 2
      else
         nbits = nbits - 4
      end if
      first = .false.
      end if

      if (lrefine < lrefine_max) then
! locaton of least sig. bit of child
        mort_child(:,:) = 0
        lbit = (lrefine_max-lrefine-1)*ndim
        j = 6-int(lbit/nbits)
        do jj = 1,j-1
          mort_child(jj,:) = mort(jj)
        end do
        lbit = lbit-(6-j)*nbits
        do ichild = 0,nchild-1
          i = ishft(ichild,lbit)
          mort_child(j,ichild+1) = mort(j) + i
        enddo
        do jj = j+1,6
          mort_child(jj,:) = mort(jj)
        end do
      endif

      return
      end subroutine morton_children




      subroutine morton_neighbors(xmin,ymin,zmin,xmax,ymax,zmax,
     .                            lperiodicx,lperiodicy,lperiodicz,
     .                            coord,bbsize,ndim,
     .                            lrefine,lrefine_max,mort_neigh)

!---------------------------------------------------------------------
!
! This routine computes the morton numbers which all possible
! neighbors of the current block, at the same refinement level
! would have. If a neighbor would lie below the minimun coordinate
! value for any axis, the morton number returned defaults to -1.
!
! Arguments:
! xmin          real           lower x range of computational domain
! ymin          real           lower y range of computational domain
! zmin          real           lower z range of computational domain
! xmax          real           upper x range of computational domain
! ymax          real           upper y range of computational domain
! zmax          real           upper z range of computational domain
! lperiodicx    logical        .true. if periodic in x-direction, otherwise
!                                set to .false.
! lperiodicy    logical        .true. if periodic in y-direction, otherwise
!                                set to .false.
! lperiodicz    logical        .true. if periodic in z-direction, otherwise
!                                set to .false.
! coord         real array     array storing coordinates of current block
! bbsize        real array     array storing size of current block
! ndim          integer        physical dimension
! lrefine       integer        refinement level of current block
! lrefine_max   integer        highest refinement level currently in use
!                                in the calculation
! mort_neigh    integer array  array storing morton numbers of all
!                                this blocks possible neighbors at its
!                                own refinement level.
!
! Written :     April 1999     Peter MacNeice
!
!---------------------------------------------------------------------

      use paramesh_mpi_interfaces, only : morton_number

      implicit none

      real, intent(in)     :: xmin,ymin,zmin,xmax,ymax,zmax
      real, intent(in)     :: coord(3),bbsize(3)
      integer, intent(in)  :: ndim
      integer, intent(in)  :: lrefine,lrefine_max
      logical, intent(in)  :: lperiodicx,lperiodicy,lperiodicz
      integer, intent(out) :: mort_neigh(6,3,3,3)

      real :: bbsize_loc(3)
      integer :: i,j,k
      real    :: x0,y0,z0,xtest,ytest,ztest,test
      integer :: k2d, k3d
      integer :: mort(6)


      mort_neigh(:,:,:,:) = -1
      k2d = 1
      k3d = 1

      bbsize_loc(:) = bbsize(:)

      if(ndim.le.2) then
          k3d = 0
          bbsize_loc(3) = 1.
      endif
      if(ndim.eq.1) then
          k2d = 0
          bbsize_loc(2) = 1.
      endif


! cycle over possible neighbors
      do k = -k3d,k3d
        z0 = coord(3) + real(k)*bbsize_loc(3)  
        if(lperiodicz.and.z0.lt.zmin) z0 = z0 + (zmax-zmin)
        if(lperiodicz.and.z0.gt.zmax) z0 = z0 - (zmax-zmin)
        z0 = z0 - zmin   ! added
        if(ndim.lt.3) z0 = 0.
        do j = -k2d,k2d
          y0 = coord(2) + real(j)*bbsize_loc(2)
          if(lperiodicy.and.y0.lt.ymin) y0 = y0 + (ymax-ymin)
          if(lperiodicy.and.y0.gt.ymax) y0 = y0 - (ymax-ymin)
          y0 = y0 - ymin   ! added
          if(ndim.lt.2) y0 = 0.
          do i = -1,1
            x0 = coord(1) + real(i)*bbsize_loc(1)
            if(lperiodicx.and.x0.lt.xmin) x0 = x0 + (xmax-xmin)
            if(lperiodicx.and.x0.gt.xmax) x0 = x0 - (xmax-xmin)
            x0 = x0 - xmin   ! added

            call morton_number (
     .         x0,y0,z0,bbsize_loc,ndim,lrefine_max,lrefine,
     .         mort)

            mort_neigh(:,i+2,j+2,k+2) = mort

! if coord is below minimum coords set mort to -1 
            xtest = x0*(x0-(xmax-xmin))
            if(lperiodicx) xtest = -1.
            ytest = y0*(y0-(ymax-ymin))
            if(lperiodicy) ytest = -1.
            ztest = z0*(z0-(zmax-zmin))
            if(lperiodicz) ztest = -1.
            test = max(xtest,ytest)
            test = max( test,ztest)
            
            if(test.gt.0.) mort_neigh(:,i+2,j+2,k+2) = -1

          enddo
        enddo
      enddo

      return
      end subroutine morton_neighbors

      subroutine morton_child_neighbors(xmin,ymin,zmin,xmax,ymax,zmax,
     .                            lperiodicx,lperiodicy,lperiodicz,
     .                            coord,bbsize,ndim,
     .                            lrefine,lrefine_max,mort_neigh)

!---------------------------------------------------------------------
!
! This routine computes the morton numbers of all possible
! blocks one level more refined than the current block and
! bordering the current block.
!
! Arguments:
! xmin          real           lower x range of computational domain
! ymin          real           lower y range of computational domain
! zmin          real           lower z range of computational domain
! xmax          real           upper x range of computational domain
! ymax          real           upper y range of computational domain
! zmax          real           upper z range of computational domain
! lperiodicx    logical        .true. if periodic in x-direction, otherwise
!                                set to .false.
! lperiodicy    logical        .true. if periodic in y-direction, otherwise
!                                set to .false.
! lperiodicz    logical        .true. if periodic in z-direction, otherwise
!                                set to .false.
! coord         real array     array storing coordinates of current block
! bbsize        real array     array storing size of current block
! ndim          integer        physical dimension
! lrefine       integer        refinement level of current block
! lrefine_max   integer        highest refinement level currently in use
!                                in the calculation
! mort_neigh    integer array  array storing morton numbers of all
!                                this blocks possible neighbors at its
!                                own refinement level.
!
! Written :     November 2002     Peter MacNeice
!
!---------------------------------------------------------------------

      use paramesh_mpi_interfaces, only : morton_number

      implicit none

      real, intent(in)     :: xmin,ymin,zmin,xmax,ymax,zmax
      real, intent(in)     :: coord(3),bbsize(3)
      integer, intent(in)  :: ndim
      integer, intent(in)  :: lrefine,lrefine_max
      logical, intent(in)  :: lperiodicx,lperiodicy,lperiodicz
      integer, intent(out) :: mort_neigh(6,4,4,4)

      real    :: x0,y0,z0,xtest,ytest,ztest,test
      real    :: cbsize(3)
      integer :: i,j,k, k2d, k3d
      integer :: mort(6)
      integer :: lrefine1

      mort_neigh(:,:,:,:) = -1
      k2d = 1
      k3d = 1

      cbsize = bbsize*.5
      lrefine1 = lrefine + 1
      if(ndim.le.2) then
          k3d = 0
          cbsize(3) = 1.
      endif
      if(ndim.eq.1) then
          k2d = 0
          cbsize(2) = 1.
      endif

! cycle over possible neighbors
      do k = -2*k3d,1*k3d
        z0 = coord(3) + real(k3d)*(real(k+2)-1.5)*cbsize(3)  
        if(lperiodicz.and.z0.lt.zmin) z0 = z0 + (zmax-zmin)
        if(lperiodicz.and.z0.gt.zmax) z0 = z0 - (zmax-zmin)
        z0 = z0 - zmin   ! added
        if(ndim.lt.3) z0 = 0.
        do j = -2*k2d,1*k2d
          y0 = coord(2) + real(k2d)*(real(j+2)-1.5)*cbsize(2)  
          if(lperiodicy.and.y0.lt.ymin) y0 = y0 + (ymax-ymin)
          if(lperiodicy.and.y0.gt.ymax) y0 = y0 - (ymax-ymin)
          y0 = y0 - ymin   ! added
          if(ndim.lt.2) y0 = 0.
          do i = -2,1
            x0 = coord(1) + (real(i+2)-1.5)*cbsize(1)  
            if(lperiodicx.and.x0.lt.xmin) x0 = x0 + (xmax-xmin)
            if(lperiodicx.and.x0.gt.xmax) x0 = x0 - (xmax-xmin)
            x0 = x0 - xmin   ! added

            call morton_number (
     .         x0,y0,z0,cbsize,ndim,lrefine_max,lrefine1,
     .         mort)

            mort_neigh(:,i+3,1+(j+2)*k2d,1+(k+2)*k3d) = mort(:)

! if coord is below minimum coords set mort to -1 
            xtest = x0*(x0-(xmax-xmin))
            if(lperiodicx) xtest = -1.
            ytest = y0*(y0-(ymax-ymin))
            if(lperiodicy) ytest = -1.
            ztest = z0*(z0-(zmax-zmin))
            if(lperiodicz) ztest = -1.
            test = max(xtest,ytest)
            test = max( test,ztest)

            if(test.gt.0.) 
     .         mort_neigh(:,i+3,1+(j+2)*k2d,1+(k+2)*k3d) = -1

          enddo
        enddo
      enddo

! mask out center of nort_neigh
      mort_neigh(:,2:3,1+k2d:1+2*k2d,1+k3d:1+2*k3d) = -1

      return
      end subroutine morton_child_neighbors




      subroutine boundary_locator(x0,y0,z0,lboundary,ibc)




!------------------------------------------------------------------------
!
! This routine detects whether a given coordinate is on a particular
! boundary.
! The boundary regions are assumed to be volumes defined by their 
! bounding boxes, and each has a boundary condition index associated
! with it.
!
! Written :     Peter MacNeice         June 2000
!------------------------------------------------------------------------
!
! Arguments:
!
!------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree

      implicit none

      logical,intent(out) :: lboundary
      integer,intent(out) :: ibc
      real,intent(in)     :: x0,y0,z0

      real    :: xtest, ytest, ztest
      integer :: ibnd 
  

!------------------------------------
 
      lboundary = .false.
      do ibnd = 1, nboundaries

        xtest = (x0 - boundary_box(1,1,ibnd))*
     .          (x0 - boundary_box(2,1,ibnd))
        ytest = (y0 - boundary_box(1,2,ibnd))*
     .          (y0 - boundary_box(2,2,ibnd))
        ztest = (z0 - boundary_box(1,3,ibnd))*
     .          (z0 - boundary_box(2,3,ibnd))
        if(ndim.lt.3) ztest = -1.
        if(ndim.lt.2) ytest = -1.
        if(xtest.le.0. .and. ytest.le.0. and. ztest.le.0.) then
          lboundary = .true.
          ibc       = boundary_index(ibnd)
        endif

      enddo

      return
      end subroutine boundary_locator



      subroutine mpi_xchange_blocks(mype,nprocs, tag_offset,
     .                              buf_dim, S_buffer, R_buffer)

!------------------------------------------------------------------------
!
! This routine uses the general hand-shaking scheme for message passing 
! between pe's. 
! sends and receives buffers containing all requested block data.
!
!
! Written :     Maharaj Bhat & Michael Gehmeyr          March 2000
!------------------------------------------------------------------------
!
! Arguments:
!      mype           rank of local processor
!      nprocs         number of processors
!      buf_dim        dimension of buffers
!      tag_offset     offset for MPI tag
!      S_buffer       send buffer
!      R_buffer       recv buffer
!
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace
      use mpi_morton
 
      implicit none

      include 'mpif.h'

      integer, intent(in)    :: mype,nprocs,buf_dim
      integer, intent(inout) :: tag_offset
      real,    intent(inout) :: S_buffer(buf_dim), R_buffer(buf_dim)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local variables

      integer :: i,j
      integer :: istrt,ilast,isize
      integer :: isrc,idest, itag
      integer :: recvrequest(nprocs)
      integer :: recvstatus(MPI_STATUS_SIZE,nprocs)
      integer :: ierrorcode,ierr
      integer :: ij, ji


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!      R_buffer = 0

      if (commatrix_send(mype+1) > 0 
     .    .or. 
     .    commatrix_recv(mype+1) > 0) then
        write(*,*) 'Paramesh error :  error in xchange : pe ',mype,
     .      ' diagonal element of commatrix is non-zero ',
     .      commatrix_recv(mype+1),
     .      commatrix_send(mype+1)
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif

      ij = 0
      ji = 0

      do i = 1,nprocs
         isrc = i-1
         idest= mype
         itag = isrc*nprocs + idest+1 + tag_offset
                                 ! receive to pe=j
         if(commatrix_recv(i).gt.0) then
            ij = ij+1
            istrt = ir_buf(1,i)
            ilast = ir_buf(2,i)
            isize = ilast-istrt+1

#ifdef AMR_ERROR_CHECKING
            if(isize.gt.buf_dim) then
               write(*,*) 'PARAMESH ERROR : mpi_xchange_blocks 1:',
     .                    ' message is bigger than buffer space  '
     .           ,' isrc ',isrc,' idest ',idest,' i ',i,
     .           ' istrt ',istrt,' ilast ',ilast,' isize ',isize,
     .           ' buf_dim ',buf_dim
                 call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
             endif
#endif
             call Mpi_Irecv(R_buffer(istrt),isize,
     .                     MPI_REAL,
     .                     isrc ,itag,MPI_COMM_WORLD,
     .                     recvrequest(ij),ierr)
        endif
      enddo

      do j = mype+1,nprocs
         isrc = mype
         idest= j-1
         itag = isrc*nprocs + idest+1 + tag_offset
				! send from mype=i
         if(commatrix_send(j).gt.0) then
            ji = ji+1
            istrt = is_buf(1,j)
            ilast = is_buf(2,j)
            isize = ilast-istrt+1

#ifdef AMR_ERROR_CHECKING
             if(isize.gt.buf_dim) then
                 write(*,*) 'PARAMESH ERROR : ',
     .                   ' message is bigger than buffer space  '
                 call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
             endif
#endif
             call MPI_Ssend(S_buffer(istrt),isize,
     .                      MPI_REAL,
     .                      idest,itag,MPI_COMM_WORLD,ierr)
        endif
      enddo

      do j = 1,mype
         isrc = mype
         idest= j-1
         itag = isrc*nprocs + idest+1 + tag_offset
				! send from mype=i
         if(commatrix_send(j).gt.0) then
            ji = ji+1
            istrt = is_buf(1,j)
            ilast = is_buf(2,j)
            isize = ilast-istrt+1

#ifdef AMR_ERROR_CHECKING
             if(isize.gt.buf_dim) then
               write(*,*) 'PARAMESH ERROR : mpi_xchange_blocks 2:',
     .                   ' message is bigger than buffer space  '
                 call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
             endif
#endif
             call MPI_Ssend(S_buffer(istrt),isize,
     .                      MPI_REAL,
     .                      idest,itag,MPI_COMM_WORLD,ierr)
        endif
      enddo

      if(ij.gt.0)
     .   call MPI_Waitall(ij,recvrequest,recvstatus,
     .                    ierrorcode)

! reset offset with largest tag
      tag_offset = (nprocs-1)*nprocs + nprocs + tag_offset          

      return
      end subroutine mpi_xchange_blocks



      subroutine mpi_xchange_tree_info(mype,nprocs, tag_offset,
     .                              buf_dim, S_buffer, R_buffer)

!------------------------------------------------------------------------
!
! This routine uses the general hand-shaking scheme for message passing 
! between pe's. 
! sends and receives buffers containing all requested block data.
!
!
! Written :     Maharaj Bhat & Michael Gehmeyr          March 2000
!------------------------------------------------------------------------
!
! Arguments:
!      mype           rank of local processor
!      nprocs         number of processors
!      buf_dim        dimension of buffers
!      tag_offset     offset for MPI tag
!      S_buffer       send buffer
!      R_buffer       recv buffer
!
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace
      use mpi_morton

      implicit none

      include 'mpif.h'

      integer, intent(in)    ::  mype,nprocs,buf_dim
      integer, intent(inout) ::  tag_offset
      real,    intent(inout) :: S_buffer(buf_dim), R_buffer(buf_dim)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local variables

      integer :: i,j, ij, ji
      integer :: istrt,ilast,isize
      integer :: isrc,idest, itag
      integer :: ierr,errcode
      integer,dimension (:),  allocatable :: recvrequest
      integer,dimension (:,:),allocatable :: recvstatus


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      if(allocated(recvrequest)) deallocate( recvrequest )
      allocate ( recvrequest(nprocs) )
      if(allocated(recvstatus)) deallocate( recvstatus )
      allocate ( recvstatus(MPI_STATUS_SIZE,nprocs) )

!!      R_buffer = 0

      ij = 0
      ji = 0

      do i = 1,nprocs
         isrc = i-1
         idest= mype
         itag = isrc*nprocs + idest+1 + tag_offset
                                 ! receive to pe=j
         if(commatrix_recv(i).gt.0) then
            ij = ij+1
            istrt = ir_buf(1,i)
            ilast = ir_buf(2,i)
            isize = ilast-istrt+1

#ifdef AMR_ERROR_CHECKING
             if(isize.gt.buf_dim) then
               write(*,*) 'PARAMESH ERROR : mpi_xchange_blocks 3:',
     .                   ' message is bigger than buffer space  '
                 call mpi_abort(MPI_COMM_WORLD,errcode,ierr)
             endif
#endif

             call Mpi_Irecv(R_buffer(istrt),isize,MPI_REAL,
     .           isrc ,itag,MPI_COMM_WORLD,recvrequest(ij),ierr)
        endif
      enddo

      do j = 1,nprocs
         isrc = mype
         idest= j-1
         itag = isrc*nprocs + idest+1 + tag_offset
				! send from mype=i
         if(commatrix_send(j).gt.0) then
            ji = ji+1
            istrt = is_buf(1,j)
            ilast = is_buf(2,j)
            isize = ilast-istrt+1

#ifdef AMR_ERROR_CHECKING
             if(isize.gt.buf_dim) then
               write(*,*) 'PARAMESH ERROR : mpi_xchange_blocks 4:',
     .                   ' message is bigger than buffer space  '
                 call mpi_abort(MPI_COMM_WORLD,errcode,ierr)
             endif
#endif

             call MPI_Ssend(S_buffer(istrt),isize,MPI_REAL,
     .           idest,itag,MPI_COMM_WORLD,ierr)
        endif
      enddo

      if(ij.gt.0)
     .   call MPI_Waitall(ij,recvrequest,recvstatus,
     .                    ierr)

! reset offset with largest tag
      tag_offset = (nprocs-1)*nprocs + nprocs + tag_offset          


      deallocate( recvrequest )
      deallocate( recvstatus )

      return
      end subroutine mpi_xchange_tree_info

!--------------------------------------------------------------------

      logical function is_remote(mort,lrefine,mype)

      use mpi_morton

      implicit none

      integer, intent(in) :: mort(6), lrefine, mype
      logical :: morton_equal
      logical :: morton_greater_than
      logical :: morton_less_than

      is_remote = .false.
      if(morton_less_than   (mort(:),morton_limits(:,1,1,mype+1)) .or.
     .   morton_greater_than(mort(:),morton_limits(:,1,2,mype+1)))
     .   is_remote = .true.

      if( morton_equal(mort(:),morton_limits(:,1,1,mype+1))
     .         .and.
     .    (lrefine < morton_limits(6,2,1,mype+1)) )
     .     is_remote = .true.

      if( morton_equal(mort(:),morton_limits(:,1,2,mype+1)) 
     .        .and.
     .    (lrefine > morton_limits(6,2,2,mype+1)) )
     .     is_remote = .true.

      end function is_remote


!------------------------------------------------------------------


      logical function is_found(mort,lrefine,pe)

      use mpi_morton

      implicit none

      integer, intent(in) :: mort(6), lrefine, pe
      logical :: morton_equal
      logical :: morton_greater_than
      logical :: morton_less_than

      is_found = .false.
      if (morton_greater_than(mort(:),morton_limits(:,1,1,pe)) .or.
     .    morton_less_than   (mort(:),morton_limits(:,1,2,pe)))
     .   is_found = .true.

      if( morton_equal(mort(:),morton_limits(:,1,1,pe))
     .         .and.
     .    (lrefine >= morton_limits(6,2,1,pe)) )
     .     is_found = .true.

      if( morton_equal(mort(:),morton_limits(:,1,2,pe))
     .        .and.
     .    (lrefine <= morton_limits(6,2,2,pe)) )
     .     is_found = .true.

      end function is_found

!-------------------------------------------------------------------


      logical function morton_greater_than(mort1,mort2)

! Logical function to test if mort1 is > mort2

      implicit none

      integer, intent(in) :: mort1(6),mort2(6)
      integer :: jj

      morton_greater_than = .false.

      do jj = 1,6
         if (mort1(jj)  >   mort2(jj))
     .       morton_greater_than = .true.
         if (mort1(jj) .ne. mort2(jj)) exit
      end do

      end function morton_greater_than

!--------------------------------------------------------------------

      logical function morton_less_than(mort1,mort2)

! Logical function to test if mort1 is < mort2

      implicit none

      integer, intent(in) ::  mort1(6),mort2(6)
      integer :: jj

      morton_less_than = .false.

      do jj = 1,6
         if (mort1(jj)  <   mort2(jj))
     .       morton_less_than = .true.
         if (mort1(jj) .ne. mort2(jj)) exit
      end do

      end function morton_less_than

!--------------------------------------------------------------------

      logical function morton_equal(mort1,mort2)

! Logical function to test if mort1 is = mort2

      implicit none

      integer, intent(in) :: mort1(6),mort2(6)

      morton_equal = .false.
      if(all(mort1(:) == mort2(:))) morton_equal = .true.

      end function morton_equal


