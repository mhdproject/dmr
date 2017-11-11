!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"
!#define DEBUG_AMR
      subroutine amr_morton_order (lnblocks_old,nprocs,mype,
     .                             l_move_solution)


c By K. Olson (NASA/GSFC and GMU) 11/96


      use paramesh_dimensions
      use physicaldata
      use tree

      use paramesh_interfaces, only : amr_compute_morton, 
     &                                amr_migrate_tree_data,
     &                                amr_redist_blk,
     &                                amr_sort_by_work,
     &                                amr_sort_morton 

      implicit none

      include 'mpif.h'

      integer, intent(in)    :: lnblocks_old,nprocs,mype
      logical, intent(in)    :: l_move_solution

! local variables and arrays ----------------------------------------------

      integer :: new_loc(2,maxblocks_tr)
      integer :: tot_blocks
      integer :: ierr
      integer :: mort_no(6,2*maxblocks_tr)
      integer :: ireduce_datain(1),ireduce_dataout(1)

      logical,save :: first = .TRUE.

! -------------------------------------------------------------------------

c compute morton numbers for each cell

      call amr_compute_morton (mort_no)
 
c Sort these morton numbers into order. The subroutine amr_sort_morton
c returns the array new_loc which gives the new locations that 
c the cells are to move to (local address is the first arguement
c and processor number is the second).

      new_loc(:,:) = -1

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,file='amr_log',status='unknown',
     &        position='append')
         write (30,*) ' calling AMR_SORT_MORTON '
         print *, ' calling AMR_SORT_MORTON '
         close(30)
      end if
#endif

      call amr_sort_morton (mort_no,new_loc,nprocs)

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,file='amr_log',status='unknown',
     &        position='append')
         close(30)
      end if
#endif
      first = .FALSE.
         
c The following call to sort_by_work attempts to realign the 
c sorted list returned by sort_morton such that the work load is 
c balanced across processors.

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,file='amr_log',status='unknown',
     &        position='append')
         write (30,*) ' calling SORT_BY_WORK '
         print *, ' calling SORT_BY_WORK '
         close(30)
      end if
#endif

      ireduce_datain(1) = lnblocks
      call mpi_int_allreduce (ireduce_datain(1),ireduce_dataout(1),
     &                    1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,ierr)
      tot_blocks = ireduce_dataout(1)

      if (tot_blocks.gt.2*nprocs) then
         call amr_sort_by_work (new_loc,nprocs,mype)
      end if

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,file='amr_log',status='unknown',
     &        position='append')
         close(30)
      end if

      if (mype.eq.0) then
         open (unit=30,file='amr_log',status='unknown',
     &        position='append')
         write (30,*) ' calling MIGRATE_TREE '
         print *, ' calling MIGRATE_TREE '
         close(30)
      end if
#endif

      call amr_migrate_tree_data (new_loc,nprocs,mype)

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,file='amr_log',status='unknown',
     &        position='append')
         close(30)
      end if
#endif

c 2) move blocks of data to new locations

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,file='amr_log',status='unknown',
     &        position='append')
         write (30,*) ' CALLING redist_blk '
         print *,' CALLING redist_blk '
         close(30)
      end if
#endif

      if(l_move_solution)
     .  call amr_redist_blk(new_loc,nprocs,mype,lnblocks_old)

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,file='amr_log',status='unknown',
     &        position='append')
         close(30)
      end if
#endif

      lnblocks = new_lnblocks

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,file='amr_log',status='unknown',
     &        position='append')
         write (30,*) ' DONE MORTON ORDERING'
         write (30,*) ' '
         print *, ' DONE MORTON ORDERING '
         print *,' '
         close(30)
      end if
#endif

      return
      end subroutine amr_morton_order

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine amr_compute_morton (mort_no)


! This subroutine computes the morton numbers of each cell by interleaving 
! bits in x, y, then z order

! Returns -> mort_no

c By K. Olson (NASA/GSFC and GMU) 12/96

      use paramesh_dimensions
      use physicaldata
      use tree
      use paramesh_mpi_interfaces, only : morton_number


      implicit none

      include 'mpif.h'

      integer, intent(out) ::  mort_no(:,:)

! local variables and arrays ----------------------------------------------

      real :: xmin,ymin,zmin,xmin_loc,ymin_loc,zmin_loc
      real :: xyz_loc_vector(3), xyz_min_loc_vector(3)
      real :: x0,y0,z0

      integer :: i
      integer :: ierr
      integer :: mype
      integer, external :: shmem_my_pe

#ifdef DEBUG_MORTON
      integer inxt_bits(3),inxt_bit
#endif

! -------------------------------------------------------------------------

      mype = shmem_my_pe()

c find local minimum + maximum values of x, y, and z

      xmin_loc = 1.e10
      ymin_loc = 1.e10
      zmin_loc = 1.e10
      do i = 1,lnblocks

         if (nodetype(i).eq.1) then

            xmin_loc = min(coord(1,i)-(bsize(1,i)/2.),xmin_loc)
            if (ndim.ge.2) then
               ymin_loc = min(coord(2,i)-(bsize(2,i)/2.),ymin_loc)
            end if
            if (ndim.eq.3) then
               zmin_loc = min(coord(3,i)-(bsize(3,i)/2.),zmin_loc)
            end if

         end if

      end do

c find global min^s across processors

! Changed by M. Zingale and J. Dursi ---------------------------------------

c pack the 3 allreduces into a single allreduce with a vector of the
c minimum in each coordinate
      xyz_loc_vector(1) = xmin_loc
      xyz_loc_vector(2) = ymin_loc
      xyz_loc_vector(3) = zmin_loc

c reduce in all ndim dimensions
      call mpi_real_allreduce(xyz_loc_vector, xyz_min_loc_vector, ndim,
     +     MPI_REAL,
     +     MPI_MIN, MPI_COMM_WORLD, ierr)

c unpack the minimums
      xmin = xyz_min_loc_vector(1)
      if (ndim .ge. 2) ymin = xyz_min_loc_vector(2)
      if (ndim .ge. 3) zmin = xyz_min_loc_vector(3)

!---------------------------------------------------------------------------

      do i = 1,lnblocks
         
         x0 = coord(1,i)-xmin
         y0 = coord(2,i)-ymin
         z0 = coord(3,i)-zmin

         call morton_number(x0,y0,z0,bsize(:,i),ndim,
     .                      lrefine_max,lrefine(i),mort_no(:,i))
 
      end do

      return
      end subroutine amr_compute_morton

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine amr_sort_morton (mort_no,new_loc,nprocs)

c Subroutine to sort morton numbers

c Input -> vector of morton numbers, nprocs (no. of processors)

c Output -> new locations that cells are to migrate to
c           new_loc(1,i) is local id to move cell i to
c           new_loc(2,i) is processor id to move cell i to

c Sorting is done without regard to work here.  The new_loc^s returned
c are computed assuming equal (or nearly so) numbers of cells per
c processor.

c By K. Olson (NASA/GSFC and GMU) 11/96

      use paramesh_dimensions
      use physicaldata
      use tree

      use paramesh_interfaces, only : morton_sort

      implicit none

      include 'mpif.h'

      integer, intent(inout) :: mort_no(:,:)
      integer, intent(inout) :: new_loc(:,:)
      integer, intent(in)    :: nprocs

! local variables and arrays ----------------------------------------------

      integer :: itemp
      integer :: lnblocks2,tot_blocks,no_per_proc,idi,idp
      integer :: i,j
      integer :: excess,nprocs_y,nprocs_x,irnkg_s
      integer :: ierr
      integer :: lnblocks_left
      integer :: ix(2*maxblocks_tr)
      integer :: ix2(2*maxblocks_tr)
      integer :: irnkg(2*maxblocks_tr)
      integer :: lreflevel(2*maxblocks_tr)
      integer :: lreflevel_old(2*maxblocks_tr)
      integer :: ireduce_datain(1),ireduce_dataout(1)
      integer :: idatain(1),idataout(1)
      integer :: max_ref_lev 

      logical :: lswap

! -------------------------------------------------------------------------

! sort morton number array mort_no and the local list of refinement levels
! and also return the index associated with this permutation
      max_ref_lev = -1
      do i = 1,lnblocks
         ix(i) = i
         ix2(i) = i
         lreflevel(i) = lrefine(i)
         max_ref_lev = max(lrefine(i),max_ref_lev)
      end do

! sort top level
      if (lnblocks > 0) then

      call morton_sort(mort_no(:,1:lnblocks),ix(1:lnblocks),lnblocks)

      lreflevel_old(:) = lreflevel(:)
      do i = 1,lnblocks
         lreflevel(i) = lreflevel_old(ix(i))
      end do

      end if

! order segments with same morton number in order of increasing
! refinement level
      lswap = .true.
      do while (lswap)
        lswap = .false.
        do i = 1,lnblocks-1
          if(mort_no(1,i).eq.mort_no(1,i+1).and.
     .       mort_no(2,i).eq.mort_no(2,i+1).and.
     .       mort_no(3,i).eq.mort_no(3,i+1).and.
     .       mort_no(4,i).eq.mort_no(4,i+1).and.
     .       mort_no(5,i).eq.mort_no(5,i+1).and.
     .       mort_no(6,i).eq.mort_no(6,i+1).and.
     .       lreflevel(i).gt.lreflevel(i+1) ) then
            lswap = .true.
            itemp = ix(i)
            ix(i) = ix(i+1)
            ix(i+1) = itemp
            itemp = lreflevel(i)
            lreflevel(i) = lreflevel(i+1)
            lreflevel(i+1) = itemp
          endif
        enddo
      enddo                           ! end do while

      do i = 1,lnblocks
         irnkg(ix(i)) = i
      end do

      lnblocks_left = 0
      idatain(1) = lnblocks
      call mpi_int_scan (idatain(1),idataout(1),1,
     &               MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,
     &               ierr)
      lnblocks_left = idataout(1)
      lnblocks_left = lnblocks_left - lnblocks
      do j = 1,lnblocks
        irnkg(j) = irnkg(j) + lnblocks_left
      end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


c 1) Compute total list length.

! I copy lnblocks to lnblocks2 since lnblocks2 can be put in a save statement.
      lnblocks2 = lnblocks 
      ireduce_datain(1) = lnblocks2
      call mpi_int_allreduce (ireduce_datain,ireduce_dataout,
     &                    1,MPI_INTEGER,
     &                    MPI_SUM,MPI_COMM_WORLD,ierr)
      tot_blocks = ireduce_dataout(1)

      no_per_proc = tot_blocks/nprocs

      excess = tot_blocks - no_per_proc*nprocs
      nprocs_y = (no_per_proc+1)*nprocs - tot_blocks
! no. of processors which will get no_per_proc + 1 blocks
      nprocs_x = nprocs - nprocs_y
! rank in list which divides those which go on processor with one number
! of blocks from those which go on another set of blocks w. a different
! no. of blocks
      irnkg_s = nprocs_x*(no_per_proc+1)

c 2) Compute new_locs from rankings (irnkg) returned by amr_bi_sort.
c    The following divides blocks evenly among processors without regard to
c    work.

      do i = 1,lnblocks

         idp = (irnkg(i)-1)/(no_per_proc+1) ! processor to send to
         if (irnkg(i).le.irnkg_s) then
            idi = mod((irnkg(i)-1),no_per_proc+1) + 1 ! rank inside 
                                                      ! local array
                                                      ! to write to
         else
            idp = (irnkg(i)-irnkg_s-1)/(no_per_proc) ! processor to send to
            idp = idp + nprocs_x
            idi = mod((irnkg(i)-irnkg_s-1),no_per_proc) + 1 ! rank inside 
                                                            ! local array
                                                            ! to write to
         end if

         new_loc(1,i) = idi
         new_loc(2,i) = idp
         
      end do

      return
      end subroutine amr_sort_morton

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine amr_sort_by_work (new_loc,nprocs,mype)


c Subroutine to balance work load

c on input takes list sorted by morton ordering
c on output returns new values of new_loc

c By K. Olson (NASA/GSFC and GMU) 11/96

      use paramesh_dimensions
      use physicaldata
      use tree
      use paramesh_interfaces, only : fill_old_loc

      implicit none

      include 'mpif.h'

      integer, intent(inout) :: new_loc(:,:)
      integer, intent(in)    :: nprocs,mype

! local variables and arrays ------------------------------------------------

      real :: work_per_proc,work_left,loc_work,tot_work
      real :: work(maxblocks_tr),workt(maxblocks_tr)
      real :: reduce_datain(1),reduce_dataout(1)
      real :: rdatain(1),rdataout(1)

      integer :: ierr,errorcode
      integer :: lnblocks2
      integer :: i,j
      integer :: pidt,lidt,lid_old
      integer :: left,right
      integer :: nsend,nrecv
      integer :: reqr(2*maxblocks_tr)
      integer :: pid(maxblocks_tr),lid(maxblocks_tr),lid2(maxblocks_tr)
      integer :: stat (MPI_STATUS_SIZE)
      integer :: new_loc_temp(2,maxblocks_tr)
      integer :: old_loc(2,maxblocks_tr)
      integer :: statr(MPI_STATUS_SIZE,2*maxblocks_tr)

      logical :: repeat,repeatt
      logical :: lreduce_datain(1),lreduce_dataout(1)

! ---------------------------------------------------------------------------

c initialize work arrary and temp work arrary

      do i = 1,maxblocks_tr

         work(i) = -1.
         workt(i) = -1.

      end do

c assign values to work array

      do i = 1,maxblocks_tr

         work(i) = 0.
         work(i) = work_block(i)
!         if (nodetype(i).eq.1) work(i) = 2.
!         if (nodetype(i).eq.2) work(i) = 1.

      end do

c move work to temp work array workt

      call fill_old_loc(new_loc,old_loc,nprocs,mype)

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,
     &        file='amr_log',
     &        position='append',
     &        status='unknown',
     &        form='formatted')
         write (30,*) ' DONE fill_old_loc '
         close(30)
      end if
#endif

      lnblocks2 = 0
      do i = 1,maxblocks_tr
        if (old_loc(1,i).gt.-1) then
          lnblocks2 = lnblocks2 + 1
        end if
      end do

      nrecv = 0
      do i = 1,lnblocks2
        if (old_loc(2,i).ne.mype) then
          nrecv = nrecv + 1
          call MPI_real_IRECV(workt(i),1,MPI_REAL,
     &         old_loc(2,i),old_loc(1,i),MPI_COMM_WORLD,
     &         reqr(nrecv),ierr)
        end if
      end do

      nsend = 0
      do i = 1,lnblocks

        if (new_loc(2,i).ne.mype) then
          nsend = nsend + 1
          call MPI_real_SSEND(work(i),1,MPI_REAL,
     &         new_loc(2,i),i,MPI_COMM_WORLD,ierr)
        else
          workt(new_loc(1,i)) = work(i)
        end if

      end do
         
      if (nrecv.gt.0) then
        call MPI_WAITALL(nrecv,reqr,statr,ierr)
      end if



      do i = 1,lnblocks2
         work(i) = workt(i)
      end do

c SUM total work within each processosr

      if (lnblocks2.gt.0) then
        workt(1) = work(1)
      else
        workt(1) = 0
      end if
      do i = 2,lnblocks2

         workt(i) = workt(i-1) + work(i)

      end do

c SUM work across processors

      if (lnblocks2.gt.0) then
        loc_work = workt(lnblocks2)
      else
        loc_work = 0
      end if

      reduce_datain(1) = loc_work
      call mpi_real_allreduce (reduce_datain(1),reduce_dataout(1),
     &     1,MPI_REAL,
     &     MPI_SUM,MPI_COMM_WORLD,ierr)
      tot_work = reduce_dataout(1)

c Compute work per processor

      work_per_proc = tot_work/nprocs

c Compute final work by looking left

      work_left = 0.

      rdatain(1) = loc_work
!      call mpi_real_scan (loc_work,work_left,1,
      call mpi_real_scan (rdatain(1),rdataout(1),1,
     &               MPI_REAL,MPI_SUM,MPI_COMM_WORLD,
     &               ierr)
      work_left = rdataout(1)

      work_left = work_left - loc_work
      do j = 1,lnblocks2
         workt(j) = workt(j) + work_left
      end do

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,
     &        file='amr_log',
     &        position='append',
     &        status='unknown',
     &        form='formatted')
!         write (30,*) ' DONE compute work_left ',mype
         close(30)
      end if
#endif
c compute processor ids

      do i = 1,maxblocks_tr

         pid(i) = 0
         lid(i) = 0

      end do

      do i = 1,lnblocks2

         pid(i) = int((workt(i)-1.)/work_per_proc)
         if (pid(i).lt.0) pid(i) = 0
         if (pid(i).gt.nprocs-1) pid(i) = nprocs-1

      end do

c compute local ids
      
      lid(1) = 1
      do i = 2,lnblocks2

         lid(i) = lid(i-1) + 1
         if (pid(i-1).lt.pid(i)) lid(i) = 1  ! start a new group

      end do

      do i = 1,maxblocks_tr
         lid2(i) = lid(i)
      end do

      left = mype - 1
      right = mype + 1
      if (mype.eq.0) left = MPI_PROC_NULL
      if (mype.eq.nprocs-1) right = MPI_PROC_NULL

      call MPI_SENDRECV 
     &     (pid(lnblocks2),1,MPI_INTEGER,right,1,
     &      pidt,          1,MPI_INTEGER,left, 1,
     &      MPI_COMM_WORLD,stat,ierr)

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,
     &        file='amr_log',
     &        position='append',
     &        status='unknown',
     &        form='formatted')
         write (30,*) ' STARTING loop 27 '
         close(30)
      end if
#endif

      lidt = 0
 27   lid_old = lidt ! lid_old stores last fetched value of lid to left

      lidt = 0

      call MPI_SENDRECV 
     &     (lid(lnblocks2),1,MPI_INTEGER,right,1,
     &      lidt,          1,MPI_INTEGER,left, 1,
     &      MPI_COMM_WORLD,stat,ierr)

      do j = 1,lnblocks2
         
         if (pidt.eq.pid(j)) then ! if pidt (which was fetched)
                                  ! equals local pid then the list
                                  ! has been split across processors
               
            lid(j) = lid2(j) + lidt
            
         end if

      end do
      
      repeat = .FALSE.
      if (lidt.ne.lid_old) repeat = .TRUE.
      lreduce_datain(1) = repeat
      call mpi_logical_allreduce (lreduce_datain(1),lreduce_dataout(1),
     &                    1,MPI_LOGICAL,MPI_LOR,MPI_COMM_WORLD,ierr)
      repeatt = lreduce_dataout(1)
      if (repeatt) go to 27

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,
     &        file='amr_log',
     &        position='append',
     &        status='unknown',
     &        form='formatted')
         write (30,*) ' DONE loop 27 '
         close(30)
      end if
#endif

c now reorder according to new pid and lid numbers

      nrecv = 0
      do i = 1,lnblocks
         if (new_loc(2,i).ne.mype) then
            nrecv = nrecv + 1
            call MPI_int_IRECV(new_loc_temp(1,i),1,MPI_INTEGER,
     &           new_loc(2,i),new_loc(1,i),MPI_COMM_WORLD,
     &           reqr(nrecv),ierr)
         else
            new_loc_temp(1,i) = lid(new_loc(1,i))
         end if
      end do

      nsend = 0
      do i = 1,lnblocks2
        if (old_loc(2,i).ne.mype) then
           nsend = nsend + 1
           call MPI_int_SSEND(lid(i),1,MPI_INTEGER,
     &          old_loc(2,i),i,
     &          MPI_COMM_WORLD,ierr)
        end if
      end do

      if (nrecv.gt.0) then
        call MPI_WAITALL(nrecv,reqr,statr,ierr)
      end if

      nrecv = 0
      do i = 1,lnblocks
         if (new_loc(2,i).ne.mype) then
            nrecv = nrecv + 1
            call MPI_int_IRECV(new_loc_temp(2,i),1,MPI_INTEGER,
     &           new_loc(2,i),new_loc(1,i),
     &           MPI_COMM_WORLD,
     &           reqr(nrecv),ierr)
         else
            new_loc_temp(2,i) = pid(new_loc(1,i))
         end if
      end do
      
      nsend = 0
      do i = 1,lnblocks2
         if (old_loc(2,i).ne.mype) then
           nsend = nsend + 1
           call MPI_int_SSEND(pid(i),1,MPI_INTEGER,
     &          old_loc(2,i),i,
     &          MPI_COMM_WORLD,ierr)
        end if
      end do
         
      if (nrecv.gt.0) then
        call MPI_WAITALL(nrecv,reqr,statr,ierr)
      end if

      do i = 1,lnblocks

         new_loc(:,i) = new_loc_temp(:,i)

         if(new_loc(1,i).gt.maxblocks) then
            open (unit=30,
     &           file='amr_log',
     &           position='append',
     &           status='unknown',
     &           form='formatted')
            write(30,*) 'PARAMESH ERROR !'
            write(30,*) mype,i
            write(30,*) ' new_loc(1 = ',new_loc(1,i),new_loc(2,i)
            write(30,*) 'New block location exceeds MAXBLOCKS limit'
            write(30,*) 'Suggestion: increase MAXBLOCKS or modify',
     .           ' refinement criteria'
            close(30)
            call MPI_ABORT (MPI_COMM_WORLD,errorcode,ierr)
         endif         
         
         if(new_loc(2,i).gt.nprocs.or.new_loc(2,i).lt.0) then
            open (unit=30,
     &           file='amr_log',
     &           position='append',
     &           status='unknown',
     &           form='formatted')
            write(30,*) 'PARAMESH ERROR !'
            write(30,*) mype,i
            write(30,*) 'New block location out of bounds',new_loc(2,i)
            write(30,*) 'Suggestion: increase MAXBLOCKS or modify',
     .           ' refinement criteria'
            close(30)
            call MPI_ABORT (MPI_COMM_WORLD,errorcode,ierr)
         endif         

      end do
      new_loc(:,lnblocks+1:maxblocks_tr) = -1

#ifdef DEBUG_AMR
      if (mype.eq.0) then
         open (unit=30,
     &        file='amr_log',
     &        position='append',
     &        status='unknown',
     &        form='formatted')
         write (30,*) ' DONE sort_by_work '
         close(30)
      end if
#endif

      return
      end subroutine amr_sort_by_work

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine amr_migrate_tree_data (new_loc,nprocs,mype)


c Subroutine move tree data and reconnect all pointers given new_loc

c By K. Olson (NASA/GSFC and GMU) 11/96


      use paramesh_dimensions
      use physicaldata
      use tree
      use paramesh_interfaces, only : fill_old_loc

      implicit none

      include 'mpif.h'

      integer, intent(inout) :: new_loc(:,:)
      integer, intent(in)    :: nprocs,mype

! local variables and arrays ------------------------------------------------

      integer, parameter :: buf_size = mdim+mdim+2*mdim
      integer, parameter :: ibuf_size = 2*mfaces+2*mchild+2+5+mflags

      real :: buffer(buf_size)
      real :: buffert(buf_size,maxblocks_tr)

      integer :: ibuffer(ibuf_size)
      integer :: ibuffert(ibuf_size,maxblocks_tr)
      integer :: neight(2,mfaces,maxblocks_tr)
      integer :: childt(2,mchild,maxblocks_tr)
      integer :: parentt(2,maxblocks_tr)
      integer :: i,j,k,jj
      integer :: old_loc(2,maxblocks_tr)
      integer :: statr(MPI_STATUS_SIZE,maxblocks_tr)
      integer :: reqr(maxblocks_tr)
      integer :: ierr,nsend,nrecv

      logical newchildt(maxblocks_tr)


#ifdef SAVE_MORTS
      integer :: isize_sm

!???? Need to rework so we don't need such a large temporary array for surr_mortst ????
      integer :: surr_mortst(6,3,1+2*k2d,1+2*k3d,maxblocks_tr)
#endif
! ----------------------------------------------------------------------------

      call fill_old_loc(new_loc,old_loc,nprocs,mype)


c count no. of new blocks

      new_lnblocks = 0
      do i = 1,maxblocks_tr
        if (old_loc(1,i).gt.-1) then
          new_lnblocks = new_lnblocks + 1
        end if
      end do

c update pointers to parents, children and neighbors

      parentt(:,1:lnblocks) = parent(:,1:lnblocks)
      childt(:,:,1:lnblocks) = child(:,:,1:lnblocks)
      neight(:,:,1:lnblocks) = neigh(:,:,1:lnblocks)

      nrecv = 0
      do i = 1,lnblocks
         if (parent(1,i).gt.0) then
           if (parent(2,i).ne.mype) then
             nrecv = nrecv + 1
             call MPI_int_IRECV(parentt(1,i),2,MPI_INTEGER,
     &            parent(2,i),i,MPI_COMM_WORLD,
     &            reqr(nrecv),ierr)
           else
             parentt(:,i) = new_loc(:,parent(1,i))
           end if
         end if
       end do
       
       nsend = 0
       do i = 1,lnblocks
         do j = 1,nchild
           if (child(1,j,i).gt.0) then
             if (child(2,j,i).ne.mype) then
               ! parent is sending to all its children
               nsend = nsend + 1
               call MPI_int_SSEND (new_loc(1,i),2,MPI_INTEGER,
     &              child(2,j,i),child(1,j,i),MPI_COMM_WORLD,
     &              ierr)
             end if
           end if
         end do
       end do

      if (nrecv.gt.0) then
        call MPI_WAITALL(nrecv,reqr,statr,ierr)
      end if

      nrecv = 0
      do i = 1,lnblocks
        do j = 1,nchild
          if (child(1,j,i).gt.0) then
            if (child(2,j,i).ne.mype) then
              nrecv = nrecv + 1
              call MPI_int_IRECV(childt(1,j,i),2,MPI_INTEGER,
     &             child(2,j,i),child(1,j,i),MPI_COMM_WORLD,
     &             reqr(nrecv),ierr)
            else
              childt(:,j,i) = new_loc(:,child(1,j,i))
            end if
          end if
        end do
       end do
       
       nsend = 0
       do i = 1,lnblocks
         if (parent(1,i).gt.0) then
           if (parent(2,i).ne.mype) then
! child is sending to its parent
             nsend = nsend + 1
             call MPI_int_SSEND (new_loc(1,i),2,MPI_INTEGER,
     &            parent(2,i),i,MPI_COMM_WORLD,
     &            ierr)
           end if
         end if
       end do

      if (nrecv.gt.0) then
        call MPI_WAITALL(nrecv,reqr,statr,ierr)
      end if

      nrecv = 0
      do i = 1,lnblocks
        do j = 1,nfaces
          if (neigh(1,j,i).gt.0) then
            if (neigh(2,j,i).ne.mype) then
              nrecv = nrecv + 1
              call MPI_int_IRECV(neight(1,j,i),2,MPI_INTEGER,
     &             neigh(2,j,i),neigh(1,j,i),MPI_COMM_WORLD,
     &             reqr(nrecv),ierr)
            else
              neight(:,j,i) = new_loc(:,neigh(1,j,i))
            end if
          end if
        end do
      end do
      
      nsend = 0
      do i = 1,lnblocks
        do j = 1,nfaces
          if (neigh(1,j,i).gt.0) then
            if (neigh(2,j,i).ne.mype) then
              nsend = nsend + 1
              call MPI_int_SSEND (new_loc(1,i),2,MPI_INTEGER,
     &             neigh(2,j,i),i,MPI_COMM_WORLD,
     &             ierr)
            end if
          end if
        end do
      end do

      if (nrecv.gt.0) then
        call MPI_WAITALL(nrecv,reqr,statr,ierr)
      end if

#ifdef SAVE_MORTS
! exchange surr_morts data
      write(*,*) 'pe ',mype,' lnblocks ',lnblocks
      write(*,*) 'pe ',mype,' new_lnblocks ',new_lnblocks
      write(*,*) 'pe ',mype,' old_loc ',old_loc(:,1:lnblocks)
      write(*,*) 'pe ',mype,' new_loc ',new_loc(:,1:lnblocks)

      isize_sm = 6*3*(1+2*k2d)*(1+2*k3d)
      nrecv = 0
      do i = 1,new_lnblocks
        if (old_loc(2,i).ne.mype.and.old_loc(2,i).gt.-1) then
          nrecv = nrecv + 1
          call MPI_int_IRECV(surr_mortst(1,1,1,1,i),
     &                  isize_sm,MPI_INTEGER,
     &                  old_loc(2,i),old_loc(1,i),MPI_COMM_WORLD,
     &                  reqr(nrecv),ierr)
!        write(*,*) 'pe ',mype,' posting recv in blk ',i,' from ',
!     .        old_loc(2,i),old_loc(1,i),' tag ',old_loc(1,i)
        else
          if(i.ne.old_loc(1,i)) then
            surr_mortst(:,:,:,:,i) = surr_morts(:,:,:,:,old_loc(1,i))
!        write(*,*) 'pe ',mype,' local copy into blk ',i,' from ',
!     .        old_loc(2,i),old_loc(1,i), 
!     .      ' surr_mortst ',surr_mortst(6,:,:,1,i)
          end if
        end if
      end do

      nsend = 0
      do i = 1,lnblocks
        if (new_loc(2,i).ne.mype.and.new_loc(2,i).gt.-1) then
          nsend = nsend + 1
          call MPI_int_SSEND (surr_morts(1,1,1,1,i),
     &             isize_sm,MPI_INTEGER,
     &             new_loc(2,i),i,MPI_COMM_WORLD,
     &             ierr)
!        write(*,*) 'pe ',mype,' posting send from blk ',i,' to ',
!     .        new_loc(2,i),new_loc(1,i),' tag ',i,
!     .   ' surr_morts ',surr_morts(6,:,:,1,i)
        end if
      end do

      if (nrecv.gt.0) then
        call MPI_WAITALL(nrecv,reqr,statr,ierr)
      end if
    
 
      do i = 1,new_lnblocks
      if(i.ne.old_loc(1,i).or.mype.ne.old_loc(2,i)) then
      surr_morts(:,:,:,:,i) = 
     .                   surr_mortst(:,:,:,:,i) 
      endif
      enddo
      
!      surr_morts(:,:,:,:,1:new_lnblocks) = 
!     .                   surr_mortst(:,:,:,:,1:new_lnblocks) 
#endif /* SAVE_MORTS */

      parent(:,1:lnblocks) = parentt(:,1:lnblocks)
      child(:,:,1:lnblocks) = childt(:,:,1:lnblocks)
      neigh(:,:,1:lnblocks) = neight(:,:,1:lnblocks)

c initialize temp buffer array

      do i = 1,maxblocks_tr
         buffert(:,i) = -1.
         ibuffert(:,i) = -1
         newchildt(i) = .FALSE.
      end do

      nrecv = 0
      do i = 1,new_lnblocks
        if (old_loc(2,i).ne.mype) then
          nrecv = nrecv + 1
          call MPI_real_IRECV(buffert(1,i),buf_size,
     &         MPI_REAL,
     &         old_loc(2,i),i,MPI_COMM_WORLD,
     &         reqr(nrecv),ierr)
          nrecv = nrecv + 1
          call MPI_int_IRECV(ibuffert(1,i),ibuf_size,MPI_INTEGER,
     &         old_loc(2,i),i+maxblocks_tr,MPI_COMM_WORLD,
     &         reqr(nrecv),ierr)
          nrecv = nrecv + 1
          call MPI_logical_IRECV(newchildt(i),1,MPI_LOGICAL,
     &         old_loc(2,i),i+2*maxblocks_tr,MPI_COMM_WORLD,
     &         reqr(nrecv),ierr)
        end if
      end do

      nsend = 0

      do i = 1,lnblocks

! pack buffer for sending

        k = 0
        do j = 1,mdim
          k = k + 1
          buffer(k) = coord(j,i)
        end do
        do j = 1,mdim
          do jj = 1,2
            k = k + 1
            buffer(k) = bnd_box(jj,j,i)
          end do
        end do
        do j = 1,mdim
          k = k + 1
          buffer(k) = bsize(j,i)
        end do

        k = 0
        do j = 1,mchild
          do jj = 1,2
            k = k + 1
            ibuffer(k) = child(jj,j,i)
          end do
        end do
        do j = 1,mfaces
          do jj = 1,2
            k = k + 1
            ibuffer(k) = neigh(jj,j,i)
          end do
        end do
        do j = 1,2
          k = k + 1
          ibuffer(k) = parent(j,i)
        end do
        k = k + 1
        ibuffer(k) = lrefine(i)
        k = k + 1
        ibuffer(k) = nodetype(i)
        k = k + 1
        ibuffer(k) = which_child(i)
        k = k + 1
        ibuffer(k) = empty(i)
        k = k + 1
        ibuffer(k) = work_block(i)
        do j=1,mflags
        k = k + 1
        ibuffer(k) = bflags(j,i)
        enddo

        if (new_loc(2,i).ne.mype) then
          nsend = nsend + 1
          call MPI_real_SSEND(buffer(1),buf_size,MPI_REAL,
     &         new_loc(2,i),new_loc(1,i),
     &         MPI_COMM_WORLD,ierr)
          nsend = nsend + 1
          call MPI_int_SSEND(ibuffer(1),ibuf_size,MPI_INTEGER,
     &         new_loc(2,i),new_loc(1,i)+maxblocks_tr,
     &         MPI_COMM_WORLD,ierr)
          nsend = nsend + 1
          call MPI_logical_SSEND(newchild(i),1,MPI_LOGICAL,
     &         new_loc(2,i),new_loc(1,i)+2*maxblocks_tr,
     &         MPI_COMM_WORLD,ierr)
        else
          buffert(1:buf_size,new_loc(1,i)) = buffer(1:buf_size)
          ibuffert(1:ibuf_size,new_loc(1,i)) = ibuffer(1:ibuf_size)
          newchildt(new_loc(1,i)) = newchild(i)
        end if

      end do

      if (nrecv.gt.0) then
        call MPI_WAITALL(nrecv,reqr,statr,ierr)
      end if

 ! now unpack the buffer
      do i = 1,maxblocks_tr

        k = 0
        do j = 1,mdim
          k = k + 1
          coord(j,i) = buffert(k,i)
        end do
        do j = 1,mdim
          do jj = 1,2
            k = k + 1
            bnd_box(jj,j,i) = buffert(k,i)
          end do
        end do
        do j = 1,mdim
          k = k + 1
          bsize(j,i) = buffert(k,i)
        end do

        k = 0
        do j = 1,mchild
          do jj = 1,2
            k = k + 1
            child(jj,j,i) = ibuffert(k,i)
          end do
        end do
        do j = 1,mfaces
          do jj = 1,2
            k = k + 1
            neigh(jj,j,i) = ibuffert(k,i)
          end do
        end do
        do j = 1,2
          k = k + 1
          parent(j,i) = ibuffert(k,i)
        end do
        k = k + 1
        lrefine(i) = ibuffert(k,i)
        k = k + 1
        nodetype(i) = ibuffert(k,i)
        k = k + 1
        which_child(i) = ibuffert(k,i)
        k = k + 1
        empty(i) = max(ibuffert(k,i),0)     ! empty must be either 1 or 0.
                                            ! It cannot be -1.
        k = k + 1
        work_block(i) = ibuffert(k,i)
        do j=1,mflags
        k = k + 1
        bflags(j,i) = ibuffert(k,i)
        enddo

        newchild(i) = newchildt(i)

      end do
      
!!!      call MPI_BARRIER (MPI_COMM_WORLD,ierr)

      return
      end subroutine amr_migrate_tree_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine fill_old_loc(new_loc,old_loc,nprocs,mype)

      use paramesh_dimensions
      use physicaldata
      use tree

      implicit none

      include 'mpif.h'

      integer, intent(inout) :: new_loc(:,:)
      integer, intent(out)   :: old_loc(:,:)
      integer, intent(in)    :: nprocs,mype

! local variables and arrays ------------------------------------------------

      integer :: nrecv,nsend
      integer :: nsend_to_proc(0:16384), nrecv_pack(0:16384)
      integer :: ierr
      integer :: i
      integer :: statr(MPI_STATUS_SIZE,maxblocks_tr)
      integer :: reqr(maxblocks_tr)
      integer :: kk(maxblocks_tr)
      integer :: idummy_array(1)

! --------------------------------------------------------------------------

! fill `old_loc' (pointer from new block location back to
! its old, unsorted location)

! count no. of receives to post

! 1) count no. of sends on each proc. to all other procs

      nsend_to_proc(:) = 0
      do i = 1,maxblocks_tr
        if (new_loc(1,i).gt.0) then
        if (new_loc(2,i).ne.mype) then ! its a send
          nsend_to_proc(new_loc(2,i)) = 
     &    nsend_to_proc(new_loc(2,i)) + 1
        end if
        end if
      end do

! 2) collect data for `this' proc from other procs so that
!    the total no. of receives to post can be computed

! Changed by M. Zingale and J. Dursi --------------------------------------

      nrecv = 0
      call MPI_AllReduce(nsend_to_proc, nrecv_pack, nprocs,
     +     MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
      nrecv = nrecv_pack(mype)

!---------------------------------------------------------------------------

      old_loc(:,:) = -1
      do i = 1,nrecv
        call MPI_int_IRECV(kk(i),
     .       1,
     .       MPI_INTEGER,
     .       MPI_ANY_SOURCE,
     .       MPI_ANY_TAG,
     .       MPI_COMM_WORLD,
     .       reqr(i),
     .       ierr)
      end do

      nsend = 0
      do i = 1,maxblocks_tr
        if (new_loc(1,i).gt.0) then
        if (new_loc(2,i).ne.mype) then
          nsend = nsend + 1
! new_loc(2,i) - PE TO SEND TO
! new_loc(1,i) - THIS IS THE TAG
          idummy_array(1) = i
          call MPI_int_SSEND(idummy_array(1),
     &         1,
     &         MPI_INTEGER,
     &         new_loc(2,i),
     &         new_loc(1,i),
     &         MPI_COMM_WORLD,
     &         ierr)
        else
          old_loc(1,new_loc(1,i)) = i
          old_loc(2,new_loc(1,i)) = mype
        end if
        end if
      end do

      if (nrecv.gt.0) then
         call MPI_WAITALL (nrecv, reqr, statr, ierr)
         do i = 1,nrecv
           old_loc(1,statr(MPI_TAG,i)) = kk(i)
           old_loc(2,statr(MPI_TAG,i)) = statr(MPI_SOURCE,i)
         end do
      end if

      call MPI_BARRIER (MPI_COMM_WORLD,ierr) ! NEEDED ????

      return
      end subroutine fill_old_loc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine morton_sort(mort_no,ix,iend)

      use paramesh_interfaces, only : amr_q_sort

      implicit none

      integer, intent(in) :: iend
      integer, intent(inout) :: mort_no(6,iend), ix(iend)

! local variables and arrays ------------------------------------------------

      integer :: ix2(iend)
      integer :: i, is, ie, level
      integer :: mort_no_old(6,iend)

      logical :: greater

! ---------------------------------------------------------------------------

      if (iend > 0) then

      do i = 1,iend
         ix2(i) = i
      end do

      is = 1
      ie = iend
      mort_no_old(:,is:ie) = mort_no(:,is:ie)

      greater = .false.
      do i = is, ie-1
         if (mort_no(1,i) > mort_no(1,i+1)) then
            greater = .true.
            exit
         end if
      end do
      if (greater) then
         call amr_q_sort(mort_no(1,is:ie),
     .                   ie-is+1,
     .                   ix2(is:ie),
     .                   ix(is:ie))
         do i = is,ie
            mort_no(:,i) = mort_no_old(:,ix2(i))
         end do
      end if

      do level = 2,6
         ie = 0
         do while (ie < iend)
         is = ie + 1
         ie = is
         do
            if (any(mort_no(1:level-1,is) .ne. 
     .              mort_no(1:level-1,ie))) then
               ie = ie - 1
               exit
            end if
            ie = ie + 1
            if (ie > iend) then
               ie = ie - 1
               exit
            end if
         end do

         if (ie > is) then
         mort_no_old(:,is:ie) = mort_no(:,is:ie)

         greater = .false.
         do i = is, ie-1
            if (mort_no(level,i) > mort_no(level,i+1)) then
               greater = .true.
               exit
            end if
         end do
         if (greater) then
            do i = is,ie
               ix2(i) = i
            end do
            call amr_q_sort(mort_no(level,is:ie),
     .                      ie-is+1,
     .                      ix2(is:ie),
     .                      ix(is:ie))
            do i = is,ie
               mort_no(:,i) = mort_no_old(:,ix2(i))
            end do
         end if
         
         end if

         end do
      end do

      end if

      end subroutine morton_sort
