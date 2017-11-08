!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"


#define DEBUG_MORTON
      subroutine amr_morton_order (nprocs,lnblocks_old,
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

#include "amr_shmem.fh"

      integer, intent(in) ::  nprocs,lnblocks_old
      logical, intent(in) ::  l_move_solution

      integer shmem_my_pe,shmem_n_pes,mype

      real tot_work

      integer i
      integer iii,lb      
      integer mort_no(2*maxblocks_tr),new_loc(2,maxblocks_tr)



      mype = shmem_my_pe()

c compute morton numbers for each cell

      call amr_compute_morton (mort_no,nprocs)

c Sort these morton numbers into order. The subroutine amr_sort_morton
c returns the array new_loc which gives the new locations that 
c the cells are to move to (local address is the first arguement
c and processor number is the second).

      new_loc(:,:) = -1


      call amr_sort_morton (mort_no,new_loc,nprocs)


c The following call to sort_by_work attempts to realign the 
c sorted list returned by sort_morton such that the work load is 
c balanced across processors.
      call amr_sort_by_work (new_loc,nprocs)

      call amr_migrate_tree_data (new_loc,nprocs)



c 2) move blocks of data to new locations
c 
      if(l_move_solution)
     .      call amr_redist_blk(mype,new_loc,lnblocks_old)

      lnblocks = new_lnblocks

      call shmem_barrier_all()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c DEBUGGING: print out total work on each processor

!      tot_work = 0.
!      do i = 1,lnblocks

!         if (nodetype(i).eq.1) tot_work = tot_work + 2.
!         if (nodetype(i).eq.2) tot_work = tot_work + 1.

!      end do

!      print *,' tot_work after sort = ',tot_work,mype

      call shmem_barrier_all()
      return
      end subroutine amr_morton_order

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine amr_compute_morton (mort_no,nprocs)


! This subroutine computes the morton numbers of each cell by interleaving 
! bits in x, y, then z order

! Returns -> mort_no

c By K. Olson (NASA/GSFC and GMU) 12/96

      use paramesh_dimensions
      use physicaldata
      use tree

      use paramesh_interfaces, only : comm_real_min_to_all,
     &                                comm_int_max_to_all
      implicit none

#include "amr_shmem.fh"

      include 'shmem_reduce.fh'

!      integer, intent(out) ::  mort_no(2*maxblocks_tr)
      integer, intent(out) ::  mort_no(:)
      integer, intent(in)  ::  nprocs


      integer shmem_my_pe,mype
      
      real    coord_norm
      real    size_min,size_min_loc,temp
      real    xmin,ymin,zmin,xmin_loc,ymin_loc,zmin_loc
      integer ix(maxblocks_tr),iy(maxblocks_tr),iz(maxblocks_tr)
      integer ipos,ipos2,i,j,k
      integer max_level,max_level_loc,nbits,nbitshft

      save size_min,size_min_loc
      save xmin,ymin,zmin,xmin_loc,ymin_loc,zmin_loc
      save max_level,max_level_loc


#ifdef DEBUG_MORTON
      integer inxt_bits(3),inxt_bit
#endif

      mype = shmem_my_pe()

c 1) find minimum size of meshes

      size_min_loc = 1.e10
      do i = 1,lnblocks

         if (nodetype(i).eq.1) then

            do j = 1,ndim
               size_min_loc = min(size_min_loc,bsize(j,i))
            end do
            
         end if

      end do

      call shmem_barrier_all()

c 1.1) find global size_min across processors
      
      call comm_real_min_to_all (size_min,size_min_loc)

      call shmem_barrier_all()

c 2) find local minimum values of x, y, and z

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

      call shmem_barrier_all()

c 2.1) find global min^s across processors

      call comm_real_min_to_all (xmin,xmin_loc)

      call shmem_barrier_all()

      if (ndim.ge.2) then
         call comm_real_min_to_all (ymin,ymin_loc)
      end if

      call shmem_barrier_all()

      if (ndim.eq.3) then
         call comm_real_min_to_all (zmin,zmin_loc)
      end if

      call shmem_barrier_all()

c 3) compute ix,iy, and iz

      do i = 1,lnblocks
         
c         ix(i) = int((coord(1,i)-xmin)/size_min)
c         iy(i) = int((coord(2,i)-ymin)/size_min)
c         if (ndim.eq.3) then
c          iz(i) = int((coord(3,i)-zmin)/size_min)
c         else
c          iz(i) = 0
c         end if

c  CODE THAT FOLLOWS gives morton ordering according to a true tree 
c  structure. I.e. child morton numbrs are at least equal to or
c  greater than their parents number.

         coord_norm = (coord(1,i)-xmin)/bsize(1,i)
         ix(i) = int( coord_norm )
!        ix(i) = int( ((coord(1,i)-xmin)/bsize(1,i)) )
         if (ndim.ge.2) then
            coord_norm = (coord(2,i)-ymin)/bsize(2,i)
            iy(i) = int( coord_norm )
!           iy(i) = int( ((coord(2,i)-ymin)/bsize(2,i)) )
         else
            iy(i) = 0
         end if
         if (ndim.eq.3) then
            coord_norm = (coord(3,i)-zmin)/bsize(3,i)
            iz(i) = int( coord_norm )
!           iz(i) = int( ((coord(3,i)-zmin)/bsize(3,i)) )
         else
            iz(i) = 0
         end if
           
      end do
      
c 4) now interleave bits of ix, iy, and iz to get the morton numbers

      do i = 1,maxblocks_tr*2
         mort_no(i) = -1
      enddo
      
      do i = 1,lnblocks
        
 
         ipos = 0
         ipos2 = 0
         mort_no(i) = 0

!        do while (ipos.lt. 31-3)
! determine bit size of a typical morton number on this machine        
         nbits = bit_size(mort_no(1)) 
         do while (ipos.lt. nbits-1-3)
              
            call mvbits (ix(i),ipos2,1,mort_no(i),ipos)
            ipos = ipos + 1
            if(ndim.ge.2) then
               call mvbits (iy(i),ipos2,1,mort_no(i),ipos)
               ipos = ipos + 1
            endif
            if(ndim.eq.3) then
               call mvbits (iz(i),ipos2,1,mort_no(i),ipos)
               ipos = ipos + 1
            endif
            
            ipos2 = ipos2 + 1

         end do

#ifdef DEBUG_MORTON
            inxt_bits(1) = ibits(ix(i),ipos2,1)
            inxt_bits(2) = ibits(iy(i),ipos2,1)
            inxt_bits(3) = ibits(iz(i),ipos2,1)
            inxt_bit = maxval(inxt_bits)
            if(inxt_bit.gt.0) then
                write(*,*) 'Error 1 in amr_morton_compute',
     .                     'mype,ipos2,inxt_bits ix iy iz : ',
     .                      mype,ipos2,inxt_bits,ix(i),iy(i),iz(i)
                call amr_abort()
            endif
#endif
              

      end do

      call shmem_barrier_all()


c  CODE THAT FOLLOWS gives morton ordering according to a true tree 
c  structure. I.e. child morton numbers are at least equal to or
c  greater than their parents number.
      
c 5.1) find maximum level

      max_level_loc = 0
      do i = 1,lnblocks
            max_level_loc = max(lrefine(i),max_level_loc)
            
      end do

      call shmem_barrier_all()

      call comm_int_max_to_all (max_level,max_level_loc)

      call shmem_barrier_all()

c 5) now shift bits to the left by max_levels - level

      do i = 1,lnblocks

         nbitshft = ndim*(max_level-lrefine(i))
#ifdef DEBUG_MORTON
! Test to see if this bitshift will move any non-zero bits beyond
! the leftmost bit.
            inxt_bit = ibits(mort_no(i),nbits-nbitshft,
     .                       nbitshft)
            if(inxt_bit.gt.0.and.nbitshft.gt.0) then
!            if(inxt_bit.gt.0) then
                write(*,*) 'Error 2 in amr_morton_compute',
     .               mort_no(i),inxt_bit,nbits,nbitshft,i,
     .               lrefine(i),max_level
                call amr_abort()
            endif
#endif
         if(nbitshft.gt.0) then
           mort_no(i) = ishft(mort_no(i),nbitshft)

! testing 
! add refine -1 to morton numbers to break degeneracy. Warning there
! may not always been enough space in the morton space for this adjustment.
          mort_no(i) = mort_no(i) + lrefine(i) - 1
! testing 



!#ifdef DEBUG_MORTON
!            nbitshft = ndim*(max_level-lrefine(i))
!            inxt_bit = ibits(mort_no(i),nbits-nbitshft,
!     .                       nbitshft)
!            if(inxt_bit.gt.0) then
!                write(*,*) 'Error 2 in amr_morton_compute',
!     .               mort_no(i),inxt_bit,nbits,nbitshft,i,
!     .               lrefine(i),max_level
!                call amr_abort()
!            endif
!#endif
         endif
      end do

 2    call shmem_barrier_all()

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

      use paramesh_interfaces, only : comm_int_sum_to_all

      implicit none

#include "amr_shmem.fh"

      include 'shmem_reduce.fh'

!      integer, intent(inout??) ::  mort_no(2*maxblocks_tr)
!      integer, intent(inout??) ::  new_loc(2,maxblocks_tr)
      integer, intent(inout) ::  mort_no(:)
      integer, intent(out) ::  new_loc(:,:)
      integer, intent(in) ::  nprocs


      integer shmem_my_pe,mype
      
      integer lnblocks2,tot_blocks,no_per_proc,idi,idp


      integer i,j,k
      integer irnkg(2*maxblocks_tr)
      integer excess,nprocs_y,nprocs_x,irnkg_s


      save tot_blocks,lnblocks2



      real    :: log2,logn
      integer :: ipow

#ifdef SHMEM_USED
      if(nprocs.ne.1) then
        log2 = log(2.)
        logn = log(real(nprocs))
        ipow = int(logn/log2)
        if(nprocs.gt.2**ipow) then
          write(*,*) 'Paramesh error : ',
     .      'you must choose NPROCS to be a power of 2 ',
     .      'when using the SHMEM version'
          call amr_abort()
        endif
      endif
#endif /* SHMEM_USED */

      mype = shmem_my_pe()

      call shmem_barrier_all()

      call amr_bi_sort(mort_no,irnkg,lnblocks)


      call shmem_barrier_all()


c 1) Compute total list length.

! I copy lnblocks to lnblocks2 since lnblocks2 can be put in a save statement.
      lnblocks2 = lnblocks 
      call comm_int_sum_to_all (tot_blocks,lnblocks2)


      no_per_proc = tot_blocks/nprocs

      excess = tot_blocks - no_per_proc*nprocs
      nprocs_y = (no_per_proc+1)*nprocs - tot_blocks
! no. of processors which will get no_per_proc + 1 blocks
      nprocs_x = nprocs - nprocs_y
! rank in list which divides those which go on processor with one number
! of blocks from those which go on another set of blocks w. a different
! no. of blocks
      irnkg_s = nprocs_x*(no_per_proc+1)

      call shmem_barrier_all()
      
c 2) Compute new_locs from rankings (irnkg) returned by amr_bi_sort.
c    The following divides blocks evenly among processors without regard to
c    work.

      if(lnblocks.gt.0) then
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
      endif

      call shmem_barrier_all()


      return
      end subroutine amr_sort_morton

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine amr_sort_by_work (new_loc,nprocs)


c Subroutine to balance work load

c on input takes list sorted by morton ordering
c on output returns new values of new_loc

c By K. Olson (NASA/GSFC and GMU) 11/96

      use paramesh_dimensions
      use physicaldata
      use tree

      use paramesh_interfaces, only : comm_real_sum_to_all

      implicit none

#include "amr_shmem.fh"

      include 'shmem_reduce.fh'

!      integer, intent(inout??) ::  new_loc(2,maxblocks_tr)
      integer, intent(inout) ::  new_loc(:,:)
      integer, intent(in) ::  nprocs



      integer shmem_my_pe,shmem_n_pes,mype
      
      integer lnblocks2,lnblocksl

      real temp

      real work(maxblocks_tr),workt(maxblocks_tr),loc_work,tot_work
      real work_per_proc,work_left(4096) ! dim. set to give max. no. of
                                         ! procs expected

      real wratio
      integer i,j,k,idi
      integer itemp,itemp2
      integer pid(maxblocks_tr),lid(maxblocks_tr),lid2(maxblocks_tr)
      integer pidt,lidt,lid_old
     
      logical,save :: ldatain(1),ldataout(1) 

      logical repeat,repeatt

      save temp,itemp,itemp2,work,workt,work_left
      save idi,lnblocks2,lnblocksl
      save pid,pidt,lid,lidt
      save tot_work,loc_work
      save repeat,repeatt

      integer,save :: idataout(1)


      mype = shmem_my_pe()

c initialize work arrary and temp work arrary

      do i = 1,maxblocks_tr

         work(i) = -1.
         workt(i) = -1.

      end do

c assign values to work array

      do i = 1,maxblocks_tr

         work(i) = 0.
!        work(i) = work_block(i)

#ifdef ADVANCE_ALL_LEVELS
         work(i) = 1.
#else
         if (nodetype(i).eq.1) work(i) = 2.              !<<< USER EDIT
         if (nodetype(i).eq.2) work(i) = 1.              !<<< USER EDIT
#endif

         work_block(i) = work(i)
      end do

      call shmem_barrier_all()

c move work to temp work array workt

      do i = 1,lnblocks

!         temp = work(i)
         call shmem_real_put (workt(new_loc(1,i)),work(i),1,
     &        new_loc(2,i))

      end do

      call shmem_barrier_all()

c copy temp work back to work

      lnblocks2 = 0
      do i = 1,maxblocks_tr
         work(i) = workt(i)
         if (workt(i).gt.-1.) lnblocks2 = lnblocks2 + 1

      end do

      call shmem_barrier_all()

c SUM total work within each processor

      workt(1) = work(1)
      do i = 2,lnblocks2

         workt(i) = workt(i-1) + work(i)

      end do

c SUM work across processors

      loc_work = workt(lnblocks2)
      call shmem_barrier_all()
      call comm_real_sum_to_all (tot_work,loc_work)
      call shmem_barrier_all()                                          ! NEW

c Compute work per processor

      work_per_proc = tot_work/nprocs

c Compute final work by looking left

      do i = 1,nprocs-1
         work_left(i) = 0.
      end do

      call shmem_barrier_all()                                          ! NEW

      do i = 1,nprocs-1

         call shmem_integer_get (idi,lnblocks2,1,i-1)
         call shmem_real_get (work_left(i),workt(idi),1,i-1)
            
      end do

      call shmem_barrier_all()

      do i = 1,mype

         do j = 1,lnblocks2

            workt(j) = workt(j) + work_left(i)

         end do

      end do

      call shmem_barrier_all()

c compute processor ids

      do i = 1,maxblocks_tr

         pid(i) = 0
         lid(i) = 0

      end do

      do i = 1,lnblocks2

         wratio = (workt(i)-1.)/work_per_proc
         pid(i) = int(wratio)
!        pid(i) = int( ((workt(i)-1.)/work_per_proc) )
         if (pid(i).lt.0) pid(i) = 0
         if (pid(i).gt.nprocs-1) pid(i) = nprocs-1

      end do

      call shmem_barrier_all()

c compute local ids
      
      lid(1) = 1
      do i = 2,lnblocks2

         lid(i) = lid(i-1) + 1
         if (pid(i-1).lt.pid(i)) lid(i) = 1  ! start a new group

      end do

      call shmem_barrier_all()


      do i = 1,maxblocks_tr
         lid2(i) = lid(i)
      end do

      pidt = 0
      if (mype.gt.0) then

         ! fetch last proc no from the left
         call shmem_integer_get (lnblocksl,lnblocks2,1,mype-1)
         call shmem_integer_get (pidt,pid(lnblocksl),1,mype-1)

      end if

      call shmem_barrier_all()
      
      lidt = 0
 27   lid_old = lidt ! lid_old stores last fetched value of lid to left

      lidt = 0

      call shmem_barrier_all() ! NEW

      if (mype.gt.0) then
         
         ! fetch last id from the left
         call shmem_integer_get (lnblocksl,lnblocks2,1,mype-1)
         call shmem_integer_get (lidt,lid(lnblocksl),1,mype-1)
         
      end if
      
      call shmem_barrier_all()
      
      do j = 1,lnblocks2
         
         if (pidt.eq.pid(j)) then ! if pidt (which was fetched)
                                  ! equals local pid then the list
                                  ! has been split across processors
               
            lid(j) = lid2(j) + lidt
            
         end if

      end do
      
      call shmem_barrier_all()
      
      repeat = .FALSE.
      if (lidt.ne.lid_old) repeat = .TRUE.
      call shmem_barrier_all() ! NEW
      do j = 0,nprocs-1
         ldatain(1) = repeat
         call shmem_logical_get (ldataout,ldatain,1,j)
         repeatt = ldataout(1)
         if (repeatt) go to 26
      end do
 26   call shmem_barrier_all()
      if (repeatt) go to 27
      
      call shmem_barrier_all()
      
c now reorder according to new pid and lid numberss

      do i = 1,lnblocks

         call shmem_integer_get (idataout(1),lid(new_loc(1,i)),1,
     &        new_loc(2,i))
         new_loc(1,i) = idataout(1)
         
         call shmem_integer_get (idataout(1),pid(new_loc(1,i)),1,
     &        new_loc(2,i))
         new_loc(2,i) = idataout(1)

      if(new_loc(1,i).gt.maxblocks) then
       write(*,*) 'PARAMESH ERROR !'
       write(*,*) 'New block location exceeds MAXBLOCKS limit'
       write(*,*) 'Suggestion: increase MAXBLOCKS or modify',
     .             ' refinement criteria'
       call amr_abort
      endif         
      end do

      call shmem_barrier_all()

      return
      end subroutine amr_sort_by_work

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine amr_migrate_tree_data (new_loc,nprocs)


c Subroutine move tree data and reconnect all pointers given new_loc

c By K. Olson (NASA/GSFC and GMU) 11/96


      use paramesh_dimensions
      use physicaldata
      use tree

      implicit none

#include "amr_shmem.fh"

!      integer, intent(inout) ::  new_loc(2,maxblocks_tr)
      integer, intent(inout) ::  new_loc(:,:)
      integer, intent(in) ::  nprocs


      integer shmem_my_pe,mype
      
      integer neight(2,mfaces,maxblocks_tr)
      integer childt(2,mchild,maxblocks_tr)
      integer which_childt(maxblocks_tr)
      integer parentt(2,maxblocks_tr),lrefinet(maxblocks_tr)
      integer nodetypet(maxblocks_tr)
      integer emptyt(maxblocks_tr),bflagst(mflags,maxblocks_tr)

      real coordt(mdim,maxblocks_tr),sizet(mdim,maxblocks_tr)
      real bnd_boxt(2,mdim,maxblocks_tr)
      real work_blockt(maxblocks_tr)

      real temp,tempb(2)

      integer i,j,k,idi,idp,lb
      integer itemp,cempty,remote_block,remote_pe
      integer ibtemp(mflags)
 
      integer new_loc2(2,maxblocks_tr)

      logical newchildt(maxblocks_tr),ltemp

      save coordt,sizet,neight,childt,parentt,nodetypet,lrefinet
      save newchildt,emptyt,cempty,remote_block,remote_pe
      save new_loc2,temp,tempb,itemp,ltemp,idi,idp
      save bnd_boxt,bflagst,ibtemp,work_blockt,which_childt

      integer,save :: idatain(1)
      real,save    :: rdatain(1)


      mype = shmem_my_pe()

c initialize temp arrays

      do i = 1,maxblocks_tr

         do j = 1,ndim
            coordt(j,i) = -1.
            bnd_boxt(1,j,i) = -1.
            bnd_boxt(2,j,i) = -1.
         end do
         do j = 1,nchild
            childt(1,j,i) = -1
            childt(2,j,i) = -1
         end do
         do j = 1,nfaces
            neight(1,j,i) = -1
            neight(2,j,i) = -1
         end do
         parentt(1,i) = -1
         parentt(2,i) = -1
         do j=1,ndim
            sizet(j,i) = -1.
         enddo
         which_childt(i) = -1
         lrefinet(i) = -1
         nodetypet(i) = -1
         emptyt(i) = 0
         bflagst(:,i) = -1
         work_blockt(i) = 0.
         newchildt(i) = .FALSE.
! new_loc2 is an exact copy of new_loc and is used for later data 
! transferal between processors
         new_loc2(1,i) = new_loc(1,i)
         new_loc2(2,i) = new_loc(2,i)

      end do

      call shmem_barrier_all()

      if(lnblocks.gt.0) then
      do i = 1,lnblocks

         do j = 1,ndim
            rdatain(1) = coord(j,i)
            call shmem_real_put (coordt(j,new_loc(1,i)),rdatain(1),1,
     &           new_loc(2,i))
         end do

         do j = 1,ndim
            tempb(1) = bnd_box(1,j,i)
            tempb(2) = bnd_box(2,j,i)
            call shmem_real_put (bnd_boxt(1,j,new_loc(1,i)),tempb(1)
     &           ,2,new_loc(2,i))
         end do

         do j = 1,nchild
            call shmem_integer_put (childt(1,j,new_loc(1,i)),
     &                      child(1,j,i),2,
     &                      new_loc(2,i))
         end do
         
         do j = 1,nfaces
            call shmem_integer_put (neight(1,j,new_loc(1,i)),
     &           neigh(1,j,i),1,
     &           new_loc(2,i))
            call shmem_udcflush()
            call shmem_integer_put (neight(2,j,new_loc(1,i)),
     &           neigh(2,j,i),1,
     &           new_loc(2,i))
         end do
         
         call shmem_integer_put (parentt(1,new_loc(1,i)),
     &        parent(1,i),1,
     &        new_loc(2,i))
         call shmem_udcflush()
         call shmem_integer_put (parentt(2,new_loc(1,i)),
     &        parent(2,i),1,
     &        new_loc(2,i))

         do j=1,ndim
            rdatain(1) = bsize(j,i)
            call shmem_real_put (sizet(j,new_loc(1,i)),rdatain(1),1,
     &        new_loc(2,i))
         enddo
         
         call shmem_integer_put (lrefinet(new_loc(1,i)),
     &        lrefine(i),1,
     &        new_loc(2,i))
         
         call shmem_integer_put (nodetypet(new_loc(1,i)),
     &        nodetype(i),1,
     &        new_loc(2,i))
         
         call shmem_integer_put (which_childt(new_loc(1,i)),
     &        which_child(i),1,
     &        new_loc(2,i))
         
         call shmem_integer_put (emptyt(new_loc(1,i)),
     &        empty(i),1,
     &        new_loc(2,i))
         
         ibtemp(:) = bflags(:,i)
         call shmem_integer_put (bflagst(1,new_loc(1,i)),ibtemp,
     &        mflags,new_loc(2,i))
         
!         ltemp = newchild(i)
         call shmem_logical_put (newchildt(new_loc(1,i)),
     &        newchild(i),1,
     &        new_loc(2,i))
         
         call shmem_real_put (work_blockt(new_loc(1,i)),
     &        work_block(i),1,new_loc(2,i))



      end do
      endif
      call shmem_barrier_all()

c copy temp arrays back to main arrays

      new_lnblocks = 0
      do i = 1,maxblocks_tr

         do j = 1,ndim
            coord(j,i) = coordt(j,i)
            bnd_box(1,j,i) = bnd_boxt(1,j,i)
            bnd_box(2,j,i) = bnd_boxt(2,j,i)
         end do
         do j = 1,nchild
            child(1,j,i) = childt(1,j,i)
            child(2,j,i) = childt(2,j,i)
         end do
         do j = 1,nfaces
            neigh(1,j,i) = neight(1,j,i)
            neigh(2,j,i) = neight(2,j,i)
         end do
         parent(1,i) = parentt(1,i)
         parent(2,i) = parentt(2,i)
         do j = 1,ndim
            bsize(j,i) = sizet(j,i)
         enddo
         which_child(i) = which_childt(i)
         lrefine(i) = lrefinet(i)
         newchild(i) = newchildt(i)
         nodetype(i) = nodetypet(i)
         empty(i) = emptyt(i)
         bflags(:,i) = bflagst(:,i)
         work_block(i) = work_blockt(i)

         ! the number of blocks on each processor has now changed
         ! due to reordering, so update lnblocks
         if (sizet(1,i).gt.0.) new_lnblocks = new_lnblocks + 1

      end do

      call shmem_barrier_all()

c update pointers to parents, children and neighbors

      do i = 1,new_lnblocks

         ! follow pointers to old locations and fetch back the
         ! the new location of the block that was stored there

         idi = parent(1,i)
         idp = parent(2,i)
         
         if (idi.ge.1) then
            
            call shmem_integer_get (parent(1,i),new_loc2(1,idi),2,idp)
            
         end if

         do k = 1,nchild

            idi = child(1,k,i)
            idp = child(2,k,i)
            
            if (idi.ge.1) then
               
               call shmem_integer_get (child(1,k,i),new_loc2(1,idi),2,
     .       idp)
               
            end if
            
         end do
         
         do k = 1,nfaces

            idi = neigh(1,k,i)
            idp = neigh(2,k,i)
            
            if (idi.ge.1) then
               
               call shmem_integer_get (neigh(1,k,i),new_loc2(1,idi),2,
     .       idp)
               
            end if

         end do

      end do

      call shmem_barrier_all()

! make sure that children of a parent which is marked as empty are also marked
! as empty.
        if(new_lnblocks.gt.0) then
        do lb = 1,new_lnblocks
        if(nodetype(lb).eq.2) then

        if(empty(lb).eq.1) then
       empty(lb) = 2
                cempty = 1

! loop over this blocks children
                do i=1,nchild

! find the location of this child
                remote_block = child(1,i,lb)
                remote_pe    = child(2,i,lb)
                idatain(1) = cempty
                call shmem_integer_put(empty(remote_block),
     .                                 idatain(1),1,remote_pe)

                enddo
        endif

        endif
        enddo
        endif

      call shmem_barrier_all()

#ifdef NEVER
! make sure that children inherit the block flags from their parents
        if(new_lnblocks.gt.0) then
        do lb = 1,new_lnblocks
        if(nodetype(lb).eq.2) then

! loop over this blocks children
                do i=1,nchild

! find the location of this child
                remote_block = child(1,i,lb)
                remote_pe    = child(2,i,lb)
                call shmem_integer_put(bflags(1,remote_block),
     .                           bflags(1,lb),mflags,remote_pe)

                enddo

        endif
        enddo
        endif

      call shmem_barrier_all()
#endif

      return
      end subroutine amr_migrate_tree_data



