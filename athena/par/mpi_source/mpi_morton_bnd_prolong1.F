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

      subroutine mpi_morton_bnd_prolong1(
     &                        mype,nprocs,tag_offset)


!------------------------------------------------------------------------
!
! This routine calculates the morton number for each block on mype.
! It stores the result along with the refinement level of each block into
! the array mortonbnd, and distributes this array among all processors.
!
!
! Written :     Peter MacNeice  and Michael Gehmeyr          February 2000
!------------------------------------------------------------------------
!
! Arguments:
!      mype           rank of local processor
!
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use timings
      use mpi_morton

      use paramesh_mpi_interfaces, only : mpi_amr_write_prol_comm

      implicit none

      include 'mpif.h'

      integer, intent(in)    ::  mype,nprocs
      integer, intent(inout) ::  tag_offset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local variables

      real    :: pbsize(3),pcoord(3)
      real    :: xmin, ymin, zmin
      real    :: xmax, ymax, zmax

      integer :: lb,i,j
      integer :: level,jstack,iprocs
      integer :: morton(6)
      integer :: mort_neigh(6,3,3,3)
      integer :: pmort_neigh(6,3,3,3)
      integer :: neigh_morts(6,3,npts_neigh)
      integer :: t1neigh_morts(6,npts_neigh)
      integer :: t2neigh_morts(6,npts_neigh)
      integer :: istart,iend,indx(npts_neigh)
      integer :: i_pe,j_pe,rem_block,rem_pe
      integer :: no_of_comm_procs
      integer :: ierrorcode,ierr,allocation_status,ierror
      integer :: no_of_remote_neighs
      integer :: max_no_to_be_received
      integer :: no_of_comms_to_send
      integer :: max_no_of_blocks
      integer :: no_of_comms_to_receive
      integer :: istack, ioff, joff, koff, k, itemp, kstack
      integer :: ll, isize, isrc, idest, itag, kk, jj, jp, ip, ii
      integer,dimension (:),  allocatable :: recvrequest
      integer,dimension (:,:),allocatable :: recvstatus
      integer, parameter :: nguarda = max(nguard,nguard_work)

      logical :: lremote,lswap,lfound
      logical :: is_remote,is_found
      logical :: morton_greater_than
      logical :: morton_equal

#ifdef TIMING_MPI
      double precision :: time1
#endif /* TIMING_MPI */
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!
! This routine assumes that the grid blocks are ordered by morton
! number and that any blocks with different refinement levels but
! the same morton number are ordered from coarse to fine.
!
!

#ifdef TIMING_MPI
      time1 = mpi_wtime()
#endif /* TIMING_MPI */
!--------------------------------------------------

! Compute xmin,ymin,zmin,xmax,ymax,zmax or get them from storage
      xmin = grid_xmin
      ymin = grid_ymin
      zmin = grid_zmin
      xmax = grid_xmax
      ymax = grid_ymax
      zmax = grid_zmax

! Initializations
      no_of_comm_procs = 0
      no_of_remote_neighs = 0
      max_no_to_be_received = 0
      max_no_to_send = 0
      commatrix_recv = 0
      commatrix_send = 0
      pe_source = -1
      pe_destination = -1
!     neigh_morts = -1
      no_of_comms_to_send = 0

!--------------------------------------------------
!
! Step 3.
! Construct a list of potential neighbors of all blocks on this
! processor, and potential neighbors of their parents.
! Exclude any which are on this processor.

      istack = 0

#ifdef DEBUG
      write(*,*) 'xmin,ymin,zmin,xmax,ymax,zmax ',
     . xmin,ymin,zmin,xmax,ymax,zmax
#endif /* DEBUG */



      do lb=1,lnblocks

#ifndef ADVANCE_ALL_LEVELS
      if(nodetype(lb).eq.1) then
#endif

!-------------

      pmort_neigh = -1
      mort_neigh  = -1

! Get the possible neighbors of the current block^s parent if it
! is a newchild
      if(newchild(lb)) then
        pbsize(:) = bsize(:,lb)*2.               ! size of parent block
        ioff = mod(which_child(lb)-1,2)        ! coord for parent block
        joff = mod((which_child(lb)-1)/2,2)
        koff = mod((which_child(lb)-1)/4,2)
        if(ioff.eq.0) then
          pcoord(1) = bnd_box(2,1,lb)
        else
          pcoord(1) = bnd_box(1,1,lb)
        endif
        if(joff.eq.0) then
          pcoord(2) = bnd_box(2,2,lb)
        else
          pcoord(2) = bnd_box(1,2,lb)
        endif
        if(ndim.lt.2) pcoord(2) = coord(2,lb)
        if(koff.eq.0) then
          pcoord(3) = bnd_box(2,3,lb)
        else
          pcoord(3) = bnd_box(1,3,lb)
        endif
        if(ndim.lt.3) pcoord(3) = coord(3,lb)

        call morton_neighbors(xmin,ymin,zmin,xmax,ymax,zmax,
     .                        lperiodicx,lperiodicy,lperiodicz,
     .                        pcoord(1),pbsize(1),ndim,
     .                        lrefine(lb)-1,lrefine_max,pmort_neigh)
      endif

!-------------

! If parent is a remote block then puts its address on the list of
! remote blocks which are required.
! The array neigh_morts stores this list. For each block the first
! word is it^s morton number, the second it^s refinement level
! and the third is a flag which tells us why it has been included
! in the list.
      if(newchild(lb).and.parent(2,lb).ne.mype) then
            istack = istack+1
#ifdef DEBUG
            if(istack.gt.npts_neigh) then
              write(*,*) 'morton_bnd_prolong : ',
     .                   'istack exceeds npts_neigh : ',
     .                   'possible solution - increase npts_neigh'
              call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
            endif
#endif /* DEBUG */
            neigh_morts(:,1,istack) = pmort_neigh(:,2,2,2)
            neigh_morts(6,2,istack) = lrefine(lb)-1
            neigh_morts(6,3,istack) = 14
#ifdef DEBUG
      write(*,*) 'parent pmort_neigh(2,2,2) ',pmort_neigh(2,2,2),
     .       ' istack ',istack,' pe ',mype,' parent ',parent(:,lb)
#endif /* DEBUG */
      endif

! Now cycle through the list of parent^s neighbors adding these
! to the list.
      do k = 2-k3d,2+k3d
      do j = 2-k2d,2+k2d
      do i = 1,3
        if(i.ne.2.or.j.ne.2.or.k.ne.2) then

        if(newchild(lb)) then
          if(pmort_neigh(6,i,j,k).gt.-1) then

          lremote = is_remote(pmort_neigh(:,i,j,k),lrefine(lb)-1,mype)

          if(lremote) then 
            istack = istack+1
            neigh_morts(:,1,istack) = pmort_neigh(:,i,j,k)
            neigh_morts(6,2,istack) = lrefine(lb)-1
            neigh_morts(6,3,istack) = (4-i)+((4-j)-1)*3+((4-k)-1)*9
            if(nguarda.gt.nmax_lays) neigh_morts(:,3,istack) = 14
#ifdef DEBUG
            write(*,*) 'pe ',mype,' blk ',lb,' ijk ',i,j,k,
     .             ' neigh_morts ',
     .              neigh_morts(:,istack),' istack ',istack
#endif /* DEBUG */
          else
#ifdef DEBUG
            write(*,*) 'pe ',mype,' blk ',lb,' ijk ',
     .         i,j,k,' parent neigh local'
#endif /* DEBUG */
          endif
          endif   
        endif


        endif
      enddo
      enddo
      enddo


!-------------

! Now get the possible neighbors of the current block.
!
! We will only need to add these to the final list of remote blocks
! if they are not new children, in which case they will may be
! needed for div B consistency checking. However it is difficult to
! determine whether they are old or new, and so we do not currently
! implement this test.
! 
      if(newchild(lb)) then

      mort_neigh = -1
      call morton_neighbors(xmin,ymin,zmin,xmax,ymax,zmax,
     .                      lperiodicx,lperiodicy,lperiodicz,
     .                      coord(1,lb),bsize(1,lb),ndim,
     .                      lrefine(lb),lrefine_max,mort_neigh)

      do k = 2-k3d,2+k3d
      do j = 2-k2d,2+k2d
      do i = 1,3
        if(i.ne.2.or.j.ne.2.or.k.ne.2) then

! if neighbor block exists at this refinement level
        if(mort_neigh(6,i,j,k).gt.-1) then
         
          lremote = is_remote(mort_neigh(:,i,j,k),lrefine(lb),mype)
          
          if(lremote) then 
            istack = istack+1
            neigh_morts(:,1,istack) = mort_neigh(:,i,j,k)
            neigh_morts(6,2,istack) = lrefine(lb)
            neigh_morts(6,3,istack) = (4-i)+((4-j)-1)*3+((4-k)-1)*9
            if(nguarda.gt.nmax_lays) neigh_morts(6,3,istack) = 14
#ifdef DEBUG
            write(*,*) 'pe ',mype,' blk ',lb,' ijk ',i,j,k,
     .             ' neigh_morts ',
     .              neigh_morts(:,istack),' istack ',istack
#endif /* DEBUG */
          else
#ifdef DEBUG
            write(*,*) 'pe ',mype,' blk ',lb,' ijk ',i,j,k,' local'
#endif /* DEBUG */
          endif
        endif


        endif
      enddo
      enddo
      enddo

      endif
!-------------

#ifndef ADVANCE_ALL_LEVELS
      endif
#endif

      enddo


!--------------------------------------------------
       if(istack.gt.0) then
!--------------------------------------------------
!
! Step 4.
! Compress this list by removing any redundancies. Do this by
! sorting the list in order of increasing morton number, and then
! sorting each sub-list with the same morton number in order of
! increasing refinement level. Then remove any identical entries.


! set indx so neigh_morts(3,:) can be permuted.
       do i=1,istack
         indx(i) = i
       enddo

! sort the neigh_morts list according to their morton numbers
       t1neigh_morts(:,1:istack) = neigh_morts(:,1,1:istack)
       call morton_sort(t1neigh_morts(:,1:istack),indx(1:istack),istack)
! now reorder the morton number part of neigh_morts
       neigh_morts(:,1,1:istack) = t1neigh_morts(:,1:istack)

       t2neigh_morts(6,1:istack) = neigh_morts(6,2,1:istack)
       do i=1,istack
         neigh_morts(6,2,i) = t2neigh_morts(6,indx(i))
       enddo

       t2neigh_morts(6,1:istack) = neigh_morts(6,3,1:istack)
       do i=1,istack
         neigh_morts(6,3,i) = t2neigh_morts(6,indx(i))
       enddo

! now scan the list identifying segments with the same morton number
! and sort each segment based on refinement level
! order segments with same morton number in order of increasing
! refinement level
      lswap = .true.
      do while (lswap)
        lswap = .false.
        do i = 1,istack-1
          if(morton_equal(neigh_morts(:,1,i),neigh_morts(:,1,i+1))
     .                          .and.
     .       neigh_morts(6,2,i).gt.neigh_morts(6,2,i+1) ) then
            lswap = .true.
            itemp = neigh_morts(6,2,i)
            neigh_morts(6,2,i) = neigh_morts(6,2,i+1)
            neigh_morts(6,2,i+1) = itemp
            itemp = neigh_morts(6,3,i)
            neigh_morts(6,3,i) = neigh_morts(6,3,i+1)
            neigh_morts(6,3,i+1) = itemp
#ifdef DEBUG 
            write(*,*) 'pe ',mype,' swapping ',i,i+1
#endif /* DEBUG  */
          endif
        enddo
      enddo                           ! end do while

!
! If any entries have the same morton number and refinement level
! but multiple data request types then mark them to fetch the complete
! blocks.
      istart = 1
      iend = 1
      i = 2
      do while(i.le.istack)
        morton(:) = neigh_morts(:,1,i-1)
        level  = neigh_morts(6,2,i-1)

        do while(morton_equal(neigh_morts(:,1,i),morton(:)) .and.
     .           neigh_morts(6,2,i).eq.level   .and.
     .           i.le.istack)
          i = i+1
        enddo
        iend = i-1
        if(istart.lt.iend) then

          call rationalize_list(neigh_morts,istart,iend)

!          lfullblock = .false.
!          jmin = npts_neigh+1       ! set to a bigger number than is possible
!          jmax = -1
!          do j = istart,iend
!            if(neigh_morts(3,j).eq.14) lfullblock = .true.
!            jmin = min(jmin,neigh_morts(3,j))
!            jmax = max(jmax,neigh_morts(3,j))
!          enddo
!          if(jmin.lt.jmax) lfullblock = .true.
!          if(lfullblock) neigh_morts(3,istart:iend) = 14
        endif
        istart = iend+1
        i = istart+1

      enddo


!
! Remove any entries which do not fall between the morton limits on
! any processor.
! Any non-existant blocks are marked with neigh_morts=(-1,-1).
        do i = 1,istack
          j_pe = 0
          lfound = .false.
          do while(.not.lfound.and.j_pe.lt.nprocs)
            j_pe = j_pe + 1
! skip the local processor while searching since any morton nos. found
! on the local proc are not required in the list of remote blocks to b
! fetched.
            if(j_pe.ne.mype+1) then

            lfound = 
     .        is_found(neigh_morts(:,1,i),neigh_morts(6,2,i),j_pe)

            endif
          enddo
          if(.not.lfound) neigh_morts(:,:,i) = -1
        enddo

! Remove any non-existant entries.
      indx = 0
      jstack = 0
      do i=1,istack
        if(neigh_morts(6,1,i).gt.-1) then
          jstack = jstack+1
          indx(jstack) = i
        endif
      enddo
      do j=1,jstack
        neigh_morts(:,:,j) = neigh_morts(:,:,indx(j))
      enddo
      neigh_morts(:,:,jstack+1:istack) = -1
      istack = jstack


#ifdef DEBUG
      do i=1,istack
      write(*,*) 'pe ',mype,' stack ',i,' neigh_morts ',
     .            neigh_morts(:,i)
      enddo
#endif /* DEBUG  */


! Finally remove any repetition of elements in this list.
       jstack = 0
       morton = -1
       level  = -1
       do i = 1,istack
         if( any(neigh_morts(:,1,i).ne.morton(:)) .or.
     .           neigh_morts(6,2,i).ne.level ) then
           jstack = jstack + 1
           neigh_morts(:,:,jstack) = neigh_morts(:,:,i)
           morton(:) = neigh_morts(:,1,i)
           level  = neigh_morts(6,2,i)
#ifdef DEBUG
         write(*,*) 'compress: pe ',mype,' jstack ',jstack,
     .      ' morton ',morton,' level ',level
#endif /* DEBUG  */
         endif
       enddo

       no_of_remote_neighs = jstack

#ifdef DEBUG
      do i=1,jstack
      write(*,*) 'pe ',mype,' jstack ',i,' neigh_morts ',
     .            neigh_morts(:,i)
      enddo
#endif /* DEBUG  */
!--------------------------------------------------
!
! Step 5.
! Construct a list of all processors from which the local processor should
! request morton number information.


! non-zero elements of COMMATRIX define which processor pairs need to 
! exchange morton number lists.

        do i = 1,no_of_remote_neighs
          i_pe = 1
          j_pe = -1
          do while( 
     .       ( morton_greater_than(
     .          neigh_morts(:,1,i),morton_limits(:,1,2,i_pe))
     .                               .or.
     .       (morton_equal(neigh_morts(:,1,i),morton_limits(:,1,2,i_pe))
     .                              .and.
     .        neigh_morts(6,2,i).gt.morton_limits(6,2,2,i_pe)  )  )
     .        .and. (i_pe.le.nprocs)
     .            )
             i_pe = i_pe + 1
             if (i_pe > nprocs) exit
          enddo
          if(i_pe.le.nprocs) j_pe = i_pe
!
! If block has been located then update commatrix
          if(j_pe.ne.-1) 
     .      commatrix_recv(j_pe) =  commatrix_recv(j_pe) + 1

        enddo

#ifdef DEBUG
        write(*,*) 'pe ',mype,' commatrix bef gather ',
     .             commatrix(1:nprocs,1:nprocs)
#endif /* DEBUG  */


! record the number of processors which will communicate with the
! local processor.
       no_of_comms_to_send = 0
       kstack = 0
       do i = 1,nprocs
         no_of_comms_to_send = no_of_comms_to_send +
     .                          min( 1, commatrix_recv(i) )
         if(commatrix_recv(i).gt.0) then
           kstack = kstack+1
           pe_source(kstack) = i
         endif
       enddo
#ifdef DEBUG
       write(*,*) 'pe ',mype,' no_of_comms_to_send ',
     .           no_of_comms_to_send
#endif /* DEBUG  */

!--------------------------------------------------
       endif                     ! end of istack if test
!--------------------------------------------------
!
! Step 6.
! provide the complete COMMATRIX to all processors


      call MPI_AlltoAll (commatrix_recv,       1,MPI_INTEGER,
     .                   commatrix_send,       1,MPI_INTEGER,
     .                   MPI_COMM_WORLD,ierror)

#ifdef DEBUG
        write(*,*) 'pe ',mype,' commatrix ',
     .             commatrix(1:nprocs,1:nprocs)
 
      if(mype.eq.0) then
         write(*,'(" ")')
         write(*,'(" COMMUNICATION MATRIX: morton_bnd_prolong")')
         write(*,'(" ")')
         do ipe=1,nprocs
         write(*,'(" ",8i3)') (commatrix(i,ipe),i=1,nprocs)
         enddo
         write(*,'(" ")')
      endif

#endif /* DEBUG  */
!--------------------------------------------------
!
! Step 7.
! Compute the maximum amount of morton information which any processor
! is going to receive.\

      max_no_to_be_received = 0
 
      iprocs = 0
      do j = 1,nprocs
         iprocs = iprocs + min(1,commatrix_recv(j))
      enddo
      max_no_to_be_received = max(1,iprocs)

#ifdef DEBUG
       write(*,*) 'pe ',mype,' max_no_to_be_received ',
     .           max_no_to_be_received
#endif /* DEBUG  */


!--------------------------------------------------
!
! Step 8.
! Dynamically allocate memory to store the remote morton information.


      call MPI_ALLREDUCE(lnblocks,
     .                   max_no_of_blocks,
     .                   1,
     .                   MPI_INTEGER,
     .                   MPI_MAX,
     .                   MPI_COMM_WORLD,
     .                   ierror)

       if(allocated(r_mortonbnd)) deallocate(r_mortonbnd)
       allocate( r_mortonbnd(6,3,max_no_of_blocks,
     .                       max(1,max_no_to_be_received) ),
     .           stat = allocation_status)
       if(allocation_status > 0) then
          write(*,*) 'morton_bnd : allocation error'
          call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
       endif

!--------------------------------------------------

       if(allocated(recvrequest)) deallocate( recvrequest )
       allocate ( recvrequest(nprocs) )
       if(allocated(recvstatus)) deallocate( recvstatus )
       allocate ( recvstatus(MPI_STATUS_SIZE,nprocs) )

!
! Step 9.
! Exchange morton information between processors.


      ll = 0

      r_mortonbnd = -1
      pe_source   = -1
      isize = 3*max_no_of_blocks*6
      k = 0

      do i = 1,nprocs
         isrc = i-1
         idest= mype
         itag = isrc*nprocs + idest+1 + tag_offset
                                 ! receive to pe=j
         if(commatrix_recv(i).gt.0) then
            k = k+1
            pe_source(k) = isrc+1
            call Mpi_Irecv(r_mortonbnd(1,1,1,k),isize,MPI_INTEGER,
     .           isrc ,itag,MPI_COMM_WORLD,recvrequest(k),ierr)
         endif
      enddo
      
      do j = 1,nprocs
         isrc = mype
         idest= j-1
         itag = isrc*nprocs + idest+1 + tag_offset
                                ! send from mype=i
         if(commatrix_send(j).gt.0) then
            ll = ll+1
            call MPI_Ssend(mortonbnd(1,1,1),isize,MPI_INTEGER,
     .           idest,itag,MPI_COMM_WORLD,ierr)
         endif
      enddo

      no_of_mortonbnds_received = k

      tag_offset = (nprocs-1)*nprocs + nprocs + tag_offset

      if(k.gt.0)
     .    call MPI_Waitall(k,recvrequest,recvstatus,
     .                     ierror)


#ifdef DEBUG
      write(*,*) 'pe ',mype,' no_of_mortonbnds_received ',
     .          no_of_mortonbnds_received
      write(*,*) 'pe ',mype,' r_mortonbnd(:,1:15,1) ',
     .          r_mortonbnd(:,1:15,1)
#endif /* DEBUG  */
        
!--------------------------------------------------
!
! Step 10.
! Loop over this processor^s list of required neighbor blocks,
! identifying their remote location from the morton information received
! in step 9.


        do i = 1,no_of_remote_neighs
          i_pe = 1
          j_pe = -1
          do while( 
     .      (morton_greater_than(
     .         neigh_morts(:,1,i),morton_limits(:,1,2,i_pe))
     .                             .or.
     .      (morton_equal(neigh_morts(:,1,i),morton_limits(:,1,2,i_pe))
     .                             .and.
     .       neigh_morts(6,2,i).gt.morton_limits(6,2,2,i_pe)  )  )
     .       .and. (i_pe.le.nprocs)
     .            )
            i_pe = i_pe + 1
            if (i_pe > nprocs) exit
          enddo
          if(i_pe.le.nprocs) j_pe = i_pe

          rem_block = -1
          rem_pe = j_pe

          kk = -1
          do k=1,no_of_mortonbnds_received
            if(pe_source(k).eq.rem_pe) kk = k 
          enddo
          if(kk.gt.0) then

          do j=1,max_no_of_blocks
            if( morton_equal(r_mortonbnd(:,1,j,kk),neigh_morts(:,1,i))
     .                              .and.
     .          r_mortonbnd(6,2,j,kk).eq.neigh_morts(6,2,i) )
     .          rem_block = j
          enddo
          endif
          if(rem_block.eq.-1) rem_pe = -1

#ifdef DEBUG 
          write(*,*) 'pe ',mype,' neigh i ',i,' rem_pe ',
     .            rem_pe,' kk ',kk,' rem_block ',rem_block
#endif /* DEBUG  */

! neigh_morts(1:2,no_of_remote_neighs) is now being used to store 
! the remote addresses of the required neighbors.
! Here proc nos. run from 1 to nprocs.
          neigh_morts(6,1,i) = rem_block
          if(neigh_morts(6,1,i).lt.0) then
            neigh_morts(6,2,i) = -1
          else
            neigh_morts(6,2,i) = rem_pe
          endif

#ifdef DEBUG 
          write(*,*) 'pe ',mype,' neigh i ',i,' address ',
     .            neigh_morts(:,i)
#endif /* DEBUG */
        enddo

!--------------------------------------------------
!
! Step 11.
! Check for any non-existent blocks in the neigh_morts list
! and remove them. Then reset commatrix.

      indx = 0
      jstack = 0
      do i=1,no_of_remote_neighs
        if(neigh_morts(6,1,i).gt.-1) then
          jstack = jstack+1
          indx(jstack) = i
        endif
      enddo
      do j=1,jstack
        neigh_morts(6,:,j) = neigh_morts(6,:,indx(j))
      enddo
      neigh_morts(6,:,jstack+1:istack) = -1
      istack = jstack
      no_of_remote_neighs = istack

!--------------------------------------------------

! Step 12.
! Reconstruct commatrix.


! non-zero elements of COMMATRIX define which processor pairs need to 
! exchange morton number lists. 
#ifdef DEBUG 
          write(*,*) 'pe ',mype,' revised no_of_remote_neighs ',
     .        no_of_remote_neighs
#endif /* DEBUG */
        commatrix_recv = 0
        commatrix_send = 0
        do i = 1,no_of_remote_neighs
          i_pe = neigh_morts(6,2,i)
#ifdef DEBUG 
          write(*,*) 'pe ',mype,' add item to commatrix(',i_pe,
     .     mype+1,')'
#endif /* DEBUG */
          commatrix_recv(i_pe) =  commatrix_recv(i_pe) + 1
        enddo

#ifdef DEBUG
       write(*,*) 'pe ',mype,' no_of_comms_to_send-xx ',
     .           no_of_comms_to_send
#endif /* DEBUG  */
!
! Eliminate any r_mortonbnds layers which are no longer required.
        jstack = 0
        do i = 1,no_of_comms_to_send
          i_pe = pe_source(i)
          if(commatrix_recv(i_pe).gt.0) then
            jstack = jstack+1
            indx(jstack) = i
          endif
        enddo
        do j=1,jstack
          r_mortonbnd(:,:,:,j) = r_mortonbnd(:,:,:,indx(j))
        enddo
        no_of_mortonbnds_received = jstack            
#ifdef DEBUG
      write(*,*) 'pe ',mype,' revised no_of_mortonbnds_received ',
     .          no_of_mortonbnds_received
#endif /* DEBUG  */

! record the number of processors which will communicate with the
! local processor.
       pe_source = -1
       no_of_comms_to_send = 0
       kstack = 0
       do i = 1,nprocs
         no_of_comms_to_send = no_of_comms_to_send +
     .                          min( 1, commatrix_recv(i) )
         if(commatrix_recv(i).gt.0) then
           kstack = kstack+1
           pe_source(kstack) = i
         endif
       enddo
#ifdef DEBUG
       write(*,*) 'pe ',mype,' no_of_comms_to_send ',
     .           no_of_comms_to_send
#endif /* DEBUG  */

!--------------------------------------------------
!
! Step 13.
! Repeat Step 6.
! provide the complete COMMATRIX to all processors


      call MPI_AlltoAll (commatrix_recv,       1,MPI_INTEGER,
     .                   commatrix_send,       1,MPI_INTEGER,
     .                   MPI_COMM_WORLD,ierror)

#ifdef DEBUG
        write(*,*) 'pe ',mype,' commatrix ',
     .             commatrix(1:nprocs,1:nprocs)

      if(mype.eq.0) then
         write(*,'(" ")')
         write(*,'(" COMMUNICATION MATRIX2: morton_bnd_prolong")')
         write(*,'(" ")')
         do ipe=1,nprocs
         write(*,'(" ",8i3)') (commatrix(i,ipe),i=1,nprocs)
         enddo
         write(*,'(" ")')
      endif
#endif /* DEBUG  */


!--------------------------------------------------

! Step 14.
! record the number of processors to which the local processor
! will send messages.

       no_of_comms_to_receive = 0
       kstack = 0
       do i = 1,nprocs
         no_of_comms_to_receive = no_of_comms_to_receive +
     .                          min( 1, commatrix_send(i) )
         if(commatrix_send(i).gt.0) then
           kstack = kstack+1
           pe_destination(kstack) = i
         endif
       enddo
#ifdef DEBUG
       write(*,*) 'pe ',mype,' no_of_comms_to_receive ',
     .           no_of_comms_to_receive
#endif /* DEBUG  */


!--------------------------------------------------
!
! Step 15.
! Compute the maximum amount of morton information which any processor
! is going to receive.

       max_no_to_be_received = 0
 
       iprocs = 0
       do j = 1,nprocs
          iprocs = iprocs + min(1,commatrix_recv(j))
       enddo
       max_no_to_be_received = max(1,iprocs)

#ifdef DEBUG
       write(*,*) 'pe ',mype,' max_no_to_be_received ',
     .           max_no_to_be_received
#endif /* DEBUG  */

!--------------------------------------------------
!
! Step 16.
! Compute the maximum amount of information which any processor
! is going to receive.

       iprocs = 0
       do j = 1,nprocs
          iprocs = iprocs + min(1,commatrix_send(j))
       enddo
       max_no_to_send = max(1,iprocs)

#ifdef DEBUG
       write(*,*) 'pe ',mype,' max_no_to_send ',
     .           max_no_to_send
#endif /* DEBUG  */

!--------------------------------------------------
!
! Step 17.
! evaluate smallest guard block starting index over all pe
! store this into variable strt_buffer which is used in amr_1blk_guardcell

      last_buffer = maxblocks_alloc

      k = last_buffer
      do i=0,nprocs-1
      k = k - commatrix_recv(i+1)
      enddo
      strt_buffer = k + 1


      if (strt_buffer.le.lnblocks) then
        write(*,*) 
     .  'ERROR in mpi_morton_bnd_prolong1:',
     .  ' guard block starting index',
     .  strt_buffer,' not larger than lnblocks',lnblocks
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif


!--------------------------------------------------
!
! Step 18.
! Dynamically allocate memory to store the lists of blocks to be
! sent and received.

      iprocs = max(maxval(commatrix_recv),maxval(commatrix_send))
      call MPI_ALLREDUCE(iprocs,
     .                   largest_no_of_blocks,
     .                   1,
     .                   MPI_INTEGER,
     .                   MPI_MAX,
     .                   MPI_COMM_WORLD,
     .                   ierror)

       if(allocated(to_be_sent)) deallocate(to_be_sent)
       if(allocated(to_be_received)) deallocate(to_be_received)
       allocate( to_be_sent(3,
     .                           max(1,largest_no_of_blocks),
     .                           max(1,max_no_to_send) ) )
       allocate( to_be_received(3,
     .                          max(1,largest_no_of_blocks),
     .                          max(1,max_no_to_be_received) ) )

!--------------------------------------------------
!
! Step 19.
! Construct arrays to_be_sent and to_be_received which contain
! the lists of blocks to be packaged.

        to_be_sent = -1
        to_be_received = -1
        laddress = 0

! First set up the array to_be_received on each processor
        if(no_of_remote_neighs.gt.0) then

          jj = 0
          do jp = 1,no_of_mortonbnds_received
            ip = pe_source(jp)
            if(commatrix_recv(ip).gt.0) then   ! this is a needless check
              do ii = 1,commatrix_recv(ip)
                jj = jj+1
                if(neigh_morts(6,2,jj).eq.ip) then
                  if(ii.gt.largest_no_of_blocks) then
          write(*,*) 'pe ',mype,' ii too large ',ii
                  endif
                  to_be_received(:,ii,jp) = neigh_morts(6,:,jj)
                endif
              enddo 
            endif
          enddo

          laddress(1,strt_buffer:strt_buffer+jj-1) =
     .          neigh_morts(6,1,1:jj)
          laddress(2,strt_buffer:strt_buffer+jj-1) =
     .          neigh_morts(6,2,1:jj)-1

        endif

#ifdef DEBUG
        do jp = 1,no_of_mortonbnds_received
        write(*,*) 'pe ',mype,' jreceive ',jp,' to_be_received ',
     .    to_be_received(:,:,jp)
        enddo
        write(*,*) 'pe ',mype,' laddress to_be_received ',
     .    laddress
#endif /* DEBUG */


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Now exchange info in to_be_received with the sending processors
! to construct the equivalent to_be_sent arrays

! Post receives 
        isize = 3*largest_no_of_blocks
        kk = 0
        do i = 1,nprocs
          isrc = i-1
          idest= mype
          itag = isrc*nprocs + idest+1 + tag_offset
                                 ! receive to pe=j
          if(commatrix_send(i).gt.0) then
            kk = kk+1
#ifdef DEBUG
            write(*,*) 'post receive on pe ',mype,' from ',isrc
#endif /* DEBUG */
            call Mpi_Irecv(to_be_sent(1,1,kk),isize,MPI_INTEGER,
     .           isrc ,itag,MPI_COMM_WORLD,recvrequest(kk),ierr)
#ifdef DEBUG
            write(*,*) 'receive posted on pe ',mype,' from ',isrc,
     .      ' into layer ',kk,' tag ',itag
#endif /* DEBUG */
          endif
        enddo

! Post sends

        ll = 0
        do j = 1,nprocs
          isrc = mype
          idest= j-1
          itag = isrc*nprocs + idest+1 + tag_offset
                                 ! send from mype=i
          if(commatrix_recv(j).gt.0) then
            ll = ll+1
            call MPI_Ssend(to_be_received(1,1,ll),isize,MPI_INTEGER,
     .           idest,itag,MPI_COMM_WORLD,ierr)
#ifdef DEBUG
            write(*,*) 'post send from pe ',mype,' to ',idest,
     .      ' from layer ',ll,' tag ',itag
#endif /* DEBUG */
          endif
        enddo

        tag_offset = (nprocs-1)*nprocs + nprocs + tag_offset

        if(kk.gt.0)
     .    call MPI_Waitall(kk,recvrequest,recvstatus,
     .                     ierror)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef DEBUG
        do jp = 1,kk
        write(*,*) 'pe ',mype,' jsend ',jp,' to_be_sent ',
     .    to_be_sent(:,:,jp)
        enddo
#endif /* DEBUG */
!--------------------------------------------------
!
! Step 20.
! Deallocate any memory which was dynamically allocated for local use in this
! routine.


       if(allocated(recvrequest)) deallocate( recvrequest )
       if(allocated(recvstatus)) deallocate( recvstatus )


!--------------------------------------------------

! Mark morton data up to date
       morton_limits_set = .true.


! Store communication info for future use
       call mpi_amr_write_prol_comm(nprocs)


#ifdef DEBUG
      write(*,*) 'pe ',mype,' exiting mpi_morton_bnd_prol'
#endif /* DEBUG */


#ifdef TIMING_MPI
              timer_mpi_morton_bnd_prolong1 =  
     .                          timer_mpi_morton_bnd_prolong1
     .                          + mpi_wtime() - time1
#endif /* TIMING_MPI */
      return
      end subroutine mpi_morton_bnd_prolong1
