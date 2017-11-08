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
!#define DEBUGX

      subroutine mpi_morton_bnd_fluxcon(
     &                        mype,nprocs,tag_offset)


!------------------------------------------------------------------------
!
! This routine calculates the morton number for each block on mype.
! It stores the result along with the refinement level of each block into
! the array mortonbnd, and distributes this array among all processors.
!
!
! Written :     Peter MacNeice  and Michael Gehmeyr          July 2000
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

      use paramesh_mpi_interfaces, only : mpi_amr_write_flux_comm

      implicit none

      include 'mpif.h'

      integer, intent(in)    ::  mype,nprocs
      integer, intent(inout) ::  tag_offset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local variables

      real    :: xmin, ymin, zmin
      real    :: xmax, ymax, zmax

      integer :: lb,i,j,ipe
      integer :: level,jstack,iprocs
      integer :: morton(6)
      integer :: mort_neigh(6,3,3,3)
      integer :: neigh_morts(6,3,npts_neigh)
      integer :: t1neigh_morts(6,npts_neigh)
      integer :: t2neigh_morts(6,npts_neigh)
      integer :: istart,iend,indx(npts_neigh)
      integer :: i_pe,j_pe,rem_block,rem_pe
      integer :: no_of_comm_procs
      integer :: ierrorcode,ierr,allocation_status,ierror
      integer :: l_nodetype(3,3,3)
      integer :: no_of_remote_neighs
      integer :: max_no_to_be_received
      integer :: no_of_comms_to_send
      integer :: max_no_of_blocks
      integer :: no_of_comms_to_receive
      integer :: istack, k, itemp, kstack, isize, isrc, idest
      integer :: itag, ll, kk, ij, ie, if1, jf1, kf1, if2, jf2, kf2
      integer :: nodetype_1, nodetype_2, jj, jp, ip, ii
      integer,dimension (:),  allocatable :: recvrequest
      integer,dimension (:,:),allocatable :: recvstatus

      logical :: lswap,lfound
      logical :: is_found 
      logical :: morton_greater_than
      logical :: morton_less_than
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


#ifdef TIMING_MPI
      time1 = mpi_wtime()
#endif /* TIMING_MPI */
!--------------------------------------------------
!
! Step 1.
!

! Compute xmin,ymin,zmin,xmax,ymax,zmax or get them from storage
      xmin = grid_xmin
      ymin = grid_ymin
      zmin = grid_zmin
      xmax = grid_xmax
      ymax = grid_ymax
      zmax = grid_zmax

!
! Initializations
      no_of_comm_procs = 0
      no_of_remote_neighs = 0
      max_no_to_be_received = 0
      max_no_to_send = 0
      commatrix_send = 0
      commatrix_recv = 0
      pe_source = -1
      pe_destination = -1
!     neigh_morts = -1
      no_of_comms_to_send = 0


!--------------------------------------------------
!
! Step 2.

!
! We want to construct a list of neighbors of leaf blocks which
! are refined. If we compute the morton number for a neighbor
! of a leaf block, it will be the same as the morton number of
! its leftmost child. Therefore if we look for the (mort,level)
! combination for this lower left child, finding it tells us that
! the leaf blocks neighbor is refined, and so we should add the
! neighbors (mort,level) to the list.

! Construct a list of the bottom left child of all neighbors of
! parents of leaf blocks. This list should include both remote
! and local children.

      istack = 0

#ifdef DEBUG
      write(*,*) 'fluxcon: xmin,ymin,zmin,xmax,ymax,zmax ',
     . xmin,ymin,zmin,xmax,ymax,zmax
#endif /* DEBUG */

      do lb=1,lnblocks

      if(nodetype(lb).eq.1) then

!-------------

!
! Get morton numbers for neighbors of any leaf blocks
      mort_neigh  = -1
      call morton_neighbors(xmin,ymin,zmin,xmax,ymax,zmax,
     .                      lperiodicx,lperiodicy,lperiodicz,
     .                      coord(1,lb),bsize(1,lb),ndim,
     .                      lrefine(lb),lrefine_max,mort_neigh)

!-------------

! Now cycle through this list of neighbors constructing the
! (morton no., refinement level) for their bottom left child.
      do k = 2-k3d,2+k3d
      do j = 2-k2d,2+k2d
      do i = 1,3
        if(i.ne.2.or.j.ne.2.or.k.ne.2) then

          if(mort_neigh(6,i,j,k).gt.-1) then
            istack = istack+1
#ifdef DEBUG
            if(istack.gt.npts_neigh) then
              write(32,*) 'morton_bnd_fluxcon : ',
     .                   'istack exceeds npts_neigh : ',
     .                   'possible solution - increase npts_neigh'
              call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
            endif
#endif /* DEBUG */
            neigh_morts(:,1,istack) = mort_neigh(:,i,j,k)
            neigh_morts(6,2,istack) = lrefine(lb)+1
            neigh_morts(6,3,istack) = (4-i)+((4-j)-1)*3+((4-k)-1)*9
            if(nguard.gt.nmax_lays) neigh_morts(6,3,istack) = 14
#ifdef DEBUGY
            write(*,*) 'flx:pe ',mype,' blk ',lb,' ijk ',i,j,k,
     .             ' neigh_morts ',
     .              neigh_morts(:,istack),' istack ',istack
#endif /* DEBUG */
          endif   

        endif
      enddo
      enddo
      enddo


      endif

      enddo


!--------------------------------------------------
       if(istack.gt.0) then
!--------------------------------------------------
!
! Step 3.
! Compress this list by removing any redundancies. Do this by
! sorting the list in order of increasing morton number, and then
! sorting each sub-list with the same morton number in order of
! increasing refinement level. Then remove any identical entries.

! set indx so neigh_morts(3,:) can be permuted.
       do i=1,istack
         indx(i) = i
       enddo

! sort the neigh_morts list according to their morton numbers
       do i = 1,6
          t1neigh_morts(i,1:istack) = neigh_morts(i,1,1:istack)
       end do
       call morton_sort(t1neigh_morts(:,1:istack),indx(1:istack),istack)
! now reorder the morton number part of neigh_morts
       do i = 1,6
          neigh_morts(i,1,1:istack) = t1neigh_morts(i,1:istack)
       end do

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
     .                         .and.
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

! If any entries have the same morton number and refinement level
! but multiple data request types then mark them to fetch the complete
! blocks.
      istart = 1
      iend = 1
      i = 2
      do while(i.le.istack)

        morton(:) = neigh_morts(:,1,i-1)
        level  = neigh_morts(6,2,i-1)

        do while(morton_equal(neigh_morts(:,1,i),morton(:)) 
     .                     .and.
     .           neigh_morts(6,2,i).eq.level   .and.
     .           i.le.istack)
          i = i+1
        enddo
        iend = i-1
        if(istart.lt.iend) then

          call rationalize_list(neigh_morts,istart,iend)

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
! Do not skip the local processor while searching since we will need this
! info to construct edge_mark. Local entries in the neigh_morts list
! will be removed once edge_mark is set up.
           lfound = 
     .        is_found(neigh_morts(:,1,i),neigh_morts(6,2,i),j_pe)
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


#ifdef DEBUGY
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

#ifdef DEBUGY
      do i=1,jstack
      write(*,*) 'flux_bnd pe ',mype,' jstack ',i,' neigh_morts ',
     .            neigh_morts(:,i)
      enddo
#endif /* DEBUG  */
!--------------------------------------------------
!
! Step 4.
! Construct a list of all processors from which the local processor should
! request morton number information.

! non-zero elements of COMMATRIX define which processor pairs need to 
! exchange morton number lists.

        do i = 1,no_of_remote_neighs
          i_pe = 1
          j_pe = -1

          do while( 
     .       ( morton_greater_than(neigh_morts(:,1,i),
     .                             morton_limits(:,1,2,i_pe))
     .                               .or.
     .         (morton_equal(neigh_morts(:,1,i),
     .                       morton_limits(:,1,2,i_pe)).and.
     .          neigh_morts(6,2,i).gt.morton_limits(6,2,2,i_pe)  )  )
     .          .and. (i_pe.le.nprocs)
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
!
! In this case some elements in the list might be local blocks. However
! we do not wish these to contribute to commatrix.
        commatrix_recv(mype+1) = 0

#ifdef DEBUGY
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
! Step 5.
! provide the complete COMMATRIX to all processors

      call MPI_AlltoAll (commatrix_recv,       1,MPI_INTEGER,
     .                   commatrix_send,       1,MPI_INTEGER,
     .                   MPI_COMM_WORLD,ierror)
       
#ifdef DEBUGY
        write(*,*) 'pe ',mype,' commatrix ',
     .             commatrix(1:nprocs,1:nprocs)
 
      if(mype.eq.0) then
         write(*,'(" ")')
         write(*,'(" COMMUNICATION MATRIX: morton_bnd_fluxcon")')
         write(*,'(" ")')
         do ipe=1,nprocs
         write(*,'(" ",8i3)') (commatrix(i,ipe),i=1,nprocs)
         enddo
         write(*,'(" ")')
      endif

#endif /* DEBUG  */
!--------------------------------------------------
!
! Step 6.
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
! Step 7.
! Dynamically allocate memory to store the remote morton information.

       call MPI_ALLREDUCE(lnblocks,
     .                    max_no_of_blocks,
     .                    1,
     .                    MPI_INTEGER,
     .                    MPI_MAX,
     .                    MPI_COMM_WORLD,
     .                    ierror)

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

!--------------------------------------------------
!
! Step 8.
! Exchange morton information between processors.

      pe_source   = -1
      isize = 3*max_no_of_blocks*6
      k = 0
      r_mortonbnd = -1

      do i = 1,nprocs
         isrc = i-1
         idest= mype
         itag = isrc*nprocs + idest+1 + tag_offset
                                ! receive to pe=j
          if(commatrix_recv(i).gt.0) then
             k = k+1
             pe_source(k) = isrc+1
             call Mpi_Irecv(r_mortonbnd(1,1,1,k),isize,MPI_INTEGER,
     .            isrc ,itag,MPI_COMM_WORLD,recvrequest(k),ierr)
          endif
       enddo

      ll = 0
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


#ifdef DEBUGY
      write(*,*) 'pe ',mype,' no_of_mortonbnds_received ',
     .          no_of_mortonbnds_received
      write(*,*) 'pe ',mype,' r_mortonbnd(:,1:15,1) ',
     .          r_mortonbnd(:,1:15,1)
#endif /* DEBUG  */
        
!--------------------------------------------------
!
! Step 9.
! Loop over this processor^s list of required neighbor blocks,
! identifying whether they exist from the morton information received
! in step 8.

        do i = 1,no_of_remote_neighs

! first test whether this block can be found locally
          rem_block = -1
          rem_pe = mype+1
          do j=1,max_no_of_blocks
            if( morton_equal(mortonbnd(:,1,j),neigh_morts(:,1,i)) 
     .                    .and.
     .          mortonbnd(6,2,j).eq.neigh_morts(6,2,i) )
     .          rem_block = j
          enddo
! If not found locally, test to see if it is in any of the remote lists
! which have been received.
          if(rem_block.eq.-1) then

            i_pe = 1
            j_pe = -1
            do while( 
     .        (  morton_greater_than(neigh_morts(:,1,i),
     .                               morton_limits(:,1,2,i_pe))
     .                             .or.
     .          (morton_equal(neigh_morts(:,1,i),
     .                        morton_limits(:,1,2,i_pe)).and.
     .           neigh_morts(6,2,i).gt.morton_limits(6,2,2,i_pe)  )  )
     .           .and. (i_pe.le.nprocs)
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
              if( morton_equal(r_mortonbnd(:,1,j,kk),
     .                         neigh_morts(:,1,i)) .and.
     .            r_mortonbnd(6,2,j,kk).eq.neigh_morts(6,2,i) ) then
                 rem_block = j

              endif
            enddo
            endif
            if(rem_block.eq.-1) rem_pe = -1

          else

          endif

!
! If the child block exists then record the (mort no., level) pair 
! of it^s parent.

          if(rem_block.lt.0) then
            neigh_morts(:,:,i) = -1
          else
            neigh_morts(6,2,i) = neigh_morts(6,2,i)-1
          endif

! At this point neigh_morts has been modified to
! contain a list of (morton no, refinement level) pairs associated
! with the parents of any children which were successfully located
! during the first sweep.

        enddo

!--------------------------------------------------
!
! Step 10.
! Check for any non-existent blocks in the neigh_morts list
! and remove them, compressing the list.

      indx = 0
      jstack = 0

      do i=1,no_of_remote_neighs
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
      no_of_remote_neighs = istack

!--------------------------------------------------

#if N_EDGE_VAR > 0

!--------------------------------------------------
!
! Step 11.
! If edge averaging is required then we need to identify leaf blocks
! which have refinement off their diagonal edges but not in the neighbors
! bounding that diagonal element. In this case we will have to make
! sure that the edge integral on the local block match those of the
! refined diagonal neighbor, if we wish to preserve a div B constraint.
!

      jstack = 0
      do lb = 1,lnblocks
        if(nodetype(lb).eq.1) then

! Compute the morton numbers of the children on diagonal edges

! Get morton numbers for neighbors of any leaf blocks
           mort_neigh  = -1
           call morton_neighbors(xmin,ymin,zmin,xmax,ymax,zmax,
     .                      lperiodicx,lperiodicy,lperiodicz,
     .                      coord(1,lb),bsize(1,lb),ndim,
     .                      lrefine(lb),lrefine_max,mort_neigh)

           l_nodetype = 1
           do k = 2-k3d,2+k3d
           do j = 2-k2d,2+k2d
           do i = 1,3
             if(i.ne.2.or.j.ne.2.or.k.ne.2) then
               lfound = .false.
               if(mort_neigh(6,i,j,k).gt.-1) then
                 morton(:) = mort_neigh(:,i,j,k)
                 level= lrefine(lb)
                 do ij = 1,no_of_remote_neighs
                   if( morton_equal(morton(:),neigh_morts(:,1,ij))
     .                      .and.
     .                 level.eq.neigh_morts(6,2,ij) ) then
                     lfound = .true.
                   endif
                 enddo
               endif
               if(lfound) l_nodetype(i,j,k) = 2
             endif
           enddo
           enddo
           enddo

! Now cycle through this list of neighbors constructing the
! (morton no., refinement level) for their bottom left child.
      do ie = 1,12
        i = 2
        j = 2
        k = 2
        if1 = 2
        jf1 = 2
        kf1 = 2
        if2 = 2
        jf2 = 2
        kf2 = 2
        if(ie.eq.1) then
          i = 1
          j = 1
          if1 = 1
          jf2 = 1
        elseif(ie.eq.2) then
          i = 1
          j = 3
          if1 = 1
          jf2 = 3
        elseif(ie.eq.3) then
          i = 3
          j = 1
          if1 = 3
          jf2 = 1
        elseif(ie.eq.4) then
          i = 3
          j = 3
          if1 = 3
          jf2 = 3
        elseif(ie.eq.5) then
          j = 1
          k = 1
          jf1 = 1
          kf2 = 1
        elseif(ie.eq.6) then
          j = 3
          k = 1
          jf1 = 3
          kf2 = 1
        elseif(ie.eq.7) then
          j = 1
          k = 3
          jf1 = 1
          kf2 = 3
        elseif(ie.eq.8) then
          j = 3
          k = 3
          jf1 = 3
          kf2 = 3
        elseif(ie.eq.9) then
          i = 1
          k = 1
          if1 = 1
          kf2 = 1
        elseif(ie.eq.10) then
          i = 1
          k = 3
          if1 = 1
          kf2 = 3
        elseif(ie.eq.11) then
          i = 3
          k = 1
          if1 = 3
          kf2 = 1
        elseif(ie.eq.12) then
          i = 3
          k = 3
          if1 = 3
          kf2 = 3
        endif

           lfound = .false.
           if(mort_neigh(6,i,j,k).gt.-1) then
             morton(:) = mort_neigh(:,i,j,k)
             level= lrefine(lb)
! find (mort,level) in the list of diagonal edges
             do istack = 1,no_of_remote_neighs
               if( morton_equal(morton(:),neigh_morts(:,1,istack))
     .                  .and.
     .             level.eq.neigh_morts(6,2,istack) ) then
                 lfound = .true.
               endif
             enddo
           endif
           if(lfound) then
! Are the neighbors which bound this edge refined?
             nodetype_1 = l_nodetype(if1,jf1,kf1)
             nodetype_2 = l_nodetype(if2,jf2,kf2)
             if(nodetype_1.eq.1.and.nodetype_2.eq.1) then
               jstack = jstack+1
               edge_mark(6,1,jstack) = ie
               edge_mark(6,2,jstack) = lb
               edge_mark(:,3,jstack) = morton(:)
               edge_mark(6,4,jstack) = level
             endif
           endif

      enddo


        endif
      enddo
      no_of_diagonal_edges = jstack

#endif /* N_EDGE_VAR */


!--------------------------------------------------
!
! Step 12.
! Now that we have identified the edges which will need diagonal correction,
! we no longer need any local entries in the neigh_morts list.

!
! Mark any entries in the neigh_morts list which are on this processor
! for removal from the list.

          do i=1,no_of_remote_neighs

            if(
     .      (  morton_greater_than(neigh_morts(:,1,i),
     .                             morton_limits(:,1,1,mype+1))
     .                             .or.
     .        (morton_equal(neigh_morts(:,1,i),
     .                      morton_limits(:,1,1,mype+1)).and.
     .         neigh_morts(6,2,i).ge.morton_limits(6,2,1,mype+1)  )  )
     .    .and.
     .      (  morton_less_than(neigh_morts(:,1,i),
     .                          morton_limits(:,1,2,mype+1))
     .                             .or.
     .        (morton_equal(neigh_morts(:,1,i),
     .                      morton_limits(:,1,2,mype+1)).and.
     .         neigh_morts(6,2,i).le.morton_limits(6,2,2,mype+1)  )  )
     .      ) then
              neigh_morts(:,:,i) = -1
            endif

          enddo

!--------------------------------------------------
!
! Step 13.
! Check for any non-existent blocks in the neigh_morts list
! and remove them.

      indx = 0
      jstack = 0
      do i=1,no_of_remote_neighs
        if(neigh_morts(6,1,i).gt.-1) then
          jstack = jstack+1
          indx(jstack) = i
        endif
      enddo
      do j=1,jstack
        neigh_morts(:,:,j) = neigh_morts(:,:,indx(j))
      enddo
      if(no_of_remote_neighs.gt.jstack)
     .      neigh_morts(:,:,jstack+1:istack) = -1
      istack = jstack
      no_of_remote_neighs = istack

!--------------------------------------------------
!
! If no such blocks were detected exit this routine.

!     if(no_of_remote_neighs.eq.0) return


!--------------------------------------------------
!
! The next stage is to locate the neighbors of the leaf
! blocks which were identified in the previous section.
! All these blocks will exist but they may be on or off
! processor.

!
! ReInitializations
      no_of_comm_procs = 0
      max_no_to_be_received = 0
      max_no_to_send = 0
      commatrix_recv = 0
      commatrix_send = 0
      pe_source = -1
      pe_destination = -1
!
! Step 14.
! Construct a list of all processors from which the local processor should
! request morton number information.

! non-zero elements of COMMATRIX define which processor pairs need to
! exchange morton number lists.

        do i = 1,no_of_remote_neighs
          i_pe = 1
          j_pe = -1
          do while(
     .       (  morton_greater_than(neigh_morts(:,1,i),
     .                              morton_limits(:,1,2,i_pe))
     .                               .or.
     .         (morton_equal(neigh_morts(:,1,i),
     .                       morton_limits(:,1,2,i_pe)).and.
     .          neigh_morts(6,2,i).gt.morton_limits(6,2,2,i_pe)  )  )
     .          .and. (i_pe.le.nprocs)
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

#ifdef DEBUGY
        write(*,*) 'pe ',mype,' commatrix bef gather2 ',
     .             commatrix(1:nprocs,1:nprocs)
        write(*,*) 'pe ',mype,' no_of_remote_neighs ',
     .              no_of_remote_neighs
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
!
! Step 12.
! provide the complete COMMATRIX to all processors

      call MPI_AlltoAll (commatrix_recv,       1,MPI_INTEGER,
     .                   commatrix_send,       1,MPI_INTEGER,
     .                   MPI_COMM_WORLD,ierror)

#ifdef DEBUGY
        write(*,*) 'pe ',mype,' commatrix ',
     .             commatrix(1:nprocs,1:nprocs)

      if(mype.eq.0) then
         write(*,'(" ")')
         write(*,'(" COMMUNICATION MATRIX2: morton_bnd_fluxcon")')
         write(*,'(" ")')
         do ipe=1,nprocs
         write(*,'(" ",8i3)') (commatrix(i,ipe),i=1,nprocs)
         enddo
         write(*,'(" ")')
      endif

#endif /* DEBUG  */
!--------------------------------------------------
!
! Step 13.
! Compute the maximum amount of morton information which any processor
! is going to receive.

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
! Step 14.
! Dynamically allocate memory to store the remote morton information.

       if(allocated(r_mortonbnd)) deallocate(r_mortonbnd)
       allocate( r_mortonbnd(6,3,max_no_of_blocks,
     .                       max(1,max_no_to_be_received) ) )

!--------------------------------------------------
!
! Step 15.
! Exchange morton information between processors.

! Exchange morton information between processors.

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

      ll = 0
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


#ifdef DEBUGY
      write(*,*) 'pe ',mype,' no_of_mortonbnds_received ',
     .          no_of_mortonbnds_received
      write(*,*) 'pe ',mype,' r_mortonbnd(:,1:15,1) ',
     .          r_mortonbnd(:,1:15,1)
#endif /* DEBUG  */
         
!--------------------------------------------------
!
! Step 16.
! Loop over this processor^s list of required neighbor blocks,
! identifying their remote location from the morton information received
! in step ??.

        do i = 1,no_of_remote_neighs
          i_pe = 1
          j_pe = -1
          do while( 
     .      (  morton_greater_than(neigh_morts(:,1,i),
     .                             morton_limits(:,1,2,i_pe))
     .                             .or.
     .        (morton_equal(neigh_morts(:,1,i),
     .                      morton_limits(:,1,2,i_pe)).and.
     .         neigh_morts(6,2,i).gt.morton_limits(6,2,2,i_pe)  )  )
     .         .and. (i_pe.le.nprocs)
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
            if( morton_equal(r_mortonbnd(:,1,j,kk),
     .                       neigh_morts(:,1,i)) .and.
     .          r_mortonbnd(6,2,j,kk).eq.neigh_morts(6,2,i) )
     .          rem_block = j
          enddo
          endif
          if(rem_block.eq.-1) rem_pe = -1

#ifdef DEBUG 
          write(*,*) 'pe ',mype,' neigh i ',i,' rem_pe ',
     .            rem_pe,' kk ',kk,' rem_block ',rem_block
#endif /* DEBUG  */

            if(rem_block.lt.0.or.rem_pe.eq.mype+1) then
              neigh_morts(:,:,i) = -1
            else
              neigh_morts(:,1,i) = rem_block
              neigh_morts(:,2,i) = rem_pe
            endif

! At this point neigh_morts has been modified to
! contain the addresses of neighbors of local leaf
! blocks, provided the neighbors are parents and are located
! off processor.

        enddo


#if N_EDGE_VAR > 0
!
! Step 16a.
! Repeat step 16 but this time replace edge_mark(3:4,..) with
! the appropriate block address.

        do i = 1,no_of_diagonal_edges

          i_pe = 1
          j_pe = -1

          do while( 
     .      (  morton_greater_than(edge_mark(:,3,i),
     .                             morton_limits(:,1,2,i_pe))
     .                             .or.
     .        (morton_equal(edge_mark(:,3,i),
     .                      morton_limits(:,1,2,i_pe)).and.
     .         edge_mark(6,4,i).gt.morton_limits(6,2,2,i_pe)  )  )
     .         .and. (i_pe.le.nprocs)
     .            )
            i_pe = i_pe + 1
            if (i_pe > nprocs) exit
          enddo
          if(i_pe.le.nprocs) j_pe = i_pe

          rem_block = -1
          rem_pe = j_pe

          if(rem_pe.eq.mype+1) then

            do j=1,max_no_of_blocks
              if( morton_equal(mortonbnd(:,1,j),
     .                         edge_mark(:,3,i)) .and.
     .            mortonbnd(6,2,j).eq.edge_mark(6,4,i) )
     .          rem_block = j
            enddo

          else

            kk = -1
            do k=1,no_of_mortonbnds_received
              if(pe_source(k).eq.rem_pe) kk = k 
            enddo
            if(kk.gt.0) then
            do j=1,max_no_of_blocks
              if( morton_equal(r_mortonbnd(:,1,j,kk),
     .                         edge_mark(:,3,i)) .and.
     .            r_mortonbnd(6,2,j,kk).eq.edge_mark(6,4,i) )
     .            rem_block = j
            enddo
            endif

          endif
          if(rem_block.eq.-1) rem_pe = -1

#ifdef DEBUG 
          write(*,*) 'pe ',mype,' neigh i ',i,' rem_pe ',
     .            rem_pe,' kk ',kk,' rem_block ',rem_block
#endif /* DEBUG  */

            if(rem_block.lt.0) then
              edge_mark(:,3:4,i) = -1
            else
              edge_mark(6,3,i) = rem_block
              edge_mark(6,4,i) = rem_pe-1
            endif

        enddo
#endif /* N_EDGE_VAR */


!--------------------------------------------------
!
! Step 17.
! Check for any non-existent blocks in the neigh_morts list
! and remove them.

      indx = 0
      jstack = 0
      do i=1,no_of_remote_neighs
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
      no_of_remote_neighs = istack




!--------------------------------------------------
! Step 17.1.
! Locate the blocks which remain in this list.



!--------------------------------------------------

! Step 18.
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
! Step 19.
! Repeat Step ??.
! provide the complete COMMATRIX to all processors

      call MPI_AlltoAll (commatrix_recv,       1,MPI_INTEGER,
     .                   commatrix_send,       1,MPI_INTEGER,
     .                   MPI_COMM_WORLD,ierror)
 
#ifdef DEBUGY
        write(*,*) 'pe ',mype,' commatrix ',
     .             commatrix(1:nprocs,1:nprocs)

      if(mype.eq.0) then
         write(*,'(" ")')
         write(*,'(" COMMUNICATION MATRIX3: morton_bnd_fluxcon")')
         write(*,'(" ")')
         do ipe=1,nprocs
         write(*,'(" ",8i3)') (commatrix(i,ipe),i=1,nprocs)
         enddo
         write(*,'(" ")')
      endif
#endif /* DEBUG  */


!--------------------------------------------------

! Step 20.
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
! Step 21.
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
! Step 22.
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
! Step 23.
! evaluate smallest guard block starting index over all pe
! store this into variable strt_buffer.

      last_buffer = maxblocks_alloc

      k = last_buffer
      do i=0,nprocs-1
      k = k - commatrix_recv(i+1)
      enddo
      strt_buffer = k + 1


      if (strt_buffer.le.lnblocks) then
        write(*,*) 
     .  'ERROR in mpi_morton_bnd_fluxcon: ',
     .  'guard block starting index',
     .  strt_buffer,' not larger than lnblocks',lnblocks
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif


!--------------------------------------------------
!
! Step 24.
! Dynamically allocate memory to store the lists of blocks to be
! sent and received.

!      iprocs = sum(commatrix_send)
      iprocs = max(maxval(commatrix_send),maxval(commatrix_recv))
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
     .                          max(1,largest_no_of_blocks),
     .                          max(1,max_no_to_send) ) )
      allocate( to_be_received(3,
     .                         max(1,largest_no_of_blocks),
     .                         max(1,max_no_to_be_received) ) )

!--------------------------------------------------
!
! Step 25.
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

#ifdef DEBUGX
        do jp = 1,no_of_mortonbnds_received
        write(*,*) 'fluxbnd pe ',mype,' jreceive ',jp,
     .     ' to_be_received ',to_be_received(:,:,jp)
        enddo
        write(*,*) 'pe ',mype,' laddress to_be_received ',
     .    laddress(:,strt_buffer:last_buffer)
#endif /* DEBUGX */


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
#ifdef DEBUGX
        do jp = 1,kk
        write(*,*) 'flxbnd pe ',mype,' jsend ',jp,
     .       ' to_be_sent ',to_be_sent(:,:,jp)
        enddo
#endif /* DEBUGX */
!--------------------------------------------------
!
! Step 27.
! Deallocate any memory which was dynamically allocated for local use in this
! routine.

       if(allocated(recvrequest)) deallocate( recvrequest )
       if(allocated(recvstatus)) deallocate( recvstatus )


!--------------------------------------------------

! Mark morton data up to date
       morton_limits_set = .true.


! Store communication info for future use
       call mpi_amr_write_flux_comm(nprocs)



#ifdef DEBUG
      write(*,*) 'pe ',mype,' exiting mpi_morton_bnd_fluxcon'
#endif /* DEBUG */


#ifdef TIMING_MPI
              timer_mpi_morton_bnd_fluxcon =  
     .                          timer_mpi_morton_bnd_fluxcon
     .                          + mpi_wtime() - time1
#endif /* TIMING_MPI */
      return
      end subroutine mpi_morton_bnd_fluxcon
