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

      subroutine mpi_mort_comm_for_surrblks
     .                   (mype,nprocs,tag_offset)



!
! DESIGN ISSUES :
!  
!  Some neighbors of parents are requested unnecessarily because we
!  cannot verify that the corresponding child neighbors actually
!  exist until after we have received morton lists.
!  Can we improve list of requested blocks by identifying neighbors
!  of parents which may not be needed?
!
!
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
      use mpi_morton

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
      integer :: morton(6),level,jstack
      integer ::  lbfirst,lblast
      integer :: mort_neigh(6,3,3,3)
      integer :: pmort_neigh(6,3,3,3)
      integer :: neigh_morts(6,3,npts_neigh)
      integer :: t1neigh_morts(6,npts_neigh)
      integer :: t2neigh_morts(6,npts_neigh)
      integer :: istart,iend,indx(npts_neigh)
      integer :: i_pe,j_pe,rem_block,rem_pe
      integer :: no_of_comm_procs
      integer :: ierrorcode,ierr,allocation_status,ierror
      integer :: interp_max_orderf,interp_max_ordere
      integer :: interp_max_order ,interp_max_ordern,interx
      integer :: no_of_remote_neighs
      integer :: max_no_to_be_received
      integer :: max_no_of_blocks
      integer :: no_of_comms_to_send
      integer :: istack, ioff, joff, koff, k, itag, ll, kk
      integer :: itemp, kstack, iprocs, isize, isrc, idest
      integer,dimension (:),  allocatable :: recvrequest
      integer,dimension (:,:),allocatable :: recvstatus
      integer, parameter :: nguarda = max(nguard,nguard_work)
#ifdef ADVANCE_ALL_LEVELS
      integer :: ii, jj
#endif

      logical :: lremote,lswap,lfound
      logical :: is_remote,is_found
      logical :: morton_greater_than
      logical :: morton_equal
#ifdef ADVANCE_ALL_LEVELS
      logical :: morton_less_than
#endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       lbfirst = 1
       lblast  = lnblocks
!
!
! This routine assumes that the grid blocks are ordered by morton
! number and that any blocks with different refinement levels but
! the same morton number are ordered from coarse to fine.

! mark morton data out of date.
!      morton_limits_set = .false.

! Find highest order interpolant used in prolongation, if using face, edge
! or corner data. If Muscl is used then this will be interpreted as order 1.
       interp_max_orderf = 0
       interp_max_ordere = 0
       interp_max_ordern = 0
       do i = 1,nfacevar
          interx = max(interp_mask_facex(i),
     .                 interp_mask_facey(i),
     .                 interp_mask_facez(i))
          if(interx.eq.20) interx = 1
          interp_max_orderf = max(interp_max_orderf,interx)
       enddo
       do i = 1,nvaredge
          interx = interp_mask_ec(i)
          if(interx.eq.20) interx = 1
          interp_max_ordere = max(interp_max_ordere,interx)
       enddo
       do i = 1,nvarcorn
          interx = interp_mask_nc(i)
          if(interx.eq.20) interx = 1
          interp_max_ordern = max(interp_max_ordern,interx)
       enddo
       interp_max_order = max(interp_max_orderf,interp_max_ordere,
     .                        interp_max_ordern)


!--------------------------------------------------
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
! Step 3.
! Construct a list of potential neighbors of all blocks on this
! processor, and potential neighbors of their parents.
! Exclude any which are on this processor.

      istack = 0

#ifdef DEBUG
      write(*,*) 'xmin,ymin,zmin,xmax,ymax,zmax ',
     . xmin,ymin,zmin,xmax,ymax,zmax
#endif /* DEBUG */


!     do lb=1,lnblocks
      do lb=lbfirst,lblast

!-------------

! First get the possible neighbors of the current block
      mort_neigh = -1
      pmort_neigh = -1
      call morton_neighbors(xmin,ymin,zmin,xmax,ymax,zmax,
     .                      lperiodicx,lperiodicy,lperiodicz,
     .                      coord(:,lb),bsize(:,lb),ndim,
     .                      lrefine(lb),lrefine_max,mort_neigh)


! Now get the possible neighbors of the current block^s parent
      if(parent(1,lb).gt.0) then
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
     .                        pcoord(:),pbsize(:),ndim,
     .                        lrefine(lb)-1,lrefine_max,pmort_neigh)

      endif                      ! end of parent if test


!-------------

! If parent is a remote block then puts its address on the list of
! remote blocks which are required.
      if(parent(1,lb).gt.0.and.parent(2,lb).ne.mype) then
            istack = istack+1
#ifdef DEBUG
            if(istack.gt.npts_neigh) then
              write(*,*) 'morton_bnd : ',
     .                   'istack exceeds npts_neigh : ',
     .                   'possible solution - increase npts_neigh'
              call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
            endif
#endif /* DEBUG */
            neigh_morts(:,1,istack) = pmort_neigh(:,2,2,2)
            neigh_morts(6,2,istack) = lrefine(lb)-1
            neigh_morts(6,3,istack) = 14      ! marks as a full block request
#ifdef DEBUG
      write(*,*) 'parent pmort_neigh(2,2,2) ',pmort_neigh(2,2,2),
     .       ' istack ',istack,' pe ',mype,' parent ',parent(:,lb)
     .  ,' of block ',lb
#endif /* DEBUG */
      endif                     ! end of parent if test

!-------------

! Now start to build the array neigh_morts which is a list of possible
! remote blocks which will be needed.

! First add any neighbors of the current block, eliminating
! any local blocks on the list


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
! compute message type - note this index is computed to reflect the part
! of the remote block to be acquired, not the part of the local blocks
! guardcells which will be filled.
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
          endif                         ! if(lremote)
        endif                           ! if(mort_neigh(i,j,k).gt.-1)

!#ifdef NEWCODE
#ifdef ADVANCE_ALL_LEVELS
!#endif /* NEWCODE */


! Now consider required neighbors of the current blocks parent
        if(parent(1,lb).gt.0) then

! If interpolation order or number of guard cells is small enough, we can limit
! this list to neighbors of the parent relevant to the child block under consideration.
! Otherwise we consider all parent neighbors.
          ioff = mod(which_child(lb)-1,2) 
          joff = mod((which_child(lb)-1)/2,2)
          koff = mod((which_child(lb)-1)/4,2)
          ii = i
          jj = j
          kk = k
          if(interp_max_order.lt.nmax_lays) then
            if(ioff.eq.0.and.i.eq.3) ii = 2
            if(ioff.eq.1.and.i.eq.1) ii = 2
            if(joff.eq.0.and.j.eq.3) jj = 2
            if(joff.eq.1.and.j.eq.1) jj = 2
            if(koff.eq.0.and.k.eq.3) kk = 2
            if(koff.eq.1.and.k.eq.1) kk = 2
          endif
          lremote = .false.
          if (pmort_neigh(6,i,j,k) > -1) then

          if( morton_less_than(pmort_neigh(:,ii,jj,kk),
     .                         morton_limits(:,1,1,mype+1)) )
     .         lremote = .true.
          if( morton_greater_than(pmort_neigh(:,ii,jj,kk),
     .                            morton_limits(:,1,2,mype+1)) )
     .         lremote = .true.
          if( (morton_equal(pmort_neigh(:,ii,jj,kk),
     .                      morton_limits(:,1,1,mype+1)))
     .             .and.
     .        (lrefine(lb)-1.lt.morton_limits(6,2,1,mype+1)) )
     .         lremote = .true.
          if( (morton_equal(pmort_neigh(:,ii,jj,kk),
     .                      morton_limits(:,1,2,mype+1)))
     .             .and.
     .        (lrefine(lb)-1.gt.morton_limits(6,2,2,mype+1)) )
     .         lremote = .true.


          if(lremote) then 
            istack = istack+1
            neigh_morts(:,1,istack) = pmort_neigh(:,ii,jj,kk)
            neigh_morts(6,2,istack) = lrefine(lb)-1
            neigh_morts(6,3,istack) = (4-ii)+((4-jj)-1)*3+((4-kk)-1)*9
            if(nguarda.gt.nmax_lays) neigh_morts(6,3,istack) = 14
#ifdef DEBUG 
            write(*,*) 'pe ',mype,' blk ',lb,' iijjkk ',ii,jj,kk,
     .             ' neigh_morts ',
     .              neigh_morts(:,istack),' istack ',istack
#endif /* DEBUG  */
          else
#ifdef DEBUG 
            write(*,*) 'pe ',mype,' blk ',lb,' iijjkk ',
     .         ii,jj,kk,' parent neigh local'
#endif /* DEBUG  */
          endif                      ! if(lremote)
          endif                      ! if (pmort_neigh > -1)
        endif                        ! if(parent(1,lb).gt.0)

!#ifdef NEWCODE
#endif /* ADVANCE_ALL_LEVELS */
!#endif /* NEWCODE */


        endif                   ! if(i.ne.2.or.j.ne.2.or.k.ne.2)
      enddo
      enddo
      enddo



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

#ifdef DEBUG 
       write(*,*) 'pe ',mype,' after sorting neigh_morts ',
     .         neigh_morts
#endif /* DEBUG  */

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
     .       neigh_morts(6,2,i) > neigh_morts(6,2,i+1) ) then
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

        do while(morton_equal(neigh_morts(:,1,i),morton(:)) 
     .                             .and.
     .           neigh_morts(6,2,i) == level   .and.
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



#ifdef DEBUG
      write(*,*) 'mpi_morton_bnd : proc ',mype,' istack ',istack
      do i=1,istack
       write(*,*) 'pe ',mype,i,' neigh_morts ',neigh_morts(:,i)
      enddo
#endif /* DEBUG */

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
        if(neigh_morts(6,1,i) > -1) then
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
           level     = neigh_morts(6,2,i)
#ifdef DEBUG
         write(*,*) 'compress: pe ',mype,' jstack ',jstack,
     .      ' morton ',morton,' level ',level
#endif /* DEBUG  */
         endif
       enddo

       no_of_remote_neighs = jstack

!      call report_comm_saving(mype,jstack,neigh_morts)


#ifdef DEBUG
      do i=1,jstack
      write(*,*) 'pe ',mype,' jstack ',i,' neigh_morts ',
     .            neigh_morts(:,i)
      enddo
!      pause
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
         write(*,'(" COMMUNICATION MATRIX1: morton_bnd")')
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
! Step 8.

       call MPI_ALLREDUCE(lnblocks, 
     .                    max_no_of_blocks,
     .                    1,
     .                    MPI_INTEGER,
     .                    MPI_MAX,
     .                    MPI_COMM_WORLD,
     .                    ierror)

! Dynamically allocate memory to store the remote morton information.

       if(allocated(r_mortonbnd)) deallocate(r_mortonbnd)
       allocate( r_mortonbnd(6,3,max_no_of_blocks,
     .           max(1,max_no_to_be_received) ),
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

      pe_source   = -1
      isize = 3*max_no_of_blocks*6
      k = 0
      r_mortonbnd = -1

      do i = 1,nprocs
         isrc = i-1
         idest= mype
         itag = isrc*nprocs + idest+1 + tag_offset

                                ! receive to pe=j
         if((commatrix_recv(i).gt.0)) then
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
     .            idest,itag,MPI_COMM_WORLD,ierr)
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

! neigh_morts(1:2,no_of_remote_neighs) is now being used to store 
! the remote addresses of the required neighbors.
! Here proc nos. run from 1 to nprocs.

          neigh_morts(:,1,i) = rem_block
          neigh_morts(:,2,i) = rem_pe

#ifdef DEBUG 
          write(*,*) 'pe ',mype,' neigh i ',i,' address ',
     .            neigh_morts(:,i)
#endif /* DEBUG  */
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
#ifdef DEBUG 
          write(*,*) 'pe ',mype,' stack entry ',neigh_morts(:,i),
     .     ' does exists - not to be removed '
#endif /* DEBUG  */
          jstack = jstack+1
          indx(jstack) = i
        endif
      enddo
      do j=1,jstack
        neigh_morts(6,:,j) = neigh_morts(6,:,indx(j))
#ifdef DEBUG 
        write(*,*) 'pe ',mype,' remaining stack entry ',j,
     . ' neigh_morts(:,j) ',neigh_morts(:,j)
#endif /* DEBUG  */
      enddo
      if(no_of_remote_neighs.gt.jstack)
     .      neigh_morts(6,:,jstack+1:no_of_remote_neighs) = -1
#ifdef DEBUG 
      write(*,*) 'pe ',mype,' removed stack items ',jstack+1,
     .       ' to ',no_of_remote_neighs
#endif /* DEBUG  */
      istack = jstack
      no_of_remote_neighs = istack


!--------------------------------------------------
! Step 12.
! Reconstruct commatrix.


! non-zero elements of COMMATRIX define which processor pairs need to 
! exchange morton number lists. 
        commatrix_send = 0
        commatrix_recv = 0
        do i = 1,no_of_remote_neighs
          i_pe = neigh_morts(6,2,i)
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
         write(*,'(" COMMUNICATION MATRIX2: morton_bnd")')
         write(*,'(" ")')
         do ipe=1,nprocs
         write(*,'(" ",8i3)') (commatrix(i,ipe),i=1,nprocs)
         enddo
         write(*,'(" ")')
      endif
#endif /* DEBUG  */


!--------------------------------------------------
!
! Step 20.
! Deallocate any memory which was dynamically allocated for local use in this
! routine.

       if(allocated(recvrequest)) deallocate( recvrequest )
       if(allocated(recvstatus)) deallocate( recvstatus )


!--------------------------------------------------

      return
      end subroutine mpi_mort_comm_for_surrblks



