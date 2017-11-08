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

      subroutine amr_mirror_blks(mype,nprocs,mirror_blks)



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
! Written :     Peter MacNeice         December 2002
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
      use paramesh_mpi_interfaces, only : morton_number

      implicit none
      include 'mpif.h'


      integer, intent(in)    ::  mype,nprocs
      integer, intent(out)   ::  mirror_blks(2,maxblocks)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local variables

      integer :: lb,mb,i,j,k,kk
      integer :: mort(6),istack


!??? Is 1000 the best value ?
      integer :: neigh_morts(6,2,npts_neigh)

      integer :: t1neigh_morts(6,npts_neigh)
      integer :: t2neigh_morts(6,npts_neigh)
      integer :: indx(npts_neigh),itemp

      integer :: i_pe,j_pe,rem_block,rem_pe
      integer :: no_of_remote_neighs

      real    :: xmin,xmax,ymin,ymax,zmin,zmax
      real    :: coordm(3)
      logical :: lswap
      logical :: morton_less_than
      logical :: morton_greater_than
      logical :: morton_equal


      integer :: no_of_comms_to_send,kstack,max_no_to_be_received
      integer :: iprocs,isize,isrc,idest,itag,tag_offset,ll
      integer :: ierr,ierror,ierrorcode
      integer :: allocation_status
      integer :: max_no_of_blocks

      integer,dimension (:),  allocatable :: recvrequest
      integer,dimension (:,:),allocatable :: recvstatus


      common/track/ icomm_indicat
      integer :: icomm_indicat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!
!
! This routine assumes that the grid blocks are ordered by morton
! number and that any blocks with different refinement levels but
! the same morton number are ordered from coarse to fine.


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
      no_of_remote_neighs = 0


!--------------------------------------------------
!
! Construct a list of morton numbers of any off-processor 
! mirror blocks.

      istack = 0

#ifdef DEBUGX
      write(*,*) 'xmin,ymin,zmin,xmax,ymax,zmax ',
     . xmin,ymin,zmin,xmax,ymax,zmax
#endif /* DEBUG */

!-------------
      do lb=1,lnblocks
      if(nodetype(lb).eq.1) then


! Mirror coordinate location offset from (xmin,ymin,zmin)

! equatorial symmetry
      coordm(1) = coord(1,lb)
      coordm(2) = grid_ymax - (coord(2,lb) - grid_ymin)
      coordm(3) = coord(3,lb)


      coordm(1) = coordm(1)  - xmin
      if(ndim.ge.2) coordm(2) = coordm(2)  - ymin
      if(ndim.eq.3) coordm(3) = coordm(3)  - zmin


      
! compute morton no for mirror block
      call morton_number(coordm(1),coordm(2),coordm(3),
     .                   bsize(1,lb),ndim,
     .                   lrefine_max,lrefine(lb),
     .                   mort)
      if( ( morton_less_than(mort(:),
     .      morton_limits(:,1,1,mype+1))  .or.
     .    ( morton_equal(mort(:),
     .      morton_limits(:,1,1,mype+1)) .and.
     .     (lrefine(lb).lt.morton_limits(6,2,1,mype+1)) ) )
     .                    .or.
     .    ( morton_greater_than(mort(:),
     .                          morton_limits(:,1,2,mype+1))  .or.
     .    ( morton_equal(mort(:),
     .                    morton_limits(:,1,2,mype+1)) .and.
     .      (lrefine(lb).gt.morton_limits(6,2,2,mype+1)) ) )
     .  ) then
        istack = istack + 1
        neigh_morts(:,1,istack) = mort(:)
        neigh_morts(:,2,istack) = lrefine(lb)
        mirror_blks(1,lb) = istack
        mirror_blks(2,lb) = -1000
#ifdef DEBUG
        write(*,*) 'pe ',mype,lb,
     .                  ' Mirror blk located off proc '
#endif /* DEBUG */
      else
        do mb = 1,lnblocks
          if( morton_equal(mort(:),mortonbnd(:,1,mb)).and.
     .        lrefine(lb).eq.mortonbnd(6,2,mb) )
     .     then
            mirror_blks(1,lb) = mb
            mirror_blks(2,lb) = mype
#ifdef DEBUG
        write(*,*) 'pe ',mype,lb,
     .       'Mirror blk located on proc at ',mb
#endif /* DEBUG */
          endif
        enddo
      endif
      endif
      enddo

      no_of_remote_neighs = istack
      write(*,*) 'mirror_blk : pe ',mype,
     .          ' no_of_remote_neighs ',istack

!-------------
      if(no_of_remote_neighs.gt.0) then

!
! Sort this local list by morton no and refinement level


!
! Sort the list in order of increasing morton number, and then
! sort each sub-list with the same morton number in order of
! increasing refinement level.

! set indx to enable sort permutation to be stored in case we need
! it again.
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

#ifdef DEBUG
       write(*,*) 'mirror pe ',mype,' after sorting neigh_morts ',
     .         neigh_morts(:,:istack)
#endif /* DEBUG  */

! now scan the list identifying segments with the same morton number
! and sort each segment based on refinement level
! order segments with same morton number in order of increasing
! refinement level
      lswap = .false.
      do while (lswap)
        lswap = .false.
        do i = 1,istack-1
          if(morton_equal(neigh_morts(:,1,i),neigh_morts(:,1,i+1))
     .                        .and.
     .       neigh_morts(6,2,i).gt.neigh_morts(6,2,i+1) ) then
            lswap = .true.
            itemp = neigh_morts(6,2,i)
            neigh_morts(6,2,i) = neigh_morts(6,2,i+1)
            neigh_morts(6,2,i+1) = itemp
#ifdef DEBUGX
            write(*,*) 'pe ',mype,' swapping ',i,i+1
#endif /* DEBUG  */
          endif
        enddo
      enddo                           ! end do while


#ifdef DEBUGX
       write(*,*) 'mirror pe ',mype,' after sorting neigh_morts ',
     .         neigh_morts(:,:istack)
#endif /* DEBUG  */



!----------------------------------------------
!----------------------------------------------

!
! Initializations
      commatrix_send = 0
      commatrix_recv = 0
      pe_source = -1

!--------------------------------------------------
!
! Step 5.
! Construct a list of all processors from which the local processor should
! request morton number information.


! non-zero elements of COMMATRIX define which processor pairs need to 
! exchange morton number lists.

        commatrix_recv = 0
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
          enddo
          if(i_pe.le.nprocs) j_pe = i_pe
!
! If block has been located then update commatrix
          if(j_pe.ne.-1) 
     .      commatrix_recv(j_pe) =  commatrix_recv(j_pe) + 1

        enddo

#ifdef DEBUGY
        write(*,*) 'pe ',mype,' mirror commatrix bef gather ',
     .             commatrix_recv(1:nprocs,1:nprocs)
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
#ifdef DEBUGY
       write(*,*) 'pe ',mype,' no_of_comms_to_send ',
     .           no_of_comms_to_send
#endif /* DEBUG  */

!--------------------------------------------------
!       endif                     ! end of istack if test
!--------------------------------------------------
!
! Step 6.
! provide the complete COMMATRIX to all processors

      call MPI_AlltoAll (commatrix_recv,       1,MPI_INTEGER,
     .                   commatrix_send,       1,MPI_INTEGER,
     .                   MPI_COMM_WORLD,ierror)

#ifdef DEBUGX
        write(*,*) 'pe ',mype,' commatrix_send ',
     .             commatrix_send(1:nprocs)
 
      if(mype.eq.0) then
         write(*,'(" ")')
         write(*,'(" COMMUNICATION MATRIX: mirror")')
         write(*,'(" ")')
         write(*,'(" ",8i3)') (commatrix_send(i),i=1,nprocs)
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
       call MPI_ALLREDUCE(iprocs, 
     .                    max_no_to_be_received,
     .                    1,
     .                    MPI_INTEGER,
     .                    MPI_MAX,
     .                    MPI_COMM_WORLD,
     .                    ierror)

#ifdef DEBUGX
       write(*,*) 'pe ',mype,' max_no_to_be_received ',
     .           max_no_to_be_received
#endif /* DEBUG  */


!--------------------------------------------------
!
! Step 8.
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

      tag_offset = 0

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
#ifdef DEBUGX
      write(*,*) 'pe ',mype,' mortonbnd receive of size ',isize,
     .          ' from ',isrc,' tag ',itag,' posted'
     .          ,' nprocs ',nprocs,' tag_offset ',tag_offset
#endif /* DEBUG  */
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
#ifdef DEBUGX
      write(*,*) 'pe ',mype,' mortonbnd of size ',isize,
     .          ' posting send to ',idest,' tag ',itag
     .          ,' nprocs ',nprocs,' tag_offset ',tag_offset
#endif /* DEBUG  */
             call MPI_Ssend(mortonbnd(1,1,1),isize,MPI_INTEGER,
     .            idest,itag,MPI_COMM_WORLD,ierr)
#ifdef DEBUGX
      write(*,*) 'pe ',mype,' mortonbnd of size ',isize,
     .          ' sent to ',idest,' tag ',itag
#endif /* DEBUG  */
          endif
      enddo

      no_of_mortonbnds_received = k

      tag_offset = (nprocs-1)*nprocs + nprocs + tag_offset

      if(k.gt.0)
     .    call MPI_Waitall(k,recvrequest,recvstatus,
     .                     ierror)


#ifdef DEBUGX
      write(*,*) 'pe ',mype,' no_of_mortonbnds_received ',
     .          no_of_mortonbnds_received
      write(*,*) 'pe ',mype,' r_mortonbnd(:,1:15,1) ',
     .          r_mortonbnd(:,1:15,1)
#endif /* DEBUG  */
        
!--------------------------------------------------




!----------------------------------------------
!----------------------------------------------


        
!--------------------------------------------------
!
! Loop over this processor^s list of required neighbor blocks,
! identifying their remote location from the morton information received

#ifdef DEBUGX 
          write(*,*) 'pe ',mype,' morton_limits ',morton_limits
#endif /* DEBUG  */
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
          enddo
          if(i_pe.le.nprocs) j_pe = i_pe

          rem_block = -1
          rem_pe = j_pe

          if(mype.eq.0.and.i.eq.1) write(*,*) ' neigh_morts ',
     .       neigh_morts(:,:,i),' rem_pe ',rem_pe,
     .      ' pe_source ',pe_source

          kk = -1
          do k=1,no_of_mortonbnds_received
            if(pe_source(k).eq.rem_pe) kk = k 
          enddo
          if(mype.eq.0.and.i.eq.1) write(*,*) ' neigh_morts ',
     .       neigh_morts(:,:,i),' kk ',kk
#ifdef DEBUG
          write(*,*) 'pe ',mype,' neigh i ',i,' i_pe ',
     .            i_pe,' neigh_morts ',neigh_morts(:,i),
     .             ' kk ',kk,' no_of_mortonbnds_received ',
     .              no_of_mortonbnds_received,
     .             ' rem_pe ',rem_pe,
     .             ' icomm_indicat ',icomm_indicat
#endif /* DEBUG  */

          if(kk.gt.0) then
          do j=1,maxblocks
            if( morton_equal(r_mortonbnd(:,1,j,kk),neigh_morts(:,1,i)) 
     .                    .and.
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

          neigh_morts(6,1,i) = rem_block
          neigh_morts(6,2,i) = rem_pe-1

#ifdef DEBUG 
          write(*,*) 'pe ',mype,' neigh i ',i,' address ',
     .            neigh_morts(:,i)
#endif /* DEBUG  */
        enddo

!--------------------------------------------------

      do lb = 1,lnblocks
      if( (nodetype(lb).eq.1) .and.
     .    (mirror_blks(2,lb).eq.-1000) ) then
        istack = mirror_blks(1,lb)
        mirror_blks(1:2,lb) = neigh_morts(6,1:2,istack)
      endif
      enddo

#ifdef DEBUGX
      do lb = 1,lnblocks
      if(nodetype(lb).eq.1) then
          write(*,*) 'pe ',mype,' lb ',lb,' mirror_blk ',
     .            mirror_blks(:,lb)
      endif
      enddo

#endif /* DEBUG  */


      endif

      return
      end subroutine amr_mirror_blks
