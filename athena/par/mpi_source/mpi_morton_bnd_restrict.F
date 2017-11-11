!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!#define ADD_CHILD_NEIGHS
#include "paramesh_preprocessor.fh"
#define TESTING

! NOTE - THIS ROUTINE HAS NOT BEEN DEBUGGED FOR VAR_DT
!        SINCE IT WAS MODIFIED FOR MPI    - July 14, 2000.


!#define DEBUG

      subroutine mpi_morton_bnd_restrict(mype,nprocs,lfulltree,lec,lnc,
     .                                   tag_offset)



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
! Written :     Peter MacNeice                  July 2000
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
      use paramesh_mpi_interfaces, only : mpi_amr_write_restrict_comm,
     .                                    morton_number

      implicit none

      include 'mpif.h'

      integer, intent(in)    ::  mype,nprocs
      logical, intent(in)    ::  lfulltree,lec,lnc
      integer, intent(inout) ::  tag_offset

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local variables

      real    :: dx,dy,dz, x0,y0,z0
      real    :: pbsize(3),pcoord(3)
      real    :: xmin, ymin, zmin
      real    :: xmax, ymax, zmax

      integer :: lb,i,j
      integer :: level,jstack,iprocs
      integer :: type_max
      integer :: morton(6)
      integer :: mort_neigh(6,3,3,3)
      integer :: pmort_neigh(6,3,3,3)
      integer, parameter :: nguarda = max(nguard,nguard_work)
      integer :: neigh_morts(6,3,npts_neigh)
      integer :: t1neigh_morts(6,npts_neigh)
      integer :: t2neigh_morts(6,npts_neigh)
      integer :: indx(npts_neigh)
      integer :: i_pe,j_pe,rem_block,rem_pe
      integer :: no_of_comm_procs
      integer :: ierrorcode,ierr,allocation_status,ierror
      integer :: mort_child(6,mchild),mort(6)
      integer :: no_of_remote_neighs
      integer :: max_no_to_be_received
      integer :: no_of_comms_to_send
      integer :: max_no_of_blocks
      integer :: no_of_comms_to_receive
      integer :: istack, ich, k, ioff, joff, koff, itemp, ll, isize
      integer :: isrc, idest, itag, kk, jj, jp, ip, ii, kstack
      integer,dimension (:),  allocatable :: recvrequest
      integer,dimension (:,:),allocatable :: recvstatus
#ifdef DEBUG_MORTON
      integer :: inxt_bit,nbitshft,nbits
#endif
#ifdef ADD_CHILD_NEIGHS
      integer :: cmort_neigh(6,4,4,4)
      integer :: i0,j0,k0,lrefine_ch
      real    :: coord_ch(3),bsize_ch(3)
      integer :: lrefine_fl
#endif /* ADD_CHILD_NEIGHS */

      logical :: lremote,lswap,lfound
      logical :: is_remote
      logical :: morton_greater_than
      logical :: morton_less_than
      logical :: morton_equal


#ifdef TIMING_MPI
      double precision :: time1
#endif /* TIMING_MPI */
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef TIMING_MPI
      time1 = mpi_wtime()
#endif /* TIMING_MPI */
#ifdef DEBUG
      write(*,*) 'pe ',mype,' entering mpi_morton_bnd_restrict'
#endif /* DEBUG */

! This routine assumes that the grid blocks are ordered by morton
! number and that any blocks with different refinement levels but
! the same morton number are ordered from coarse to fine.


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
      commatrix_send = 0
      commatrix_recv = 0
      pe_source = -1
      pe_destination = -1
!     neigh_morts = -1
      no_of_comms_to_send = 0

! set message types to whole block. May want to tailor this later
!      neigh_morts(:,3,:) = 14


! Select the nodetypes of parent blocks which will be considered for
! receiving restricted data.
      if(lfulltree) then
        type_max = 3
      else
        type_max = 2
      endif
!--------------------------------------------------
!
! Step 3.
! Compute morton numbers of the children of any parents of leaf blocks.
! We could use the child(:,:) addresses directly, but this would make it
! difficult to remove redundancies in the list of required remote
! blocks when these may be needed for odd nxb or for higher order 
! interpolation.

      istack = 0

#ifdef DEBUG
      write(*,*) 'xmin,ymin,zmin,xmax,ymax,zmax ',
     . xmin,ymin,zmin,xmax,ymax,zmax
#endif /* DEBUG */

      do lb=1,lnblocks

! Loop over parent blocks only, identifying their children
! and putting them on the list of needed remote blocks
      if(nodetype(lb) == 2) then

        dx = bsize(1,lb)
        dy = bsize(2,lb)
        dz = bsize(3,lb)
! set the block center
        x0 = coord(1,lb) - xmin
        y0 = coord(2,lb) - ymin
        z0 = coord(3,lb) - zmin

! compute the morton number for block lb
        call morton_number(x0,y0,z0,bsize(:,lb),ndim,
     .                     lrefine_max,lrefine(lb),mort)
!
! shift the morton number by the refinement level of the block
! and store the morton number and refinement level
#ifdef DEBUG_MORTON
            nbitshft = ndim*(lrefine_max-lrefine(lb))
            nbits = bit_size(mort)
            if(nbitshft.gt.0) then
            inxt_bit = ibits(mort,nbits-1-(nbitshft-1),
     .                       nbitshft-1)
            if(inxt_bit.gt.0) then
                write(*,*) 'Error 2 in mpi_morton_bnd_restrict'
                call amr_abort()
            endif
            endif
#endif


        call morton_children(mort,lrefine(lb),lrefine_max,ndim,
     .                           mort_child,nchild)

        do ich = 1,nchild

          lremote = is_remote(mort_child(:,ich),lrefine(lb)+1,mype)

          if(lremote) then
            istack = istack+1
#ifdef DEBUG
            if(istack.gt.npts_neigh) then
              write(*,*) 'morton_bnd_restrict : ',
     .                   'istack exceeds npts_neigh : ',
     .                   'possible solution - increase npts_neigh'
              call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
            endif
#endif /* DEBUG */
            neigh_morts(:,1,istack) = mort_child(:,ich)
            neigh_morts(6,2,istack) = lrefine(lb)+1
!!!
            neigh_morts(6,3,istack) = 14
!!!
          endif
        enddo

#ifndef VAR_DT
      if(lec.or.lnc) then
#endif

!
! add neighbors (only those in the neigh arrays) of current parent
! to list if using variable timestepping.
        mort_neigh = -1
        call morton_neighbors(xmin,ymin,zmin,xmax,ymax,zmax,
     .                          lperiodicx,lperiodicy,lperiodicz,
     .                          coord(:,lb),bsize(:,lb),ndim,
     .                          lrefine(lb),lrefine_max,mort_neigh)

! if neighbor block exists at this refinement level
        do k = 1,1+2*k3d
        do j = 1,1+2*k2d
        do i = 1,3
          if(surr_blks(1,i,j,k,lb).gt.0.and.
     .       surr_blks(2,i,j,k,lb).ne.mype) then
            istack = istack+1
            neigh_morts(:,1,istack) = mort_neigh(:,i,j,k)
            neigh_morts(6,2,istack) = lrefine(lb)
!!!
!            jj = j
!            kk = k
!            if (j == 1 .and. k2d == 0) jj = 2
!            if (k == 1 .and. k3d == 0) kk = 2
!            neigh_morts(6,3,istack) = (4-i)+((4-jj)-1)*3+((4-kk)-1)*9
!            if(nguarda.gt.nmax_lays) neigh_morts(6,3,istack) = 14
            neigh_morts(6,3,istack) = 14
!!!
          endif
        enddo
        enddo
        enddo

#ifndef VAR_DT
      endif            ! lec.or.lnc iftest
#endif

      endif
      enddo


      if(lec.or.lnc) then

!
! If there is an odd number of data points in each block then identify
! neighbors of the leaf node children of local parents

      do lb=1,lnblocks

! Loop over parent blocks only, identifying their neighbors
! and their children and putting them on the list of needed remote blocks
      if(nodetype(lb) == 2) then

!-------------
! First get this parents neighbors
      mort_neigh = -1
      call morton_neighbors(xmin,ymin,zmin,xmax,ymax,zmax,
     .                      lperiodicx,lperiodicy,lperiodicz,
     .                      coord(:,lb),bsize(:,lb),ndim,
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
! compute message type - note this index is computed to reflect the part
! of the remote block to be acquired, not the part of the local blocks
! guardcells which will be filled.
            neigh_morts(6,3,istack) = 14
!!!
!            neigh_morts(6,3,istack) = (4-i)+((4-j)-1)*3+((4-k)-1)*9
!            if(nguarda.gt.nmax_lays) neigh_morts(6,3,istack) = 14
!!!
          else
          endif                         ! if(lremote)
        endif                           ! if(mort_neigh(i,j,k).gt.-1)

        endif                   ! if(i.ne.2.or.j.ne.2.or.k.ne.2)
      enddo
      enddo
      enddo

!-------------
! Now get this parents children
        do ich = 1, nchild

! Now get the possible neighbors of the current block^s child ich
          pbsize(:) = bsize(:,lb)/2.          ! size of child block
          ioff = mod(ich-1,2)        
          joff = mod((ich-1)/2,2)
          koff = mod((ich-1)/4,2)
          if(ioff.eq.0) then
            pcoord(1) = 0.5*(coord(1,lb)+bnd_box(1,1,lb))
          else
            pcoord(1) = 0.5*(coord(1,lb)+bnd_box(2,1,lb))
          endif
          if(joff.eq.0) then
            pcoord(2) = 0.5*(coord(2,lb)+bnd_box(1,2,lb))
          else
            pcoord(2) = 0.5*(coord(2,lb)+bnd_box(2,2,lb))
          endif
          if(ndim.lt.2) pcoord(2) = coord(2,lb)
          if(koff.eq.0) then
            pcoord(3) = 0.5*(coord(3,lb)+bnd_box(1,3,lb))
          else
            pcoord(3) = 0.5*(coord(3,lb)+bnd_box(2,3,lb))
          endif
          if(ndim.lt.3) pcoord(3) = coord(3,lb)

          pmort_neigh = -1
          call morton_neighbors(xmin,ymin,zmin,xmax,ymax,zmax,
     .                          lperiodicx,lperiodicy,lperiodicz,
     .                          pcoord(:),pbsize(:),ndim,
     .                          lrefine(lb)+1,lrefine_max,pmort_neigh)

!-------------

      do k = 2-k3d,2+k3d
      do j = 2-k2d,2+k2d
      do i = 1,3
        if(i.ne.2.or.j.ne.2.or.k.ne.2) then

! if neighbor block exists at this refinement level

        if(pmort_neigh(6,i,j,k).gt.-1) then
          
          lremote = is_remote(pmort_neigh(:,i,j,k),lrefine(lb),mype)

          if(lremote) then 
            istack = istack+1
            neigh_morts(:,1,istack) = pmort_neigh(:,i,j,k)
            neigh_morts(6,2,istack) = lrefine(lb)+1
            neigh_morts(6,3,istack) = 14
!!!
!            neigh_morts(6,3,istack) = (4-i)+((4-j)-1)*3+((4-k)-1)*9
!            if(nguarda.gt.nmax_lays) neigh_morts(6,3,istack) = 14
!!!
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


      enddo                ! end of loop over children

      endif                ! end of nodetype if test



#ifdef ADD_CHILD_NEIGHS

! Now get the possible neighbors of the current blocks children,
! if the current block is a leaf block, and even if those children do not exist.
! This is needed to support amr_1blk_guardcell_f_to_c called during guardcell
! filling for restriction of edge and node centered data

      if(nodetype(lb).eq.1) then

      cmort_neigh = -1
      call morton_child_neighbors(xmin,ymin,zmin,xmax,ymax,zmax,
     .                      lperiodicx,lperiodicy,lperiodicz,
     .                      coord(:,lb),bsize(:,lb),ndim,
     .                      lrefine(lb),lrefine_max,cmort_neigh)

! add any neighbors of the children of the current block which are off proc
! to the neigh_morts array

      do k = 1,1+3*k3d
      do j = 1,1+3*k2d
      do i = 1,4

        lrefine_fl = lrefine(lb)+1

! if neighbor block exists at this refinement level

        if(cmort_neigh(6,i,j,k).gt.-1) then
          lremote = .false.

          lremote = is_remote(cmort_neigh(:,i,j,k),lrefine_fl,mype)

          if(lremote) then
            istack = istack+1
            neigh_morts(:,1,istack) = cmort_neigh(:,i,j,k)
            neigh_morts(6,2,istack) = lrefine_fl
! compute message type - note this index is computed to reflect the part
! of the remote block to be acquired, not the part of the local blocks
! guardcells which will be filled.
            i0 = i/2+1
            j0 = 2 + k2d*( j/2+1 -2)
            k0 = 2 + k3d*( k/2+1 -2)
!            neigh_morts(6,3,istack) = (4-i0)+((4-j0)-1)*3+((4-k0)-1)*9 
!     .                              + 27
            neigh_morts(6,3,istack) = 14
            if(nguarda.gt.nmax_lays) neigh_morts(:,3,istack) = 14
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

      enddo
      enddo
      enddo

      endif




! Now treat the case where a local parent block, needs its children
! with guardcells filled, but one or more of the children have more
! refined neighbors which are off processor. If we do not add these
! morton numbers to the list, the child will have bad guardcell data
! after amr_1blk_guardcell_f_to_c is called, and this will give the
! wrong result for node data during restriction.

      if(nodetype(lb).eq.2) then

      do ich = 1,nchild

       bsize_ch(1) = bsize(1,lb)*.5
       bsize_ch(2) = bsize(2,lb)/real(1+k2d)
       bsize_ch(3) = bsize(3,lb)/real(1+k3d)
       coord_ch(1) = coord(1,lb)+ bsize_ch(1)*(.5-real(mod(ich,2)))
       coord_ch(2) = coord(2,lb)+ 
     .          bsize_ch(2)*(.5-real(mod((ich+1)/2,2)))*real(k2d)
       coord_ch(3) = coord(3,lb)+ 
     .          bsize_ch(3)*(.5-real(mod((ich-1)/4+1,2)))*real(k3d)
       lrefine_ch = lrefine(lb)+1

      cmort_neigh = -1
      call morton_child_neighbors(xmin,ymin,zmin,xmax,ymax,zmax,
     .                      lperiodicx,lperiodicy,lperiodicz,
     .                      coord_ch(:),bsize_ch(:),ndim,
     .                      lrefine_ch,lrefine_max,cmort_neigh)

! add any neighbors of the children of the current block which are off proc
! to the neigh_morts array

      do k = 1,1+3*k3d
      do j = 1,1+3*k2d
      do i = 1,4

! if neighbor block exists at this refinement level

        if(cmort_neigh(6,i,j,k).gt.-1) then

          lremote = is_remote(cmort_neigh(:,i,j,k),lrefine(lb),mype)

          if(lremote) then
            istack = istack+1
            neigh_morts(:,1,istack) = cmort_neigh(:,i,j,k)
            neigh_morts(6,2,istack) = lrefine(lb)+2
! compute message type - note this index is computed to reflect the part
! of the remote block to be acquired, not the part of the local blocks
! guardcells which will be filled.
            i0 = i/2+1
            j0 = 2 + k2d*( j/2+1 -2)
            k0 = 2 + k3d*( k/2+1 -2)
!           neigh_morts(6,3,istack) = (4-i0)+((4-j0)-1)*3+((4-k0)-1)*9 
!     .                              + 27
            neigh_morts(6,3,istack) = 14
            if(nguarda.gt.nmax_lays) neigh_morts(:,3,istack) = 14
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

      enddo
      enddo
      enddo

      enddo                  ! end of loop over children

      endif

#endif /* ADD_CHILD_NEIGHS */

      enddo                ! end of loop over blocks

      endif               ! if (lec 


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
          if( (morton_equal(neigh_morts(:,1,i),neigh_morts(:,1,i+1)))
     .                             .and.
     .         neigh_morts(6,2,i).gt.neigh_morts(6,2,i+1) ) then
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

            if(  
     .    ( (morton_greater_than(neigh_morts(:,1,i),
     .                           morton_limits(:,1,1,j_pe)))
     .                               .or.
     .      ( (morton_equal(neigh_morts(:,1,i),
     .                      morton_limits(:,1,1,j_pe))).and.
     .        (neigh_morts(6,2,i).ge.morton_limits(6,2,1,j_pe))    )
     .    ) .and.
     .    ( (morton_less_than(neigh_morts(:,1,i),
     .                        morton_limits(:,1,2,j_pe)))
     .                               .or.
     .      ( (morton_equal(neigh_morts(:,1,i),
     .                      morton_limits(:,1,2,j_pe))).and.
     .        neigh_morts(6,2,i).le.morton_limits(6,2,2,j_pe)    )
     .    ) 
     .        ) lfound = .true.

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
       morton(:) = -1
       level  = -1
       do i = 1,istack
         if(any((neigh_morts(:,1,i).ne.morton(:))) .or.
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
     .       (  morton_greater_than(neigh_morts(:,1,i),
     .                              morton_limits(:,1,2,i_pe))
     .                               .or.
     .         (morton_equal(neigh_morts(:,1,i),
     .                       morton_limits(:,1,2,i_pe))
     .                              .and.
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
         write(*,'(" COMMUNICATION MATRIX1: morton_restrict")')
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
            if( (morton_equal(r_mortonbnd(:,1,j,kk),
     .                        neigh_morts(:,1,i))) 
     .                      .and.
     .          r_mortonbnd(6,2,j,kk).eq.neigh_morts(6,2,i) )
     .          rem_block = j
          enddo
          endif
          if(rem_block.eq.-1) rem_pe = -1

#ifdef DEBUG 
          write(*,*) 'pe ',mype,' neigh i ',i,' rem_pe ',
     .            rem_pe,' kk ',kk,' rem_block ',rem_block
     .            ,' no_of_remote_neighs ',no_of_remote_neighs
#endif /* DEBUG  */

! neigh_morts(1:2,no_of_remote_neighs) is now being used to store 
! the remote addresses of the required neighbors.
! Here proc nos. run from 1 to nprocs.
          neigh_morts(6,1,i) = rem_block
          neigh_morts(6,2,i) = rem_pe

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
     . ' neigh_morts(6,j) ',neigh_morts(6,j)
#endif /* DEBUG  */
      enddo
      if(no_of_remote_neighs.gt.jstack)
     .      neigh_morts(:,:,jstack+1:no_of_remote_neighs) = -1
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
        commatrix_recv = 0
        commatrix_send = 0
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
        call shmem_barrier_all()
        write(*,*) 'pe ',mype,' commatrix ',
     .             commatrix_recv(1:nprocs)


         write(*,'(" ")')
         write(*,'(" COMMUNICATION MATRIX2: morton_restrict")')
         write(*,'(" ")')
         write(*,'("pe  ",i3," commatrix_send ",
     .       2i3)') mype,(commatrix_send(i),i=1,nprocs)
         write(*,'(" ")')
        call shmem_barrier_all()
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
     .  'ERROR in mpi_morton_bnd_restrict: ',
     .  'guard block starting index',
     .  strt_buffer,' not larger than lnblocks',lnblocks
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif

#ifdef DEBUG
       write(*,*) 'pe ',mype,' end step 17 '
#endif /* DEBUG  */

!--------------------------------------------------
!
! Step 18.
! Dynamically allocate memory to store the lists of blocks to be
! sent and received.


      iprocs = max(maxval(commatrix_send),maxval(commatrix_recv))
      call MPI_ALLREDUCE(iprocs,
     .                   largest_no_of_blocks,
     .                   1,
     .                   MPI_INTEGER,
     .                   MPI_MAX,
     .                   MPI_COMM_WORLD,
     .                   ierror)

       if(allocated(to_be_sent)) 
     .    deallocate(to_be_sent)
       if(allocated(to_be_received)) deallocate(to_be_received)
       allocate( to_be_sent(3,
     .                      max(1,largest_no_of_blocks),
     .                      max(1,max_no_to_send) ) )
       allocate( to_be_received(3,
     .                          max(1,largest_no_of_blocks),
     .                          max(1,max_no_to_be_received) ) )      


#ifdef DEBUG
       write(*,*) 'pe ',mype,' end step 18 '
#endif /* DEBUG  */
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
                  to_be_received(1:2,ii,jp) = neigh_morts(6,1:2,jj)
                  to_be_received(  3,ii,jp) = neigh_morts(6,3,jj)
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
#endif /* DEBUG  */

#ifdef DEBUG
       write(*,*) 'pe ',mype,' end step 19a '
#endif /* DEBUG  */
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
            call Mpi_Irecv(to_be_sent(1,1,kk),isize,
     .           MPI_INTEGER,
     .           isrc ,itag,MPI_COMM_WORLD,recvrequest(kk),ierr)
#ifdef DEBUG
        write(*,*) 'pe ',mype,' posted receive with tag ',
     .          itag,' source ',isrc
#endif /* DEBUG  */
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
        write(*,*) 'pe ',mype,' posted send with tag ',
     .          itag,' destination ',idest
#endif /* DEBUG  */
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

#ifdef DEBUG
      write(*,*) 'pe ',mype,' calling mpi_amr_write_restrict_comm'
#endif /* DEBUG */

! Store communication info for future use
       call mpi_amr_write_restrict_comm(nprocs)

#ifdef DEBUG
      write(*,*) 'pe ',mype,' exiting mpi_morton_bnd_restrict'
#endif /* DEBUG */

#ifdef TIMING_MPI
              timer_mpi_morton_bnd_restrict =  
     .                          timer_mpi_morton_bnd_restrict
     .                          + mpi_wtime() - time1
#endif /* TIMING_MPI */
      return
      end subroutine mpi_morton_bnd_restrict
