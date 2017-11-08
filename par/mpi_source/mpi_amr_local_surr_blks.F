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
!#define DEBUG

      subroutine mpi_amr_local_surr_blks(mype,lb,nprocs,
     &                          max_no_of_blocks,
     &                          surrblks,l_parent,psurrblks)


!------------------------------------------------------------------------
!
! This routine calculates the addresses of surrounding blocks of 
! the block lb on the local processor, and of its parent.
!
!
! Written :     Peter MacNeice          June 2000
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

      include 'mpif.h'

      integer, intent(in)    ::  mype,lb,nprocs,max_no_of_blocks
      integer, intent(out)   ::  surrblks(:,:,:,:)
      integer, intent(out)   ::  psurrblks(:,:,:,:)
      logical, intent(in)    ::  l_parent

!------------------------------------------------------------------------
! local variables

      integer :: i,j,k
      real    :: dx,dy,dz, x0,y0,z0
      integer :: level
      integer :: morton(6)

      real    :: pbsize(3),pcoord(3)
      integer :: mort_neigh(6,3,3,3),pmort_neigh(6,3,3,3)
      integer :: mortref(6,2),pmortref(6,2)
      integer :: mort(6)

      integer :: i_pe,rem_block,rem_pe,rem_type
      integer :: ierrorcode,ierr,ibc
      logical :: lboundary
      logical :: lfound
      logical :: morton_greater_than
      logical :: morton_less_than
      logical :: morton_equal


#ifdef DEBUG_MORTON
      integer :: nbitshft,nbits,inxt_bit
#endif

!------------------------------------------------------------------------
!
!
! This routine assumes that the grid blocks are ordered by morton
! number and that any blocks with different refinement levels but
! the same morton number are ordered from coarse to fine.

!------------------------------------------------------------------------

      if(.not.morton_limits_set) then
        write(*,*) 'Error : mpi_amr_local_surr_blks : morton info '
     &            ,'is out of date or was never set up. Make sure '
     &            ,'there is a call to mpi_morton_limits before '
     &            ,'this routine is called. '
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif

!----------------------------

      surrblks  = -1
      psurrblks = -1

!
! Compute xmin,ymin,zmin,xmax,ymax,zmax or get them from storage
      xmin = grid_xmin
      ymin = grid_ymin
      zmin = grid_zmin
      xmax = grid_xmax
      ymax = grid_ymax
      zmax = grid_zmax


      dx = bsize(1,lb)
      dy = bsize(2,lb)
      dz = bsize(3,lb)

! set the block center
      x0 = coord(1,lb) - xmin
      y0 = coord(2,lb) - ymin
      z0 = coord(3,lb) - zmin

! compute the morton number for block lb 

!#ifdef SAVE_MORTS
!      morton = surr_morts(:,2,1+k2d,1+k3d,lb)
!#else
      call morton_number (x0,y0,z0,bsize(:,lb),ndim,lrefine_max,
     .                   lrefine(lb),morton)
!#endif /* SAVE_MORTS */

#ifdef DEBUG
      write(*,*) mype,' local_surr} ',
     .         'current block ',lb,mype,
     .         ' morton no before shift ',morton,
     .         ' x0 y0 z0 ',z0,y0,z0
#endif /* DEBUG */
      mortref(1:6,1) = morton(:)
      mortref(1:6,2) = lrefine(lb)

#ifdef DEBUG
      write(*,*) mype,' local_surr} ',
     .         'current block ',lb,mype,
     .         ' morton no after shift ',morton,
     .         ' x0 y0 z0 ',z0,y0,z0
#endif /* DEBUG */

      if(l_parent) then
      if(parent(1,lb).gt.0) then
        call morton_parent(mortref(:,1),mortref(1,2),lrefine_max,
     &                     ndim,pmortref(:,1))
        pmortref(:,2) = lrefine(lb)-1 
      endif
      endif

!----------------------------
!
! Construct a list of neighbors of this block

! First get the possible neighbors of the current block
!      mort_neigh = -1
!#ifdef SAVE_MORTS
!      mort_neigh(:,:,2-k2d:2+k2d,2-k3d:2+k3d) =
!     .      surr_morts(:,:,1:1+2*k2d,1:1+2*k3d,lb)
!#else /*  SAVE_MORTS */
      call morton_neighbors(xmin,ymin,zmin,xmax,ymax,zmax,
     .                      lperiodicx,lperiodicy,lperiodicz,
     .                      coord(1,lb),bsize(1,lb),ndim,
     .                      lrefine(lb),lrefine_max,mort_neigh)
!#endif /* SAVE_MORTS */

#ifdef DEBUGX
        if(mype.eq.1.and.lb.eq.9)
     . write(*,*) mype,' local_surr ',
     .         'current block ',lb,mype,
     .         ' mort_neigh ',mort_neigh
#endif /* DEBUG */
      if(l_parent) then
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
        
!        pmort_neigh = -1
!#ifdef SAVE_MORTS
!      if(parent(2,lb).eq.mype) then
!        pmort_neigh(:,:,2-k2d:2+k2d,2-k3d:2+k3d) =
!     .      surr_morts(:,:,1:1+2*k2d,1:1+2*k3d,parent(1,lb))
!      else
!        call morton_neighbors(xmin,ymin,zmin,xmax,ymax,zmax,
!     .                        lperiodicx,lperiodicy,lperiodicz,
!     .                        pcoord(1),pbsize(1),ndim,
!     .                        lrefine(lb)-1,lrefine_max,pmort_neigh)
!      endif
!#else /*  SAVE_MORTS */
        call morton_neighbors(xmin,ymin,zmin,xmax,ymax,zmax,
     .                        lperiodicx,lperiodicy,lperiodicz,
     .                        pcoord(1),pbsize(1),ndim,
     .                        lrefine(lb)-1,lrefine_max,pmort_neigh)
!#endif /* SAVE_MORTS */
      endif
      endif

          
!----------------------------
!
! Loop over this blocks^s list of required neighbor blocks,
! identifying their remote location from the morton information received
! above

! First get the block addresses of all elements in the mort_neigh array

        do k = 2-k3d,2+k3d
        do j = 2-k2d,2+k2d
        do i = 1,3
        mort(:)  = mort_neigh(:,i,j,k)
        level = lrefine(lb)

#ifdef DEBUGX
        if(mype.eq.1.and.lb.eq.9)
     .   write(*,*) mype,' local_surr} ',
     .           'searching for mort/level ',mort,level
#endif /* DEBUG */

        if( mort(6).gt.-1 .and. (i.ne.2.or.j.ne.2.or.k.ne.2) ) then

          j_pe = -1
          do i_pe = 1,nprocs
            if( 
     .          ( ( morton_greater_than(mort(:),
     .                            morton_limits(:,1,1,i_pe)) ) .or. 
     .            ( morton_equal(mort(:),
     .                           morton_limits(:,1,1,i_pe)) .and. 
     .              level.ge.morton_limits(6,2,1,i_pe) ) )
     .                      .and.
     .          ( ( morton_less_than(mort(:),
     .                               morton_limits(:,1,2,i_pe)) ) .or. 
     .            ( morton_equal(mort(:),
     .                           morton_limits(:,1,2,i_pe)) .and. 
     .              level.le.morton_limits(6,2,2,i_pe) ) ) 
     .        ) then
            j_pe = i_pe
#ifdef DEBUGX
        if(mype.eq.1.and.lb.eq.9)
     .   write(*,*) mype,' local_surr} ',
     .   'found mort/level ',mort,level,' on pe ',j_pe
#endif /* DEBUG */
          endif
          enddo


          rem_type = -1
          rem_block = -1
          rem_pe = j_pe

          if(rem_pe.ne.-1) then
          if(rem_pe.ne.mype+1) then
            kk = -1
            do k1=1,no_of_mortonbnds_received
              if(pe_source(k1).eq.rem_pe) kk = k1 
            enddo
            if(kk.gt.0) then

            lfound = .false.
            j1 = 0
            do while (.not.lfound.and.j1.lt.max_no_of_blocks)
              j1 = j1+1
              if( morton_equal(r_mortonbnd(:,1,j1,kk),mort(:)) 
     .                        .and.
     .            r_mortonbnd(6,2,j1,kk).eq.level ) then
                  rem_block = j1
                  rem_type = r_mortonbnd(6,3,j1,kk)
                  lfound = .true.
              endif

#ifdef DEBUGX
        if(mype.eq.1.and.lb.eq.9)
     .   write(*,*) mype,' local_surr} ',
     .    'address of mort/level ',mort,level,' is ',
     .              rem_block,rem_pe
#endif /* DEBUG */
            enddo
            endif
          else
            lfound = .false.
            j1 = 0
            do while (.not.lfound.and.j1.lt.lnblocks)
              j1 = j1+1
              if( morton_equal(mortonbnd(:,1,j1),mort(:)) 
     .                          .and.
     .            mortonbnd(6,2,j1).eq.level ) then
                  rem_block = j1
                  rem_type = mortonbnd(6,3,j1)
                  lfound = .true.
              endif

#ifdef DEBUGX
        if(mype.eq.1.and.lb.eq.9)
     .   write(*,*) mype,' local_surr} ',
     .   'address of mort/level ',mort,level,' is ',
     .              rem_block,rem_pe
#endif /* DEBUG */
            enddo
          endif
          endif
          if(rem_block.eq.-1) rem_pe = -1

! surrblks is used to store 
! the addresses of the existing neighbors.
! Here proc nos. run from 0 to nprocs-1
          surrblks(1,i,j,k) = rem_block
          if(rem_pe.gt.-1) surrblks(2,i,j,k) = rem_pe-1
          surrblks(3,i,j,k) = rem_type

!        else
!          surrblks(:,i,j,k) = -1
!          x0 = coord(1,lb)+bsize(1,lb)*real(i-2)
!          y0 = coord(2,lb)+bsize(2,lb)*real(j-2)
!          z0 = coord(3,lb)+bsize(3,lb)*real(k-2)
!          call boundary_locator(x0,y0,z0,lboundary,ibc)
!          if(lboundary) surrblks(:,i,j,k) = ibc
        endif
        enddo
        enddo
        enddo

! Search for boundaries

        do k = 2-k3d,2+k3d
        do j = 2-k2d,2+k2d
        do i = 1,3
        mort(:)  = mort_neigh(:,i,j,k)
        level = lrefine(lb)
        if (surrblks(1,i,j,k).eq.-1.or.mort(6).eq.-1) then
          x0 = coord(1,lb)+bsize(1,lb)*real(i-2)
          y0 = coord(2,lb)+bsize(2,lb)*real(j-2)
          z0 = coord(3,lb)+bsize(3,lb)*real(k-2)
          call boundary_locator(x0,y0,z0,lboundary,ibc)
          if(lboundary) then
            surrblks(1:2,i,j,k) = ibc
            surrblks(3,i,j,k) = -1
          endif
        end if
        end do
        end do
        end do

        surrblks(1,2,2,2) = lb
        surrblks(2,2,2,2) = mype
        surrblks(3,2,2,2) = nodetype(lb)


      if(l_parent) then
! Now get the block addresses of all elements in the pmort_neigh array
        if(parent(1,lb).gt.0) then

        do k = 2-k3d,2+k3d
        do j = 2-k2d,2+k2d
        do i = 1,3
        mort(:)  = pmort_neigh(:,i,j,k)
        level = lrefine(lb)-1

        if( mort(6).gt.-1 .and. (i.ne.2.or.j.ne.2.or.k.ne.2) ) then
          i_pe = 1
          do while( 
     .          morton_greater_than(mort(:),morton_limits(:,1,2,i_pe))
     .                               .or.
     .         (morton_equal(mort(:),morton_limits(:,1,2,i_pe)).and.
     .          level.gt.morton_limits(6,2,2,i_pe)    ) 
     .            )
             i_pe = i_pe + 1
          enddo

          rem_type = -1
          rem_block = -1
          rem_pe = i_pe

          if(rem_pe.ne.mype+1) then
            do k1=1,no_of_mortonbnds_received
              if(pe_source(k1).eq.rem_pe) kk = k1 
            enddo
            do j1=1,max_no_of_blocks
              if( morton_equal(r_mortonbnd(:,1,j1,kk),mort(:))
     .                          .and.
     .            r_mortonbnd(6,2,j1,kk).eq.level ) then
                  rem_block = j1
                  rem_type = r_mortonbnd(6,3,j1,kk)
              endif
            enddo
          else
            do j1=1,lnblocks
              if( morton_equal(mortonbnd(:,1,j1),mort(:))
     .                         .and.
     .            mortonbnd(6,2,j1).eq.level ) then
                  rem_block = j1
                  rem_type = r_mortonbnd(6,3,j1,kk)
              endif
            enddo
          endif

! surrblks is used to store 
! the addresses of the existing neighbors.
! Here proc nos. run from 0 to nprocs-1
          psurrblks(1,i,j,k) = rem_block
          psurrblks(2,i,j,k) = rem_pe-1
          psurrblks(3,i,j,k) = rem_type

        endif

        enddo
        enddo
        enddo

! Search for boundaries

        do k = 2-k3d,2+k3d
        do j = 2-k2d,2+k2d
        do i = 1,3
        mort(:)  = pmort_neigh(:,i,j,k)
        level = lrefine(lb)-1
        if( mort(6).eq.-1 .or. psurrblks(1,i,j,k) .eq. -1) then
          psurrblks(:,i,j,k) = -1
          x0 = pcoord(1)+pbsize(1)*real(i-2)
          y0 = pcoord(2)+pbsize(2)*real(j-2)
          z0 = pcoord(3)+pbsize(3)*real(k-2)
          call boundary_locator(x0,y0,z0,lboundary,ibc)
          if(lboundary) then
            psurrblks(1:2,i,j,k) = ibc
            psurrblks(3,i,j,k) = -1
          endif
        end if
        end do
        end do
        end do

        psurrblks(1:2,2,2,2) = parent(1:2,lb)
        psurrblks(3  ,2,2,2) = 2               ! Not set right yet. Cannot detect
                                               ! parent type 3 in this way

        endif
      endif
!----------------------------

      return
      end subroutine mpi_amr_local_surr_blks



      subroutine mpi_amr_local_surr_blks_lkup(mype,lb,
     &                          surrblks,l_parent,psurrblks)

!#define DEBUG


!------------------------------------------------------------------------
!
! This routine calculates the addresses of surrounding blocks of 
! the block lb on the local processor, and of its parent.
!
!
! Written :     Peter MacNeice          June 2000
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

      include 'mpif.h'


      integer, intent(in)    ::  mype,lb
      integer, intent(out)   ::  surrblks(:,:,:,:)
      integer, intent(out)   ::  psurrblks(:,:,:,:)
      logical, intent(in)    ::  l_parent

!------------------------------------------------------------------------
! local variables

      integer :: iblk
      integer :: remote_block,remote_pe
      integer :: ierrorcode,ierr
      logical :: lfound

!------------------------------------------------------------------------
!
!
! This routine assumes that the grid blocks are ordered by morton
! number and that any blocks with different refinement levels but
! the same morton number are ordered from coarse to fine.

!------------------------------------------------------------------------

      surrblks = -1
      psurrblks = -1

      surrblks(:,:,2-k2d:2+k2d,2-k3d:2+k3d) =
     .       surr_blks(:,:,1:1+2*k2d,1:1+2*k3d,lb)


      if( l_parent .and. (parent(1,lb).gt.0) ) then

          lfound = .false.
          if(parent(2,lb).ne.mype) then

            iblk = ladd_strt(parent(2,lb))
            do while(.not.lfound.and.
     .               iblk.le.ladd_end(parent(2,lb)))

!            iblk = strt_buffer
!            do while(.not.lfound.and.iblk.le.last_buffer)
#ifdef DEBUG
             write(*,*) mype,' local_surr} ',
     .            ' searching buffer for ',
     .            parent(:,lb),' current buffer entry ',
     .          ' iblk ',iblk,' laddress ',laddress(:,iblk)
#endif /* DEBUG */
              if(parent(1,lb).eq.laddress(1,iblk).and.
     .           parent(2,lb).eq.laddress(2,iblk) ) then
                remote_block = iblk
                remote_pe    = mype
                lfound = .true.
#ifdef DEBUG
             write(*,*) mype,' local_surr} ',
     .        ' remote block ',
     .          parent(:,lb),' located in buffer slot ',
     .          iblk,' has surr_blks ',surr_blks(:,:,:,:,iblk)
#endif /* DEBUG */
              else
                iblk = iblk+1
              endif
            enddo

          else

            remote_block = parent(1,lb)
            remote_pe    = mype
            if(remote_block.le.lnblocks) lfound = .true.

          endif

        if(.not.lfound) then
          write(*,*) 'Error in mpi_amr_local_surr_blks_lkup : ',
     .               'remote block ',parent(:,lb),
     .               ' not located on pe ',mype,
     .               ' while processing blk ',lb,mype
          call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
        endif

        if(lfound) then
          psurrblks(:,:,2-k2d:2+k2d,2-k3d:2+k3d) =
     .         surr_blks(:,:,1:1+2*k2d,1:1+2*k3d,remote_block)
        endif

      endif

!----------------------------

      return
      end subroutine mpi_amr_local_surr_blks_lkup
