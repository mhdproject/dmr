!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"

      subroutine amr_refine_blocks (nprocs,mype)





c By K. Olson (NASA/GSFC and GMU), 11/96
c modified 2/2001 By M. Zingale and J. Dursi (U. of Chicago)

      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton
#ifdef TIMING_MPI
      use timings
#endif /* TIMING_MPI */

      implicit none

      include 'mpif.h'

      integer, intent(in)    :: nprocs,mype

! local variables -------------------------------------------------------

      real :: hx,hy,hz
      real :: hh(mdim)      

      integer :: i,j,k,jj,ipar,ipar_proc
      integer :: ineigh,ineigh_proc
      integer :: lnblocks2,ichi
      integer :: ierr,jr,n3
      integer :: nsend,nrecv
      integer :: istart,iend,istart2,iend2,istep
      integer :: nsend_to_proc(0:16384), nrecv_pack(0:16384)
      integer :: errorcode
      integer :: reqr(maxblocks_tr)
      integer :: kk(maxblocks_tr)
      integer :: nodetype_chi(nchild,maxblocks_tr)
      integer :: statr(MPI_STATUS_SIZE,maxblocks_tr)
      integer :: tneigh(2,mfaces,maxblocks_tr)
      integer :: idummy_array(1)
      integer :: iii,ii,jjj,kkk

      logical :: refine_neigh(maxblocks_tr)
#ifdef TIMING_MPI
      double precision :: time1
#endif /* TIMING_MPI */

! ------------------------------------------------------------------------
#ifdef SAVE_MORTS
      real   :: xmin,ymin,zmin,xmax,ymax,zmax
      integer:: mort_neigh(6,3,3,3)

! Find the coordinate ranges (TO DO : This particular call need only be made at the start of the run,
! since it is always updated at the end of refine_derefine.)
#ifdef TIMING_MPI
      time1 = mpi_wtime()
#endif /* TIMING_MPI */
      call mpi_amr_global_domain_limits
#ifdef TIMING_MPI
      timer_amr_global_domain_limits =
     .     timer_amr_global_domain_limits + mpi_wtime() - time1
#endif /* TIMING_MPI */

! Compute xmin,ymin,zmin,xmax,ymax,zmax or get them from storage
      xmin = grid_xmin
      ymin = grid_ymin
      zmin = grid_zmin
      xmax = grid_xmax
      ymax = grid_ymax
      zmax = grid_zmax
#endif /* SAVE_MORTS */

c 1) refine blocks, create their children, the children's parents,
c    turn the children on, and turn the parents off

      do i = 1,maxblocks_tr
         newchild(i) = .FALSE.
      end do

!!!!!!!!!!!!!!

      lnblocks2 = lnblocks

      do i = 1,lnblocks

         if (refine(i)) then ! refine block 'i' on this processor

            do j = 1,nchild ! create 'nchild' child blocks

               lnblocks2 = lnblocks2 + 1

      if(lnblocks2.gt.maxblocks_tr) then
       open(unit=30,file='amr_log',status='unknown',
     &        position='append')
       write(30,*) 'PARAMESH ERROR !'
       write(30,*) 'Too many blocks created! '
       write(30,*) 'Increase MAXBLOCKS_TR and rerun! '
       close(30)
       write(*,*) 'PARAMESH ERROR !'
       write(*,*) 'Too many blocks created! '
       write(*,*) 'Increase MAXBLOCKS_TR and rerun! '
       call MPI_ABORT(MPI_COMM_WORLD,errorcode,ierr)
      endif

               child(1,j,i) = lnblocks2 ! child j's on-processor id
               child(2,j,i) = mype      ! child j's processor no.

               lrefine(lnblocks2) = lrefine(i) + 1 ! refinement level of child
               newchild(lnblocks2) = .TRUE. ! this is a new child

               parent(1,lnblocks2) = i ! child j's parent
               parent(2,lnblocks2) = mype ! child j's parent's processor

               bnd_box(:,:,lnblocks2) = bnd_box(:,:,i)

               ! ordering of children in space is morton ordered !

               hh(:) = bsize(:,i)/4.
               if (j.eq.1) then
                  hx = -1.
                  hy = -1.
                  hz = -1.
                  bnd_box(2,1,lnblocks2) = coord(1,i)
                  bnd_box(2,2,lnblocks2) = coord(2,i)
                  bnd_box(2,3,lnblocks2) = coord(3,i)
               else if (j.eq.2) then
                  hx = 1.
                  hy = -1.
                  hz = -1.
                  bnd_box(1,1,lnblocks2) = coord(1,i)
                  bnd_box(2,2,lnblocks2) = coord(2,i)
                  bnd_box(2,3,lnblocks2) = coord(3,i)
               else if (j.eq.3) then
                  hx = -1.
                  hy = 1.
                  hz = -1.
                  bnd_box(2,1,lnblocks2) = coord(1,i)
                  bnd_box(1,2,lnblocks2) = coord(2,i)
                  bnd_box(2,3,lnblocks2) = coord(3,i)
               else if (j.eq.4) then
                  hx = 1.
                  hy = 1.
                  hz = -1.
                  bnd_box(1,1,lnblocks2) = coord(1,i)
                  bnd_box(1,2,lnblocks2) = coord(2,i)
                  bnd_box(2,3,lnblocks2) = coord(3,i)
               else if (j.eq.5) then
                  hx = -1.
                  hy = -1.
                  hz = 1.
                  bnd_box(2,1,lnblocks2) = coord(1,i)
                  bnd_box(2,2,lnblocks2) = coord(2,i)
                  bnd_box(1,3,lnblocks2) = coord(3,i)
               else if (j.eq.6) then
                  hx = 1.
                  hy = -1.
                  hz = 1.
                  bnd_box(1,1,lnblocks2) = coord(1,i)
                  bnd_box(2,2,lnblocks2) = coord(2,i)
                  bnd_box(1,3,lnblocks2) = coord(3,i)
               else if (j.eq.7) then
                  hx = -1.
                  hy = 1.
                  hz = 1.
                  bnd_box(2,1,lnblocks2) = coord(1,i)
                  bnd_box(1,2,lnblocks2) = coord(2,i)
                  bnd_box(1,3,lnblocks2) = coord(3,i)
               else if (j.eq.8) then
                  hx = 1.
                  hy = 1.
                  hz = 1.
                  bnd_box(1,1,lnblocks2) = coord(1,i)
                  bnd_box(1,2,lnblocks2) = coord(2,i)
                  bnd_box(1,3,lnblocks2) = coord(3,i)
               end if

               if(ndim.lt.3) bnd_box(:,3,lnblocks2) = bnd_box(:,3,i)
               if(ndim.lt.2) bnd_box(:,2,lnblocks2) = bnd_box(:,2,i)

               which_child(lnblocks2) = j
                     
               coord(1,lnblocks2) = coord(1,i) + hx*hh(1)
               if (ndim.ge.2) then
                  coord(2,lnblocks2) = coord(2,i) + hy*hh(2)
               end if
               if (ndim.eq.3) then
                  coord(3,lnblocks2) = coord(3,i) + hz*hh(3)
               end if

               bsize(:,lnblocks2) = bsize(:,i)/2.

               bflags(:,lnblocks2) = bflags(:,i)

#ifdef EMPTY_CELLS
       if(empty(i).eq.1) empty(lnblocks2)=1
#endif
#ifdef SAVE_MORTS
      call morton_neighbors(xmin,ymin,zmin,xmax,ymax,zmax,
     .                      lperiodicx,lperiodicy,lperiodicz,
     .                      coord(:,lnblocks2),bsize(:,lnblocks2),ndim,
     .                      lrefine(lnblocks2),lrefine_max,
     .                      mort_neigh)
      surr_morts(:,:,1:1+2*k2d,1:1+2*k3d,lnblocks2) = 
     .                           mort_neigh(:,:,2-k2d:2+k2d,2-k3d:2+k3d)
      write(*,*) 'pe ',mype,' added surr_morts for new block ',
     .  lnblocks2,
     . ' mort_neigh(6,:,2-k2d:2+k2d,2-k3d:2+k3d) ',
     .   mort_neigh(6,:,2-k2d:2+k2d,2-k3d:2+k3d)
#endif /* SAVE_MORTS */

            end do

#ifdef EMPTY_CELLS
       if(empty(i).eq.1) empty(i)=0
#endif
         end if

       end do

c Connect neighbors of newly created sub-blocks
c
c                    4              6
c                    |              |
c  in x-y plane  1 - i - 2 ; in z   i
c                    |              |
c                    3              5
c

! In some cases, old data in the parts of the neigh array
! about to be used by the new children can confuse the algorithm
! which searches for off processor neighbors. In particular, if
! an external boundary is not properly detected, then a spurious
! neighbor link can be set up in section XXX below.
! Reinitializing these parts of neigh avoids the confusion.
      neigh(:,:,lnblocks+1:lnblocks2) = -1
      child(:,:,lnblocks+1:lnblocks2) = -1


c 1) connect with siblings (which at this point are all on processor)
      
      do i = 1,lnblocks
         if (refine(i)) then
            do j = 1,nchild ! cycle through children
               ichi = child(1,j,i)
               do k = 1,nfaces
                  neigh(2,k,ichi) = mype
               end do
            end do

            ! connect in x direction
               
            do j = 2,nchild,2
               ichi = child(1,j,i)
               k = j - 1
               ineigh = child(1,k,i)
               neigh(1,1,ichi) = ineigh
            end do
            do j = 1,nchild-1,2
               ichi = child(1,j,i)
               k = j + 1
               ineigh = child(1,k,i)
               neigh(1,2,ichi) = ineigh
            end do

            ! connect in y direction

            if(ndim.ge.2) then
            do j = 3,4
               ichi = child(1,j,i)
               k = j - 2
               ineigh = child(1,k,i)
               neigh(1,3,ichi) = ineigh
            end do
            if (ndim.eq.3) then
            do j = 7,8
               ichi = child(1,j,i)
               k = j - 2
               ineigh = child(1,k,i)
               neigh(1,3,ichi) = ineigh
            end do
            end if
            do j = 1,2
               ichi = child(1,j,i)
               k = j + 2
               ineigh = child(1,k,i)
               neigh(1,4,ichi) = ineigh
            end do
            if (ndim.eq.3) then
            do j = 5,6
               ichi = child(1,j,i)
               k = j + 2
               ineigh = child(1,k,i)
               neigh(1,4,ichi) = ineigh
            end do
            end if
            end if

!!!#if N_DIM == 3
            if (ndim == 3) then

            ! connect in z direction

            if (ndim.eq.3) then
            do j = 5,8
               ichi = child(1,j,i)
               k = j - 4
               ineigh = child(1,k,i)
               neigh(1,5,ichi) = ineigh
            end do
            do j = 1,4
               ichi = child(1,j,i)
               k = j + 4
               ineigh = child(1,k,i)
               neigh(1,6,ichi) = ineigh
            end do
            end if
!!!#endif /*N_DIM*/
            end if
         end if
      end do


c 2) connect with off-processor neighbors by looking at neighbors of parent

      tneigh(:,:,:) = 0

! Send refine flags to neighbors

      do j = 1,nfaces
         if (j.eq.1) then
            jr = 2
            istart = 1
            iend = nchild-1
            istart2 = 2
            iend2 = nchild
            istep = 2
         elseif (j.eq.2) then
            jr = 1
            istart = 2
            iend = nchild
            istart2 = 1
            iend2 = nchild-1
            istep = 2
         elseif (j.eq.3) then
            jr = 4
            istart = 1
            iend = 2
            istart2 = 3
            iend2 = 4
            istep = 1
         elseif (j.eq.4) then
            jr = 3
            istart = 3
            iend = 4
            istart2 = 1
            iend2 = 2
            istep = 1
         elseif (j.eq.5) then
            jr = 6
            istart = 1
            iend = 4
            istart2 = 5
            iend2 = 8
            istep = 1
         elseif (j.eq.6) then
            jr = 5
            istart = 5
            iend = 8
            istart2 = 1
            iend2 = 4
            istep = 1
         end if

         refine_neigh(:) = .FALSE.
         refine(lnblocks+1:maxblocks_tr) = .FALSE.

         nrecv = 0
         do i = 1,lnblocks
            if (neigh(1,jr,i).gt.0) then
               if (neigh(2,jr,i).ne.mype) then
                  nrecv = nrecv + 1
                  call MPI_logical_IRECV(refine_neigh(i),
     &                           1,
     &                           MPI_LOGICAL,
     &                           neigh(2,jr,i),
     &                           neigh(1,jr,i),
     &                           MPI_COMM_WORLD,
     &                           reqr(nrecv),
     &                           ierr)
               end if
            end if
         end do  

         nsend = 0
         do i = 1,lnblocks
            ineigh = neigh(1,j,i)
            ineigh_proc = neigh(2,j,i)
            if (ineigh.ge.1) then
               if (ineigh_proc.ne.mype) then
                  nsend = nsend + 1
                  call MPI_logical_SSEND(refine(i),
     &                           1,
     &                           MPI_LOGICAL,
     &                           ineigh_proc,
     &                           i,
     &                           MPI_COMM_WORLD,
     &                           ierr)
               else
                  refine_neigh(ineigh) = refine(i)
               end if
            end if
         end do

         if (nrecv.gt.0) then
            call MPI_WAITALL(nrecv,reqr,statr,ierr)
         end if
 100     continue

         do jj = istart,iend,istep
            
            if (j.eq.1) then
               k = jj+1         ! child no. of neighbor which lies on border
            elseif (j.eq.2) then
               k = jj-1
            elseif (j.eq.3) then
               k = jj+2
            elseif (j.eq.4) then
               k = jj-2
            elseif (j.eq.5) then
               k = jj+4
            elseif (j.eq.6) then
               k = jj-4
            end if

         nrecv = 0
         do i = 1,lnblocks
            if (refine(i)) then
               ineigh = neigh(1,j,i) ! neighbor of parent
               ineigh_proc = neigh(2,j,i) ! neighbor of parent
               if (ineigh.gt.0) then
                                            ! 'i' on left side
                  ichi = child(1,jj,i) ! child 'jj' of parent 'i'

                     ! fetch child no. k of ineigh from ineigh_proc
                     ! and store in neighbor 1 of child jj (ichi)

                  if (ineigh_proc.ne.mype) then
                     nrecv = nrecv + 1
                     call MPI_int_IRECV(neigh(1,j,ichi),
     &                              2,
     &                              MPI_INTEGER,
     &                              neigh(2,j,i),
     &                              neigh(1,j,i),
     &                              MPI_COMM_WORLD,
     &                              reqr(nrecv),
     &                              ierr)
                  else
                     neigh(:,j,ichi) = child(:,k,ineigh)
                  end if
               end if
            end if
         end do


         nsend = 0
         do i = 1,lnblocks
            if (refine_neigh(i)) then
               ineigh = neigh(1,jr,i) ! right neighbor of parent
               ineigh_proc = neigh(2,jr,i)
               if (ineigh.gt.0) then
                  ! send child no. jj of block i to neigh on right
                  if (ineigh_proc.ne.mype) then
                     nsend = nsend + 1
                     call MPI_int_SSEND(child(1,k,i),
     &                    2,
     &                    MPI_INTEGER,
     &                    neigh(2,jr,i),
     &                    i,
     &                    MPI_COMM_WORLD,
     &                    ierr)
                  end if
               end if
            end if
         end do
         
         if (nrecv.gt.0) then
            call MPI_WAITALL(nrecv,reqr,statr,ierr)
         end if

         end do ! end loop over jj

         if (j.eq.3.and.istart.eq.1.and.ndim.eq.3) then
            istart = 5
            iend = 6
            istart2 = 7
            iend2 = 8
            go to 100
         elseif (j.eq.4.and.istart.eq.3.and.ndim.eq.3) then
            istart = 7
            iend = 8
            istart2 = 5
            iend2 = 6
            go to 100
         end if

! 1) count no. of sends on each proc. to all other procs

         nsend_to_proc(:) = 0
         do i = lnblocks+1,lnblocks2
            if (newchild(i)) then
               ineigh = neigh(1,j,i)
               ineigh_proc = neigh(2,j,i)
               if (ineigh.ge.1) then
                  if (ineigh_proc.ne.mype) then
                     nsend_to_proc(ineigh_proc) = 
     &                    nsend_to_proc(ineigh_proc) + 1
                  end if
               end if
            end if
         end do

! 2) collect data for `this' proc from other procs so that
!    the total no. of receives to post can be computed

! Change made by M. Zingale and J. Dursi -------------------------------------
         nrecv = 0
         call MPI_int_AllReduce(nsend_to_proc, nrecv_pack, nprocs,
     +        MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
         nrecv = nrecv_pack(MyPE)

!-----------------------------------------------------------------------------

         do i = 1,nrecv
            call MPI_int_IRECV(kk(i),
     &                     1,
     &                     MPI_INTEGER,
     &                     MPI_ANY_SOURCE,
     &                     MPI_ANY_TAG,
     &                     MPI_COMM_WORLD,
     &                     reqr(i),
     &                     ierr)
         end do  

! section XXX
! NOW new children send to their new neighbors so that they can be set
                     
         nsend = 0
         do i = lnblocks+1,lnblocks2
            if (newchild(i)) then
               ineigh = neigh(1,j,i)
               ineigh_proc = neigh(2,j,i)
               if (ineigh.ge.1) then
                  if (ineigh_proc.ne.mype) then
                     nsend = nsend + 1
                     idummy_array(1) = ineigh
                     call MPI_int_SSEND(idummy_array(1),
     &                              1,
     &                              MPI_INTEGER,
     &                              ineigh_proc,
     &                              i,
     &                              MPI_COMM_WORLD,
     &                              ierr)
                  else
                     tneigh(1,jr,ineigh) = i
                     tneigh(2,jr,ineigh) = mype
                  end if
               end if
            end if
         end do

! it seems to work only if tneigh is used ????

         if (nrecv.gt.0) then
            call MPI_WAITALL(nrecv,reqr,statr,ierr)
            do i = 1,nrecv
               tneigh(1,jr,kk(i)) = statr(MPI_TAG,i)
               tneigh(2,jr,kk(i)) = statr(MPI_SOURCE,i)
            end do
         end if
         call MPI_BARRIER(MPI_COMM_WORLD,ierr) ! NEEDED ????

      end do ! end loop over faces

                     
      do i = 1,lnblocks2
         do j = 1,nfaces
            if (tneigh(1,j,i).gt.0.and..not.newchild(i).and.
     &           tneigh(1,j,i).ne.neigh(1,j,i)) then
               neigh(1,j,i) = tneigh(1,j,i)
               neigh(2,j,i) = tneigh(2,j,i)
          
            end if
         end do
      end do
      

      lnblocks = lnblocks2

c reset node types
! Do we really need this since it gets done completely in amr_derefine_blocks
      
      do i = 1,lnblocks
         if (newchild(i)) nodetype(i) = 1
         if (refine(i)) then
            nodetype(i) = 2
         end if
      end do

      nrecv = 0
      nodetype_chi(:,1:lnblocks) = 0
      do i = 1,lnblocks
         do j = 1,nchild
            if (child(1,j,i).gt.0) then
            if (child(2,j,i).ne.mype) then
               nrecv = nrecv + 1
               call MPI_int_IRECV(nodetype_chi(j,i),
     &                        1,
     &                        MPI_INTEGER,
     &                        child(2,j,i),
     &                        child(1,j,i),
     &                        MPI_COMM_WORLD,
     &                        reqr(nrecv),
     &                        ierr)
            else
               nodetype_chi(j,i) = nodetype(child(1,j,i))
            end if
            end if
         end do
      end do  

      nsend = 0
      do i = 1,lnblocks
         if (parent(1,i).gt.0) then
         if (parent(2,i).ne.mype) then
            nsend = nsend + 1
            call MPI_int_SSEND(nodetype(i),
     &                     1,
     &                     MPI_INTEGER,
     &                     parent(2,i),
     &                     i,
     &                     MPI_COMM_WORLD,
     &                     ierr)
         end if
         end if
      end do

      if (nrecv.gt.0) then
         call MPI_WAITALL(nrecv,reqr,statr,ierr)
      end if
      do i = 1,lnblocks
         n3 = 0
         do j = 1,nchild
            if (nodetype_chi(j,i).ne.1) then
               n3 = n3 + 1
            end if
            if (n3.eq.nchild) nodetype(i) = 3
         end do
      end do 

c Now set neighbor pointers of new children if they lie on a boundary

      do i = 1,lnblocks

         if (newchild(i)) then

            ! fetch i's parent

            ipar = parent(1,i)
            ipar_proc = parent(2,i)

            do j = 1,nfaces

               if (neigh(1,j,i).le.-1) then ! this neighbor may be on a border

                  ! fetch i's parent's neighbor j

! Here we know that all new children on the same processor as
! the parent so we don't need any communications !!!

                  ineigh = neigh(1,j,ipar)
               
                  ! if parent's neighbor is lt -1 then i's neighbor is
                  ! also on the domain border and is set to the parent's
                  ! value

                  if (ineigh.le.-20) neigh(1,j,i) = ineigh

               end if

            end do

         end if

      end do

 56   continue

c reset refine flags

      do i = 1,maxblocks_tr
         refine(i) = .FALSE.
      end do

      return
      end subroutine amr_refine_blocks

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine amr_check_refine(nprocs,mype,icontinue)

      use paramesh_dimensions
      use physicaldata
      use tree

      implicit none

      integer, intent(in)    :: nprocs,mype
      integer, intent(in)    :: icontinue

      include 'mpif.h'

! local variables ----------------------------------------------------------

      integer :: i,j,ipar,ipar_proc
      integer :: ineigh,ineigh_proc
      integer :: ix, iy, iz, ix_n, iy_n, iz_n
      integer :: ierr,jr
      integer :: nsend,nrecv
      integer :: reqr(maxblocks_tr)
      integer :: statr(MPI_STATUS_SIZE,maxblocks_tr)

      logical :: repeat,repeat_t
      logical :: ref_testt(2,1+k2d,1+k3d,maxblocks_tr)
      logical :: ref_test(2,1+k2d,1+k3d,maxblocks_tr)
      logical :: lreduce_datain(1), lreduce_dataout(1)

! --------------------------------------------------------------------------



! Turn off refine flags of non-leaf blocks
      do i = 1,lnblocks
         if (nodetype(i).gt.1.and.refine(i)) refine(i) = .false.
      end do


! If no leaf blocks were marked for refinement then return.
 21   if(icontinue.eq.0) return


!!!!!!!!!!      
! CHECK FOR neighboring blocks which differ by more than 1 level of refinement
      
      do i = 1,lnblocks
        do iz = 1,1+k3d
          do iy = 1,1+k2d
            do ix = 1,2
              ref_test(ix,iy,iz,i) = .FALSE.
            end do
          end do
        end do
      end do

! Step 1 - set up an array, ref_test, which contains the 
!          logical refine flag for each child of every local block

      nrecv = 0
      do i = 1,lnblocks
         do j = 1,nchild
            if (child(1,j,i).gt.0) then

              if (j.eq.1) then
                 ix = 1
                 iy = 1
                 iz = 1
              elseif (j.eq.2) then
                 ix = 2
                 iy = 1
                 iz = 1
              elseif (j.eq.3) then
                 ix = 1
                 iy = 2
                 iz = 1
              elseif (j.eq.4) then
                 ix = 2
                 iy = 2
                 iz = 1
              elseif (j.eq.5) then
                 ix = 1
                 iy = 1
                 iz = 2
              elseif (j.eq.6) then
                 ix = 2
                 iy = 1
                 iz = 2
              elseif (j.eq.7) then
                 ix = 1
                 iy = 2
                 iz = 2
              elseif (j.eq.8) then
                 ix = 2
                 iy = 2
                 iz = 2
              endif

               if (child(2,j,i).ne.mype) then
                  nrecv = nrecv + 1
                  call MPI_logical_IRECV(ref_test(ix,iy,iz,i),
     &                           1,
     &                           MPI_LOGICAL,
     &                           child(2,j,i),
     &                           child(1,j,i),
     &                           MPI_COMM_WORLD,
     &                           reqr(nrecv),
     &                           ierr)
               else
                  ref_test(ix,iy,iz,i) = refine(child(1,j,i))
               end if
            end if
         end do
      end do
           
      nsend = 0
      do i = 1,lnblocks
!         if (refine(i)) then
            ipar = parent(1,i)
            ipar_proc = parent(2,i)
            if (ipar.ge.1) then
               if (ipar_proc.ne.mype) then
                  nsend = nsend + 1
                  call MPI_logical_SSEND(refine(i),
     &                           1,
     &                           MPI_LOGICAL,
     &                           ipar_proc,
     &                           i,
     &                           MPI_COMM_WORLD,
     &                           ierr)
               end if
            end if
!         end if
      end do

      if (nrecv.gt.0) then
        call MPI_WAITALL(nrecv,reqr,statr,ierr)
      end if


! Step 2 - now cycle over all faces, and for each set up the 
!          array, ref_testt, which contains the logical refine flag
!          for each child of the neighbor block across this face.


      do j = 1,nfaces     

         do i = 1,lnblocks
            ref_testt(:,:,:,i) = .FALSE.
         end do

         if (j.eq.1) then
            jr = 2
         elseif (j.eq.2) then
            jr = 1
         elseif (j.eq.3) then
            jr = 4
         elseif (j.eq.4) then
            jr = 3
         elseif (j.eq.5) then
            jr = 6
         elseif (j.eq.6) then
            jr = 5
         end if
        
 
         nrecv = 0
         do i = 1,lnblocks
            if (neigh(1,jr,i).gt.0) then
               if (neigh(2,jr,i).ne.mype) then
                  nrecv = nrecv + 1
                  call MPI_logical_IRECV(ref_testt(1,1,1,i),
     &                           2**ndim,
     &                           MPI_LOGICAL,
     &                           neigh(2,jr,i),
     &                           neigh(1,jr,i),
     &                           MPI_COMM_WORLD,
     &                           reqr(nrecv),
     &                           ierr)
               end if
            end if
         end do
         
         nsend = 0
         do i = 1,lnblocks
!            if (any(ref_test(:,:,:,i))) then
               ineigh = neigh(1,j,i)
               ineigh_proc = neigh(2,j,i)
               if (ineigh.ge.1) then
                  if (ineigh_proc.ne.mype) then
                     nsend = nsend + 1
                     call MPI_logical_SSEND(ref_test(1,1,1,i),
     &                              2**ndim,
     &                              MPI_LOGICAL,
     &                              ineigh_proc,
     &                              i,
     &                              MPI_COMM_WORLD,
     &                              ierr)
                  else
                     ref_testt(:,:,:,ineigh) = ref_test(:,:,:,i)
                  end if
               end if
!            end if
         end do
         
         if (nrecv.gt.0) then
            call MPI_WAITALL(nrecv,reqr,statr,ierr)
         end if

! Step 3 - if any elements of ref_testt are true, then make sure
!          that the elements of ref_test across face j are also
!          set to true. (Note: j here is the face as indexed by
!          the neighbor of block i.)
         
         do i = 1,lnblocks
            
            if (any(ref_testt(:,:,:,i))) then
               
               if (j.eq.1) then
                  
                  ix = 2
                  ix_n = 1
                  do iz = 1,1+k3d
                     do iy = 1,1+k2d
                        ref_test(ix,iy,iz,i) = 
     &                       ref_test(ix,iy,iz,i) .or. 
     &                       ref_testt(ix_n,iy,iz,i)
                     end do
                  end do
                  
               elseif (j.eq.2) then
                  
                  ix = 1
                  ix_n = 2
                  do iz = 1,1+k3d
                     do iy = 1,1+k2d
                        ref_test(ix,iy,iz,i) = 
     &                       ref_test(ix,iy,iz,i) .or. 
     &                       ref_testt(ix_n,iy,iz,i)
                     end do
                  end do
              
               elseif (j.eq.3) then
                  
                  iy = 2
                  iy_n = 1
                  do iz = 1,1+k3d
                     do ix = 1,2
                        ref_test(ix,iy,iz,i) = 
     &                       ref_test(ix,iy,iz,i) .or. 
     &                       ref_testt(ix,iy_n,iz,i)
                     end do
                  end do
                  
               elseif (j.eq.4) then
                  
                  iy = 1
                  iy_n = 2
                  do iz = 1,1+k3d
                     do ix = 1,2
                        ref_test(ix,iy,iz,i) = 
     &                       ref_test(ix,iy,iz,i) .or. 
     &                       ref_testt(ix,iy_n,iz,i)
                     end do
                  end do
                  
               elseif (j.eq.5) then
                  
                  iz = 2
                  iz_n = 1
                  do iy = 1,1+k2d
                     do ix = 1,2
                        ref_test(ix,iy,iz,i) = 
     &                       ref_test(ix,iy,iz,i) .or. 
     &                       ref_testt(ix,iy,iz_n,i)
                     end do
                  end do
                  
               elseif (j.eq.6) then
                  
                  iz = 1
                  iz_n = 2
                  do iy = 1,1+k2d
                     do ix = 1,2
                        ref_test(ix,iy,iz,i) = 
     &                       ref_test(ix,iy,iz,i) .or. 
     &                       ref_testt(ix,iy,iz_n,i)
                     end do
                  end do
                  
               end if
            end if
            
         end do
      end do

! Step 4      
! SET REFINE FLAGS BASED ON ref_test flags and repeat refinement process
! if necessary

      repeat = .FALSE.
      do i = 1,lnblocks

         do iz = 1,1+k3d
           do iy = 1,1+k2d
             do ix = 1,2
           
               if (ref_test(ix,iy,iz,i).and.nodetype(i).eq.1
     &               .and..not.refine(i)) then
                 repeat = .TRUE.
                 refine(i) = .TRUE.
                 derefine(i) = .FALSE.
               end if

             end do
           end do
         end do
               
      end do

c cycle through all processors to see if any repeat, if so all repeat
c this should be done via a scan function

      lreduce_datain(1) = repeat
      call mpi_logical_allreduce (lreduce_datain,lreduce_dataout,
     &                    1,MPI_LOGICAL,
     &                    MPI_LOR,MPI_COMM_WORLD,ierr)
      repeat_t = lreduce_dataout(1)
      repeat = repeat_t

      if (repeat) then
         go to 21
      end if

      return
      end subroutine amr_check_refine

















