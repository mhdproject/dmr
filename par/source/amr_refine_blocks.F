!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"


      subroutine amr_refine_blocks





c By K. Olson (NASA/GSFC and GMU), 11/96


      use paramesh_dimensions
      use physicaldata
      use tree

      implicit none

#include "amr_shmem.fh"


      integer nodetype2(6,maxblocks_tr)
      integer i,j,k,ipar,ipar_proc
      integer ineigh,ineigh_proc
      integer ix,iy,iz
      integer nprocs,mype,lnblocks2,ichi
      integer shmem_my_pe,shmem_n_pes
      integer nref1,nref2,nref3
      integer tneigh(2,mfaces,maxblocks_tr)
      integer,save :: i_arr(1)
      logical,save :: l_arr(1),l_arr_in(1)

      real h,hh(mdim),hx,hy,hz

      logical refinet(maxblocks_tr),
     &     refinet2(maxblocks_tr),refinet3(maxblocks_tr),
     &     repeat,repeat_t,lt

      save repeat,repeat_t,ichi,k,nodetype2,ineigh,ineigh_proc,
     &     refinet,refinet2,refinet3,lt,ipar,ipar_proc,mype,
     &     tneigh


c 1) refine blocks, create their children, the children's parents,
c    turn the children on, and turn the parents off

      mype = shmem_my_pe()
      nprocs = shmem_n_pes()

      do i = 1,maxblocks_tr
         newchild(i) = .FALSE.
      end do

c SAFETY checking, if the block is marked for refinement and it is not
c on then do not refine it !

 21   do i = 1,lnblocks

         if (refine(i).and.nodetype(i).ne.1) refine(i) = .FALSE.

      end do

      call shmem_barrier_all()

      lnblocks2 = lnblocks

      do i = 1,lnblocks

         if (refine(i)) then ! refine block 'i' on this processor

            do j = 1,nchild ! create 'nchild' child blocks

               lnblocks2 = lnblocks2 + 1
      if(lnblocks2.gt.maxblocks_tr) then
       write(*,*) 'PARAMESH ERROR !'
       write(*,*) 'Too many blocks created! '
       write(*,*) 'Increase MAXBLOCKS_TR and rerun! '
       call amr_abort
      endif
               
               child(1,j,i) = lnblocks2 ! child j's on-processor id
               child(2,j,i) = mype  ! child j's processor no.
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
                  bnd_box(2,1:ndim,lnblocks2) = coord(1:ndim,i)
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
            if (ndim.eq.3) then

            ! connect in z direction
            
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

      end do

      call shmem_barrier_all()

c 2) connect with off-processor neighbors by looking at neighbors of parent

      tneigh(:,:,:) = 0

      call shmem_barrier_all()

      do i = 1,lnblocks

         if (refine(i)) then

            ! connect left side (i.e. neighbor no. 1)
            
            ineigh = neigh(1,1,i) ! left neighbor of parent
            ineigh_proc = neigh(2,1,i) ! left neigbor's processor

            if (ineigh.gt.0) then

               do j = 1,nchild-1,2 ! cycle through children of 'i' on left side

                  k = j+1 ! child no. of neighbor which lies on border

                  ichi = child(1,j,i) ! child 'j' of parent 'i'

                  ! fetch child no. k of ineigh from ineigh_proc
                  ! and store in neighbor 1 of child j (ichi)

                  call shmem_integer_get (neigh(1,1,ichi),
     .       child(1,k,ineigh),2,ineigh_proc)

                  ! write yourself to other neighbor

                  if (neigh(1,1,ichi).gt.0) then

                     i_arr(1) = ichi
                     call shmem_integer_put 
     &                    (tneigh(1,2,neigh(1,1,ichi)),
     &                    i_arr(1),1,neigh(2,1,ichi))
!     &                    ichi,1,neigh(2,1,ichi))
                     i_arr(1) = mype
                     call shmem_integer_put 
     &                    (tneigh(2,2,neigh(1,1,ichi)),
     &                    i_arr(1),1,neigh(2,1,ichi))
!     &                    mype,1,neigh(2,1,ichi))

                  end if

               end do

            end if

            ! connect right side (i.e. neighbor no. 2)

            ineigh = neigh(1,2,i) ! right neighbor of parent
            ineigh_proc = neigh(2,2,i) ! right neighbor's processor

            if (ineigh.gt.0) then

               do j = 2,nchild,2 ! cycle through children of 'i' on left side

                  k = j-1 ! child no. of neighbor which lies on border

                  ichi = child(1,j,i) ! child 'j' of parent 'i'

                  ! fetch child no. k data (processor and id) of ineigh 
                  ! from ineigh_proc and store in neighbor 2 of child j (ichi)

                  call shmem_integer_get (neigh(1,2,ichi),
     .       child(1,k,ineigh),2,ineigh_proc)

                  ! write yourself to other neignbor

                  if (neigh(1,2,ichi).gt.0) then

                     i_arr(1) = ichi
                     call shmem_integer_put 
     &                    (tneigh(1,1,neigh(1,2,ichi)),
     &                    i_arr(1),1,neigh(2,2,ichi))
!     &                    ichi,1,neigh(2,2,ichi))
                     i_arr(1) = mype
                     call shmem_integer_put 
     &                    (tneigh(2,1,neigh(1,2,ichi)),
     &                    i_arr(1),1,neigh(2,2,ichi))
!     &                    mype,1,neigh(2,2,ichi))

                  end if

               end do

            end if

           ! connect in -y direction (i.e. neighbor no. 3)
           if (ndim.ge.2) then

            ineigh = neigh(1,3,i) ! bottom neighbor of parent
            ineigh_proc = neigh(2,3,i) ! bottom neigbor's processor

            if (ineigh.gt.0) then

               do j = 1,2 ! cycle through children of 'i' on bottom

                  k = j+2 ! child no. of neighbor which lies on border

                  ichi = child(1,j,i) ! child 'j' of parent 'i'

                  ! fetch child no. k of ineigh from ineigh_proc
                  ! and store in neighbor 3 of child j (ichi)

                  call shmem_integer_get (neigh(1,3,ichi),
     .                child(1,k,ineigh),2,ineigh_proc)

                  ! write yourself to other neignbor

                  if (neigh(1,3,ichi).gt.0) then

                     i_arr(1) = ichi
                     call shmem_integer_put 
     &                    (tneigh(1,4,neigh(1,3,ichi)),
     &                    i_arr(1),1,neigh(2,3,ichi))
!     &                    ichi,1,neigh(2,3,ichi))
                     i_arr(1) = mype
                     call shmem_integer_put 
     &                    (tneigh(2,4,neigh(1,3,ichi)),
     &                    i_arr(1),1,neigh(2,3,ichi))
!     &                    mype,1,neigh(2,3,ichi))

                  end if

               end do

               if (ndim.eq.3) then

               do j = 5,6 ! cycle through children of 'i' on bottom

                  k = j+2 ! child no. of neighbor which lies on border

                  ichi = child(1,j,i) ! child 'j' of parent 'i'

                  ! fetch child no. k of ineigh from ineigh_proc
                  ! and store in neighbor 3 of child j (ichi)

                  call shmem_integer_get (neigh(1,3,ichi),
     .                 child(1,k,ineigh),2,ineigh_proc)

                  ! write yourself to other neignbor


                  if (neigh(1,3,ichi).gt.0) then

                     i_arr(1) = ichi
                     call shmem_integer_put 
     &                    (tneigh(1,4,neigh(1,3,ichi)),
     &                    i_arr(1),1,neigh(2,3,ichi))
!     &                    ichi,1,neigh(2,3,ichi))
                     i_arr(1) = mype
                     call shmem_integer_put 
     &                    (tneigh(2,4,neigh(1,3,ichi)),
     &                    i_arr(1),1,neigh(2,3,ichi))
!     &                    mype,1,neigh(2,3,ichi))

                 end if

               end do

               end if

            end if

            ! connect in +y direction (i.e. neighbor no. 4)

            ineigh = neigh(1,4,i) ! bottom neighbor of parent
            ineigh_proc = neigh(2,4,i) ! bottom neigbor's processor

            if (ineigh.gt.0) then

               do j = 3,4 ! cycle through children of 'i' on bottom

                  k = j-2 ! child no. of neighbor which lies on border

                  ichi = child(1,j,i) ! child 'j' of parent 'i'

                  ! fetch child no. k of ineigh from ineigh_proc
                  ! and store in neighbor 4 of child j (ichi)

                  call shmem_integer_get (neigh(1,4,ichi),
     .       child(1,k,ineigh),2,ineigh_proc)

                  ! write yourself to other neignbor

                  if (neigh(1,4,ichi).gt.0) then

                     i_arr(1) = ichi
                     call shmem_integer_put 
     &                    (tneigh(1,3,neigh(1,4,ichi)),
     &                    i_arr(1),1,neigh(2,4,ichi))
!     &                    ichi,1,neigh(2,4,ichi))
                     i_arr(1) = mype
                     call shmem_integer_put 
     &                    (tneigh(2,3,neigh(1,4,ichi)),
     &                    i_arr(1),1,neigh(2,4,ichi))
!     &                    mype,1,neigh(2,4,ichi))

                  end if

               end do
!!!#if N_DIM == 3
               if (ndim.eq.3) then

               do j = 7,8 ! cycle through children of 'i' on bottom

                  k = j-2 ! child no. of neighbor which lies on border

                  ichi = child(1,j,i) ! child 'j' of parent 'i'

                  ! fetch child no. k of ineigh from ineigh_proc
                  ! and store in neighbor 4 of child j (ichi)

                  call shmem_integer_get (neigh(1,4,ichi),
     .       child(1,k,ineigh),2,ineigh_proc)

                  ! write yourself to other neignbor

                  if (neigh(1,4,ichi).gt.0) then

                     i_arr(1) = ichi
                     call shmem_integer_put 
     &                    (tneigh(1,3,neigh(1,4,ichi)),
     &                    i_arr(1),1,neigh(2,4,ichi))
!     &                    ichi,1,neigh(2,4,ichi))
                     i_arr(1) = mype
                     call shmem_integer_put 
     &                    (tneigh(2,3,neigh(1,4,ichi)),
     &                    i_arr(1),1,neigh(2,4,ichi))
!     &                    mype,1,neigh(2,4,ichi))

                  end if

               end do

               end if
!!!#endif

            end if
           end if

!!!#if N_DIM == 3
            if (ndim.eq.3) then

            ! connect -z (i.e. neighbor no. 5)

            ineigh = neigh(1,5,i) ! right neighbor of parent
            ineigh_proc = neigh(2,5,i) ! right neigbor's processor

            if (ineigh.gt.0) then

               do j = 1,4 ! cycle through children of 'i' on left side

                  k = j+4 ! child no. of neighbor which lies on border

                  ichi = child(1,j,i) ! child 'j' of parent 'i'

                  ! fetch child no. k of ineigh from ineigh_proc
                  ! and store in neighbor 5 of child j (ich)

                  call shmem_integer_get (neigh(1,5,ichi),
     .       child(1,k,ineigh),2,ineigh_proc)

                  ! write yourself to other neignbor


                  if (neigh(1,5,ichi).gt.0) then

                     i_arr(1) = ichi
                     call shmem_integer_put 
     &                    (tneigh(1,6,neigh(1,5,ichi)),
     &                    i_arr(1),1,neigh(2,5,ichi))
!     &                    ichi,1,neigh(2,5,ichi))
                     i_arr(1) = mype
                     call shmem_integer_put 
     &                    (tneigh(2,6,neigh(1,5,ichi)),
     &                    i_arr(1),1,neigh(2,5,ichi))
!     &                    mype,1,neigh(2,5,ichi))

                  end if

               end do

            end if

            ! connect +z (i.e. neighbor no. 6)

            ineigh = neigh(1,6,i) ! right neighbor of parent
            ineigh_proc = neigh(2,6,i) ! right neigbor's processor

            if (ineigh.gt.0) then

               do j = 5,8 ! cycle through children of 'i' on left side

                  k = j-4 ! child no. of neighbor which lies on border

                  ichi = child(1,j,i) ! child 'j' of parent 'i'

                  ! fetch child no. k of ineigh from ineigh_proc
                  ! and store in neighbor 6 of child j (ich)

                  call shmem_integer_get (neigh(1,6,ichi),
     .       child(1,k,ineigh),2,ineigh_proc)

                  ! write yourself to other neignbor

                  if (neigh(1,6,ichi).gt.0) then

                     i_arr(1) = ichi
                     call shmem_integer_put 
     &                    (tneigh(1,5,neigh(1,6,ichi)),
     &                    i_arr(1),1,neigh(2,6,ichi))
!     &                    ichi,1,neigh(2,6,ichi))
                     i_arr(1) = mype
                     call shmem_integer_put 
     &                    (tneigh(2,5,neigh(1,6,ichi)),
     &                    i_arr(1),1,neigh(2,6,ichi))
!     &                    mype,1,neigh(2,6,ichi))

                  end if

               end do

            end if

            end if
!!!#endif /*N_DIM*/

         end if

      end do

      call shmem_barrier_all()

      do i = 1,lnblocks2
         do j = 1,nfaces
            if (tneigh(1,j,i).gt.0.and.
     &           tneigh(1,j,i).ne.neigh(1,j,i)) then
               neigh(1,j,i) = tneigh(1,j,i)
               neigh(2,j,i) = tneigh(2,j,i)
            end if
         end do
      end do
      
      call shmem_barrier_all()

      lnblocks = lnblocks2

      call shmem_barrier_all()

c reset node types (2)
!
! nodetype2 is set for each face on each block. It stores information on the
! maximum refinement level to be found immediately interior to that face of the
! block. If a block is a leaf block then nodetype2 is 1 for all faces. If it is
! a parent block then nodetype2 is 2. If it is a grandparent than for any faces 
! which also bound a grand-child, nodetype2 is set to 3, otherwise they are set to 2.
! By testing the value of nodetype2 on a neighbor block, a block can decide if
! it needs to refine in order to prevent a jump in refinement level of more than
! 1 level from developing.
!
      do i = 1,maxblocks_tr
         do j = 1,nfaces
            nodetype2(j,i) = 1
         end do
      end do

      do i = 1,lnblocks

         do j = 1,nfaces
            nodetype2(j,i) = 1
         end do

         ! 1) cycle through children of i

         do j = 1,nchild

            if (child(1,j,i).ge.1) then

               do k = 1,nfaces
                  if (nodetype2(k,i).ne.3) nodetype2(k,i) = 2
               end do
               
               ! 2) cycle through children of child(1,j,i)

               do k = 1,nchild
                  
                  ! fetch child of child
                 call shmem_integer_get(i_arr(1),child(1,k,child(1,j,i))
     &                 ,1,child(2,j,i))
                  ichi = i_arr(1)

                  if (ichi.ge.1) then

                     if (j.eq.1) then
                        nodetype2(2,i) = 3
                        if (ndim.ge.2) nodetype2(4,i) = 3
                        if (ndim.eq.3) nodetype2(6,i) = 3
                     else if (j.eq.2) then
                        nodetype2(1,i) = 3
                        if (ndim.ge.2) nodetype2(4,i) = 3
                        if (ndim.eq.3) nodetype2(6,i) = 3
                     else if (j.eq.3) then
                        nodetype2(2,i) = 3
                        if (ndim.ge.2) nodetype2(3,i) = 3
                        if (ndim.eq.3) nodetype2(6,i) = 3
                     else if (j.eq.4) then
                        nodetype2(1,i) = 3
                        if (ndim.ge.2) nodetype2(3,i) = 3
                        if (ndim.eq.3) nodetype2(6,i) = 3
                     else if (j.eq.5) then
                        nodetype2(2,i) = 3
                        if (ndim.ge.2) nodetype2(4,i) = 3
                        if (ndim.eq.3) nodetype2(5,i) = 3
                     else if (j.eq.6) then
                        nodetype2(1,i) = 3
                        if (ndim.ge.2) nodetype2(4,i) = 3
                        if (ndim.eq.3) nodetype2(5,i) = 3
                     else if (j.eq.7) then
                        nodetype2(2,i) = 3
                        if (ndim.ge.2) nodetype2(3,i) = 3
                        if (ndim.eq.3) nodetype2(5,i) = 3
                     else if (j.eq.8) then
                        nodetype2(1,i) = 3
                        if (ndim.ge.2) nodetype2(3,i) = 3
                        if (ndim.eq.3) nodetype2(5,i) = 3
                     end if

                  end if

               end do
               
            end if
            
         end do

      end do

c reset node types

      do i = 1,lnblocks

         nodetype(i) = 3
         
         ! 1) cycle through children of i
         
         do j = 1,nchild
            
            if (child(1,j,i).le.-1) then
                  
               nodetype(i) = 1
                  
            else
                  
               ! 2) cycle through children of child(1,j,i)
               
               do k = 1,nchild
                  
                  ! fetch child of child
                 call shmem_integer_get(i_arr(1),child(1,k,child(1,j,i))
     &                 ,1,child(2,j,i))
                  ichi = i_arr(1)
                  
                  if (ichi.le.-1) nodetype(i) = 2
                  
               end do
               
            end if
            
         end do
         
      end do

      call shmem_barrier_all()

c Check for neighboring blocks which are more than one level of refinement
c different

      do i = 1,maxblocks_tr
         refinet(i) = .FALSE.
         refinet2(i) = .FALSE.
         refinet3(i) = .FALSE.
         refine(i) = .FALSE.
      end do
      repeat = .FALSE.

! FIRST STEP, check neighbors and set refinet

      nref1 = 0
      nref2 = 0
      nref3 = 0

      do i = 1,lnblocks

         do j = 1,nfaces     ! cycle through neighbors of i in i's x-y plane
            
            ineigh = neigh(1,j,i)
            ineigh_proc = neigh(2,j,i)
            
            if (ineigh.ge.1) then

               ! fetch nodetype2 of neighbor
               call shmem_integer_get (i_arr(1),nodetype2(j,ineigh),1,
     .       ineigh_proc)
               k = i_arr(1)

               
               if (k.eq.3) then
                  refinet(i) = .TRUE.
                  nref1 = nref1 + 1
                  go to 20
               end if

            end if
               
         end do

 20      continue ! continue and go to next block 'i'

      end do

      call shmem_barrier_all()

! CHECK EDGES

      if (ndim.ge.2) then

      do i = 1,lnblocks

         ix = 0
         iy = 0
         iz = 0

         do j = 1,nfaces  ! cycle through neighbors of i in i's x-y plane
            
            ineigh = neigh(1,j,i)
            ineigh_proc = neigh(2,j,i)
            
            if (ineigh.ge.1) then ! there exits a neighbor of i
                                     ! so examine its children
               
                                ! fetch refinet of neighbor
               
               call shmem_logical_get (l_arr(1),refinet(ineigh),
     &                                             1,ineigh_proc)
               lt = l_arr(1)

               if (lt.and.(j.eq.1.or.j.eq.2)) then
                  ix = ix + 1
               else if (lt.and.(j.eq.3.or.j.eq.4)) then
                  iy = iy + 1
               else if (lt.and.(j.eq.5.or.j.eq.6)) then
                  iz = iz + 1
               end if
               
            end if
                  
         end do

         if ((ix.ge.1.and.iy.ge.1).or.(ix.ge.1.and.iz.ge.1).or.
     &        (iy.ge.1.and.iz.ge.1)) then
            refinet2(i) = .TRUE.
            nref2 = nref2 + 1
         end if

      end do

      do i = 1,lnblocks

         if (refinet2(i).or.refinet(i)) refinet2(i) = .TRUE.

      end do

      endif

      call shmem_barrier_all()

! CHECK CORNERS

      if (ndim.eq.3) then

      do i = 1,lnblocks

         ix = 0
         iy = 0
         iz = 0

         do j = 1,nfaces   ! cycle through neighbors of i in i's x-y plane
            
            ineigh = neigh(1,j,i)
            ineigh_proc = neigh(2,j,i)
            
            if (ineigh.ge.1) then ! there exits a neighbor of i
                                     ! so examine its children
               
                                ! fetch refinet of neighbor
               
              call shmem_logical_get (l_arr(1),refinet2(ineigh),
     &                                             1,ineigh_proc)
              lt = l_arr(1)

               if (lt.and.(j.eq.1.or.j.eq.2)) then
                  ix = ix + 1
               else if (lt.and.(j.eq.3.or.j.eq.4)) then
                  iy = iy + 1
               else if (lt.and.(j.eq.5.or.j.eq.6)) then
                  iz = iz + 1
               end if

            end if
                  
         end do

         if (ix.ge.1.and.iy.ge.1.and.iz.ge.1) then
            refinet3(i) = .TRUE.
            nref3 = nref3 + 1
         end if

      end do

      end if

      call shmem_barrier_all()

! SET REFINE FLAGS BASED ON REFINET

      do i = 1,lnblocks

         refine(i) = .FALSE.

         if ((refinet(i).or.refinet2(i).or.refinet3(i))
     &        .and.nodetype(i).eq.1) then
            repeat = .TRUE.
            refine(i) = .TRUE.
         end if

      end do

      call shmem_barrier_all()

c cycle through all processors to see if any repeat, if so all repeat
c this should be done via a scan function

      ix = mype
      do i = 0,nprocs-1

         l_arr_in(1) = repeat
         call shmem_logical_get(l_arr(1),l_arr_in(1),1,i)
         repeat_t = l_arr(1)
         if (repeat_t) repeat = .TRUE.
         ix = ix + 1
         if (ix.ge.nprocs-1) ix = 0

      end do

      if (repeat) then
         go to 21
      end if

      call shmem_barrier_all()

c Now set neighbor pointers of new children if they lie on a boundry

      do i = 1,lnblocks

         if (newchild(i)) then

            ! fetch i's parent

            ipar = parent(1,i)
            ipar_proc = parent(2,i)

            do j = 1,nfaces

               if (neigh(1,j,i).le.-1) then ! this neighbor may be on a border

                  ! fetch i's parent's neighbor j

                  call shmem_integer_get (i_arr(1),neigh(1,j,ipar),1,
     .                                    ipar_proc)
                  ineigh = i_arr(1)
               
                  ! if parent's neighbor is lt -1 then i's neighbor is
                  ! also on the domain border and is set to the parent's
                  ! value

                  if (ineigh.le.-20) neigh(1,j,i) = ineigh

               end if

            end do

         end if

      end do

      call shmem_barrier_all()

c reset refine flags

      do i = 1,maxblocks_tr
         refine(i) = .FALSE.
      end do

      call shmem_barrier_all


      return
      end subroutine amr_refine_blocks
