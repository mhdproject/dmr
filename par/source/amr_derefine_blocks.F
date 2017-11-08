!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"


      subroutine amr_derefine_blocks(lnblocks_old)




c By K. Olson (NASA/GSFC and GMU), 11/96


      use paramesh_dimensions
      use physicaldata
      use tree

      use paramesh_interfaces, only : amr_check_derefine

      implicit none

#include "amr_shmem.fh"

      integer, intent(inout) :: lnblocks_old

      integer new_loc(maxblocks_tr),new_loct
      integer ichi
      integer i,j,k
      integer mype,lnblocks2
      integer shmem_my_pe,shmem_n_pes
      integer cempty(1)


      save new_loc,new_loct,ichi
      save k,cempty

      integer,save :: idataout(1)


      mype = shmem_my_pe()

c 1) Check if it really is OK to derefine blocks marked for derefinement

      call amr_check_derefine


c 2) remove blocks marked for derefinement by packing the data

      do i = 1,maxblocks_tr
         new_loc(i) = -1
      end do

#ifdef EMPTY_CELLS
        do i = 1,lnblocks
        if( (derefine(i)) .and. empty(i).eq.1) then
                cempty(1)=1
                empty(i) = 0
                call shmem_integer_put(empty(parent(1,i)),cempty(1),1,
     .       parent(2,i))
        endif
        enddo
        call shmem_barrier_all()
#endif

        do i = 1,lnblocks
        if( derefine(i) ) bflags(:,i) = -1
        enddo
        call shmem_barrier_all()


      k = 1
      lnblocks2 = lnblocks
      do i = 1,lnblocks
         
         if (.not.derefine(i)) then

            do j = 1,nchild
               child(1,j,k) = child(1,j,i)
               child(2,j,k) = child(2,j,i)
            end do
            parent(1,k) = parent(1,i)
            parent(2,k) = parent(2,i)
            do j = 1,nfaces
               neigh(1,j,k) = neigh(1,j,i)
               neigh(2,j,k) = neigh(2,j,i)
            end do
            do j = 1,ndim
               coord(j,k) = coord(j,i)
               bnd_box(1,j,k) = bnd_box(1,j,i)
               bnd_box(2,j,k) = bnd_box(2,j,i)
            end do
            bsize(:,k) = bsize(:,i)
            newchild(k) = newchild(i)
            which_child(k) = which_child(i)
            lrefine(k) = lrefine(i)
            bflags(:,k) = bflags(:,i)
            work_block(k) = work_block(i)
#ifdef EMPTY_CELLS
                empty(k) = empty(i)
#endif

      if(i.le.maxblocks) then
            if (nvar.gt.0) unk(:,:,:,:,k) = unk(:,:,:,:,i)
            if (nfacevar.gt.0) then
               facevarx(:,:,:,:,k) = facevarx(:,:,:,:,i)
               facevary(:,:,:,:,k) = facevary(:,:,:,:,i)
               facevarz(:,:,:,:,k) = facevarz(:,:,:,:,i)
            end if
            if (nvaredge.gt.0) then
               unk_e_x(:,:,:,:,k) = unk_e_x(:,:,:,:,i)
               unk_e_y(:,:,:,:,k) = unk_e_y(:,:,:,:,i)
               unk_e_z(:,:,:,:,k) = unk_e_z(:,:,:,:,i)
            end if
            if (nvarcorn.gt.0) unk_n(:,:,:,:,k) = unk_n(:,:,:,:,i)
        end if
                  
            
            new_loc(i) = k
            
            k = k + 1
            
         else
            
            lnblocks2 = lnblocks2 - 1
            lnblocks_old = lnblocks_old - 1
            
         end if
         
      end do
      
      call shmem_barrier_all()

c 3) overwrite old locations

      do i = lnblocks2+1,lnblocks
         
         derefine(i) = .FALSE.
         do j = 1,nchild
            child(1,j,i) = -1
            child(2,j,i) = -1
         end do
         parent(1,i) = -1
         parent(2,i) = -1
         do j = 1,nfaces
            neigh(1,j,i) = -1
            neigh(2,j,i) = -1
         end do
         do j = 1,ndim
            coord(j,i) = -1.
            bnd_box(1,j,i) = -1.
            bnd_box(2,j,i) = -1.
         end do
         bsize(:,i) = -1.
         nodetype(i) = -1
         which_child(i) = -1
         newchild(i) = .FALSE.
         lrefine(i) = -1
         bflags(:,i) = -1
         work_block(i) = 0.
#ifdef EMPTY_CELLS
         empty(i) = 0
#endif


      if(i.le.maxblocks) then
         if (nvar.gt.0) unk(:,:,:,:,i) = -1.
         if (nfacevar.gt.0) then
            facevarx(:,:,:,:,i) = -1.
            facevary(:,:,:,:,i) = -1.
            facevarz(:,:,:,:,i) = -1.
         end if
         if (nvaredge.gt.0) then
           unk_e_x(:,:,:,:,i) = -1.
           unk_e_y(:,:,:,:,i) = -1.
           unk_e_z(:,:,:,:,i) = -1.
         end if
         if (nvarcorn.gt.0) unk_n(:,:,:,:,i) = -1.
        end if
         
      end do
      
      lnblocks = lnblocks2

      call shmem_barrier_all()
 
c 4) reconnect all pointers

      do i = 1,lnblocks

         if (parent(1,i).gt.0) then
               
            ! fetch new_loc at old location of parent
            call shmem_integer_get(idataout(1),new_loc(parent(1,i)),
     $           1,parent(2,i))
            new_loct = idataout(1)
            parent(1,i) = new_loct
            
         end if
         
         do j = 1,nchild
            
            if (child(1,j,i).gt.0) then
               
               ! fetch new_loc at old location of child
               call shmem_integer_get(idataout(1),
     &               new_loc(child(1,j,i)),1
     &              ,child(2,j,i))
               new_loct = idataout(1)
               child(1,j,i) = new_loct
               
            end if

         end do
         
         do j = 1,nfaces
            
            if (neigh(1,j,i).gt.0) then
               
               ! fetch new_loc at old location of neighbor
               call shmem_integer_get(idataout(1),
     &               new_loc(neigh(1,j,i)),1
     &              ,neigh(2,j,i))
               new_loct = idataout(1)
               neigh(1,j,i) = new_loct
               
            end if
            
         end do

      end do

      call shmem_barrier_all()

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
                  call shmem_integer_get(idataout(1),
     &                 child(1,k,child(1,j,i)),
     &                 1,child(2,j,i))
                  ichi = idataout(1)
                  
                  if (ichi.le.-1) nodetype(i) = 2
                  
               end do
               
            end if
            
         end do
         
      end do

      call shmem_barrier_all()

c reset derefine flags

      do i = 1,maxblocks_tr
         derefine(i) = .FALSE.
      end do

      call shmem_barrier_all()

      return
      end subroutine amr_derefine_blocks

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine amr_check_derefine

c By K. Olson (NASA/GSFC and GMU), 4/97


      use paramesh_dimensions
      use physicaldata
      use tree

      implicit none

#include "amr_shmem.fh"

      integer ipar,ipar_proc,ichi,ichi_proc
      integer i,j,k,ineigh,ineigh_proc,childt(2),childt2(2,mchild)
      integer shmem_my_pe,shmem_n_pes
      integer cempty(1)
      integer nodetype2(maxblocks_tr)

      logical derefinet
      logical refinet(maxblocks_tr),refinet2(maxblocks_tr)
      logical refinet3(maxblocks_tr)
      logical lt

      save childt,childt2,derefinet,ichi
      save k,lt,refinet,refinet2,nodetype2,cempty

      integer,save :: idataout(1)
      logical,save :: ldataout(1)

c SAFETY checking, if the block is marked for derefinement and it is not
c a leaf block then do not derefine it !

      call shmem_barrier_all()


      do i = 1,lnblocks

         if (derefine(i).and.nodetype(i).ne.1) derefine(i) = .FALSE.

      end do

      call shmem_barrier_all()


#ifdef EMPTY_CELLS
! if an empty=1 block is marked for derefinement but it has at least
! one non-empty sibling then cancel its derefinement.
        do i = 1,lnblocks
        if(empty(i).eq.1 .and. derefine(i)) then

! get list of siblings
           call shmem_integer_get(childt2(1,1),child(1,1,parent(1,i)),
     .                                          2*mchild,parent(2,i))
! loop over siblings
                do j=1,nchild
                        call shmem_integer_get(cempty(1),
     .       empty(childt2(1,j)),1,childt2(2,j))
                        if(cempty(1).eq.0) derefine(i)=.false.
                enddo

        endif
        enddo
        call shmem_barrier_all()
#endif

c 1.1) Do a check of neighbors to check if OK to derefine

c set node types (2)

      do i = 1,lnblocks

         do j = 1,nfaces
            nodetype2(i) = 1
         end do

         do j = 1,nchild

            if (child(1,j,i).ge.1) then ! this node has children, then its
                                        ! type is 2

                  nodetype2(i) = 2
               
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
      end do

! FIRST STEP, check neighbors and set refinet

      do i = 1,lnblocks

         do j = 1,nfaces                     ! cycle through neighbors of i
            
            ineigh = neigh(1,j,i)
            ineigh_proc = neigh(2,j,i)
            
            if (ineigh.ge.1) then

               ! fetch nodetype2 of neighbor
               call shmem_integer_get (idataout(1),nodetype2(ineigh),
     $              1,ineigh_proc)
               k = idataout(1)
 
               if (k.eq.2) then
                  refinet(i) = .TRUE.
                  go to 20
               end if

            end if
               
         end do

 20      continue                            ! continue and go to next block i

      end do

      call shmem_barrier_all()

! CHECK EDGES
      if(ndim.ge.2) then
      do i = 1,lnblocks


         do j = 1,nfaces                     ! cycle through neighbors of i
            
            ineigh = neigh(1,j,i)
            ineigh_proc = neigh(2,j,i)
            
            if (ineigh.ge.1) then            ! there exits a neighbor of i
                                             ! so examine its children
               
               ! fetch refinet of neighbor
               call SHMEM_LOGICAL_GET(ldataout(1),
     &                                refinet(ineigh),1,ineigh_proc)
               lt = ldataout(1)

               if (lt) then
                  refinet2(i) = .TRUE.
               end if
               
            end if
                  
         end do

      end do

      do i = 1,lnblocks

         if (refinet2(i).or.refinet(i)) refinet2(i) = .TRUE.

      end do

      call shmem_barrier_all()
      endif

! CHECK CORNERS

!!!#if N_DIM == 3
      if (ndim == 3) then
      do i = 1,lnblocks


         do j = 1,nfaces   ! cycle through neighbors of i
            
            ineigh = neigh(1,j,i)
            ineigh_proc = neigh(2,j,i)
            
            if (ineigh.ge.1) then ! there exits a neighbor of i
                                  ! so examine its children
               
               ! fetch refinet of neighbor
               call SHMEM_LOGICAL_GET(ldataout(1),
     &                                refinet2(ineigh),1,ineigh_proc)
               lt = ldataout(1)

               if (lt) then
                  refinet3(i) = .TRUE.
               end if

            end if
                  
         end do

      end do
!!!#endif
      end if

      call shmem_barrier_all()

! SET DEREFINE FLAGS BASED ON REFINET

      do i = 1,lnblocks

         if ((refinet(i).or.refinet2(i).or.refinet3(i))
     &        .and.nodetype(i).eq.1) then
            derefine(i) = .FALSE.
         end if

      end do

      call shmem_barrier_all()

c 1.2) If a block does not have a parent (i.e. = -1) then you can^t refine
c      it further so if it is marked for derefinement turn derefine off

      do i = 1,lnblocks

         if (derefine(i).and.parent(1,i).lt.0) derefine(i) = .FALSE.

      end do

      call shmem_barrier_all()

c 1.3) check if all siblings are also marked for derefinement, if not then
c      don^t derefine this block

      do i = 1,lnblocks

         if (derefine(i)) then

            ipar = parent(1,i) ! parent of i
            ipar_proc = parent(2,i) ! processor parent is stored on

            ! fetch child ids and procno^s of ipar
            call shmem_integer_get(childt2(1,1),child(1,1,ipar),
     $           2*mchild,ipar_proc)

            do j = 1,nchild

               ichi = childt2(1,j)
               ichi_proc = childt2(2,j)

               ! fetch derefinement value of this child
               call SHMEM_LOGICAL_GET(ldataout(1),derefine(ichi),
     $              1,ichi_proc)
               derefinet = ldataout(1)

               if (.not.derefinet) derefine(i) = .FALSE.

            end do

         end if
            
      end do

      call shmem_barrier_all()


      return
      end subroutine amr_check_derefine

