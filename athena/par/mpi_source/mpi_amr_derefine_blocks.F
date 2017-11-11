!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"


      subroutine amr_derefine_blocks(lnblocks_old,mype)




c By K. Olson (NASA/GSFC and GMU), 11/96

      use paramesh_dimensions
      use physicaldata
      use tree

      implicit none

      integer, intent(inout) :: lnblocks_old
      integer, intent(in)    :: mype

      include 'mpif.h'

! local variables --------------------------------------------------------

      integer :: new_loc(maxblocks_tr)
      integer :: i,j,k,jsend
      integer :: lnblocks2
      integer :: neight(2,mfaces,maxblocks_tr)
      integer :: childt(2,mchild,maxblocks_tr)
      integer :: parentt(2,maxblocks_tr)
      integer :: statr(MPI_STATUS_SIZE,maxblocks_tr)
      integer :: reqr(maxblocks_tr)
      integer :: ierr,nsend,nrecv
      integer :: nodetype_chi(nchild,maxblocks_tr)

! ------------------------------------------------------------------------

c remove blocks marked for derefinement by packing the data

      do i = 1,maxblocks_tr
         new_loc(i) = -1
      end do

! Compute new_loc, new_loc marks where each block will end up after the
! derefinement is done

      k = 1
      do i = 1,lnblocks
         if (.not.derefine(i)) then
            new_loc(i) = k
            k = k + 1
          endif
      end do

c 4) reconnect all pointers

      parentt(:,1:lnblocks) = parent(:,1:lnblocks)
      childt(:,:,1:lnblocks) = child(:,:,1:lnblocks)
      neight(:,:,1:lnblocks) = neigh(:,:,1:lnblocks)

      nrecv = 0
      do i = 1,lnblocks
         if (parent(1,i).gt.0) then
           if (parent(2,i).ne.mype) then
             nrecv = nrecv + 1
             call MPI_int_IRECV(parentt(1,i),1,MPI_INTEGER,
     &            parent(2,i),i,MPI_COMM_WORLD,
     &            reqr(nrecv),ierr)
           else
             parentt(1,i) = new_loc(parent(1,i))
           end if
         end if
       end do
       
       nsend = 0
       do i = 1,lnblocks
         do j = 1,nchild
           if (child(1,j,i).gt.0) then
             if (child(2,j,i).ne.mype) then
               ! parent is sending to all its children
               nsend = nsend + 1
               call MPI_int_SSEND (new_loc(i),1,MPI_INTEGER,
     &              child(2,j,i),child(1,j,i),MPI_COMM_WORLD,
     &              ierr)
             end if
           end if
         end do
       end do

      if (nrecv.gt.0) then
        call MPI_WAITALL(nrecv,reqr,statr,ierr)
      end if

      nrecv = 0
      do i = 1,lnblocks
        do j = 1,nchild
          if (child(1,j,i).gt.0) then
            if (child(2,j,i).ne.mype) then
              nrecv = nrecv + 1
              call MPI_int_IRECV(childt(1,j,i),1,MPI_INTEGER,
     &             child(2,j,i),child(1,j,i),MPI_COMM_WORLD,
     &             reqr(nrecv),ierr)
            else
              childt(1,j,i) = new_loc(child(1,j,i))
            end if
          end if
        end do
       end do
       
       nsend = 0
       do i = 1,lnblocks
         if (parent(1,i).gt.0) then
           if (parent(2,i).ne.mype) then
! child is sending to its parent
             nsend = nsend + 1
             call MPI_int_SSEND (new_loc(i),1,MPI_INTEGER,
     &            parent(2,i),i,MPI_COMM_WORLD,
     &            ierr)
           end if
         end if
       end do

      if (nrecv.gt.0) then
        call MPI_WAITALL(nrecv,reqr,statr,ierr)
      end if

      do j = 1,nfaces

         if (mod(j,2).eq.0) then
            jsend = j - 1
         else
            jsend = j + 1
         end if
            
         nrecv = 0
         do i = 1,lnblocks
            if (neigh(1,j,i).gt.0) then
               if (neigh(2,j,i).ne.mype) then
                  nrecv = nrecv + 1
                  call MPI_int_IRECV(neight(1,j,i),1,MPI_INTEGER,
     &                 neigh(2,j,i),neigh(1,j,i),MPI_COMM_WORLD,
     &                 reqr(nrecv),ierr)
               else
                  neight(1,j,i) = new_loc(neigh(1,j,i))
               end if
            end if
         end do
      
         nsend = 0
         do i = 1,lnblocks
            if (neigh(1,jsend,i).gt.0) then
               if (neigh(2,jsend,i).ne.mype) then
                  nsend = nsend + 1
                  call MPI_int_SSEND (new_loc(i),1,MPI_INTEGER,
     &                 neigh(2,jsend,i),i,MPI_COMM_WORLD,
     &                 ierr)
               end if
            end if
         end do

         if (nrecv.gt.0) then
            call MPI_WAITALL(nrecv,reqr,statr,ierr)
         end if

      end do

      do i = 1,lnblocks_old
        if (new_loc(i).ne.i.and.new_loc(i).gt.0) then
          if(nvar.gt.0) unk(:,:,:,:,new_loc(i)) = unk(:,:,:,:,i)
          if (nfacevar.gt.0) then
             facevarx(:,:,:,:,new_loc(i)) = facevarx(:,:,:,:,i)
             facevary(:,:,:,:,new_loc(i)) = facevary(:,:,:,:,i)
             facevarz(:,:,:,:,new_loc(i)) = facevarz(:,:,:,:,i)
          end if
          if (nvaredge.gt.0) then
             unk_e_x(:,:,:,:,new_loc(i)) = unk_e_x(:,:,:,:,i)
             unk_e_y(:,:,:,:,new_loc(i)) = unk_e_y(:,:,:,:,i)
             unk_e_z(:,:,:,:,new_loc(i)) = unk_e_z(:,:,:,:,i)
          end if
          if (nvarcorn.gt.0) unk_n(:,:,:,:,new_loc(i)) = 
     .                                       unk_n(:,:,:,:,i)
        end if
      end do

      parent(1,1:lnblocks) = parentt(1,1:lnblocks)
      child(1,:,1:lnblocks) = childt(1,:,1:lnblocks)
      neigh(1,:,1:lnblocks) = neight(1,:,1:lnblocks)

      k = 1
      lnblocks2 = lnblocks
      do i = 1,lnblocks
         
         if (.not.derefine(i)) then
            
            if (k.ne.i) then
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

               
            end if

            k = k + 1
            
         else
            
            lnblocks2 = lnblocks2 - 1
            lnblocks_old = lnblocks_old - 1
            
         end if
         
      end do

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

         
      end do
      
      lnblocks = lnblocks2

c reset node types

      do i = 1,lnblocks
         nodetype(i) = 3
         if (child(1,1,i).le.-1) then
            nodetype(i) = 1
         end if
      end do
      nrecv = 0
      do i = 1,lnblocks
         do j = 1,nchild
            nodetype_chi(j,i) = -1 
            if (child(1,j,i).gt.-1) then
            if (child(2,j,i).ne.mype) then
               nrecv = nrecv + 1
               call MPI_int_IRECV(nodetype_chi(j,i),
     .                        1,
     .                        MPI_INTEGER,
     .                        child(2,j,i),
     .                        child(1,j,i),
     .                        MPI_COMM_WORLD,
     .                        reqr(nrecv),
     .                        ierr)
            else
               nodetype_chi(j,i) = nodetype(child(1,j,i))
            end if
            end if
         end do
      end do

      nsend = 0
      do i = 1,lnblocks
! send nodetype to your parent
         if (parent(1,i).ge.1) then
         if (parent(2,i).ne.mype) then
            nsend = nsend + 1
 ! parent(2,i) - PE TO SEND TO
 ! i           - THIS IS THE TAG
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
         call MPI_WAITALL (nrecv, reqr, statr, ierr)
      end if

      do i = 1,lnblocks
         do j = 1,nchild
            if (nodetype_chi(j,i).eq.1) nodetype(i) = 2
         end do
      end do

c reset derefine flags

      do i = 1,maxblocks_tr
         derefine(i) = .FALSE.
      end do



      return
      end subroutine amr_derefine_blocks

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine amr_check_derefine (mype)

c By K. Olson (NASA/GSFC and GMU), 4/97

      use paramesh_dimensions
      use physicaldata
      use tree
      use timings

      implicit none

      integer, intent(in) :: mype

      include 'mpif.h'

! local variables ----------------------------------------------------------

      integer :: i,j
      integer :: nodetype2(maxblocks_tr)
      integer :: nodetype_recv(maxblocks_tr)
      integer :: nodetype_send(maxblocks_tr)
      integer :: ipar,ipar_proc
      integer :: reqr(maxblocks)
      integer :: statr(MPI_STATUS_SIZE,maxblocks)
      integer :: isg,ierr,neighs,neighr,jsend

      logical :: derefine_chi(nchild,maxblocks_tr)
      logical :: refine_par(maxblocks_tr)

      double precision :: time1

! --------------------------------------------------------------------------

      no_of_calls_check_derefine = no_of_calls_check_derefine + 1

c If the block is marked for derefinement and it is not
c a leaf block then do not derefine it

#ifdef TIMING_MPI 
      time1 = mpi_wtime()
#endif /* TIMING_MPI */

      do i = 1,lnblocks
         if (derefine(i).and.nodetype(i).ne.1) derefine(i) = .FALSE.
      end do

c Allow the child blocks to only derefine if their parent is NOT
c marked for refinement.

      neighr = 0
      do i = 1,lnblocks
         refine_par(i) = .FALSE.
         if (parent(1,i).gt.0) then
            if (parent(2,i).ne.mype) then
               neighr = neighr + 1
               call MPI_logical_IRECV(refine_par(i),
     .                        1,
     .                        MPI_LOGICAL,
     .                        parent(2,i),
     .                        i,
     .                        MPI_COMM_WORLD,
     .                        reqr(neighr),
     .                        ierr)
            else
               refine_par(i) = refine(parent(1,i))
            end if
         end if
      end do

      neighs = 0
      do i = 1,lnblocks
         do j = 1,nchild
          if(child(1,j,i).gt.0) then
             if (child(2,j,i).ne.mype) then
                neighs = neighs + 1
! child(2,j,i) - PE TO SEND TO                
! child(1,j,i) - THIS IS THE TAG
                call MPI_logical_SSEND(refine(i),
     &                         1,
     &                         MPI_LOGICAL,
     &                         child(2,j,i),
     &                         child(1,j,i),
     &                         MPI_COMM_WORLD,
     &                         ierr)
            end if
          end if
         end do
      end do

      if (neighr.gt.0) then
         call MPI_WAITALL (neighr, reqr, statr, ierr)
      end if

#ifdef TIMING_MPI 
      timer_amr_check_derefine(1) = timer_amr_check_derefine(1)
     .                + mpi_wtime() - time1
      time1 = mpi_wtime()
#endif /* TIMING_MPI */

      do i = 1,lnblocks
         if(nodetype(i).eq.1.and.derefine(i)) then
            if(refine_par(i)) derefine(i)=.false.
         endif
      enddo

c Turn off refine flags of non-leaf blocks

      do i = 1,lnblocks
         if (nodetype(i).gt.1.and.refine(i)) refine(i) = .false.
      end do

c Check neighbors to check if OK to derefine

c set nodetype2 = 2 if it either has children or it is marked for
c refinement

      do i = 1,lnblocks
        nodetype2(i) = 1
        if (child(1,1,i).ge.1.or.refine(i)) then ! this node has children 
                                                 ! or it is marked for 
                                                 ! refinement then its
                                                 ! type is 2
          nodetype2(i) = 2
        end if
      end do

c Check for neighboring blocks which are more than one level of refinement
c different

! cycle through block faces
      
      do j = 1,nfaces

         if (j.eq.1) jsend = 2
         if (j.eq.2) jsend = 1
         if (j.eq.3) jsend = 4
         if (j.eq.4) jsend = 3
         if (j.eq.5) jsend = 6
         if (j.eq.6) jsend = 5
         
         neighr = 0
         do isg = 1,lnblocks
            nodetype_recv(isg) = 0
            if(neigh(1,j,isg).gt.-1) then
               if(neigh(2,j,isg).ne.mype) then
                  neighr = neighr + 1
                  call MPI_int_IRECV(nodetype_recv(isg),
     .                           1,
     .                           MPI_INTEGER,
     .                           neigh(2,j,isg),
     .                           neigh(1,j,isg),
     .                           MPI_COMM_WORLD,
     .                           reqr(neighr),
     .                           ierr)
                  mess_counter_chk_deref = mess_counter_chk_deref + 1
               else
                  nodetype_recv(isg) = nodetype2(neigh(1,j,isg))
               end if
            end if
         end do

! send nodetype2 to neigh if neighbor is off processor and nodetype2 = 2

         neighs = 0
         do isg = 1,lnblocks
!            if (nodetype2(isg).eq.2) then
               if(neigh(1,jsend,isg).gt.-1) then
                  if(neigh(2,jsend,isg).ne.mype) then
                     neighs = neighs + 1
! neigh(2,jsend,isg) - PE TO SEND TO
! isg                - THIS IS THE TAG
                     call MPI_int_SSEND(nodetype2(isg),
     &                              1,
     &                              MPI_INTEGER,
     &                              neigh(2,jsend,isg),
     &                              isg,
     &                              MPI_COMM_WORLD,
     &                              ierr)
                  end if
               end if
!            end if
         end do
         
         if (neighr.gt.0) then
            call MPI_WAITALL (neighr, reqr, statr, ierr)
         end if
      
#ifdef TIMING_MPI 
      timer_amr_check_derefine(2) = timer_amr_check_derefine(2)
     .                + mpi_wtime() - time1
      time1 = mpi_wtime()
#endif /* TIMING_MPI */

         do i = 1,lnblocks
            if (nodetype_recv(i).eq.2) nodetype2(i) = 2
         end do

      end do

! Now reset derefine flags based on value of nodetype2

      do i = 1,lnblocks

         if (nodetype2(i).eq.2 .and. derefine(i)) then
            derefine(i) = .FALSE.
         end if
         
      end do

c 1.2) If a block does not have a parent (i.e. = -1) then you can^t derefine
c      it further so if it is marked for derefinement turn derefine off

      do i = 1,lnblocks

         if (derefine(i).and.parent(1,i).lt.0) derefine(i) = .FALSE.

      end do

c 1.3) check if all siblings are also marked for derefinement, if not then
c      don^t derefine this block

! parents collect messages from children and count the number of children
! marked for derefinement (stored in nodetype_recv).

      neighr = 0
      do isg = 1,lnblocks
         do j = 1,nchild
            derefine_chi(j,isg) = .FALSE.
            if (child(1,j,isg).gt.-1) then
            if (child(2,j,isg).ne.mype) then
               neighr = neighr + 1
! derefine_chi(j,isg) - this is just junk
               call MPI_logical_IRECV(derefine_chi(j,isg),
     .                        1,
     .                        MPI_LOGICAL,
     .                        child(2,j,isg),
     .                        child(1,j,isg),
     .                        MPI_COMM_WORLD,
     .                        reqr(neighr),
     .                        ierr)
            else
               derefine_chi(j,isg) = derefine(child(1,j,isg))
            end if
            end if
         end do
      end do

! Children send a message to parent if they ar marked for derefinement

      neighs = 0
      nodetype_recv(:) = 0    ! using this variable as a counter here

      do i = 1,lnblocks
!         if (derefine(i)) then
            ipar = parent(1,i) ! parent of i
            ipar_proc = parent(2,i) ! processor parent is stored on
            if (ipar.gt.-1) then
            if (ipar_proc.ne.mype) then
               neighs = neighs + 1
! ipar_proc - PE TO SEND TO
! i         - THIS IS THE TAG
               call MPI_logical_SSEND(derefine(i),
     &                        1,
     &                        MPI_LOGICAL,
     &                        ipar_proc,
     &                        i,
     &                        MPI_COMM_WORLD,
     &                        ierr)
            end if
            end if
!         end if
      end do

      if (neighr.gt.0) then
         call MPI_WAITALL (neighr, reqr, statr, ierr)
      end if
      
#ifdef TIMING_MPI 
      timer_amr_check_derefine(3) = timer_amr_check_derefine(3)
     .                + mpi_wtime() - time1
      time1 = mpi_wtime()
#endif /* TIMING_MPI */

      do i = 1,lnblocks
         do j = 1,nchild
            if (derefine_chi(j,i)) then
               nodetype_recv(i) = nodetype_recv(i) + 1
            end if
         end do
      end do
      nodetype_send(1:lnblocks) = nodetype_recv(1:lnblocks)

! Now parent sends nodetype_recv to its children if nodetype_recv = nchild

! child blocks post recieves      

      neighr = 0
      do isg = 1,lnblocks
         if(parent(1,isg).gt.-1) then
            if(parent(2,isg).ne.mype) then
               neighr = neighr + 1
               call MPI_int_IRECV(nodetype_recv(isg),
     .                        1,
     .                        MPI_INTEGER,
     .                        parent(2,isg),
     .                        isg,
     .                        MPI_COMM_WORLD,
     .                        reqr(neighr),
     .                        ierr)
            end if
         end if
      end do

      neighs = 0
      do isg = 1,lnblocks
!         if (nodetype_recv(isg).eq.nchild) then
            do j = 1,nchild
               if (child(1,j,isg).ge.1) then
                  if (child(2,j,isg).ne.mype) then
                     neighs = neighs + 1
! child(2,j,isg) - PE TO SEND TO
! child(1,j,isg) - THIS IS THE TAG
                     call MPI_int_SSEND(nodetype_send(isg),
     &                              1,
     &                              MPI_INTEGER,
     &                              child(2,j,isg),
     &                              child(1,j,isg),
     &                              MPI_COMM_WORLD,
     &                              ierr)
                  else
                     nodetype_recv(child(1,j,isg)) = nodetype_send(isg)
                  end if
               end if
            end do
!         end if
      end do
      
      if (neighr.gt.0) then
         call MPI_WAITALL (neighr, reqr, statr, ierr)
      end if

! Now loop though the blocks one final time and if nodetype_recv .ne. nchild
!  and
! derefine = .TRUE. then don't derefine

      do isg = 1,lnblocks
         if (derefine(isg).and.nodetype_recv(isg).ne.nchild) then
            derefine(isg) = .FALSE.
         end if
      end do
      
#ifdef TIMING_MPI 
      timer_amr_check_derefine(4) = timer_amr_check_derefine(4)
     .                + mpi_wtime() - time1
#endif /* TIMING_MPI */


      return
      end subroutine amr_check_derefine



