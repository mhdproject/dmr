!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"


      subroutine amr_refine_derefine





      use paramesh_dimensions
      use physicaldata
      use tree

      use paramesh_interfaces, only : amr_check_derefine,
     &                                amr_refine_blocks,
     &                                amr_derefine_blocks,
     &                                amr_morton_order,
     &                                amr_ser_distribute,
     &                                amr_gsurrounding_blks,
     &                                comm_int_max_to_all,
     &                                comm_int_sum_to_all,
     &                                comm_real_min_to_all,
     &                                comm_real_max_to_all

      implicit none

#include "amr_shmem.fh"
      include 'shmem_reduce.fh'


#ifdef TIMINGS
#include "timer.fh"
#endif


      integer nprocs,mype
      integer shmem_my_pe,shmem_n_pes,i,j,k,l
      integer lnblocks_old

      logical refinet(maxblocks_tr)
      logical lparref

      integer,save :: icontinue,icontinuet
      integer,save :: nrefs,nrefst

      integer lnblocks2,tot_blocks,max_blocks 
      integer tot_ref,tot_deref,lb
      integer parent_blk,parent_pe,cempty(1)

      save    max_blocks
      save    lparref,parent_blk,parent_pe
      save    tot_blocks,lnblocks2,refinet,cempty

      real :: xmin1,xmax1,ymin1,ymax1,zmin1,zmax1
      real :: xmin,xmax,ymin,ymax,zmin,zmax

      logical      :: l_move_solution
      logical,save :: l_arr(1)
      logical,save :: first_call = .true.

#ifdef TIMINGS
      call shmem_barrier_all()
      itimer1 = irtc()
#endif

      nprocs = shmem_n_pes()
      mype = shmem_my_pe()

! reset gsurrblks_set flag
      gsurrblks_set = -1


! error trap for lrefine_max and lrefine_min
      if(lrefine_max.lt.1.or.lrefine_max.gt.100) then
        write(*,*) 'PARAMESH error : lrefine_max has a bad value'
        call amr_abort()
      endif
      if(lrefine_min.lt.1.or.lrefine_min.gt.100.or.
     .   lrefine_min.gt.lrefine_max) then
        write(*,*) 'PARAMESH error : lrefine_min or lrefine_max ',
     .   'has a bad value'
        call amr_abort()
      endif

! enforce refinement level limits
! first upper limit
      where(lrefine==lrefine_max) refine = .false.
! then lower limit
      where(lrefine==lrefine_min) derefine = .false.
! finally force grid to refine toward base level if too coarse
      where( (lrefine<lrefine_min) .and.
     .       (nodetype==1) ) refine = .true.



      if(lnblocks.gt.0) refinet(1:lnblocks) = refine(1:lnblocks)
      newchild(:) = .FALSE.


#ifdef EMPTY_CELLS
! Cancel all empty=2 flags and mark all empty=1 blocks for derefinement
      if(lnblocks.gt.0) then
      do l=1,lnblocks
       if(empty(l).eq.2) empty(l)=0
       if(empty(l).eq.1) derefine(l)=.true.
      enddo
      endif
#endif


! CHECK derefinements
      call amr_check_derefine


! If a leaf block is marked for derefinement but its parent is
! marked for refinement then cancel the leaf block derefinement.
       do l=1,lnblocks
       if(nodetype(l).eq.1.and.derefine(l)) then
         parent_blk = parent(1,l)
         parent_pe = parent(2,l)
         call shmem_logical_get(l_arr(1),refinet(parent_blk),
     .       1,parent_pe)
         lparref = l_arr(1)
         if(lparref) derefine(l)=.false.
       endif
       enddo
 22     call shmem_barrier_all()


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! test to see if any refinement or derefinement has been requested.
        icontinue=0
        if(lnblocks.gt.0) then
! BUG FIXED, added test of nodetype .eq. 1
        do l=1,lnblocks
                if(nodetype(l).eq.1.and.
     &          (refine(l).or.derefine(l))) then
                        icontinue=1
                        goto 10
                endif
        enddo
        endif
10      continue

        do j = 1,shmem_reduce_sync_size
           ipsync(j) = shmem_sync_value
        end do
        icontinuet = icontinue             
        call shmem_barrier_all()
        call comm_int_max_to_all (icontinue,icontinuet)
        call shmem_barrier_all()
!        if(icontinue.eq.0) then
          if(first_call) then

!pmn bug fix to be tested
            lnblocks_old = lnblocks
!pmn end bug fix to be tested
            lnblocks2 = lnblocks 
            call shmem_barrier_all()
            call comm_int_sum_to_all (tot_blocks,lnblocks2)
            call shmem_barrier_all()
            call comm_int_max_to_all (max_blocks,lnblocks2)
            call shmem_barrier_all()


            if (tot_blocks.gt.nprocs) then
! set work values
              work_block(:) = 0.
              if(lnblocks.gt.0) then
              do i = 1,lnblocks
                if (nodetype(i).eq.1) work_block(i) = 2.  !<<< USER EDIT
                if (nodetype(i).ge.2) work_block(i) = 1.  !<<< USER EDIT
              end do
              endif

              l_move_solution = .true.
              call amr_morton_order (nprocs,lnblocks_old,
     .                               l_move_solution)
            else
              call amr_ser_distribute (nprocs,mype,lnblocks_old)
            end if

            first_call = .false.

          else    
            if(icontinue.eq.0) go to 11
          endif                ! end of first_call if test
!        endif                  ! end of icontinue if test

      call shmem_barrier_all()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! compute total no. of requested refinements

      nrefs = 0
      do i = 1,lnblocks
        if (nodetype(i).eq.1) then
          if (refine(i)) nrefs = nrefs + 1
        end if
      end do

      nrefst = nrefs
      call shmem_barrier_all()
      call comm_int_sum_to_all(nrefs,nrefst)
      call shmem_barrier_all()
          
      lnblocks_old = lnblocks


      if (nrefs.gt.0) call amr_refine_blocks
      call amr_derefine_blocks(lnblocks_old)

!      do j = 1,shmem_reduce_sync_size
!         ipsync(j) = shmem_sync_value
!      end do             
! I copy lnblocks to lnblocks2 since lnblocks2 can be put in a save statement.
      lnblocks2 = lnblocks 
      call shmem_barrier_all()
      call comm_int_sum_to_all (tot_blocks,lnblocks2)
      call shmem_barrier_all()
      call comm_int_max_to_all (max_blocks,lnblocks2)
      call shmem_barrier_all()

      if (mype.eq.0) print *,' tot_blocks ',tot_blocks
      if (mype.eq.0) print *,' max_blocks ',max_blocks

      call shmem_barrier_all()


      if (tot_blocks.gt.nprocs) then
! set work values
         work_block(:) = 0.
         if(lnblocks.gt.0) then
         do i = 1,lnblocks
            if (nodetype(i).eq.1) work_block(i) = 2.  !<<< USER EDIT
            if (nodetype(i).ge.2) work_block(i) = 1.  !<<< USER EDIT
         end do
         endif

         l_move_solution = .true.
         call amr_morton_order (nprocs,lnblocks_old,
     .                          l_move_solution)
      else
         call amr_ser_distribute (nprocs,mype,lnblocks_old)
      end if

#ifdef EMPTY_CELLS
        if(lnblocks.gt.0) then
        do l = 1,lnblocks
       cempty(1)=2
       if(empty(l).eq.1) call shmem_integer_put(empty(parent(1,l)),
     .                       cempty(1),1,parent(2,l))
      enddo
      endif
#endif

!
! Fill the surr_blks arrays
 11   call amr_gsurrounding_blks(mype,.true.)

! Find the coordinate ranges
      xmin1 = minval(bnd_box(1,1,1:lnblocks))
      ymin1 = minval(bnd_box(1,2,1:lnblocks))
      zmin1 = minval(bnd_box(1,3,1:lnblocks))
      xmax1 = maxval(bnd_box(2,1,1:lnblocks))
      ymax1 = maxval(bnd_box(2,2,1:lnblocks))
      zmax1 = maxval(bnd_box(2,3,1:lnblocks))
      xmin = min(1.e30, xmin1)
      ymin = min(1.e30, ymin1)
      zmin = min(1.e30, zmin1)
      xmax = max(-1.e30, xmax1)
      ymax = max(-1.e30, ymax1)
      zmax = max(-1.e30, zmax1)
      call comm_real_min_to_all(grid_xmin,xmin)
      call comm_real_min_to_all(grid_ymin,ymin)
      call comm_real_min_to_all(grid_zmin,zmin)
      call comm_real_max_to_all(grid_xmax,xmax)
      call comm_real_max_to_all(grid_ymax,ymax)
      call comm_real_max_to_all(grid_zmax,zmax)
      if(mype.eq.0) then
      write(*,*) ' grid_xmin ',grid_xmin
      write(*,*) ' grid_xmax ',grid_xmax
      write(*,*) ' grid_ymin ',grid_ymin
      write(*,*) ' grid_ymax ',grid_ymax
      write(*,*) ' grid_zmin ',grid_zmin
      write(*,*) ' grid_zmax ',grid_zmax
      endif


! 
! St up an array of cell sizes for each grid refinement level.
! These can be used to minimize variation due to roundoff, but
! should ONLY be used with a uniformly spaced grid.
      level_cell_sizes = 0.
      level_cell_sizes(1,1) = (grid_xmax-grid_xmin)/real(nxb)
      if(ndim.gt.1)
     .  level_cell_sizes(2,1) = (grid_ymax-grid_ymin)/real(nyb)
      if(ndim.eq.3)
     .  level_cell_sizes(3,1) = (grid_zmax-grid_zmin)/real(nzb)
      do i=2,lrefine_max
        level_cell_sizes(1:ndim,i) = .5*level_cell_sizes(1:ndim,i-1)
      enddo

!
! set grid modification flag
      grid_changed = 1



      call shmem_barrier_all()


#ifdef TIMINGS
      itimer2 = irtc()
      call shmem_barrier_all()
      irtc_refderef = itimer2-itimer1+irtc_refderef
#endif

      return
      end subroutine amr_refine_derefine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine amr_ser_distribute (nprocs,mype,lnblocks_old)


      use paramesh_dimensions
      use physicaldata
      use tree

      use paramesh_interfaces, only : amr_redist_blk,
     &                                amr_migrate_tree_data


      implicit none

      integer, intent(in)  ::  nprocs,mype,lnblocks_old
       
      integer new_loc(2,maxblocks_tr)
      integer i,j,k
      logical mark, markt
      save mark,markt
     
      call shmem_barrier_all()

! find new_locs

      mark = .FALSE.
      if (lnblocks.ge.1) mark = .TRUE.

      if (lnblocks.ge.1) then
         new_loc(1,1) = 1
         new_loc(2,1) = mype
      end if

      do i = 0,nprocs-1

         if (i.eq.mype) then

            if (lnblocks.gt.1) then

               do j = 2,lnblocks

                  do k = 0,nprocs-1

                     call shmem_logical_get(markt,mark,1,k)
                     if (.not.markt) then
                        markt = .TRUE.
                        call shmem_logical_put (mark,markt,1,k)
                        new_loc(1,j) = 1
                        new_loc(2,j) = k
                        go to 20
                     end if

                  end do

 20               continue

               end do

            end if

         end if

         call shmem_barrier_all()

      end do

      call shmem_barrier_all()

! distribute tree data
      call amr_migrate_tree_data (new_loc,nprocs)

! distribute block data
      call amr_redist_blk (mype,new_loc,lnblocks_old)

      lnblocks = new_lnblocks

      call shmem_barrier_all()

      return
      end subroutine amr_ser_distribute 
