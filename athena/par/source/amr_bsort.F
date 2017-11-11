!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

      subroutine amr_bi_sort (list,gid,npp)




      implicit none
      
      integer mpp
      parameter (mpp=100000)

      integer list(mpp),dlist(2*mpp)
      integer npp,nppt,nppt2
      integer gid(mpp)
      integer gidt(mpp),gidt2(mpp)
      integer ix(mpp),iproc(mpp)
      integer dix(2*mpp),diproc(2*mpp)
      integer shmem_n_pes,shmem_my_pe
      integer npp_max,near_p_2
      integer i,j,k

      save nppt,nppt2,gidt,gidt2

! find max local list length among processors

!      nppt = npp
!      npp_max = npp
      call shmem_barrier_all()
!      do i = 0,shmem_n_pes()-1

!         call shmem_integer_get(nppt2,nppt,1,i)
!         npp_max = max(nppt2,npp_max)

!      end do
      call comm_int_max_to_all (npp_max,npp)
      call shmem_barrier_all()

! compute nearest power of 2

      near_p_2 = 1
      do while (near_p_2.lt.npp_max)

         near_p_2 = 2*near_p_2

      end do

! load dummy sorting list

      do i = 1,npp
         dlist(i) = list(i)
      end do
      do i = npp+1,mpp
         dlist(i) = 1000000000
      end do

      call shmem_barrier_all()

! sort dummy list using a parallel bitonic sort

      call amr_bsort(dlist,dix,diproc,near_p_2)

! compute global id for each item in the sorted list

      do i = 1,near_p_2

         gidt(i) = i + near_p_2*shmem_my_pe()

      end do

! valid items in sorted list send global ids back to their original
! locations 

      call shmem_barrier_all()

      do i = 1,near_p_2
         if (dlist(i).lt.1000000000) then
            call shmem_integer_put(gidt2(dix(i)),gidt(i),1,diproc(i))
         end if
      end do
      
      call shmem_barrier_all()

      do i = 1,near_p_2
         list(i) = dlist(i)
         gid(i) = gidt2(i)
      end do

      call shmem_barrier_all()

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine amr_bsort (list,ix,iproc,npp)

! Written K. Olson, 6/97

      implicit none

      integer npp
      integer list(npp)
      integer n,nprocs
      integer n_fold,fold_len,next_fold(100000)
      integer ix(npp),iproc(npp)
      integer shmem_n_pes,shmem_my_pe
      integer i,j,k

      interface

      subroutine amr_q_sort (ix,n,ia,ib)

      integer, intent(in) :: n
      integer, dimension(n), intent(inout) :: ix
      integer, optional, dimension(n), intent(inout) :: ia, ib

      end subroutine amr_q_sort

      end interface

      if (npp.gt.100000) print *,' ERROR in bsort: npp > 100000 '

      n_fold = 0
      fold_len = 1
      nprocs = shmem_n_pes()

      do i = 1,npp
         ix(i) = i
         iproc(i) = shmem_my_pe()
      end do

! local sort to produce bitonic sequence on processor
      call amr_q_sort(list,npp,ix,iproc)

      do while (fold_len.le.nprocs/2)

! apply bitonic merge on blocks of procs of length fold_len
         call amr_fold_and_sort (list,ix,iproc,fold_len,npp)
         
         n_fold = n_fold + 1
         next_fold(n_fold) = fold_len

         do j = 1,n_fold-1

            call amr_fold_and_sort (list,ix,iproc,next_fold(j),npp)

            n_fold = n_fold + 1
            next_fold(n_fold) = next_fold(j)

            call shmem_barrier_all()

         end do

         fold_len = 2*fold_len

         call shmem_barrier_all()

      end do

      call shmem_barrier_all()

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine amr_fold_and_sort (list,ix,iproc,fold_len,npp)

      implicit none

      integer nx
      parameter (nx = 4*16384)

      integer npp
      integer list(npp),ix(npp),iproc(npp)
      integer list_temp(npp),ix_temp(npp),iproc_temp(npp)
      integer fold_len,list_len,mype,list_no,left_most,pivot
      integer delta,fetch_from
      integer i,j,k
      integer shmem_my_pe
      integer list_com(nx),ix_com(nx),iproc_com(nx)
      integer list_com2(nx),ix_com2(nx),iproc_com2(nx)

      common /comm_fold/ list_com,ix_com,iproc_com,list_com2,ix_com2,
     &     iproc_com2

      interface

      subroutine amr_q_sort (ix,n,ia,ib)

      integer, intent(in) :: n
      integer, dimension(n), intent(inout) :: ix
      integer, optional, dimension(n), intent(inout) :: ia, ib

      end subroutine amr_q_sort

      end interface

      if (npp.gt.nx) print *,' ERROR IN fold_and_sort: npp > nx '

      mype = shmem_my_pe()

      list_len = fold_len*2
      list_no = int(mype/list_len)
      left_most = list_len*list_no
      pivot = left_most + fold_len - 1
      
      if (mype.le.pivot) then
         delta = pivot - mype + 1
         fetch_from = pivot + delta
      else
         delta = mype - pivot - 1
         fetch_from = pivot - delta
      end if

      do i = 1,npp
         list_com(i) = list(i)
         ix_com(i) = ix(i)
         iproc_com(i) = iproc(i)
      end do

      call shmem_barrier_all()

      call shmem_integer_get (list_com2(1),list_com(1),npp,fetch_from)
      call shmem_integer_get (ix_com2(1),ix_com(1),npp,fetch_from)
      call shmem_integer_get (iproc_com2(1),iproc_com(1),npp,fetch_from)
      
      call shmem_barrier_all()

! change fetched list from monotonically increasing to decreasing
      do i = 1,npp
         k = npp - i + 1
         list_temp(i) = list_com2(k)
         ix_temp(i) = ix_com2(k)
         iproc_temp(i) = iproc_com2(k)
      end do

! compare 2 halves of current list, ipoint by point,
! putting minimum list in first half
! and maximum list in second half.
      if (mype.le.pivot) then

         do i = 1,npp

            if (list_temp(i).lt.list(i)) then
               list(i) = list_temp(i)
               ix(i) = ix_temp(i)
               iproc(i) = iproc_temp(i)
            end if

         end do

      else

         do i = 1,npp

            if (list_temp(i).gt.list(i)) then
               list(i) = list_temp(i)
               ix(i) = ix_temp(i)
               iproc(i) = iproc_temp(i)
            end if

         end do

      end if

      call amr_q_sort(list,npp,ix,iproc)

      call shmem_barrier_all()

      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine amr_q_sort (ix,n,ia,ib)

! quik_sort subroutine
! Written: K. Olson, April, 2003

      implicit none

      integer, intent(in) :: n
      integer, dimension(n),  intent(inout) :: ix
      integer, optional, dimension(n), intent(inout) :: ia, ib
      integer, allocatable, dimension(:) :: ia_temp
      integer, allocatable, dimension(:) :: ib_temp
      integer, allocatable, dimension(:) :: iperm
      integer :: i

      allocate(iperm(size(ix)))
      do i = 1,size(ix)
         iperm(i) = i
      end do

      call q_sort_1(1,size(ix))

      if (present(ia)) then
      allocate(ia_temp(size(ix)))
      ia_temp(:) = ia(:)
      do i = 1,size(ix)
         ia(i) = ia_temp(iperm(i))
      end do
      deallocate(ia_temp)
      end if

      if (present(ib)) then
      allocate(ib_temp(size(ix)))
      ib_temp(:) = ib(:)
      do i = 1,size(ix)
         ib(i) = ib_temp(iperm(i))
      end do
      deallocate(ib_temp)
      end if

      deallocate(iperm)

      contains

      recursive subroutine q_sort_1 (ismall,ibig)

      integer, intent(in) :: ismall, ibig
      integer :: i, j
      integer :: pivot, temp
      integer, parameter :: max_qsort_size = 10

      if (ibig < ismall + max_qsort_size) then

         call simple_sort(ismall, ibig)

      else

         pivot = ix((ismall + ibig)/2)
         i = ismall-1
         j = ibig+1

         do
            do
               i = i + 1
               if (ix(i) >= pivot) exit
            end do

            do
               j = j - 1
               if (ix(j) <= pivot) exit
            end do

            if (i < j) then
               temp = ix(i)
               ix(i) = ix(j)
               ix(j) = temp
               temp = iperm(i)
               iperm(i) = iperm(j)
               iperm(j) = temp
            else if (i == j) then
               i = i + 1
               exit
            else
               exit
            end if
         end do

         if (ismall < j) call q_sort_1(ismall, j)
         if (i < ibig) call q_sort_1(i, ibig)

      end if

      end subroutine q_sort_1

      subroutine simple_sort(ismall, ibig)

      integer, intent(in) :: ismall, ibig
      integer :: i, j
      integer :: temp

      do i = ismall, ibig-1
         do j = i + 1, ibig
            if (ix(i) > ix(j)) then
               temp = ix(i)
               ix(i) = ix(j)
               ix(j) = temp
               temp = iperm(i)
               iperm(i) = iperm(j)
               iperm(j) = temp
            end if
         end do
      end do

      end subroutine simple_sort

      end subroutine amr_q_sort

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
