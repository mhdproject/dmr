!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"

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
