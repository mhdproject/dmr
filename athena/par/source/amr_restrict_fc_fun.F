!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"

#define FIRST_ORDER_RESTRICTION

      subroutine amr_restrict_fc_fun(recv,temp,icoord)




!------------------------------------------------------------------------
!
! This routine performs a user defined reduction operation on the 
! array recv and returns the result in temp.
!
! Note that this does not update guard cell elements of temp.
!
! Also note that we use stride 2 along each dimension when computing
! reduced data values on cell faces, so not all values of temp
! have been updated.
!
!
! This particular version is only appropriate for 2nd order schemes 
! using linear interpolation with even number of mesh points along 
! each block axis.
!
! Written :     Peter MacNeice          July 1997
!------------------------------------------------------------------------


      use paramesh_dimensions
      use physicaldata

      implicit none

      real,    intent(in)    :: recv(:,:,:,:)
      real,    intent(inout) :: temp(:,:,:,:)
      integer, intent(in)    :: icoord

      real    :: fact
      real,save :: weight(2,-1:1)
      real    :: xi, xj, www

      integer :: i,j,k,ivar
      integer :: iparmin,jparmin,kparmin
      integer :: iparmax,jparmax,kparmax
      integer :: iii,jjj,kkk
      integer :: is,js,ks,iw,jw,kw
      integer :: i0, j0, k0

      logical,save :: first = .true.

!------------------------------------
      fact = .25
#ifdef CURVILINEAR
      fact = 1.
      if(ndim.eq.2) fact = .5
#endif /* CURVILINEAR */

      if (first) then

      first = .false.
      xi = -.5
      do i = -1,1
         weight(1,i) = 1.
         xj = -.5
         do j = -1,1
            if (i .ne. j) then
               weight(1,i) =
     &              weight(1,i)*(0.-xj)/(xi-xj)
            end if
            xj = xj + 1.
         end do
         xi = xi + 1.
      end do

      xi = -1.5
      do i = -1,1
         weight(2,i) = 1.
         xj = -1.5
         do j = -1,1
            if (i .ne. j) then
               weight(2,i) =
     &              weight(2,i)*(0.-xj)/(xi-xj)
            end if
            xj = xj + 1.
         end do
         xi = xi + 1.
      end do
      end if

      if(icoord.eq.1) then                         ! x-face variables

#ifdef FIRST_ORDER_RESTRICTION
       do k=1+nguard*k3d,nzb+nguard*k3d,2
       do j=1+nguard*k2d,nyb+nguard*k2d,2
       do i=1+nguard,nxb+nguard+1,2
       do ivar=1,nbndvar
       temp(ivar,i,j,k) = ( 
     .       recv(ivar,i,j,k) +
     .       recv(ivar,i,j+k2d,k) +
     .       recv(ivar,i,j,k+k3d) +
     .       recv(ivar,i,j+k2d,k+k3d)) *fact
       enddo
       enddo
       enddo
       enddo
#else
        iparmin = 1+nguard
        iparmax = nxb+nguard+1
        jparmin = 1+nguard*k2d
        jparmax = nyb+nguard*k2d
        kparmin = 1+nguard*k3d
        kparmax = nzb+nguard*k3d

        do k0 = kparmin,kparmax,2
        do j0 = jparmin,jparmax,2
        do i0 = iparmin,iparmax,2

           if (k0 <= kparmin + (kparmax-kparmin)/2) then
              ks = k0
              kw = 1
           else
              ks = k0-1
              kw = 2
           end if
           if (j0 <= jparmin + (jparmax-jparmin)/2) then
              js = j0
              jw = 1
           else
              js = j0-1
              jw = 2
           end if
           i = i0

!           if (i >= imin .and. i <= imax) then

              temp(1:nbndvar,i0,j0,k0) = 0.

              k = ks
              do kkk = -k3d,k3d
              j = js
              do jjj = -k2d,k2d

                 if (ndim == 1) then
                    www = 1.
                 elseif (ndim == 2) then
                    www = 1.*
     .                   weight(jw,jjj)
                 elseif (ndim == 3) then
                    www = 1.*
     .                   weight(jw,jjj)*
     .                   weight(kw,kkk)
                 end if

                 temp(1:nbndvar,i0,j0,k0) =
     .                temp(1:nbndvar,i0,j0,k0) +
     .                (www*recv(1:nbndvar,i,j,k))
#ifdef CURVILINEAR
     .           *fact
#endif /* CURVILINEAR */


              j = j + 1
              end do
              k = k + 1
              end do

!           end if

        enddo
        enddo
        enddo
#endif   /* not FIRST_ORDER_RESTRICTION*/

      elseif(icoord.eq.2) then                     ! y-face variables
#ifdef FIRST_ORDER_RESTRICTION
       do k=1+nguard*k3d,nzb+nguard*k3d,2
       do j=1+nguard*k2d,nyb+(nguard+1)*k2d,2
       do i=1+nguard,nxb+nguard,2
       do ivar=1,nbndvar
       temp(ivar,i,j,k) = ( 
     .       recv(ivar,i,j,k) +
     .       recv(ivar,i+1,j,k) +
     .       recv(ivar,i,j,k+k3d) +
     .       recv(ivar,i+1,j,k+k3d)) *fact
       enddo
       enddo
       enddo
       enddo
#else
        iparmin = 1+nguard
        iparmax = nxb+nguard
        jparmin = 1+nguard*k2d
        jparmax = nyb+(nguard+1)*k2d
        kparmin = 1+nguard*k3d
        kparmax = nzb+nguard*k3d

        do k0 = kparmin,kparmax,2
        do j0 = jparmin,jparmax,2
        do i0 = iparmin,iparmax,2

           if (k0 <= kparmin + (kparmax-kparmin)/2) then
              ks = k0
              kw = 1
           else
              ks = k0-1
              kw = 2
           end if
           if (i0 <= iparmin + (iparmax-iparmin)/2) then
              is = i0
              iw = 1
           else
              is = i0-1
              iw = 2
           end if
           j = j0

!           if (j >= jmin .and. j <= jmax) then

              temp(1:nbndvar,i0,j0,k0) = 0.

              k = ks
              do kkk = -k3d,k3d
              i = is
              do iii = -1,1

              if (ndim == 1) then
                 www = 1.
              elseif (ndim == 2) then
                 www = 1.*
     .                 weight(iw,iii)
              elseif (ndim == 3) then
                 www = 1.*
     .                 weight(iw,iii)*
     .                 weight(kw,kkk)
              end if

              temp(1:nbndvar,i0,j0,k0) =
     .           temp(1:nbndvar,i0,j0,k0) +
     .           (www*recv(1:nbndvar,i,j,k))
#ifdef CURVILINEAR
     .           *fact
#endif /* CURVILINEAR */

           i = i + 1
           end do
           k = k + 1
           end do

!           end if

        enddo
        enddo
        enddo
#endif /* not FIRST_ORDER_RESTRICTION */
      elseif(icoord.eq.3) then                     ! z-face variables
#ifdef FIRST_ORDER_RESTRICTION
       do k=1+nguard*k3d,nzb+(nguard+1)*k3d,2
       do j=1+nguard*k2d,nyb+nguard*k2d,2
       do i=1+nguard,nxb+nguard,2
       do ivar=1,nbndvar
       temp(ivar,i,j,k) = ( 
     .       recv(ivar,i,j,k) +
     .       recv(ivar,i+1,j,k) +
     .       recv(ivar,i,j+k2d,k) +
     .       recv(ivar,i+1,j+k2d,k)) *fact
       enddo
       enddo
       enddo
       enddo
#else
        iparmin = 1+nguard
        iparmax = nxb+nguard
        jparmin = 1+nguard*k2d
        jparmax = nyb+nguard*k2d
        kparmin = 1+nguard*k3d
        kparmax = nzb+(nguard+1)*k3d

        do k0 = kparmin,kparmax,2
        do j0 = jparmin,jparmax,2
        do i0 = iparmin,iparmax,2

           if (j0 <= jparmin + (jparmax-jparmin)/2) then
              js = j0
              jw = 1
           else
              js = j0-1
              jw = 2
           end if
           if (i0 <= iparmin + (iparmax-iparmin)/2) then
              is = i0
              iw = 1
           else
              is = i0-1
              iw = 2
           end if
           k = k0

!           if (k >= kmin .and. k <= kmax) then

              temp(1:nbndvar,i0,j0,k0) = 0.

              j = js
              do jjj = -k2d,k2d
              i = is
              do iii = -1,1

              www = 1.*
     .              weight(iw,iii)*
     .              weight(jw,jjj)

              temp(1:nbndvar,i0,j0,k0) =
     .           temp(1:nbndvar,i0,j0,k0) +
     .           (www*recv(1:nbndvar,i,j,k))
#ifdef CURVILINEAR
     .           *fact
#endif /* CURVILINEAR */

              i = i + 1
              end do
              j = j + 1
              end do

!           end if

        enddo
        enddo
        enddo

#endif /* FIRST_ORDER_RESTRICTION */
      endif

      return
      end subroutine amr_restrict_fc_fun
