!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"


      subroutine amr_restrict_work_fun(datain,dataout)




!------------------------------------------------------------------------
!
! This routine performs restriction on the array datain and
! returns the result in dataout. Note that this does not updata
! guard cell elements of dataout.
!
! This particular version is only appropriate for 2nd order schemes 
! using linear interpolation with even number of mesh points along 
! each block axis.
!
! Written :     Peter MacNeice          January 1997
!------------------------------------------------------------------------


      use paramesh_dimensions
      use physicaldata
      use workspace

      implicit none

      integer nguard_work0
      parameter(nguard_work0 = nguard_work*npgs)

      real, intent(in)    :: datain(:,:,:)
      real, intent(inout) :: dataout(:,:,:)

!------------------------------------
! local arrays
!      real datain(ilw1:iuw1,jlw1:juw1,klw1:kuw1)
!      real dataout(ilw1:iuw1,jlw1:juw1,klw1:kuw1)

      real    :: fact, xi, xj, www
      real,    save :: weight(3,-1:1)

      integer :: i,j,k, i0, j0, k0, is, js, ks
      integer :: iw, jw, kw, iii, jjj, kkk
      integer, save :: iparmin,iparmax
      integer, save :: jparmin,jparmax
      integer, save :: kparmin,kparmax

      logical, save :: first = .true.

!------------------------------------

       fact = .125
#ifdef CURVILINEAR
       fact = 1.
       if(ndim.eq.2) fact = .5
       if(ndim.eq.1) fact = .25
#endif /* CURVILINEAR */

!       do k=1+nguard_work*k3d,nzb+nguard_work*k3d,2
!       do j=1+nguard_work*k2d,nyb+nguard_work*k2d,2
!       do i=1+nguard_work,nxb+nguard_work,2
!       dataout(i,j,k) = (
!     .       datain(i,j,k) +
!     .       datain(i+1,j,k) +
!     .       datain(i,j+k2d,k) +
!     .       datain(i+1,j+k2d,k) +
!     .       datain(i,j,k+k3d) +
!     .       datain(i+1,j,k+k3d) +
!     .       datain(i,j+k2d,k+k3d) +
!     .       datain(i+1,j+k2d,k+k3d)) *fact
!       enddo
!       enddo
!       enddo

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

      iparmin = 1+nguard_work
      iparmax = nxb+nguard_work
      jparmin = 1+nguard_work*k2d
      jparmax = nyb+nguard_work*k2d
      kparmin = 1+nguard_work*k3d
      kparmax = nzb+nguard_work*k3d

      end if

      do k0 = kparmin,kparmax,2
      do j0 = jparmin,jparmax,2
      do i0 = iparmin,iparmax,2

        dataout(i0,j0,k0) = 0.

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
        if (i0 <= iparmin + (iparmax-iparmin)/2) then
           is = i0
           iw = 1
        else
           is = i0-1
           iw = 2
        end if

       k = ks
        do kkk = -k3d,k3d
           j = js
           do jjj = -k2d,k2d
              i = is
              do iii = -1,1

                 if (ndim == 1) then
                    www = weight(iw,iii)
                 elseif (ndim == 2) then
                    www = weight(iw,iii)*
     .                    weight(jw,jjj)
                 elseif (ndim == 3) then
                    www = weight(iw,iii)*
     .                    weight(jw,jjj)*
     .                    weight(kw,kkk)
                 end if

                 dataout(i0,j0,k0) =
     .                dataout(i0,j0,k0) +
     .                (www*datain(i,j,k))
#ifdef CURVILINEAR
     .                 *fact
#endif /* CURVILINEAR */

                 i = i + 1
              end do
              j = j + 1
           end do
           k = k + 1
        end do

      enddo
      enddo
      enddo

      return
      end subroutine amr_restrict_work_fun
