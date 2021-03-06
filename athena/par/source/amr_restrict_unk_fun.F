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

      subroutine amr_restrict_unk_fun(datain,dataout,lb)




!------------------------------------------------------------------------
!
! This routine performs restriction on the array datain and
! returns the result in dataout. Note that this does not update
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

      implicit none

      real, intent(in)    :: datain(:,:,:,:)
      real, intent(inout) :: dataout(:,:,:,:)
      integer, intent(in) :: lb

!------------------------------------
! local arrays
!      real datain(nvar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1)
!      real dataout(nvar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1)

      real    :: fact, xi, xj, www
      real,    save :: weight(3,-1:1)

      integer :: i,j,k,ivar
      integer :: i0, j0, k0, is, js, ks, iw, jw, kw
      integer :: iii, jjj, kkk
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


#ifdef FIRST_ORDER_RESTRICTION
       do k=1+nguard*k3d,nzb+nguard*k3d,2
       do j=1+nguard*k2d,nyb+nguard*k2d,2
       do i=1+nguard,nxb+nguard,2
       do ivar=1,nvar
       dataout(ivar,i,j,k) = (
     .       datain(ivar,i,j,k) +
     .       datain(ivar,i+1,j,k) +
     .       datain(ivar,i,j+k2d,k) +
     .       datain(ivar,i+1,j+k2d,k) +
     .       datain(ivar,i,j,k+k3d) +
     .       datain(ivar,i+1,j,k+k3d) +
     .       datain(ivar,i,j+k2d,k+k3d) +
     .       datain(ivar,i+1,j+k2d,k+k3d)) *fact
       enddo
       enddo
       enddo
       enddo

       go to 2
#endif

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

      iparmin = 1+nguard
      iparmax = nxb+nguard
      jparmin = 1+nguard*k2d
      jparmax = nyb+nguard*k2d
      kparmin = 1+nguard*k3d
      kparmax = nzb+nguard*k3d

      end if

      do k0 = kparmin,kparmax,2
      do j0 = jparmin,jparmax,2
      do i0 = iparmin,iparmax,2

        dataout(1:nvar,i0,j0,k0) = 0.
        
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

                 dataout(1:nvar,i0,j0,k0) =
     .                dataout(1:nvar,i0,j0,k0) +
     .                (www*datain(1:nvar,i,j,k))
#ifdef CURVILINEAR
     .                *fact
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
        
 2    return
      end subroutine amr_restrict_unk_fun
