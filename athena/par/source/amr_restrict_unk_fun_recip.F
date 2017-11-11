!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"

      subroutine amr_restrict_unk_fun_recip(datain,dataout)




!------------------------------------------------------------------------
!
! This routine performs restriction on the array datain and
! returns the result in dataout. Note that this does not update
! guard cell elements of dataout.
!
! This particular version applies the 3D generalization of the
! restriction operator in eqn (19.6.17) of the 2nd edition of
! Numerical recipes.
! The 2D case is
!                  | 1  2  1 |
!                  | 2  4  2 | /16.
!                  | 1  2  1 |
!
!
! Written :     Peter MacNeice          January 1997
!------------------------------------------------------------------------


      use paramesh_dimensions
      use physicaldata

      implicit none

      real, intent(in)    :: datain(:,:,:,:)
      real, intent(inout) :: dataout(:,:,:,:)

      integer :: i,j,k,ivar

!------------------------------------

      if (nvarcorn > 0) then
       do k=1+nguard*k3d,nzb+nguard*k3d+k3d
       do j=1+nguard*k2d,nyb+nguard*k2d+k2d
       do i=1+nguard,nxb+nguard+1
       do ivar=1,nvarcorn

!       dataout(ivar,i,j,k) = (
!     .   ( datain(ivar,i-1,j-k2d,k-k3d) + 2.*datain(ivar,i,j-k2d,k-k3d)
!     .              + datain(ivar,i+1,j-k2d,k-k3d) ) +
!     .  2.*( datain(ivar,i-1,j,k-k3d) + 2.*datain(ivar,i,j,k-k3d) +
!     .              datain(ivar,i+1,j,k-k3d) ) +
!     .   ( datain(ivar,i-1,j+k2d,k-k3d) + 2.*datain(ivar,i,j+k2d,k-k3d)
!     .              + datain(ivar,i+1,j+k2d,k-k3d) ) +
!     .  2.*( datain(ivar,i-1,j-k2d,k) + 2.*datain(ivar,i,j-k2d,k) +
!     .              datain(ivar,i+1,j-k2d,k) ) +
!     .  4.*( datain(ivar,i-1,j,k) + 2.*datain(ivar,i,j,k) +
!     .              datain(ivar,i+1,j,k) ) +
!     .  2.*( datain(ivar,i-1,j+k2d,k) + 2.*datain(ivar,i,j+k2d,k) +
!     .              datain(ivar,i+1,j+k2d,k) ) +
!     .   ( datain(ivar,i-1,j-k2d,k+k3d) + 2.*datain(ivar,i,j-k2d,k+k3d)
!     .              + datain(ivar,i+1,j-k2d,k+k3d) ) +
!     .  2.*( datain(ivar,i-1,j,k+k3d) + 2.*datain(ivar,i,j,k+k3d) +
!     .              datain(ivar,i+1,j,k+k3d) ) +
!     .   ( datain(ivar,i-1,j+k2d,k+k3d) + 2.*datain(ivar,i,j+k2d,k+k3d)
!     .           +   datain(ivar,i+1,j+k2d,k+k3d) )      )/64.

          dataout(ivar,i,j,k) = datain(ivar,i,j,k)

       enddo
       enddo
       enddo
       enddo
       end if

      return
      end subroutine amr_restrict_unk_fun_recip
