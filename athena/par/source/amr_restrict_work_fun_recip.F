!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

      subroutine amr_restrict_work_fun_recip(datain,dataout)




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
      use workspace

      implicit none

      real, intent(in)    :: datain(:,:,:)
      real, intent(inout) :: dataout(:,:,:)

!------------------------------------
! local arrays
!      real datain(ilw1:iuw1,jlw1:juw1,klw1:kuw1)
!      real dataout(ilw1:iuw1,jlw1:juw1,klw1:kuw1)

      integer :: i,j,k

!------------------------------------

       do k=1+nguard_work*k3d,nzb+nguard_work*k3d
       do j=1+nguard_work*k2d,nyb+nguard_work*k2d
       do i=1+nguard_work,nxb+nguard_work
       dataout(i,j,k) = (
     .   ( datain(i-1,j-k2d,k-k3d) + 2.*datain(i,j-k2d,k-k3d) +
     .              datain(i+1,j-k2d,k-k3d)   ) +
     .  2.*(   datain(i-1,j,k-k3d) + 2.*datain(i,j,k-k3d) +
     .              datain(i+1,j,k-k3d)   ) +
     .   (   datain(i-1,j+k2d,k-k3d) + 2.*datain(i,j+k2d,k-k3d) +
     .              datain(i+1,j+k2d,k-k3d)   ) +
     .  2.*(   datain(i-1,j-k2d,k) + 2.*datain(i,j-k2d,k) +
     .              datain(i+1,j-k2d,k)   ) +
     .  4.*(   datain(i-1,j,k) +  2.*datain(i,j,k) +
     .              datain(i+1,j,k)   ) +
     .  2.*(   datain(i-1,j+k2d,k) + 2.*datain(i,j+k2d,k) +
     .              datain(i+1,j+k2d,k)   ) +
     .   (   datain(i-1,j-k2d,k+k3d) + 2.*datain(i,j-k2d,k+k3d) +
     .              datain(i+1,j-k2d,k+k3d)   ) +
     .  2.*(   datain(i-1,j,k+k3d) + 2.*datain(i,j,k+k3d) +
     .              datain(i+1,j,k+k3d)   ) +
     .   (   datain(i-1,j+k2d,k+k3d) + 2.*datain(i,j+k2d,k+k3d) +
     .              datain(i+1,j+k2d,k+k3d)   )    )/64.

       enddo
       enddo
       enddo

      return
      end subroutine amr_restrict_work_fun_recip
