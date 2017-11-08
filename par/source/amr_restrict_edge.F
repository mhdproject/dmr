!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

      subroutine amr_restrict_edge(icoord)





!------------------------------------------------------------------------
!
! This routine performs a user defined reduction operation on the 
! arrays recvarx(y)(z)[1,2] and returns the result in the same arrays.
! These data arrays are defined on block boundaries only.
!
! Note that this does not update guard cell elements of recvarx(y)(z)[1,2].
!
! Also note that we use stride 2 along each dimension when computing
! reduced data values on block faces, so not all values of dataout
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

      integer, intent(in)    :: icoord

      real :: fact

      integer, parameter ::  nguard0 = nguard*npgs
      integer :: i,j,k

!------------------------------------
      fact = .5
!#ifdef CURVILINEAR
#ifdef EDGE_VALUE_INTEG
      fact = 1.
#endif /* EDGE_VALUE_INTEG */
!#endif /* CURVILINEAR */

      if(icoord.eq.1) then                              ! edges on x-face
        do k=1+nguard0*k3d,nzb+(nguard0+1)*k3d,2
          do j=1+nguard0,nyb+nguard0,2
            do i=1,2
                  ! y pointing edge first
              recvarx1e(:,i,j,k) = (
     .                  recvarx1e(:,i,j,k) +
     .                  recvarx1e(:,i,j+1,k) )*fact
            enddo
          enddo
        enddo

        if(ndim.eq.3) then
          do k=1+nguard0*k3d,nzb+nguard0*k3d,2
            do j=1+nguard0,nyb+nguard0+1,2
              do i=1,2
                    ! z pointing edge 
                recvarx2e(:,i,j,k) = (
     .                    recvarx2e(:,i,j,k) +
     .                    recvarx2e(:,i,j,k+k3d) )*fact
              enddo
            enddo
          enddo
        endif

      elseif(icoord.eq.2) then                          ! edges on y-face
        if(ndim.eq.3) then
          do k=1+nguard0*k3d,nzb+nguard0*k3d,2
            do j=1,2
              do i=1+nguard0,nxb+nguard0+1,2
                  ! z pointing edge first
                recvary1e(:,i,j,k) = (
     .                    recvary1e(:,i,j,k) +
     .                    recvary1e(:,i,j,k+k3d) )*fact
              enddo
            enddo
          enddo
        endif

        do k=1+nguard0*k3d,nzb+(nguard0+1)*k3d,2
          do j=1,2
            do i=1+nguard0,nxb+nguard0,2
                    ! x pointing edge
              recvary2e(:,i,j,k) = (
     .                  recvary2e(:,i,j,k) +
     .                  recvary2e(:,i+1,j,k) )*fact
            enddo
          enddo
        enddo

      elseif(icoord.eq.3) then                          ! edges on z-face
        do k=1,2
          do j=1+nguard0,nyb+nguard0+1,2
            do i=1+nguard0,nxb+nguard0,2
                    ! x pointing edge first
              recvarz1e(:,i,j,k) = (
     .                  recvarz1e(:,i,j,k) +
     .                  recvarz1e(:,i+1,j,k) )*fact
            enddo
          enddo
        enddo

        do k=1,2
          do j=1+nguard0,nyb+nguard0,2
            do i=1+nguard0,nxb+nguard0+1,2
                    ! y pointing edge
              recvarz2e(:,i,j,k) = (
     .                  recvarz2e(:,i,j,k) +
     .                  recvarz2e(:,i,j+1,k) )*fact
            enddo
          enddo
        enddo

      endif

      return
      end subroutine amr_restrict_edge
