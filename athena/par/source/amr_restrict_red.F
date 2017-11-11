!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

      subroutine amr_restrict_red(icoord)




!------------------------------------------------------------------------
!
! This routine performs a user defined reduction operation on the 
! array recvarx(y)(z) and returns the result in bndtempx(y)(z).
! These data arrays are defined on block boundaries only.
!
! Note that this does not update guard cell elements of bndtempx(y)(z).
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
! Written :     Peter MacNeice          February 1997
!------------------------------------------------------------------------


      use paramesh_dimensions
      use physicaldata

      implicit none

      integer, intent(in)    :: icoord


      integer :: i,j,k,ivar

      integer nguard0
      parameter(nguard0 = nguard*npgs)

!------------------------------------

      if(icoord.eq.1) then                         ! x-comp of flux
         do k=1+nguard0*k3d,nzb+nguard0*k3d,2
           do j=1+nguard0*k2d,nyb+nguard0*k2d,2
             do i=1,2
               do ivar=1,nfluxes
                 bndtempx1(ivar,i,j,k) = ( 
     .                        recvarxf(ivar,i,j,k) +
     .                        recvarxf(ivar,i,j+k2d,k) +
     .                        recvarxf(ivar,i,j,k+k3d) +
     .                        recvarxf(ivar,i,j+k2d,k+k3d))
     .                       *red_f
               enddo
             enddo
           enddo
         enddo

      elseif(icoord.eq.2) then                     ! y-comp of flux
        do k=1+nguard0*k3d,nzb+nguard0*k3d,2
          do j=1,2
            do i=1+nguard0,nxb+nguard0,2
              do ivar=1,nfluxes
                bndtempy1(ivar,i,j,k) = ( 
     .                       recvaryf(ivar,i,j,k) +
     .                       recvaryf(ivar,i+1,j,k) +
     .                       recvaryf(ivar,i,j,k+k3d) +
     .                       recvaryf(ivar,i+1,j,k+k3d))
     .                      *red_f
               enddo
             enddo
           enddo
         enddo

      elseif(icoord.eq.3) then                      ! z-comp of flux
        do k=1,2
          do j=1+nguard0,nyb+nguard0,2
            do i=1+nguard0,nxb+nguard0,2
              do ivar=1,nfluxes
                bndtempz1(ivar,i,j,k) = ( 
     .                       recvarzf(ivar,i,j,k) +
     .                       recvarzf(ivar,i+1,j,k) +
     .                       recvarzf(ivar,i,j+1,k) +
     .                       recvarzf(ivar,i+1,j+1,k))
     .                      *red_f
               enddo
             enddo
           enddo
         enddo

      endif

      return
      end subroutine amr_restrict_red
