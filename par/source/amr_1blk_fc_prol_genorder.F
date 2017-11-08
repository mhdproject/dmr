!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"


      subroutine amr_1blk_fc_prol_genorder
     .  (recv,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff,
     .   mype,iface,ivar,order)


!
!------------------------------------------------------------------------
!
! This routine takes data from the array recv, originally extracted 
! from the solution array unk, and performs a prolongation operation 
! on it, between the bounds ranges ia to ib, ja to jb, and ka to kb. 
! The data in recv is from a parent block and the
! result of the prolongation operation is written directly into one
! layer of the working block array unk1(...,idest).
! The position of the child within the parent block is specified by 
! the ioff, joff and koff arguments.
!
! This particular prolongation uses a more general interpolation proceedure then
! some of the other routines provided with PARAMESH.  Any 'order' if interpolation
! can be selected for any variable (as described below).
! It does this by explicitly computing the necessary Taylor expansions out 
! to the specified order.  The interpolations are performed first in the `x' 
! direction.  `Y' interapolations follow, but use the interpolated
! data from the `x' sweep.  The 'Z' sweep is similarly performed.  
!
! To select the `order' (we use the term order here loosely) of interpolation 
! the array interp_mask must have data in it that is >= 0.  
! Since the interpolation scheme is general, one can select
! different orders of interpolation for different variables as.  For instance,
! if,
! interp_mask(1) = 0
! interp_mask(2) = 1
! interp_mask(3) = 2
! then variable 1 will be prolongated used simple direct injection, variable 2
! will be prolongated using linear interpolation and variable 3 will be prolongated
! using quadratic interpolation.
!
! Finally, the `order' of interpolation must be equal or less than nguard.
!
! It is applied to all UNK variables whose corresponding element
! of interp_mask is set to 0.
!
! NOTE: This routine may not be as effcient as some of the other, similar routines
!       provided for prolongation. So, if you don't need the flexibility 
!       of this routine, you might want to consider using another or writing 
!       another yourself.
!
! NOTE2:  This routine does NOT guarantee conservative prologation at refinement
!         jumps.  This is described in the documentation.
!
! Written :     Kevin Olson,  March 2002 and based on similar routines 
!               by Peter MacNeice.
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree

      implicit none

      integer, intent(in) :: ia,ib,ja,jb,ka,kb
      integer, intent(in) :: idest,ioff,joff,koff,mype
      integer, intent(in) :: ivar,iface,order
      real,    intent(inout) :: recv(:,:,:,:)

!------------------------------------
! local arrays

      integer,parameter :: maxorder = 4

      real,save :: weight_right(2,-maxorder-maxorder/2:
     .                             maxorder+maxorder/2,0:maxorder)
      real,save :: weight_left(2,-maxorder-maxorder/2:
     .                            maxorder+maxorder/2,0:maxorder)
      real,save :: weight_half(2,-maxorder-maxorder/2:
     .                            maxorder+maxorder/2,0:maxorder)
      real :: f_intx(iu_bnd1+1,
     .               ju_bnd1+k2d,
     .               ku_bnd1+k3d)
      real :: f_inty(iu_bnd1+1,
     .               ju_bnd1+k2d,
     .               ku_bnd1+k3d)
      real :: f_intz(iu_bnd1+1,
     .               ju_bnd1+k2d,
     .               ku_bnd1+k3d)

      integer :: i,j,k
      integer :: offi,offj,offk
      integer :: ii,jj,kk,iorder
      integer :: icmin,icmax,jcmin,jcmax,kcmin,kcmax
      integer :: ifmin,ifmax,jfmin,jfmax,kfmin,kfmax
      integer :: imin, imax, jmin, jmax, kmin, kmax
      integer,save :: iminl(2,0:maxorder)
      integer,save :: imaxl(2,0:maxorder)
      integer,save :: iminr(2,0:maxorder)
      integer,save :: imaxr(2,0:maxorder)
      integer,save :: iminh(2,0:maxorder)
      integer,save :: imaxh(2,0:maxorder)
      integer :: iw, jw, kw, ipar, jpar, kpar
      integer,parameter :: largei = 100

      logical,save :: first_call = .true.
      
!------------------------------------


      if (first_call) then
         first_call = .false.

         do iorder = 0,maxorder
! LEFT

            iminl(1,iorder) = -iorder/2
            imaxl(1,iorder) = -iorder/2 + iorder

            do ipar = iminl(1,iorder),imaxl(1,iorder)
               weight_left(1,ipar,iorder) = 1.
               do jpar = iminl(1,iorder),imaxl(1,iorder)
                  if (jpar.ne.ipar) then
                     weight_left(1,ipar,iorder) =
     .                    weight_left(1,ipar,iorder)*
     .                    (-.25-jpar)/(ipar-jpar)
                  end if
               end do
            end do

            iminl(2,iorder) = -iorder
            imaxl(2,iorder) = 0

            do ipar = iminl(2,iorder),imaxl(2,iorder)
               weight_left(2,ipar,iorder) = 1.
               do jpar = iminl(2,iorder),imaxl(2,iorder)
                  if (jpar.ne.ipar) then
                     weight_left(2,ipar,iorder) =
     .                    weight_left(2,ipar,iorder)*
     .                    (-.25-jpar)/(ipar-jpar)
                  end if
               end do
            end do
! RIGHT

            iminr(1,iorder) = 0
            imaxr(1,iorder) = iorder

            do ipar = iminr(1,iorder),imaxr(1,iorder)
               weight_right(1,ipar,iorder) = 1.
               do jpar = iminr(1,iorder),imaxr(1,iorder)
                  if (jpar.ne.ipar) then
                     weight_right(1,ipar,iorder) =
     .                    weight_right(1,ipar,iorder)*
     .                    (.25-jpar)/(ipar-jpar)
                  end if
               end do
            end do

            iminr(2,iorder) = iorder/2 - iorder
            imaxr(2,iorder) = iorder/2

            do ipar = iminr(2,iorder),imaxr(2,iorder)
               weight_right(2,ipar,iorder) = 1.
               do jpar = iminr(2,iorder),imaxr(2,iorder)
                  if (jpar.ne.ipar) then
                     weight_right(2,ipar,iorder) =
     .                    weight_right(2,ipar,iorder)*
     .                    (.25-jpar)/(ipar-jpar)
                  end if
               end do
            end do
! HALF CELL TO THE RIGHT of CELL FACE

            iminh(1,iorder) = 0
            imaxh(1,iorder) = iorder

            do ipar = iminh(1,iorder),imaxh(1,iorder)
               weight_half(1,ipar,iorder) = 1.
               do jpar = iminh(1,iorder),imaxh(1,iorder)
                  if (jpar.ne.ipar) then
                     weight_half(1,ipar,iorder) =
     .                    weight_half(1,ipar,iorder)*
     .                    (.5-jpar)/(ipar-jpar)
                  end if
               end do
            end do

            iminh(2,iorder) = iorder/2 - iorder
            imaxh(2,iorder) = iorder/2

            do ipar = iminh(2,iorder),imaxh(2,iorder)
               weight_half(2,ipar,iorder) = 1.
               do jpar = iminh(2,iorder),imaxh(2,iorder)
                  if (jpar.ne.ipar) then
                     weight_half(2,ipar,iorder) =
     .                    weight_half(2,ipar,iorder)*
     .                    (.5-jpar)/(ipar-jpar)
                  end if
               end do
            end do

         end do

      end if                    ! end if (first_call


! Set the bounds on the loop controlling the interpolation.
      ifmin=ia
      ifmax=ib
      jfmin=ja
      jfmax=jb
      kfmin=ka
      kfmax=kb


      offi = 0
      offj = 0
      offk = 0
      if(ioff.gt.0) offi = nxb/2
      if(joff.gt.0) offj = nyb*k2d/2
      if(koff.gt.0) offk = nzb*k3d/2

      kcmin = ((kfmin-nguard-1+largei)/2 +
     .                nguard - largei/2 )*k3d + 
     .                1 + offk
      kcmax = ((kfmax-nguard-1+largei)/2 +
     .                nguard - largei/2 )*k3d + 
     .                1 + offk
      jcmin = ((jfmin-nguard-1+largei)/2 +
     .                nguard - largei/2 )*k2d + 
     .                1 + offj
      jcmax = ((jfmax-nguard-1+largei)/2 +
     .                nguard - largei/2 )*k2d + 
     .                1 + offj
      icmin = ((ifmin-nguard-1+largei)/2 +
     .                nguard - largei/2 ) + 
     .                1 + offi
      icmax = ((ifmax-nguard-1+largei)/2 +
     .                nguard - largei/2 ) + 
     .                1 + offi





! Main Interpolation loop.





      if (iface == 1) then ! interpolate facevarx


! Interpolate in x direction





      if (ndim >= 1) then


      f_intx(:,:,:) = 0.

      kmin = kcmin-nguard*k3d
      if (kmin < 1) kmin = 1
      kmax = kcmax+nguard*k3d 
      if (kmax > nguard*2+nzb) kmax = nguard*2+nzb
      jmin = jcmin-nguard*k2d
      if (jmin < 1) jmin = 1
      jmax = jcmax+nguard*k2d
      if (jmax > nguard*2+nyb) jmax = nguard*2+nyb 

      do k = kmin,kmax
      do j = jmin,jmax

         ! 1) now interpolate to half points

         ! starting parent index
         i = icmin
         
         do ii = ifmin,ifmax

            if ((mod(ii,2) .ne. 0 .and. mod(nguard,2)  ==  0) .or.
     .          (mod(ii,2)  ==  0 .and. mod(nguard,2) .ne. 0)) then
                                       ! this point is on parent's face
                                       ! and does not need to be interpolated
               f_intx(ii,j,k) = recv(ivar,i,j,k)
            else
               if (ii < nguard + nxb/2) then
                  iw = 1
               else
                  iw = 2
               end if

               imin = i + iminh(iw,order)
               imax = i + imaxh(iw,order)

               do ipar = imin,imax
                  f_intx(ii,j,k) = f_intx(ii,j,k) +
     .                 weight_half(iw,ipar-i,order)*
     .                 recv(ivar,ipar,j,k)
               end do
               ! update parent index
               i = i + 1
            end if

         end do                 ! end loop over ii
      end do                    ! end loop over j
      end do                    ! end loop over k

      if (ndim == 1) then
         do k = kfmin,kfmax
            do j = jfmin,jfmax
               do i = ifmin,ifmax
                  facevarx1(ivar,i,j,k,idest) = f_intx(i,j,k)
               end do
            end do
         end do
      end if

      end if                    ! end if (ndim




! Interpolate in y direction




      if (ndim >= 2) then
      

      f_inty(:,:,:) = 0.

      kmin = kcmin-nguard*k3d
      if (kmin < 1) kmin = 1
      kmax = kcmax+nguard*k3d
      if (kmax > nguard*2+nzb) kmax = nguard*2+nzb

      do k = kmin,kmax
      do i = ifmin,ifmax

         ! 1) interpolate to half points

         ! starting parent index
         j = jcmin
            
         do jj = jfmin,jfmax
            
            if ((mod(jj,2) .ne. 0 .and. mod(nguard,2) .ne. 0) .or.
     .          (mod(jj,2)  ==  0 .and. mod(nguard,2)  ==  0)) then

                                ! right point

               if (jj < nguard + nyb/2) then
                  jw = 1
               else
                  jw = 2
               end if

               jmin = iminr(jw,order) + j
               jmax = imaxr(jw,order) + j

               do jpar = jmin,jmax
                  f_inty(i,jj,k) = f_inty(i,jj,k) +
     .                 weight_right(jw,jpar-j,order)*f_intx(i,jpar,k)
               end do
               ! update parent index
               j = j + 1

            else

                                ! left point

               if (jj < nguard + nyb/2) then
                  jw = 1
               else
                  jw = 2
               end if

               jmin = iminl(jw,order) + j
               jmax = imaxl(jw,order) + j

               do jpar = jmin,jmax
                  f_inty(i,jj,k) = f_inty(i,jj,k) +
     .                 weight_left(jw,jpar-j,order)*f_intx(i,jpar,k)
               end do

            end if
            
         end do                 ! end loop over jj
      end do                    ! end loop over i
      end do                    ! end loop over k

      if (ndim == 2) then
         do k = kfmin,kfmax
            do j = jfmin,jfmax
               do i = ifmin,ifmax
                  facevarx1(ivar,i,j,k,idest) = f_inty(i,j,k)
               end do
            end do
         end do
      end if

      end if                    ! end if (ndim




      
! Interpolate in z direction





      if (ndim == 3) then


      do j = jfmin,jfmax
      do i = ifmin,ifmax

         ! 1) interpolate to half points

         ! starting parent index
         k = kcmin

         do kk = kfmin,kfmax
            
            facevarx1(ivar,i,j,kk,idest) = 0.
            if ((mod(kk,2) .ne. 0 .and. mod(nguard,2) .ne. 0) .or.
     .          (mod(kk,2)  ==  0 .and. mod(nguard,2)  ==  0)) then

                                ! right point

               if (kk < nguard + nzb/2) then
                  kw = 1
               else
                  kw = 2
               end if

               kmin = iminr(kw,order) + k
               kmax = imaxr(kw,order) + k

               do kpar = kmin,kmax
                  facevarx1(ivar,i,j,kk,idest) = 
     .                 facevarx1(ivar,i,j,kk,idest) +
     .                 weight_right(kw,kpar-k,order)*f_inty(i,j,kpar)
               end do
               ! update parent index
               k = k + 1

            else

                                ! left point

               if (kk < nguard + nzb/2) then
                  kw = 1
               else
                  kw = 2
               end if

               kmin = iminl(kw,order) + k
               kmax = imaxl(kw,order) + k

               do kpar = kmin,kmax
                  facevarx1(ivar,i,j,kk,idest) = 
     .                 facevarx1(ivar,i,j,kk,idest) +
     .                 weight_left(kw,kpar-k,order)*f_inty(i,j,kpar)
               end do

            end if
            
         end do                 ! end loop over kk
      end do                    ! end loop over j
      end do                    ! end loop over i

      end if                    ! end if (ndim





      
      elseif (iface == 2) then ! interpolate facevary


! Interpolate in y direction





      if (ndim >= 1) then


      f_inty(:,:,:) = 0.

      kmin = kcmin-nguard*k3d
      if (kmin < 1) kmin = 1
      kmax = kcmax+nguard*k3d
      if (kmax > nguard*2+nzb) kmax = nguard*2+nzb
      imin = icmin-nguard
      if (imin < 1) imin = 1
      imax = icmax+nguard
      if (imax > nguard*2+nxb) imax = nguard*2+nxb

      do k = kmin,kmax
      do i = imin,imax

         ! 1) interpolate to half points

         ! starting parent index
         j = jcmin
         
         do jj = jfmin,jfmax

            if ((mod(jj,2) .ne. 0 .and. mod(nguard,2)  ==  0) .or.
     .          (mod(jj,2)  ==  0 .and. mod(nguard,2) .ne. 0)) then
                                       ! this point is on the parent's face
                                       ! and does not need to be interpolated
               f_inty(i,jj,k) = recv(ivar,i,j,k)
            else

               if (jj < nguard + nyb/2) then
                  jw = 1
               else
                  jw = 2
               end if

               jmin = iminh(jw,order) + j
               jmax = imaxh(jw,order) + j

               do jpar = jmin,jmax
                  f_inty(i,jj,k) = f_inty(i,jj,k) +
     .                 weight_half(jw,jpar-j,order)*
     .                 recv(ivar,i,jpar,k)
               end do
               ! update parent index
               j = j + 1
            end if

         end do                 ! end loop over jj
      end do                    ! end loop over i
      end do                    ! end loop over k

      if (ndim == 1) then
         do k = kfmin,kfmax
            do j = jfmin,jfmax
               do i = ifmin,ifmax
                  facevary1(ivar,i,j,k,idest) = f_inty(i,j,k)
               end do
            end do
         end do
      end if

      end if                    ! end if (ndim




! Interpolate in x direction




      if (ndim >= 2) then
      

      f_intx(:,:,:) = 0.

      kmin = kcmin-nguard*k3d
      if (kmin < 1) kmin = 1
      kmax = kcmax+nguard*k3d
      if (kmax > nguard*2+nzb) kmax = nguard*2+nzb

      do k = kmin,kmax
      do j = jfmin,jfmax

         ! 1) interpolate to half points

         ! starting parent index
         i = icmin
            
         do ii = ifmin,ifmax
            
            if ((mod(ii,2) .ne. 0 .and. mod(nguard,2) .ne. 0) .or.
     .          (mod(ii,2)  ==  0 .and. mod(nguard,2)  ==  0)) then

                                ! right point

               if (ii < nguard + nxb/2) then
                  iw = 1
               else
                  iw = 2
               end if

               imin = iminr(iw,order) + i
               imax = imaxr(iw,order) + i

               do ipar = imin,imax
                  f_intx(ii,j,k) = f_intx(ii,j,k) +
     .                 weight_right(iw,ipar-i,order)*f_inty(ipar,j,k)
               end do
               ! update parent index
               i = i + 1

            else

                                ! left point

               if (ii < nguard + nxb/2) then
                  iw = 1
               else
                  iw = 2
               end if

               imin = iminl(iw,order) + i
               imax = imaxl(iw,order) + i

               do ipar = imin,imax
                  f_intx(ii,j,k) = f_intx(ii,j,k) +
     .                 weight_left(iw,ipar-i,order)*f_inty(ipar,j,k)
               end do

            end if
            
         end do                 ! end loop over ii
      end do                    ! end loop over j
      end do                    ! end loop over k

      if (ndim == 2) then
         do k = kfmin,kfmax
            do j = jfmin,jfmax
               do i = ifmin,ifmax
                  facevary1(ivar,i,j,k,idest) = f_intx(i,j,k)
               end do
            end do
         end do
      end if

      end if                    ! end if (ndim




      
! Interpolate in z direction





      if (ndim == 3) then


      do j = jfmin,jfmax
      do i = ifmin,ifmax

         ! 2) now interpolate to half points

         ! starting parent index
         k = kcmin

         do kk = kfmin,kfmax
            
            facevary1(ivar,i,j,kk,idest) = 0.
            if ((mod(kk,2) .ne. 0 .and. mod(nguard,2) .ne. 0) .or.
     .          (mod(kk,2)  ==  0 .and. mod(nguard,2)  ==  0)) then

                                ! right point

               if (kk < nguard + nzb/2) then
                  kw = 1
               else
                  kw = 2
               end if

               kmin = iminr(kw,order) + k
               kmax = imaxr(kw,order) + k

               do kpar = kmin,kmax
                  facevary1(ivar,i,j,kk,idest) =
     .                 facevary1(ivar,i,j,kk,idest) +
     .                 weight_right(kw,kpar-k,order)*f_intx(i,j,kpar)
               end do
               ! update parent index
               k = k + 1

            else

                                ! left point

               if (kk < nguard + nzb/2) then
                  kw = 1
               else
                  kw = 2
               end if
               
               kmin = iminl(kw,order) + k
               kmax = imaxl(kw,order) + k

               do kpar = kmin,kmax
                  facevary1(ivar,i,j,kk,idest) =
     .                 facevary1(ivar,i,j,kk,idest) +
     .                 weight_left(kw,kpar-k,order)*f_intx(i,j,kpar)
               end do

            end if
            
         end do                 ! end loop over kk
      end do                    ! end loop over j
      end do                    ! end loop over i

      end if                    ! end if (ndim






      
      elseif (iface == 3) then ! interpolate facevarz




! Interpolate in z direction





      if (ndim >= 1) then


      f_intz(:,:,:) = 0.

      imin = icmin-nguard
      if (imin < 1) imin = 1
      imax = icmax+nguard
      if (imax > nguard*2+nxb) imax = nguard*2+nxb
      jmin = jcmin-nguard*k2d
      if (jmin < 1) jmin = 1
      jmax = jcmax+nguard*k2d
      if (jmax > nguard*2+nyb) jmax = nguard*2+nyb

      do j = jmin,jmax
      do i = imin,imax

         ! 1) interpolate to half points

         ! starting parent index
         k = kcmin
         
         do kk = kfmin,kfmax

            if ((mod(kk,2) .ne. 0 .and. mod(nguard,2)  ==  0) .or.
     .          (mod(kk,2)  ==  0 .and. mod(nguard,2) .ne. 0)) then
                                       ! and does not need to be interpolated
               f_intz(i,j,kk) = recv(ivar,i,j,k) 
            else

               if (kk < nguard + nzb/2) then
                  kw = 1
               else
                  kw = 2
               end if

               kmin = iminh(kw,order) + k
               kmax = imaxh(kw,order) + k

               do kpar = kmin,kmax
                  f_intz(i,j,kk) = f_intz(i,j,kk) +
     .                 weight_half(kw,kpar-k,order)*
     .                 recv(ivar,i,j,kpar)
               end do
               ! update parent index
               k = k + 1
            end if

         end do                 ! end loop over ii
      end do                    ! end loop over j
      end do                    ! end loop over k

      if (ndim == 1) then
         do k = kfmin,kfmax
            do j = jfmin,jfmax
               do i = ifmin,ifmax
                  facevarz1(ivar,i,j,k,idest) = f_intz(i,j,k)
               end do
            end do
         end do
      end if

      end if                    ! end if (ndim




! Interpolate in y direction




      if (ndim >= 2) then
      

      f_inty(:,:,:) = 0.

      imin = icmin-nguard
      if (imin < 1) imin = 1
      imax = icmax+nguard
      if (imax > nguard*2+nxb) imax = nguard*2+nxb

      do k = kfmin,kfmax
      do i = imin,imax

         ! 1) interpolate to half points

         ! starting parent index
         j = jcmin
            
         do jj = jfmin,jfmax
            
            if ((mod(jj,2) .ne. 0 .and. mod(nguard,2) .ne. 0) .or.
     .          (mod(jj,2)  ==  0 .and. mod(nguard,2)  ==  0)) then

                                ! right point

               if (jj < nguard + nyb/2) then
                  jw = 1
               else
                  jw = 2
               end if

               jmin = iminr(jw,order) + j
               jmax = imaxr(jw,order) + j

               do jpar = jmin,jmax
                  f_inty(i,jj,k) = f_inty(i,jj,k) +
     .                 weight_right(jw,jpar-j,order)*f_intz(i,jpar,k)
               end do
               ! update parent index
               j = j + 1

            else
                                ! left point

               if (jj < nguard + nyb/2) then
                  jw = 1
               else
                  jw = 2
               end if

               jmin = iminl(jw,order) + j
               jmax = imaxl(jw,order) + j

               do jpar = jmin,jmax
                  f_inty(i,jj,k) = f_inty(i,jj,k) +
     .                 weight_left(jw,jpar-j,order)*f_intz(i,jpar,k)
               end do

            end if
            
         end do                 ! end loop over jj
      end do                    ! end loop over i
      end do                    ! end loop over k

      if (ndim == 2) then
         do k = kfmin,kfmax
            do j = jfmin,jfmax
               do i = ifmin,ifmax
                  facevarz1(ivar,i,j,k,idest) = f_inty(i,j,k)
               end do
            end do
         end do
      end if

      end if                    ! end if (ndim




      
! Interpolate in x direction





      if (ndim == 3) then


      do k = kfmin,kfmax
      do j = jfmin,jfmax

         ! 1) interpolate to half points

         ! starting parent index
         i = icmin

         do ii = ifmin,ifmax
            
            facevarz1(ivar,ii,j,k,idest) = 0.
            if ((mod(ii,2) .ne. 0 .and. mod(nguard,2) .ne. 0) .or.
     .          (mod(ii,2)  ==  0 .and. mod(nguard,2)  ==  0)) then

                                ! right point

               if (ii < nguard + nxb/2) then
                  iw = 1
               else
                  iw = 2
               end if

               imin = iminr(iw,order) + i
               imax = imaxr(iw,order) + i

               do ipar = imin,imax
                  facevarz1(ivar,ii,j,k,idest)= 
     .                 facevarz1(ivar,ii,j,k,idest) +
     .                 weight_right(iw,ipar-i,order)*f_inty(ipar,j,k)
               end do
               ! update parent index
               i = i + 1

            else

                                ! left point

               if (ii < nguard + nxb/2) then
                  iw = 1
               else
                  iw = 2
               end if

               imin = iminl(iw,order) + i
               imax = imaxl(iw,order) + i

               do ipar = imin,imax
                  facevarz1(ivar,ii,j,k,idest) = 
     .                 facevarz1(ivar,ii,j,k,idest) +
     .                 weight_left(iw,ipar-i,order)*f_inty(ipar,j,k)
               end do

            end if
            
         end do                 ! end loop over ii
      end do                    ! end loop over j
      end do                    ! end loop over k

      end if                    ! end if (ndim

      end if                    ! if (iface
      
      end subroutine amr_1blk_fc_prol_genorder






