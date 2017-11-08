!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"


      subroutine amr_1blk_cc_prol_genorder
     .  (recv,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff,
     .   mype,ivar,order)


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
! This particular prolongation uses a more general interpolation proceedure
! than some of the other routines provided with PARAMESH.  Any 'order' of
! interpolation can be selected for any variable (as described below).
! This routine works by using a general polynomial fitting algorithm.
! The interpolations are performed first in the `x'
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
! will be prolongated using linear interpolation and variable 3 will be
! prolongated using quadratic interpolation.
!
! Finally, the `order' of interpolation must be equal or less than nguard.
! This ensures that enough guardcells space is available to compute
! the interpolation weights for the polynomial fits.
!
! It is applied to all UNK variables whose corresponding element
! of interp_mask is set to 0.
!
! Finally, the `order' of interpolation must be equal or less than nguard.
! This ensures that enough guardcells space is available to compute
! the interpolation weights for the polynomial fits.
!
! It is applied to all UNK variables whose corresponding element
! of interp_mask is set to 0.
!
! NOTE: This routine may not be as effcient as some of the other, similar
!       routine provided for prolongation. So, if you don't need the
!       flexibility of this routine, you might want to consider using another
!       or writing another yourself.
!
! NOTE2:  This routine does NOT guarantee conservative prologation at
!         refinement jumps.  This is described in the documentation.
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
      integer, intent(in) :: ivar,order
      real,    intent(inout) :: recv(:,:,:,:)


!------------------------------------
! local arrays

      real :: weight_right
      real :: weight_left
      real :: tempy, tempx

      integer :: i,j,k,ii
      integer :: offi,offj,offk
      integer,parameter :: largei = 100
      integer,parameter :: maxorder = 4
      integer :: iorder
      integer :: icmin,jcmin,kcmin
      integer :: ifmin,ifmax,jfmin,jfmax,kfmin,kfmax
      integer :: ipar, jpar, kpar

#ifndef LIBRARY
      integer,save :: imina(iu_bnd1,0:maxorder),
     .                imaxa(iu_bnd1,0:maxorder)
      real,save    :: weightx(0:iu_bnd1,iu_bnd1,0:maxorder)

      integer,save :: jmina(ju_bnd1,0:maxorder),
     .                jmaxa(ju_bnd1,0:maxorder)
      real,save    :: weighty(0:ju_bnd1,ju_bnd1,0:maxorder)

      integer,save :: kmina(ku_bnd1,0:maxorder),
     .                kmaxa(ku_bnd1,0:maxorder)
      real,save    :: weightz(0:ku_bnd1,ku_bnd1,0:maxorder)
#else
      integer,save,allocatable :: imina(:,:),
     .                            imaxa(:,:)
      real,save,allocatable    :: weightx(:,:,:)

      integer,save,allocatable :: jmina(:,:),
     .                            jmaxa(:,:)
      real,save,allocatable    :: weighty(:,:,:)

      integer,save,allocatable :: kmina(:,:),
     .                            kmaxa(:,:)
      real,save,allocatable    :: weightz(:,:,:)
#endif
      
      logical,save :: first_call = .true.

!------------------------------------

      if (first_call) then
#ifdef LIBRARY
         allocate (imina(iu_bnd1,0:maxorder))
         allocate (imaxa(iu_bnd1,0:maxorder))
         allocate (weightx(iu_bnd1,iu_bnd1,0:maxorder))

         allocate (jmina(ju_bnd1,0:maxorder))
         allocate (jmaxa(ju_bnd1,0:maxorder))
         allocate (weighty(ju_bnd1,ju_bnd1,0:maxorder))

         allocate (kmina(ku_bnd1,0:maxorder))
         allocate (kmaxa(ku_bnd1,0:maxorder))
         allocate (weightz(ku_bnd1,ku_bnd1,0:maxorder))
#endif
         first_call = .false.

         do iorder = 0,maxorder
         
!!! XXXXX !!!

         i = ((1-nguard-1+largei)/2 +
     .         nguard - largei/2 ) + 1 
         do ii = 1,iu_bnd1
            
            if ((mod(ii,2) .ne. 0 .and. mod(nguard,2) .ne. 0) .or.
     .          (mod(ii,2)  ==  0 .and. mod(nguard,2)  ==  0)) then
               
                                ! right point
               
               if (ii < nguard + nxb/2) then
                  imina(ii,iorder) = i
                  imaxa(ii,iorder) = i + iorder
               else
                  imina(ii,iorder) = i + iorder/2 - iorder
                  imaxa(ii,iorder) = i + iorder/2
               end if

               do ipar = imina(ii,iorder),imaxa(ii,iorder)
                  weight_right = 1.
                  do jpar = imina(ii,iorder),imaxa(ii,iorder)
                     if (jpar.ne.ipar) then
                        weight_right = 
     .                   weight_right*(.25-(jpar-i))/(ipar-jpar)
                     end if
                  end do
                  weightx(ipar,ii,iorder) = weight_right
               end do
                                ! update parent index
               i = i + 1
               
            else
                                ! left point
               
               if (ii < nguard + nxb/2) then
                  imina(ii,iorder) = i - iorder/2
                  imaxa(ii,iorder) = i - iorder/2 + iorder
               else
                  imina(ii,iorder) = i - iorder
                  imaxa(ii,iorder) = i
               end if

               do ipar = imina(ii,iorder),imaxa(ii,iorder)
                  weight_left = 1.
                  do jpar = imina(ii,iorder),imaxa(ii,iorder)
                     if (jpar.ne.ipar) then
                        weight_left = 
     .                   weight_left*(-.25-(jpar-i))/(ipar-jpar)
                     end if
                  end do
                  weightx(ipar,ii,iorder) = weight_left
               end do
            
            end if

         end do                 ! end loop over ii

!!! YYYYY !!!

         if (ndim >= 2) then

         i = ((1-nguard-1+largei)/2 +
     .         nguard - largei/2 ) + 1 
         do ii = 1,ju_bnd1
            
            if ((mod(ii,2) .ne. 0 .and. mod(nguard,2) .ne. 0) .or.
     .          (mod(ii,2)  ==  0 .and. mod(nguard,2)  ==  0)) then
               
                                ! right point
               
               if (ii < nguard + nyb/2) then
                  jmina(ii,iorder) = i
                  jmaxa(ii,iorder) = i + iorder
               else
                  jmina(ii,iorder) = i + iorder/2 - iorder
                  jmaxa(ii,iorder) = i + iorder/2
               end if

               do ipar = jmina(ii,iorder),jmaxa(ii,iorder)
                  weight_right = 1.
                  do jpar = jmina(ii,iorder),jmaxa(ii,iorder)
                     if (jpar.ne.ipar) then
                        weight_right = 
     .                   weight_right*(.25-(jpar-i))/(ipar-jpar)
                     end if
                  end do
                  weighty(ipar,ii,iorder) = weight_right
               end do
                                ! update parent index
               i = i + 1
               
            else
                                ! left point
               
               if (ii < nguard + nyb/2) then
                  jmina(ii,iorder) = i - iorder/2
                  jmaxa(ii,iorder) = i - iorder/2 + iorder
               else
                  jmina(ii,iorder) = i - iorder
                  jmaxa(ii,iorder) = i
               end if

               do ipar = jmina(ii,iorder),jmaxa(ii,iorder)
                  weight_left = 1.
                  do jpar = jmina(ii,iorder),jmaxa(ii,iorder)
                     if (jpar.ne.ipar) then
                        weight_left = 
     .                   weight_left*(-.25-(jpar-i))/(ipar-jpar)
                     end if
                  end do
                  weighty(ipar,ii,iorder) = weight_left
               end do
            
            end if

         end do                 ! end loop over ii

         end if                 ! if (ndim

!!! ZZZZZ !!!

         if (ndim == 3) then

         i = ((1-nguard-1+largei)/2 +
     .         nguard - largei/2 ) + 1 
         do ii = 1,ku_bnd1
            
            if ((mod(ii,2) .ne. 0 .and. mod(nguard,2) .ne. 0) .or.
     .          (mod(ii,2)  ==  0 .and. mod(nguard,2)  ==  0)) then
               
                                ! right point
               
               if (ii < nguard + nzb/2) then
                  kmina(ii,iorder) = i
                  kmaxa(ii,iorder) = i + iorder
               else
                  kmina(ii,iorder) = i + iorder/2 - iorder
                  kmaxa(ii,iorder) = i + iorder/2
               end if

               do ipar = kmina(ii,iorder),kmaxa(ii,iorder)
                  weight_right = 1.
                  do jpar = kmina(ii,iorder),kmaxa(ii,iorder)
                     if (jpar.ne.ipar) then
                        weight_right = 
     .                   weight_right*(.25-(jpar-i))/(ipar-jpar)
                     end if
                  end do
                  weightz(ipar,ii,iorder) = weight_right
               end do
                                ! update parent index
               i = i + 1
               
            else
                                ! left point
               
               if (ii < nguard + nzb/2) then
                  kmina(ii,iorder) = i - iorder/2
                  kmaxa(ii,iorder) = i - iorder/2 + iorder
               else
                  kmina(ii,iorder) = i - iorder
                  kmaxa(ii,iorder) = i
               end if

               do ipar = kmina(ii,iorder),kmaxa(ii,iorder)
                  weight_left = 1.
                  do jpar = kmina(ii,iorder),kmaxa(ii,iorder)
                     if (jpar.ne.ipar) then
                        weight_left = 
     .                   weight_left*(-.25-(jpar-i))/(ipar-jpar)
                     end if
                  end do
                  weightz(ipar,ii,iorder) = weight_left
               end do
            
            end if

         end do                 ! end loop over ii

         end if                 ! end if (ndim

         end do                 ! end loop over iorder

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
      jcmin = ((jfmin-nguard-1+largei)/2 +
     .                nguard - largei/2 )*k2d + 
     .                1 + offj
      icmin = ((ifmin-nguard-1+largei)/2 +
     .                nguard - largei/2 ) + 
     .                1 + offi





! Main Interpolation loop.






! Interpolate !





      do k = kfmin,kfmax
      do j = jfmin,jfmax
      do i = ifmin,ifmax
         
         unk1(ivar,i,j,k,idest) = 0.

         if (ndim == 3) then
               
            do kpar = kmina(k,order),kmaxa(k,order)
            tempy = 0.
            do jpar = jmina(j,order),jmaxa(j,order)
            tempx = 0.
            do ipar = imina(i,order),imaxa(i,order)
               tempx = tempx +
     .              weightx(ipar,i,order)*
     .              recv(ivar,ipar+offi,jpar+offj,kpar+offk)
            end do
               tempy = tempy +
     .              weighty(jpar,j,order)*tempx
            end do
               unk1(ivar,i,j,k,idest) = unk1(ivar,i,j,k,idest) + 
     .              weightz(kpar,k,order)*tempy
            end do

         elseif (ndim == 2) then

            kpar = 1
            do jpar = jmina(j,order),jmaxa(j,order)
            do ipar = imina(i,order),imaxa(i,order)
               unk1(ivar,i,j,k,idest) = unk1(ivar,i,j,k,idest) + 
     .              weightx(ipar,i,order)*
     .              weighty(jpar,j,order)*
     .              recv(ivar,ipar+offi,jpar+offj,kpar+offk)
            end do
            end do

         elseif (ndim == 1) then

            kpar = 1
            jpar = 1
            do ipar = imina(i,order),imaxa(i,order)
               unk1(ivar,i,j,k,idest) = unk1(ivar,i,j,k,idest) + 
     .              weightx(ipar,i,order)*
     .              recv(ivar,ipar+offi,jpar+offj,kpar+offk)
            end do

         end if

      end do                    ! end loop over i
      end do                    ! end loop over j
      end do                    ! end loop over k


 2    return
      end subroutine amr_1blk_cc_prol_genorder

