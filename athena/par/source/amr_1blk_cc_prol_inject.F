!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"


      subroutine amr_1blk_cc_prol_inject
     .  (recv,ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff,
     .   mype,ivar)


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
! This particular prolongation is simple injection.
!
! It is applied to all UNK variables whose corresponding element
! of interp_mask is set to 0.
!
! Conservative prolongation. Special treatment for the  cells immediately
! adjacent to a boundary (ie i=nguard,nguard+1,iu_bnd1-nguard,iu_bnd1-nguard+1
! and likewise for j and k indeces) if using an even number of grid cells
! per block along that axis. No special treatment is required when the number
! of cells is odd.
!
! Note: before using this routine in your program, make sure that the
! routine prolong_unk_fun_init has been called.
!
!
! Written :     Peter MacNeice          January 1997
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use prolong_arrays

      implicit none

      integer, intent(in) :: ia,ib,ja,jb,ka,kb
      integer, intent(in) :: idest,ioff,joff,koff,mype
      integer, intent(in) :: ivar
      real,    intent(inout) :: recv(:,:,:,:)

!------------------------------------
! local arrays

      real :: dx,dy,dz,cx,cy,cz
      integer :: icl,icu,jcl,jcu,kcl,kcu
      integer :: i,j,k,i1,j1,k1,i1p,j1p,k1p
      integer :: offi,offj,offk

      integer,parameter :: largei = 100

!------------------------------------


      if(prol_init.ne.100) then
       write(*,*) 'PARAMESH ERROR !'
       write(*,*) 'Error : prolong_gen_unk_fun. ',
     .       'You must call amr_initialize ',
     .       'before you can use this routine!'
       call amr_abort
      endif


! Set the bounds on the loop controlling the interpolation.
        icl=ia
        icu=ib
        jcl=ja
        jcu=jb
        kcl=ka
        kcu=kb


        offi = 0
        offj = 0
        offk = 0
        if(ioff.gt.0) offi = nxb/2
        if(joff.gt.0) offj = nyb*k2d/2
        if(koff.gt.0) offk = nzb*k3d/2

! Interpolation loop.

        do k=kcl,kcu
             k1 = ((k-nguard-1+largei)/2 + 
     .                nguard - largei/2 )*k3d + 1 +offk
             k1p= k1
             dz = 1.
             cz = 0.
             do j=jcl,jcu
                   j1 = ((j-nguard-1+largei)/2 + 
     .                      nguard - largei/2 )*k2d + 1 + offj
                   j1p= j1
                   dy = 1.
                   cy = 0.
                   do i=icl,icu
                         i1 = (i-nguard-1+largei)/2 + 
     .                           nguard - largei/2 + 1 + offi
                         i1p = i1
                         dx = 1.
                         cx = 0.

! compute interpolated values at location (i,j,k)
                             unk1(ivar,i,j,k,idest) =
     .                          dz*( dy*( dx*recv(ivar,i1,j1,k1) +
     .                          cx*recv(ivar,i1p,j1,k1))  +
     .                          cy*( dx*recv(ivar,i1,j1p,k1) +
     .                          cx*recv(ivar,i1p,j1p,k1) ) ) +
     .                          cz*( dy*( dx*recv(ivar,i1,j1,k1p) +
     .                          cx*recv(ivar,i1p,j1,k1p))  +
     .                          cy*( dx*recv(ivar,i1,j1p,k1p) +
     .                          cx*recv(ivar,i1p,j1p,k1p) ) )

                    enddo
             enddo
        enddo


      return
      end subroutine amr_1blk_cc_prol_inject
