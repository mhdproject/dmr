!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

      subroutine amr_1blk_cc_prol_work_inject(recv,
     .       ia,ib,ja,jb,ka,kb,
     .       idest,ioff,joff,koff,mype)


!------------------------------------------------------------------------
!
! This routine takes data from the array recvw1, originally extracted 
! from the workspace array work on some block, 
! and performs a prolongation operation on it, between the bounds ranges 
! ia to ib, ja to jb, and ka to kb. The data in recvw1 is from a parent 
! block and the result of the prolongation operation is returned in
! the working block `work' array work1.
! The position of the child within the 
! parent block is specified by the ioff, joff and koff arguments.
!
! This particular prolongation is simple injection.
!
! Note: before using this routine in your program, make sure that the
! routine prolong_fun_init has been called.
!
!
! Written :     Peter MacNeice          January 1997
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace
      use prolong_arrays

      implicit none

      integer, intent(in) :: ia,ib,ja,jb,ka,kb
      integer, intent(in) :: idest,ioff,joff,koff,mype
      real,    intent(inout) :: recv(:,:,:)

      real :: dx,dy,dz,cx,cy,cz
      integer :: icl,icu,jcl,jcu,kcl,kcu
      integer :: i,j,k,i1,j1,k1,i1p,j1p,k1p
      integer :: offi,offj,offk

      integer,parameter :: largei = 100

!------------------------------------

      if(prolw_init.ne.100) then
       write(*,*) 'PARAMESH ERROR !'
       write(*,*) 'Error : prolong_work_fun. ',
     .       'You must call amr_prolong_fun_init ',
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
             k1 = ((k-nguard_work-1+largei)/2 + 
     .                nguard_work - largei/2 )*k3d + 1 +offk
             k1p= k1
             dz = 1.
             cz = 0.
             do j=jcl,jcu
                   j1 = ((j-nguard_work-1+largei)/2 + 
     .                      nguard_work - largei/2 )*k2d + 1 + offj
                   j1p= j1
                   dy = 1.
                   cy = 0.
                   do i=icl,icu
                         i1 = (i-nguard_work-1+largei)/2 + 
     .                           nguard_work - largei/2 + 1 + offi
                         i1p = i1
                         dx = 1.
                         cx = 0.

! compute interpolated values at location (i,j,k)
                             work1(i,j,k,idest) =
     .                          dz*( dy*( dx*recv(i1,j1,k1) +
     .                          cx*recv(i1p,j1,k1))  +
     .                          cy*( dx*recv(i1,j1p,k1) +
     .                          cx*recv(i1p,j1p,k1) ) ) +
     .                          cz*( dy*( dx*recv(i1,j1,k1p) +
     .                          cx*recv(i1p,j1,k1p))  +
     .                          cy*( dx*recv(i1,j1p,k1p) +
     .                          cx*recv(i1p,j1p,k1p) ) )



                    enddo
             enddo
        enddo


      return
      end subroutine amr_1blk_cc_prol_work_inject
