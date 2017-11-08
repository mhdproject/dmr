!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"

      subroutine amr_1blk_nc_prol_linear (recv,ia,ib,ja,jb,ka,kb,idest,
     .       ioff,joff,koff,mype,ivar)



      use paramesh_dimensions
      use physicaldata
      use tree
      use prolong_arrays

      implicit none

!------------------------------------

      integer, intent(in)    :: ia,ib,ja,jb,ka,kb,idest
      integer, intent(in)    :: ioff,joff,koff,mype
      integer, intent(in)    :: ivar
      real,    intent(inout) :: recv(:,:,:,:)

      integer :: icl,icu,jcl,jcu,kcl,kcu,i_ind,j_ind,k_ind
      integer :: i,j,k,i1,j1,k1,i1p,j1p,k1p
      
      real :: dx,dy,dz,cx,cy,cz

!------------------------------------

      if(prol_init.ne.100) then
       write(*,*) 'PARAMESH ERROR !'
       write(*,*) 'Error : prolong_face_fun. ',
     .       'You must call amr_prolong_face_fun_init ',
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


      i_ind = 1
      j_ind = 1
      k_ind = 1
      if(ioff.gt.0) i_ind = 2
      if(joff.gt.0) j_ind = 2
      if(koff.gt.0) k_ind = 2


! Interpolation loop.
!
! Note that the range of indeces used in the facevar plane differs
! depending on the value of iface_off. This assumes that the face values
! corresponding to index 0 (ie nguard faces to the left of the block
! boundary) are never needed, when iface_off=-1. 

      do k=kcl,kcu+iface_off
         k1 = prol_f_indexz(1,k,k_ind)
         k1p= prol_f_indexz(2,k,k_ind)
         dz = prol_f_dz(k)
         cz = 1.-dz
         do j=jcl,jcu+iface_off
            j1 = prol_f_indexy(1,j,j_ind)
            j1p= prol_f_indexy(2,j,j_ind)
            dy = prol_f_dy(j)
            cy = 1.-dy
            do i=icl,icu+iface_off
               i1 = prol_f_indexx(1,i,i_ind)
               i1p= prol_f_indexx(2,i,i_ind)
               dx = prol_f_dx(i)
               cx = 1.-dx

! compute interpolated values at location (i,j,k)
               unk_n1(ivar,i,j,k,idest) =
     .              dz*( dy*( dx*recv(ivar,i1,j1,k1) +
     .              cx*recv(ivar,i1p,j1,k1))  +
     .              cy*( dx*recv(ivar,i1,j1p,k1) +
     .              cx*recv(ivar,i1p,j1p,k1) ) ) +
     .              cz*( dy*( dx*recv(ivar,i1,j1,k1p) +
     .              cx*recv(ivar,i1p,j1,k1p))  +
     .              cy*( dx*recv(ivar,i1,j1p,k1p) +
     .              cx*recv(ivar,i1p,j1p,k1p) ) )
               
            enddo
         enddo
      enddo


      return
      end subroutine amr_1blk_nc_prol_linear
