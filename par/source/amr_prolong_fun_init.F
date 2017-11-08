!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

      subroutine amr_prolong_fun_init





      use paramesh_interfaces, only : amr_prolong_cc_fun_init,
     &                                amr_prolong_face_fun_init

!------------------------------------------------------------------------
!
! This routine calls the routines which compute the values of dx,dy and 
! dz and some index vectors used during the interpolation process. 
! These are used inside the prolongation routines
! saving needless repetitive computation at the cost of minimal storage
! space.
!
! Written :     Peter MacNeice          July 1997
!------------------------------------------------------------------------

      implicit none

      call amr_prolong_cc_fun_init

      call amr_prolong_face_fun_init


      call shmem_barrier_all()
      return
      end
