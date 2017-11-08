!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"
        subroutine amr_set_runtime_parameters()

        use paramesh_dimensions

        implicit none

#ifdef LIBRARY
        open (unit=35,
     &        file='amr_runtime_parameters', 
     &        status='old',
     &        action='READ',
     &        form='formatted')
        read (35,*) maxblocks
        read (35,*) ndim
        read (35,*) l2p5d
        read (35,*) nxb
        read (35,*) nyb
        read (35,*) nzb
        read (35,*) nvar
        read (35,*) nfacevar
        read (35,*) nvaredge
        read (35,*) nvarcorn
        read (35,*) nvar_work
        close(35)
#endif

        return
        end subroutine amr_set_runtime_parameters
