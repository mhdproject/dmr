!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"


      subroutine amr_flux_conserve(mype,nsub,flux_dir)




!------------------------------------------------------------------------
! This is a wrapper routine which makes the appropriate call to the
! routines which manage flux conservation at the boundaries between
! grid blocks of different refinement level.
!
! 
! These routines get block boundary data from neighbors who are
! parents of leaf blocks. This is required in flux conserving schemes
! where the coarser block needs to use the same fluxes and mean pressures
! as will be used on the finer blocks across their shared boundary.
!
! The data structure used to store and pass this data is defined
! in the include file 'block_boundary_data.fh' which can be included
! in 'physicaldata.fh'.
!
!
! Written :     Peter MacNeice          February 1997
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree

      use paramesh_interfaces, only : amr_flux_conserve_udt,
     &                                amr_flux_conserve_vdt


      integer, optional, intent(in) :: flux_dir
      integer, intent(in)  ::  mype,nsub

      integer :: lb

!------------------------------------


      if(lnblocks.gt.0) then
      do lb = 1,lnblocks
#ifndef ADVANCE_ALL_LEVELS
      if(nodetype(lb).eq.1) then
#endif

! Store fluxes in temporary storage
       tflux_x(:,:,:,:,lb) = flux_x(:,:,:,:,lb)
       tflux_y(:,:,:,:,lb) = flux_y(:,:,:,:,lb)
       tflux_z(:,:,:,:,lb) = flux_z(:,:,:,:,lb)


#ifndef ADVANCE_ALL_LEVELS
      endif
#endif
      enddo
      endif

      call shmem_barrier_all()

#ifdef VAR_DT
      call amr_flux_conserve_vdt(mype,nsub) ! called if variable dt
#else
      call amr_flux_conserve_udt(mype)      ! called if uniform dt
#endif

      return
      end subroutine amr_flux_conserve
