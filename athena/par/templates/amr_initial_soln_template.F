
#include "paramesh_preprocessor.fh"



      subroutine amr_initial_soln





!
! This file is a template describing how the solution can be
! initialized on the initial grid. Modify it for your own use.
!
!--------------------------------------------------------------
! include files for amr
      use paramesh_dimensions
      use physicaldata
      use tree

#include "amr_shmem.fh"
      include 'shmem_reduce.fh'


      integer,parameter ::  nguard0 = nguard*npgs

!--------------------------------------------------------------


! loop over leaf grid blocks
      if(lnblocks.gt.0) then
      do lb=1,lnblocks
#ifndef ADVANCE_ALL_LEVELS
      if(nodetype(lb).eq.1) then
#endif


        if(nvar.gt.0) then

! set values for unk
        do k=kl_bnd+nguard0*k3d,ku_bnd-nguard0*k3d
          do j=jl_bnd+nguard0*k2d,ju_bnd-nguard0*k2d
            do i=il_bnd+nguard0,iu_bnd-nguard0
              unk(1,i,j,k,lb) = ???
              unk(2,i,j,k,lb) = ???
              .
              .
              .
            enddo
          enddo
        enddo

        endif

#if N_VAR_CORN > 0
        if(nvarcorn.gt.0) then

! set values for unk_n
        do k=kl_bnd+nguard0*k3d,ku_bnd+(nguard0+1)*k3d
          do j=jl_bnd+nguard0*k2d,ju_bnd+(nguard0+1)*k2d
            do i=il_bnd+nguard0,iu_bnd+(nguard0+1)
              unk_n(1,i,j,k,lb) = ???
              unk_n(2,i,j,k,lb) = ???
              .
              .
              .
            enddo
          enddo
        enddo

        endif

#endif /* N_VAR_CORN > 0 */

#if N_FACEVAR > 0

        if(nfacevar.gt.0) then

! set values for facevarx
        do k=kl_bnd+nguard0*k3d,ku_bnd+nguard0*k3d
          do j=jl_bnd+nguard0*k2d,ju_bnd+nguard0*k2d
            do i=il_bnd+nguard0,iu_bnd+nguard0+1
              facevarx(1,i,j,k,lb) = ???
              facevarx(2,i,j,k,lb) = ???
              .
              .
              .
            enddo
          enddo
        enddo

! set values for facevary
        do k=kl_bnd+nguard0*k3d,ku_bnd+nguard0*k3d
          do j=jl_bnd+nguard0*k2d,ju_bnd+(nguard0+1)*k2d
            do i=il_bnd+nguard0,iu_bnd+nguard0
              facevary(1,i,j,k,lb) = ???
              facevary(2,i,j,k,lb) = ???
              .
              .
              .
            enddo
          enddo
        enddo

! set values for facevarz
        do k=kl_bnd+nguard0*k3d,ku_bnd+(nguard0+1)*k3d
          do j=jl_bnd+nguard0*k2d,ju_bnd+nguard0*k2d
            do i=il_bnd+nguard0,iu_bnd+nguard0
              facevarz(1,i,j,k,lb) = ???
              facevarz(2,i,j,k,lb) = ???
              .
              .
              .
            enddo
          enddo
        enddo

        endif


#endif /* N_FACEVAR > 0 */

#if N_VAR_EDGE > 0

        if(nvaredge.gt.0) then

! set values for unk_e_x
        do k=kl_bnd+nguard0*k3d,ku_bnd+(nguard0+1)*k3d
          do j=jl_bnd+nguard0*k2d,ju_bnd+(nguard0+1)*k2d
            do i=il_bnd+nguard0,iu_bnd+nguard0
              unk_e_x(1,i,j,k,lb) = ???
              unk_e_x(2,i,j,k,lb) = ???
              .
              .
              .
            enddo
          enddo
        enddo

! set values for unk_e_y
        do k=kl_bnd+nguard0*k3d,ku_bnd+(nguard0+1)*k3d
          do j=jl_bnd+nguard0*k2d,ju_bnd+nguard0*k2d
            do i=il_bnd+nguard0,iu_bnd+(nguard0+1)
              unk_e_y(1,i,j,k,lb) = ???
              unk_e_y(2,i,j,k,lb) = ???
              .
              .
              .
            enddo
          enddo
        enddo

! set values for unk_e_z
        do k=kl_bnd+nguard0*k3d,ku_bnd+nguard0*k3d
          do j=jl_bnd+nguard0*k2d,ju_bnd+(nguard0+1)*k2d
            do i=il_bnd+nguard0,iu_bnd+(nguard0+1)
              unk_e_z(1,i,j,k,lb) = ???
              unk_e_z(2,i,j,k,lb) = ???
              .
              .
              .
            enddo
          enddo
        enddo

        endif

#endif /* N_VAREDGE > 0 */

#ifndef ADVANCE_ALL_LEVELS
      endif
#endif
      enddo ! end loop over grid blocks
      endif

      return
      end subroutine amr_initial_soln
