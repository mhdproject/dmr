!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"

      subroutine amr_guardcell(mype,iopt,nlayers,
     .                         nlayersx,nlayersy,nlayersz)



!------------------------------------------------------------------------
!
! This routine manages the exchange of guard cell information.
! It cycles through the blocks, applying the routine amr_1blk_guardcell
! to each block in turn.
!
! Written :     Peter MacNeice          January 1997
!------------------------------------------------------------------------
!
! Arguments:
!      mype         integer          local processor number
!      iopt         integer          indicates which data arrays are to
!                                    have their guardcells filled 
!                                      iopt = 1    unk, facevar[x][y][z]
!                                      iopt >= 2   work
!      nlayers      integer          no. of guardcell layers at each block
!                                    boundary
!
!------------------------------------

      use paramesh_dimensions
      use physicaldata
      use workspace
      use tree

      use paramesh_interfaces, only : amr_1blk_guardcell_reset,
     &                                amr_restrict,
     &                                amr_1blk_guardcell,
     &                                amr_1blk_to_perm,
     &                                amr_gsurrounding_blks


      integer, optional, intent(in) :: nlayersx, nlayersy, nlayersz
      integer, intent(in) :: mype,iopt,nlayers



      logical  lcc,lfc,lec,lnc,l_srl_only,ldiag
      integer  lb,pe,icoord


#ifdef TIMINGS
#include "timer.fh"
#endif

!------------------------------------
#ifndef DIAGONALS
      write(*,*) 'amr_guardcell:  diagonals off'
#endif

#ifdef NO_PERMANENT_GUARDCELLS

      if(mype.eq.0) then
        write(*,*) 'amr_guardcell call ignored!'
        write(*,*) 'NO_PERMANENT_GUARDCELLS is defined'
      endif

#else /* NO_PERMANENT_GUARDCELLS */


! make sure that nlayers and iopt are set consistently.
      if(iopt.eq.1.and.nlayers.ne.nguard) then
        if(mype.eq.0) then
          write(*,*) 'PARAMESH ERROR !'
          write(*,*) 'Error in guardcell - iopt and nlayers'
          write(*,*) 'are not consistent. For iopt=1 you must'
          write(*,*) 'set nlayers=nguard.'
        endif
        call amr_abort
      elseif(iopt.ge.2.and.nlayers.gt.nguard_work) then
        if(mype.eq.0) then
          write(*,*) 'PARAMESH ERROR !'
          write(*,*) 'Error in guardcell - iopt and nlayers'
          write(*,*) 'are not consistent. For iopt>=2 you must'
          write(*,*) 'set nlayers le nguard_work.'
        endif
        call amr_abort
      endif

! reinitialize addresses of cached parent blocks
      call amr_1blk_guardcell_reset

      call shmem_barrier_all()


#ifdef TIMINGS
      call shmem_barrier_all()
      itimer1 = irtc()
#endif



      lcc = .false.
      lfc = .false.
      lec = .false.
      lnc = .false.
      if(iopt.eq.1) then
        if(nvar.gt.0) lcc = .true.
        if(nfacevar.gt.0) lfc = .true.
        if(nvaredge.gt.0) lec = .true.
        if(nvarcorn.gt.0) lnc = .true.
      elseif(iopt.ge.2) then
        lcc = .true.
      endif
!
! Restrict solution to parent blocks
#ifndef ADVANCE_ALL_LEVELS

      iempty = 0
      call amr_restrict(mype,iopt,iempty)

      call amr_1blk_guardcell_reset
#endif

#ifdef FORCE_CONSISTENCY_AT_SRL_INTERFACES
      if(lfc) then
        do lb = 1,lnblocks
          gt_facevarx(:,1,:,:,lb) = facevarx(:,1+nguard,:,:,lb)
          gt_facevarx(:,2,:,:,lb) = facevarx(:,nxb+1+nguard,:,:,lb)
          if(ndim.ge.2) then
          gt_facevary(:,:,1,:,lb) = facevary(:,:,1+nguard*k2d,:,lb)
          gt_facevary(:,:,1+k2d,:,lb) = 
     .                        facevary(:,:,nyb+(1+nguard)*k2d,:,lb)
          endif
          if(ndim.eq.3) then
          gt_facevarz(:,:,:,1,lb) = facevarz(:,:,:,1+nguard*k3d,lb)
          gt_facevarz(:,:,:,1+k3d,lb) = 
     .                        facevarz(:,:,:,nzb+(1+nguard)*k3d,lb)
          endif
        enddo
      endif
      call shmem_barrier_all()
#endif  /* FORCE_CONSISTENCY_AT_SRL_INTERFACES */

        ldiag = .false.
#ifdef DIAGONALS
        ldiag = .true.
#endif

! If surr_blks has not been set up do it now.
      if(gsurrblks_set.ne.1) call amr_gsurrounding_blks(mype,ldiag)


      if(lnblocks.gt.0) then
      do lb = 1,lnblocks
#ifndef ADVANCE_ALL_LEVELS
      if(nodetype(lb).eq.1) then
#endif

! Copy this blocks data into the working block, and fill its guardcells
        l_srl_only = .false.                     ! fill srl and coarse
        icoord = 0                               ! fill in all coord directions
        call amr_1blk_guardcell(mype,iopt,nlayers,lb,mype,
     .                          lcc,lfc,lec,lnc,
     .                          l_srl_only,icoord,ldiag)


! Now copy the result back to permanent storage arrays
        idest = 1
        call amr_1blk_to_perm( lcc,lfc,lec,lnc,lb,iopt,idest )


#ifndef ADVANCE_ALL_LEVELS
      endif
#endif
      enddo                                 ! end of loop over blocks
      endif


! reinitialize addresses of cached parent blocks
      call amr_1blk_guardcell_reset

#ifdef TIMINGS
      itimer2 = irtc()
      call shmem_barrier_all()
      irtc_guard = itimer2-itimer1+irtc_guard
#endif

#endif /* NO_PERMANENT_GUARDCELLS */


      return
      end subroutine amr_guardcell
