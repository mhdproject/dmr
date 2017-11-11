!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

      subroutine amr_flux_conserve_udt(mype)




!------------------------------------------------------------------------
!
! This routine gets block boundary data from neighbors who are
! parents of leaf blocks. This is required in flux conserving schemes
! where the coarser block needs to use the same fluxes and mean pressures
! as will be used on the finer blocks across their shared boundary.
!
! The data structure used to store and pass this data is defined
! in the include file 'block_boundary_data.h' which can be included
! in 'physicaldata.h'.
!
! This version is called when uniform timesteps are being used across
! the blocks in the computation.
!
!
! Written :     Peter MacNeice          February 1997
!------------------------------------------------------------------------


      use paramesh_dimensions
      use physicaldata
      use tree

      use paramesh_interfaces, only : amr_restrict_bnd_data


#include "amr_shmem.fh"

      integer, intent(in)  ::  mype

!------------------------------------
! local variables

      integer :: lb,jf
      integer remote_pe,remote_block
      integer cnodetype
      save    cnodetype

!------------------------------------


! all leaf blocks provide reduced boundary data to their parents
      call amr_restrict_bnd_data(mype)

      call shmem_barrier_all()


! cycle through the grid blocks on this processor
      if(lnblocks.gt.0) then
      do lb = 1,lnblocks

! Is this a leaf block and not at the original refinement level ?
!      if(nodetype(lb).eq.1.and.lrefine(lb).gt.1) then
      if(nodetype(lb).eq.1) then


! Cycle over the blocks faces
       do jf = 1,nfaces

          remote_pe = neigh(2,jf,lb)
          remote_block  = neigh(1,jf,lb)

! Is the neighbor to this face a parent of a leaf block?
          cnodetype = 0
          if(remote_block.gt.0) call shmem_integer_get(cnodetype,
     .                      nodetype(remote_block),1,remote_pe)

          if(cnodetype.eq.2) then


! If yes then copy the appropriate layer from its boundary variable data 

            if(jf.eq.1) then
               call shmem_real_get(recvarxf(1,1,1,1),
     .               flux_x(1,1,1,1,remote_block),
     .               len_block_bndx*nfluxes,remote_pe)
               flux_x(1:nfluxes,1,:,:,lb) = recvarxf(1:nfluxes,2,:,:)
            elseif(jf.eq.2) then
               call shmem_real_get(recvarxf(1,1,1,1),
     .               flux_x(1,1,1,1,remote_block),
     .               len_block_bndx*nfluxes,remote_pe)
               flux_x(1:nfluxes,2,:,:,lb) = recvarxf(1:nfluxes,1,:,:)
            elseif(jf.eq.3) then
               call shmem_real_get(recvaryf(1,1,1,1),
     .               flux_y(1,1,1,1,remote_block),
     .               len_block_bndy*nfluxes,remote_pe)
               flux_y(1:nfluxes,:,1,:,lb) = recvaryf(1:nfluxes,:,2,:)
            elseif(jf.eq.4) then
               call shmem_real_get(recvaryf(1,1,1,1),
     .               flux_y(1,1,1,1,remote_block),
     .               len_block_bndy*nfluxes,remote_pe)
               flux_y(1:nfluxes,:,2,:,lb) = recvaryf(1:nfluxes,:,1,:)
            elseif(jf.eq.5) then
               call shmem_real_get(recvarzf(1,1,1,1),
     .               flux_z(1,1,1,1,remote_block),
     .               len_block_bndz*nfluxes,remote_pe)
               flux_z(1:nfluxes,:,:,1,lb) = recvarzf(1:nfluxes,:,:,2)
            elseif(jf.eq.6) then
               call shmem_real_get(recvarzf(1,1,1,1),
     .               flux_z(1,1,1,1,remote_block),
     .               len_block_bndz*nfluxes,remote_pe)
               flux_z(1:nfluxes,:,:,2,lb) = recvarzf(1:nfluxes,:,:,1)

            endif


          endif

       enddo

      endif
      enddo
      endif

      call shmem_barrier_all()

      return
      end subroutine amr_flux_conserve_udt
