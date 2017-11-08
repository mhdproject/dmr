!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

      subroutine amr_flux_conserve_vdt(mype,nsub)




!------------------------------------------------------------------------
!
! This routine gets block boundary data from neighbors who are
! parents of leaf blocks. This is required in flux conserving schemes
! where the coarser block needs to use the same fluxes and mean pressures
! as will be used on the finer blocks across their shared boundary.
!
! The data structure used to store and pass this data is defined
! in the include file 'block_boundary_data.fh' which can be included
! in 'physicaldata.fh'.
!
! This version is used when variable timesteps are allowed across the
! blocks in the computation.
!
! Arguments:
!      mype          processor number
!      nsub          current time subcycle. If this is 1 then this
!                     info is used to reset the temporary boundary flux
!                     arrays to 0.
!
! Written :     Peter MacNeice          February 1997
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree

      use paramesh_interfaces, only : amr_restrict_bnd_data_vdt

#include "amr_shmem.fh"

      integer, intent(in)  ::  mype,nsub

!------------------------------------
! local variables

      integer remote_pe,remote_block
      integer cnodetype,cchild(2),ich
      logical lnodetime
      save    lnodetime,cchild,cnodetype

!------------------------------------
#ifdef VAR_DT


      if(lnblocks.gt.0) then
      do lb = 1,lnblocks

! Is this a parent of at least one leaf block ?
      if(nodetype(lb).eq.2) then

! Set timestep phases for the current block, and for the next finer level.
        lcycle = loc_cycle(lrefine(lb))
        phase0 = phase_dt(lrefine(lb))
        phase1 = phase_dt(lrefine(lb)+1)

! At start of the current blocks timestep zero out the arrays used to 
! accumulate boundary fluxes from its children.
!        if(phase0.eq.1.and.phase1.eq.1) then
        if(lcycle.eq.1) then
           ttflux_x(:,:,:,:,lb) = 0.
           if(ndim.ge.2) ttflux_y(:,:,:,:,lb) = 0.
           if(ndim.eq.3) ttflux_z(:,:,:,:,lb) = 0.
        endif

      endif
      enddo
      endif
      call shmem_barrier_all()
!------------------------------------

! Leaf blocks which have completed their timestep provide reduced 
! boundary data to their parents.
! Fluxes are accumulated in the ttflux_ arrays.
      call amr_restrict_bnd_data_vdt(mype)


      call shmem_barrier_all()
!------------------------------------


! Parents who have completed their timestep and border a leaf block
! update their fluxes.
      do lb = 1,lnblocks


! Is this a parent block of at least one leaf node?
      if((nodetype(lb).eq.2).and.ldtcomplete(lb)) then

! If yes then cycle through its neighbors.
        do iface=1,nfaces

! If this neighbor is a leaf block or an external boundary then 
! replace fluxes with restricted fluxes.
          cnodetype = 1
          if(neigh(1,iface,lb).ge.1) then
            remote_pe    = neigh(2,iface,lb)
            remote_block = neigh(1,iface,lb)

            call shmem_integer_get(cnodetype,nodetype(remote_block),
     .                             1,remote_pe)
          endif
          if(cnodetype.eq.1) then
            if(iface.eq.1) flux_x(:,1,:,:,lb)=ttflux_x(:,1,:,:,lb)
            if(iface.eq.2) flux_x(:,2,:,:,lb)=ttflux_x(:,2,:,:,lb)
            if(iface.eq.3) flux_y(:,:,1,:,lb)=ttflux_y(:,:,1,:,lb)
            if(iface.eq.4) flux_y(:,:,2,:,lb)=ttflux_y(:,:,2,:,lb)
            if(iface.eq.5) flux_z(:,:,:,1,lb)=ttflux_z(:,:,:,1,lb)
            if(iface.eq.6) flux_z(:,:,:,2,lb)=ttflux_z(:,:,:,2,lb)
          endif
        enddo
      endif
      enddo


      call shmem_barrier_all()
!------------------------------------

! cycle through the grid blocks on this processor
      if(lnblocks.gt.0) then
      do lb = 1,lnblocks

! Is this a leaf block and not at the original refinement level ?
!      if((nodetype(lb).eq.1).and.(lrefine(lb).gt.1)) then
      if(nodetype(lb).eq.1) then

! Has this block completed its timestep?
      if(ldtcomplete(lb)) then

! Cycle over the blocks faces
       do jf = 1,nfaces

          remote_pe = neigh(2,jf,lb)
          remote_block  = neigh(1,jf,lb)

! Is the neighbor to this face a parent of a leaf block?
          cnodetype = 0
          if(remote_block.gt.0) then
             call shmem_integer_get(cnodetype,
     .                       nodetype(remote_block),1,remote_pe)
          endif


          if(cnodetype.eq.2) then


             if(jf.eq.1) then
               call shmem_real_get(recvarxf(1,1,1,1),
     .               flux_x(1,1,1,1,remote_block),
     .               len_block_bndx*nfluxes,remote_pe)
                 flux_x(1:nfluxes,1,:,:,lb) = 
     .                           recvarxf(1:nfluxes,2,:,:)

             elseif(jf.eq.2) then
               call shmem_real_get(recvarxf(1,1,1,1),
     .               flux_x(1,1,1,1,remote_block),
     .               len_block_bndx*nfluxes,remote_pe)
                 flux_x(1:nfluxes,2,:,:,lb) = 
     .                             recvarxf(1:nfluxes,1,:,:)

             elseif(jf.eq.3) then
               call shmem_real_get(recvaryf(1,1,1,1),
     .               flux_y(1,1,1,1,remote_block),
     .               len_block_bndy*nfluxes,remote_pe)
                 flux_y(1:nfluxes,:,1,:,lb) = 
     .                             recvaryf(1:nfluxes,:,2,:)

             elseif(jf.eq.4) then
               call shmem_real_get(recvaryf(1,1,1,1),
     .               flux_y(1,1,1,1,remote_block),
     .               len_block_bndy*nfluxes,remote_pe)
                 flux_y(1:nfluxes,:,2,:,lb) = 
     .                             recvaryf(1:nfluxes,:,1,:)

             elseif(jf.eq.5) then
               call shmem_real_get(recvarzf(1,1,1,1),
     .               flux_z(1,1,1,1,remote_block),
     .               len_block_bndz*nfluxes,remote_pe)
                 flux_z(1:nfluxes,:,:,1,lb)=
     .                             recvarzf(1:nfluxes,:,:,2)

             elseif(jf.eq.6) then
               call shmem_real_get(recvarzf(1,1,1,1),
     .               flux_z(1,1,1,1,remote_block),
     .               len_block_bndz*nfluxes,remote_pe)
                 flux_z(1:nfluxes,:,:,2,lb) = 
     .                             recvarzf(1:nfluxes,:,:,1)

             endif


          endif

        enddo

      endif                      ! end of ldtcomplete if test

      endif
      enddo
      endif

      call shmem_barrier_all()
!------------------------------------

#endif /* VAR_DT */

      return
      end subroutine amr_flux_conserve_vdt
