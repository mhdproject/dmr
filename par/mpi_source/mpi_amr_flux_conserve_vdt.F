!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"

      subroutine amr_flux_conserve_vdt(mype,nsub)

! NOTE - THIS ROUTINE HAS NOT BEEN TESTED - July 14,2000



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

#ifdef VAR_DT

      use paramesh_interfaces, only : amr_restrict_bnd_data_vdt
#ifdef MPI_USED
      use paramesh_mpi_interfaces, only : mpi_amr_comm_setup
#endif

      implicit none

#include "amr_shmem.fh"

      integer, intent(in)  ::  mype,nsub

!------------------------------------
! local variables

      integer :: remote_pe,remote_block
      integer,save ::  anodetype(1)
      integer :: cnodetype,cchild(2),ich
      logical :: lnodetime,lfound
      save       lnodetime,cchild,cnodetype

      integer :: tag_offset,nprocs
      integer, external :: shmem_n_pes

      logical :: lcc ,lfc,lec,lnc
      logical :: lguard,lprolong,lflux,ledge,lrestrict,lfulltree
      integer :: iopt
      integer :: ierrorcode,ierr

!------------------------------------

      lcc = .false.
      lfc = .false.
      lec = .false.
      lnc = .false.
      iopt = 1


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
!------------------------------------

#ifdef MPI_USED
      nprocs = shmem_n_pes()
      tag_offset = 100


! Note, both lflux and lrestrict are true so that the fluxes
! are acquired which are needed in the restriction operation.
      lguard    = .false.
      lprolong  = .false.
      lflux     = .true.
      ledge     = .false.
      lrestrict = .true.
      lrestrict = .false.
      call mpi_amr_comm_setup(mype,nprocs,lguard,lprolong,
     .                        lflux,ledge,lrestrict,lfulltree,
     .                        iopt,lcc,lfc,lec,lnc,tag_offset)


#endif /* MPI_USED */

! Leaf blocks which have completed their timestep provide reduced 
! boundary data to their parents.
! Fluxes are accumulated in the ttflux_ arrays.
      call amr_restrict_bnd_data_vdt(mype)

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

#ifdef MPI_USED
! if (remote_block,remote_pe) is not a local block then it must have a
! local copy available in the buffer space at the end of the local
! block list.
          if(remote_pe.ne.mype) then
            do iblk = strt_buffer,last_buffer
#ifdef DEBUG
             write(*,*) 'pe ',mype,' searching buffer for ',
     .            remote_block,remote_pe,' current buffer entry ',
     .          ' iblk ',iblk,' laddress ',laddress(:,iblk)
#endif /* DEBUG */
              if(remote_block.eq.laddress(1,iblk).and.
     .             remote_pe .eq.laddress(2,iblk) ) then
                remote_block = iblk
                remote_pe    = mype
#ifdef DEBUG
             write(*,*) 'pe ',mype,' remote block ',
     .          remote_block,remote_pe,' located in buffer slot ',
     .          iblk
#endif /* DEBUG */
              endif
            enddo
          endif
#endif /* MPI_USED */

            call shmem_integer_get(anodetype,nodetype(remote_block),
     .                             1,remote_pe)
            cnodetype = anodetype(1)
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

#ifdef MPI_USED
      tag_offset = 100


      lguard    = .false.
      lprolong  = .false.
      lflux     = .true.
      ledge     = .false.
      lrestrict = .false.
      lfulltree = .false.
      call mpi_amr_comm_setup(mype,nprocs,lguard,lprolong,
     .                        lflux,ledge,lrestrict,lfulltree,
     .                        iopt,lcc,lfc,lec,lnc,tag_offset)

#endif /* MPI_USED */

!------------------------------------

! cycle through the grid blocks on this processor
      if(lnblocks.gt.0) then
      do lb = 1,lnblocks

! Is this a leaf block and not at the original refinement level ?
!      if(nodetype(lb).eq.1.and.lrefine(lb).gt.1) then
      if(nodetype(lb).eq.1) then

! Has this block completed its timestep?
      if(ldtcomplete(lb)) then

! Cycle over the blocks faces
       do jf = 1,nfaces

          remote_pe = neigh(2,jf,lb)
          remote_block  = neigh(1,jf,lb)
          cnodetype = 0
          lfound = .false.


          if(remote_block.gt.0) then

#ifdef MPI_USED
! (remote_block,remote_pe) may be a local block, a remote block,
! or it may not exist.
! If it is a local block then check its nodetype.
! If it is found in the list of remote blocks stored in buffer space
! then check its nodetype.
! If it is not found in either of these places, then set its nodetype
! to 0.
          if(remote_pe.ne.mype) then
            do iblk = strt_buffer,last_buffer
              if(remote_block.eq.laddress(1,iblk).and.
     .             remote_pe .eq.laddress(2,iblk) ) then
                remote_block = iblk
                remote_pe    = mype
                lfound = .true.
              endif
            enddo
          elseif(remote_pe.eq.mype) then
            lfound = .true.
          endif
#else  /* MPI_USED */
          lfound = .true.
#endif /* MPI_USED */


! Is the neighbor to this face a parent of a leaf block?
          if(lfound) then
             call shmem_integer_get(anodetype,
     .                       nodetype(remote_block),1,remote_pe)
             cnodetype = anodetype(1)
          endif

          endif  !  end of remote_block if test

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

            endif      ! end of jf if


          endif        ! end of cnodetype if

       enddo

      endif                      ! end of ldtcomplete if test

      endif
      enddo
      endif

!------------------------------------

#endif /* VAR_DT */

      return
      end subroutine amr_flux_conserve_vdt
