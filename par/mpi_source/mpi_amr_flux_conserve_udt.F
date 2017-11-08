!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"

      subroutine amr_flux_conserve_udt(mype,flux_dir)

!#define DEBUG



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
#ifdef MPI_USED
      use paramesh_mpi_interfaces, only : mpi_amr_comm_setup
#endif

      implicit none

#include "amr_shmem.fh"

      integer, optional, intent(in)  ::  flux_dir
      integer, intent(in)  ::  mype

!------------------------------------
! local variables

      integer :: remote_pe,remote_block
      integer,save :: anodetype(1)
      integer,save :: cnodetype

      integer :: tag_offset,nprocs
      integer, external :: shmem_n_pes
      integer :: iopt, lb, jf, iblk
      integer :: face_min, face_max, flux_dirt

      logical :: lfound
      logical :: lcc ,lfc,lec,lnc
      logical :: lguard,lprolong,lflux,ledge,lrestrict,lfulltree

!------------------------------------

      if (present(flux_dir)) then
         flux_dirt = flux_dir
      else
         flux_dirt = 0
      endif

      if (flux_dirt == 1) then
         face_min = 1
         face_max = 2
      elseif (flux_dirt == 2) then
         face_min = 3
         face_max = 4
      elseif (flux_dirt == 3) then
         face_min = 5
         face_max = 6
      else
         face_min = 1
         face_max = nfaces
      end if

      lcc = .false.
      lfc = .false.
      lec = .false.
      lnc = .false.
      iopt = 1

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
      lfulltree = .false.
      call mpi_amr_comm_setup(mype,nprocs,lguard,lprolong,
     .                        lflux,ledge,lrestrict,lfulltree,
     .                        iopt,lcc,lfc,lec,lnc,tag_offset,
     .                        flux_dir=flux_dirt)


#endif /* MPI_USED */

! all leaf blocks provide reduced boundary data to their parents
      call amr_restrict_bnd_data(mype,flux_dirt)

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
     .                        iopt,lcc,lfc,lec,lnc,tag_offset,
     .                        flux_dir=flux_dirt)

#endif /* MPI_USED */

! cycle through the grid blocks on this processor
      if(lnblocks.gt.0) then
      do lb = 1,lnblocks

! Is this a leaf block and not at the original refinement level ?
!      if(nodetype(lb).eq.1.and.lrefine(lb).gt.1) then
      if(nodetype(lb).eq.1) then


! Cycle over the blocks faces
       do jf = face_min, face_max

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
     .                      nodetype(remote_block),1,remote_pe)
            cnodetype = anodetype(1)
          endif

          endif


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

      return
      end subroutine amr_flux_conserve_udt
