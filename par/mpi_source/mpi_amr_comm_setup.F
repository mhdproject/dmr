#define GCELL_ON
#define GCELL_ON_FC
#define GCELL_ON_EC
#define GCELL_ON_NC

#define ORIG
!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"

!#define DEBUG

      subroutine mpi_amr_comm_setup(mype,nprocs,
     .                              lguard,lprolong,
     .                              lflux,ledge,lrestrict,lfulltree,
     .                              iopt,lcc,lfc,lec,lnc,tag_offset,
     .                              nlayersx,nlayersy,nlayersz,
     .                              flux_dir)


!------------------------------------------------------------------------
!
! This routine is the pre-analysis operation which precedes
! a code block in which interprocessor communications may
! be required.
!
!
! Written :     Peter MacNeice          June 2000
!------------------------------------------------------------------------
!
! Arguments:
!      mype           rank of local processor
!      nprocs         number of processors
!      lguard         if true, info for guardcell filling will be
!                      communicated
!      lprolong       if true, info for prolongation will be
!                      communicate
!      lflux          if true, info for flux conservation will be
!                      communicated
!      ledge          if true, info for edge averaging will be
!                      communicated
!      ledge          if true, info for restriction will be
!                      communicated
!      iopt           a switch to control which data source is to be used
!                      iopt=1 will use 'unk'
!                      iopt>=2 will use 'work'
!      lcc            a logical switch controlling whether unk or work data
!                      is filled.
!      lfc            a logical switch controlling whether facevar data
!                      is filled.
!      lec            a logical switch controlling whether unk_e_? data
!                      is filled.
!      lnc            a logical switch controlling whether unk_n data
!                      is filled.
!      tag_offset     defines the last tag number used for an mpi message.
!
!------------------------------------------------------------------------
      use paramesh_dimensions
      use physicaldata
      use workspace
      use tree
      use timings
      use mpi_morton


      use paramesh_mpi_interfaces, only :
     &                                mpi_pack_blocks,
     &                                mpi_unpack_blocks,
     &                                mpi_pack_edges,
     &                                mpi_unpack_edges,
     &                                mpi_pack_fluxes,
     &                                mpi_unpack_fluxes,
     &                                mpi_pack_tree_info,
     &                                mpi_unpack_tree_info,
     &                                mpi_amr_read_guard_comm,
     &                                mpi_amr_read_prol_comm,
     &                                mpi_amr_read_flux_comm,
     &                                mpi_amr_read_restrict_comm,
     &                                mpi_set_message_sizes,
     &                                mpi_xchange_blocks

      implicit none

      integer, intent(in)    :: mype,nprocs,iopt
      integer, intent(inout) :: tag_offset
      logical, intent(in)    :: lcc,lfc,lec,lnc,lfulltree
      logical, intent(in)    :: lguard,lprolong,lflux,ledge,lrestrict
      integer, intent(in), optional :: nlayersx,nlayersy,nlayersz
      integer, intent(in), optional :: flux_dir

      include 'mpif.h'

#ifdef ORIG
      real, dimension (:), allocatable :: send_buf
      real, dimension (:), allocatable :: recv_buf
#endif /* ORIG */
      integer :: buffer_dim,buf_dim_bytes
      integer :: bufdim
      integer :: max_blks_sent
      integer :: len_surr_blks
      integer :: itemp
      integer :: offset_tree
      integer :: nlayerstx, nlayersty, nlayerstz
      integer :: flux_dirt
      integer :: invar, ibndvar, ivaredge, ivarcorn, ierror
#ifdef TIMING_MPI
      double precision :: time1
      double precision :: time2
#endif /* TIMING_MPI */

#ifdef GCELL_ON
      integer :: ii,jj
#endif /* GCELL_ON */

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef TIMING_MPI
       time1 = mpi_wtime()
       time2 = mpi_wtime()
#endif /* TIMING_MPI */


#ifdef DEBUG
      write(*,*) 'pe ',mype,' entered mpi_amr_comm_setup: '
     .           ,' max_blks_sent ',
     .           max_blks_sent,' tag_offset ',tag_offset,
     .           ' nprocs ',nprocs
#endif /* DEBUG */

      if (iopt == 1) then


! install user selection for guardcell variables on reset defaults.
      if(lguard) then
        int_gcell_on_cc = gcell_on_cc
        int_gcell_on_fc = gcell_on_fc
        int_gcell_on_ec = gcell_on_ec
        int_gcell_on_nc = gcell_on_nc
        lguard_in_progress = .true.

#ifdef GCELL_ON
        if(nvar.gt.0) then
          ngcell_on_cc = 0
          do ii = 1,nvar
            if(int_gcell_on_cc(ii)) then
              ngcell_on_cc =  ngcell_on_cc + 1
              gcell_on_cc_pointer(ngcell_on_cc) = ii
            endif
          enddo
        endif
#else
        ngcell_on_cc = nvar
#endif /* GCELL_ON */
#ifdef GCELL_ON_FC
        if(nfacevar.gt.0) then
          ngcell_on_fc = 0
          do ii = 1,nfacevar
          do jj = 1,3
            if(int_gcell_on_fc(jj,ii)) then
              ngcell_on_fc(jj) =  ngcell_on_fc(jj) + 1
              gcell_on_fc_pointer(jj,ngcell_on_fc(jj)) = ii
            endif
          enddo
          enddo
        endif
#else
        ngcell_on_fc = nfacevar
#endif /* GCELL_ON_FC */
#ifdef GCELL_ON_EC
        if(nvaredge.gt.0) then
          ngcell_on_ec = 0
          do ii = 1,nvaredge
          do jj = 1,3
            if(int_gcell_on_ec(jj,ii)) then
              ngcell_on_ec(jj) =  ngcell_on_ec(jj) + 1
              gcell_on_ec_pointer(jj,ngcell_on_ec(jj)) = ii
            endif
          enddo
          enddo
        endif
#else
        ngcell_on_ec = nvaredge
#endif /* GCELL_ON_EC */
#ifdef GCELL_ON_NC
        if(nvarcorn.gt.0) then
          ngcell_on_nc = 0
          do ii = 1,nvarcorn
            if(int_gcell_on_nc(ii)) then
              ngcell_on_nc =  ngcell_on_nc + 1
              gcell_on_nc_pointer(ngcell_on_nc) = ii
            endif
          enddo
        endif
#else
        ngcell_on_nc = nvarcorn
#endif /* GCELL_ON_NC */

      else                                    ! lguard

        int_gcell_on_cc(:) = .true.
        int_gcell_on_fc(:,:) = .true.
        int_gcell_on_ec(:,:) = .true.
        int_gcell_on_nc(:) = .true.
        lguard_in_progress = .false.
        ngcell_on_cc = nvar
#ifdef GCELL_ON
        do ii=1,nvar
          gcell_on_cc_pointer(ii) = ii
        enddo
#endif /* GCELL_ON */
        ngcell_on_fc = nfacevar
#ifdef GCELL_ON_FC
        do ii=1,nfacevar
          gcell_on_fc_pointer(:,ii) = ii
        enddo
#endif /* GCELL_ON_FC */
        ngcell_on_ec = nvaredge
#ifdef GCELL_ON_EC
        do ii=1,nvaredge
          gcell_on_ec_pointer(:,ii) = ii
        enddo
#endif /* GCELL_ON_EC */
        ngcell_on_nc = nvarcorn
#ifdef GCELL_ON_NC
        do ii=1,nvarcorn
          gcell_on_nc_pointer(ii) = ii
        enddo
#endif /* GCELL_ON_NC */

      endif                                    ! lguard

      if (present(nlayersx)) then
      nlayerstx = nlayersx
!      write(*,*) 'nlayersx = ',nlayerstx,' present'
      else
      nlayerstx = nguard
!      write(*,*) 'nlayerstx = ',nlayerstx,' not present',
!     .       ' set to nguard ',nguard
      end if

      if (present(nlayersy)) then
      nlayersty = nlayersy
      else
      nlayersty = nguard
      end if

      if (present(nlayersz)) then
      nlayerstz = nlayersz
      else
      nlayerstz = nguard
      end if

      else

      if (present(nlayersx)) then
      nlayerstx = nlayersx
      else
      nlayerstx = nguard_work
      end if

      if (present(nlayersy)) then
      nlayersty = nlayersy
      else
      nlayersty = nguard_work
      end if

      if (present(nlayersz)) then
      nlayerstz = nlayersz
      else
      nlayerstz = nguard_work
      end if

      endif ! end of if (iopt

      if (present(flux_dir)) then
         flux_dirt = flux_dir
      else
         flux_dirt = 0
      end if

      call mpi_set_message_sizes(iopt,
     .     nlayerstx,nlayersty,nlayerstz)

#ifdef TIMING_MPI
      timer_amr_comm_setup(1) =  timer_amr_comm_setup(1)
     .                          + mpi_wtime() - time2
      time2 = mpi_wtime()
#endif /* TIMING_MPI */

      if(lguard.and.(.not.lrestrict)) then

#ifdef DEBUG
      write(*,*) 'pe ',mype,' calling mpi_amr_read_guard_comm'
#endif /* DEBUG */

         call mpi_amr_read_guard_comm(nprocs)

#ifdef DEBUG
      write(*,*) 'pe ',mype,' exited mpi_amr_read_guard_comm'
#endif /* DEBUG */

      elseif(lprolong) then

         call mpi_amr_read_prol_comm(nprocs)

      elseif((lflux.or.ledge).and.(.not.lrestrict)) then

#ifdef DEBUG
      write(*,*) 'pe ',mype,' calling mpi_amr_read_flux_comm'
#endif /* DEBUG */

        call mpi_amr_read_flux_comm(nprocs)

      elseif(lrestrict) then

#ifdef DEBUG
      write(*,*) 'pe ',mype,' calling mpi_amr_read_restrict_comm'
#endif /* DEBUG */

        call mpi_amr_read_restrict_comm(nprocs)

#ifdef DEBUG
      write(*,*) 'pe ',mype,' exited mpi_amr_read_restrict_comm'
#endif /* DEBUG */

      endif

#ifdef TIMING_MPI
      timer_amr_comm_setup(2) =  timer_amr_comm_setup(2)
     .                          + mpi_wtime() - time2
      time2 = mpi_wtime()
#endif /* TIMING_MPI */

      itemp = max(sum(commatrix_send), sum(commatrix_recv))
      call MPI_ALLREDUCE (itemp,
     .                    max_blks_sent,
     .                    1,
     .                    MPI_INTEGER,
     .                    MPI_MAX,
     .                    MPI_COMM_WORLD,
     .                    ierror)

      invar   = 0
      ibndvar = 0
      ivaredge = 0
      ivarcorn = 0

      if(lguard.or.lprolong) then

        if(iopt.eq.1) then
          invar   = 0
          ibndvar = 0
          ivaredge = 0
          ivarcorn = 0
          if (lcc) invar = nvar
          if (lfc) ibndvar = nbndvar
          if (lec) ivaredge = nvaredge
          if (lnc) ivarcorn = nvarcorn

#ifdef GCELL_ON
          if(lcc) invar = ngcell_on_cc
#endif /* GCELL_ON */
#ifdef GCELL_ON_FC
          if(lfc) ibndvar = maxval(ngcell_on_fc(1:ndim))
#endif /* GCELL_ON_FC */
#ifdef GCELL_ON_EC
          if(lec) ivaredge = maxval(ngcell_on_ec(1:3))
#endif /* GCELL_ON_EC */
#ifdef GCELL_ON_NC
          if(lnc) ivarcorn = ngcell_on_nc
#endif /* GCELL_ON_NC */
          if(ndim.eq.1) buffer_dim = invar*nxb*nyb*nzb + ibndvar*
     &                  (nxb+1)*nyb*nzb 
     &                   + ivarcorn*(nxb+1)*(nyb+k2d)*(nzb+k3d)
     &                   + 3
          if(ndim.eq.2) buffer_dim =
     &                  invar*nxb*nyb*nzb + ibndvar*(
     &                  (nxb+1)*nyb*nzb +nxb*(nyb+k2d)*nzb )
     &                   + ivaredge*( nxb*(nyb+k2d)*(nzb+k3d) +
     &                                (nxb+1)*nyb*(nzb+k3d) )
     &                   + ivarcorn*(nxb+1)*(nyb+k2d)*(nzb+k3d)
     &                 + 3
          if(ndim.eq.3.or.l2p5d.eq.1) buffer_dim =
     &                  invar*nxb*nyb*nzb + ibndvar*(
     &                  (nxb+1)*nyb*nzb +nxb*(nyb+k2d)*nzb +
     &                  nxb*nyb*(nzb+k3d))
     &                   + ivaredge*( nxb*(nyb+k2d)*(nzb+k3d) +
     &                                (nxb+1)*nyb*(nzb+k3d)  +
     &                                (nxb+1)*(nyb+k2d)*nzb )
     &                   + ivarcorn*(nxb+1)*(nyb+k2d)*(nzb+k3d)
     &                 + 3
        elseif(iopt.gt.1) then
          buffer_dim = nxb*nyb*nzb + 3
        endif

      elseif(lflux) then

        if (flux_dirt == 1) then

        buffer_dim = nfluxes*2*(nyb*nzb+0+0)
        buffer_dim = buffer_dim + 3

        elseif (flux_dirt == 2) then

        buffer_dim = nfluxes*2*(0+nxb*nzb*k2d+0)
        buffer_dim = buffer_dim + 3

        elseif (flux_dirt == 3) then

        buffer_dim = nfluxes*2*(0+0+nxb*nyb*k3d)
        buffer_dim = buffer_dim + 3

        else

        buffer_dim = nfluxes*2*(nyb*nzb+nxb*nzb*k2d+nxb*nyb*k3d)
        buffer_dim = buffer_dim + 3

        end if

      elseif(ledge) then

        buffer_dim = nedges*2*2*( (nyb+k2d)*(nzb+k3d)
     .                           +(nxb+1)  *(nzb+k3d)*k2d
     .                           +(nxb+1)  *(nyb+k2d)*k3d )
        if(l2p5d.eq.1) buffer_dim = nedges*( 
     .                          2*(  (nyb+k2d) + (nyb+k2d) +
     .                               (nxb+1)   + (nxb+1)  ) 
     .                           +(nxb+1)*(nyb+k2d)
     .                           +(nxb+1)*(nyb+k2d)   )
        buffer_dim = buffer_dim + 3

      endif
      buffer_dim = buffer_dim*max_blks_sent


! If lrestrict is true then we only need tree information.
      if( (.not.lrestrict) .or. (lrestrict.and.lflux)
     .                     .or. (lrestrict.and.ledge)
     .                     .or. (lrestrict.and.lguard)
     .                                            ) then
!----
! Set up buffer size info, including space for tree data
      len_surr_blks = 3*3*(1+2*k2d)*(1+2*k3d)
#ifdef SAVE_MORTS
     .                + 6*3*(1+2*k2d)*(1+2*k3d)
#endif
      buffer_dim = buffer_dim +
     .             (32+16+len_surr_blks)*max_blks_sent
      offset_tree = 32+16+len_surr_blks
!      buffer_dim = buffer_dim +
!     .             (32+len_surr_blks)*max_blks_sent
!      offset_tree = 32+len_surr_blks


#ifdef DEBUG
      write(*,*) 'pe ',mype,' logic switches ',
     .                 ' lrestrict,lflux,lguard,ledge,lprolong ',
     .                 lrestrict,lflux,lguard,ledge,lprolong
      write(*,*) 'pe ',mype,' setup : buffer_dim ',
     .           buffer_dim
#endif /* DEBUG */

#ifdef ORIG
      if(allocated(send_buf)) deallocate(send_buf)
#ifdef DEBUG
      write(*,*) 'pe ',mype,' setup : send_buf deallocated '
#endif /* DEBUG */
      allocate( send_buf(buffer_dim))
#ifdef DEBUG
      write(*,*) 'pe ',mype,' setup : send_buf allocated ',
     .     ' size ',buffer_dim
#endif /* DEBUG */
      if(allocated(temprecv_buf)) deallocate(temprecv_buf)
      allocate( temprecv_buf(buffer_dim))
      bufdim = buffer_dim
#endif /* ORIG */

!      write(*,*) 'pe ',mype,' allocated(temprecv_buf) ',
!     . allocated(temprecv_buf),' allocated(send_buf) ',
!     . allocated(send_buf),' size(temprecv_buf) ',
!     . size(temprecv_buf,dim=1),' buffer_dim ',buffer_dim
#ifndef ORIG
      if(allocated(temprecv_buf).and.allocated(send_buf).and.
     .    (size(temprecv_buf,dim=1).lt.buffer_dim
     .         .or.
     .     size(temprecv_buf,dim=1).gt.2*buffer_dim) ) then
        deallocate(send_buf)
        allocate( send_buf(buffer_dim))
        deallocate(temprecv_buf)
        allocate( temprecv_buf(buffer_dim))
      endif
      if(.not.allocated(send_buf)) allocate(send_buf(buffer_dim))
      if(.not.allocated(temprecv_buf)) 
     .                         allocate(temprecv_buf(buffer_dim))
      bufdim = size(send_buf,dim=1)
#endif /* ORIG */
!      write(*,*) 'size ',size(send_buf,dim=1),
!     .     ' shape ',shape(send_buf)
      

#ifdef DEBUG
      write(*,*) 'pe ',mype,' setup : allocated '
#endif /* DEBUG */

#ifdef TIMING_MPI
      timer_amr_comm_setup(3) =  timer_amr_comm_setup(3)
     .                          + mpi_wtime() - time2
      time2 = mpi_wtime()
#endif /* TIMING_MPI */

      if(lguard.or.lprolong) then



        call mpi_pack_blocks(mype,nprocs,iopt,lcc,lfc,lec,lnc,
     .                       bufdim,send_buf,offset_tree,
     .                       nlayerstx,nlayersty,nlayerstz)
!        call mpi_pack_blocks(mype,nprocs,iopt,lcc,lfc,lec,lnc,
!     .                       buffer_dim,send_buf,offset_tree,
!     .                       nlayerstx,nlayersty,nlayerstz)

      elseif(lflux) then

!        call mpi_pack_fluxes(mype,nprocs,buffer_dim,send_buf,
        call mpi_pack_fluxes(mype,nprocs,bufdim,send_buf,
     .                       offset_tree,flux_dir)

      elseif(ledge) then

!        call mpi_pack_edges(mype,nprocs,buffer_dim,send_buf,
        call mpi_pack_edges(mype,nprocs,bufdim,send_buf,
     .                      offset_tree)

      endif

#ifdef DEBUG
      write(*,*) 'pe ',mype,' setup : exited mpi_pack_blocks'
#endif /* DEBUG */

#ifdef DEBUG
      write(*,*) 'pe ',mype,' setup : entering mpi_xchange_blocks'
     .       ,' buffer_dim ',buffer_dim
#endif /* DEBUG */
#ifdef TIMING_MPI
      timer_amr_comm_setup(4) =  timer_amr_comm_setup(4)
     .                          + mpi_wtime() - time2
      time2 = mpi_wtime()
#endif /* TIMING_MPI */



      call mpi_xchange_blocks(mype,nprocs,tag_offset,
     .                        bufdim,send_buf,temprecv_buf)
!     .                        buffer_dim,send_buf,temprecv_buf)


#ifdef TIMING_MPI
      timer_amr_comm_setup(5) =  timer_amr_comm_setup(5)
     .                          + mpi_wtime() - time2
      time2 = mpi_wtime()
#endif /* TIMING_MPI */

#ifdef DEBUG
      write(*,*) 'pe ',mype,' setup : exited mpi_xchange_blocks'
#endif /* DEBUG */

      if(lguard.or.lprolong) then

         call mpi_unpack_blocks(mype,iopt,lcc,lfc,lec,lnc,
     .                          bufdim,temprecv_buf,
     .                          nlayerstx,nlayersty,nlayerstz)
!         call mpi_unpack_blocks(mype,iopt,lcc,lfc,lec,lnc,
!     .                          buffer_dim,temprecv_buf,
!     .                          nlayerstx,nlayersty,nlayerstz)

#ifdef DEBUG
      write(*,*) 'pe ',mype,' setup : unpack blocks from buffer'
#endif /* DEBUG */
      elseif(lflux) then
         
!         call mpi_unpack_fluxes(mype,buffer_dim,temprecv_buf,
         call mpi_unpack_fluxes(mype,bufdim,temprecv_buf,
     .                          flux_dir)

#ifdef DEBUG
      write(*,*) 'pe ',mype,' setup : exited mpi_unpack_fluxes'
#endif /* DEBUG */

      elseif(ledge) then

!         call mpi_unpack_edges(mype,buffer_dim,temprecv_buf)
         call mpi_unpack_edges(mype,bufdim,temprecv_buf)

#ifdef DEBUG
      write(*,*) 'pe ',mype,' setup : exited mpi_unpack_edges'
#endif /* DEBUG */

      endif

#ifdef TIMING_MPI
      timer_amr_comm_setup(6) =  timer_amr_comm_setup(6)
     .                          + mpi_wtime() - time2
#endif /* TIMING_MPI */

!      if(lguard.or.lprolong) then
! copy recv_buf into temprecv_buf because recv_buf will be reused
! to communicate tree info. Guardcell and prolongation routines will
! read directly from temprecv_buf.
!        temprecv_buf = recv_buf
!      endif

#ifdef ORIG
      if(allocated(send_buf)) deallocate(send_buf)
#endif /* ORIG */

#ifdef DEBUG
      write(*,*) 'pe ',mype,' setup : deallocated '
#endif /* DEBUG */

!----

      else

!--------------------------------------------------------------
! 
! Establish tree info which is needed about remote blocks

#ifdef TIMING_MPI
      time2 = mpi_wtime()
#endif /* TIMING_MPI */


#ifdef DEBUG
      write(*,*) 
     .    'pe ',mype,' starting tree setup : max_blks_sent ',
     .           max_blks_sent,' offset ',tag_offset
#endif /* DEBUG */

!------
! If we change this buffer size, remember to make the same change in
! mpi_amr_tree_setup
!------
!
! Set up buffer size info for tree messages
      len_surr_blks = 3*3*(1+2*k2d)*(1+2*k3d)
#ifdef SAVE_MORTS
     .                + 6*3*(1+2*k2d)*(1+2*k3d)
#endif
      buffer_dim = (32+16+len_surr_blks)*max_blks_sent

#ifdef DEBUG
      write(*,*) 'pe ',mype,' setup : buffer_dim ',
     .           buffer_dim
#endif /* DEBUG */
      
#ifndef ORIG
      if(allocated(recv_buf).and.allocated(send_buf).and.
     .    (size(recv_buf,dim=1).lt.buffer_dim
     .         .or.
     .     size(recv_buf,dim=1).gt.2*buffer_dim) ) then
        deallocate(send_buf)
        allocate( send_buf(buffer_dim))
        deallocate(recv_buf)
        allocate( recv_buf(buffer_dim))
      endif
      if(.not.allocated(send_buf)) allocate(send_buf(buffer_dim))
      if(.not.allocated(recv_buf)) allocate(recv_buf(buffer_dim))
      bufdim = size(send_buf,dim=1)
#endif /* ORIG */

#ifdef ORIG
      if(allocated(send_buf)) deallocate(send_buf)
      allocate( send_buf(buffer_dim))
      if(allocated(recv_buf)) deallocate(recv_buf)
      allocate( recv_buf(buffer_dim))
      bufdim = buffer_dim
#endif /* ORIG */

#ifdef DEBUG
      write(*,*) 'pe ',mype,' setup : allocated '
#endif /* DEBUG */

#ifdef TIMING_MPI
      timer_amr_comm_setup(7) =  timer_amr_comm_setup(7)
     .                          + mpi_wtime() - time2
      time2 = mpi_wtime()
#endif /* TIMING_MPI */


! Pack the tree info to be sent

#ifdef DEBUG
      write(*,*) 'pe ',mype,' tree setup : entering mpi_pack_tree_info'
#endif /* DEBUG */
      call mpi_pack_tree_info(mype,nprocs,buf_dim_bytes,
     .                     bufdim,send_buf)
!     .                     buffer_dim,send_buf)

#ifdef DEBUG
      write(*,*) 'pe ',mype,' tree setup : exited mpi_pack_tree_info'
#endif /* DEBUG */

#ifdef DEBUG
      write(*,*) 
     .'pe ',mype,' tree setup : entering mpi_xchange_blocks',
     . ' args : mype,nprocs,tag_offset,buffer_dim ',
     .  mype,nprocs,tag_offset,buffer_dim
#endif /* DEBUG */
#ifdef TIMING_MPI
      timer_amr_comm_setup(8) =  timer_amr_comm_setup(8)
     .                          + mpi_wtime() - time2
      time2 = mpi_wtime()
#endif /* TIMING_MPI */

!
! Exchange the tree info to be sent


      call mpi_xchange_blocks(mype,nprocs,tag_offset,
     .                        bufdim,send_buf,recv_buf)
!     .                        buffer_dim,send_buf,recv_buf)

#ifdef TIMING_MPI
      timer_amr_comm_setup(5) =  timer_amr_comm_setup(5)
     .                          + mpi_wtime() - time2
      time2 = mpi_wtime()
#endif /* TIMING_MPI */

#ifdef DEBUG
      write(*,*) 
     .'pe ',mype,' tree setup : exited mpi_xchange_blocks'
#endif /* DEBUG */

! Unpack the tree info which has been received
#ifdef DEBUG
      write(*,*) 
     . 'pe ',mype,' comm setup : entering mpi_unpack_treeinfo'
     .  ,' lguard ',lguard,' lprolong ',lprolong,' lflux ',
     .     lflux,' ledge ',' lrestrict ',lrestrict
#endif /* DEBUG */

      call mpi_unpack_tree_info(mype,nprocs,buf_dim_bytes,
     .                       bufdim,recv_buf)
!     .                       buffer_dim,recv_buf)

#ifdef TIMING_MPI
      timer_amr_comm_setup(9) =  timer_amr_comm_setup(9)
     .                          + mpi_wtime() - time2
      time2 = mpi_wtime()
#endif /* TIMING_MPI */

#ifdef DEBUG
      write(*,*) 
     . 'pe ',mype,' tree setup : exited mpi_unpack_treeinfo'
#endif /* DEBUG */

#ifdef ORIG
      if(allocated(send_buf)) deallocate(send_buf)
      if(allocated(recv_buf)) deallocate(recv_buf)
#endif /* ORIG */


#ifdef DEBUG
      write(*,*) 'pe ',mype,' tree setup : deallocated '
      write(*,*) 'pe ',mype,' exiting mpi_amr_comm_setup ',
     .           ' l_datapacked ',l_datapacked
#endif /* DEBUG */

      end if ! end of lrestrict if test


!---------------------------------------------------
! Store laddress pe limits to help optimize searching

      if(last_buffer.gt.strt_buffer.and.nprocs.gt.1) then

         ladd_strt(1:nprocs-1) = last_buffer
         ladd_strt(0         ) = strt_buffer
         do jj= last_buffer,strt_buffer,-1
           ladd_strt(laddress(2,jj)) = jj
         enddo
         do jj= nprocs-2,0,-1
           ladd_strt(jj) = minval(ladd_strt(jj:jj+1))
         enddo

         ladd_end(0:nprocs-2) = strt_buffer
         ladd_end(nprocs-1  ) = last_buffer
         do jj= strt_buffer,last_buffer
           ladd_end(laddress(2,jj)) = jj
         enddo
         do jj= 1,nprocs-1
           ladd_end(jj) = maxval(ladd_end(jj-1:jj))
         enddo

      else
        ladd_strt = strt_buffer
        ladd_end  = strt_buffer
      endif

!---------------------------------------------------


#ifdef TIMING_MPI
        timer_amr_comm_setup(0) =  timer_amr_comm_setup(0)
     .                          + mpi_wtime() - time1
#endif /* TIMING_MPI */

      return
      end subroutine mpi_amr_comm_setup


