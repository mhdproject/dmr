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

      subroutine amr_edge_diagonal_check(mype)




!------------------------------------------------------------------------
!
! This routine checks to see if the diagonal block between two
! leaf-neighbors at the same refinement level as the current block,
! is refined. If it is then the edge-based variables along the edge
! shared with that diagonal block is given the edge values
! form the refined diagonal block, to insure conservation properties.
!
!
! Written :     Peter MacNeice          October 1997
!------------------------------------------------------------------------


      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton

      implicit none 

#include "amr_shmem.fh"
      include 'mpif.h'

      integer, intent(in)  ::  mype

!------------------------------------
! local variables

      integer :: nguard0
      parameter (nguard0 = nguard*npgs)


      integer :: klo,kup
      integer :: jlo,jup
      integer :: ilo,iup
#ifndef LIBRARY
      parameter(klo=1+nguard0*k3d,kup=klo*k3d+nzb)
      parameter(jlo=1+nguard0*k2d,jup=jlo*k2d+nyb)
      parameter(ilo=1+nguard0,iup=ilo+nxb)
#endif


      integer :: remote_pe,remote_block
      integer :: mark_edge(12,maxblocks)
      integer :: i, ie, iblk, lb, k, j, ierrcode

      logical :: lfound

      real :: receive(nedges,maxdim+2*nguard0)

#ifdef LIBRARY
      klo=1+nguard0*k3d
      kup=klo*k3d+nzb
      jlo=1+nguard0*k2d
      jup=jlo*k2d+nyb
      ilo=1+nguard0
      iup=ilo+nxb
#endif

!------------------------------------

#ifdef DEBUG
       call amr_flush(6)
       call mpi_barrier (MPI_COMM_WORLD, ierrcode)
       write(*,*) 'starting amr_edge_diagonal_check : pe ',mype
       call amr_flush(6)
       call mpi_barrier (MPI_COMM_WORLD, ierrcode)
#endif /* DEBUG */

!#ifndef NO_PERMANENT_GUARDCELLS

!!!#if N_DIM >= 2
       if (ndim >= 2) then


! Initialize array marking edges for diagonal patching.
      mark_edge(:,:) = 0

      do i = 1,no_of_diagonal_edges
        ie   = edge_mark(6,1,i)
        iblk = edge_mark(6,2,i)
        mark_edge(ie,iblk) = i
      enddo


! Loop over the blocks on this processor.
      if(lnblocks.gt.0) then
      do lb=1,lnblocks


#ifdef VAR_DT
! Is this a leaf block which has finished its current timestep?
      if(nodetype(lb).eq.1.and.ldtcomplete(lb)) then
#else
! Is this a leaf block ?
      if(nodetype(lb).eq.1) then
#endif

#ifdef DEBUG
         write(*,*) 'amr_edge_diagonal_check : ',
     .              'checking edges on blk ',lb
         write(*,*) 'nbedges ',nbedges
#endif /* DEBUG */


! Any edges on this block which are still marked need a diagonal patch.
! Note that in the shmem_gets below, we can always assume that a
! neighbor block exists, since the edge would not have been marked
! earlier if that was not so.

! Loop over the edges on this block.
       do ie=1,nbedges


       if(mark_edge(ie,lb).ge.1) then

        lfound = .false.
        remote_block = edge_mark(6,3,mark_edge(ie,lb))
        remote_pe    = edge_mark(6,4,mark_edge(ie,lb))

#ifdef DEBUG
         write(*,*) 'amr_edge_diagonal_check : ',
     .              'data source edge for ',lb, 
     .              ' is on blk ',remote_block,remote_pe
         write(*,*) 'nbedges ',nbedges
#endif /* DEBUG */

! (remote_block,remote_pe) may be a local block or a remote block.
          if(remote_pe.ne.mype) then
            do iblk = strt_buffer,last_buffer
#ifdef DEBUG
          write(*,*) 'diag-check :pe ',mype,' blk ',lb,' edge ',ie,
     .    ' searching buffer for edge_location ',
     .            remote_block,remote_pe,' current buffer entry ',
     .          ' iblk ',iblk,' laddress ',laddress(:,iblk)
#endif /* DEBUG */
              if(remote_block.eq.laddress(1,iblk).and.
     .             remote_pe .eq.laddress(2,iblk) ) then
                remote_block = iblk
                remote_pe    = mype
                lfound = .true.
#ifdef DEBUG
        write(*,*) 'diag-check :pe ',mype,' edge_location ',
     .  edge_mark(6,3:4,mark_edge(ie,lb)),' located in buffer slot ',
     .          iblk
#endif /* DEBUG */
              endif
            enddo
          elseif(remote_pe.eq.mype) then
            lfound = .true.
          endif

! The edge data on the neighboring faces can be assumed to have been averaged
! correctly from the refined diagonal blocks.


! Now copy over the edge data from one of the neighbors.
       if(ie.eq.1) then                    ! edge: x low edge, y low edge
         do k=klo,kup-k3d
           call shmem_real_get(
     .         receive(1:nedgevar,k),
     .         bedge_facex_z(1:nedgevar,2,jup,k,remote_block),
     .         nedgevar,remote_pe)
           bedge_facex_z(:,1,1+nguard0*k2d,k,lb)= receive(:,k)
         enddo
         bedge_facey_z(:,1+nguard0,1,klo:kup-k3d,lb)=
     .                bedge_facex_z(:,1,1+nguard0*k2d,klo:kup-k3d,lb)


       elseif(ie.eq.2) then               ! edge: x low edge, y high edge
         do k=klo,kup-k3d
           call shmem_real_get(
     .         receive(1:nedgevar,k),
     .         bedge_facex_z(1:nedgevar,2,jlo,k,remote_block),
     .         nedgevar,remote_pe)
           bedge_facex_z(:,1,k2d+nguard0*k2d+nyb,k,lb)= 
     .         receive(:,k)
         enddo
         bedge_facey_z(:,1+nguard0,2,klo:kup-k3d,lb)=
     .     bedge_facex_z(:,1,k2d+nguard0*k2d+nyb,klo:kup-k3d,lb)

       elseif(ie.eq.3) then               ! edge: x high edge, y low edge
         do k=klo,kup-k3d
           call shmem_real_get(
     .         receive(1:nedgevar,k),
     .         bedge_facex_z(1:nedgevar,1,jup,k,remote_block),
     .         nedgevar,remote_pe)
           bedge_facex_z(:,2,1+nguard0*k2d,k,lb)= receive(:,k)
         enddo
         bedge_facey_z(:,1+nguard0+nxb,1,klo:kup-k3d,lb)=
     .     bedge_facex_z(:,2,1+nguard0*k2d,klo:kup-k3d,lb)


       elseif(ie.eq.4) then               ! edge: x high edge, y high edge
         do k=klo,kup-k3d
           call shmem_real_get(
     .         receive(1:nedgevar,k),
     .         bedge_facex_z(1:nedgevar,1,jlo,k,remote_block),
     .         nedgevar,remote_pe)
           bedge_facex_z(:,2,k2d+nguard0*k2d+nyb,k,lb)= 
     .         receive(:,k)
         enddo
         bedge_facey_z(:,1+nguard0+nxb,2,klo:kup-k3d,lb)=
     .     bedge_facex_z(:,2,k2d+nguard0*k2d+nyb,klo:kup-k3d,lb)


       elseif(ie.eq.5) then                ! edge: y low edge, z low edge
         do i=ilo,iup-1
           call shmem_real_get(
     .         receive(1:nedgevar,i),
     .         bedge_facey_x(1:nedgevar,i,2,kup,remote_block),
     .         nedgevar,remote_pe)
           bedge_facey_x(:,i,1,klo,lb)= receive(:,i)
         enddo
         bedge_facez_x(:,ilo:iup-1,1+nguard0*k3d,1,lb)=
     .                bedge_facey_x(:,ilo:iup-1,1,klo,lb)


       elseif(ie.eq.6) then                ! edge: y high edge, z low edge
         do i=ilo,iup-1
           call shmem_real_get(
     .         receive(1:nedgevar,i),
     .         bedge_facey_x(1:nedgevar,i,1,kup,remote_block),
     .         nedgevar,remote_pe)
           bedge_facey_x(:,i,2,klo,lb)= receive(:,i)
         enddo
         bedge_facez_x(:,ilo:iup-1,k2d+nguard0*k2d+nyb,1,lb)=
     .                bedge_facey_x(:,ilo:iup-1,2,klo,lb)


       elseif(ie.eq.7) then                ! edge: y low edge, z high edge
         do i=ilo,iup-1
           call shmem_real_get(
     .         receive(1:nedgevar,i),
     .         bedge_facey_x(1:nedgevar,i,2,klo,remote_block),
     .         nedgevar,remote_pe)
           bedge_facey_x(:,i,1,kup,lb)= receive(:,i)
         enddo
         bedge_facez_x(:,ilo:iup-1,1+nguard0*k2d,2,lb)=
     .                bedge_facey_x(:,ilo:iup-1,1,kup,lb)


       elseif(ie.eq.8) then                ! edge: y high edge, z high edge
         do i=ilo,iup-1
           call shmem_real_get(
     .         receive(1:nedgevar,i),
     .         bedge_facey_x(1:nedgevar,i,1,klo,remote_block),
     .         nedgevar,remote_pe)
           bedge_facey_x(:,i,2,kup,lb)= receive(:,i)
         enddo
         bedge_facez_x(:,ilo:iup-1,k2d+nguard0*k2d+nyb,2,lb)=
     .                bedge_facey_x(:,ilo:iup-1,2,kup,lb)

       elseif(ie.eq.9) then                ! edge: x low edge, z low edge
         do j=jlo,jup-k2d
           call shmem_real_get(
     .         receive(1:nedgevar,j),
     .         bedge_facex_y(1:nedgevar,2,j,kup,remote_block),
     .         nedgevar,remote_pe)
           bedge_facex_y(:,1,j,klo,lb)= receive(:,j)
         enddo
         bedge_facez_y(:,1+nguard0,jlo:jup-k2d,1,lb)=
     .                bedge_facex_y(:,1,jlo:jup-k2d,klo,lb)


       elseif(ie.eq.10) then                ! edge: x low edge, z high edge
         do j=jlo,jup-k2d
           call shmem_real_get(
     .         receive(1:nedgevar,j),
     .         bedge_facex_y(1:nedgevar,2,j,klo,remote_block),
     .         nedgevar,remote_pe)
           bedge_facex_y(:,1,j,kup,lb)= receive(:,j)
         enddo
         bedge_facez_y(:,1+nguard0,jlo:jup-k2d,2,lb)=
     .                bedge_facex_y(:,1,jlo:jup-k2d,kup,lb)


       elseif(ie.eq.11) then                ! edge: x high edge, z low edge
         do j=jlo,jup-k2d
           call shmem_real_get(
     .         receive(1:nedgevar,j),
     .         bedge_facex_y(1:nedgevar,1,j,kup,remote_block),
     .         nedgevar,remote_pe)
           bedge_facex_y(:,2,j,klo,lb)= receive(:,j)
         enddo
         bedge_facez_y(:,1+nguard0+nxb,jlo:jup-k2d,1,lb)=
     .                bedge_facex_y(:,2,jlo:jup-k2d,klo,lb)


       elseif(ie.eq.12) then                ! edge: x high edge, z high edge
         do j=jlo,jup-k2d
           call shmem_real_get(
     .         receive(1:nedgevar,j),
     .         bedge_facex_y(1:nedgevar,1,j,klo,remote_block),
     .         nedgevar,remote_pe)
           bedge_facex_y(:,2,j,kup,lb)= receive(:,j)
         enddo
         bedge_facez_y(:,1+nguard0+nxb,jlo:jup-k2d,2,lb)=
     .                bedge_facex_y(:,2,jlo:jup-k2d,kup,lb)
       endif


       endif


       enddo ! loop over edges

      endif

      enddo
      endif


!!!#endif       
      end if ! if (ndim >= 2)


!#endif

      return
      end subroutine amr_edge_diagonal_check
