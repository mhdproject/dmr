!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"


      subroutine amr_restrict_bnd_data(mype,flux_dir)




!------------------------------------------------------------------------
!
! This routine does the data averaging required when a child block
! passes data back to its parent. The parent receives data at the
! block boundary only.
!
! This routine provides a mechanism for passing data defined at block
! boundaries from leaf blocks back to their parents.
! The averaging rules used to combine interface values on the finer
! mesh to construct interface values on the coarser parent mesh are
! specified by the user who provides a function called amr_restrict_red
! to do this.
!
! This routine is only relevant for schemes with even number of grid points.
!
!
! Written :     Peter MacNeice          February 1997
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton

      use paramesh_interfaces, only : amr_restrict_red

#include "amr_shmem.fh"

      integer, intent(in)    :: flux_dir
      integer, intent(in)    :: mype

!------------------------------------
! local arrays

      integer remote_pe,remote_block
      integer, save :: anodetype(1)
      integer cnodetype,cneigh(2,6)
      save cnodetype,cneigh

      integer :: nprocs
      integer,external :: shmem_n_pes

      integer nguard0
      parameter(nguard0 = nguard*npgs)
      integer :: minn, maxx

      logical :: lfound



!------------------------------------

      nprocs = shmem_n_pes()

      if(mpi_pattern_id.ne.40 .and. nprocs.gt.1) then
        write(*,*) 'Paramesh error : wrong pattern being',
     .' used for pre-communication for flux restriction :',
     .' Fix - insert appropriate call to mpi_amr_comm_setup'
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif


! cycle through the grid blocks on this processor
      if(lnblocks.gt.0) then
      do lb = 1,lnblocks

! Is this a parent block of at least one leaf node?
      if(nodetype(lb).eq.2) then

! If yes then cycle through its children.
        do ich=1,nchild

          jchild = ich
          remote_pe = child(2,ich,lb)
          remote_block  = child(1,ich,lb)

#ifdef MPI_USED
! if (remote_block,remote_pe) is not a local block then it must have a
! local copy available in the buffer space at the end of the local
! block list.
          if(remote_pe.ne.mype) then
             lfound = .false.
             iblk = ladd_strt(remote_pe)
             do while(.not.lfound.and.iblk.le.ladd_end(remote_pe))
!            do iblk = strt_buffer,last_buffer
#ifdef DEBUG
             write(*,*) 'pe ',mype,' searching buffer for ',
     .            remote_block,remote_pe,' current buffer entry ',
     .          ' iblk ',iblk,' laddress ',laddress(:,iblk)
#endif /* DEBUG */
              if(remote_block.eq.laddress(1,iblk).and.
     .             remote_pe .eq.laddress(2,iblk) ) then
                remote_block = iblk
                remote_pe    = mype
                lfound = .true.
#ifdef DEBUG
             write(*,*) 'pe ',mype,' remote block ',
     .          remote_block,remote_pe,' located in buffer slot ',
     .          iblk
#endif /* DEBUG */
              else
                iblk = iblk+1
              endif
            enddo
          endif
#endif /* MPI_USED */


! Is this child a leaf block(nodetype=1)?
! If it is then fetch its data.
          call shmem_integer_get(anodetype,nodetype(remote_block),
     .       1,remote_pe)
          cnodetype = anodetype(1)
          if(cnodetype.eq.1) then

! fetch child's neighbors. This info will be needed when we are
! advancing the solution on all levels.
            call shmem_integer_get(cneigh(1,1),neigh(1,1,remote_block),
     .       12,remote_pe)

            if (flux_dir == 1 .or.
     .          flux_dir == 2 .or.
     .          flux_dir == 3) then
               minn = flux_dir
               maxx = flux_dir
            else
               minn = 1
               maxx = ndim
            end if

            do icoord=minn,maxx
              if(icoord.eq.1) then
                call shmem_real_get( recvarxf(1,1,1,1),
     .                               flux_x(1,1,1,1,remote_block),
     .                               len_block_bndx*nfluxes,remote_pe)
              elseif(icoord.eq.2) then
                call shmem_real_get( recvaryf(1,1,1,1),
     .                               flux_y(1,1,1,1,remote_block),
     .                               len_block_bndy*nfluxes,remote_pe)
              elseif(icoord.eq.3) then
                call shmem_real_get( recvarzf(1,1,1,1),
     .                               flux_z(1,1,1,1,remote_block),
     .                               len_block_bndz*nfluxes,remote_pe)
              endif


! compute the offset in the parent block appropriate for this child
       iface = mod(jchild-1,2)
       jface = mod((jchild-1)/2,2)
       kface = mod((jchild-1)/4,2)
       ioff = iface*nxb/2
       joff = jface*nyb/2
       koff = kface*nzb/2

! Compute restricted data from the data in the buffer and
! update only boundary values on the parent block
       if(icoord.eq.1) then

         i = iface+1
! apply, only if the appropriate child neighbor does not exist
         if(cneigh(1,i).gt.-20.and.cneigh(1,i).lt.0) then


         call amr_restrict_red(icoord)
         do k=1+nguard0*k3d,nzb+nguard0*k3d,2
           kk = (k-nguard0*k3d)/2+nguard0*k3d+1
           do j=1+nguard0*k2d,nyb+nguard0*k2d,2
             jj = (j-nguard0*k2d)/2+nguard0*k2d+1
             do ivar=1,nfluxes
               flux_x(ivar,i,jj+joff,kk+koff,lb) = 
     .                                     bndtempx1(ivar,i,j,k)
             enddo
           enddo
         enddo
         endif


       elseif(icoord.eq.2) then

         j = jface+1
! apply, only if the appropriate child neighbor does not exist
         if(cneigh(1,j+2).gt.-20.and.cneigh(1,j+2).lt.0) then

         call amr_restrict_red(icoord)
         do k=1+nguard0*k3d,nzb+nguard0*k3d,2
           kk = (k-nguard0*k3d)/2+nguard0*k3d+1
           do i=1+nguard0,nxb+nguard0,2
             ii = (i-nguard0)/2+nguard0+1
             do ivar=1,nfluxes
               flux_y(ivar,ii+ioff,j,kk+koff,lb)=bndtempy1(ivar,i,j,k)
             enddo
           enddo
         enddo
         endif

       elseif(icoord.eq.3) then

         k = kface+1
! apply, only if the appropriate child neighbor does not exist
         if(cneigh(1,k+4).gt.-20.and.cneigh(1,k+4).lt.0) then

         call amr_restrict_red(icoord)
         do j=1+nguard0*k2d,nyb+nguard0*k2d,2
           jj = (j-nguard0*k2d)/2+nguard0*k2d+1
           do i=1+nguard0,nxb+nguard0,2
             ii = (i-nguard0)/2+nguard0+1
             do ivar=1,nfluxes
               flux_z(ivar,ii+ioff,jj+joff,k,lb)=bndtempz1(ivar,i,j,k)
             enddo
           enddo
         enddo

         endif

       endif

       enddo

       endif

       enddo

      endif

      enddo
      endif

      return
      end subroutine amr_restrict_bnd_data
