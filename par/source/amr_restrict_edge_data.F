!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

      subroutine amr_restrict_edge_data(mype)




!------------------------------------------------------------------------
!
! This routine does the data averaging on cell edges required when a 
! child block passes data back to its parent. The parent receives data 
! at the block boundary only.
!
! This routine provides a mechanism for passing data defined at block
! boundaries from leaf blocks back to their parents.
! The averaging rules used to combine interface values on the finer
! mesh to construct interface values on the coarser parent mesh are
! specified by the user who provides a function called amr_restrict_edge
! to do this.
!
! This routine is only relevant for schemes with even number of grid points.
!
!
! Written :     Peter MacNeice          July 1997
!------------------------------------------------------------------------


      use paramesh_dimensions
      use physicaldata
      use tree

      use paramesh_interfaces, only : amr_restrict_edge


#include "amr_shmem.fh"

      integer, intent(in)    :: mype

!------------------------------------
! local arrays

      integer remote_pe,remote_block
      integer cnodetype,cneigh(2,6)
      save cnodetype,cneigh

      integer nguard0
      parameter(nguard0 = nguard*npgs)

      integer :: lb,ich,jchild,icoord,iface,jface,kface
      integer :: ioff,joff,koff,i,j,k,ii,jj,kk

!------------------------------------

      call shmem_barrier_all()

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

! Is this child a leaf block(nodetype=1)? 
! If it is then fetch its data.
       call shmem_integer_get(cnodetype,nodetype(remote_block),
     .       1,remote_pe)
       if(cnodetype.eq.1) then

! fetch child's neighbors. This info will be needed when we are
! advancing the solution on all levels.
            call shmem_integer_get(cneigh,neigh(1,1,remote_block),
     .       12,remote_pe)


         do icoord=1,ndim
           if(icoord.eq.1) then
             call shmem_real_get(recvarx1e(1,1,1,1),
     .           bedge_facex_y(1,1,1,1,remote_block),
     .           len_block_ex*nedges,remote_pe)
             if((ndim.eq.3).or.(l2p5d.eq.1))
     .         call shmem_real_get(recvarx2e(1,1,1,1),
     .             bedge_facex_z(1,1,1,1,remote_block),
     .             len_block_ex*nedges,remote_pe)
           elseif(icoord.eq.2) then
             if((ndim.eq.3).or.(l2p5d.eq.1))
     .          call shmem_real_get(recvary1e(1,1,1,1),
     .              bedge_facey_z(1,1,1,1,remote_block),
     .              len_block_ey*nedges,remote_pe)
             call shmem_real_get(recvary2e(1,1,1,1),
     .           bedge_facey_x(1,1,1,1,remote_block),
     .           len_block_ey*nedges,remote_pe)
           elseif(icoord.eq.3) then
             call shmem_real_get(recvarz1e(1,1,1,1),
     .           bedge_facez_x(1,1,1,1,remote_block),
     .           len_block_ez*nedges,remote_pe)
             call shmem_real_get(recvarz2e(1,1,1,1),
     .           bedge_facez_y(1,1,1,1,remote_block),
     .           len_block_ez*nedges,remote_pe)
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

             call amr_restrict_edge(icoord)
             do k=1+nguard0*k3d,nzb+(nguard0+1)*k3d,2
               kk = (k-nguard0*k3d)/2+nguard0*k3d+1
               do j=1+nguard0,nyb+nguard0,2
                 jj = (j-nguard0)/2+nguard0+1
                 bedge_facex_y(:,i,jj+joff,kk+koff,
     .                                    lb) = recvarx1e(:,i,j,k)
               enddo
             enddo

             if((ndim.eq.3).or.(l2p5d.eq.1)) then
               do k=1+nguard0*k3d,nzb+nguard0*k3d,2
                 kk = (k-nguard0*k3d)/2+nguard0*k3d+1
                 do j=1+nguard0,nyb+nguard0+1,2
                   jj = (j-nguard0)/2+nguard0+1
                   bedge_facex_z(:,i,jj+joff,kk+koff,
     .                                    lb) = recvarx2e(:,i,j,k)
                 enddo
               enddo
             endif

             endif


           elseif(icoord.eq.2) then

             j = jface+1
! apply, only if the appropriate child neighbor does not exist
             if(cneigh(1,j+2).gt.-20.and.cneigh(1,j+2).lt.0) then

             call amr_restrict_edge(icoord)
             do k=1+nguard0*k3d,nzb+(nguard0+1)*k3d,2
               kk = (k-nguard0*k3d)/2+nguard0*k3d+1
               do i=1+nguard0,nxb+nguard0,2
                 ii = (i-nguard0)/2+nguard0+1
                 bedge_facey_x(:,ii+ioff,j,kk+koff,
     .                                    lb) = recvary2e(:,i,j,k)

               enddo
             enddo

             if((ndim.eq.3).or.(l2p5d.eq.1)) then
               do k=1+nguard0*k3d,nzb+nguard0*k3d,2
                 kk = (k-nguard0*k3d)/2+nguard0*k3d+1
                 do i=1+nguard0,nxb+nguard0+1,2
                   ii = (i-nguard0)/2+nguard0+1
                   bedge_facey_z(:,ii+ioff,j,kk+koff,
     .                                    lb) = recvary1e(:,i,j,k)
                 enddo
               enddo
             endif

             endif

           elseif(icoord.eq.3) then

             k = kface+1
! apply, only if the appropriate child neighbor does not exist
             if(cneigh(1,k+4).gt.-20.and.cneigh(1,k+4).lt.0) then

             call amr_restrict_edge(icoord)
             do j=1+nguard0,nyb+nguard0+1,2
               jj = (j-nguard0)/2+nguard0+1
               do i=1+nguard0,nxb+nguard0,2
                 ii = (i-nguard0)/2+nguard0+1
                 bedge_facez_x(:,ii+ioff,jj+joff,k,
     .                                    lb) = recvarz1e(:,i,j,k)
               enddo
             enddo
             do j=1+nguard0,nyb+nguard0,2
               jj = (j-nguard0)/2+nguard0+1
               do i=1+nguard0,nxb+nguard0+1,2
                 ii = (i-nguard0)/2+nguard0+1
                 bedge_facez_y(:,ii+ioff,jj+joff,k,
     .                                    lb) = recvarz2e(:,i,j,k)
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

      call shmem_barrier_all()

      return
      end subroutine amr_restrict_edge_data
