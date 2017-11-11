!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"


#undef DEBUG

      subroutine amr_gsurrounding_blks(mype,ldiag)




!------------------------------------------------------------------------
!
! This routine constructs a list of the addresses of the blocks which
! surround all the local blocks.
!
! Written :     Peter MacNeice          July 2000
!------------------------------------------------------------------------
!
! Arguments:
!       mype            local processor number
!       ldiag           a logical switch which controls whether blocks
!                        diagonally opposite edges and corners are included
!                        in the list of surrounding blocks.
!                                   if .true.  include diagonals 
!                                   if .false. do not include diagonals 
!
!------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace


#include "amr_shmem.fh"

      integer, intent(in)    ::  mype
      logical, intent(in)    ::  ldiag


!------------------------------------
! local arrays

        integer :: surrblks(3,3,3,3)
        integer :: psurrblks(3,3,1+2*k2d,1+2*k3d)
        integer :: coarse_surrblks(2,6,6,6)

        integer :: cneigh(2,6),dneigh(2,6),cchild(2,8)
        save    cneigh,dneigh,psurrblks,cchild

        integer :: rem_blk,rem_pe,ctype
        save ctype

        logical :: lmask(3,3,3),lupdate_diag

        integer :: pe_r,lb_r,remote_addr
        integer :: parent_blk,parent_pe,remote_blk,remote_pe

        integer :: lb,iface

!------------------------------------

      do lb = 1,lnblocks
      surr_blks(:,:,:,:,lb) = -1


        dneigh(:,:) = neigh(:,:,lb)


        if(ndim.eq.1) then
          surrblks(1:2,1,2,2) = dneigh(1:2,1)
          surrblks(1:2,3,2,2) = dneigh(1:2,2)
          surrblks(1,2,2,2) = lb
          surrblks(2,2,2,2) = mype
          surrblks(3,2,2,2) = nodetype(lb)

          surr_blks(:,:,1,1,lb) = surrblks(:,:,2,2)

        elseif(ndim.ge.1) then
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section 1.
!
! Get addresses of neighbor blocks
        surrblks(:,:,:,:) = -1

        surrblks(1:2,1,2,2) = dneigh(1:2,1)
        surrblks(1:2,3,2,2) = dneigh(1:2,2)
        if(ndim.ge.2) then
          surrblks(1:2,2,1,2) = dneigh(1:2,3)
          surrblks(1:2,2,3,2) = dneigh(1:2,4)
          if(ndim.eq.3) then
            surrblks(1:2,2,2,1) = dneigh(1:2,5)
            surrblks(1:2,2,2,3) = dneigh(1:2,6)
          endif
        endif
        surrblks(1,2,2,2) = lb
        surrblks(2,2,2,2) = mype

!
! Get neigh's neigh array
        do iface=1,nfaces
          lb_r = dneigh(1,iface)
          pe_r = dneigh(2,iface)

          if(lb_r.gt.0) then
            call shmem_integer_get(cneigh,neigh(1,1,lb_r),12,pe_r)
          elseif(lb_r.le.-20) then
            cneigh(1:2,1:6) = lb_r
          endif
          if(lb_r.gt.0.or.lb_r.le.-20) then
            if(iface.eq.1) then
              surrblks(1:2,1,1,2) = cneigh(1:2,3)
              surrblks(1:2,1,3,2) = cneigh(1:2,4)
              surrblks(1:2,1,2,1) = cneigh(1:2,5)
              surrblks(1:2,1,2,3) = cneigh(1:2,6)
            elseif(iface.eq.2) then
              surrblks(1:2,3,1,2) = cneigh(1:2,3)
              surrblks(1:2,3,3,2) = cneigh(1:2,4)
              surrblks(1:2,3,2,1) = cneigh(1:2,5)
              surrblks(1:2,3,2,3) = cneigh(1:2,6)
            elseif(iface.eq.3) then
              surrblks(1:2,1,1,2) = cneigh(1:2,1)
              surrblks(1:2,3,1,2) = cneigh(1:2,2)
              surrblks(1:2,2,1,1) = cneigh(1:2,5)
              surrblks(1:2,2,1,3) = cneigh(1:2,6)
            elseif(iface.eq.4) then
              surrblks(1:2,1,3,2) = cneigh(1:2,1)
              surrblks(1:2,3,3,2) = cneigh(1:2,2)
              surrblks(1:2,2,3,1) = cneigh(1:2,5)
              surrblks(1:2,2,3,3) = cneigh(1:2,6)
            elseif(iface.eq.5) then
              surrblks(1:2,1,2,1) = cneigh(1:2,1)
              surrblks(1:2,3,2,1) = cneigh(1:2,2)
              surrblks(1:2,2,1,1) = cneigh(1:2,3)
              surrblks(1:2,2,3,1) = cneigh(1:2,4)
            elseif(iface.eq.6) then
              surrblks(1:2,1,2,3) = cneigh(1:2,1)
              surrblks(1:2,3,2,3) = cneigh(1:2,2)
              surrblks(1:2,2,1,3) = cneigh(1:2,3)
              surrblks(1:2,2,3,3) = cneigh(1:2,4)
            endif                                  ! end of iface if
          endif                                    ! end of lb_r if
        enddo                                      ! end of loop over iface


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Now pick up addresses of corner blocks
       if(ndim.eq.3.and.ldiag) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Section 2.

! First pick up corners from left x face

! Lower corners( ie low z)
        lb_r = surrblks(1,1,2,1)
        pe_r = surrblks(2,1,2,1)
        if(lb_r.gt.0) then
            call shmem_integer_get(cneigh,neigh(1,1,lb_r),12,pe_r)
        elseif(lb_r.le.-20) then
            cneigh(1:2,3:4) = lb_r
        endif
        if(lb_r.gt.0.or.lb_r.le.-20) then
            surrblks(1:2,1,1,1) = cneigh(1:2,3)
            surrblks(1:2,1,3,1) = cneigh(1:2,4)
        endif

! Upper corners( ie high z)
        lb_r = surrblks(1,1,2,3)
        pe_r = surrblks(2,1,2,3)
        if(lb_r.gt.0) then
            call shmem_integer_get(cneigh,neigh(1,1,lb_r),12,pe_r)
        elseif(lb_r.le.-20) then
            cneigh(1:2,3:4) = lb_r
        endif
        if(lb_r.gt.0.or.lb_r.le.-20) then
            surrblks(1:2,1,1,3) = cneigh(1:2,3)
            surrblks(1:2,1,3,3) = cneigh(1:2,4)
        endif

! Finally find corners from right x face
! Lower corners( ie low z)
        lb_r = surrblks(1,3,2,1)
        pe_r = surrblks(2,3,2,1)
        if(lb_r.gt.0) then
            call shmem_integer_get(cneigh,neigh(1,1,lb_r),12,pe_r)
        elseif(lb_r.le.-20) then
            cneigh(1:2,3:4) = lb_r
        endif
        if(lb_r.gt.0.or.lb_r.le.-20) then
            surrblks(1:2,3,1,1) = cneigh(1:2,3)
            surrblks(1:2,3,3,1) = cneigh(1:2,4)
        endif

! Upper corners( ie high z)
        lb_r = surrblks(1,3,2,3)
        pe_r = surrblks(2,3,2,3)
        if(lb_r.gt.0) then
            call shmem_integer_get(cneigh,neigh(1,1,lb_r),12,pe_r)
        elseif(lb_r.le.-20) then
            cneigh(1:2,3:4) = lb_r
        endif
        if(lb_r.gt.0.or.lb_r.le.-20) then
            surrblks(1:2,3,1,3) = cneigh(1:2,3)
            surrblks(1:2,3,3,3) = cneigh(1:2,4)
        endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! If an external boundary is known to interfere with this corner search,
! repeat section 2 but for y and z faces.

! Fix corner 1 by trying 2nd path to it.
        if(surrblks(1,1,1,1).eq.-1) then
        lb_r = surrblks(1,2,1,1)
        pe_r = surrblks(2,2,1,1)
        if(lb_r.gt.0) then
          call shmem_integer_get(cneigh,neigh(1,1,lb_r),12,pe_r)
        elseif(lb_r.le.-20) then
            cneigh(1:2,1) = lb_r
        endif
        if(lb_r.gt.0.or.lb_r.le.-20) then
          surrblks(1:2,1,1,1) = cneigh(1:2,1)
        else
! If 2nd path didn't work try 3rd path.
          lb_r = surrblks(1,1,1,2)
          pe_r = surrblks(2,1,1,2)
          if(lb_r.gt.0) then
            call shmem_integer_get(cneigh,neigh(1,1,lb_r),12,pe_r)
          elseif(lb_r.le.-20) then
            cneigh(1:2,5) = lb_r
          endif
          if(lb_r.gt.0.or.lb_r.le.-20) then
            surrblks(1:2,1,1,1) = cneigh(1:2,5)
          endif
        endif
        endif


! Fix corner 2 by trying 2nd path to it.
        if(surrblks(1,1,3,1).eq.-1) then
        lb_r = surrblks(1,2,3,1)
        pe_r = surrblks(2,2,3,1)
        if(lb_r.gt.0) then
          call shmem_integer_get(cneigh,neigh(1,1,lb_r),12,pe_r)
        elseif(lb_r.le.-20) then
            cneigh(1:2,1) = lb_r
        endif
        if(lb_r.gt.0.or.lb_r.le.-20) then
          surrblks(1:2,1,3,1) = cneigh(1:2,1)
        else
! If 2nd path didn't work try 3rd path.
          lb_r = surrblks(1,1,3,2)
          pe_r = surrblks(2,1,3,2)
          if(lb_r.gt.0) then
            call shmem_integer_get(cneigh,neigh(1,1,lb_r),12,pe_r)
          elseif(lb_r.le.-20) then
            cneigh(1:2,5) = lb_r
          endif
          if(lb_r.gt.0.or.lb_r.le.-20) then
            surrblks(1:2,1,3,1) = cneigh(1:2,5)
          endif
        endif
        endif


! Fix corner 3 by trying 2nd path to it.
        if(surrblks(1,1,1,3).eq.-1) then
        lb_r = surrblks(1,2,1,3)
        pe_r = surrblks(2,2,1,3)
        if(lb_r.gt.0) then
          call shmem_integer_get(cneigh,neigh(1,1,lb_r),12,pe_r)
        elseif(lb_r.le.-20) then
            cneigh(1:2,1) = lb_r
        endif
        if(lb_r.gt.0.or.lb_r.le.-20) then
          surrblks(1:2,1,1,3) = cneigh(1:2,1)
        else
! If 2nd path didn't work try 3rd path.
          lb_r = surrblks(1,1,1,2)
          pe_r = surrblks(2,1,1,2)
          if(lb_r.gt.0) then
            call shmem_integer_get(cneigh,neigh(1,1,lb_r),12,pe_r)
          elseif(lb_r.le.-20) then
            cneigh(1:2,6) = lb_r
          endif
          if(lb_r.gt.0.or.lb_r.le.-20) then
            surrblks(1:2,1,1,3) = cneigh(1:2,6)
          endif
        endif
        endif


! Fix corner 4 by trying 2nd path to it.
        if(surrblks(1,1,3,3).eq.-1) then
        lb_r = surrblks(1,2,3,3)
        pe_r = surrblks(2,2,3,3)
        if(lb_r.gt.0) then
          call shmem_integer_get(cneigh,neigh(1,1,lb_r),12,pe_r)
        elseif(lb_r.le.-20) then
            cneigh(1:2,1) = lb_r
        endif
        if(lb_r.gt.0.or.lb_r.le.-20) then
          surrblks(1:2,1,3,3) = cneigh(1:2,1)
        else
! If 2nd path didn't work try 3rd path.
          lb_r = surrblks(1,1,3,2)
          pe_r = surrblks(2,1,3,2)
          if(lb_r.gt.0) then
            call shmem_integer_get(cneigh,neigh(1,1,lb_r),12,pe_r)
          elseif(lb_r.le.-20) then
            cneigh(1:2,6) = lb_r
          endif
          if(lb_r.gt.0.or.lb_r.le.-20) then
            surrblks(1:2,1,3,3) = cneigh(1:2,6)
          endif
        endif
        endif

! Fix corner 5 by trying 2nd path to it.
        if(surrblks(1,3,1,1).eq.-1) then
        lb_r = surrblks(1,2,1,1)
        pe_r = surrblks(2,2,1,1)
        if(lb_r.gt.0) then
          call shmem_integer_get(cneigh,neigh(1,1,lb_r),12,pe_r)
        elseif(lb_r.le.-20) then
            cneigh(1:2,2) = lb_r
        endif
        if(lb_r.gt.0.or.lb_r.le.-20) then
          surrblks(1:2,3,1,1) = cneigh(1:2,2)
        else
! If 2nd path didn't work try 3rd path.
          lb_r = surrblks(1,3,1,2)
          pe_r = surrblks(2,3,1,2)
          if(lb_r.gt.0) then
            call shmem_integer_get(cneigh,neigh(1,1,lb_r),12,pe_r)
          elseif(lb_r.le.-20) then
            cneigh(1:2,5) = lb_r
          endif
          if(lb_r.gt.0.or.lb_r.le.-20) then
            surrblks(1:2,3,1,1) = cneigh(1:2,5)
          endif
        endif
        endif


! Fix corner 6 by trying 2nd path to it.
        if(surrblks(1,3,3,1).eq.-1) then
        lb_r = surrblks(1,2,3,1)
        pe_r = surrblks(2,2,3,1)
        if(lb_r.gt.0) then
          call shmem_integer_get(cneigh,neigh(1,1,lb_r),12,pe_r)
        elseif(lb_r.le.-20) then
            cneigh(1:2,2) = lb_r
        endif
        if(lb_r.gt.0.or.lb_r.le.-20) then
          surrblks(1:2,3,3,1) = cneigh(1:2,2)
        else
! If 2nd path didn't work try 3rd path.
          lb_r = surrblks(1,3,3,2)
          pe_r = surrblks(2,3,3,2)
          if(lb_r.gt.0) then
            call shmem_integer_get(cneigh,neigh(1,1,lb_r),12,pe_r)
          elseif(lb_r.le.-20) then
            cneigh(1:2,5) = lb_r
          endif
          if(lb_r.gt.0.or.lb_r.le.-20) then
            surrblks(1:2,3,3,1) = cneigh(1:2,5)
          endif
        endif
        endif


! Fix corner 7 by trying 2nd path to it.
        if(surrblks(1,3,1,3).eq.-1) then
        lb_r = surrblks(1,2,1,3)
        pe_r = surrblks(2,2,1,3)
        if(lb_r.gt.0) then
          call shmem_integer_get(cneigh,neigh(1,1,lb_r),12,pe_r)
        elseif(lb_r.le.-20) then
            cneigh(1:2,2) = lb_r
        endif
        if(lb_r.gt.0.or.lb_r.le.-20) then
          surrblks(1:2,3,1,3) = cneigh(1:2,2)
        else
! If 2nd path didn't work try 3rd path.
          lb_r = surrblks(1,3,1,2)
          pe_r = surrblks(2,3,1,2)
          if(lb_r.gt.0) then
            call shmem_integer_get(cneigh,neigh(1,1,lb_r),12,pe_r)
          elseif(lb_r.le.-20) then
            cneigh(1:2,6) = lb_r
          endif
          if(lb_r.gt.0.or.lb_r.le.-20) then
            surrblks(1:2,3,1,3) = cneigh(1:2,6)
          endif
        endif
        endif


! Fix corner 8 by trying 2nd path to it.
        if(surrblks(1,3,3,3).eq.-1) then
        lb_r = surrblks(1,2,3,3)
        pe_r = surrblks(2,2,3,3)
        if(lb_r.gt.0) then
          call shmem_integer_get(cneigh,neigh(1,1,lb_r),12,pe_r)
        elseif(lb_r.le.-20) then
            cneigh(1:2,2) = lb_r
        endif
        if(lb_r.gt.0.or.lb_r.le.-20) then
          surrblks(1:2,3,3,3) = cneigh(1:2,2)
        else
! If 2nd path didn't work try 3rd path.
          lb_r = surrblks(1,3,3,2)
          pe_r = surrblks(2,3,3,2)
          if(lb_r.gt.0) then
            call shmem_integer_get(cneigh,neigh(1,1,lb_r),12,pe_r)
          elseif(lb_r.le.-20) then
            cneigh(1:2,6) = lb_r
          endif
          if(lb_r.gt.0.or.lb_r.le.-20) then
            surrblks(1:2,3,3,3) = cneigh(1:2,6)
          endif
        endif
        endif


        endif                                 ! end of ndim=3 if

        surr_blks(:,:,1:1+2*k2d,1:1+2*k3d,lb) =
     .         surrblks(:,:,2-k2d:2+k2d,2-k3d:2+k3d)
!          if(ndim.eq.2) then
!            surr_blks(:,:,:,1,lb) = surrblks(:,:,:,1)
!          elseif(ndim.eq.3) then
!            surr_blks(:,:,:,:,lb) = surrblks(:,:,:,:)
!          endif

        endif                                 ! end of ndim if


!---
! get nodetypes for neighbor blocks
        do k = 1,1+2*k3d
        do j = 1,1+2*k2d
        do i = 1,3
          if(surr_blks(1,i,j,k,lb).gt.0) then
          if(i.eq.2.and.j.eq.2.and.k.eq.2) then
            surr_blks(3,i,j,k,lb) = nodetype(lb)
          else 
            rem_blk = surr_blks(1,i,j,k,lb)
            rem_pe  = surr_blks(2,i,j,k,lb)
            call shmem_integer_get(ctype,nodetype(rem_blk),1,rem_pe)
            surr_blks(3,i,j,k,lb) = ctype
          endif
          else
            surr_blks(3,i,j,k,lb) = -1
          endif
        enddo
        enddo
        enddo
!---


      enddo                            ! end of loop over blocks


      call shmem_barrier_all()     ! make sure other procs are done before
                                   ! trusting their parent surr_blks which 
                                   ! may be needed in the next section

!
! At this point the only remaining neighbors which have not been
! located are edge and corner blocks which are not reachable by
! testing off the remote face. To locate any of these, we need to
! use a brute force approach.
      do lb = 1,lnblocks
      if(nodetype(lb).eq.1.and.parent(1,lb).gt.0) then

! Do we need to check this block for missing diagonal addresses?
      lupdate_diag = .false.

      lmask=.false.
      lmask(1,1,1) = .true.
      lmask(3,1,1) = .true.
      if(ndim.eq.2) then
        lmask(1,3,1) = .true.
        lmask(3,3,1) = .true.
      endif
      if(ndim.eq.3) then
        lmask(1,1,:) = .true.
        lmask(3,1,:) = .true.
        lmask(1,3,:) = .true.
        lmask(3,3,:) = .true.
        lmask(:,1,1) = .true.
        lmask(:,3,1) = .true.
        lmask(:,1,3) = .true.
        lmask(:,3,3) = .true.
        lmask(1,:,1) = .true.
        lmask(3,:,1) = .true.
        lmask(1,:,3) = .true.
        lmask(3,:,3) = .true.
      endif
 
      do k=1,1+2*k3d
      do j=1,1+2*k2d
      do i=1,3
        remote_addr = surr_blks(1,i,j,k,lb)
        if( (lmask(i,j,k)) .and.
     .      remote_addr.lt.0 ) lupdate_diag=.true.
      enddo 
      enddo 
      enddo 

      if(lupdate_diag) then

      parent_blk = parent(1,lb)
      parent_pe  = parent(2,lb)

! Get the list of blocks surrounding the parent of this block.
! This will be sufficiently complete for our current purpose,
! since, if the parent has children, it will already know any
! of its neighbors which are required to help construct the childs
! list of neighbors.
      call shmem_integer_get(psurrblks,surr_blks(1,1,1,1,parent_blk),
     .                       3*3*(1+2*k2d)*(1+2*k3d),parent_pe)


! Loop over all the neighbors of the parent.
      do k=1,1+2*k3d
      do j=1,1+2*k2d
      do i=1,3

! If this neighbor exists then get the list of its children
        remote_blk = psurrblks(1,i,j,k)
        remote_pe  = psurrblks(2,i,j,k)
        if(remote_blk.gt.0) then           ! block exists
          call shmem_integer_get(cchild,
     .                          child(1,1,remote_blk),16,remote_pe)
        elseif(remote_blk.le.-20) then     ! boundary condition to be recorded
          do ii = 1,2
          do jj = 1,8
            cchild(ii,jj) = remote_blk
          enddo
          enddo
        endif

! Now read list of surrounding children into enlarged array.
        ll = 0
        do kk = 2*k-1,2*k-1+k3d
        do jj = 2*j-1,2*j-1+k2d
        do ii = 2*i-1,2*i
          ll = ll+1
          coarse_surrblks(1:2,ii,jj,kk) = cchild(1:2,ll)
        enddo
        enddo
        enddo

      enddo
      enddo
      enddo

! coarse_surrblks now stores any edge and corner addresses needed
! by block lb. To find their addresses we need to use which_child 
! to locate lb in its parent, and then extract the 3x3x3 array
! of addresses centered on that location.
      ich = which_child(lb)
      ioff = mod(which_child(lb)-1,2)
      joff = mod((which_child(lb)-1)/2,2)
      koff = mod((which_child(lb)-1)/4,2) 
      surr_blks(1:2,1:3,1:1+2*k2d,1:1+2*k3d,lb) = 
     .         coarse_surrblks(1:2,2+ioff:4+ioff,
     .                           1+(1+joff)*k2d:1+(3+joff)*k2d,
     .                           1+(1+koff)*k3d:1+(3+koff)*k3d )


!---
! finally get nodetypes for these neighbor blocks
        do k = 1,1+2*k3d
        do j = 1,1+2*k2d
        do i = 1,3
          if(surr_blks(1,i,j,k,lb).gt.0) then
          if(i.eq.2.and.j.eq.2.and.k.eq.2) then
            surr_blks(3,i,j,k,lb) = nodetype(lb)
          else
            rem_blk = surr_blks(1,i,j,k,lb)
            rem_pe  = surr_blks(2,i,j,k,lb)
            call shmem_integer_get(ctype,nodetype(rem_blk),1,rem_pe)
            surr_blks(3,i,j,k,lb) = ctype
          endif
          else
            surr_blks(3,i,j,k,lb) = -1
          endif
        enddo
        enddo
        enddo
!---


      endif              ! end of lupdate_diag if test
      endif              ! end of nodetype if test

#ifdef DEBUG
      write(*,*) 'Neighboring blocks of block ',lb,' on pe ',mype
      do k=1,1+2*k3d
      do j=1,1+2*k2d
      write(*,10) j,k, (surrblks(1,i,j,k),surrblks(2,i,j,k),i=1,3)
10    format(i3,2x,i3,3(2x,'(',i3,'/',i3,')'))
      enddo
      enddo
#endif

      enddo                            ! end of loop over blocks

! Record that routine has been called.
      gsurrblks_set = +1


      return
      end subroutine amr_gsurrounding_blks
