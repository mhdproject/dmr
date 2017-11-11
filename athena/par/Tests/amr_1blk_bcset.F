!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"

      subroutine amr_1blk_bcset(mype,ibc,lb,pe,
     .                          idest,iopt,ibnd,jbnd,kbnd,
     .                          surrblks)






!------------------------------------------------------------------------
!
! This routine sets guardcell values at external boundaries in the case
! where a single block is having its guardcells filled.
!
! It can be assumed in writing this routine, that all guardcells for this
! block which are not across an external boundary have already been
! properly filled.
!
!
! Arguments:
!      mype             local processor
!      ibc              the integer specifying the particular boundary
!                        condition to be imposed
!      lb               block number of selected block
!      pe               processor on which block lb is located
!      idest            selects the storage space in data_1blk.fh which is to
!                        be used in this call. If the leaf node is having its
!                        guardcells filled then set this to 1, if its parent
!                        is being filled set it to 2.
!      ibnd             a selector setting designating whether the guarcells
!                        to be set are at the left, center or right sections
!                        of the i index range, eg
!                           ibnd = -1      left end
!                                =  0      middle
!                                = +1      right. For example, if ibnd=-1,
!                        the i index applied when filling unk will run
!                        from 1:nguard, if ibnd=0 from 1+nguard:nxb+nguard,
!                        and if ibnd=+1 from nxb+nguard+1:nxb+2*nguard.
!      jbnd             a selector setting designating whether the guarcells
!                        to be set are at the left, center or right sections
!                        of the j index range.
!      kbnd             a selector setting designating whether the guarcells
!                        to be set are at the left, center or right sections
!                        of the k index range.
!
!
! Written :     Peter MacNeice          August 1998
! Modified:     Peter MacNeice          January 2001
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace

      implicit real(a-h,o-z)

#include "amr_shmem.fh"

! Only necessary for programs in ./Tests
#include "test_defs.fh"

      integer, intent(in) :: mype,ibc,lb,pe
      integer, intent(in) :: idest,iopt,ibnd,jbnd,kbnd      
      integer, intent(in) :: surrblks(:,:,:,:)

      real ccoord(3),cbsize(3)
      save    ccoord,csize

#ifdef MPI_USED
      integer :: rem_block, rem_pe, ierrorcode,ierr
      logical :: lfound
#endif /* MPI_USED */

!-------------------------
! Test specific



#if defined(SHMEM_USED) | defined(SHMEM_DUMMYLIB)
      call shmem_real_get(cbsize(1),bsize(1,lb),3,pe)
      call shmem_real_get(ccoord(1),coord(1,lb),3,pe)
#endif /* SHMEM_USED */
#ifdef MPI_USED
      if(pe.eq.mype) then
        call shmem_real_get(cbsize(1),bsize(1,lb),3,pe)
        call shmem_real_get(ccoord(1),coord(1,lb),3,pe)
      else
        lfound = .false.
        iblk = strt_buffer
        do while(.not.lfound.and.iblk.le.last_buffer)
          if(lb.eq.laddress(1,iblk).and.pe.eq.laddress(2,iblk) ) then
              rem_block = iblk
              rem_pe    = mype
              lfound = .true.
          endif
          iblk = iblk+1
        enddo
        if(lfound) then
          call shmem_real_get(cbsize(1),bsize(1,rem_block),3,rem_pe)
          call shmem_real_get(ccoord(1),coord(1,rem_block),3,rem_pe)
        else
          write(*,*) 'Paramesh error: amr_1blk_bcset : ',
     .      ' remote block is not in list received on this pe'
          call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
        endif
      endif
#endif /* MPI_USED */


! default values used if ndim < 3.
      dz = 0.
      z0 = 0.


      if(mod(nxb,2).eq.0) then
       if(ndim.eq.3) dz = cbsize(3)/real(nzb)
       dy = cbsize(2)/real(nyb)
       dx = cbsize(1)/real(nxb)
      else
       if(ndim.eq.3) dz = cbsize(3)/real(nzb-k3d)
       dy = cbsize(2)/real(nyb-1)
       dx = cbsize(1)/real(nxb-1)
      endif


!-------------------------

! Set index ranges
!
! unk
      idcc1 = bc_index_i(1,ibnd+2,1)
      idcc2 = bc_index_i(2,ibnd+2,1)
      jdcc1 = bc_index_j(1,jbnd+2,1)
      jdcc2 = bc_index_j(2,jbnd+2,1)
      kdcc1 = bc_index_k(1,kbnd+2,1)
      kdcc2 = bc_index_k(2,kbnd+2,1)
! facevar^s
      idfc1 = bc_index_i(1,ibnd+2,2)
      idfc2 = bc_index_i(2,ibnd+2,2)
      jdfc1 = bc_index_j(1,jbnd+2,2)
      jdfc2 = bc_index_j(2,jbnd+2,2)
      kdfc1 = bc_index_k(1,kbnd+2,2)
      kdfc2 = bc_index_k(2,kbnd+2,2)
! unk_e^s
      idec1 = bc_index_i(1,ibnd+2,3)
      idec2 = bc_index_i(2,ibnd+2,3)
      jdec1 = bc_index_j(1,jbnd+2,3)
      jdec2 = bc_index_j(2,jbnd+2,3)
      kdec1 = bc_index_k(1,kbnd+2,3)
      kdec2 = bc_index_k(2,kbnd+2,3)
! unk_n
      idnc1 = bc_index_i(1,ibnd+2,4)
      idnc2 = bc_index_i(2,ibnd+2,4)
      jdnc1 = bc_index_j(1,jbnd+2,4)
      jdnc2 = bc_index_j(2,jbnd+2,4)
      kdnc1 = bc_index_k(1,kbnd+2,4)
      kdnc2 = bc_index_k(2,kbnd+2,4)
      if(ibnd.eq.1) idnc1=idnc1+1
      if(ibnd.ge.0) idnc2=idnc2+1
      if(jbnd.eq.1) jdnc1=jdnc1+k2d
      if(jbnd.ge.0) jdnc2=jdnc2+k2d
      if(kbnd.eq.1) kdnc1=kdnc1+k3d
      if(kbnd.ge.0) kdnc2=kdnc2+k3d
! work
      idwc1 = bc_index_i(1,ibnd+2,5)
      idwc2 = bc_index_i(2,ibnd+2,5)
      jdwc1 = bc_index_j(1,jbnd+2,5)
      jdwc2 = bc_index_j(2,jbnd+2,5)
      kdwc1 = bc_index_k(1,kbnd+2,5)
      kdwc2 = bc_index_k(2,kbnd+2,5)

!-------------------------


! Which boundary condition has been specified?
      if(ibc.eq.-21) then

        if(nfacevar.gt.0.and.iopt.eq.1) then

!
! Do cell-face-centered data

!          facevarx1(:,idfc1:idfc2+1,jdfc1:jdfc2,kdfc1:kdfc2,idest) = 
!     .            ????


       l = idest

#ifdef FACEX
       do k=kdfc1,kdfc2
         if(ndim.eq.3) z0 = ccoord(3)-.5*(cbsize(3)+dz)
         zk = z0 + dz*real(k-nguard)
         do j=jdfc1,jdfc2
           y0 = ccoord(2)-.5*(cbsize(2)+dy)
           yj = y0 + dy*real(j-nguard)
           do i=idfc1,idfc2+1
             x0 = ccoord(1)-.5*cbsize(1)-dx
             xi = x0 + dx*real(i-nguard)
             do ivar=1,nbndvar
                value = ax*xi**interp_mask_facex(ivar) + 
     .                  ay*yj**interp_mask_facex(ivar) + 
     .                  az*zk**interp_mask_facex(ivar)
                facevarx1(ivar,i,j,k,l)=value*real(ivar)
             enddo
           enddo
         enddo
       enddo
#endif

#ifdef FACEY
!          facevary1(:,idfc1:idfc2,jdfc1:jdfc2+k2d,kdfc1:kdfc2,idest) = 
!     .            ????

       if(ndim.ge.2) then
       do k=kdfc1,kdfc2
         if(ndim.eq.3) z0 = ccoord(3)-.5*(cbsize(3)+dz)
         zk = z0 + dz*real(k-nguard)
         do j=jdfc1,jdfc2+k2d
           y0 = ccoord(2)-.5*cbsize(2)-dy
           yj = y0 + dy*real(j-nguard)
           do i=idfc1,idfc2
             x0 = ccoord(1)-.5*(cbsize(1)+dx)
             xi = x0 + dx*real(i-nguard)
             do ivar=1,nbndvar
                value = ax*xi**interp_mask_facey(ivar) + 
     .                  ay*yj**interp_mask_facey(ivar) + 
     .                  az*zk**interp_mask_facey(ivar)
                facevary1(ivar,i,j,k,l)=value*real(ivar)
             enddo
           enddo
         enddo
       enddo
       endif
#endif

#ifdef FACEZ
!          facevarz1(:,idfc1:idfc2,jdfc1:jdfc2,kdfc1:kdfc2+k3d,idest) = 
!     .            ????

       if(ndim.eq.3) then
       do k=kdfc1,kdfc2+k3d
         if(ndim.eq.3) z0 = ccoord(3)-.5*cbsize(3)-dz
         zk = z0 + dz*real(k-nguard)
         do j=jdfc1,jdfc2
           y0 = ccoord(2)-.5*(cbsize(2)+dy)
           yj = y0 + dy*real(j-nguard)
           do i=idfc1,idfc2
             x0 = ccoord(1)-.5*(cbsize(1)+dx)
             xi = x0 + dx*real(i-nguard)
             do ivar=1,nbndvar
                value = ax*xi**interp_mask_facez(ivar) + 
     .                  ay*yj**interp_mask_facez(ivar) + 
     .                  az*zk**interp_mask_facez(ivar)
                facevarz1(ivar,i,j,k,l)=value*real(ivar)
             enddo
           enddo
         enddo
       enddo
       endif
#endif

        endif                            ! end of nfacevar if test


!
! Now do cell centered data

        if(nvar.gt.0.and.iopt.eq.1) then

!          unk1(:,idcc1:idcc2,jdcc1:jdcc2,kdcc1:kdcc2,idest) = 
!     .             ????

          l = idest
          do k = kdcc1,kdcc2
          do j = jdcc1,jdcc2
          do i = idcc1,idcc2

            x0 = ccoord(1)-.5*(cbsize(1)+dx)+dx*real(i-nguard)
            y0 = ccoord(2)-.5*(cbsize(2)+dy)+dy*real(j-nguard)
            if(ndim.eq.3) z0 = 
     .          ccoord(3)-.5*(cbsize(3)+dz)+dz*real(k-nguard)

            do ivar=1,nvar
               if(int_gcell_on_cc(ivar)) then
               value = ax*x0**interp_mask_unk(ivar) + 
     .                 ay*y0**interp_mask_unk(ivar) + 
     .                 az*z0**interp_mask_unk(ivar)
               unk1(ivar,i,j,k,l) = value*real(ivar)
               endif
            enddo

          enddo
          enddo
          enddo

        endif


!!!#if N_DIM > 1
        if (ndim > 1) then
 
        if(nvaredge.gt.0.and.iopt.eq.1) then
!
! Now do cell edge centered data
 
 
!          unk_e_x1(:,idec1:idec2,jdec1:jdec2+k2d,kdec1:kdec2+k3d,idest) =
!     .             ????
 
          l = idest
          do k = kdec1,kdec2+k3d
          do j = jdec1,jdec2+k2d
          do i = idec1,idec2
 
            x0 = ccoord(1)-.5*(cbsize(1)+dx)+dx*real(i-nguard)
            y0 = ccoord(2)-.5*cbsize(2)-dy+dy*real(j-nguard)
            if(ndim.eq.3) z0 =
     .          ccoord(3)-.5*cbsize(3)-dz+dz*real(k-nguard)
 
            do ivar=1,nvaredge
              value = ax*x0**interp_mask_ec(ivar) + 
     .                ay*y0**interp_mask_ec(ivar) + 
     .                az*z0**interp_mask_ec(ivar)
              unk_e_x1(ivar,i,j,k,l) = value*real(ivar)
            enddo
 
          enddo
          enddo
          enddo
                           
!          unk_e_y1(:,idec1:idec2+1,jdec1:jdec2,kdec1:kdec2+k3d,idest) =
!     .             ????
 
          l = idest
          do k = kdec1,kdec2+k3d
          do j = jdec1,jdec2
          do i = idec1,idec2+1
 
            x0 = ccoord(1)-.5*cbsize(1)-dx+dx*real(i-nguard)
            y0 = ccoord(2)-.5*(cbsize(2)+dy)+dy*real(j-nguard)
            if(ndim.eq.3) z0 =
     .          ccoord(3)-.5*cbsize(3)-dz+dz*real(k-nguard)
 
            do ivar=1,nvaredge
              value = ax*x0**interp_mask_ec(ivar) + 
     .                ay*y0**interp_mask_ec(ivar) + 
     .                az*z0**interp_mask_ec(ivar)
              unk_e_y1(ivar,i,j,k,l) = value*real(ivar)
            enddo
 
          enddo
          enddo
          enddo
                           
!!!#if N_DIM == 3
          if (ndim == 3) then
!          unk_e_z1(:,idec1:idec2+1,jdec1:jdec2+k2d,kdec1:kdec2,idest) =
!     .             ????
 
          l = idest
          do k = kdec1,kdec2
          do j = jdec1,jdec2+k2d
          do i = idec1,idec2+1
 
            x0 = ccoord(1)-.5*cbsize(1)-dx+dx*real(i-nguard)
            y0 = ccoord(2)-.5*cbsize(2)-dy+dy*real(j-nguard)
            if(ndim.eq.3) z0 =
     .          ccoord(3)-.5*(cbsize(3)+dz)+dz*real(k-nguard)
 
            do ivar=1,nvaredge
              value = ax*x0**interp_mask_ec(ivar) + 
     .                ay*y0**interp_mask_ec(ivar) + 
     .                az*z0**interp_mask_ec(ivar)
              unk_e_z1(ivar,i,j,k,l) = value*real(ivar)
            enddo
 
          enddo
          enddo
          enddo
!!!#endif /* N_DIM == 3 */
          end if
                         
        endif
!!!#endif /* N_DIM > 1 */
        endif

                                                
        if(nvarcorn.gt.0.and.iopt.eq.1) then
!
! Now do cell corner data
!          unk_n1(:,idnc1:idnc2,jdnc1:jdnc2,kdnc1:kdnc2,idest) = 
!     .             ????

          l = idest
          do k = kdnc1,kdnc2
          do j = jdnc1,jdnc2
          do i = idnc1,idnc2

            x0 = ccoord(1)-.5*cbsize(1)-dx+dx*real(i-nguard)
            y0 = ccoord(2)-.5*cbsize(2)-dy+dy*real(j-nguard)
            if(ndim.eq.3) z0 = 
     .          ccoord(3)-.5*cbsize(3)-dz+dz*real(k-nguard)

            do ivar=1,nvarcorn
               value = ax*(x0**interp_mask_nc(ivar)) + 
     .                 ay*(y0**interp_mask_nc(ivar)) + 
     .                 az*(z0**interp_mask_nc(ivar))
               unk_n1(ivar,i,j,k,l) = value*real(ivar)
            enddo

          enddo
          enddo
          enddo

        endif

        if(iopt.ge.2) then

!          work1(idwc1:idwc2,jdwc1:jdwc2,kdwc1:kdwc2,idest) = 
!     .             ????

          l = idest
          do k = kdwc1,kdwc2
          do j = jdwc1,jdwc2
          do i = idwc1,idwc2

              x0=ccoord(1)-.5*(cbsize(1)+dx)+dx*real(i-nguard_work)
              y0=ccoord(2)-.5*(cbsize(2)+dy)+dy*real(j-nguard_work)
              if(ndim.eq.3) z0 =
     .          ccoord(3)-.5*(cbsize(3)+dz)+dz*real(k-nguard_work)

              value = ax*x0**interp_mask_work(1) + 
     .                ay*y0**interp_mask_work(1) + 
     .                az*z0**interp_mask_work(1)
              work1(i,j,k,l) = value

          enddo
          enddo
          enddo


        endif
!-------------------------
#ifdef MLR
! functions fnbx,fnby,fnbz are not available unless Rilee^s
! support routines are linked in.

! Which boundary condition has been specified?
      if(ibc.eq.-210) then

        if(nfacevar.gt.0.and.iopt.eq.1) then

!
! Do cell-face-centered data

!          facevarx1(:,idfc1:idfc2+1,jdfc1:jdfc2,kdfc1:kdfc2,idest) = 
!     .            ????


       l = idest

#ifdef FACEX
       do k=kdfc1,kdfc2
         if(ndim.eq.3) z0 = ccoord(3)-.5*(cbsize(3)+dz)
         zk = z0 + dz*real(k-nguard)
         do j=jdfc1,jdfc2
           y0 = ccoord(2)-.5*(cbsize(2)+dy)
           yj = y0 + dy*real(j-nguard)
           do i=idfc1,idfc2+1
             x0 = ccoord(1)-.5*cbsize(1)-dx
             xi = x0 + dx*real(i-nguard)
             do ivar=1,nbndvar
                value = fnbx(xi,yj,zk)
                facevarx1(ivar,i,j,k,l)=value
             enddo
           enddo
         enddo
       enddo
#endif

#ifdef FACEY
!          facevary1(:,idfc1:idfc2,jdfc1:jdfc2+k2d,kdfc1:kdfc2,idest) = 
!     .            ????

       if(ndim.ge.2) then
       do k=kdfc1,kdfc2
         if(ndim.eq.3) z0 = ccoord(3)-.5*(cbsize(3)+dz)
         zk = z0 + dz*real(k-nguard)
         do j=jdfc1,jdfc2+k2d
           y0 = ccoord(2)-.5*cbsize(2)-dy
           yj = y0 + dy*real(j-nguard)
           do i=idfc1,idfc2
             x0 = ccoord(1)-.5*(cbsize(1)+dx)
             xi = x0 + dx*real(i-nguard)
             do ivar=1,nbndvar
                value = fnby(xi,yj,zk)
                facevary1(ivar,i,j,k,l)=value
             enddo
           enddo
         enddo
       enddo
       endif
#endif

#ifdef FACEZ
!          facevarz1(:,idfc1:idfc2,jdfc1:jdfc2,kdfc1:kdfc2+k3d,idest) = 
!     .            ????

       if(ndim.eq.3) then
       do k=kdfc1,kdfc2+k3d
         if(ndim.eq.3) z0 = ccoord(3)-.5*cbsize(3)-dz
         zk = z0 + dz*real(k-nguard)
         do j=jdfc1,jdfc2
           y0 = ccoord(2)-.5*(cbsize(2)+dy)
           yj = y0 + dy*real(j-nguard)
           do i=idfc1,idfc2
             x0 = ccoord(1)-.5*(cbsize(1)+dx)
             xi = x0 + dx*real(i-nguard)
             do ivar=1,nbndvar
                value = fnby(xi,yj,zk)
                facevarz1(ivar,i,j,k,l)=value
             enddo
           enddo
         enddo
       enddo
       endif
#endif

        endif                            ! end of nfacevar if test


!
! Now do cell centered data

        if(nvar.gt.0.and.iopt.eq.1) then

!          unk1(:,idcc1:idcc2,jdcc1:jdcc2,kdcc1:kdcc2,idest) = 
!     .             ????

          l = idest
          do k = kdcc1,kdcc2
          do j = jdcc1,jdcc2
          do i = idcc1,idcc2

            x0 = ccoord(1)-.5*(cbsize(1)+dx)+dx*real(i-nguard)
            y0 = ccoord(2)-.5*(cbsize(2)+dy)+dy*real(j-nguard)
            if(ndim.eq.3) z0 = 
     .          ccoord(3)-.5*(cbsize(3)+dz)+dz*real(k-nguard)

            do ivar=1,nvar
              if(int_gcell_on_cc(ivar)) then
               value = fnbx(x0,y0,z0)
               unk1(ivar,i,j,k,l) = value
              endif
            enddo

          enddo
          enddo
          enddo

        endif


!!!#if N_DIM > 1
        if (ndim > 1) then
 
        if(nvaredge.gt.0.and.iopt.eq.1) then
!
! Now do cell edge centered data
 
 
!          unk_e_x1(:,idec1:idec2,jdec1:jdec2+k2d,kdec1:kdec2+k3d,idest) =
!     .             ????
 
          l = idest
          do k = kdec1,kdec2+k3d
          do j = jdec1,jdec2+k2d
          do i = idec1,idec2
 
            x0 = ccoord(1)-.5*(cbsize(1)+dx)+dx*real(i-nguard)
            y0 = ccoord(2)-.5*cbsize(2)-dy+dy*real(j-nguard)
            if(ndim.eq.3) z0 =
     .          ccoord(3)-.5*cbsize(3)-dz+dz*real(k-nguard)
 
            do ivar=1,nvaredge
               value = fnbx(x0,y0,z0)
              unk_e_x1(ivar,i,j,k,l) = value
            enddo
 
          enddo
          enddo
          enddo
                           
!          unk_e_y1(:,idec1:idec2+1,jdec1:jdec2,kdec1:kdec2+k3d,idest) =
!     .             ????
 
          l = idest
          do k = kdec1,kdec2+k3d
          do j = jdec1,jdec2
          do i = idec1,idec2+1
 
            x0 = ccoord(1)-.5*cbsize(1)-dx+dx*real(i-nguard)
            y0 = ccoord(2)-.5*(cbsize(2)+dy)+dy*real(j-nguard)
            if(ndim.eq.3) z0 =
     .          ccoord(3)-.5*cbsize(3)-dz+dz*real(k-nguard)
 
            do ivar=1,nvaredge
              value = fnby(x0,y0,z0)
              unk_e_y1(ivar,i,j,k,l) = value
            enddo
 
          enddo
          enddo
          enddo
                           
!!!#if N_DIM == 3
          if (ndim == 3) then
!          unk_e_z1(:,idec1:idec2+1,jdec1:jdec2+k2d,kdec1:kdec2,idest) =
!     .             ????
 
          l = idest
          do k = kdec1,kdec2
          do j = jdec1,jdec2+k2d
          do i = idec1,idec2+1
 
            x0 = ccoord(1)-.5*cbsize(1)-dx+dx*real(i-nguard)
            y0 = ccoord(2)-.5*cbsize(2)-dy+dy*real(j-nguard)
            if(ndim.eq.3) z0 =
     .          ccoord(3)-.5*(cbsize(3)+dz)+dz*real(k-nguard)
 
            do ivar=1,nvaredge
              value = fnbz(x0,y0,z0)
              unk_e_z1(ivar,i,j,k,l) = value
            enddo
 
          enddo
          enddo
          enddo
!!!#endif /* N_DIM == 3 */
          end if
                         
        endif
!!!#endif /* N_DIM > 1 */
        endif

                                                
        if(nvarcorn.gt.0.and.iopt.eq.1) then
!
! Now do cell corner data
!          unk_n1(:,idnc1:idnc2,jdnc1:jdnc2,kdnc1:kdnc2,idest) = 
!     .             ????

          l = idest
          do k = kdnc1,kdnc2
          do j = jdnc1,jdnc2
          do i = idnc1,idnc2

            x0 = ccoord(1)-.5*cbsize(1)-dx+dx*real(i-nguard)
            y0 = ccoord(2)-.5*cbsize(2)-dy+dy*real(j-nguard)
            if(ndim.eq.3) z0 = 
     .          ccoord(3)-.5*cbsize(3)-dz+dz*real(k-nguard)

            do ivar=1,nvarcorn
               value = fnbx(x0,y0,z0)
               unk_n1(ivar,i,j,k,l) = value
            enddo

          enddo
          enddo
          enddo

        endif

        if(iopt.ge.2) then

!          work1(idwc1:idwc2,jdwc1:jdwc2,kdwc1:kdwc2,idest) = 
!     .             ????

          l = idest
          do k = kdwc1,kdwc2
          do j = jdwc1,jdwc2
          do i = idwc1,idwc2

              x0=ccoord(1)-.5*(cbsize(1)+dx)+dx*real(i-nguard_work)
              y0=ccoord(2)-.5*(cbsize(2)+dy)+dy*real(j-nguard_work)
              if(ndim.eq.3) z0 =
     .          ccoord(3)-.5*(cbsize(3)+dz)+dz*real(k-nguard_work)

              value = fnbx(x0,y0,z0)
              work1(i,j,k,l) = value

          enddo
          enddo
          enddo


       endif ! iopt
       endif ! ibc.eq.-210
#endif /* MLR */
!-------------------------

!       elseif(ibc.eq.??) then


       endif                            ! end of test of bc flag

      return
      end
