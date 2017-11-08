#define GCELL_ON
#define GCELL_ON_FC
#define GCELL_ON_EC
!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#define TESTING_RESTRICT
#include "paramesh_preprocessor.fh"


      subroutine amr_1blk_ec_cp_remote(mype,remote_pe,remote_block,
     .   idest,id,jd,kd,is,js,ks,ilays,jlays,klays,ip1,jp1,kp1,
     .   ip2,jp2,kp2,ip3,jp3,kp3,iface,nblk_ind)



!------------------------------------------------------------------------
!
! This routine copies guard cell information for cell edge centered
! data to layer idest of unk_e_x, unk_e_y and unk_e_z, from the appropriate
! edge data of the neighboring block.
!
! Arguments:
!      mype             local processor
!      remote_pe        remote processor
!      remote_block     local block id of the block to be copied from
!                        the remote processor
!      idest            selects the storage space in data_1blk.fh which is to
!                        be used in this call. If the leaf node is having its
!                        guardcells filled then set this to 1, if its parent
!                        is being filled set it to 2.
!      id               lower limit of index range of points in x direction
!                        on destination block
!      jd               lower limit of index range of points in y direction
!                        on destination block
!      kd               lower limit of index range of points in z direction
!                        on destination block
!      is               lower limit of index range of points in x direction
!                        on source block
!      js               lower limit of index range of points in y direction
!                        on source block
!      ks               lower limit of index range of points in z direction
!                        on source block
!      ilay             no. of mesh points in x direction to be copied
!      jlay             no. of mesh points in y direction to be copied
!      klay             no. of mesh points in z direction to be copied
!      ip1              i index offset dependent on the face being treated,
!                        1 if at low or high x range, 0 otherwise.
!      jp1              j index offset dependent on the face being treated,
!                        1 if at low or high y range, 0 otherwise.
!      kp1              k index offset dependent on the face being treated,
!                        1 if at low or high z range, 0 otherwise.
!      ip2              extend range in i coord for unk_e_x by this amount
!                        must be set to either 1 or 0 
!                        (should be 1 for face 2, otherwise 0.)
!      jp2              extend range in j coord for unk_e_y by this amount
!                        must be set to either 1 or 0
!                        (should be 1 for face 4, otherwise 0.)
!      kp2              extend range in k coord for unk_e_z by this amount
!                        must be set to either 1 or 0
!                        (should be 1 for face 6, otherwise 0.)
!      iface            set between 1 and 6 if working on a block face,
!                        otherwise set to 0.
!
!
!
! Written :     Peter MacNeice          December 2000
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
!     use workspace
#ifdef MPI_USED
      use mpi_morton

      use paramesh_interfaces, only : amr_mpi_find_blk_in_buffer
      use paramesh_mpi_interfaces, only : mpi_set_message_limits
#endif /* MPI_USED */
      implicit none

#include "amr_shmem.fh"
      

!-------------------------

      integer, intent(in) :: mype,remote_pe,remote_block
      integer, intent(in) :: idest,id,jd,kd,is,js,ks
      integer, intent(in) :: ilays,jlays,klays
      integer, intent(in) :: ip1,jp1,kp1,ip2,jp2,kp2,ip3,jp3,kp3,iface
      integer, intent(in) :: nblk_ind

!-------------------------
! local arrays
#ifndef MPI_USED
      real recvx(nbndvare,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,
     .       kl_bnd:ku_bnd+k3d)
      real recvy(nbndvare,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,
     .       kl_bnd:ku_bnd+k3d)
      real recvz(nbndvare,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+k2d,
     .       kl_bnd:ku_bnd)
      save recvx,recvy,recvz
#endif

      integer :: il,jl,kl
      integer :: ill,jll,kll
      integer :: il1, jl1, kl1, id1, jd1, kd1
      integer :: is1, js1, ks1, index
      integer :: ii, jj, kk, i, j, k
      integer :: ia, ib, ja, jb, ka, kb
      integer :: ivar, ivar_next

#ifdef MPI_USED
      logical :: lfound
      integer :: dtype
      integer :: vtype
#endif /* MPI_USED */
!-------------------------

      if (ndim > 1) then

!
! Adjust index ranges
      il = ilays - ip3
      jl = jlays*k2d - jp3*k2d
      kl = klays*k3d - kp3*k3d

      il1 = il-1
      jl1 = (jl-1)*k2d
      kl1 = (kl-1)*k3d

      id1 = id + ip1
      jd1 = jd + jp1*k2d
      kd1 = kd + kp1*k3d
      is1 = is + ip1
      js1 = js + jp1*k2d
      ks1 = ks + kp1*k3d

#ifndef MPI_USED

      recvx = 0.
      recvy = 0.
      recvz = 0.

! Copy complete remote block into a buffer block called recv.
#ifdef NO_PERMANENT_GUARDCELLS
      call shmem_real_get(recvx(1,1,1,1),
     .                   gt_unk_e_x(1,1,1,1,remote_block),
     .       len_blockex*nbndvare,remote_pe) 
      call shmem_real_get(
     .                    recvy(1,1,1,1),
     .                    gt_unk_e_y(1,1,1,1,remote_block),
     .                    len_blockey*nbndvare,remote_pe) 
      if (ndim == 3) then
      call shmem_real_get(
     .                    recvz(1,1,1,1),
     .                    gt_unk_e_z(1,1,1,1,remote_block),
     .                    len_blockez*nbndvare,remote_pe) 
      end if

#else
      call shmem_real_get(recvx(1,1,1,1),
     .                    unk_e_x(1,1,1,1,remote_block),
     .                    len_blockex*nbndvare,remote_pe) 
      call shmem_real_get(
     .                    recvy(1,1,1,1),
     .                    unk_e_y(1,1,1,1,remote_block),
     .                    len_blockey*nbndvare,remote_pe) 
      if (ndim == 3) then
      call shmem_real_get(
     .                    recvz(1,1,1,1),
     .                    unk_e_z(1,1,1,1,remote_block),
     .                    len_blockez*nbndvare,remote_pe) 
      end if
#endif /* NO_PERMANENT_GUARDCELLS */

#else /* MPI_USED */

!--
      if(remote_block.le.lnblocks.and.remote_pe.eq.mype) then
!--

#ifdef NO_PERMANENT_GUARDCELLS
       if(.not.l_f_to_c) then
!-- 
       unk_e_x1(1:nbndvare,
     .          id:id+il1+ip2,
     .          jd1:jd1+jl,
     .          kd1:kd1+kl,
     .          idest)
     .    =  gt_unk_e_x(1:nbndvare,
     .                  is:is+il1+ip2,
     .                  js1:js1+jl,
     .                  ks1:ks1+kl,remote_block)

       unk_e_y1(1:nbndvare,
     .          id1:id1+il,
     .          jd:jd+jl1+jp2*k2d,
     .          kd1:kd1+kl,
     .          idest)
     .    =  gt_unk_e_y(1:nbndvare,
     .                  is1:is1+il,
     .                  js:js+jl1+jp2*k2d,
     .                  ks1:ks1+kl,remote_block)

       if (ndim == 3) then
       unk_e_z1(1:nbndvare,
     .          id1:id1+il,
     .          jd1:jd1+jl,
     .          kd:kd+kl1+kp2*k3d,
     .          idest)
     .    =  gt_unk_e_z(1:nbndvare,
     .                  is1:is1+il,
     .                  js1:js1+jl,
     .                  ks:ks+kl1+kp2*k3d,remote_block)
       end if

!-- 
      else
!-- 
       unk_e_x1_fl(1:nbndvare,
     .          id:id+il1+ip2,
     .          jd1:jd1+jl,
     .          kd1:kd1+kl)
     .    =  gt_unk_e_x(1:nbndvare,
     .                 is:is+il1+ip2,
     .                 js1:js1+jl,
     .                 ks1:ks1+kl,remote_block)

       unk_e_y1_fl(1:nbndvare,
     .          id1:id1+il,
     .          jd:jd+jl1+jp2*k2d,
     .          kd1:kd1+kl)
     .    =  gt_unk_e_y(1:nbndvare,
     .                  is1:is1+il,
     .                  js:js+jl1+jp2*k2d,
     .                  ks1:ks1+kl,remote_block)

       if (ndim == 3) then
       unk_e_z1_fl(1:nbndvare,
     .          id1:id1+il,
     .          jd1:jd1+jl,
     .          kd:kd+kl1+kp2*k3d)
     .    =  gt_unk_e_z(1:nbndvare,
     .                  is1:is1+il,
     .                  js1:js1+jl,
     .                  ks:ks+kl1+kp2*k3d,remote_block)
       end if

      end if

#else /* NO_PERMANENT_GUARDCELLS */

      if(.not.l_f_to_c) then
!-- 
       unk_e_x1(1:nbndvare,
     .          id:id+il1+ip2,
     .          jd1:jd1+jl,
     .          kd1:kd1+kl,
     .          idest)
     .    =  unk_e_x(1:nbndvare,
     .               is:is+il1+ip2,
     .               js1:js1+jl,
     .               ks1:ks1+kl,remote_block)

       unk_e_y1(1:nbndvare,
     .          id1:id1+il,
     .          jd:jd+jl1+jp2*k2d,
     .          kd1:kd1+kl,
     .          idest)
     .    =  unk_e_y(1:nbndvare,
     .               is1:is1+il,
     .               js:js+jl1+jp2*k2d,
     .               ks1:ks1+kl,remote_block)

       if (ndim == 3) then
       unk_e_z1(1:nbndvare,
     .          id1:id1+il,
     .          jd1:jd1+jl,
     .          kd:kd+kl1+kp2*k3d,
     .          idest)
     .    =  unk_e_z(1:nbndvare,
     .               is1:is1+il,
     .               js1:js1+jl,
     .               ks:ks+kl1+kp2*k3d,remote_block)
       end if

!-- 
      else
!-- 
       unk_e_x1_fl(1:nbndvare,
     .          id:id+il1+ip2,
     .          jd1:jd1+jl,
     .          kd1:kd1+kl)
     .    =  unk_e_x(1:nbndvare,
     .               is:is+il1+ip2,
     .               js1:js1+jl,
     .               ks1:ks1+kl,remote_block)

       unk_e_y1_fl(1:nbndvare,
     .          id1:id1+il,
     .          jd:jd+jl1+jp2*k2d,
     .          kd1:kd1+kl)
     .    =  unk_e_y(1:nbndvare,
     .               is1:is1+il,
     .               js:js+jl1+jp2*k2d,
     .               ks1:ks1+kl,remote_block)

       if (ndim == 3) then
       unk_e_z1_fl(1:nbndvare,
     .          id1:id1+il,
     .          jd1:jd1+jl,
     .          kd:kd+kl1+kp2*k3d)
     .    =  unk_e_z(1:nbndvare,
     .               is1:is1+il,
     .               js1:js1+jl,
     .               ks:ks+kl1+kp2*k3d,remote_block)
       end if

      end if

#endif /* NO_PERMANENT_GUARDCELLS */

!--
      else                          ! otherwise if block is remote
!--

        call amr_mpi_find_blk_in_buffer(mype,remote_block,
     .                        remote_pe,idest,dtype,index,lfound)


! If this routine is executing a copy to fill guardcells of a
! leaf blocks^s parent, and the remote block is not found, then
! it is assumed that it is not in the list of buffered remote blocks
! because it is not really needed. Therefore in this case we
! return without copying anything.
        if(idest.eq.2.and.(.not.lfound)) return

! starting index if cell-centered data is also included in recv_buf
        if(l_datapacked(2)) index =
#ifdef GCELL_ON
     .                      index + ngcell_on_cc*message_size_cc(dtype)
#else
     .                      index + nvar*message_size_cc(dtype)
#endif
        if(l_datapacked(3)) index =
#ifdef GCELL_ON_FC
     .                      index + maxval(ngcell_on_fc(1:ndim))
     .                                   *message_size_fc(dtype)
#else
     .                      index + nfacevar*message_size_fc(dtype)
#endif

        if (l_f_to_c) then
           if (ilays == 2*nguard) then
              ill = nguard
           else
              ill = (ilays-1)/2
           end if
           if (jlays == 2*nguard) then
              jll = nguard
           else
              jll = (jlays-1)/2
           end if
           if (klays == 2*nguard) then
              kll = nguard
           else
              kll = (klays-1)/2
           end if
        else
           ill = ilays
           jll = jlays
           kll = klays
        end if

        vtype = 5
        call mpi_set_message_limits(
     .               dtype,ia,ib,ja,jb,ka,kb,vtype,
     .               ill,jll,kll)

        kk = kd1
        do k = ka,kb
        jj = jd1
        do j = ja,jb
        ii = id
        do i = ia,ib
          if (k >= ks1 .and. k <= ks1 + kl) then
          if (j >= js1 .and. j <= js1 + jl) then
          if (i >= is .and.  i <= is + il1 + ip2) then

        do ivar=1,ngcell_on_ec(1)
          ivar_next = gcell_on_ec_pointer(1,ivar)

          if (.not.l_f_to_c) then
           unk_e_x1(ivar_next,ii,jj,kk,idest) =
     .              temprecv_buf(index+ivar)
!pmn           unk_e_x1(1:nbndvare,ii,jj,kk,idest) =
!pmn     .              temprecv_buf(index+1:index+nbndvare)
          else
           unk_e_x1_fl(ivar_next,ii,jj,kk) =
     .                 temprecv_buf(index+ivar)
!pmn           unk_e_x1_fl(1:nbndvare,ii,jj,kk) =
!pmn     .                 temprecv_buf(index+1:index+nbndvare)
          endif

        enddo

          endif
          endif
          endif
          if (i >= is .and. i <= is + il1 + ip2) ii = ii + 1
#ifdef GCELL_ON_EC
          index = index+ngcell_on_ec(1)
#else
          index = index+nbndvare
#endif /* GCELL_ON_EC */
        enddo
        if (j >= js1 .and. j <= js1 + jl) jj = jj + 1
        enddo
        if (k >= ks1 .and. k <= ks1 + kl) kk = kk + 1
        enddo

       if(ndim.ge.2) then
        vtype = 6
        call mpi_set_message_limits(
     .               dtype,ia,ib,ja,jb,ka,kb,vtype,
     .               ill,jll,kll)

        kk = kd1
        do k = ka,kb
        jj = jd
        do j = ja,jb
        ii = id1
        do i = ia,ib
          if (k >= ks1 .and. k <= ks1 + kl) then
          if (j >= js  .and. j <= js  + jl1 + jp2*k2d) then
          if (i >= is1 .and. i <= is1 + il) then

        do ivar=1,ngcell_on_ec(2)
          ivar_next = gcell_on_ec_pointer(2,ivar)

          if (.not.l_f_to_c) then
           unk_e_y1(ivar_next,ii,jj,kk,idest) =
     .              temprecv_buf(index+ivar)
!pmn           unk_e_y1(1:nbndvare,ii,jj,kk,idest) =
!pmn     .              temprecv_buf(index+1:index+nbndvare)
          else
           unk_e_y1_fl(ivar_next,ii,jj,kk) =
     .                 temprecv_buf(index+ivar)
!pmn           unk_e_y1_fl(1:nbndvare,ii,jj,kk) =
!pmn     .                 temprecv_buf(index+1:index+nbndvare)
          endif

        enddo

          endif
          endif
          endif
          if (i >= is1 .and.  i <= is1 + il) ii = ii + 1
#ifdef GCELL_ON_EC
          index = index+ngcell_on_ec(2)
#else
          index = index+nbndvare
#endif /* GCELL_ON_EC */
        enddo
        if (j >= js  .and.  j <= js + jl1 + jp2*k2d) jj = jj + 1
        enddo
        if (k >= ks1 .and. k <= ks1 + kl) kk = kk + 1
        enddo
       endif

       if((ndim == 3).or.(l2p5d.eq.1)) then
        vtype =7 
        call mpi_set_message_limits(
     .               dtype,ia,ib,ja,jb,ka,kb,vtype,
     .               ill,jll,kll)

        kk = kd
        do k = ka,kb
        jj = jd1
        do j = ja,jb
        ii = id1
        do i = ia,ib
          if (k >= ks  .and.  k <= ks + kl1 + kp2*k3d) then
          if (j >= js1 .and.  j <= js1 + jl) then
          if (i >= is1 .and.  i <= is1 + il) then

        do ivar=1,ngcell_on_ec(3)
          ivar_next = gcell_on_ec_pointer(3,ivar)

          if (.not.l_f_to_c) then
           unk_e_z1(ivar_next,ii,jj,kk,idest) =
     .              temprecv_buf(index+ivar)
!pmn           unk_e_z1(1:nbndvare,ii,jj,kk,idest) =
!pmn     .              temprecv_buf(index+1:index+nbndvare)
          else
           unk_e_z1_fl(ivar_next,ii,jj,kk) =
     .                 temprecv_buf(index+ivar)
!pmn           unk_e_z1_fl(1:nbndvare,ii,jj,kk) =
!pmn     .                 temprecv_buf(index+1:index+nbndvare)
          endif

        enddo

          endif
          endif
          endif
          if (i >= is1 .and.  i <= is1 + il) ii = ii + 1
#ifdef GCELL_ON_EC
          index = index+ngcell_on_ec(3)
#else
          index = index+nbndvare
#endif /* GCELL_ON_EC */
        enddo
        if (j >= js1 .and.  j <= js1 + jl) jj = jj + 1
        enddo
        if (k >= ks .and. k <= ks + kl1 + kp2*k3d) kk = kk + 1
        enddo
       endif

       endif
#endif /* MPI_USED */

#ifndef MPI_USED

       if(.not.l_f_to_c) then

!-- 
       unk_e_x1(1:nbndvare,
     .          id:id+il1+ip2,
     .          jd1:jd1+jl,
     .          kd1:kd1+kl,
     .          idest)
     .    =  recvx(1:nbndvare,
     .             is:is+il1+ip2,
     .             js1:js1+jl,
     .             ks1:ks1+kl)

       unk_e_y1(1:nbndvare,
     .          id1:id1+il,
     .          jd:jd+jl1+jp2*k2d,
     .          kd1:kd1+kl,
     .          idest)
     .    =  recvy(1:nbndvare,
     .             is1:is1+il,
     .             js:js+jl1+jp2*k2d,
     .             ks1:ks1+kl)


       if (ndim == 3) then
       unk_e_z1(1:nbndvare,
     .          id1:id1+il,
     .          jd1:jd1+jl,
     .          kd:kd+kl1+kp2*k3d,
     .          idest)
     .    =  recvz(1:nbndvare,
     .             is1:is1+il,
     .             js1:js1+jl,
     .             ks:ks+kl1+kp2*k3d)
       end if

!-- 
      else
!-- 
       unk_e_x1_fl(1:nbndvare,
     .          id:id+il1+ip2,
     .          jd1:jd1+jl,
     .          kd1:kd1+kl)
     .    =  recvx(1:nbndvare,
     .             is:is+il1+ip2,
     .             js1:js1+jl,
     .             ks1:ks1+kl)

       unk_e_y1_fl(1:nbndvare,
     .          id1:id1+il,
     .          jd:jd+jl1+jp2*k2d,
     .          kd1:kd1+kl)
     .    =  recvy(1:nbndvare,
     .             is1:is1+il,
     .             js:js+jl1+jp2*k2d,
     .             ks1:ks1+kl)


       if (ndim == 3) then
       unk_e_z1_fl(1:nbndvare,
     .          id1:id1+il,
     .          jd1:jd1+jl,
     .          kd:kd+kl1+kp2*k3d)
     .    =  recvz(1:nbndvare,
     .             is1:is1+il,
     .             js1:js1+jl,
     .             ks:ks+kl1+kp2*k3d)
       end if

       end if ! if (.not.l_f_to_c

#endif /* MPI_USED */

#ifdef TESTING_RESTRICT

       if(l_f_to_c) then

        f2c_ind_unkex(1,1,nblk_ind) = min( id,
     .                                   f2c_ind_unkex(1,1,nblk_ind))
        f2c_ind_unkex(2,1,nblk_ind) = max( id+il1+ip2,
     .                                   f2c_ind_unkex(2,1,nblk_ind))
        f2c_ind_unkex(1,2,nblk_ind) = min( jd1,
     .                                   f2c_ind_unkex(1,2,nblk_ind))
        f2c_ind_unkex(2,2,nblk_ind) = max( jd1+jl,
     .                                   f2c_ind_unkex(2,2,nblk_ind))
        f2c_ind_unkex(1,3,nblk_ind) = min( kd1,
     .                                   f2c_ind_unkex(1,3,nblk_ind))
        f2c_ind_unkex(2,3,nblk_ind) = max( kd1+kl,
     .                                   f2c_ind_unkex(2,3,nblk_ind))

        f2c_ind_unkey(1,1,nblk_ind) = min( id1,
     .                                   f2c_ind_unkey(1,1,nblk_ind))
        f2c_ind_unkey(2,1,nblk_ind) = max( id1+il,
     .                                   f2c_ind_unkey(2,1,nblk_ind))
        f2c_ind_unkey(1,2,nblk_ind) = min( jd,
     .                                   f2c_ind_unkey(1,2,nblk_ind))
        f2c_ind_unkey(2,2,nblk_ind) = max( jd+jl1+jp2,
     .                                   f2c_ind_unkey(2,2,nblk_ind))
        f2c_ind_unkey(1,3,nblk_ind) = min( kd1,
     .                                   f2c_ind_unkey(1,3,nblk_ind))
        f2c_ind_unkey(2,3,nblk_ind) = max( kd1+kl,
     .                                   f2c_ind_unkey(2,3,nblk_ind))

        f2c_ind_unkez(1,1,nblk_ind) = min( id1,
     .                                   f2c_ind_unkez(1,1,nblk_ind))
        f2c_ind_unkez(2,1,nblk_ind) = max( id1+il,
     .                                   f2c_ind_unkez(2,1,nblk_ind))
        f2c_ind_unkez(1,2,nblk_ind) = min( jd1,
     .                                   f2c_ind_unkez(1,2,nblk_ind))
        f2c_ind_unkez(2,2,nblk_ind) = max( jd1+jl,
     .                                   f2c_ind_unkez(2,2,nblk_ind))
        f2c_ind_unkez(1,3,nblk_ind) = min( kd,
     .                                   f2c_ind_unkez(1,3,nblk_ind))
        f2c_ind_unkez(2,3,nblk_ind) = max( kd+kl1+kp2,
     .                                   f2c_ind_unkez(2,3,nblk_ind))

!-- 
      endif  ! if (l_f_to_c

#endif

#ifdef XFORCE_CONSISTENCY_AT_SRL_INTERFACES
! This is not yet complete. Need to cater for cases where 4 blocks share
! a common edge.

       if(iface.eq.1) then

        do ivar=1,ngcell_on_ec(2)
          ivar_next = gcell_on_ec_pointer(2,ivar)
         unk_e_y1(ivar_next,1+nguard,1+nguard*k2d:nyb+nguard*k2d,
     .            1+nguard*k3d:nzb+(1+nguard)*k3d,idest) = .5*(
     .   unk_e_y1(ivar_next,1+nguard,1+nguard*k2d:nyb+nguard*k2d,
     .                       1+nguard*k3d:nzb+(1+nguard)*k3d,idest)
     .  + recvy(ivar_next,nxb+1,1:nyb,1:nzb+k3d) )
        enddo

         if (ndim == 3) then
        do ivar=1,ngcell_on_ec(3)
          ivar_next = gcell_on_ec_pointer(3,ivar)
         unk_e_z1(ivar_next,1+nguard,1+nguard*k2d:nyb+(1+nguard)*k2d,
     .            1+nguard*k3d:nzb+nguard*k3d,idest) = .5*(
     .   unk_e_z1(ivar_next,1+nguard,1+nguard*k2d:nyb+(1+nguard)*k2d,
     .            1+nguard*k3d:nzb+nguard*k3d,idest) 
     .  + recvz(ivar_next,nxb+1,1:nyb+k2d,1:nzb) )
        enddo
         end if

       elseif(iface.eq.2) then

        do ivar=1,ngcell_on_ec(2)
          ivar_next = gcell_on_ec_pointer(2,ivar)
         unk_e_y1(ivar_next,nxb+1+nguard,1+nguard*k2d:nyb+nguard*k2d,
     .            1+nguard*k3d:nzb+(1+nguard)*k3d,idest) = .5*(
     .   unk_e_y1(ivar_next,nxb+1+nguard,1+nguard*k2d:nyb+nguard*k2d,
     .            1+nguard*k3d:nzb+(1+nguard)*k3d,idest)
     .  + recvy(ivar_next,1,1:nyb,1:nzb+k3d) )
        enddo

         if (ndim == 3) then
        do ivar=1,ngcell_on_ec(3)
          ivar_next = gcell_on_ec_pointer(3,ivar)
         unk_e_z1(ivar_next,nxb+1+nguard,
     .            1+nguard*k2d:nyb+(1+nguard)*k2d,
     .            1+nguard*k3d:nzb+nguard*k3d,idest) = .5*(
     .   unk_e_z1(ivar_next,nxb+1+nguard,
     .            1+nguard*k2d:nyb+(1+nguard)*k2d,
     .            1+nguard*k3d:nzb+nguard*k3d,idest) 
     .  + recvz(ivar_next,1,1:nyb+k2d,1:nzb) )
        enddo
         end if

       elseif(iface.eq.3) then
        do ivar=1,ngcell_on_ec(1)
          ivar_next = gcell_on_ec_pointer(1,ivar)
         unk_e_x1(ivar_next,1+nguard:nxb+nguard,1+nguard*k2d,
     .            1+nguard*k3d:nzb+(1+nguard)*k3d,idest) = .5*(
     .   unk_e_x1(ivar_next,1+nguard:nxb+nguard,1+nguard*k2d,
     .            1+nguard*k3d:nzb+(1+nguard)*k3d,idest)
     .  + recvx(ivar_next,1:nxb,nyb+k2d,1:nzb+k3d) )
        enddo

         if (ndim == 3) then
        do ivar=1,ngcell_on_ec(3)
          ivar_next = gcell_on_ec_pointer(3,ivar)
         unk_e_z1(ivar_next,1+nguard:nxb+nguard+1,1+nguard*k2d,
     .            1+nguard*k3d:nzb+nguard*k3d,idest) = .5*(
     .   unk_e_z1(ivar_next,1+nguard:nxb+nguard+1,1+nguard*k2d,
     .            1+nguard*k3d:nzb+nguard*k3d,idest)
     .  + recvz(ivar_next,1:nxb+1,1,1:nzb) )
        enddo
         end if

       elseif(iface.eq.4) then
        do ivar=1,ngcell_on_ec(1)
          ivar_next = gcell_on_ec_pointer(1,ivar)
         unk_e_x1(ivar_next,1+nguard:nxb+nguard,nyb+(1+nguard)*k2d,
     .            1+nguard*k3d:nzb+(1+nguard)*k3d,idest) = .5*(
     .   unk_e_x1(ivar_next,1+nguard:nxb+nguard,nyb+(1+nguard)*k2d,
     .            1+nguard*k3d:nzb+(1+nguard)*k3d,idest)
     .  + recvx(ivar_next,1:nxb,1,1:nzb+k3d) )
        enddo

         if (ndim == 3) then
        do ivar=1,ngcell_on_ec(3)
          ivar_next = gcell_on_ec_pointer(3,ivar)
         unk_e_z1(ivar_next,1+nguard:nxb+nguard+1,nyb+(1+nguard)*k2d,
     .            1+nguard*k3d:nzb+nguard*k3d,idest) = .5*(
     .   unk_e_z1(ivar_next,1+nguard:nxb+nguard+1,nyb+(1+nguard)*k2d,
     .            1+nguard*k3d:nzb+nguard*k3d,idest)
     .  + recvz(ivar_next,1:nxb+1,1,1:nzb) )
        enddo
         end if

       elseif(iface.eq.5 .and. ndim == 3) then
        do ivar=1,ngcell_on_ec(1)
          ivar_next = gcell_on_ec_pointer(1,ivar)
         unk_e_x1(ivar_next,1+nguard:nxb+nguard,
     .            1+nguard*k2d:nyb+(1+nguard)*k2d,
     .            1+nguard*k3d,idest) = .5*(
     .   unk_e_x1(ivar_next,1+nguard:nxb+nguard,
     .            1+nguard*k2d:nyb+(1+nguard)*k2d,
     .            1+nguard*k3d,idest) 
     .  + recvx(ivar_next,1:nxb,1:nyb+k2d,nzb) )
        enddo

        do ivar=1,ngcell_on_ec(2)
          ivar_next = gcell_on_ec_pointer(2,ivar)
         unk_e_y1(ivar_next,1+nguard:nxb+nguard+1,
     .            1+nguard*k2d:nyb+nguard*k2d,
     .            1+nguard*k3d,idest) = .5*(
     .   unk_e_y1(ivar_next,1+nguard:nxb+nguard+1,
     .            1+nguard*k2d:nyb+nguard*k2d,
     .            1+nguard*k3d,idest) 
     .  + recvy(ivar_next,1:nxb+1,1:nyb,nzb) )
        enddo

       elseif(iface.eq.6 .and. ndim == 3) then
        do ivar=1,ngcell_on_ec(1)
          ivar_next = gcell_on_ec_pointer(1,ivar)
         unk_e_x1(ivar_next,1+nguard:nxb+nguard,
     .            1+nguard*k2d:nyb+(1+nguard)*k2d,
     .            nzb+(1+nguard)*k3d,idest) = .5*(
     .   unk_e_x1(ivar_next,1+nguard:nxb+nguard,
     .            1+nguard*k2d:nyb+(1+nguard)*k2d,
     .            nzb+(1+nguard)*k3d,idest)
     .  + recvx(ivar_next,1:nxb,1:nyb+k2d,1) )
        enddo

        do ivar=1,ngcell_on_ec(2)
          ivar_next = gcell_on_ec_pointer(2,ivar)
         unk_e_y1(ivar_next,1+nguard:nxb+nguard+1,
     .            1+nguard*k2d:nyb+nguard*k2d,
     .            nzb+(1+nguard)*k3d,idest) = .5*(
     .   unk_e_y1(ivar_next,1+nguard:nxb+nguard+1,
     .            1+nguard*k2d:nyb+nguard*k2d,
     .            nzb+(1+nguard)*k3d,idest)
     .  + recvy(ivar_next,1:nxb+1,1:nyb,1) )
        enddo

       endif
#endif /* FORCE_CONSISTENCY_AT_SRL_INTERFACES */

       end if ! if (ndim > 1

      return
      end subroutine amr_1blk_ec_cp_remote


