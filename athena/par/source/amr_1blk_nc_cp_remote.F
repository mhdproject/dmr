#define GCELL_ON
#define GCELL_ON_FC
#define GCELL_ON_EC
#define GCELL_ON_NC
!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"
#define TESTING_RESTRICT
!#define DEBUG

      subroutine amr_1blk_nc_cp_remote(mype,remote_pe,remote_block,
     .   idest,id,jd,kd,is,js,ks,ilays,jlays,klays,ip1,jp1,kp1,
     .   ip3,jp3,kp3,nblk_ind)



!------------------------------------------------------------------------
!
! This routine copies guard cell information for cell corner
! data to layer idest of unk_n, from the appropriate
! corner data of the neighboring block.
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
!      ip1              offset added to index range defined by (id,ilay)
!                        0 if guardcells are at lower end of i index
!                        1 if guardcells are at upper end of i index
!      jp1              offset added to index range defined by (jd,jlay)
!                        0 if guardcells are at lower end of j index
!                        1 if guardcells are at upper end of j index
!      kp1              offset added to index range defined by (kd,klay)
!                        0 if guardcells are at lower end of k index
!                        1 if guardcells are at upper end of k index
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
      integer, intent(in) :: ip1,jp1,kp1,ip3,jp3,kp3
      integer, intent(in) :: nblk_ind


!-------------------------
#ifndef MPI_USED
! local arrays
      real recv(nbndvarc,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+k2d,
     .       kl_bnd:ku_bnd+k3d)
      save recv
#endif

      integer :: il,jl,kl,id1,jd1,kd1,is1,js1,ks1
      integer :: ill,jll,kll
      integer :: index, dtype 
      integer :: ii, jj, kk, i, j, k
      integer :: ia, ib, ja, jb, ka, kb
      integer :: ivar, ivar_next

#ifdef MPI_USED
      logical :: lfound
      integer :: vtype
#endif /* MPI_USED */
!-------------------------

#ifdef DEBUG
       if(l_f_to_c)
     . write(*,*) 'amr_1blk_nc_cp_remote : args ',
     .      mype,remote_pe,remote_block,
     .   idest,id,jd,kd,is,js,ks,ilays,jlays,klays,ip1,jp1,kp1
#endif /* DEBUG */

!
! Adjust index ranges
      il = ilays - ip3
      jl = jlays*k2d - jp3*k2d
      kl = klays*k3d - kp3*k3d

      id1 = id + ip1
      jd1 = jd + jp1*k2d
      kd1 = kd + kp1*k3d
      is1 = is + ip1
      js1 = js + jp1*k2d
      ks1 = ks + kp1*k3d

#ifndef MPI_USED

      recv = 0.

! Copy complete remote block into a buffer block called recv.
#ifdef NO_PERMANENT_GUARDCELLS
      call shmem_real_get(recv(1,1,1,1),
     .                   gt_unk_n(1,1,1,1,remote_block),
     .       len_blockn*nvarcorn,remote_pe) 
#else
      call shmem_real_get(recv(1,1,1,1),
     .                    unk_n(1,1,1,1,remote_block),
     .                    len_blockn*nvarcorn,remote_pe) 
#endif

#else /* MPI_USED */

!--
      if(remote_block.le.lnblocks.and.remote_pe.eq.mype) then
!--

#ifdef NO_PERMANENT_GUARDCELLS 

       if(.not.l_f_to_c) then
       unk_n1(1:nvarcorn,id1:id1+il,jd1:jd1+jl,
     .                                  kd1:kd1+kl,idest)
     .    =  gt_unk_n(1:nvarcorn,is1:is1+il,js1:js1+jl,
     .                           ks1:ks1+kl,
     .                           remote_block)
       else
       unk_n1_fl(1:nvarcorn,id1:id1+il,jd1:jd1+jl,
     .                                  kd1:kd1+kl)
     .    =  gt_unk_n(1:nvarcorn,is1:is1+il,js1:js1+jl,
     .                           ks1:ks1+kl,
     .                           remote_block)
       end if

#else

       if(.not.l_f_to_c) then
       unk_n1(1:nvarcorn,id1:id1+il,jd1:jd1+jl,
     .                                  kd1:kd1+kl,idest)
     .    =  unk_n(1:nvarcorn,is1:is1+il,js1:js1+jl,
     .                        ks1:ks1+kl,
     .                        remote_block)
       else
       unk_n1_fl(1:nvarcorn,id1:id1+il,jd1:jd1+jl,
     .                                  kd1:kd1+kl)
     .    =  unk_n(1:nvarcorn,is1:is1+il,js1:js1+jl,
     .                        ks1:ks1+kl,
     .                        remote_block)
       end if

#endif /* NO_PERMANENT_GUARDCELLS */

!--
      else                          ! otherwise if block is remote
!--

        call amr_mpi_find_blk_in_buffer(mype,remote_block,
     .                        remote_pe,idest,dtype,index,lfound)

#ifdef DEBUG
         write(*,*) 'pe ',mype,' find blk in buff ',remote_block,
     .                  remote_pe,' lfound ',lfound,' dtype ',dtype,
     .                  ' index ',index,' original address ',
     .                  laddress(:,remote_block)
#endif /* DEBUG */

! If this routine is executing a copy to fill guardcells of a
! leaf blocks^s parent, and the remote block is not found, then
! it is assumed that it is not in the list of buffered remote blocks
! because it is not really needed. Therefore in this case we
! return without copying anything.
        if(idest.eq.2.and.(.not.lfound)) return

! starting index if cell-centered data is also included in recv_buf
        if(l_datapacked(2)) index =
#ifdef GCELL_ON
     .               index + ngcell_on_cc*message_size_cc(dtype)
#else
     .                      index + nvar*message_size_cc(dtype)
#endif
        if(l_datapacked(3)) index =
#ifdef GCELL_ON_FC
     .                      index + maxval(ngcell_on_fc(1:ndim))
     .                              *message_size_fc(dtype)
#else
     .                      index + nfacevar*message_size_fc(dtype)
#endif
        if(l_datapacked(4)) index =
#ifdef GCELL_ON_EC
     .                      index + maxval(ngcell_on_ec(1:3))
     .                              *message_size_ec(dtype)
#else
     .                      index + nvaredge*message_size_ec(dtype)
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

        vtype = 8
        call mpi_set_message_limits(
     .               dtype,ia,ib,ja,jb,ka,kb,vtype,
     .               ill,jll,kll)



!       unk_n1(1:nvarcorn,id1:id1+il,jd1:jd1+jl,
!     .                                  kd1:kd1+kl,idest)
!     .    =  recv(1:nvarcorn,is1:is1+il,js1:js1+jl,
!     .                                       ks1:ks1+kl)


!       unk_n1_fl(1:nvarcorn,id1:id1+il,jd1:jd1+jl,
!     .                                  kd1:kd1+kl)
!     .    =  recv(1:nvarcorn,is1:is1+il,js1:js1+jl,
!     .                                       ks1:ks1+kl)

        kk = kd1
        do k = ka,kb
        jj = jd1
        do j = ja,jb
        ii = id1
        do i = ia,ib
          if (k >= ks1 .and. k <= ks1 + kl) then
          if (j >= js1 .and. j <= js1 + jl) then
          if (i >= is1 .and. i <= is1 + il) then

        do ivar=1,ngcell_on_nc
          ivar_next = gcell_on_nc_pointer(ivar)

          if (.not.l_f_to_c) then
          unk_n1(ivar_next,ii,jj,kk,idest) =
     .              temprecv_buf(index+ivar)
!pmn          unk_n1(1:nvarcorn,ii,jj,kk,idest) =
!pmn     .              temprecv_buf(index+1:index+nvarcorn)
          else
          unk_n1_fl(ivar_next,ii,jj,kk) =
     .              temprecv_buf(index+ivar)
!pmn          unk_n1_fl(1:nvarcorn,ii,jj,kk) =
!pmn     .              temprecv_buf(index+1:index+nvarcorn)
          end if
        enddo

          end if
          end if
          end if
          if (i >= is1 .and. i <= is1 + il) ii = ii + 1
#ifdef GCELL_ON_NC
          index = index+ngcell_on_nc
#else
          index = index+nvarcorn
#endif /* GCELL_ON_NC */`
#ifdef DEBUG
         write(*,*) 'pe ',mype,' unpacked recv for unk_n1 ',i,j,k,
     .      unk_n1(1:nvarcorn,ii,jj,kk,idest),' index ',index+1
#endif /* DEBUG */
        enddo
        if (j >= js1 .and. j <= js1 + jl) jj = jj + 1
        enddo
        if (k >= ks1 .and. k <= ks1 + kl) kk = kk + 1
        enddo


       endif
#endif /* MPI_USED */

#ifndef MPI_USED
       if(.not.l_f_to_c) then

       unk_n1(1:nvarcorn,id1:id1+il,jd1:jd1+jl,
     .                                  kd1:kd1+kl,idest)
     .    =  recv(1:nvarcorn,is1:is1+il,js1:js1+jl,
     .                                       ks1:ks1+kl)

       else

       unk_n1_fl(1:nvarcorn,id1:id1+il,jd1:jd1+jl,
     .                                  kd1:kd1+kl)
     .    =  recv(1:nvarcorn,is1:is1+il,js1:js1+jl,
     .                                       ks1:ks1+kl)

       end if
#endif /* MPI_USED */


#ifdef TESTING_RESTRICT

        if (l_f_to_c) then

        f2c_ind_unkn(1,1,nblk_ind) = min( id1,
     .                                   f2c_ind_unkn(1,1,nblk_ind))
        f2c_ind_unkn(2,1,nblk_ind) = max( id1+il,
     .                                   f2c_ind_unkn(2,1,nblk_ind))
        f2c_ind_unkn(1,2,nblk_ind) = min( jd1,
     .                                   f2c_ind_unkn(1,2,nblk_ind))
        f2c_ind_unkn(2,2,nblk_ind) = max( jd1+jl,
     .                                   f2c_ind_unkn(2,2,nblk_ind))
        f2c_ind_unkn(1,3,nblk_ind) = min( kd1,
     .                                   f2c_ind_unkn(1,3,nblk_ind))
        f2c_ind_unkn(2,3,nblk_ind) = max( kd1+kl,
     .                                   f2c_ind_unkn(2,3,nblk_ind))

       endif

#endif

      return
      end subroutine amr_1blk_nc_cp_remote
