!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"


      subroutine amr_perm_to_1blk( lcc,lfc,lec,lnc,lb,pe,iopt,idest)



!------------------------------------------------------------------------
!
! This routine copies data to the 1-block working arrays with guardcells
! from the permanent data arrays, which may or may not have permanent
! guardcells, depending on whether NO_PERMANENT_GUARDCELLS is defined 
! in physicaldata.fh.
!
!
! Arguments :
!      lcc          logical       copies cell centered data if true
!      lfc          logical       copies cell face-centered data if true
!      lec          logical       copies cell edge-centered data if true
!      lnc          logical       copies cell corner data if true
!      lb           integer       block from which data is to be copied
!      pe           integer       processor from which data is to be copied
!      iopt         integer       data structure to be copied
!      idest        integer       sets value for last dimension index
!                                  in the 1-blk data arrays
!
!
! Written :     Peter MacNeice          February 1999
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace
#ifdef MPI_USED
      use mpi_morton
      use paramesh_mpi_interfaces, only : mpi_set_message_limits
      use paramesh_interfaces, only : amr_mpi_find_blk_in_buffer
#endif /* MPI_USED */

      implicit none

#include "amr_shmem.fh"

      integer, intent(in) ::  lb,pe,iopt,idest
      logical, intent(in) ::  lcc,lfc,lec,lnc


!------------------------------------
! local variables

#ifndef MPI_USED
      real recv(nvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
      real recvx(nbndvar,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,
     .       kl_bnd:ku_bnd)
      real recvy(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,
     .       kl_bnd:ku_bnd)
      real recvz(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,
     .       kl_bnd:ku_bnd+k3d)
      real recvex(nbndvare,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,
     .       kl_bnd:ku_bnd+k3d)
      real recvey(nbndvare,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,
     .       kl_bnd:ku_bnd+k3d)
      real recvez(nbndvare,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+k2d,
     .       kl_bnd:ku_bnd)
      real recvn(nbndvarc,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+k2d,
     .       kl_bnd:ku_bnd+k3d)
      save recv,recvn
      save recvx,recvy,recvz
      save recvex,recvey,recvez
#endif

      integer :: iopt0
      integer :: nguard0,nguard_work0
      integer :: ia, ib, ja, jb, ka, kb
      integer :: i, j, k, ivar, ivar_next
      parameter(nguard0 = nguard*npgs, nguard_work0 = nguard_work*npgs)
#ifdef MPI_USED
      integer :: vtype,dtype,rem_blk,rem_pe,mype,shmem_my_pe
      integer :: index,index0
      logical :: lfound
#endif /* MPI_USED */

!------------------------------------

#ifdef MPI_USED
         if(lb.gt.lnblocks) then
           mype = shmem_my_pe()
           rem_blk = lb
           rem_pe  = mype
           call amr_mpi_find_blk_in_buffer(mype,rem_blk,
     .                       rem_pe,idest,dtype,index0,lfound)

!           if(.not.lfound) 
!      write(*,*) 
!     .              'perm to 1blk reporting blk not found'
         endif


#endif /* MPI_USED */

!
! Put block lb's data into the data_1blk.fh datastructures, with the
! appropriate guardcell padding.
          if(iopt.eq.1) then


          if(lcc) then

#ifdef MPI_USED
          if(lb.le.lnblocks) then

            do ivar=1,nvar
              if(int_gcell_on_cc(ivar)) then
#ifdef NO_PERMANENT_GUARDCELLS
            unk1(ivar,1+nguard:nxb+nguard,1+nguard*k2d:nyb+nguard*k2d,
     .            1+nguard*k3d:nzb+nguard*k3d,idest) =
     .  unk(ivar,1+nguard0:nxb+nguard0,1+nguard0*k2d:nyb+nguard0*k2d,
     .            1+nguard0*k3d:nzb+nguard0*k3d,lb)
#else
            unk1(ivar,1+nguard:nxb+nguard,1+nguard*k2d:nyb+nguard*k2d,
     .            1+nguard*k3d:nzb+nguard*k3d,idest) =
     .  unk(ivar,1+nguard:nxb+nguard,1+nguard*k2d:nyb+nguard*k2d,
     .            1+nguard*k3d:nzb+nguard*k3d,lb)
#endif /* NO_PERMANENT_GUARDCELLS */
              endif
            enddo

          elseif(lb.gt.lnblocks) then

          vtype = 1
          call mpi_set_message_limits(
     .                 dtype,ia,ib,ja,jb,ka,kb,vtype)
          index = index0

#ifdef NO_PERMANENT_GUARDCELLS
          ia = ia + nguard
          ib = ib + nguard
          ja = ja + nguard*k2d
          jb = jb + nguard*k2d
          ka = ka + nguard*k3d
          kb = kb + nguard*k3d
#endif /* NO_PERMANENT_GUARDCELLS */

          do k = ka,kb
          do j = ja,jb
          do i = ia,ib
#ifdef OLD
            do ivar=1,nvar
              if(int_gcell_on_cc(ivar)) then
                unk1(ivar,i,j,k,idest) = 
     .              temprecv_buf(index+ivar)
              endif
            enddo
            index = index+nvar
#else
            do ivar=1,ngcell_on_cc
              ivar_next = gcell_on_cc_pointer(ivar)
              unk1(ivar_next,i,j,k,idest) = temprecv_buf(index+ivar)
            enddo
            index = index+ngcell_on_cc
#endif /* OLD */
          enddo
          enddo
          enddo

          endif

#else /* MPI_USED */

          call shmem_real_get(recv(1,1,1,1),
     .                        unk(1,1,1,1,lb),len_block,pe)
          do ivar=1,nvar
            if(int_gcell_on_cc(ivar)) then
              unk1(ivar,:,:,:,idest) = 0.
#ifdef NO_PERMANENT_GUARDCELLS
          unk1(ivar,1+nguard:nxb+nguard,1+nguard*k2d:nyb+nguard*k2d,
     .           1+nguard*k3d:nzb+nguard*k3d,idest) =
     .  recv(ivar,1+nguard0:nxb+nguard0,1+nguard0*k2d:nyb+nguard0*k2d,
     .            1+nguard0*k3d:nzb+nguard0*k3d)
#else
          unk1(ivar,:,:,:,idest) = recv(ivar,:,:,:)
#endif
            endif
          enddo

#endif /* MPI_USED */

          endif ! end if (lcc

          if(lfc) then

#ifdef MPI_USED

            if(lb.le.lnblocks) then

            do ivar=1,nfacevar
              if(int_gcell_on_fc(1,ivar)) then
#ifdef NO_PERMANENT_GUARDCELLS
            facevarx1(ivar,1+nguard:nxb+nguard+1,
     .                  1+nguard*k2d:nyb+nguard*k2d,
     .                  1+nguard*k3d:nzb+nguard*k3d,idest) =
     .       facevarx(ivar,1+nguard0:nxb+nguard0+1,
     .                  1+nguard0*k2d:nyb+nguard0*k2d,
     .                  1+nguard0*k3d:nzb+nguard0*k3d,lb)
#else
            facevarx1(ivar,1+nguard:nxb+nguard+1,
     .                  1+nguard*k2d:nyb+nguard*k2d,
     .                  1+nguard*k3d:nzb+nguard*k3d,idest) =
     .       facevarx(ivar,1+nguard:nxb+nguard+1,
     .                  1+nguard*k2d:nyb+nguard*k2d,
     .                  1+nguard*k3d:nzb+nguard*k3d,lb)
#endif
              endif
            enddo

          elseif(lb.gt.lnblocks) then

! starting index if cell-centered data is also included in recv_buf
#ifdef OLD
          index = index0 + nvar*message_size_cc(dtype)
#else
          index = index0 + ngcell_on_cc*message_size_cc(dtype)
#endif /* OLD */

          vtype = 2
          call mpi_set_message_limits(
     .                 dtype,ia,ib,ja,jb,ka,kb,vtype)

#ifdef NO_PERMANENT_GUARDCELLS
          ia = ia + nguard
          ib = ib + nguard
          ja = ja + nguard*k2d
          jb = jb + nguard*k2d
          ka = ka + nguard*k3d
          kb = kb + nguard*k3d
#endif
          do k = ka,kb
          do j = ja,jb
          do i = ia,ib
#ifdef OLD
            do ivar=1,nfacevar
              if(int_gcell_on_fc(1,ivar)) then
            facevarx1(ivar,i,j,k,idest) = 
     .                 temprecv_buf(index+ivar)
              endif
            enddo
            index = index+nfacevar
#else
            do ivar=1,ngcell_on_fc(1)
              ivar_next = gcell_on_fc_pointer(1,ivar)
              facevarx1(ivar_next,i,j,k,idest) = 
     .                 temprecv_buf(index+ivar)
            enddo
            index = index+ngcell_on_fc(1)
#endif /* OLD */
          enddo
          enddo
          enddo

          endif

#else /* MPI_USED */

            call shmem_real_get(recvx(1,1,1,1),
     .                            facevarx(1,1,1,1,lb),
     .                            nbndvar*len_blockfx,pe)

            do ivar=1,nfacevar
              if(int_gcell_on_fc(1,ivar)) then
            facevarx1(ivar,:,:,:,idest) = 0.
#ifdef NO_PERMANENT_GUARDCELLS
            facevarx1(ivar,1+nguard:nxb+nguard+1,
     .                  1+nguard*k2d:nyb+nguard*k2d,
     .                  1+nguard*k3d:nzb+nguard*k3d,idest) =
     .       recvx(ivar,1+nguard0:nxb+nguard0+1,
     .                  1+nguard0*k2d:nyb+nguard0*k2d,
     .                  1+nguard0*k3d:nzb+nguard0*k3d)
#else
            facevarx1(ivar,:,:,:,idest) = recvx(ivar,:,:,:)
#endif
              endif
            enddo

#endif /* MPI_USED */

            if(ndim.ge.2) then

#ifdef MPI_USED

            if(lb.le.lnblocks) then

            do ivar=1,nfacevar
              if(int_gcell_on_fc(2,ivar)) then
#ifdef NO_PERMANENT_GUARDCELLS
              facevary1(ivar,1+nguard:nxb+nguard,
     .                   1+nguard*k2d:nyb+(nguard+1)*k2d,
     .                   1+nguard*k3d:nzb+nguard*k3d,idest) =
     .        facevary(ivar,1+nguard0:nxb+nguard0,
     .                   1+nguard0*k2d:nyb+(nguard0+1)*k2d,
     .                   1+nguard0*k3d:nzb+nguard0*k3d,lb)
#else
              facevary1(ivar,1+nguard:nxb+nguard,
     .                   1+nguard*k2d:nyb+(nguard+1)*k2d,
     .                   1+nguard*k3d:nzb+nguard*k3d,idest) =
     .        facevary(ivar,1+nguard:nxb+nguard,
     .                   1+nguard*k2d:nyb+(nguard+1)*k2d,
     .                   1+nguard*k3d:nzb+nguard*k3d,lb)
#endif
              endif
            enddo

          elseif(lb.gt.lnblocks) then

          vtype = 3
          call mpi_set_message_limits(
     .                 dtype,ia,ib,ja,jb,ka,kb,vtype)

#ifdef NO_PERMANENT_GUARDCELLS
          ia = ia + nguard
          ib = ib + nguard
          ja = ja + nguard*k2d
          jb = jb + nguard*k2d
          ka = ka + nguard*k3d
          kb = kb + nguard*k3d
#endif

          do k = ka,kb
          do j = ja,jb
          do i = ia,ib
#ifdef OLD
            do ivar=1,nfacevar
              if(int_gcell_on_fc(2,ivar)) then
            facevary1(ivar,i,j,k,idest) =
     .                 temprecv_buf(index+ivar)
              endif
            enddo
            index = index+nfacevar
#else
            do ivar=1,ngcell_on_fc(2)
              ivar_next = gcell_on_fc_pointer(2,ivar)
              facevary1(ivar_next,i,j,k,idest) = 
     .                 temprecv_buf(index+ivar)
            enddo
            index = index+ngcell_on_fc(2)
#endif /* OLD */
          enddo
          enddo
          enddo

          endif

#else /* MPI_USED */

            call shmem_real_get(recvy(1,1,1,1),
     .                            facevary(1,1,1,1,lb),
     .                            nbndvar*len_blockfy,pe)

            do ivar=1,nfacevar
              if(int_gcell_on_fc(2,ivar)) then
              facevary1(ivar,:,:,:,idest) = 0.
#ifdef NO_PERMANENT_GUARDCELLS
              facevary1(ivar,1+nguard:nxb+nguard,
     .                   1+nguard*k2d:nyb+(nguard+1)*k2d,
     .                   1+nguard*k3d:nzb+nguard*k3d,idest) =
     .        recvy(ivar,1+nguard0:nxb+nguard0,
     .                   1+nguard0*k2d:nyb+(nguard0+1)*k2d,
     .                   1+nguard0*k3d:nzb+nguard0*k3d)
#else
            facevary1(ivar,:,:,:,idest) = recvy(ivar,:,:,:)
#endif
              endif
            enddo

#endif /* MPI_USED */

            endif  ! end (if ndim

            if(ndim.eq.3) then

#ifdef MPI_USED

            if(lb.le.lnblocks) then


            do ivar=1,nfacevar
              if(int_gcell_on_fc(3,ivar)) then
#ifdef NO_PERMANENT_GUARDCELLS
              facevarz1(ivar,1+nguard:nxb+nguard,
     .                    1+nguard*k2d:nyb+nguard*k2d,
     .                    1+nguard*k3d:nzb+(nguard+1)*k3d,idest) =
     .         facevarz(ivar,1+nguard0:nxb+nguard0,
     .                    1+nguard0*k2d:nyb+nguard0*k2d,
     .                    1+nguard0*k3d:nzb+(nguard0+1)*k3d,lb)
#else
              facevarz1(ivar,1+nguard:nxb+nguard,
     .                    1+nguard*k2d:nyb+nguard*k2d,
     .                    1+nguard*k3d:nzb+(nguard+1)*k3d,idest) =
     .         facevarz(ivar,1+nguard:nxb+nguard,
     .                    1+nguard*k2d:nyb+nguard*k2d,
     .                    1+nguard*k3d:nzb+(nguard+1)*k3d,lb)
#endif
              endif
            enddo

            elseif(lb.gt.lnblocks) then

          vtype = 4
          call mpi_set_message_limits(
     .                 dtype,ia,ib,ja,jb,ka,kb,vtype)

#ifdef NO_PERMANENT_GUARDCELLS
          ia = ia + nguard
          ib = ib + nguard
          ja = ja + nguard*k2d
          jb = jb + nguard*k2d
          ka = ka + nguard*k3d
          kb = kb + nguard*k3d
#endif

          do k = ka,kb
          do j = ja,jb
          do i = ia,ib
#ifdef OLD
            do ivar=1,nfacevar
              if(int_gcell_on_fc(3,ivar)) then
            facevarz1(ivar,i,j,k,idest) =
     .                 temprecv_buf(index+ivar)
              endif
            enddo
            index = index+nfacevar
#else
            do ivar=1,ngcell_on_fc(3)
              ivar_next = gcell_on_fc_pointer(3,ivar)
              facevarz1(ivar_next,i,j,k,idest) = 
     .                 temprecv_buf(index+ivar)
            enddo
            index = index+ngcell_on_fc(3)
#endif /* OLD */
          enddo
          enddo
          enddo

          endif

#else /* MPI_USED */

              call shmem_real_get(recvz(1,1,1,1),
     .                            facevarz(1,1,1,1,lb),
     .                            nbndvar*len_blockfz,pe)

            do ivar=1,nfacevar
              if(int_gcell_on_fc(3,ivar)) then
              facevarz1(ivar,:,:,:,idest) = 0.
#ifdef NO_PERMANENT_GUARDCELLS
              facevarz1(ivar,1+nguard:nxb+nguard,
     .                   1+nguard*k2d:nyb+nguard*k2d,
     .                   1+nguard*k3d:nzb+(nguard+1)*k3d,idest) =
     .         recvz(ivar,1+nguard0:nxb+nguard0,
     .                   1+nguard0*k2d:nyb+nguard0*k2d,
     .                   1+nguard0*k3d:nzb+(nguard0+1)*k3d)
#else
            facevarz1(ivar,:,:,:,idest) = recvz(ivar,:,:,:)
#endif
              endif
            enddo

#endif /* MPI_USED */

            endif                 ! end if (ndim

          endif                   ! end of lfc if test

          if (ndim > 1) then
          if(lec) then

#ifdef MPI_USED

            if(lb.le.lnblocks) then

            do ivar=1,nvaredge
              if(int_gcell_on_ec(1,ivar)) then
#ifdef NO_PERMANENT_GUARDCELLS
            unk_e_x1(ivar,1+nguard:nxb+nguard,
     .                 1+nguard*k2d:nyb+nguard*k2d+k2d,
     .                 1+nguard*k3d:nzb+nguard*k3d+k3d,idest) =
     .        unk_e_x(ivar,1+nguard0:nxb+nguard0,
     .                  1+nguard0*k2d:nyb+nguard0*k2d+k2d,
     .                  1+nguard0*k3d:nzb+nguard0*k3d+k3d,lb)
#else
            unk_e_x1(ivar,1+nguard:nxb+nguard,
     .                 1+nguard*k2d:nyb+nguard*k2d+k2d,
     .                 1+nguard*k3d:nzb+nguard*k3d+k3d,idest) =
     .        unk_e_x(ivar,1+nguard:nxb+nguard,
     .                  1+nguard*k2d:nyb+nguard*k2d+k2d,
     .                  1+nguard*k3d:nzb+nguard*k3d+k3d,lb)
#endif
              endif
            enddo

            elseif(lb.gt.lnblocks) then

! starting index if cell-centered data is also included in recv_buf
#ifdef OLD
         index = index0 + nvar*message_size_cc(dtype)
     .                  + nfacevar*message_size_fc(dtype)
#else
         index = index0 + ngcell_on_cc*message_size_cc(dtype)
     .                  + maxval(ngcell_on_fc(1:ndim))*
     .                             message_size_fc(dtype)
#endif /* OLD */

          vtype = 5
          call mpi_set_message_limits(
     .                 dtype,ia,ib,ja,jb,ka,kb,vtype)

#ifdef NO_PERMANENT_GUARDCELLS
          ia = ia + nguard
          ib = ib + nguard
          ja = ja + nguard*k2d
          jb = jb + nguard*k2d
          ka = ka + nguard*k3d
          kb = kb + nguard*k3d
#endif

          do k = ka,kb
          do j = ja,jb
          do i = ia,ib
#ifdef OLD
            do ivar=1,nvaredge
              if(int_gcell_on_ec(1,ivar)) then
            unk_e_x1(ivar,i,j,k,idest) =
     .                 temprecv_buf(index+ivar)
              endif
            enddo
            index = index+nvaredge
#else
            do ivar=1,ngcell_on_ec(1)
              ivar_next = gcell_on_ec_pointer(1,ivar)
              unk_e_x1(ivar_next,i,j,k,idest) =
     .                 temprecv_buf(index+ivar)
            enddo
            index = index+ngcell_on_ec(1)
#endif /* OLD */
          enddo
          enddo
          enddo

          endif

#else /* MPI_USED */

            call shmem_real_get(recvex(1,1,1,1),
     .                          unk_e_x(1,1,1,1,lb),
     .                          nbndvare*len_blockex,pe)

            do ivar=1,nvaredge
              if(int_gcell_on_ec(1,ivar)) then
            unk_e_x1(ivar,:,:,:,idest) = 0.
#ifdef NO_PERMANENT_GUARDCELLS
            unk_e_x1(ivar,1+nguard:nxb+nguard,
     .                 1+nguard*k2d:nyb+nguard*k2d+k2d,
     .                 1+nguard*k3d:nzb+nguard*k3d+k3d,idest) =
     .        recvex(ivar,1+nguard0:nxb+nguard0,
     .                 1+nguard0*k2d:nyb+nguard0*k2d+k2d,
     .                 1+nguard0*k3d:nzb+nguard0*k3d+k3d)
#else
            unk_e_x1(ivar,:,:,:,idest) = recvex(ivar,:,:,:)
#endif
              endif
            enddo

#endif /* MPI_USED */


#ifdef MPI_USED
            if(lb.le.lnblocks) then

            do ivar=1,nvaredge
              if(int_gcell_on_ec(2,ivar)) then
#ifdef NO_PERMANENT_GUARDCELLS
            unk_e_y1(ivar,1+nguard:nxb+nguard+1,
     .                 1+nguard*k2d:nyb+nguard*k2d,
     .                 1+nguard*k3d:nzb+(nguard+1)*k3d,idest) =
     .      unk_e_y(ivar,1+nguard0:nxb+nguard0+1,
     .                1+nguard0*k2d:nyb+nguard0*k2d,
     .                1+nguard0*k3d:nzb+(nguard0+1)*k3d,lb)
#else
            unk_e_y1(ivar,1+nguard:nxb+nguard+1,
     .                 1+nguard*k2d:nyb+nguard*k2d,
     .                 1+nguard*k3d:nzb+(nguard+1)*k3d,idest) =
     .      unk_e_y(ivar,1+nguard:nxb+nguard+1,
     .                1+nguard*k2d:nyb+nguard*k2d,
     .                1+nguard*k3d:nzb+(nguard+1)*k3d,lb)
#endif
              endif
            enddo

            elseif(lb.gt.lnblocks) then

          vtype = 6
          call mpi_set_message_limits(
     .                 dtype,ia,ib,ja,jb,ka,kb,vtype)

#ifdef NO_PERMANENT_GUARDCELLS
          ia = ia + nguard
          ib = ib + nguard
          ja = ja + nguard*k2d
          jb = jb + nguard*k2d
          ka = ka + nguard*k3d
          kb = kb + nguard*k3d
#endif

          do k = ka,kb
          do j = ja,jb
          do i = ia,ib
#ifdef OLD
            do ivar=1,nvaredge
              if(int_gcell_on_ec(2,ivar)) then
            unk_e_y1(ivar,i,j,k,idest) =
     .                 temprecv_buf(index+ivar)
              endif
            enddo
            index = index+nvaredge
#else
            do ivar=1,ngcell_on_ec(2)
              ivar_next = gcell_on_ec_pointer(2,ivar)
              unk_e_y1(ivar_next,i,j,k,idest) =
     .                 temprecv_buf(index+ivar)
            enddo
            index = index+ngcell_on_ec(2)
#endif /* OLD */
          enddo
          enddo
          enddo

          endif

#else /* MPI USED */

            call shmem_real_get(recvey(1,1,1,1),
     .                          unk_e_y(1,1,1,1,lb),
     .                          nbndvare*len_blockey,pe)

            do ivar=1,nvaredge
              if(int_gcell_on_ec(2,ivar)) then
            unk_e_y1(ivar,:,:,:,idest) = 0.
#ifdef NO_PERMANENT_GUARDCELLS
            unk_e_y1(ivar,1+nguard:nxb+nguard+1,
     .                 1+nguard*k2d:nyb+nguard*k2d,
     .                 1+nguard*k3d:nzb+(nguard+1)*k3d,idest) =
     .      recvey(ivar,1+nguard0:nxb+nguard0+1,
     .               1+nguard0*k2d:nyb+nguard0*k2d,
     .               1+nguard0*k3d:nzb+(nguard0+1)*k3d)
#else
            unk_e_y1(ivar,:,:,:,idest) = recvey(ivar,:,:,:)
#endif
              endif
            enddo

#endif /* MPI_USED */

            if (ndim == 3) then

#ifdef MPI_USED
            if(lb.le.lnblocks) then

            do ivar=1,nvaredge
              if(int_gcell_on_ec(3,ivar)) then
#ifdef NO_PERMANENT_GUARDCELLS
            unk_e_z1(ivar,1+nguard:nxb+nguard+1,
     .                 1+nguard*k2d:nyb+(nguard+1)*k2d,
     .                 1+nguard*k3d:nzb+nguard*k3d,idest) =
     .        unk_e_z(ivar,1+nguard0:nxb+nguard0+1,
     .                  1+nguard0*k2d:nyb+(nguard0+1)*k2d,
     .                  1+nguard0*k3d:nzb+nguard0*k3d,lb)
#else
            unk_e_z1(ivar,1+nguard:nxb+nguard+1,
     .                 1+nguard*k2d:nyb+(nguard+1)*k2d,
     .                 1+nguard*k3d:nzb+nguard*k3d,idest) =
     .        unk_e_z(ivar,1+nguard:nxb+nguard+1,
     .                  1+nguard*k2d:nyb+(nguard+1)*k2d,
     .                  1+nguard*k3d:nzb+nguard*k3d,lb)
#endif
              endif
            enddo

            elseif(lb.gt.lnblocks) then

          vtype = 7
          call mpi_set_message_limits(
     .                 dtype,ia,ib,ja,jb,ka,kb,vtype)

#ifdef NO_PERMANENT_GUARDCELLS
          ia = ia + nguard
          ib = ib + nguard
          ja = ja + nguard*k2d
          jb = jb + nguard*k2d
          ka = ka + nguard*k3d
          kb = kb + nguard*k3d
#endif

          do k = ka,kb
          do j = ja,jb
          do i = ia,ib
#ifdef OLD
            do ivar=1,nvaredge
              if(int_gcell_on_ec(3,ivar)) then
            unk_e_z1(ivar,i,j,k,idest) =
     .                 temprecv_buf(index+ivar)
              endif
            enddo
            index = index+nvaredge
#else
            do ivar=1,ngcell_on_ec(3)
              ivar_next = gcell_on_ec_pointer(3,ivar)
              unk_e_z1(ivar_next,i,j,k,idest) =
     .                 temprecv_buf(index+ivar)
            enddo
            index = index+ngcell_on_ec(3)
#endif /* OLD */
          enddo
          enddo
          enddo

          endif

#else /* MPI_USED */

            call shmem_real_get(recvez(1,1,1,1),
     .                          unk_e_z(1,1,1,1,lb),
     .                          nbndvare*len_blockez,pe)

            do ivar=1,nvaredge
              if(int_gcell_on_ec(3,ivar)) then
            unk_e_z1(ivar,:,:,:,idest) = 0.
#ifdef NO_PERMANENT_GUARDCELLS
            unk_e_z1(ivar,1+nguard:nxb+nguard+1,
     .                 1+nguard*k2d:nyb+(nguard+1)*k2d,
     .                 1+nguard*k3d:nzb+nguard*k3d,idest) =
     .        recvez(ivar,1+nguard0:nxb+nguard0+1,
     .                 1+nguard0*k2d:nyb+(nguard0+1)*k2d,
     .                 1+nguard0*k3d:nzb+nguard0*k3d)
#else
            unk_e_z1(ivar,:,:,:,idest) = recvez(ivar,:,:,:)
#endif
              endif
            enddo

#endif /* MPI_USED */

            end if                ! end if (ndim == 3

          endif                   ! end of lec if test
          end if

          if(lnc) then

#ifdef MPI_USED
            if(lb.le.lnblocks) then

            do ivar=1,nvarcorn
              if(int_gcell_on_nc(ivar)) then
#ifdef NO_PERMANENT_GUARDCELLS
            unk_n1(ivar,1+nguard:nxb+nguard+1,
     .               1+nguard*k2d:nyb+(nguard+1)*k2d,
     .               1+nguard*k3d:nzb+(nguard+1)*k3d,idest) =
     .        unk_n(ivar,1+nguard0:nxb+nguard0+1,
     .                1+nguard0*k2d:nyb+(nguard0+1)*k2d,
     .                1+nguard0*k3d:nzb+(nguard0+1)*k3d,lb)
#else
            unk_n1(ivar,1+nguard:nxb+nguard+1,
     .               1+nguard*k2d:nyb+(nguard+1)*k2d,
     .               1+nguard*k3d:nzb+(nguard+1)*k3d,idest) =
     .        unk_n(ivar,1+nguard:nxb+nguard+1,
     .                1+nguard*k2d:nyb+(nguard+1)*k2d,
     .                1+nguard*k3d:nzb+(nguard+1)*k3d,lb)
#endif
              endif
            enddo

            elseif(lb.gt.lnblocks) then

! starting index if cell-centered data is also included in recv_buf
#ifdef OLD
         index = index0 + nvar*message_size_cc(dtype)
     .                 + nfacevar*message_size_fc(dtype)
     .                 + nvaredge*message_size_ec(dtype)
#else
         index = index0 + ngcell_on_cc*message_size_cc(dtype)
     .                  + maxval(ngcell_on_fc(1:ndim))*
     .                             message_size_fc(dtype)
     .                  + maxval(ngcell_on_ec(1:ndim))*
     .                             message_size_ec(dtype)
#endif /* OLD */

          vtype = 8
          call mpi_set_message_limits(
     .                 dtype,ia,ib,ja,jb,ka,kb,vtype)

#ifdef NO_PERMANENT_GUARDCELLS
          ia = ia + nguard
          ib = ib + nguard
          ja = ja + nguard*k2d
          jb = jb + nguard*k2d
          ka = ka + nguard*k3d
          kb = kb + nguard*k3d
#endif

          do k = ka,kb
          do j = ja,jb
          do i = ia,ib
#ifdef OLD
            do ivar=1,nvarcorn
              if(int_gcell_on_nc(ivar)) then
            unk_n1(ivar,i,j,k,idest) =
     .                 temprecv_buf(index+ivar)
              endif
            enddo
            index = index+nvarcorn
#else
            do ivar=1,ngcell_on_nc
              ivar_next = gcell_on_nc_pointer(ivar)
              unk_n1(ivar_next,i,j,k,idest) =
     .                 temprecv_buf(index+ivar)
            enddo
            index = index+ngcell_on_nc
#endif /* OLD */
          enddo
          enddo
          enddo

          endif

#else /* MPI_USED */

            call shmem_real_get(recvn(1,1,1,1),
     .                          unk_n(1,1,1,1,lb),
     .                          len_blockn*nvarcorn,pe)

            do ivar=1,nvarcorn
              if(int_gcell_on_nc(ivar)) then
            unk_n1(ivar,:,:,:,idest) = 0.
#ifdef NO_PERMANENT_GUARDCELLS
            unk_n1(ivar,1+nguard:nxb+nguard+1,
     .               1+nguard*k2d:nyb+(nguard+1)*k2d,
     .               1+nguard*k3d:nzb+(nguard+1)*k3d,idest) =
     .        recvn(ivar,1+nguard0:nxb+nguard0+1,
     .                1+nguard0*k2d:nyb+(nguard0+1)*k2d,
     .                1+nguard0*k3d:nzb+(nguard0+1)*k3d)
#else
            unk_n1(ivar,:,:,:,idest) = recvn(ivar,:,:,:)
#endif
              endif
            enddo

#endif /* MPI_USED */

          endif                   ! end of lnc if test

          elseif(iopt.ge.2) then
            iopt0 = iopt-1

#ifdef MPI_USED

            if(lb.le.lnblocks) then

#ifdef NO_PERMANENT_GUARDCELLS
            work1(1+nguard_work:nxb+nguard_work,
     .            1+nguard_work*k2d:nyb+nguard_work*k2d,
     .            1+nguard_work*k3d:nzb+nguard_work*k3d,idest) =
     .  work(1+nguard_work0:nxb+nguard_work0,
     .       1+nguard_work0*k2d:nyb+nguard_work0*k2d,
     .       1+nguard_work0*k3d:nzb+nguard_work0*k3d,lb,iopt0)
#else
            work1(1+nguard_work:nxb+nguard_work,
     .            1+nguard_work*k2d:nyb+nguard_work*k2d,
     .            1+nguard_work*k3d:nzb+nguard_work*k3d,idest) =
     .  work(1+nguard_work:nxb+nguard_work,
     .       1+nguard_work*k2d:nyb+nguard_work*k2d,
     .       1+nguard_work*k3d:nzb+nguard_work*k3d,lb,iopt0)
#endif

            elseif(lb.gt.lnblocks) then

          vtype = 0
          index = index0
          call mpi_set_message_limits(
     .                 dtype,ia,ib,ja,jb,ka,kb,vtype)

#ifdef NO_PERMANENT_GUARDCELLS
          ia = ia + nguard_work
          ib = ib + nguard_work
          ja = ja + nguard_work*k2d
          jb = jb + nguard_work*k2d
          ka = ka + nguard_work*k3d
          kb = kb + nguard_work*k3d
#endif

          do k = ka,kb
          do j = ja,jb
          do i = ia,ib
            work1(i,j,k,idest) = 
     .              temprecv_buf(index+1)
            index = index+1
          enddo
          enddo
          enddo

          endif

#else /* MPI_USED */

            call shmem_real_get(recvw(1,1,1),
     .                          work(1,1,1,lb,iopt0),
     .                          len_wblock,pe)
            work1(:,:,:,idest) = 0.
#ifdef NO_PERMANENT_GUARDCELLS
            work1(1+nguard_work:nxb+nguard_work,
     .            1+nguard_work*k2d:nyb+nguard_work*k2d,
     .            1+nguard_work*k3d:nzb+nguard_work*k3d,idest) =
     .         recvw(1+nguard_work0:nxb+nguard_work0,
     .               1+nguard_work0*k2d:nyb+nguard_work0*k2d,
     .               1+nguard_work0*k3d:nzb+nguard_work0*k3d)
#else
            work1(:,:,:,idest) = recvw(:,:,:)
#endif
#endif /* MPI_USED */

          endif                 ! end of iopt if test



      return
      end subroutine amr_perm_to_1blk
