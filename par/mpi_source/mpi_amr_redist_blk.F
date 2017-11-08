!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2002
! By: Peter J. MacNeice, Drexel University.
!     Kevin M. Olson, Univ. of MD Baltimore Campus.
! 
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
! 
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
! 
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
! USA
!----------------------------------------------------------------------
#include "paramesh_preprocessor.fh"

      subroutine amr_redist_blk(new_loc,nprocs,mype,lnblocks_old)



      use paramesh_dimensions
      use physicaldata
      use tree

      use paramesh_interfaces, only : fill_old_loc

      implicit none

      include 'mpif.h'


      integer, intent(inout) :: new_loc(:,:)
      integer, intent(in)    :: nprocs,mype,lnblocks_old

      integer :: lb,ierr,errorcode
      logical :: free(maxblocks), moved(maxblocks), sent(maxblocks)
      logical :: repeat, repeatt
      integer :: old_loc(2,maxblocks_tr)
      integer :: nsend, nrecv
      integer :: reqr(maxblocks_tr)
      integer :: reqs(maxblocks_tr)
      integer :: statr(MPI_STATUS_SIZE,maxblocks_tr)
      integer :: stats(MPI_STATUS_SIZE,maxblocks_tr)
      integer :: nmoved, nit
      integer :: test(maxblocks), point_to(maxblocks)
      integer :: nm, nm2, nm2_old

      integer :: ireduce_datain(1),ireduce_dataout(1)
      logical :: lreduce_datain(1),lreduce_dataout(1)
      integer, save :: myblockint, block_int2d, block_int3d
      logical, save :: first = .true.
      integer,save ::  istart,jstart,kstart,iend,jend,kend

      if (first) then
      first = .false.

      istart = nguard*npgs+1
      jstart = nguard*k2d*npgs+1
      kstart = nguard*k3d*npgs+1
      iend   = nguard*npgs+nxb
      jend   = nguard*k2d*npgs+nyb
      kend   = nguard*k3d*npgs+nzb

      ! DEFINE BLOCK INTERIOR
      call MPI_TYPE_VECTOR (nyb, 
     &                      nvar*nxb, 
     &                      nvar*iu_bnd, 
     &                      MPI_DOUBLE_PRECISION, 
     &                      block_int2d, 
     &                      ierr)
      myblockint = block_int2d
      if (ndim.eq.3) then
         call MPI_TYPE_HVECTOR (nzb, 
     &                          1, 
     &                          nvar*iu_bnd*ju_bnd*8, 
     &                          block_int2d, 
     &                          block_int3d, 
     &                          ierr)
         myblockint = block_int3d
      end if
      if (ndim.eq.3) call MPI_TYPE_COMMIT(block_int2d,ierr)
      call MPI_TYPE_COMMIT(myblockint,ierr)
      end if


! 1) compute old_loc
      call fill_old_loc (new_loc,old_loc,nprocs,mype)

      nrecv = 0
      nsend = 0


!--------------
! treat unk

      if(nvar.gt.0) then

! Post all receives for unk
         do lb = 1,new_lnblocks
            if (.not.newchild(lb)) then
               if (old_loc(2,lb).ne.mype) then
                  nrecv = nrecv + 1
                  call MPI_IRECV (unk(1,istart,jstart,kstart,lb),
     &                            1,
     &                            myblockint,
     &                            old_loc(2,lb),
     &                            lb,
     &                            MPI_COMM_WORLD,
     &                            reqr(nrecv),
     &                            ierr)

               end if
            end if
         end do

      end if

!--------------------
! Treat Facevariables
!--------------------

      if (nfacevar.gt.0) then

!---------------
! Treat facevarx
!---------------

         do lb = 1,new_lnblocks
            if (.not.newchild(lb)) then
               if (old_loc(2,lb).ne.mype) then
                  nrecv = nrecv + 1
                  call MPI_IRECV (facevarx(1,1,1,1,lb),
     &                            nbndvar*len_blockfx,
     &                            MPI_REAL,
     &                            old_loc(2,lb),
     &                            lb+2*maxblocks,
     &                            MPI_COMM_WORLD,
     &                            reqr(nrecv),
     &                            ierr)
               end if
            end if
         end do

!---------------
! Treat facevary
!---------------

         if (ndim >= 2) then
         do lb = 1,new_lnblocks
            if (.not.newchild(lb)) then
               if (old_loc(2,lb).ne.mype) then
                  nrecv = nrecv + 1
                  call MPI_IRECV (facevary(1,1,1,1,lb),
     &                            nbndvar*len_blockfy,
     &                            MPI_REAL,
     &                            old_loc(2,lb),
     &                            lb+3*maxblocks,
     &                            MPI_COMM_WORLD,
     &                            reqr(nrecv),
     &                            ierr)
               end if
            end if
         end do
         end if

!---------------
! Treat Facevarz
!---------------

         if (ndim == 3) then
         do lb = 1,new_lnblocks
            if (.not.newchild(lb)) then
               if (old_loc(2,lb).ne.mype) then
                  nrecv = nrecv + 1
                  call MPI_IRECV (facevarz(1,1,1,1,lb),
     &                            nbndvar*len_blockfz,
     &                            MPI_REAL,
     &                            old_loc(2,lb),
     &                            lb+4*maxblocks,
     &                            MPI_COMM_WORLD,
     &                            reqr(nrecv),
     &                            ierr)
               end if
            end if
         end do
         end if

      end if

!--------------------
! Treat Edge variables
!--------------------

      if (nvaredge.gt.0) then

!---------------
! Treat unk_e_x
!---------------

         do lb = 1,new_lnblocks
            if (.not.newchild(lb)) then
               if (old_loc(2,lb).ne.mype) then
                  nrecv = nrecv + 1
                  call MPI_IRECV (unk_e_x(1,1,1,1,lb),
     &                            nbndvare*len_blockex,
     &                            MPI_REAL,
     &                            old_loc(2,lb),
     &                            lb+5*maxblocks,
     &                            MPI_COMM_WORLD,
     &                            reqr(nrecv),
     &                            ierr)
               end if
            end if
         end do

!---------------
! Treat unk_e_y
!---------------

         if (ndim >= 2) then
         do lb = 1,new_lnblocks
            if (.not.newchild(lb)) then
               if (old_loc(2,lb).ne.mype) then
                  nrecv = nrecv + 1
                  call MPI_IRECV (unk_e_y(1,1,1,1,lb),
     &                            nbndvare*len_blockey,
     &                            MPI_REAL,
     &                            old_loc(2,lb),
     &                            lb+6*maxblocks,
     &                            MPI_COMM_WORLD,
     &                            reqr(nrecv),
     &                            ierr)
               end if
            end if
         end do
         end if
!---------------
! Treat unk_e_z
!---------------

         if (ndim == 3) then
         do lb = 1,new_lnblocks
            if (.not.newchild(lb)) then
               if (old_loc(2,lb).ne.mype) then
                  nrecv = nrecv + 1
                  call MPI_IRECV (unk_e_z(1,1,1,1,lb),
     &                            nbndvare*len_blockez,
     &                            MPI_REAL,
     &                            old_loc(2,lb),
     &                            lb+7*maxblocks,
     &                            MPI_COMM_WORLD,
     &                            reqr(nrecv),
     &                            ierr)
               end if
            end if
         end do
         end if

      endif

!--------------
! treat unk_n

      if(nvarcorn.gt.0) then

! Post all receives for unk_n
         do lb = 1,new_lnblocks
            if (.not.newchild(lb)) then
               if (old_loc(2,lb).ne.mype) then
                  nrecv = nrecv + 1
                  call MPI_IRECV (unk_n(1,1,1,1,lb),
     &                            nbndvarc*len_blockn,
     &                            MPI_REAL,
     &                            old_loc(2,lb),
     &                            lb+8*maxblocks,
     &                            MPI_COMM_WORLD,
     &                            reqr(nrecv),
     &                            ierr)
               end if
            end if
         end do

      end if

!--------------

      
      moved(:) = .false.
      moved(lnblocks_old+1:maxblocks) = .true.
      free(:) = .false.
      free(lnblocks_old+1:maxblocks) = .true.
      sent(:) = .false.
      repeat = .TRUE.
      nmoved = 0 
      test(:) = 0
      point_to(:) = 0
      
      nit = 0
      nm2 = 0
      nm2_old = 1
      do while (repeat.and.nit<=100) 
         
         do lb = 1,max(lnblocks_old,new_lnblocks)
            call send_block_data (lb, new_loc, old_loc, free, 
     &                            moved, sent,
     &                            lnblocks_old, mype, nmoved,
     &                            test, point_to,
     &                            reqs, nsend, myblockint)
         end do
         repeat = any(.not.moved(:))
         lreduce_datain(1) = repeat
         call mpi_logical_allreduce(
     &           lreduce_datain(1),lreduce_dataout(1),
     &           1,MPI_LOGICAL,
     &           MPI_LOR,MPI_COMM_WORLD,ierr)
         repeatt = lreduce_dataout(1)
         repeat = repeatt
         
         nm2_old = nm2
         nm = count(.not.moved(:))
         ireduce_datain(1) = nm
         call mpi_int_allreduce(
     &        ireduce_datain(1),ireduce_dataout(1),
     &        1,MPI_INTEGER,
     &        MPI_SUM,MPI_COMM_WORLD,ierr)
         nm2 = ireduce_dataout(1)
         if (mype.eq.0) then
            print *,' interation, no. not moved = ',nit,nm2
         end if
         
         nit = nit + 1
         
      end do
      
      if (nm2_old.eq.nm2.and.nm2.ne.0.and.nit>=100) then
         if (mype.eq.0) then
          print *,' ERROR: could not move all blocks in amr_redist_blk'
          print *,' Try increasing maxblocks or use more processors'
          print *,' nm2_old, nm2 = ',nm2_old,nm2
          print *,' ABORTING !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
         end if
         call MPI_ABORT(MPI_COMM_WORLD,errorcode,ierr)
      end if
      
      if (nrecv.gt.0) then
         call MPI_WAITALL(nrecv,reqr,statr,ierr)
      end if

! NOTE: Bug fix by Paul Ricker and Marcus Gross (5/2003).  Added Waitall 
! to isends so that SGI MPI buffers do not overflow.

      if (nsend.gt.0) then
         call MPI_WAITALL(nsend,reqs,stats,ierr)
      end if
      
      return
      end subroutine amr_redist_blk
      
!----------------------------------------------------------------------------
      
      subroutine send_block_data (lb, new_loc, old_loc, free,
     &                            moved, sent,
     &                            lnblocks_old, mype, nmoved,
     &                            test, point_to,
     &                            reqs, nsend, myblockint)
         
      use paramesh_dimensions
      use physicaldata
      use tree

      implicit none

      include 'mpif.h'


      integer :: new_loc(2,maxblocks_tr), old_loc(2,maxblocks_tr)
      logical :: free(maxblocks), moved(maxblocks), sent(maxblocks)
      integer :: lb, lnblocks_old, mype
      logical :: success
      integer :: status(MPI_STATUS_SIZE)
      integer :: reqs(maxblocks_tr), nsend
      integer :: ierr, nmoved, lb2
      integer :: point_to(maxblocks),test(maxblocks)
      integer :: myblockint
      logical, save :: first = .true.
      integer,save ::  istart,jstart,kstart,iend,jend,kend

      if (first) then

      first = .false.
      istart = nguard*npgs+1
      jstart = nguard*k2d*npgs+1
      kstart = nguard*k3d*npgs+1
      iend   = nguard*npgs+nxb
      jend   = nguard*k2d*npgs+nyb
      kend   = nguard*k3d*npgs+nzb

      end if

      if (new_loc(1,lb).eq.lb.and.new_loc(2,lb).eq.mype) then
         if (.not.moved(lb)) moved(lb) = .true.
         return
      end if

      if (lb.le.max(lnblocks_old,new_lnblocks)) then

         if (lb.le.lnblocks_old) then
           if (new_loc(2,lb).ne.mype) then
            success = .false.
            call MPI_IPROBE (new_loc(2,lb),
     &                       maxblocks+new_loc(1,lb),
     &                       MPI_COMM_WORLD,
     &                       success,
     &                       status,
     &                       ierr)
            if (.not.moved(lb).and.success) then
               call MPI_RECV (success,
     &                        1,
     &                        MPI_LOGICAL,
     &                        new_loc(2,lb),
     &                        maxblocks+new_loc(1,lb),
     &                        MPI_COMM_WORLD,
     &                        status,
     &                        ierr)
               if (free(lb)) then
                if (nvar.gt.0) then
                call MPI_SSEND (
     &   unk(1,istart,jstart,kstart,point_to(lb)),
     &                          1,
     &                          myblockint,
     &                          new_loc(2,lb),
     &                          new_loc(1,lb),
     &                          MPI_COMM_WORLD,
     &                          ierr)
                endif
! send facevariables
                if (nfacevar.gt.0) then
                call MPI_SSEND (facevarx(1,1,1,1,point_to(lb)),
     &                          nbndvar*len_blockfx,
     &                          MPI_REAL,
     &                          new_loc(2,lb),
     &                          new_loc(1,lb)+2*maxblocks,
     &                          MPI_COMM_WORLD,
     &                          ierr)
                if (ndim >= 2) then
                call MPI_SSEND (facevary(1,1,1,1,point_to(lb)),
     &                          nbndvar*len_blockfy,
     &                          MPI_REAL,
     &                          new_loc(2,lb),
     &                          new_loc(1,lb)+3*maxblocks,
     &                          MPI_COMM_WORLD,
     &                          ierr)
                end if
                if (ndim == 3) then
                call MPI_SSEND (facevarz(1,1,1,1,point_to(lb)),
     &                          nbndvar*len_blockfz,
     &                          MPI_REAL,
     &                          new_loc(2,lb),
     &                          new_loc(1,lb)+4*maxblocks,
     &                          MPI_COMM_WORLD,
     &                          ierr)
                end if

                end if

! send edge variables
                if (nvaredge.gt.0) then
                call MPI_SSEND (unk_e_x(1,1,1,1,point_to(lb)),
     &                          nbndvare*len_blockex,
     &                          MPI_REAL,
     &                          new_loc(2,lb),
     &                          new_loc(1,lb)+5*maxblocks,
     &                          MPI_COMM_WORLD,
     &                          ierr)
                if (ndim >= 2) then
                call MPI_SSEND (unk_e_y(1,1,1,1,point_to(lb)),
     &                          nbndvare*len_blockey,
     &                          MPI_REAL,
     &                          new_loc(2,lb),
     &                          new_loc(1,lb)+6*maxblocks,
     &                          MPI_COMM_WORLD,
     &                          ierr)
                end if
                if (ndim == 3) then
                call MPI_SSEND (unk_e_z(1,1,1,1,point_to(lb)),
     &                          nbndvare*len_blockez,
     &                          MPI_REAL,
     &                          new_loc(2,lb),
     &                          new_loc(1,lb)+7*maxblocks,
     &                          MPI_COMM_WORLD,
     &                          ierr)
                end if
                end if

! send corner variables
                if (nvarcorn.gt.0) then
                call MPI_SSEND (unk_n(1,1,1,1,point_to(lb)),
     &                          nbndvarc*len_blockn,
     &                          MPI_REAL,
     &                          new_loc(2,lb),
     &                          new_loc(1,lb)+8*maxblocks,
     &                          MPI_COMM_WORLD,
     &                          ierr)
                end if

                test(point_to(lb)) = -1
               else
                if (nvar.gt.0) then
                call MPI_SSEND (unk(1,istart,jstart,kstart,lb),
     &                          1,
     &                          myblockint,
     &                          new_loc(2,lb),
     &                          new_loc(1,lb),
     &                          MPI_COMM_WORLD,
     &                          ierr)
                end if
! send facevariables
                if (nfacevar.gt.0) then
                call MPI_SSEND (facevarx(1,1,1,1,lb),
     &                          nbndvar*len_blockfx,
     &                          MPI_REAL,
     &                          new_loc(2,lb),
     &                          new_loc(1,lb)+2*maxblocks,
     &                          MPI_COMM_WORLD,
     &                          ierr)
                if (ndim >= 2) then
                call MPI_SSEND (facevary(1,1,1,1,lb),
     &                          nbndvar*len_blockfy,
     &                          MPI_REAL,
     &                          new_loc(2,lb),
     &                          new_loc(1,lb)+3*maxblocks,
     &                          MPI_COMM_WORLD,
     &                          ierr)
                end if
                if (ndim == 3) then
                call MPI_SSEND (facevarz(1,1,1,1,lb),
     &                          nbndvar*len_blockfz,
     &                          MPI_REAL,
     &                          new_loc(2,lb),
     &                          new_loc(1,lb)+4*maxblocks,
     &                          MPI_COMM_WORLD,
     &                          ierr)
                end if
                end if

! send edge variables
                if (nvaredge.gt.0) then
                call MPI_SSEND (unk_e_x(1,1,1,1,lb),
     &                          nbndvare*len_blockex,
     &                          MPI_REAL,
     &                          new_loc(2,lb),
     &                          new_loc(1,lb)+5*maxblocks,
     &                          MPI_COMM_WORLD,
     &                          ierr)
                if (ndim >= 2) then
                call MPI_SSEND (unk_e_y(1,1,1,1,lb),
     &                          nbndvare*len_blockey,
     &                          MPI_REAL,
     &                          new_loc(2,lb),
     &                          new_loc(1,lb)+6*maxblocks,
     &                          MPI_COMM_WORLD,
     &                          ierr)
                end if
                if (ndim == 3) then
                call MPI_SSEND (unk_e_z(1,1,1,1,lb),
     &                          nbndvare*len_blockez,
     &                          MPI_REAL,
     &                          new_loc(2,lb),
     &                          new_loc(1,lb)+7*maxblocks,
     &                          MPI_COMM_WORLD,
     &                          ierr)
                end if
                end if

! send corner variables
                if (nvarcorn.gt.0) then
                call MPI_SSEND (unk_n(1,1,1,1,lb),
     &                          nbndvarc*len_blockn,
     &                          MPI_REAL,
     &                          new_loc(2,lb),
     &                          new_loc(1,lb)+8*maxblocks,
     &                          MPI_COMM_WORLD,
     &                          ierr)
                end if

                free(lb) = .true.
               end if
               moved(lb) = .true.
            end if
           else
            if (.not.moved(lb).and.free(new_loc(1,lb))) then
             if (free(lb)) then

               if (nvar.gt.0) then
         unk(:,istart:iend,jstart:jend,kstart:kend,new_loc(1,lb)) = 
     &    unk(:,istart:iend,jstart:jend,kstart:kend,point_to(lb))
               end if

! move facevars
               if (nfacevar.gt.0) then
                  facevarx(:,:,:,:,new_loc(1,lb)) = 
     &                 facevarx(:,:,:,:,point_to(lb))
                  if (ndim >= 2) then
                  facevary(:,:,:,:,new_loc(1,lb)) = 
     &                 facevary(:,:,:,:,point_to(lb))
                  end if
                  if (ndim == 3) then
                  facevarz(:,:,:,:,new_loc(1,lb)) = 
     &                 facevarz(:,:,:,:,point_to(lb))
                  end if
               end if

! move edgevars
               if (nvaredge.gt.0) then
                  unk_e_x(:,:,:,:,new_loc(1,lb)) =
     &                 unk_e_x(:,:,:,:,point_to(lb))
                  if (ndim >= 2) then
                  unk_e_y(:,:,:,:,new_loc(1,lb)) =
     &                 unk_e_y(:,:,:,:,point_to(lb))
                  end if
                  if (ndim == 3) then
                  unk_e_z(:,:,:,:,new_loc(1,lb)) =
     &                 unk_e_z(:,:,:,:,point_to(lb))
                  end if
               end if

               if (nvarcorn.gt.0) then
                  unk_n(:,:,:,:,new_loc(1,lb)) = 
     &                 unk_n(:,:,:,:,point_to(lb))
               end if

               test(point_to(lb)) = -1
             else
               if (nvar.gt.0) then
        unk(:,istart:iend,jstart:jend,kstart:kend,new_loc(1,lb)) = 
     &   unk(:,istart:iend,jstart:jend,kstart:kend,lb)
               end if
! move facevars
               if (nfacevar.gt.0) then
                  facevarx(:,:,:,:,new_loc(1,lb)) = 
     &                 facevarx(:,:,:,:,lb)
                  if (ndim >= 2) then
                  facevary(:,:,:,:,new_loc(1,lb)) = 
     &                 facevary(:,:,:,:,lb)
                  end if
                  if (ndim == 3) then
                  facevarz(:,:,:,:,new_loc(1,lb)) = 
     &                 facevarz(:,:,:,:,lb)
                  end if
               end if
! move edgevars
               if (nvaredge.gt.0) then
                  unk_e_x(:,:,:,:,new_loc(1,lb)) =
     &                 unk_e_x(:,:,:,:,lb)
                  if (ndim >= 2) then
                  unk_e_y(:,:,:,:,new_loc(1,lb)) =
     &                 unk_e_y(:,:,:,:,lb)
                  end if
                  if (ndim == 3) then
                  unk_e_z(:,:,:,:,new_loc(1,lb)) =
     &                 unk_e_z(:,:,:,:,lb)
                  end if
               end if
               if (nvarcorn.gt.0) then
               unk_n(:,:,:,:,new_loc(1,lb)) = unk_n(:,:,:,:,lb)
               end if

               free(lb) = .true.
             end if
             moved(lb) = .true.
            end if
           end if
         end if

         if (lb.le.new_lnblocks) then
            if (free(lb).and..not.sent(lb)) then
               sent(lb) = .true.
               if (.not.newchild(lb)) then
                  if (old_loc(2,lb).ne.mype) then
                     nsend = nsend + 1
                     call MPI_ISEND (free(lb),
     &                               1,
     &                               MPI_LOGICAL,
     &                               old_loc(2,lb),
     &                               maxblocks+lb,
     &                               MPI_COMM_WORLD,
     &                               reqs(nsend),
     &                               ierr)
                  end if
               end if
            end if
         end if

         if (lb.le.lnblocks_old.and..not.free(lb)) then
            nmoved = nmoved + 1
            point_to(lb) = max(lnblocks_old,new_lnblocks)+nmoved
            if (point_to(lb).gt.maxblocks) then
               do lb2 = max(lnblocks_old,new_lnblocks)+1,maxblocks
                  if (test(lb2).eq.-1) then
                     point_to(lb) = lb2
                     go to 22
                  end if
               end do
            end if
 22         if (point_to(lb).le.maxblocks) then
               test(point_to(lb)) = 1
               if (nvar.gt.0) then
               unk(:,istart:iend,jstart:jend,kstart:kend,point_to(lb)) = 
     &          unk(:,istart:iend,jstart:jend,kstart:kend,lb)
               end if
! move facevars
               if (nfacevar.gt.0) then
                  facevarx(:,:,:,:,point_to(lb)) = 
     &                 facevarx(:,:,:,:,lb)
                  if (ndim >= 2) then
                  facevary(:,:,:,:,point_to(lb)) = 
     &                 facevary(:,:,:,:,lb)
                  end if
                  if (ndim == 3) then
                  facevarz(:,:,:,:,point_to(lb)) = 
     &                 facevarz(:,:,:,:,lb)
                  end if
               end if
! move edgevars
               if (nvaredge.gt.0) then
                  unk_e_x(:,:,:,:,point_to(lb)) =
     &                 unk_e_x(:,:,:,:,lb)
                  if (ndim >= 2) then
                  unk_e_y(:,:,:,:,point_to(lb)) =
     &                 unk_e_y(:,:,:,:,lb)
                  end if
                  if (ndim == 3) then
                  unk_e_z(:,:,:,:,point_to(lb)) =
     &                 unk_e_z(:,:,:,:,lb)
                  end if
               end if
               if (nvarcorn.gt.0) then
               unk_n(:,:,:,:,point_to(lb)) = unk_n(:,:,:,:,lb)
               end if
               free(lb) = .TRUE.
            end if
         end if

         return

      else

         return

      end if

      return
      end subroutine send_block_data
