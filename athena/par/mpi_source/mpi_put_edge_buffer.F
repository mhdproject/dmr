!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"


      subroutine mpi_put_edge_buffer(lb,offset,
     .                          buffer_size,R_buffer)

!------------------------------------------------------------------------
!
! This subroutine unpacks the buffer R_buffer which has been received on 
! this processor into the bedge_faceX_X arrays for the local block lb.
!
!
! Written :     Maharaj Bhat & Michael Gehmeyr          March 2000
!------------------------------------------------------------------------
!
! Arguments:
!      lb             guard block id to be unpacked from receiving
!      offset         offset for buffer index
!      buffer_size    size of receive buffer
!      R_buffer       receive buffer
!
!------------------------------------------------------------------------
      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace
      use mpi_morton

      use paramesh_mpi_interfaces, only : mpi_set_message_limits

      implicit none

      include 'mpif.h'

      integer, intent(in)    :: lb,buffer_size
      integer, intent(inout) :: offset
      real,    intent(inout) :: R_buffer(buffer_size)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local variables
      integer, parameter :: nguard0 = nguard*npgs
      integer, parameter :: nguard_work0 = nguard_work*npgs
      integer :: index
      integer :: ia,ib,ja,jb,ka,kb
      integer :: ia0,ib0,ja0,jb0,ka0,kb0
      integer :: dtype,vtype
      integer :: i, j, k, n, ii, jj, iii

      integer :: mype
      integer,external :: shmem_my_pe

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      mype = shmem_my_pe()

! take over incoming index offset for block lb to be sent to remote pe

      index = offset          

      dtype = anint(R_buffer(index+2))
      index = index+3

      vtype = 8
      call mpi_set_message_limits(dtype,
     .                            ia0,ib0,ja0,jb0,ka0,kb0,vtype)


! unpack the bedge_facex_y and bedge_facex_z arrays for block lb
      if(dtype.eq.13.or.dtype.eq.15.or.dtype.eq.14) then

      ia = ia0
      ib = ib0
      ja = ja0
      jb = jb0
      ka = ka0
      kb = kb0

      if(dtype.eq.13) then
        ia = 1
        ib = 1
      elseif(dtype.eq.15) then
        ia = 2
        ib = 2
      elseif(dtype.eq.14) then
        ia = 1
        ib = 2
      endif

      do k = ka,kb
      do j = ja,jb
      do i = ia,ib
        do n=1,nedges
          bedge_facex_y(n,i,j,k,lb) = R_buffer(index)
          index  = index + 1
        enddo
      enddo
      enddo
      enddo
      if(dtype.eq.13) then
        bedge_facez_y(1:nedges,1+nguard0,ja:jb,1,lb) =
     .      bedge_facex_y(1:nedges,1,ja:jb,1+nguard0*k3d,lb)
        bedge_facez_y(1:nedges,1+nguard0,ja:jb,1+k3d,lb) =
     .      bedge_facex_y(1:nedges,1,ja:jb,nzb+(1+nguard0)*k3d,lb)
      elseif(dtype.eq.15) then
        bedge_facez_y(1:nedges,nxb+1+nguard0,ja:jb,1,lb) =
     .      bedge_facex_y(1:nedges,2,ja:jb,1+nguard0*k3d,lb)
        bedge_facez_y(1:nedges,nxb+1+nguard0,ja:jb,1+k3d,lb) =
     .      bedge_facex_y(1:nedges,2,ja:jb,nzb+(1+nguard0)*k3d,lb)
      endif

      if(ndim.eq.3.or.l2p5d.eq.1) then
      do k = ka , kb
      do j = ja , jb
      do i = ia , ib
        do n=1,nedges
          bedge_facex_z(n,i,j,k,lb) = R_buffer(index)
          index  = index + 1
        enddo
      enddo
      enddo
      enddo
      if(dtype.eq.13) then
        bedge_facey_z(1:nedges,1+nguard0,1,ka:kb,lb) =
     .      bedge_facex_z(1:nedges,1,1+nguard0*k2d,ka:kb,lb)
        bedge_facey_z(1:nedges,1+nguard0,2,ka:kb,lb) =
     .      bedge_facex_z(1:nedges,1,nyb+(1+nguard0)*k2d,ka:kb,lb)
      elseif(dtype.eq.15) then
        bedge_facey_z(1:nedges,nxb+1+nguard0,1,ka:kb,lb) =
     .      bedge_facex_z(1:nedges,2,1+nguard0*k2d,ka:kb,lb)
        bedge_facey_z(1:nedges,nxb+1+nguard0,2,ka:kb,lb) =
     .      bedge_facex_z(1:nedges,2,nyb+(1+nguard0)*k2d,ka:kb,lb)
      endif
      endif

      endif


! pack the bedge_facey_x and bedge_facey_z arrays for block lb

      if(ndim.ge.2) then
      if(dtype.eq.11.or.dtype.eq.17.or.dtype.eq.14) then

      ia = ia0
      ib = ib0
      ja = ja0
      jb = jb0
      ka = ka0
      kb = kb0

      if(dtype.eq.11) then
        ja = 1
        jb = 1
      elseif(dtype.eq.17) then
        ja = 2
        jb = 2
      elseif(dtype.eq.14) then
        ja = 1
        jb = 2
      endif

      if(mype.eq.0.and.lb.eq.10)write(*,*) 'updating bedge blk 10',dtype

      do k = ka,kb
      do j = ja,jb
      do i = ia,ib
        do n=1,nedges
          bedge_facey_x(n,i,j,k,lb) = R_buffer(index)
          index  = index + 1
        enddo
      enddo
      enddo
      enddo
      if(dtype.eq.11) then
        bedge_facez_x(1:nedges,ia:ib,1+nguard0*k2d,1,lb) =
     .      bedge_facey_x(1:nedges,ia:ib,1,1+nguard0*k3d,lb)
        bedge_facez_x(1:nedges,ia:ib,1+nguard0*k2d,1+k3d,lb) =
     .      bedge_facey_x(1:nedges,ia:ib,1,nzb+(1+nguard0)*k3d,lb)
      elseif(dtype.eq.17) then
        bedge_facez_x(1:nedges,ia:ib,nyb+(1+nguard0)*k2d,1,lb) =
     .      bedge_facey_x(1:nedges,ia:ib,2,1+nguard0*k3d,lb)
        bedge_facez_x(1:nedges,ia:ib,nyb+(1+nguard0)*k2d,1+k3d,lb) =
     .      bedge_facey_x(1:nedges,ia:ib,2,nzb+(1+nguard0)*k3d,lb)
      endif

      if(ndim.eq.3.or.l2p5d.eq.1) then
      do k = ka , kb
      do j = ja , jb
      do i = ia , ib
        do n=1,nedges
          bedge_facey_z(n,i,j,k,lb) = R_buffer(index)
          index  = index + 1
        enddo
      enddo
      enddo
      enddo
      if(dtype.eq.11) then
        bedge_facex_z(1:nedges,1,1+nguard0*k2d,ka:kb,lb) =
     .      bedge_facey_z(1:nedges,1+nguard0,1,ka:kb,lb)
        bedge_facex_z(1:nedges,2,1+nguard0*k2d,ka:kb,lb) =
     .      bedge_facey_z(1:nedges,nxb+1+nguard0,1,ka:kb,lb)
      elseif(dtype.eq.17) then
        bedge_facex_z(1:nedges,1,nyb+(1+nguard0)*k2d,ka:kb,lb) =
     .      bedge_facey_z(1:nedges,1+nguard0,2,ka:kb,lb)
        bedge_facex_z(1:nedges,2,nyb+(1+nguard0)*k2d,ka:kb,lb) =
     .      bedge_facey_z(1:nedges,nxb+1+nguard0,2,ka:kb,lb)
      endif
      endif


      endif
      endif

! pack the bedge_facez_x array for block lb

!pmn      if(ndim.eq.3.or.l2p5d.eq.1) then
      if(ndim.eq.3) then
      if(dtype.eq.5.or.dtype.eq.23.or.dtype.eq.14) then

      ia = ia0
      ib = ib0
      ja = ja0
      jb = jb0
      ka = ka0
      kb = kb0

      if(dtype.eq.5) then
        ka = 1
        kb = 1
      elseif(dtype.eq.23) then
        ka = 2
        kb = 2
      elseif(dtype.eq.14) then
        ka = 1
        kb = 1+k3d
      endif


      do k = ka,kb
      do j = ja,jb
      do i = ia,ib
        do n=1,nedges
          bedge_facez_x(n,i,j,k,lb) = R_buffer(index)
          index  = index + 1
        enddo
      enddo
      enddo
      enddo
      if(dtype.eq.5) then
      if(mype.eq.0.and.lb.eq.10)write(*,*) 'updating bedge blk 10',dtype
        bedge_facey_x(1:nedges,ia:ib,1,1+nguard0*k3d,lb) = 
     .        bedge_facez_x(1:nedges,ia:ib,1+nguard0*k2d,1,lb)
        bedge_facey_x(1:nedges,ia:ib,2,1+nguard0*k3d,lb) = 
     .        bedge_facez_x(1:nedges,ia:ib,nyb+(1+nguard0)*k2d,1,lb)
      elseif(dtype.eq.23) then
        bedge_facey_x(1:nedges,ia:ib,1,nzb+(1+nguard0)*k3d,lb) = 
     .        bedge_facez_x(1:nedges,ia:ib,1+nguard0*k2d,2,lb)
        bedge_facey_x(1:nedges,ia:ib,2,nzb+(1+nguard0)*k3d,lb) = 
     .        bedge_facez_x(1:nedges,ia:ib,nyb+(1+nguard0)*k2d,2,lb)
      endif

      do k = ka,kb
      do j = ja,jb
      do i = ia,ib
        do n=1,nedges
          bedge_facez_y(n,i,j,k,lb) = R_buffer(index)
          index  = index + 1
        enddo
      enddo
      enddo
      enddo
      if(dtype.eq.5) then
        bedge_facex_y(1:nedges,1,ja:jb,1+nguard0*k3d,lb) =
     .        bedge_facez_y(1:nedges,1+nguard0,ja:jb,1,lb)
        bedge_facex_y(1:nedges,2,ja:jb,1+nguard0*k3d,lb) =
     .        bedge_facez_y(1:nedges,nxb+1+nguard0,ja:jb,1,lb)
      elseif(dtype.eq.23) then
        bedge_facex_y(1:nedges,1,ja:jb,nzb+(1+nguard0)*k3d,lb) =
     .        bedge_facez_y(1:nedges,1+nguard0,ja:jb,2,lb)
        bedge_facex_y(1:nedges,2,ja:jb,nzb+(1+nguard0)*k3d,lb) =
     .        bedge_facez_y(1:nedges,nxb+1+nguard0,ja:jb,2,lb)
      endif


      endif
      endif


! Now unpack single edges


! first x edges
      if(ndim.eq.3) then
      if(dtype.eq.2.or.dtype.eq.8.or.dtype.eq.20.or.dtype.eq.26) then

      ia = ia0
      ib = ib0
      j = 1+nguard0
      k = 1

      if(dtype.eq.8) then
        j = nyb+k2d*(1+nguard0)
      elseif(dtype.eq.20) then
        k = 2
      elseif(dtype.eq.26) then
        k = 2
        j = nyb+k2d*(1+nguard0)
      endif
      do i = ia , ib
        do n=1,nedges
          bedge_facez_x(n,i,j,k,lb) = R_buffer(index)
          index  = index + 1
        enddo
      enddo

      ia = ia0
      ib = ib0
      j = 1
      k = 1+nguard0*k3d

      if(dtype.eq.8) then
        j = 2
      elseif(dtype.eq.20) then
        k = nzb+k3d*(1+nguard0)
      elseif(dtype.eq.26) then
        k = nzb+k3d*(1+nguard0)
        j = 2
      endif
      if(mype.eq.0.and.lb.eq.10)write(*,*) 'updating bedge blk 10',dtype
      do i = ia , ib
        do n=1,nedges
          bedge_facey_x(n,i,j,k,lb) = R_buffer(index)
          index  = index + 1
        enddo
      enddo


      endif
      endif


! now y edges

      if(ndim.eq.3) then
      if(dtype.eq.4.or.dtype.eq.6.or.dtype.eq.22.or.dtype.eq.24) then

      ja = ja0
      jb = jb0
      i = 1
      k = 1+nguard0*k3d

      if(dtype.eq.6) then
        i = 2
      elseif(dtype.eq.22) then
        k = nzb+k3d*(1+nguard0)
      elseif(dtype.eq.24) then
        k = nzb+k3d*(1+nguard0)
        i = 2
      endif
      do j = ja , jb
        do n=1,nedges
          bedge_facex_y(n,i,j,k,lb) = R_buffer(index)
          index  = index + 1
        enddo
      enddo

      ja = ja0
      jb = jb0
      i = 1+nguard0
      k = 1

      if(dtype.eq.6) then
        i = nxb+1+nguard0
      elseif(dtype.eq.22) then
        k = 2
      elseif(dtype.eq.24) then
        i = nxb+1+nguard0
        k = 2
      endif
      do j = ja , jb
        do n=1,nedges
          bedge_facez_y(n,i,j,k,lb) = R_buffer(index)
          index  = index + 1
        enddo
      enddo


      endif
      endif


! finally z edges

      if(ndim.eq.3.or.l2p5d.eq.1) then
      if(dtype.eq.10.or.dtype.eq.12.or.dtype.eq.16.or.dtype.eq.18) then

      ka = ka0
      kb = kb0
      i = 1
      j = 1+nguard0*k2d

      if(dtype.eq.12) then
        i = 2
      elseif(dtype.eq.16) then
        j = nyb+k2d*(1+nguard0)
      elseif(dtype.eq.18) then
        i = 2
        j = nyb+k2d*(1+nguard0)
      endif
      do k = ka , kb
        do n=1,nedges
          bedge_facex_z(n,i,j,k,lb) = R_buffer(index)
          index  = index + 1
        enddo
      enddo

      ka = ka0
      kb = kb0
      i = 1+nguard0
      j = 1

      if(dtype.eq.12) then
        i = nxb+1+nguard0
      elseif(dtype.eq.16) then
        j = 2
      elseif(dtype.eq.18) then
        i = nxb+1+nguard0
        j = 2
      endif
      do k = ka , kb
        do n=1,nedges
          bedge_facey_z(n,i,j,k,lb) = R_buffer(index)
          index  = index + 1
        enddo
      enddo


      endif
      endif

! unpack tree info

      do ii = 1,mdim
         coord(ii,lb) = R_buffer(index)
         index = index + 1
      end do

      do ii = 1,mdim
         bsize(ii,lb) = R_buffer(index)
         index = index + 1
      end do

      do ii = 1,2
         do jj = 1,mdim
            bnd_box(ii,jj,lb) = R_buffer(index)
            index = index + 1
         end do
      end do

      do ii = 1,2
         parent(ii,lb) = int(R_buffer(index))
         index = index + 1
      end do

      do ii = 1,2
         do jj = 1,mchild
            child(ii,jj,lb) = R_buffer(index)
            index = index + 1
         end do
      end do

      which_child(lb) = int(R_buffer(index))
      index = index + 1

      if (int(R_buffer(index)) == 1) then
         newchild(lb) = .true.
      else
         newchild(lb) = .false.
      end if
      index = index + 1

      do ii = 1,2
         do jj = 1,mfaces
            neigh(ii,jj,lb) = int(R_buffer(index))
            index = index + 1
         end do
      end do

      lrefine(lb) = int(R_buffer(index))
      index = index + 1

      nodetype(lb) = int(R_buffer(index))
      index = index + 1

      do ii = 1,3
         do k = 1,1+2*k3d
            do jj = 1,1+2*k2d
               do iii = 1,3
                  surr_blks(ii,iii,jj,k,lb) =
     &                 int(R_buffer(index))
                  index = index + 1
               end do
            end do
         end do
      end do


#ifdef SAVE_MORTS
      do ii = 1,6
         do k = 1,1+2*k3d
            do jj = 1,1+2*k2d
               do iii = 1,3
                  surr_morts(ii,iii,jj,k,lb) =
     &                 int(R_buffer(index))
                  index = index + 1
               end do
            end do
         end do
      end do
#endif /* SAVE_MORTS */

      empty(lb) = int(R_buffer(index))
      index = index + 1

      if (int(R_buffer(index)) == 1) then
         ldtcomplete(lb) = .true.
      else
         ldtcomplete(lb) = .false.
      end if
      index = index + 1


! overwrite outgoing index offset

      offset = index

      return
      end subroutine mpi_put_edge_buffer
