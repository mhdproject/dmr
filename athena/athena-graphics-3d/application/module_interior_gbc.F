#include "paramesh_preprocessor.fh"
#include "athena_preprocessor.fh"
      module interior_gbc
      use paramesh_dimensions
      private
c      save
      public earth_boundary
      public earth_cells
      public interior_bc_states_n
      public interior_bc_states_nplus_1half
      public interior_bc_fl
      public internal_hardwall_bc
      public bubble_sort
      
      integer,parameter, public              :: maxzones = 4096
      integer,dimension(3,maxzones,maxblocks),public::p_zone,p_rzone 
      integer,dimension(4,maxzones,maxblocks),public::p_face,p_rface 
      integer,dimension(3,maxzones,maxblocks),public::p_efields 
      integer,dimension(maxblocks),public::idzones,irdzones 
      integer,dimension(maxblocks),public::idfaces,irdfaces
      integer,dimension(maxblocks),public::idefields
      integer,dimension(maxblocks),public::boundary_blocks
      integer,dimension(maxblocks),public::domain_blocks
      logical,dimension(maxblocks),public::earth_blocks
      logical,dimension(iu_bnd1,ju_bnd1,ku_bnd1,
     &                               maxblocks_alloc),public::learth_g
      logical,dimension(iu_bnd1,ju_bnd1,ku_bnd1),public::learth
      integer, public:: iearth,idomain,iloop
      real, public:: ebox_xmin,ebox_xmax,ebox_ymin,ebox_ymax
      real, public:: ebox_zmin,ebox_zmax
      real, save,public:: roa,arfact2
!     real, public:: roa,arfact2,iearth,idomain
      real,dimension(il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1,      
     &                  maxblocks_alloc),public::bmask
      real,dimension(il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,      
     &                  maxblocks_alloc),public::bmask_rb
      integer,public:: izones,irzones,ifaces,irfaces,iedges
      integer, dimension(4,1000),public:: iface,irface
      integer, dimension(3,1000),public:: izone
      integer, dimension(3,1000),public:: iedge
      integer, dimension(2,1000),public:: iface_1d
      integer, dimension(  1000),public:: izone_1d
      integer                   ,public:: ifaces_1d,izones_1d
      integer                   ,public:: izones_min,izones_max
      integer, dimension(3,1000),public:: irzone
      contains
!---------------------------------------------------------------
!
!---------------------------------------------------------------
      subroutine earth_boundary(mype)
      use physcons
      use paramesh_dimensions
      use physicaldata
      use tree
      use grid
      use strongfields
! include files for amr


#include "amr_shmem.fh"
#include "shmem_reduce.fh"
#include "ionosphere.fh"
!--------------------------------------------------------------

      integer, parameter:: ng0= nguard*npgs
      integer, dimension(maxblocks):: bpoint,rbpoint
c     implicit none
      integer :: ix,iy,iz,ibpoint,irbpoint,mype
      real    :: xpos,ypos,zpos,rpos2,rpos
      real    :: xlo,xhi,ylo,yhi,zlo,zhi
      real    :: xll,xhh,yll,yhh,zll,zhh
      real    :: xmin,xmax,ymin,ymax,zmin,zmax,factor
      real    :: xmincc,xmaxcc,ymincc,ymaxcc,zmincc,zmaxcc
      real    :: xoff_set,yoff_set,zoff_set
      logical :: lblock,lrblock
      integer :: errcode,nprocs

#ifdef MPI_USED
      include "mpif.h"
!     integer :: errcode
#endif /* MPI_USED */
!
#ifdef MPI_USED
        call mpi_barrier (MPI_COMM_WORLD, errcode)
#else
        call shmem_barrier_all()
#endif /*MPI_USED*/
      xhi = -1.e30
      xlo =  1.e30
      yhi = -1.e30
      ylo =  1.e30
      zhi = -1.e30
      zlo =  1.e30
      ebox_xmax = -1.e30
      ebox_xmin =  1.e30
      ebox_ymax = -1.e30
      ebox_ymin =  1.e30
      ebox_zmax = -1.e30
      ebox_zmin =  1.e30
      ibpoint       = 0
      irbpoint      = 0
c      write(*,*) 'irbpoint', irbpoint
      bflags(1,:)   =-1.0
      bmask(:,:,:,:)= 1.0
      bmask_rb(:,:,:,:)= 1.0


#ifdef SPHERE
      if(lnblocks.gt.0) then
      do lb=1,lnblocks
!      write(*,*) 'SPHERE: starting block ',lb,' lnblocks ',lnblocks
#ifndef ADVANCE_ALL_LEVELS
      if(nodetype(lb).eq.1) then
#endif
!         write(*,*)'earth_boundary:sphere: lb=',lb
        call grid_variables(lb,0.0)

        xmax  = coord(1,lb) + 0.5*bsize(1,lb)+nguard*delx
        xmin  = coord(1,lb) - 0.5*bsize(1,lb)-nguard*delx

        ymax  = coord(2,lb) + 0.5*bsize(2,lb)+nguard*dely
        ymin  = coord(2,lb) - 0.5*bsize(2,lb)-nguard*dely
#if N_DIM == 3
        zmax  = coord(3,lb) + 0.5*bsize(3,lb)+nguard*delz
        zmin  = coord(3,lb) - 0.5*bsize(3,lb)-nguard*delz
#else
        zmax  = 0.
        zmin  = 0.
        delz  = 0.
#endif
!      write(*,*) 'SPHERE: block ',lb, roa,
!    &                xmin,xmax,ymin,ymax,zmin,zmax
!
        lblock = .true.
        do iz=kl_bnd1,ku_bnd1
         do iy=jl_bnd1,ju_bnd1
          do ix=il_bnd1,iu_bnd1
!xpos = xmin + delx*(ix-0.5)+x0
           xpos = xmin + delx*(ix-0.5)
           ypos = ymin + dely*(iy-0.5)

#if N_DIM == 3
           zpos = zmin + delz*(iz-0.5)
#else
           zpos = 0.
#endif              
           rpos2 =(xpos**2) + (ypos**2) + (zpos**2) 
           rpos  = sqrt(rpos2)
           if(rpos.le.roa)then
            bmask(ix,iy,iz,lb)=0.0
            if(lblock)then
             ibpoint = ibpoint +1
             bpoint(ibpoint) = lb
              write(*,*)'SPHERE: ibpoint,lb',ibpoint,lb
             lblock = .false.
            endif
           endif        
          enddo
         enddo
        enddo

        xoff_set = nguard*delx*float(k1d)
        yoff_set = nguard*dely*float(k2d)
        zoff_set = nguard*delz*float(k3d)
        xmincc = xmin+xoff_set
        xmaxcc = xmax-xoff_set
        ymincc = ymin+yoff_set
        ymaxcc = ymax-yoff_set
        zmincc = zmin+zoff_set
        zmaxcc = zmax-zoff_set
        lrblock = .true.
        do iz=1,nzb
         do iy=1,nyb
          do ix=1,nxb
           xpos = xmincc + delx*(ix-0.5)
           ypos = ymincc + dely*(iy-0.5)

#if N_DIM == 3
           zpos = zmincc + delz*(iz-0.5)
#else
           zpos = 0.
#endif              
           rpos2 =(xpos**2) + (ypos**2) + (zpos**2) 
           rpos  = sqrt(rpos2)
           if(rpos.le.roa)then
            bmask_rb(ix,iy,iz,lb)=0.0
            if(lrblock)then
c               write(*,*) 'irbpoint', irbpoint
             irbpoint = irbpoint +1
c               write(*,*) 'irbpoint', irbpoint
             rbpoint(irbpoint) = lb
              write(*,*)'SPHERE: irbpoint,lb',irbpoint,lb
             lrblock = .false.
            endif
           endif        
          enddo
         enddo
        enddo

#ifndef ADVANCE_ALL_LEVELS
      endif
#endif
c      write(*,*) 'init: finished block ',lb
      enddo ! end loop over grid blocks
      endif


#endif /*SPHERE*/

#ifdef CUBE
     write(*,*)'earth_boundary-----------------roa = ',roa,mype

!e-box
      xhi = -1.e30; xlo =  1.e30;
      yhi = -1.e30; ylo =  1.e30;
      zhi = -1.e30; zlo =  1.e30;
      xhh = -1.e30; xll =  1.e30;
      yhh = -1.e30; yll =  1.e30;
      zhh = -1.e30; zll =  1.e30;

!-----------------------------------------------------------------
      if(lnblocks.gt.0) then
         do lb=1,lnblocks

#ifndef ADVANCE_ALL_LEVELS
           if(nodetype(lb).eq.1) then
#endif
              call grid_variables(lb,0.0)
              xmax  = coord(1,lb) + 0.5*bsize(1,lb)+nguard*delx
              xmin  = coord(1,lb) - 0.5*bsize(1,lb)-nguard*delx

              ymax  = coord(2,lb) + 0.5*bsize(2,lb)+nguard*dely
              ymin  = coord(2,lb) - 0.5*bsize(2,lb)-nguard*dely
#if N_DIM == 3
              zmax  = coord(3,lb) + 0.5*bsize(3,lb)+nguard*delz
              zmin  = coord(3,lb) - 0.5*bsize(3,lb)-nguard*delz
#else
              zmax  = 0.
              zmin  = 0.
              delz  = 0.
#endif

              do iz=kl_bnd1,ku_bnd1
                 do iy=jl_bnd1,ju_bnd1
                    do ix=il_bnd1,iu_bnd1

                       xpos = xmin + delx*(ix-0.5)
                       ypos = ymin + dely*(iy-0.5)

#if N_DIM == 3
                       zpos = zmin + delz*(iz-0.5)
#else
                       zpos = 0.
#endif 
#if N_DIM < 3
                 if(abs(xpos).le.roa.and.abs(ypos).le.roa)then
                       xlo = min(xlo, (xpos-.5*delx))
                       xhi = max(xhi, (xpos+.5*delx))
                       ylo = min(ylo, (ypos-.5*dely))
                       yhi = max(yhi, (ypos+.5*dely))
                       zlo = min(zlo, (zpos-.5*delz))
                       zhi = max(zhi, (zpos+.5*delz))
                 endif 
             
#else
                 if(abs(xpos).le.roa.and.abs(ypos).le.roa.and.
     &           abs(zpos).le.roa)then
                       xlo = min(xlo, (xpos-.5*delx))
                       xhi = max(xhi, (xpos+.5*delx))
                       ylo = min(ylo, (ypos-.5*dely))
                       yhi = max(yhi, (ypos+.5*dely))
                       zlo = min(zlo, (zpos-.5*delz))
                       zhi = max(zhi, (zpos+.5*delz))
                 endif

#endif /* NDIM3 */

                    enddo
                 enddo
              enddo

              xoff_set = nguard*delx*float(k1d)
              yoff_set = nguard*dely*float(k2d)
              zoff_set = nguard*delz*float(k3d)
              xmincc   = xmin+xoff_set
              xmaxcc   = xmax-xoff_set
              ymincc   = ymin+yoff_set
              ymaxcc   = ymax-yoff_set
              zmincc   = zmin+zoff_set
              zmaxcc   = zmax-zoff_set

              do iz=kl_bnd,ku_bnd
                 do iy=jl_bnd,ju_bnd
                    do ix=il_bnd,iu_bnd

                       xpos = xmincc + delx*(ix-0.5)
                       ypos = ymincc + dely*(iy-0.5)

#if N_DIM == 3
                       zpos = zmincc + delz*(iz-0.5)
#else
                       zpos = 0.
#endif 
#if N_DIM < 3
                 if(abs(xpos).le.roa.and.abs(ypos).le.roa)then
                       xll = min(xll, (xpos-.5*delx))
                       xhh = max(xhh, (xpos+.5*delx))
                       yll = min(yll, (ypos-.5*dely))
                       yhh = max(yhh, (ypos+.5*dely))
                       zll = min(zll, (zpos-.5*delz))
                       zhh = max(zhh, (zpos+.5*delz))
                 endif 
             
#else
                 if(abs(xpos).le.roa.and.abs(ypos).le.roa.and.
     &           abs(zpos).le.roa)then
                       xll = min(xll, (xpos-.5*delx))
                       xhh = max(xhh, (xpos+.5*delx))
                       yll = min(yll, (ypos-.5*dely))
                       yhh = max(yhh, (ypos+.5*dely))
                       zll = min(zll, (zpos-.5*delz))
                       zhh = max(zhh, (zpos+.5*delz))
                 endif

#endif /* NDIM3 */

                    enddo
                 enddo
              enddo
#
#ifndef ADVANCE_ALL_LEVELS
           endif
#endif
         enddo ! end loop over grid blocks
      endif

!----------------------------------------------------------------------

      ebox_xmin = xlo; ebox_xmax = xhi;
      ebox_ymin = ylo; ebox_ymax = yhi;
      ebox_zmin = zlo; ebox_zmax = zhi;


      call shmem_barrier_all()
      call comm_real_min_to_all(ebox_xmin,ebox_xmin)
      call shmem_barrier_all()
      call comm_real_min_to_all(ebox_ymin,ebox_ymin)
      call shmem_barrier_all()
      call comm_real_min_to_all(ebox_zmin,ebox_zmin)
      call shmem_barrier_all()
      call comm_real_max_to_all(ebox_xmax,ebox_xmax)
      call shmem_barrier_all()
      call comm_real_max_to_all(ebox_ymax,ebox_ymax)
      call shmem_barrier_all()
      call comm_real_max_to_all(ebox_zmax,ebox_zmax)

      call shmem_barrier_all()

      if(mype.eq.0) then
         write(*,*)'earth: ebox_xmin,ebox_xmax=',ebox_xmin,ebox_xmax
         write(*,*)'earth: ebox_ymin,ebox_ymax=',ebox_ymin,ebox_ymax
         write(*,*)'earth: ebox_zmin,ebox_zmax=',ebox_zmin,ebox_zmax
      endif

!----------------------------------------------------------------------

!set bmax and bpoint

      if(lnblocks.gt.0) then
         do lb=1,lnblocks

#ifndef ADVANCE_ALL_LEVELS
           if(nodetype(lb).eq.1) then
#endif
              call grid_variables(lb,0.0)
              xmax  = coord(1,lb) + 0.5*bsize(1,lb)+nguard*delx
              xmin  = coord(1,lb) - 0.5*bsize(1,lb)-nguard*delx

              ymax  = coord(2,lb) + 0.5*bsize(2,lb)+nguard*dely
              ymin  = coord(2,lb) - 0.5*bsize(2,lb)-nguard*dely
#if N_DIM == 3
              zmax  = coord(3,lb) + 0.5*bsize(3,lb)+nguard*delz
              zmin  = coord(3,lb) - 0.5*bsize(3,lb)-nguard*delz
#else
              zmax  = 0.
              zmin  = 0.
              delz  = 0.
#endif

              lblock = .true.
              do iz=kl_bnd1,ku_bnd1
                 do iy=jl_bnd1,ju_bnd1
                    do ix=il_bnd1,iu_bnd1

                       xpos = xmin + delx*(ix-0.5)
                       ypos = ymin + dely*(iy-0.5)

#if N_DIM == 3
                       zpos = zmin + delz*(iz-0.5)*k3d
#else
                       zpos = 0.
#endif 
#if N_DIM < 3
                       if(xpos.ge.xlo.and.xpos.le.xhi.and.
     &                    ypos.ge.ylo.and.ypos.le.yhi    )then

                              bmask(ix,iy,iz,lb)=0.0
                              if(lblock)then
                                 ibpoint = ibpoint +1
                                 bpoint(ibpoint) = lb
                                 lblock = .false.
                              endif
                       endif        
             
#else
                       if(xpos.ge.xlo.and.xpos.le.xhi.and.
     &                 ypos.ge.ylo.and.ypos.le.yhi.and.    
     &                 zpos.ge.zlo.and.zpos.le.zhi    )then
                              bmask(ix,iy,iz,lb)=0.0
                              if(lblock)then
                                 ibpoint = ibpoint +1
                                 bpoint(ibpoint) = lb
                                 lblock = .false.
                              endif
                       endif        

#endif /* NDIM3 */
                    enddo
                 enddo
              enddo

             xoff_set = nguard*delx*float(k1d)
             yoff_set = nguard*dely*float(k2d)
             zoff_set = nguard*delz*float(k3d)
             xmincc   = xmin+xoff_set
             xmaxcc   = xmax-xoff_set
             ymincc   = ymin+yoff_set
             ymaxcc   = ymax-yoff_set
             zmincc   = zmin+zoff_set
             zmaxcc   = zmax-zoff_set

              lrblock = .true.
              do iz=kl_bnd,ku_bnd
                 do iy=jl_bnd,ju_bnd
                    do ix=il_bnd,iu_bnd

                       xpos = xmincc + delx*(ix-0.5)
                       ypos = ymincc + dely*(iy-0.5)

#if N_DIM == 3
                       zpos = zmincc + delz*(iz-0.5)*k3d
#else
                       zpos = 0.
#endif 
#if N_DIM < 3
                       if(xpos.ge.xll.and.xpos.le.xhh.and.
     &                    ypos.ge.yll.and.ypos.le.yhh    )then

                              bmask_rb(ix,iy,iz,lb)=0.0
                              if(lrblock)then
                                 irbpoint = irbpoint +1
                                 rbpoint(irbpoint) = lb
                                 lrblock = .false.
                              endif
                       endif        
             
#else
                       if(xpos.ge.xll.and.xpos.le.xhh.and.
     &                 ypos.ge.yll.and.ypos.le.yhh.and.    
     &                 zpos.ge.zll.and.zpos.le.zhh    )then
                              bmask_rb(ix,iy,iz,lb)=0.0
                              if(lrblock)then
                                 irbpoint = irbpoint +1
                                 rbpoint(irbpoint) = lb
                                 lrblock = .false.
                              endif
                       endif        

#endif /* NDIM3 */
                    enddo
                 enddo
              enddo


#ifndef ADVANCE_ALL_LEVELS
           endif
#endif
         enddo ! end loop over grid blocks
      endif
#endif /*CUBE*/

!----------------------------------------------------------------------


      idfaces(:)   = 0
      idzones(:)   = 0
      irdfaces(:)  = 0
      irdzones(:)  = 0
      earth_blocks(:) = .false.

      do ii=1,ibpoint
      lb = bpoint(ii)
      earth_blocks(lb) = .true.


      do k = kl_bnd1,ku_bnd1
       do j = jl_bnd1,ju_bnd1
        do i = il_bnd1+1,iu_bnd1
         if(bmask(i,j,k,lb).eq.0.0.and.bmask(i-1,j,k,lb).eq.1.0)then
          idfaces(lb) = idfaces(lb)+1
          p_face(1,idfaces(lb),lb) = i-1
          p_face(2,idfaces(lb),lb) = j
          p_face(3,idfaces(lb),lb) = k
          p_face(4,idfaces(lb),lb) = 2
c        write(*,*)'i,j,k,p_face(4,idfaces(lb),lb)=',i,j,k,
c     &p_face(4,idfaces(lb),lb)
         elseif(bmask(i,j,k,lb).eq.1.0.and.bmask(i-1,j,k,lb).eq.0.0)then
          idfaces(lb) = idfaces(lb)+1
          p_face(1,idfaces(lb),lb) = i
          p_face(2,idfaces(lb),lb) = j
          p_face(3,idfaces(lb),lb) = k
          p_face(4,idfaces(lb),lb) = 1
c       write(*,*)'i,j,k,p_face(4,idfaces(lb),lb)=',i,j,k,
c     &p_face(4,idfaces(lb),lb)
         endif
          
        enddo
       enddo
      enddo

      do k = kl_bnd1,ku_bnd1
       do i = il_bnd1,iu_bnd1
        do j = jl_bnd1+1,ju_bnd1
         if(bmask(i,j,k,lb).eq.0.0.and.bmask(i,j-1,k,lb).eq.1.0)then
          idfaces(lb) = idfaces(lb)+1
          p_face(1,idfaces(lb),lb) = i
          p_face(2,idfaces(lb),lb) = j-1
          p_face(3,idfaces(lb),lb) = k
          p_face(4,idfaces(lb),lb) = 4
c       write(*,*)'i,j,k,p_face(4,idfaces(lb),lb)=',i,j,k,
c     &p_face(4,idfaces(lb),lb)
         elseif(bmask(i,j,k,lb).eq.1.0.and.bmask(i,j-1,k,lb).eq.0.0)then
          idfaces(lb) = idfaces(lb)+1
          p_face(1,idfaces(lb),lb) = i
          p_face(2,idfaces(lb),lb) = j
          p_face(3,idfaces(lb),lb) = k
          p_face(4,idfaces(lb),lb) = 3
c       write(*,*)'i,j,k,p_face(4,idfaces(lb),lb)=',i,j,k,
c     &p_face(4,idfaces(lb),lb)
         endif
          
        enddo
       enddo
      enddo

      do i = il_bnd1,iu_bnd1
       do j = jl_bnd1,ju_bnd1
        do k = kl_bnd1+k3d,ku_bnd1
         if(bmask(i,j,k,lb).eq.0.0.and.bmask(i,j,k-k3d,lb).eq.1.0)then
          idfaces(lb) = idfaces(lb)+1
          p_face(1,idfaces(lb),lb) = i
          p_face(2,idfaces(lb),lb) = j
          p_face(3,idfaces(lb),lb) = k-k3d
          p_face(4,idfaces(lb),lb) = 6
c      write(*,*)'i,j,k,p_face(4,idfaces(lb),lb)=',i,j,k,
c     &p_face(4,idfaces(lb),lb)
      elseif(bmask(i,j,k,lb).eq.1.0.and.bmask(i,j,k-k3d,lb).eq.0.0)then
          idfaces(lb) = idfaces(lb)+1
          p_face(1,idfaces(lb),lb) = i
          p_face(2,idfaces(lb),lb) = j
          p_face(3,idfaces(lb),lb) = k
          p_face(4,idfaces(lb),lb) = 5
c      write(*,*)'i,j,k,p_face(4,idfaces(lb),lb)=',i,j,k,
c     &iface(4,idfaces(lb),lb)
         endif          
        enddo
       enddo
      enddo
      learth_g(:,:,:,lb)=.false.
      do i = il_bnd1,iu_bnd1
       do j = jl_bnd1,ju_bnd1
        do k = kl_bnd1+k3d,ku_bnd1
         if(bmask(i,j,k,lb).eq.0.0)learth_g(i,j,k,lb)=.true.
        enddo
       enddo
      enddo

c
c-----real blocks
c
      do k = kl_bnd,ku_bnd
       do j = jl_bnd,ju_bnd
        do i = il_bnd+1,iu_bnd
         if(bmask_rb(i,j,k,lb).eq.0.0.and.bmask_rb(i-1,j,k,lb).eq.1.0)then
          irdfaces(lb) = irdfaces(lb)+1
          p_rface(1,irdfaces(lb),lb) = i-1
          p_rface(2,irdfaces(lb),lb) = j
          p_rface(3,irdfaces(lb),lb) = k
          p_rface(4,irdfaces(lb),lb) = 2
        write(*,*)'i,j,k,p_rface(4,irdfaces(lb),lb)=',i,j,k,
     &p_rface(4,irdfaces(lb),lb)
         elseif(bmask_rb(i,j,k,lb).eq.1.0.and.bmask_rb(i-1,j,k,lb).eq.0.0)then
          irdfaces(lb) = irdfaces(lb)+1
          p_rface(1,irdfaces(lb),lb) = i
          p_rface(2,irdfaces(lb),lb) = j
          p_rface(3,irdfaces(lb),lb) = k
          p_rface(4,irdfaces(lb),lb) = 1
       write(*,*)'i,j,k,p_rface(4,irdfaces(lb),lb)=',i,j,k,
     &p_rface(4,irdfaces(lb),lb)
         endif
          
        enddo
       enddo
      enddo

      do k = kl_bnd,ku_bnd
       do i = il_bnd,iu_bnd
        do j = jl_bnd+1,ju_bnd
         if(bmask_rb(i,j,k,lb).eq.0.0.and.bmask_rb(i,j-1,k,lb).eq.1.0)then
          irdfaces(lb) = irdfaces(lb)+1
          p_rface(1,irdfaces(lb),lb) = i
          p_rface(2,irdfaces(lb),lb) = j-1
          p_rface(3,irdfaces(lb),lb) = k
          p_rface(4,irdfaces(lb),lb) = 4
c       write(*,*)'i,j,k,p_rface(4,irdfaces(lb),lb)=',i,j,k,
c     &p_rface(4,irdfaces(lb),lb)
         elseif(bmask_rb(i,j,k,lb).eq.1.0.and.bmask_rb(i,j-1,k,lb).eq.0.0)then
          irdfaces(lb) = irdfaces(lb)+1
          p_rface(1,irdfaces(lb),lb) = i
          p_rface(2,irdfaces(lb),lb) = j
          p_rface(3,irdfaces(lb),lb) = k
          p_rface(4,irdfaces(lb),lb) = 3
c       write(*,*)'i,j,k,p_rface(4,irdfaces(lb),lb)=',i,j,k,
c     &p_rface(4,irdfaces(lb),lb)
         endif
          
        enddo
       enddo
      enddo

      do i = il_bnd,iu_bnd
       do j = jl_bnd,ju_bnd
        do k = kl_bnd+k3d,ku_bnd
         if(bmask_rb(i,j,k,lb).eq.0.0.and.bmask_rb(i,j,k-k3d,lb).eq.1.0)then
          irdfaces(lb) = irdfaces(lb)+1
          p_rface(1,irdfaces(lb),lb) = i
          p_rface(2,irdfaces(lb),lb) = j
          p_rface(3,irdfaces(lb),lb) = k-k3d
          p_rface(4,irdfaces(lb),lb) = 6
c      write(*,*)'i,j,k,p_rface(4,irdfaces(lb),lb)=',i,j,k,
c     &p_rface(4,irdfaces(lb),lb)
      elseif(bmask_rb(i,j,k,lb).eq.1.0.and.bmask_rb(i,j,k-k3d,lb).eq.0.0)then
          irdfaces(lb) = irdfaces(lb)+1
          p_rface(1,irdfaces(lb),lb) = i
          p_rface(2,irdfaces(lb),lb) = j
          p_rface(3,irdfaces(lb),lb) = k
          p_rface(4,irdfaces(lb),lb) = 5
c      write(*,*)'i,j,k,p_rface(4,irdfaces(lb),lb)=',i,j,k,
c     &iface(4,irdfaces(lb),lb)
         endif          
        enddo
       enddo
      enddo


      do k = kl_bnd1,ku_bnd1
       do j = jl_bnd1,ju_bnd1
        do i = il_bnd1,iu_bnd1
         if(bmask(i,j,k,lb).eq.0.0) then
          idzones(lb) = idzones(lb)+1
          p_zone(1,idzones(lb),lb) = i
          p_zone(2,idzones(lb),lb) = j 
          p_zone(3,idzones(lb),lb) = k
c          write(*,*)'iz,jz,kz=',i,j,k
         endif
        enddo
       enddo
      enddo
      if(mype.eq.0)
     &write(*,*)'idzones(lb),idfaces(lb)=',idzones(lb),idfaces(lb)

      do k = kl_bnd,ku_bnd
       do j = jl_bnd,ju_bnd
        do i = il_bnd,iu_bnd
         if(bmask_rb(i,j,k,lb).eq.0.0) then
          irdzones(lb) = irdzones(lb)+1
          p_rzone(1,irdzones(lb),lb) = i
          p_rzone(2,irdzones(lb),lb) = j 
          p_rzone(3,irdzones(lb),lb) = k
c          write(*,*)'iz,jz,kz=',i,j,k
         endif
        enddo
       enddo
      enddo
      if(mype.eq.0)
     &write(*,*)'irdzones(lb)=',irdzones(lb)

      enddo



       do ii = 1,ibpoint
        if(mype.eq.0)
     &  write(*,*)'ii,bpoint(ii)=',ii,bpoint(ii)
       enddo
       do ii=1,ibpoint
        lb = bpoint(ii)
        if(lb.eq.175)then
         do jj = 1,idfaces(lb)
         i  = p_face(1,jj,lb)
         j  = p_face(2,jj,lb)
         k  = p_face(3,jj,lb)
         if = p_face(4,jj,lb)
         write(*,*)'lb,i,j,k,if=',lb,i,j,k,if
         enddo
c         stop
        endif
       enddo

       do ii=1,ibpoint
        lb = bpoint(ii)
        bflags(1,lb)=1.0
       enddo


#ifdef MPI_USED
        call mpi_barrier (MPI_COMM_WORLD, errcode)
#else
        call shmem_barrier_all()
#endif /*MPI_USED*/

c      deallocate(bmask)
      return
      end subroutine earth_boundary
!---------------------------------------------------------------
!
!---------------------------------------------------------------
!----------------------------------------------------------
!
!----------------------------------------------------------
      subroutine earth_cells(lb)


!----------------------------------------------------------
!
! This routine identifies any cells in block lb which lie
! on the square earth bounding box, and whether these are
! cut, and provides the outward normal vector to the earth
! surface in that cell. The terms used to describe cut cells
! are as follows:
!   earth cell  - a cell entirely inside the earth boundary
!   edge cell   - a cell which is sliced along a diagonal
!   corner cell - a cell which is sliced along all 3 diagonals
!   wall cell   - a cell outside the earth boundary which shares
!                 1 face with the earth boundary
!
!

!----------------------------------------------------------
      use paramesh_dimensions
      use physicaldata
      use tree
      implicit none

!----------------------------------------------------------

      integer ixnorm(iu_bnd1,ju_bnd1,ku_bnd1)
      integer iynorm(iu_bnd1,ju_bnd1,ku_bnd1)
      integer iznorm(iu_bnd1,ju_bnd1,ku_bnd1)
      logical lcorner(iu_bnd1,ju_bnd1,ku_bnd1)
      logical ledge(iu_bnd1,ju_bnd1,ku_bnd1)
      logical lwall(iu_bnd1,ju_bnd1,ku_bnd1)
      integer lb

!----------------------------------------------------------

      real    dmask(il_bnd1-1:iu_bnd1+1,jl_bnd1-k2d:ju_bnd1+k2d,
     .              kl_bnd1-k3d:ku_bnd1+k3d)
      integer nrwall,i,j,k
      real    rho_max,xpos,ypos,zpos,delx,dely,delz
      real    xtest,ytest,ztest,xearth,yearth,zearth
      real    xmin,ymin,zmin,xmax,ymax,zmax
      real    w_ebox_x,w_ebox_y,w_ebox_z,w_ebox
      real    x0,y0,z0,dmask_tot
      integer ix,iy,iz,ixn,iyn,izn
      integer ncells
      real    ncells_buffer
!----------------------------------------------------------

      lcorner(:,:,:) = .false.
      ledge(:,:,:)   = .false.
      lwall(:,:,:)   = .false.
      ixnorm(:,:,:)  = 0.
      iynorm(:,:,:)  = 0.
      iznorm(:,:,:)  = 0.

! Initialize counters for each type of cell pointer
      idzones(lb)    = 0
      idfaces(lb)    = 0
      irdfaces(lb)   = 0
      idefields(lb)  = 0
      learth(:,:,:)  = .false.

      xmin = bnd_box(1,1,lb)
      xmax = bnd_box(2,1,lb)
      ymin = bnd_box(1,2,lb)
      ymax = bnd_box(2,2,lb)
      zmin = bnd_box(1,3,lb)
      zmax = bnd_box(2,3,lb)

      write(*,*)'xmin,xmax=',xmin,xmax
      write(*,*)'ymin,ymax=',ymin,ymax
      write(*,*)'zmin,zmax=',zmin,zmax
      write(*,*)'ebox_xmax-ebox_xmin=',ebox_xmax-ebox_xmin
      write(*,*)'ebox_ymax-ebox_ymin=',ebox_ymax-ebox_ymin
      write(*,*)'ebox_zmax-ebox_zmin=',ebox_zmax-ebox_zmin


! Step 1
! Fill mask array, including guardcells
      dmask(:,:,:) = 0.


! Get grid cell sizes
      delx = bsize(1,lb)/real(nxb)
      dely = bsize(2,lb)/real(nyb)
      if(ndim.eq.3)then
       delz = bsize(3,lb)/real(nzb)
      else
       delz = 0.
      endif
      write(*,*)'delx,dely,delz=',delx,dely,delz
c      pause
!---
! Error trapping. Test size of earth box.
!---
      xearth = 0.
      yearth = 0.
      zearth = 0.
      xtest = (xearth-xmin)*(xearth-xmax)
      ytest = (yearth-ymin)*(yearth-ymax)
      if(ndim.eq.3)then
        ztest=(zearth-zmin)*(zearth-zmax)
      else
       ztest = -1.
      endif
      if(xtest.lt.0.and.ytest.lt.0.and.ztest.lt.0.) then
        w_ebox_x = abs(ebox_xmax-ebox_xmin)/delx
        w_ebox_y = abs(ebox_ymax-ebox_ymin)/dely
        if(ndim.eq.3) then
          w_ebox_z = abs(ebox_zmax-ebox_zmin)/delz
        else
          w_ebox_z = 3.
        endif
        w_ebox = min( min(w_ebox_x,w_ebox_y), w_ebox_z)
        if(w_ebox.lt.1.9) then
          write(*,*) 'Earth box too small!!'
          write(*,*) 'lb dels w_ebox ',lb,delx,dely,delz,w_ebox_x,
     .                               w_ebox_y,w_ebox_z
          write(*,*) 'bnd_box ',bnd_box(:,:,lb)
          stop
         endif
       endif
!---


!---
! set dmask to 1.0 inside earth box
!
      dmask_tot = 0.
      do k = kl_bnd1-k3d,ku_bnd1+k3d
       do j = jl_bnd1-k2d,ju_bnd1+k2d
        do i = il_bnd1-1  ,iu_bnd1+1
         xpos = xmin+delx*(real(i-nguard    )-.5)
         ypos = ymin+dely*(real(j-nguard*k2d)-.5)
         if(ndim.eq.3)then
          zpos = zmin+delz*(real(k-nguard*k3d)-.5)
         else
          zpos = 0.
         endif

         xtest = (xpos-ebox_xmin)*(xpos-ebox_xmax)
         ytest = (ypos-ebox_ymin)*(ypos-ebox_ymax)
         if(ndim.eq.3)then
          ztest = (zpos-ebox_zmin)*(zpos-ebox_zmax)
         else
          ztest = -1.
         endif
        if(xtest.le.0..and.ytest.le.0..and.ztest.le.0.) then
          dmask(i,j,k) = 1.
          dmask_tot    = dmask_tot + 1.0
!        write(*,*) 'dmask ',i,j,k,dmask(i,j,k),lb
         endif
        enddo
       enddo
      enddo
!
      write(*,*) 'dmask_tot ',dmask_tot

!---
! Step 2
! Flag corners and edges of earth box, cells entirely within the
! earth boundary, and cells in flow domain which abut the earth boundary
! on a single face.
      lcorner(:,:,:)    = .false.
      ledge(:,:,:)      = .false.
      lwall(:,:,:)      = .false.
      learth(:,:,:)  = .false.

! detect edges, corners or cells inside earth
      do k = kl_bnd1,ku_bnd1
      do j = jl_bnd1,ju_bnd1
      do i = il_bnd1,iu_bnd1
        nrwall =
     .  int( dmask(i-1,j    ,k    ) + dmask(i+1,j    ,k    )
     .     + dmask(i  ,j-k2d,k    ) + dmask(i  ,j+k2d,k    )
     .     + dmask(i  ,j    ,k-k3d) + dmask(i  ,j    ,k+k3d) + .1 )

      if(dmask_tot.gt.1.5) then
       if(nrwall.ge.3) learth(i,j,k)   = .true.
       if(nrwall.eq.1) lwall(i,j,k    )   = .true.
      elseif(dmask_tot.eq.1.) then
       if(dmask(i,j,k).eq.1.) learth(i,j,k) = .true.
       if(nrwall.eq.1) lwall(i,j,k)   = .true.
c       if(nrwall.eq.1) write(*,*)'wall detected at ',i,j,k
      endif

      enddo
      enddo
      enddo
!---


!---
! Step 5
! Compute direction of normal of cut surface (directed outward from earth)
! for wall cells
      do k = kl_bnd1,ku_bnd1
      do j = jl_bnd1,ju_bnd1
      do i = il_bnd1,iu_bnd1
      if(lwall(i,j,k)) then
        ixnorm(i,j,k) = int( dmask(i-1,j    ,k    )+.1 )
     .                - int( dmask(i+1,j    ,k    )+.1 )
        iynorm(i,j,k) = int( dmask(i  ,j-k2d,k    )+.1 )
     .                - int( dmask(i  ,j+k2d,k    )+.1 )
        iznorm(i,j,k) = int( dmask(i  ,j    ,k-k3d)+.1 )
     .                - int( dmask(i  ,j    ,k+k3d)+.1 )
      endif
      enddo
      enddo
      enddo


! Set up pointer arrays

      ixn = 0
      iyn = 0
      izn = 0

      do k = kl_bnd1,ku_bnd1
       do j = jl_bnd1,ju_bnd1
        do i = il_bnd1,iu_bnd1
        ix = ixnorm(i,j,k)
        iy = iynorm(i,j,k)
        iz = iznorm(i,j,k)
        if(learth(i,j,k)) then
          idzones(lb)              = idzones(lb)+1
          p_zone(1,idzones(lb),lb) = i
          p_zone(2,idzones(lb),lb) = j 
          p_zone(3,idzones(lb),lb) = k
c          write(*,*)'iz,jz,kz=',i,j,k
        endif
        if(lwall(i,j,k)) then
          idfaces(lb)              = idfaces(lb)+1
          p_face(1,idfaces(lb),lb) = i
          p_face(2,idfaces(lb),lb) = j
          p_face(3,idfaces(lb),lb) = k
          p_face(4,idfaces(lb),lb) = 
     .        abs(ix)*(1 - (ix-1)/2) +
     .        abs(iy)*(3 - (iy-1)/2) +
     .        abs(iz)*(5 - (iz-1)/2)
!
!
          if(ix.ne.0)ixn=ixn+1
          if(iy.ne.0)iyn=iyn+1
          if(iz.ne.0)izn=izn+1
c      write(*,*)'if,jf,kf=',i,j,k
        endif
      enddo
      enddo
      enddo

c      write(*,*)'idzones(lb)=',idzones(lb)
c      write(*,*)'idfaces(lb)=',idfaces(lb)
c      write(*,*)'ixn,iyn,izn=',ixn,iyn,izn
c      pause


      ncells         = -1
      ncells_buffer  = real(ncells)
      do k = kl_bnd1,ku_bnd1
        z0 = zmin+delz*(real(k-nguard)-.5)
        do j = jl_bnd1,ju_bnd1
          y0 = ymin+dely*(real(j-nguard)-.5)
          do i = il_bnd1,iu_bnd1
            x0 = xmin+delx*(real(i-nguard)-.5)


! Test size of modified earth box.
            xtest = (x0-(ebox_xmin-ncells_buffer*delx))
     .             *(x0-(ebox_xmax+ncells_buffer*delx))
            ytest = (y0-(ebox_ymin-ncells_buffer*dely))
     .             *(y0-(ebox_ymax+ncells_buffer*dely))
            ztest = -1.
            if(ndim.eq.3) ztest=
     .               (z0-(ebox_zmin-ncells_buffer*delz))*
     .               (z0-(ebox_zmax+ncells_buffer*delz))

            if(xtest.lt.0.and.ytest.lt.0.and.ztest.lt.0.) then
              idefields(lb)                 = idefields(lb) + 1
              p_efields(1,idefields(lb),lb) = i
              p_efields(2,idefields(lb),lb) = j
              p_efields(3,idefields(lb),lb) = k
            endif
          enddo
        enddo
      enddo

!
!-----now set up pointers for global blocks and determine
!     number of domain blocks and boundary blocks

      return
      end subroutine earth_cells
!---------------------------------------------------------------
!
!---------------------------------------------------------------
       subroutine interior_bc_fl(fx,fy,fz,vfx,vfy,vfz,v2,lb)
!=======================================================================
      use physcons
      use paramesh_dimensions
      use physicaldata
      use tree
      use strongfields
      use athena_specific_parameters
      implicit none
!--------------------------------------------------------------
! include files for amr


#include "amr_shmem.fh"
#include "shmem_reduce.fh"
#define face_x  real,dimension(ne,iu_bnd1p,ju_bnd1p,ku_bnd1p)
#define face_f  real,dimension(nf,iu_bnd1p,ju_bnd1p,ku_bnd1p)
#define cell_v  real,dimension(ne,iu_bnd1 ,ju_bnd1 ,ku_bnd1 )
      integer,intent(in  ) ::  lb
      face_f ,intent(inout) :: fx
      face_f ,intent(inout) :: fy
      face_f ,intent(inout) :: fz
      face_x ,intent(inout) ::vfx
      face_x ,intent(inout) ::vfy
      face_x ,intent(inout) ::vfz
      cell_v ,intent(in   ) :: v2

      integer i,j,k,n
!
      real :: bx,by,bz,bxt,byt,bzt,bs
!
!-----------------------------------------------------------------------
!
!-----Magnetic Field pointers for face based parameters
!
c      real, pointer :: bxf(:,:,:),byf(:,:,:),bzf(:,:,:)
!
!-----------------------------------------------------------------------
      integer :: ix,iy,iz,f,ii,jj,kk,if,l
      real    :: dt,time
      logical :: lflax,lflay
!---------------------------------------------------------------------------
!
      real :: rh, ph, eh
      real :: vx,vy,vz,bx1,by1,bz1,bxs,bys,bzs

c      bxf = facevarx1(1,:,:,:,1)
c      byf = facevary1(1,:,:,:,1)
c      bzf = facevarz1(1,:,:,:,1)

#ifdef LORENTZ_FORCE
       do if =1, ifaces
        i = iface(1,if)
        j = iface(2,if)
        k = iface(3,if)
        f = iface(4,if)
        if(f.eq.2)then  ! boundary face to your right
            ii = i+1
            fx(:5,ii,j,k) = 0.0
            fx( 2,ii,j,k) = v2(5,i,j,k)
            bxs           = bxf_sf(1,ii,j,k)
            bys           = bxf_sf(2,ii,j,k)
            bzs           = bxf_sf(3,ii,j,k)
            bx            = facevarx1(1,ii,j,k,1)
            by            = v2(7,i,j,k)
            bz            = v2(8,i,j,k)
            vx            = 0.0
            vy            = 0.
            vz            = 0.
            bxt           = bxs+bx
            byt           = bys+by
            bzt           = bzs+bz
            fx( 6,ii,j,k) = 0.
            fx( 7,ii,j,k) =-vy*bxt
            fx( 8,ii,j,k) =-vz*bxt
            bs            = (by*by+bz*bz-bx*bx)*pie8i
            fx( 9,ii,j,k) = bs+pie4i*(bys*by+bzs*bz-bxs*bx)
            fx(10,ii,j,k) =-(bxt*by+bx*bys)*pie4i
            fx(11,ii,j,k) =-(bxt*bz+bx*bzs)*pie4i
            vfx(2,ii,j,k) = 0.0
            vfx(3,ii,j,k) = 0.0
            vfx(4,ii,j,k) = 0.0
            vfx(6,ii,j,k) = 0.0
        elseif(f.eq.1)then   ! wall is on the left
            fx(:5,i,j,k)  = 0.0
            fx( 2,i,j,k)  = v2(5,i,j,k)
            bxs           = bxf_sf(1,i,j,k)
            bys           = bxf_sf(2,i,j,k)
            bzs           = bxf_sf(3,i,j,k)
            bx            = facevarx1(1,i,j,k,1)
            by            = v2(7,i,j,k)
            bz            = v2(8,i,j,k)
            vy            = 0.
            vz            = 0.
            vx            = 0.0
            bxt           = bxs+bx
            byt           = bys+by
            bzt           = bzs+bz
            bs            = (by*by+bz*bz-bx*bx)*pie8i
            fx( 6,i  ,j,k)= 0.0
            fx( 7,i  ,j,k)=-vy*bxt
            fx( 8,i  ,j,k)=-vz*bxt
            fx( 9,i  ,j,k) = bs+pie4i*(bys*by+bzs*bz-bxs*bx)
            fx(10,i  ,j,k) =-(bxt*by+bx*bys)*pie4i
            fx(11,i  ,j,k) =-(bxt*bz+bx*bzs)*pie4i
            vfx(2,i  ,j,k) = 0.0
            vfx(3,i  ,j,k) = 0.0
            vfx(4,i  ,j,k) = 0.0
            vfx(6,i  ,j,k) = 0.0
        elseif(f.eq.4)then
            jj = j+1
            fy(:5,i,jj,k) = 0.0
            fy( 3,i,jj,k) = v2(5,i,j,k)
            bxs           = byf_sf(3,i,jj,k)
            bys           = byf_sf(1,i,jj,k)
            bzs           = byf_sf(2,i,jj,k)
            bx            = v2(6,i,j,k)
            by            = facevary1(1,i,jj,k,1)
            bz            = v2(8,i,j,k)
            vx            = 0.
            vz            = 0.
            vy            = 0.0
            bxt           = bxs+bx
            byt           = bys+by
            bzt           = bzs+bz
            bs            = (bx*bx+bz*bz-by*by)*pie8i
            fy( 6,i,jj,k) =-vx*byt
            fy( 7,i,jj,k) = 0.0
            fy( 8,i,jj,k) =-vz*byt
            fy( 9,i,jj,k) =-(byt*bx+by*bxs)*pie4i
            fy(10,i,jj,k) = bs+pie4i*(bxs*bx+bzs*bz-bys*by)
            fy(11,i,jj,k) =-(byt*bz+by*bzs)*pie4i
            vfy(2,i,jj,k) = 0.0
            vfy(3,i,jj,k) = 0.0
            vfy(4,i,jj,k) = 0.0
            vfy(7,i,jj,k) = 0.0
        elseif(f.eq.3)then
            fy(:5,i,j,k)   = 0.0
            fy( 3,i,j,k)   = v2(5,i,j,k)
            bxs            = byf_sf(3,i,j,k)
            bys            = byf_sf(1,i,j,k)
            bzs            = byf_sf(2,i,j,k)
            bx             = v2(6,i,j,k)
            by             = facevary1(1,i,j,k,1)
            bz             = v2(8,i,j,k)
            vx             = 0.
            vz             = 0.
            vy             = 0.0
            bxt            = bxs+bx
            byt            = bys+by
            bzt            = bzs+bz
            bs             = (bx*bx+bz*bz-by*by)*pie8i
            fy( 6,i,j  ,k) =-vx*byt
            fy( 7,i,j  ,k) = 0.0
            fy( 8,i,j  ,k) =-vz*byt
            fy( 9,i,j  ,k) =-(byt*bx+by*bxs)*pie4i
            fy(10,i,j  ,k) = bs+pie4i*(bxs*bx+bzs*bz-bys*by)
            fy(11,i,j  ,k) =-(byt*bz+by*bzs)*pie4i
            vfy(2,i,j  ,k) = 0.0
            vfy(3,i,j  ,k) = 0.0
            vfy(4,i,j  ,k) = 0.0
            vfy(7,i,j  ,k) = 0.0
        elseif(f.eq.6)then
            kk             = k+1
            fz(:5,i,j,kk)  = 0.0
            fz( 4,i,j,kk)  = v2(5,i,j,k)
            bxs            = bzf_sf(2,i,j,kk)
            bys            = bzf_sf(3,i,j,kk)
            bzs            = bzf_sf(1,i,j,kk)
            bx             = v2(6,i,j,k)
            by             = v2(7,i,j,k)
            bz             = facevarz1(1,i,j,kk,1)
            vx             = 0.0
            vy             = 0.0
            vz             = 0.0
            bxt            = bxs+bx
            byt            = bys+by
            bzt            = bzs+bz
            bs             = (bx*bx+by*by-bz*bz)*pie8i

            fz( 6,i,j,kk) =-vx*bzt
            fz( 7,i,j,kk) =-vy*bzt
            fz( 8,i,j,kk) = 0.0
            fz( 9,i,j,kk) =-(bzt*bx+bz*bxs)*pie4i
            fz(10,i,j,kk) =-(bzt*by+bz*bys)*pie4i
            fz(11,i,j,kk) = bs+pie4i*(bxs*bx+bys*by-bzs*bz)
            vfz(2,i,j,kk) = 0.0
            vfz(3,i,j,kk) = 0.0
            vfz(4,i,j,kk) = 0.0
            vfz(8,i,j,kk) = 0.0
        elseif(f.eq.5)then
            fz(:5,i,j,k  ) = 0.0
            fz( 4,i,j,k  ) = v2(5,i,j,k)
            bxs            = bzf_sf(2,i,j,k)
            bys            = bzf_sf(3,i,j,k)
            bzs            = bzf_sf(1,i,j,k)
            bx             = v2(6,i,j,k)
            by             = v2(7,i,j,k)
            bz             = facevarz1(1,i,j,k,1)
            vx             = 0.
            vy             = 0.
            vz             = 0.0
            bxt            = bxs+bx
            byt            = bys+by
            bzt            = bzs+bz
            bs             = (bx*bx+by*by-bz*bz)*pie8i
            fz( 6,i,j,k  ) =-vx*bzt
            fz( 7,i,j,k  ) =-vy*bzt
            fz( 8,i,j,k  ) = 0.0
            fz( 9,i,j,k  ) =-(bzt*bx+bz*bxs)*pie4i
            fz(10,i,j,k  ) =-(bzt*by+bz*bys)*pie4i
            fz(11,i,j,k  ) = bs+pie4i*(bxs*bx+bys*by-bzs*bz)
            vfz(2,i,j,k  ) = 0.0
            vfz(3,i,j,k  ) = 0.0
            vfz(4,i,j,k  ) = 0.0
            vfz(8,i,j,k  ) = 0.0
        endif
       enddo
c      endif

#else
       do if =1, ifaces
        i = iface(1,if)
        j = iface(2,if)
        k = iface(3,if)
        f = iface(4,if)
        if(f.eq.2)then  ! boundary face to your right
            fx(:ne,i+1,j,k)= 0.0
            fx(  2,i+1,j,k)= v2(5,i,j,k)
            fx( 9,i+1,j,k) = 0.0
            fx(10,i+1,j,k) = 0.0
            fx(11,i+1,j,k) = 0.0

        elseif(f.eq.1)then   ! wall is on the left
            fx(:ne,i  ,j,k)= 0.0
            fx(  2,i  ,j,k)= v2(5,i,j,k)
            fx( 9,i  ,j,k) = 0.0
            fx(10,i  ,j,k) = 0.0
            fx(11,i  ,j,k) = 0.0

        elseif(f.eq.4)then
            fy(:ne,i,j+1,k)= 0.0
            fy(  3,i,j+1,k)= v2(5,i,j,k)
            fy(  9,i,j+1,k)= 0.0
            fy( 10,i,j+1,k)= 0.0
            fy( 11,i,j+1,k)= 0.0

        elseif(f.eq.3)then
            fy(:ne,i,j  ,k)= 0.0
            fy(  3,i,j  ,k)= v2(5,i,j,k)
            fy(  9,i,j  ,k)= 0.0
            fy( 10,i,j  ,k)= 0.0
            fy( 11,i,j  ,k)= 0.0

        elseif(f.eq.6)then
            fz(:ne,i,j,k+1)= 0.0
            fz(  4,i,j,k+1)= v2(5,i,j,k)
            fz( 9,i,j,k+1) = 0.0
            fz(10,i,j,k+1) = 0.0
            fz(11,i,j,k+1) = 0.0
        elseif(f.eq.5)then
            fz(:ne,i,j,k  )= 0.0
            fz( 4,i,j,k  ) = v2(5,i,j,k)
            fz( 9,i,j,k  ) = 0.0
            fz(10,i,j,k  ) = 0.0
            fz(11,i,j,k  ) = 0.0

        endif
       enddo


#endif /*LORENTZ_FORCE*/
      return
      end  subroutine interior_bc_fl
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      subroutine interior_bc_states_nplus_1half(vlx,vrx,idir,lb)
      use paramesh_dimensions
      use physicaldata
      use tree
      use athena_specific_parameters
      implicit none
#define cell_vp real,dimension(ne,iu_bnd1p,ju_bnd1p,ku_bnd1p)
      integer,intent(in   )::idir,lb
      cell_vp,intent(inout)::vlx
      cell_vp,intent(inout)::vrx

      integer, parameter :: imin = 1
      integer, parameter :: imax = iu_bnd1-1
      integer, parameter :: jmin = 1
      integer, parameter :: jmax = ju_bnd1-k2d
      integer, parameter :: kmin = 1
      integer, parameter :: kmax = ku_bnd1-k3d

      integer ::il,jl,kl,i,j,k,in
      integer :: f,ii,iif
c       write(*,*)'il_bnd1,nguard,k1d=',il_bnd1,nguard,k1d
c       write(*,*)'jl_bnd1,nguard,k2d=',jl_bnd1,nguard,k2d
c       write(*,*)'imin,imax,jmin,jmax=',imin,imax,jmin,jmax
      do iif =1, ifaces
        i = iface(1,iif)
        j = iface(2,iif)
        k = iface(3,iif)
        f = iface(4,iif)
        if(idir.eq.1.and.f.eq.2)then
c         if(i.le.imax.and.i.ge.imin)then
#ifdef NO_SLIP
            vrx(1,i+1,j,k) = vlx(1,i+1,j,k)
            vrx(2,i+1,j,k) = 0.0
            vrx(3,i+1,j,k) = 0.0
            vrx(4,i+1,j,k) = 0.0
            vlx(2,i+1,j,k) = 0.0
            vlx(3,i+1,j,k) = 0.0
            vlx(4,i+1,j,k) = 0.0

            vrx(5,i+1,j,k) = vlx(5,i+1,j,k)
            vrx(6,i+1,j,k) = 0.0
            vlx(6,i+1,j,k) = 0.0
            vrx(7,i+1,j,k) = vlx(7,i+1,j,k)
            vrx(8,i+1,j,k) = vlx(8,i+1,j,k)
#endif
#ifdef SLIP
            vrx(1,i+1,j,k) = vlx(1,i+1,j,k)
            vrx(2,i+1,j,k) = 0.0
            vlx(2,i+1,j,k) = 0.0
            vrx(3,i+1,j,k) = vlx(3,i+1,j,k)
            vrx(4,i+1,j,k) = vlx(4,i+1,j,k)
            vrx(5,i+1,j,k) = vlx(5,i+1,j,k)
            vrx(6,i+1,j,k) = 0.0
            vlx(6,i+1,j,k) = 0.0
            vrx(7,i+1,j,k) = vlx(7,i+1,j,k)
            vrx(8,i+1,j,k) = vlx(8,i+1,j,k)
#endif
c         endif
        endif
        if(idir.eq.1.and.f.eq.1)then   ! wall is on the left
c         if(i.le.imax.and.i.ge.imin)then
#ifdef NO_SLIP
            vlx(1,i,j,k  ) = vrx(1,i,j,k)
            vlx(2,i,j,k  ) = 0.0
            vlx(3,i,j,k  ) = 0.0
            vlx(4,i,j,k  ) = 0.0
            vrx(2,i,j,k  ) = 0.0
            vrx(3,i,j,k  ) = 0.0
            vrx(4,i,j,k  ) = 0.0
            vlx(5,i,j,k  ) = vrx(5,i,j,k)
            vlx(6,i,j,k  ) = 0.0
            vrx(6,i,j,k  ) = 0.0
            vlx(7,i,j,k  ) = vrx(7,i,j,k)
            vlx(8,i,j,k  ) = vrx(8,i,j,k)
#endif
#ifdef SLIP
            vlx(1,i,j,k  ) = vrx(1,i,j,k)
            vlx(2,i,j,k  ) = 0.0
            vrx(2,i,j,k  ) = 0.0
            vlx(3,i,j,k  ) = vrx(3,i,j,k)
            vlx(4,i,j,k  ) = vrx(4,i,j,k)
            vlx(5,i,j,k  ) = vrx(5,i,j,k)
            vlx(6,i,j,k  ) = 0.0
            vrx(6,i,j,k  ) = 0.0
            vlx(7,i,j,k  ) = vrx(7,i,j,k)
            vlx(8,i,j,k  ) = vrx(8,i,j,k)
#endif
c         endif
        endif
        if(idir.eq.2.and.f.eq.4)then
         if(j.le.jmax.and.j.ge.jmin)then
#ifdef NO_SLIP
            vrx(1,i,j+1,k) = vlx(1,i,j+1,k)
            vrx(2,i,j+1,k) = 0.0
            vrx(3,i,j+1,k) = 0.0
            vrx(4,i,j+1,k) = 0.0
            vlx(2,i,j+1,k) = 0.0
            vlx(3,i,j+1,k) = 0.0
            vlx(4,i,j+1,k) = 0.0
            vrx(5,i,j+1,k) = vlx(5,i,j+1,k)
            vrx(6,i,j+1,k) = 0.0
            vlx(6,i,j+1,k) = 0.0
            vrx(7,i,j+1,k) = vlx(7,i,j+1,k)
            vrx(8,i,j+1,k) = vlx(8,i,j+1,k)
#endif
#ifdef SLIP
            vrx(1,i,j+1,k) = vlx(1,i,j+1,k)
            vrx(2,i,j+1,k) = 0.0
            vlx(2,i,j+1,k) = 0.0
            vrx(3,i,j+1,k) = vlx(3,i,j+1,k)
            vrx(4,i,j+1,k) = vlx(4,i,j+1,k)
            vrx(5,i,j+1,k) = vlx(5,i,j+1,k)
            vrx(6,i,j+1,k) = 0.0
            vlx(6,i,j+1,k) = 0.0
            vrx(7,i,j+1,k) = vlx(7,i,j+1,k)
            vrx(8,i,j+1,k) = vlx(8,i,j+1,k)
#endif
c         endif
        endif
        if(idir.eq.2.and.f.eq.3)then
c         if(j.le.jmax.and.j.ge.jmin)then
#ifdef NO_SLIP
            vlx(1,i,j,k  ) = vrx(1,i,j,k)
            vlx(2,i,j,k  ) = 0.0
            vlx(3,i,j,k  ) = 0.0
            vlx(4,i,j,k  ) = 0.0
            vrx(2,i,j,k  ) = 0.0
            vrx(3,i,j,k  ) = 0.0
            vrx(4,i,j,k  ) = 0.0
            vlx(5,i,j,k  ) = vrx(5,i,j,k)
            vlx(6,i,j,k  ) = 0.0
            vrx(6,i,j,k  ) = 0.0
            vlx(7,i,j,k  ) = vrx(7,i,j,k)
            vlx(8,i,j,k  ) = vrx(8,i,j,k)
#endif
#ifdef SLIP
            vlx(1,i,j,k  ) = vrx(1,i,j,k)
            vlx(2,i,j,k  ) = 0.0
            vrx(2,i,j,k  ) = 0.0
            vlx(3,i,j,k  ) = vrx(3,i,j,k)
            vlx(4,i,j,k  ) = vrx(4,i,j,k)
            vlx(5,i,j,k  ) = vrx(5,i,j,k)
            vlx(6,i,j,k  ) = 0.0
            vrx(6,i,j,k  ) = 0.0
            vlx(7,i,j,k  ) = vrx(7,i,j,k)
            vlx(8,i,j,k  ) = vrx(8,i,j,k)
#endif
         endif
        endif
        if(idir.eq.3.and.f.eq.6)then
c         if(k.le.kmax.and.k.ge.kmin)then
#ifdef NO_SLIP
            vrx(1,i,j,k+1) = vlx(1,i,j,k+1)
            vrx(2,i,j,k+1) = 0.0
            vrx(3,i,j,k+1) = 0.0
            vrx(4,i,j,k+1) = 0.0
            vlx(2,i,j,k+1) = 0.0
            vlx(3,i,j,k+1) = 0.0
            vlx(4,i,j,k+1) = 0.0
            vrx(5,i,j,k+1) = vlx(5,i,j,k+1)
            vrx(6,i,j,k+1) = 0.0
            vlx(6,i,j,k+1) = 0.0
            vrx(7,i,j,k+1) = vlx(7,i,j,k+1)
            vrx(8,i,j,k+1) = vlx(8,i,j,k+1)
#endif
#ifdef SLIP
            vrx(1,i,j,k+1) = vlx(1,i,j,k+1)
            vrx(2,i,j,k+1) = 0.0
            vlx(2,i,j,k+1) = 0.0
            vrx(3,i,j,k+1) = vlx(3,i,j,k+1)
            vrx(4,i,j,k+1) = vlx(4,i,j,k+1)
            vrx(5,i,j,k+1) = vlx(5,i,j,k+1)
            vrx(6,i,j,k+1) = 0.0
            vlx(6,i,j,k+1) = 0.0
            vrx(7,i,j,k+1) = vlx(7,i,j,k+1)
            vrx(8,i,j,k+1) = vlx(8,i,j,k+1)
#endif
c         endif
        endif
        if(idir.eq.3.and.f.eq.5)then
c         if(k.le.kmax.and.k.ge.kmin)then
#ifdef NO_SLIP
            vlx(1,i,j,k  ) = vrx(1,i,j,k)
            vlx(2,i,j,k  ) = 0.0
            vlx(3,i,j,k  ) = 0.0
            vlx(4,i,j,k  ) = 0.0
            vrx(2,i,j,k  ) = 0.0
            vrx(3,i,j,k  ) = 0.0
            vrx(4,i,j,k  ) = 0.0
            vlx(5,i,j,k  ) = vrx(5,i,j,k)
            vlx(6,i,j,k  ) = 0.0
            vrx(6,i,j,k  ) = 0.0
            vlx(7,i,j,k  ) = vrx(7,i,j,k)
            vlx(8,i,j,k  ) = vrx(8,i,j,k)
#endif
#ifdef SLIP
            vlx(1,i,j,k  ) = vrx(1,i,j,k)
            vlx(2,i,j,k  ) = 0.0
            vrx(2,i,j,k  ) = 0.0
            vlx(3,i,j,k  ) = vrx(3,i,j,k)
            vlx(4,i,j,k  ) = vrx(4,i,j,k)
            vlx(5,i,j,k  ) = vrx(5,i,j,k)
            vlx(6,i,j,k  ) = 0.0
            vrx(6,i,j,k  ) = 0.0
            vlx(7,i,j,k  ) = vrx(7,i,j,k)
            vlx(8,i,j,k  ) = vrx(8,i,j,k)
#endif
c         endif
        endif
      enddo

      return
      end subroutine interior_bc_states_nplus_1half
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
      subroutine interior_bc_states_n(v,idir,lb)
      use paramesh_dimensions
      use physicaldata
      use tree
      use athena_specific_parameters
      implicit none
      integer,intent(in   )::idir,lb
      cell_v ,intent(inout)::v
      integer, parameter :: imin = 1
      integer, parameter :: imax = iu_bnd1-1
      integer, parameter :: jmin = 1
      integer, parameter :: jmax = ju_bnd1-k2d
      integer, parameter :: kmin = 1
      integer, parameter :: kmax = ku_bnd1-k3d

      integer ::il,jl,kl,i,j,k,in
      integer :: f,ii,iif
c       write(*,*)'il_bnd1,nguard,k1d=',il_bnd1,nguard,k1d
c       write(*,*)'jl_bnd1,nguard,k2d=',jl_bnd1,nguard,k2d
c       write(*,*)'imin,imax,jmin,jmax=',imin,imax,jmin,jmax
      do iif =1, ifaces
        i = iface(1,iif)
        j = iface(2,iif)
        k = iface(3,iif)
        f = iface(4,iif)
        if(idir.eq.1.and.f.eq.2)then
c         if(i.le.imax.and.i.ge.imin)then
#ifdef NO_SLIP
            v(1,i+1,j,k) = v(1,i,j,k)
            v(2,i+1,j,k) =-v(2,i,j,k)
            v(3,i+1,j,k) =-v(3,i,j,k)
            v(4,i+1,j,k) =-v(4,i,j,k)
            v(5,i+1,j,k) = v(5,i,j,k)
            v(6,i+1,j,k) =-v(6,i,j,k)
            v(7,i+1,j,k) = v(7,i,j,k)
            v(8,i+1,j,k) = v(8,i,j,k)
#endif
#ifdef SLIP
            v(1,i+1,j,k) = v(1,i,j,k)
            v(2,i+1,j,k) =-v(2,i,j,k)
            v(3,i+1,j,k) = v(3,i,j,k)
            v(4,i+1,j,k) = v(4,i,j,k)
            v(5,i+1,j,k) = v(5,i,j,k)
            v(6,i+1,j,k) =-v(6,i,j,k)
            v(7,i+1,j,k) = v(7,i,j,k)
            v(8,i+1,j,k) = v(8,i,j,k)
#endif
c         endif
        endif
        if(idir.eq.1.and.f.eq.1)then   ! wall is on the left
c         if(i.le.imax.and.i.ge.imin)then
#ifdef NO_SLIP
            v(1,i-1,j,k) = v(1,i,j,k)
            v(2,i-1,j,k) =-v(2,i,j,k)
            v(3,i-1,j,k) =-v(3,i,j,k)
            v(4,i-1,j,k) =-v(4,i,j,k)
            v(5,i-1,j,k) = v(5,i,j,k)
            v(6,i-1,j,k) =-v(6,i,j,k)
            v(7,i-1,j,k) = v(7,i,j,k)
            v(8,i-1,j,k) = v(8,i,j,k)
#endif
#ifdef SLIP
            v(1,i-1,j,k) = v(1,i,j,k)
            v(2,i-1,j,k) =-v(2,i,j,k)
            v(3,i-1,j,k) = v(3,i,j,k)
            v(4,i-1,j,k) = v(4,i,j,k)
            v(5,i-1,j,k) = v(5,i,j,k)
            v(6,i-1,j,k) =-v(6,i,j,k)
            v(7,i-1,j,k) = v(7,i,j,k)
            v(8,i-1,j,k) = v(8,i,j,k)
#endif
c         endif
        endif
        if(idir.eq.2.and.f.eq.4)then
c         if(j.le.jmax.and.j.ge.jmin)then
#ifdef NO_SLIP
            v(1,i,j+1,k) = v(1,i,j,k)
            v(2,i,j+1,k) =-v(2,i,j,k)
            v(3,i,j+1,k) =-v(3,i,j,k)
            v(4,i,j+1,k) =-v(4,i,j,k)
            v(5,i,j+1,k) = v(5,i,j,k)
            v(6,i,j+1,k) = v(6,i,j,k)
            v(7,i,j+1,k) =-v(7,i,j,k)
            v(8,i,j+1,k) = v(8,i,j,k)
#endif
#ifdef SLIP
            v(1,i,j+1,k) = v(1,i,j,k)
            v(2,i,j+1,k) = v(2,i,j,k)
            v(3,i,j+1,k) =-v(3,i,j,k)
            v(4,i,j+1,k) = v(4,i,j,k)
            v(5,i,j+1,k) = v(5,i,j,k)
            v(6,i,j+1,k) = v(6,i,j,k)
            v(7,i,j+1,k) =-v(7,i,j,k)
            v(8,i,j+1,k) = v(8,i,j,k)
#endif
c         endif
        endif
        if(idir.eq.2.and.f.eq.3)then
c         if(j.le.jmax.and.j.ge.jmin)then
#ifdef NO_SLIP
            v(1,i,j-1,k) = v(1,i,j,k)
            v(2,i,j-1,k) =-v(2,i,j,k)
            v(3,i,j-1,k) =-v(3,i,j,k)
            v(4,i,j-1,k) =-v(4,i,j,k)
            v(5,i,j-1,k) = v(5,i,j,k)
            v(6,i,j-1,k) = v(6,i,j,k)
            v(7,i,j-1,k) =-v(7,i,j,k)
            v(8,i,j-1,k) = v(8,i,j,k)
#endif
#ifdef SLIP
            v(1,i,j-1,k) = v(1,i,j,k)
            v(2,i,j-1,k) = v(2,i,j,k)
            v(3,i,j-1,k) =-v(3,i,j,k)
            v(4,i,j-1,k) = v(4,i,j,k)
            v(5,i,j-1,k) = v(5,i,j,k)
            v(6,i,j-1,k) = v(6,i,j,k)
            v(7,i,j-1,k) =-v(7,i,j,k)
            v(8,i,j-1,k) = v(8,i,j,k)
#endif
c         endif
        endif
        if(idir.eq.3.and.f.eq.6)then
c         if(k.le.kmax.and.k.ge.kmin)then
#ifdef NO_SLIP
            v(1,i,j,k+1) = v(1,i,j,k)
            v(2,i,j,k+1) =-v(2,i,j,k)
            v(3,i,j,k+1) =-v(3,i,j,k)
            v(4,i,j,k+1) =-v(4,i,j,k)
            v(5,i,j,k+1) = v(5,i,j,k)
            v(6,i,j,k+1) = v(6,i,j,k)
            v(7,i,j,k+1) = v(7,i,j,k)
            v(8,i,j,k+1) =-v(8,i,j,k)
#endif
#ifdef SLIP
            v(1,i,j,k+1) = v(1,i,j,k)
            v(2,i,j,k+1) = v(2,i,j,k)
            v(3,i,j,k+1) = v(3,i,j,k)
            v(4,i,j,k+1) =-v(4,i,j,k)
            v(5,i,j,k+1) = v(5,i,j,k)
            v(6,i,j,k+1) = v(6,i,j,k)
            v(7,i,j,k+1) = v(7,i,j,k)
            v(8,i,j,k+1) =-v(8,i,j,k)
#endif
c         endif
        endif
        if(idir.eq.3.and.f.eq.5)then
c         if(k.le.kmax.and.k.ge.kmin)then
#ifdef NO_SLIP
            v(1,i,j,k-1) = v(1,i,j,k)
            v(2,i,j,k-1) =-v(2,i,j,k)
            v(3,i,j,k-1) =-v(3,i,j,k)
            v(4,i,j,k-1) =-v(4,i,j,k)
            v(5,i,j,k-1) = v(5,i,j,k)
            v(6,i,j,k-1) = v(6,i,j,k)
            v(7,i,j,k-1) = v(7,i,j,k)
            v(8,i,j,k-1) =-v(8,i,j,k)
#endif
#ifdef SLIP
            v(1,i,j,k-1) = v(1,i,j,k)
            v(2,i,j,k-1) = v(2,i,j,k)
            v(3,i,j,k-1) = v(3,i,j,k)
            v(4,i,j,k-1) =-v(4,i,j,k)
            v(5,i,j,k-1) = v(5,i,j,k)
            v(6,i,j,k-1) = v(6,i,j,k)
            v(7,i,j,k-1) = v(7,i,j,k)
            v(8,i,j,k-1) =-v(8,i,j,k)
#endif
c         endif
        endif
      enddo

      return
      end subroutine interior_bc_states_n
!---------------------------------------------------------------
!
!---------------------------------------------------------------
      subroutine internal_hardwall_bc(i,j,k,idir)
      use paramesh_dimensions
      use physicaldata
      use physcons
      use grid
      use athena_specific_parameters
      implicit none
      integer,intent(in  ) :: i,j,k,idir
      integer :: ii,m,jj,kk,lb
      integer, dimension(1000) :: unsorted
      integer, dimension(1000) :: sorted
      select case (idir) 
      case (1) ! x-direction

      ifaces_1d = 0
      do ii =1, ifaces
       if(iface(2,ii).eq.j.and.iface(3,ii).eq.k ) then	! boundary on cell
        m=iface(4,ii)
        if(m.eq.1 .or. m.eq.2) then ! and in current direction
         ifaces_1d             = ifaces_1d+1
         iface_1d(1,ifaces_1d) = iface(1,ii)
         iface_1d(2,ifaces_1d) = m
        endif
       endif
      enddo
	    
      izones_1d=0
      do ii =1, izones
       if(izone(2,ii).eq.j .and. izone(3,ii).eq.k) then
        izones_1d           = izones_1d+1
        izone_1d(izones_1d) = izone(1,ii)
       endif
      enddo

      case (2) ! y-direction

      ifaces_1d=0
      do ii =1, ifaces
       if(iface(3,ii).eq.k .and. iface(1,ii).eq.i) then  ! boundary on cell
        m=iface(4,ii)
        if(m.eq.3 .or. m.eq.4) then                    ! and in current direction
         ifaces_1d             = ifaces_1d+1
         iface_1d(1,ifaces_1d) = iface(2,ii)
         iface_1d(2,ifaces_1d) = m-2
        endif
       endif
      enddo

      izones_1d=0
      do ii =1, izones
       if(izone(3,ii).eq.k .and. izone(1,ii).eq.i) then
        izones_1d           = izones_1d+1
        izone_1d(izones_1d) = izone(2,ii)
       endif
      enddo

      case (3) ! z-direction

      ifaces_1d=0
      do ii =1, ifaces
       if(iface(2,ii).eq.j .and. iface(1,ii).eq.i) then  ! boundary on cell
        m=iface(4,ii)
        if(m.eq.5 .or. m.eq.6) then                    ! and in current direction
         ifaces_1d             = ifaces_1d+1
         iface_1d(1,ifaces_1d) = iface(3,ii)
         iface_1d(2,ifaces_1d) = m-4
        endif
       endif
      enddo

      izones_1d=0
      do ii =1, izones
       if(izone(2,ii).eq.k .and. izone(1,ii).eq.i) then
        izones_1d           = izones_1d+1
        izone_1d(izones_1d) = izone(3,ii)
       endif
      enddo

      end select

      if(izones_1d>1)then
       unsorted(:izones_1d) = izone_1d(:izones_1d)
       call bubble_sort(unsorted,sorted,izones_1d)
       izone_1d(:izones_1d) = sorted(:izones_1d)
       izones_min = izone_1d(1        )+1
       izones_max = izone_1d(izones_1d)
      endif 

      return
      end subroutine internal_hardwall_bc
!---------------------------------------------------------------
!
!---------------------------------------------------------------
      subroutine bubble_sort(unsorted,sorted,n)
      use paramesh_dimensions
      use physicaldata
      integer :: i,j,n,temporary,nmax
      integer, dimension(1000),intent(in ) :: unsorted
      integer, dimension(1000),intent(out) :: sorted
      sorted = unsorted
      do i = 1, n-1,1
       do j = i+1,n,1
        if(sorted(i)>sorted(j))then
         temporary = sorted(j)
         sorted(j) = sorted(i)
         sorted(i) = temporary
        endif
       enddo
      enddo
      return 
      end subroutine bubble_sort
!---------------------------------------------------------------
!
!---------------------------------------------------------------

      end module interior_gbc
