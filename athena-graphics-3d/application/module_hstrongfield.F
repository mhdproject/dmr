#include "paramesh_preprocessor.fh"
#include "athena_preprocessor.fh"
#define NEVER30
      module strongfields
      use paramesh_dimensions
      implicit none 
      private
      save
      public strongfields_wb
      public rotate_dipole
!
!-----Strong Field Arrays
!
      public bxf_sf,byf_sf,bzf_sf,b_sf
      real, dimension(3,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,
     &                  kl_bnd1:ku_bnd1     )::bxf_sf
      real, dimension(3,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,
     &                  kl_bnd1:ku_bnd1     )::byf_sf
      real, dimension(3,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,
     &                  kl_bnd1:ku_bnd1+k3d )::bzf_sf
      real, dimension(3,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1
     &                 ,kl_bnd1:ku_bnd1     )::  b_sf
      real,public :: xsym,ysym,zsym
      real,public :: adipol
      real,public :: adipolx,adipoly,adipolz,xx0,yy0,zz0
      real,public :: adipolxt,adipolyt,adipolzt
      contains
!================================================================
!
!================================================================
      subroutine strongfields_wb(pe,lb)
      use physcons
      use paramesh_dimensions
      use physicaldata

! include file defining the tree
      use tree
c      use strongfields
!-----------------------------------------------------------------------
!
! include files for amr


#include "amr_shmem.fh"
#include "shmem_reduce.fh"
      integer :: pe,lb
      real, dimension(  il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,
     &                  kl_bnd1:ku_bnd1+k3d) ::ax,ay,az
      real, parameter  :: tiny=(1e-1)**2
      real :: xmin,xmax,ymin,ymax,zmin,zmax
      real :: bx1,by1,bz1,bxa,bya,bza
      real :: delx,dely,delz
      real :: r2,r1i,r2i,r3i,rdi,temp,x1,y1,z1,rdx,rdy,rdz
      real :: time
      common/evolution/time
      integer, parameter:: ng0= nguard*npgs
      real :: bxmax,bxmin,bymax,bymin,bzmax,bzmin
c      real, dimension(        iu_bnd1+1)   :: xpp,ypp,zpp
c      real :: xpos,ypos,zpos
      integer :: nx, nxp, ny, nyp, nz, nzp
      real, save :: ccbsize(3),bbox(2,3),ccoord(3) 
      integer :: icounter
      integer :: ilb, pe0, lb0
      integer :: ix,iy,iz

#ifdef MPI_USED
      include "mpif.h"
      integer,external :: shmem_my_pe
      integer          :: mype,ierrorcode,ierr
#endif /* MPI_USED */

#ifdef MPI_USED
      mype = shmem_my_pe()

!     if(pe.ne.mype) then
!       write(*,*) 'Error: strongfields_wb : calling shmem',
!    .             ' for a non-local block'
!       call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
!     endif
#endif /* MPI_USED */

      icounter = 0
!
!--------coordinates a working block 
!
!
c      write(*,*)'entering strongfields_wb:1'
      nx  = nxb+2*nguard        !<<< EDITED
      nxp = nx+1
      ny  = nyb+2*nguard*k2d         !<<< EDITED
      nyp = ny+k2d
      nz  = nzb+2*nguard*k3d
      nzp = nz+k3d

      pe0=pe; lb0=lb;
#ifdef MPI_USED
      if(pe.ne.mype) then
!write(*,*)'--------entering finding ilb'
         lb0=0;  ilb=strt_buffer;

         do while( (ilb.ge.strt_buffer.and.
     &              ilb.le.last_buffer).and.
     &              lb0.eq.0) 
            if(lb.eq.laddress(1,ilb).and.pe.eq.laddress(2,ilb)) then
                lb0=ilb
                pe0=mype
            endif
            ilb=ilb+1
         end do

         if(lb0.eq.0) then
            write(*,*) 'Error: strongfields_wb : could not find block',
     .             'and pe in laddress ------ stopping'
            call amr_close
            stop
         endif
!write(*,*)'--------exiting finding ilb'

      endif
#endif
      call shmem_real_get(bbox(1,1),bnd_box(1,1,lb0),6,mype)
      call shmem_real_get(ccoord(1),coord(1,lb0),3,mype)



      ccbsize(:) = bbox(2,:)-bbox(1,:)
      xmax  = ccoord(1) + 0.5*ccbsize(1)
      xmin  = ccoord(1) - 0.5*ccbsize(1)
      xmin  = bbox(1,1)
      xmax  = bbox(2,1)
c      xmax  = bnd_box(1,1,lb)
c      xmin  = bnd_box(2,1,lb)
c      write(*,*)'xmin,xmax,nxb=',xmin,xmax,nxb
      delx  = (xmax - xmin)/float(nxb)
      rdx   = 1./delx

      ymax  = ccoord(2) + 0.5*ccbsize(2)
      ymin  = ccoord(2) - 0.5*ccbsize(2)
      ymin  = bbox(1,2)
      ymax  = bbox(2,2)
c      write(*,*)'ymin,ymax,nyb=',ymin,ymax,nyb
      dely  = (ymax - ymin)/float(nyb)
      rdy   = 1./dely
#if N_DIM == 3
c      zmax  = ccoord(3) + 0.5*ccbsize(3)
c      zmin  = ccoord(3) - 0.5*ccbsize(3)
c      zmax  = bnd_box(1,3,lb)
c      zmin  = bnd_box(2,3,lb)
c      write(*,*)'zmin,zmax,nzb=',zmin,zmax,nzb
      zmin  = bbox(1,3)
      zmax  = bbox(2,3)
c      write(*,*)'zmin,zmax,nzb=',zmin,zmax,nzb
c      pause
      delz  = (zmax - zmin)/float(nzb)
      rdz   = 1./delz
#else
      zmax  = 0.
      zmin  = 0.
      delz  = 0.
      rdz   = 0.                ! The generalized inverse.
#endif /*NDIM3*/
c      write(*,*)'entering strongfields_wb:5'

#ifdef DEBUGG
      write(*,*)'after new calculation of xmax,xmin'
      write(*,*)'xmax,xmin=',xmax,xmin
      write(*,*)'ymax,ymin=',ymax,ymin
      write(*,*)'zmax,zmin=',zmax,zmin
#endif /*DEBUGG*/
!
!--------coordinates of block (not working block)
!



#ifdef MAGNETOSPHERE
c      write(*,*)'before a0_magnetosphere_strongfield call'
c      write(*,*)'xmin,delx  =',xmin,delx
c      if(delx .lt. 0.0)pause
c             icounter = icounter+1

       call rotate_dipole (adipolx ,adipoly ,adipolz ,0.0 )
       call rotate_dipole (adipolxt,adipolyt,adipolzt,time)

       call a0_magnetosphere_strongfield( ax,ay,az,pe,lb )
c      write(*,*)'after a0 call'
c      write(*,*)'xmin,delx  =',xmin,delx
c      write(*,*)'ax(5,1,1)=',ax(5,1,1)
c      write(*,*)'ay(5,1,1)=',ay(5,1,1)
c      write(*,*)'az(5,1,1)=',az(5,1,1)

#ifdef DEBUGGGG
       nzc = (ku_bnd1+k3d)/2
       write(*,*)'nzc=',nzc
       write(*,*)'ax '
      write(*,*)((ax(ix,iy,nzc),ix=il_bnd1,iu_bnd1+1  ),
     &                           iy=jl_bnd1,ju_bnd1+k2d)
       write(*,*)'ay '
      write(*,*)((ay(ix,iy,nzc),ix=il_bnd1,iu_bnd1+1),
     &                            iy=jl_bnd1,ju_bnd1+k2d)
       write(*,*)'az '
      write(*,*)((az(ix,iy,nzc),ix=il_bnd1,iu_bnd1+1),
     &                            iy=jl_bnd1,ju_bnd1+k2d)
#endif /*DEBUGG*/

! Compute the line integrated vector potentials
! from the vector potential on points.
      do iz=kl_bnd1,ku_bnd1+k3d
       do iy=jl_bnd1,ju_bnd1+k2d
         do ix=il_bnd1,iu_bnd1
           ax(ix,iy,iz) = 0.5*(ax(ix,iy,iz)+ax(ix+1 ,iy,iz)) 
         enddo
       enddo
      enddo
      do iz=kl_bnd1,ku_bnd1+k3d
       do iy=jl_bnd1,ju_bnd1
        do ix=il_bnd1,iu_bnd1+1
         ay(ix,iy,iz) = 0.5*(ay(ix,iy,iz)+ay(ix,iy+k2d,iz)) 
        enddo
       enddo
      enddo
      do iz=kl_bnd1,ku_bnd1
       do iy=jl_bnd1,ju_bnd1+k2d
        do ix=il_bnd1,iu_bnd1+1
         az(ix,iy,iz) = 0.5*(az(ix,iy,iz)+az(ix,iy,iz+k3d)) 
        enddo
       enddo
      enddo
#ifdef DEBUGG
      write(*,*)'ax(5,1,1)=',ax(5,1,1)
      write(*,*)'ay(5,1,1)=',ay(5,1,1)
      write(*,*)'az(5,1,1)=',az(5,1,1)
#endif
!
!------Overwrite magnetic fields at faces from vector potential on edges.
!
!      B = curl(A)
!
!      Bx = Az,y - Ay,z
!      By = Ax,z - Az,x
!      Bz = Ay,x - Ax,y
!
            do iz= kl_bnd1,ku_bnd1
             do iy= jl_bnd1,ju_bnd1
              do ix= il_bnd1,iu_bnd1+1
               bxf_sf(1,ix,iy,iz) = 
     &           ( az(ix,iy+k2d,iz    ) 
     &           - az(ix,iy    ,iz    ))*rdy
     &           -(ay(ix,iy    ,iz+k3d) 
     &           - ay(ix,iy    ,iz    ))*rdz
              enddo
             enddo
            enddo
!pmn
            if(neigh(1,1,lb).eq.-30) then
c              bxf_sf(1,1:nguard,:,:) = 0.
            endif
!pmn
            do iz= kl_bnd1,ku_bnd1
             do iy= jl_bnd1,ju_bnd1+k2d
              do ix= il_bnd1,iu_bnd1
               byf_sf(1,ix,iy,iz) =
     &           ( ax(ix    ,iy,iz+k3d) 
     &           - ax(ix    ,iy,iz    ))*rdz
     &           -(az(ix+1,iy,iz    )
     &           - az(ix    ,iy,iz    ))*rdx
            enddo
           enddo
          enddo
!pmn
            if(neigh(1,1,lb).eq.-30) then
c              byf_sf(1,1:nguard,:,:) = 0.
            endif
!pmn
          do iz= kl_bnd1,ku_bnd1+k3d
           do iy= jl_bnd1,ju_bnd1
            do ix= il_bnd1,iu_bnd1
             bzf_sf(1,ix,iy,iz) =
     &         ( ay(ix+1,iy    ,iz)
     &         - ay(ix    ,iy    ,iz))*rdx
     &         -(ax(ix    ,iy+k2d,iz)
     &         - ax(ix    ,iy    ,iz))*rdy
            enddo
           enddo
          enddo
!pmn
            if(neigh(1,1,lb).eq.-30) then
c              bzf_sf(1,1:nguard,:,:) = 0.
            endif
!pmn
#ifdef DEBUGG
           bxmax = maxval( bxf_sf(1,:,:,:) )
           bxmin = minval( bxf_sf(1,:,:,:) )
           bymax = maxval( byf_sf(1,:,:,:) )
           bymin = minval( byf_sf(1,:,:,:) )
           bzmax = maxval( bzf_sf(1,:,:,:) )
           bzmin = minval( bzf_sf(1,:,:,:) )
           write(*,*)'within strongfields_wb: 1' 
           write(*,*)'bxmax,bxmin = ',bxmax,bxmin
           write(*,*)'bymax,bymin = ',bymax,bymin
           write(*,*)'bzmax,bzmin = ',bzmax,bzmin
#endif /*DEBUGG*/
!
!-----Now compute face strong fields and then average the
!     result to the cell centers
c      write(*,*)'xmin,delx  =',xmin,delx
c      write(*,*)'entering strongfields_wb:2'

!     Loop over the x-faces of the cells.
      do iz=1,nz
         do iy=1,ny
            do ix=1,nxp
               
!  compute the coordinates of the lower-left corner                     

c               x1 = xmin + delx*(ix-1  )
c               y1 = ymin + dely*(iy-0.5)*k2d
c               z1 = zmin + delz*(iz-0.5)*k3d
               x1 = xmin + delx*(ix-nguard-1  )
               y1 = ymin + dely*(iy-nguard-0.5)*k2d
c              if(iy.eq.65)then
c               write(*,*)'ix,iy,x1,y1=',ix,iy,x1,y1
c              endif
#if N_DIM == 3
               z1 = zmin + delz*(iz-nguard-0.5)*k3d
#else
               z1 = 0.0
#endif               
               r2  = x1*x1 + y1*y1 + z1*z1
!  There is roundoff error in determining the coordinates
!  of the grid points from the coord() and bsize() arrays.
!  We use the SVD like trick of setting the inverse to zero.
               if( r2 .gt. tiny) then
                  r1i = 1./sqrt(r2)
               else
                  r1i = 0.
               endif
               r2i = r1i*r1i
               r3i = r2i*r1i
#if N_DIM < 3
               rdi = r2i*r2i
               z1  = 0.
#else /*NDIM3*/
               rdi = r3i*r2i
#endif /*NDIM3*/
               
               temp = adipolx*x1+adipoly*y1+adipolz*z1
               bx1 = (3.*x1*temp - r2*adipolx)*rdi
               by1 = (3.*y1*temp - r2*adipoly)*rdi
               bz1 = (3.*z1*temp - r2*adipolz)*rdi
               
               bxa = bx1
               bya = by1
               bza = bz1

#ifdef DEBUGG
               if(ix.eq.5.and.iy.eq.1)then
                write(*,*)'before image'
                write(*,*)'xmin,delx  =',xmin,delx
                write(*,*)'x1,y1,z1   =',x1,y1,z1
                write(*,*)'bxa,bya,bza=',bxa,bya,bza
               endif
#endif /*DEBUGG*/               
#ifdef IMAGE_DIPOLE
c        Now add symmetry dipole of same sign, but located at 
c        (2.0*xsym,0,0) instead of (0,0,0).

               x1 = x1 - 2.*xsym
               y1 = y1 - 2.*ysym
               z1 = z1 - 2.*zsym
               
               r2  = x1*x1 + y1*y1 + z1*z1
               if( r2 .gt. tiny ) then
                  r1i = 1./sqrt(r2)
               else
                  r1i = 0.
               endif
               r2i = r1i*r1i
               r3i = r2i*r1i
#if N_DIM < 3
               rdi = r2i*r2i
               z1 = 0.
#else /*NDIM3*/
               rdi = r3i*r2i
#endif /*NDIM3*/
               
               temp = adipolx*x1+adipoly*y1+adipolz*z1
               bx1 = (3.*x1*temp - r2*adipolx)*rdi
               by1 = (3.*y1*temp - r2*adipoly)*rdi
               bz1 = (3.*z1*temp - r2*adipolz)*rdi
               
! reverse direction of image dipole
!               bx1 = -bx1
!               by1 = -by1
!               bz1 = -bz1
               
               bxa = bxa+bx1
               bya = bya+by1
               bza = bza+bz1
#ifdef DEBUGG
               if(ix.eq.5.and.iy.eq.1)then
                write(*,*)'after image'
                write(*,*)'xsym,ysym  =',xsym,ysym
                write(*,*)'x1,y1,z1   =',x1,y1,z1
                write(*,*)'bxa,bya,bza=',bxa,bya,bza
               endif
#endif /*DEBUGG*/               

#endif /*IMAGE_DIPOLE*/
#if N_DIM < 3
               bza = 0.0
#endif
c               bxf_sf(1,ix,iy,iz) = bxa
               bxf_sf(2,ix,iy,iz) = bya
               bxf_sf(3,ix,iy,iz) = bza
        enddo               ! do ix
       enddo                  ! do iy
      enddo                     ! do iz
#ifdef DEBUGG
           bxmax = maxval( bxf_sf(1,:,:,:) )
           bxmin = minval( bxf_sf(1,:,:,:) )
           bymax = maxval( byf_sf(1,:,:,:) )
           bymin = minval( byf_sf(1,:,:,:) )
           bzmax = maxval( bzf_sf(1,:,:,:) )
           bzmin = minval( bzf_sf(1,:,:,:) )

           write(*,*)'x-direction'
           write(*,*)'bxmin,bxmax=',bxmin,bxmax
           write(*,*)'bymin,bymax=',bymin,bymax
           write(*,*)'bzmin,bzmax=',bzmin,bzmax

           write(*,*)'adipolx,adipoly,adipolz:1=',
     &                adipolx,adipoly,adipolz
      write(*,*)'bxf_sf(1,5,1,1)=',bxf_sf(1,5,1,1)
      write(*,*)'bxf_sf(2,5,1,1)=',bxf_sf(2,5,1,1)
      write(*,*)'bxf_sf(3,5,1,1)=',bxf_sf(3,5,1,1)
#endif /*DEBUGG*/               

!pmn
            if(neigh(1,1,lb).eq.-30) then
c              bxf_sf(1:3,1:nguard,:,:) = 0.
            endif
!pmn

!     Loop over the y-faces of the cells.
      do iz=1,nz
         do iy=1,nyp
            do ix=1,nx
               
!  compute the coordinates of the lower-left corner                     

c               x1 = xmin + delx*(ix-0.5)
c               y1 = ymin + dely*(iy-1  )*k2d
c               z1 = zmin + delz*(iz-0.5)*k3d
               x1 = xmin + delx*(ix-nguard-0.5)
               y1 = ymin + dely*(iy-nguard-1  )*k2d
#if N_DIM == 3
               z1 = zmin + delz*(iz-nguard-0.5)*k3d
#else
               z1 = 0.
#endif
               
               r2  = x1*x1 + y1*y1 + z1*z1
!  There is roundoff error in determining the coordinates
!  of the grid points from the coord() and bsize() arrays.
!  We use the SVD like trick of setting the inverse to zero.
               if( r2 .gt. tiny) then
                  r1i = 1./sqrt(r2)
               else
                  r1i = 0.
               endif
               r2i = r1i*r1i
               r3i = r2i*r1i
#if N_DIM < 3
               rdi = r2i*r2i
               z1 = 0.               !pmn
#else /*NDIM3*/
               rdi = r3i*r2i
!pmn               z1 = 0.
#endif /*NDIM3*/
               
               temp = adipolx*x1+adipoly*y1+adipolz*z1
               bx1 = (3.*x1*temp - r2*adipolx)*rdi
               by1 = (3.*y1*temp - r2*adipoly)*rdi
               bz1 = (3.*z1*temp - r2*adipolz)*rdi

               bxa = bx1
               bya = by1
               bza = bz1
               
#ifdef IMAGE_DIPOLE
c        Now add symmetry dipole of same sign, but located at 
c        (2.0*xsym,0,0) instead of (0,0,0).

               x1 = x1 - 2.*xsym
               y1 = y1 - 2.*ysym
               z1 = z1 - 2.*zsym
               
               r2  = x1*x1 + y1*y1 + z1*z1
               if( r2 .gt. tiny ) then
                  r1i = 1./sqrt(r2)
               else
                  r1i = 0.
               endif
               r2i = r1i*r1i
               r3i = r2i*r1i
#if N_DIM < 3
               rdi = r2i*r2i
               z1 = 0.                !pmn
#else /*NDIM3*/
               rdi = r3i*r2i
!pmn               z1 = 0.
#endif /*NDIM3*/
               
               temp = adipolx*x1+adipoly*y1+adipolz*z1
               bx1 = (3.*x1*temp - r2*adipolx)*rdi
               by1 = (3.*y1*temp - r2*adipoly)*rdi
               bz1 = (3.*z1*temp - r2*adipolz)*rdi
               
! reverse direction of image dipole
!               bx1 = -bx1
!               by1 = -by1
!               bz1 = -bz1
               
               bxa = bxa+bx1
               bya = bya+by1
               bza = bza+bz1

#endif /*IMAGE_DIPOLE*/
#if N_DIM < 3
               bza = 0.0
#endif
c               byf_sf(1,ix,iy,iz) = bya
               byf_sf(2,ix,iy,iz) = bza
               byf_sf(3,ix,iy,iz) = bxa
        enddo               ! do ix=1,nxp
       enddo                  ! do iy=1,nyp
      enddo                     ! do iz=1,nzp
!pmn
#ifdef DEBUGG
           write(*,*)'y-direction'
           bymax = maxval( byf_sf(1,:,:,:) )
           bymin = minval( byf_sf(1,:,:,:) )
           bzmax = maxval( byf_sf(2,:,:,:) )
           bzmin = minval( byf_sf(2,:,:,:) )
           bxmax = maxval( byf_sf(3,:,:,:) )
           bxmin = minval( byf_sf(3,:,:,:) )
           write(*,*)'bxmin,bxmax=',bxmin,bxmax
           write(*,*)'bymin,bymax=',bymin,bymax
           write(*,*)'bzmin,bzmax=',bzmin,bzmax
      write(*,*)'byf_sf(1,5,1,1)=',byf_sf(1,5,1,1)
      write(*,*)'byf_sf(2,5,1,1)=',byf_sf(2,5,1,1)
      write(*,*)'byf_sf(3,5,1,1)=',byf_sf(3,5,1,1)
#endif /*DEBUGG*/               

            if(neigh(1,1,lb).eq.-30) then
c              byf_sf(1:3,1:nguard,:,:) = 0.
            endif
!pmn

c     Loop over the z-faces of the cells.
      do iz=1,nzp
         do iy=1,ny
            do ix=1,nx
               
c  compute the coordinates of the lower-left corner                     

c               x1 = xmin + delx*(ix-0.5)
c               y1 = ymin + dely*(iy-0.5)*k2d
c               z1 = zmin + delz*(iz-1  )*k3d
               x1 = xmin + delx*(ix-nguard-0.5)
               y1 = ymin + dely*(iy-nguard-0.5)*k2d
#if N_DIM == 3
               z1 = zmin + delz*(iz-nguard-1  )*k3d
#else
               z1 = 0.
#endif

               r2 = x1*x1 + y1*y1 + z1*z1
c  There is roundoff error in determining the coordinates
c  of the grid points from the coord() and bsize() arrays.
c  We use the SVD like trick of setting the inverse to zero.
               if( r2 .gt. tiny) then
                  r1i = 1./sqrt(r2)
               else
                  r1i = 0.
               endif
               r2i = r1i*r1i
               r3i = r2i*r1i
#if N_DIM < 3
               rdi = r2i*r2i
               z1 = 0.            !pmn
#else /*NDIM3*/
               rdi = r3i*r2i
!pmn               z1 = 0.
#endif /*NDIM3*/
               
               temp = adipolx*x1+adipoly*y1+adipolz*z1
               bx1 = (3.*x1*temp - r2*adipolx)*rdi
               by1 = (3.*y1*temp - r2*adipoly)*rdi
               bz1 = (3.*z1*temp - r2*adipolz)*rdi
               
               bxa = bx1
               bya = by1
               bza = bz1
               
#ifdef IMAGE_DIPOLE
c        Now add symmetry dipole of same sign, but located at 
c        (2.0*xsym,0,0) instead of (0,0,0).

               x1 = x1 - 2.*xsym
               y1 = y1 - 2.*ysym
               z1 = z1 - 2.*zsym
               
               r2  = x1*x1 + y1*y1 + z1*z1
               if( r2 .gt. tiny ) then
                  r1i = 1./sqrt(r2)
               else
                  r1i = 0.
               endif
               r2i = r1i*r1i
               r3i = r2i*r1i
#if N_DIM < 3
               rdi = r2i*r2i
               z1 = 0.                !pmn
#else /*NDIM3*/
               rdi = r3i*r2i
!pmn               z1 = 0.
#endif /*NDIM3*/
               
               temp = adipolx*x1+adipoly*y1+adipolz*z1
               bx1 = (3.*x1*temp - r2*adipolx)*rdi
               by1 = (3.*y1*temp - r2*adipoly)*rdi
               bz1 = (3.*z1*temp - r2*adipolz)*rdi
               
! reverse direction of image dipole
!               bx1 = -bx1
!               by1 = -by1
!               bz1 = -bz1

               bxa = bxa+bx1
               bya = bya+by1
               bza = bza+bz1

#endif /*IMAGE_DIPOLE*/
c               bzf_sf(1,ix,iy,iz) = bza
               bzf_sf(2,ix,iy,iz) = bxa
               bzf_sf(3,ix,iy,iz) = bya
        enddo               ! do ix=1,nxp
       enddo                  ! do iy=1,nyp
      enddo                     ! do iz=1,nzp1
c           write(*,*)'z-direction'
           bzmax = maxval( bzf_sf(1,:,:,:) )
           bzmin = minval( bzf_sf(1,:,:,:) )
           bxmax = maxval( bzf_sf(2,:,:,:) )
           bxmin = minval( bzf_sf(2,:,:,:) )
           bymax = maxval( bzf_sf(3,:,:,:) )
           bymin = minval( bzf_sf(3,:,:,:) )
#ifdef DEBUGG
           write(*,*)'bxmin,bxmax=',bxmin,bxmax
           write(*,*)'bymin,bymax=',bymin,bymax
           write(*,*)'bzmin,bzmax=',bzmin,bzmax
      write(*,*)'bzf_sf(1,5,1,1)=',bzf_sf(1,5,1,1)
      write(*,*)'bzf_sf(2,5,1,1)=',bzf_sf(2,5,1,1)
      write(*,*)'bzf_sf(3,5,1,1)=',bzf_sf(3,5,1,1)
#endif /*DEBUGG*/               
!pmn
            if(neigh(1,1,lb).eq.-30) then
c              bzf_sf(1:3,1:nguard,:,:) = 0.
            endif
!pmn
      

! Overwrite magnetic fields at cell centers from those on faces
      do iz=kl_bnd1,ku_bnd1
       do iy=jl_bnd1,ju_bnd1
        do ix=il_bnd1,iu_bnd1
         b_sf(1,ix,iy,iz)= 0.5*( bxf_sf(1,ix+1,iy    ,iz    )
     &                         + bxf_sf(1,ix  ,iy    ,iz    ))
         b_sf(2,ix,iy,iz)= 0.5*( byf_sf(1,ix  ,iy+k2d,iz    )+ 
     &                           byf_sf(1,ix  ,iy    ,iz    ))
         b_sf(3,ix,iy,iz)= 0.5*( bzf_sf(1,ix  ,iy    ,iz+k3d)+ 
     &                           bzf_sf(1,ix  ,iy    ,iz    ))
!-----------------------------------------------------------------------
        end do
       end do
      end do
      if(neigh(1,1,lb).eq.-30) then
c       b_sf(1:3,1:nguard,:,:) = 0.
      endif
#ifdef DEBUGGG1
      write(*,*)'b_sf(1,4,45,1)=',b_sf(1,4,45,1)
      write(*,*)'b_sf(2,4,45,1)=',b_sf(2,4,45,1)
      write(*,*)'b_sf(3,4,45,1)=',b_sf(3,4,45,1)
      write(*,*)'bxf_sf(1,5,68,1)=',bxf_sf(1,5,68,1)
      write(*,*)'bxf_sf(2,5,68,1)=',bxf_sf(2,5,68,1)
      write(*,*)'bxf_sf(1,5,95,1)=',bxf_sf(1,5,95,1)
      write(*,*)'bxf_sf(2,5,95,1)=',bxf_sf(2,5,95,1)
c      pause
#endif

#ifdef DEBUGGFF
      write(*,*)'b_sf(1,4,45,1)=',b_sf(1,4,45,1)
      write(*,*)'b_sf(2,4,45,1)=',b_sf(2,4,45,1)
      write(*,*)'b_sf(3,4,45,1)=',b_sf(3,4,45,1)
      write(*,*)'bxf_sf(1,4,45,1)=',bxf_sf(1,4,45,1)
      write(*,*)'bxf_sf(1,5,45,1)=',bxf_sf(1,5,45,1)
      pause
#endif /*DEBUGG*/               



#endif /*MAGNETOSPHERE*/

#ifdef EXPLOSION
      do iz=kl_bnd1,ku_bnd1
       do iy=jl_bnd1,ju_bnd1
        do ix=il_bnd1,iu_bnd1
         b_sf(1,ix,iy,iz) = 0.0
         b_sf(2,ix,iy,iz) = 0.0 !<==1.1
         b_sf(3,ix,iy,iz) = 0.0
        enddo               ! do ix=1,nxp
       enddo                  ! do iy=1,nyp
      enddo                     ! do iz=1,nzp
      do k=kl_bnd1,ku_bnd1
       do j=jl_bnd1,ju_bnd1
        do i=il_bnd1+1,iu_bnd1
         bxf_sf(1,i,j,k)=0.5*(b_sf(1,i,j,k)+b_sf(1,i-1,j,k))
         bxf_sf(2,i,j,k)=0.5*(b_sf(2,i,j,k)+b_sf(2,i-1,j,k))
         bxf_sf(3,i,j,k)=0.5*(b_sf(3,i,j,k)+b_sf(3,i-1,j,k))
        enddo
       enddo
      enddo
      do k=kl_bnd1,ku_bnd1
       do j=jl_bnd1+k2d,ju_bnd1
        do i=il_bnd1,iu_bnd1
         byf_sf(1,i,j,k)=0.5*(b_sf(2,i,j,k)+b_sf(2,i,j-1,k))
         byf_sf(2,i,j,k)=0.5*(b_sf(3,i,j,k)+b_sf(3,i,j-1,k))
         byf_sf(3,i,j,k)=0.5*(b_sf(1,i,j,k)+b_sf(1,i,j-1,k))
        enddo
       enddo
      enddo
#if N_DIM == 3
      do k=kl_bnd1+k3d,ku_bnd1
       do j=jl_bnd1,ju_bnd1
        do i=il_bnd1,iu_bnd1
         bzf_sf(1,i,j,k)=0.5*(b_sf(3,i,j,k)+b_sf(3,i,j,k-1))
         bzf_sf(2,i,j,k)=0.5*(b_sf(1,i,j,k)+b_sf(1,i,j,k-1))
         bzf_sf(3,i,j,k)=0.5*(b_sf(2,i,j,k)+b_sf(2,i,j,k-1))
        enddo
       enddo
      enddo
#else
      bzf_sf(:,:,:,:)=0.
#endif

c      bxf_sf(:,:,:,:)=0.
c      byf_sf(:,:,:,:)=0.
c      bzf_sf(:,:,:,:)=0.

#endif /*EXPLOSION*/
#ifdef DIVERGENCE_TEST
      do iz=kl_bnd1,ku_bnd1
       do iy=jl_bnd1,ju_bnd1
        do ix=il_bnd1,iu_bnd1
         if(coord(1,lb)<0.0)then
           b_sf(1,ix,iy,iz) = 0.0
           b_sf(2,ix,iy,iz) = 0.0
           b_sf(3,ix,iy,iz) = 0.0
         else
           b_sf(1,ix,iy,iz) = 0.0
           b_sf(2,ix,iy,iz) = 0.0
           b_sf(3,ix,iy,iz) = 0.0
         endif
        enddo               
       enddo                 
      enddo                     
      do k=kl_bnd1,ku_bnd1
       do j=jl_bnd1,ju_bnd1
        do i=il_bnd1+1,iu_bnd1
         bxf_sf(1,i,j,k)=0.5*(b_sf(1,i,j,k)+b_sf(1,i-1,j,k))
         bxf_sf(2,i,j,k)=0.5*(b_sf(2,i,j,k)+b_sf(2,i-1,j,k))
         bxf_sf(3,i,j,k)=0.5*(b_sf(3,i,j,k)+b_sf(3,i-1,j,k))
        enddo
       enddo
      enddo
      do k=kl_bnd1,ku_bnd1
       do j=jl_bnd1+k2d,ju_bnd1
        do i=il_bnd1,iu_bnd1
         byf_sf(1,i,j,k)=0.5*(b_sf(2,i,j,k)+b_sf(2,i,j-1,k))
         byf_sf(2,i,j,k)=0.5*(b_sf(3,i,j,k)+b_sf(3,i,j-1,k))
         byf_sf(3,i,j,k)=0.5*(b_sf(1,i,j,k)+b_sf(1,i,j-1,k))
        enddo
       enddo
      enddo
#if N_DIM == 3
      do k=kl_bnd1+k3d,ku_bnd1
       do j=jl_bnd1,ju_bnd1
        do i=il_bnd1,iu_bnd1
         bzf_sf(1,i,j,k)=0.5*(b_sf(3,i,j,k)+b_sf(3,i,j,k-1))
         bzf_sf(2,i,j,k)=0.5*(b_sf(1,i,j,k)+b_sf(1,i,j,k-1))
         bzf_sf(3,i,j,k)=0.5*(b_sf(2,i,j,k)+b_sf(2,i,j,k-1))
        enddo
       enddo
      enddo
#else
      bzf_sf(:,:,:,:)=0.
#endif

c      bxf_sf(:,:,:,:)=0.
c      byf_sf(:,:,:,:)=0.
c      bzf_sf(:,:,:,:)=0.


#endif /*DIVERGENCE_TEST*/
#ifdef DEBUGG
           bxmax = maxval( b_sf(1,:,:,:) )
           bxmin = minval( b_sf(1,:,:,:) )
           bymax = maxval( b_sf(2,:,:,:) )
           bymin = minval( b_sf(2,:,:,:) )
           bzmax = maxval( b_sf(3,:,:,:) )
           bzmin = minval( b_sf(3,:,:,:) )
           write(*,*)'b_sf: within strongfields_wb: 1' 
           write(*,*)'bxmin,bxmax = ',bxmin,bxmax
           write(*,*)'bymin,bymax = ',bymin,bymax
           write(*,*)'bzmin,bzmax = ',bzmin,bzmax
#endif /*DEBUG*/
!     b_sf(:,:,:,:)=0.
!     bxf_sf(:,:,:,:)=0.
!     byf_sf(:,:,:,:)=0.
!     bzf_sf(:,:,:,:)=0.

!=======================================================================

      return
      end subroutine strongfields_wb
!=======================================================================
!
!=======================================================================

      subroutine a0_magnetosphere_strongfield(ax ,ay ,az,pe,lb )
      use physcons
      use paramesh_dimensions
      use physicaldata

! include file defining the tree
      use tree
c      use strongfields
      implicit none


#include "amr_shmem.fh"
#include "shmem_reduce.fh"
      real, dimension(il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,
     &                kl_bnd1:ku_bnd1+k3d) ::ax,ay,az

      real, dimension(         iu_bnd1 +1) :: xpp,ypp,zpp
      integer :: ix,iy,iz,lb
      integer :: nx, nxp, ny, nyp, nz, nzp


      real :: xpos,ypos,zpos

! Local Variables
      real ax1,ay1,az1
      real psym
      real adipole
      real adipol,alpha1
      real xp,yp,zp,phi,theta,omega,phi2
      real x1,y1,z1,xr,yr,zr,r2,r1i,r2i,r3i,rdi


! Mesh Variables
      real    :: xmin,xmax, delx, rdx
      real    :: ymin,ymax, dely, rdy
      real    :: zmin,zmax, delz, rdz

      real, parameter  :: tiny=(1e-1)**2
      real, save :: ccbsize(3),bbox(2,3),ccoord(3) 
      integer :: pe
      integer ilb, pe0, lb0

#ifdef MPI_USED
      include "mpif.h"
      integer,external :: shmem_my_pe
      integer          :: mype,ierrorcode,ierr
#endif /* MPI_USED */

#ifdef MPI_USED
      mype = shmem_my_pe()
!     if(pe.ne.mype) then
!       write(*,*) 'Error: a0 : calling shmem',
!    .             ' for a non-local block'
!       call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
!     endif
#endif /* MPI_USED */


c      call shmem_real_get(ccbsize(1),bsize(1,lb),3,pe)

      pe0=pe; lb0=lb;
#ifdef MPI_USED
      if(pe.ne.mype) then
!write(*,*)'--------entering finding ilb'
         lb0=0;  ilb=strt_buffer;

         do while( (ilb.ge.strt_buffer.and.
     &              ilb.le.last_buffer).and.
     &              lb0.eq.0)
            if(lb.eq.laddress(1,ilb).and.pe.eq.laddress(2,ilb)) then
                lb0=ilb
                pe0=mype
            endif
            ilb=ilb+1
         end do

         if(lb0.eq.0) then
            write(*,*) 'Error: strongfields_wb : could not find block',
     .             'and pe in laddress ------ stopping'
            call amr_close
            stop
         endif
!write(*,*)'--------exiting finding ilb'

      endif
#endif

      call shmem_real_get(bbox(1,1),bnd_box(1,1,lb0),6,mype)
      call shmem_real_get(ccoord(1),coord(1,lb0),3,mype)

      ccbsize(:) = bbox(2,:)-bbox(1,:)
      nx  = nxb+2*nguard
      nxp = nx+1
      ny  = nyb+2*nguard*k2d
      nyp = ny+k2d
      nz  = nzb+2*nguard*k3d
      nzp = nz+k3d


      ax(:,:,:)=0.
      ay(:,:,:)=0.
      az(:,:,:)=0.

c      xmax  = ccoord(1) + 0.5*ccbsize(1)
c      xmin  = ccoord(1) - 0.5*ccbsize(1)
      xmax  = bbox(1,1)
      xmin  = bbox(2,1)
      delx  = (xmax - xmin)/float(nxb)
      rdx   = 1./delx

c      ymax  = ccoord(2) + 0.5*ccbsize(2)
c      ymin  = ccoord(2) - 0.5*ccbsize(2)
      ymax  = bbox(1,2)
      ymin  = bbox(2,2)
      dely  = (ymax - ymin)/float(nyb)
      rdy   = 1./dely
#if N_DIM == 3
c      zmax  = ccoord(3) + 0.5*ccbsize(3)
c      zmin  = ccoord(3) - 0.5*ccbsize(3)
      zmax  = bbox(1,3)
      zmin  = bbox(2,3)
      delz  = (zmax - zmin)/float(nzb)
      rdz   = 1./delz
#else
      zmax  = 0.
      zmin  = 0.
      delz  = 0.
      rdz   = 0.                ! The generalized inverse.
#endif /*NDIM3*/
!
!--------coordinates of block (not working block)
!
      xmin = bbox(1,1) 
      xmax = bbox(2,1) 
      ymin = bbox(1,2) 
      ymax = bbox(2,2)
      delx = (xmax-xmin)/real(nxb) 
      dely = (ymax-ymin)/real(nyb) 
      rdx  = 1./delx
      rdy  = 1./dely
c      xmin = xmin-nguard*delx
c      xmax = xmax+nguard*delx 
c      ymin = ymin-nguard*dely 
c      ymax = ymax+nguard*dely
#if N_DIM == 3
      zmin = bbox(1,3) 
      zmax = bbox(2,3) 
      delz = (zmax-zmin)/real(nzb)
      rdz  = 1./delz
c      zmin = zmin-nguard*delz
c      zmax = zmax+nguard*delz 
#else
      zmax  = 0.
      zmin  = 0.
      delz  = 0.
      rdz   = 0.                ! The generalized inverse.
#endif
      
!      write(*,*)'delx,dely,delz=',delx,dely,delz
!      write(*,*)'rdx,rdy,rdz   =',rdx,rdy,rdz   
!      write(*,*)'adipolx,adipoly,adipolz   =',adipolx,adipoly,adipolz
!      write(*,*)'xsym,ysym,zsym=',xsym,ysym,zsym

! COMPUTE THE POINT VECTOR POTENTIALS
! Loop over the corner points of the cells.
c      write(*,*)'xmin,ymin,zmin=',xmin,ymin,zmin
      do iz=1,nzp
         do iy=1,nyp
            do ix=1,nxp
               
! compute the coordinates of the lower-left corner                     
c               xpos    = xmin+delx*(ix-1)
c               ypos    = ymin+dely*(iy-1)
c               zpos    = zmin+delz*(iz-1)
               xpos    = xmin + delx*(ix-nguard-1)
               ypos    = ymin + dely*(iy-nguard-1)
#if N_DIM == 3
               zpos    = zmin + delz*(iz-nguard-1)
#else
               zpos    = 0.
#endif 
               xpp(ix) = xpos
               ypp(ix) = ypos
               zpp(ix) = zpos*k3d
               
         enddo
            
        do ix=1,nxp
               
               x1 = xpp(ix)
               y1 = ypp(ix)
               z1 = zpp(ix)
               
               r2  = x1*x1 + y1*y1 + z1*z1
c  There is roundoff error in determining the coordinates
c  of the grid points from the coord() and bsize() arrays.
c  We use the SVD like trick of setting the inverse to zero.
               if( r2 .gt. tiny) then
                  r1i = 1./sqrt(r2)
               else
                  r1i = 0.
               endif
               r2i = r1i*r1i
               r3i = r2i*r1i
#if N_DIM < 3
               rdi = r2i
               z1  = 0.
#else /*NDIM3*/
               rdi = r3i
#endif /*NDIM3*/

               ax1 = (z1*adipolyt - y1*adipolzt)*rdi
               ay1 = (x1*adipolzt - z1*adipolxt)*rdi
               az1 = (y1*adipolxt - x1*adipolyt)*rdi
               
               ax(ix,iy,iz) = ax1
               ay(ix,iy,iz) = ay1
               az(ix,iy,iz) = az1
#ifdef DEBUGG
         if(ix.eq.5.and.iy.eq.1) then
            write(*,*) 'before image'
            write(*,*) 'adipolx,adipoly=',adipolx,adipoly
            write(*,*) 'x1,y1,z1       =',x1,y1,z1
            write(*,*) 'ax1,ay1,az1    =',ax1,ay1,az1
            write(*,*) 'xmin,delx      =',xmin,delx
          endif
#endif /*DEBUGG*/               
              
#ifdef IMAGE_DIPOLE
! Now add symmetry dipole of same sign, but located at 
! (2.0*xsym,2.*ysym,2.*zsym) instead of (0,0,0).

               x1 = x1 - 2.*xsym
               y1 = y1 - 2.*ysym
               z1 = z1 - 2.*zsym

               
               r2  = x1*x1 + y1*y1 + z1*z1
               if( r2 .gt. tiny ) then
                  r1i = 1./sqrt(r2)
               else
                  r1i = 0.
               endif
               r2i = r1i*r1i
               r3i = r2i*r1i
#if N_DIM < 3
               rdi = r2i
               z1  = 0.
#else /*NDIM3*/
               rdi = r3i
#endif /*NDIM3*/
               ax1 = (z1*adipoly - y1*adipolz)*rdi
               ay1 = (x1*adipolz - z1*adipolx)*rdi
               az1 = (y1*adipolx - x1*adipoly)*rdi

! reverse sign of image dipole
!               ax1 = -ax1
!               ay1 = -ay1
!               az1 = -az1

               ax(ix,iy,iz) = ax(ix,iy,iz) + ax1
               ay(ix,iy,iz) = ay(ix,iy,iz) + ay1
               az(ix,iy,iz) = az(ix,iy,iz) + az1
               
#ifdef DEBUGG
         if(ix.eq.5.and.iy.eq.1) then
            write(*,*) 'after image'
            write(*,*) 'adipolx,adipoly=',adipolx,adipoly
            write(*,*) 'x1,y1,z1       =',x1,y1,z1
            write(*,*) 'xsym,ysym      =',xsym,ysym
            write(*,*) 'ax(ix,iy,iz)   =',ax(ix,iy,iz) 
            write(*,*) 'ay(ix,iy,iz)   =',ay(ix,iy,iz) 
            write(*,*) 'az(ix,iy,iz)   =',az(ix,iy,iz) 
          endif
#endif /*DEBUGG*/               

#endif /*IMAGE_DIPOLE*/
            enddo               ! do ix=1,nxp1
         enddo                  ! do iy=1,nyp1
      enddo                     ! do iz=1,nzp1

      return
      end subroutine a0_magnetosphere_strongfield
!=======================================================================
!
c------------------------------------------------------------------
!
c------------------------------------------------------------------
      subroutine rotate_dipole(adipolxx,adipolyy,adipolzz,time)
      use physcons
c      use strongfields
      implicit none

      real :: adipolxx,adipolyy,adipolzz
      real :: time

c------------------------------------------------------------------
c
c This routine computes the electric field in the rest frame used
c by ODIN, which is introduced by the time dependence of the Earth's
c dipole field. This time dependence can include daily rotation
c and annual orbit of the Earth about the sun.
c Written :  Peter MacNeice        October 1998

c
c arguments :
c     adipol            magnitude of earth's dipole moment
c     adipolxx           x-component of earth's dipole moment in ODIN's frame
c     adipolyy           y-component of earth's dipole moment in ODIN's frame
c     adipolzz           z-component of earth's dipole moment in ODIN's frame
c
c------------------------------------------------------------------



      real ::  alpha,omega,phi,phi0,theta,xi,xi0
      real ::  sintheta,costheta,sinxi,cosxi
      real ::  mx_2,my_2,mz_2
      real ::  coef

c------------------------------------------------------- 
c
c Set location of dipole
      xx0 = 0.
      yy0 = 0.
      zz0 = 0.
c#include "earth_rotation.fh"
c
c     coef  = 0.0, no tilt, coef =1.0 full tilt
c
      coef  = 0.0
c
c  Set earth's daily rotational angular velocity
      omega = 0.*2.*pie/86400.
c      omega = 1000.*2.*pie/86400.

c
c Earth's rotational axis is tilted by 23.5 degrees.
c alpha = tilt of earth's rotational axis
      alpha = coef*23.5*pie/180.

c
c theta = tilt of dipole relative to rotation axis of earth
      theta = coef*11.*pie/180.

c
c phi = angle of rotation associated with annual variation in the
c       plane containing the rotation axis of the earth
c       When phi=0 the earth is at winter solstice
c       ie azimuthal position of earth relative to winter solstice.

c Set initial time of year
c      phi0 = pie/2.                 ! vernal equinox
c      phi0 = pie                    ! summer solstice
c      phi0 = 3.*pie/2.              ! autumnal equinox
       phi0 = 0.                      ! winter solstice

c
c include earth's orbit about sun
      phi = phi0 + 0.*omega*time/365.0

c
c xi = angle associated with daily rotation of earth
c Set initial time of day
      xi0 = 0.

c
c xi measures the angle between the plane containing the earth's
c rotation axis and the magnetic moment vector at time t, and the
c equivalent plane at t=0. At t=0 this plane is assumed to be
c the night-side at winter solstice.
      xi = xi0 + omega*time

c Given a vector in the cartesian coordinate frame rotating with the earth,
c this routine computes the vector in the frame of reference of the ODIN
c calculation, which has x axis always pointing radially away from the 
c sun and z axis in the direction of the solar rotation vector.
c 

    
c------------------------------------------------------- 
c
c Earth's dipole moment in the rotation frame of earth 
      sintheta = sin(theta)
      costheta = cos(theta)
      sinxi    = sin(xi)
      cosxi    = cos(xi)
#if N_DIM == 2
      sinxi = 0.
      cosxi = 1.
#endif
      mx_2 =  adipol*sintheta*cosxi 
      my_2 =  adipol*sintheta*sinxi 
      mz_2 =  adipol*costheta

c      write(*,*)'within rotate_dipole:adipol,xi,theta=',adipol,xi,theta
c
c In ODIN frame (adds tilt of rotation axis and annual orbit)
      adipolxx=  cos(alpha)*cos(phi)*mx_2 + sin(phi)*my_2
     .         + sin(alpha)*cos(phi)*mz_2
      adipolyy= -cos(alpha)*sin(phi)*mx_2 + cos(phi)*my_2
     .         - sin(alpha)*sin(phi)*mz_2
      adipolzz= -sin(alpha)*mx_2 + cos(alpha)*mz_2


#if N_DIM == 2
! Put dipole into xy plane instead of xz plane
      adipolyy = adipolzz
      adipolzz = 0.
      adipolxx = 0.


! testing for symmetry
cdss      adipolx = 0.
cdss      adipoly = -adipol
! testing for symmetry


#endif

c      adipolxx = adipolx
c      adipolyy = adipoly
c      adipolzz = adipolz

      return
      end subroutine rotate_dipole
      end module strongfields
