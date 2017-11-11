#include "paramesh_preprocessor.fh"
#include "athena_preprocessor.fh"
#include "jet_conditions.fh"

      subroutine amr_initial_grid_jet
     &      (mype,xminb,yminb,zminb,xmaxb,ymaxb,zmaxb)
      use physcons
      use paramesh_dimensions
      use physicaldata


      use paramesh_interfaces
#ifdef MPI_USED
      use paramesh_mpi_interfaces
#endif /* MPI_USED */



! include file defining the tree
      use tree
      use strongfields
      implicit none
      integer :: mype

      real :: xminb,yminb,zminb,xmaxb,ymaxb,zmaxb




!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! $RCSfile: amr_initial_jet.F,v $
! $Revision: 1.8 $
! $Date: 2003/11/24 17:06:29 $


!
! This file is a template describing how the solution can be
! initialized on the initial grid. Modify it for your own use.
!
!--------------------------------------------------------------
! include files for amr

#include "amr_shmem.fh"
#include "shmem_reduce.fh"
!--------------------------------------------------------------
      real xmin,xmax,ymin,ymax,zmin,zmax,factor
      real :: time
      common/evolution/time
      integer :: nxc, nyc, nzc
      common/ncentral/nxc,nyc,nzc
     
      integer, parameter:: ng0= nguard*npgs
      
      integer :: ifaces


! set up initial grid aspect ratios.
      xmin =   JET_XMIN
      xmax =   JET_XMAX
      ymin =   JET_YMIN
      ymax =   JET_YMAX
      zmin =  -100.
      zmax =   100.


cggc change factor to 1
#if N_DIM < 3
      factor = 2.0
      factor = 1.0
#else
      factor = 1.0
#endif

      xmin = xmin * factor
      xmax = xmax * factor
      ymin = ymin * factor
      ymax = ymax * factor
      zmin = zmin * factor
      zmax = zmax * factor
      xminb= xmin
      xmaxb= xmax
      yminb= ymin
      ymaxb= ymax
      zminb= zmin
      zmaxb= zmax


      nxc  = iu_bnd/2+1
      nyc  = ju_bnd/2+1
#if N_DIM == 3
      nzc  = ku_bnd/2+1
#else
      nzc  = 1
#endif      

! set up step function on a single block covering the whole cubic domain
      lnblocks = 0
      if(mype .eq. 0 ) then
         lnblocks       = 1
         nodetype(1)    = 1
         lrefine(1)     = 1
         refine(1)      = .true.
         coord(1,1)     = 0.5*(xmax+xmin)
         coord(2,1)     = 0.5*(ymax+ymin)
         bsize(1,1)     = (xmax-xmin)
         bsize(2,1)     = (ymax-ymin)
         bnd_box(1,1,1) = xmin
         bnd_box(2,1,1) = xmax
         bnd_box(1,2,1) = ymin
         bnd_box(2,2,1) = ymax
#if N_DIM == 3
         coord(3,1)     = 0.5*(zmax+zmin)
         bsize(3,1)     = (zmax-zmin)
         bnd_box(1,3,1) = zmin
         bnd_box(2,3,1) = zmax
#endif

! first make all boundaries periodic
         do ifaces=1,nfaces
            neigh(1,ifaces,1) = 1 ! the initial block is its own
            neigh(2,ifaces,1) = 0 ! neighbor,  that is,  periodic BCs
         enddo


! Now set x-boundaries to be free

                                ! Now over-ride the default in the x-direction
         neigh(1,1,1) = -38     ! initial block is not its own
         neigh(2,1,1) = -38     ! neighbor. (static boundary?)
         neigh(1,2,1) = -21     ! initial block is not its own
         neigh(2,2,1) = -21     ! neighbor. (free boundary)

! Now set y-boundaries to be free

         neigh(1,3,1) = -40     ! initial block is not its own
         neigh(2,3,1) = -40     ! neighbor. (static boundary?)

         neigh(1,4,1) = -50     ! initial block is not its own
         neigh(2,4,1) = -50     ! neighbor. (free boundary)
#if N_DIM == 3
         neigh(1,5,1) = -21     ! initial block is not its own
         neigh(2,5,1) = -21     ! neighbor. (free boundary)
         neigh(1,6,1) = -21     ! initial block is not its own
         neigh(2,6,1) = -21     ! neighbor. (free boundary)
#endif
      endif
c------NOTE: ALL of what follows must be seen by all processors
       boundary_index(1:6)     = -21

!       boundary_index(1  )     = -22
       boundary_index(1  )     = -38
       boundary_index(3  )     = -40
       boundary_index(4  )     = -50

! x boundaries
       boundary_box(1,2:3,1:2) = -1.e30
       boundary_box(2,2:3,1:2) =  1.e30
       boundary_box(1,1,1)     = -1.e30
       boundary_box(2,1,1)     =  xmin
       boundary_box(1,1,2)     =  xmax
       boundary_box(2,1,2)     =  1.e30
! y boundaries
      if(ndim.ge.2) then
        boundary_box(1,1,3:4) =  xmin
        boundary_box(2,1,3:4) =  xmax
        boundary_box(1,3,3:4) = -1.e30
        boundary_box(2,3,3:4) =  1.e30
        boundary_box(1,2,3)   = -1.e30
        boundary_box(2,2,3)   =  ymin
        boundary_box(1,2,4)   =  ymax
        boundary_box(2,2,4)   =  1.e30
      endif
! z boundaries
      if(ndim.eq.3) then
        boundary_box(1,1,5:6) =  xmin
        boundary_box(2,1,5:6) =  xmax
        boundary_box(1,2,5:6) =  ymin
        boundary_box(2,2,5:6) =  ymax
        boundary_box(1,3,5)   = -1.e30
        boundary_box(2,3,5)   =  zmin
        boundary_box(1,3,6)   =  zmax
        boundary_box(2,3,6)   =  1.e30
      endif


c      write(*,*)'initial jet: boundary_index ',boundary_index(1:6)
      return
      end

!-----------------------------------------------------
      subroutine amr_initial_jet(mype,xmin,ymin,zmin,xmax,ymax,zmax)
      use physcons
      use paramesh_dimensions
      use physicaldata
      use strongfields


      use paramesh_interfaces
#ifdef MPI_USED
      use paramesh_mpi_interfaces
#endif /* MPI_USED */



! include file defining the tree
      use tree
      implicit none

      integer :: mype
c
#ifdef MPI_USED
       include 'mpif.h'
#endif /* MPI_USED */
c
#ifdef MPI_USED
      integer :: my_id
      integer :: ierror
#endif /*MPI_USED*/

! $RCSfile: amr_initial_jet.F,v $
! $Revision: 1.8 $
! $Date: 2003/11/24 17:06:29 $


!
! This file is a template describing how the solution can be
! initialized on the initial grid. Modify it for your own use.
!
!--------------------------------------------------------------

#include "amr_shmem.fh"
#include "shmem_reduce.fh"
!--------------------------------------------------------------
      real xmin,xmax,ymin,ymax,zmin,zmax,factor
      real :: time
      common/evolution/time
      integer :: nxc,nyc,nzc

      common/ncentral/nxc,nyc,nzc
      integer, parameter:: ng0= nguard*npgs
      real, save :: ccbsize(3),bbox(2,3),ccoord(3) 
      real :: bx,by,bz,b2      



      real :: jet_den, jet_p, jet_vx,jet_vy,jet_vz, jet_radius, jet_beta,
     .     ambient_den, ambient_vx, ambient_vy, ambient_p
      common/jet_state/jet_den, jet_p, jet_vx,jet_vy,jet_vz,jet_radius, jet_beta,
     .     ambient_den, ambient_vx, ambient_vy, ambient_p

      real :: xpos, ypos, zpos, rpos, rpos2, distance,rposb, rpos2b, distanceb,  ek, eth, rr,
     .     shear, rdx, rdy, rdz, delx, dely, delz, 
     .     correction, p
      
      integer :: ib, i, j, k



      bx = 0.0
      by = 0.0
      bz = 0.0

! loop over leaf grid blocks
      if(lnblocks.gt.0) then
      do ib=1,lnblocks
#ifndef ADVANCE_ALL_LEVELS
      if(nodetype(ib).eq.1) then
#endif
         bbox(:,:)=bnd_box(:,:,ib)
         xmin  = bbox(1,1)
         xmax  = bbox(2,1)
         delx  = (xmax - xmin)/float(nxb)

         ymin  = bbox(1,2)
         ymax  = bbox(2,2)
         dely  = (ymax - ymin)/float(nyb)
#if N_DIM == 3
         zmin  = bbox(1,3)
         zmax  = bbox(2,3)
         delz  = (zmax - zmin)/float(nzb)
#else
         zmax  = 0.
         zmin  = 0.
         delz  = 0.
#endif

!      
!----initial jet parameters
!
!************************* Stuff for initialising **************************



C     Let jet_beta*R be the distance at which the shear first becomes non-1.  Then
C     the volume-averaged value of the shear is given by the equation
C     <shear> = (1-jet_beta^4)/(4*(1-jet_beta))
C     This is used to correct the initden and initp values so that the total
C     density and energy released are correct. */


c Turlough's suggested Jet Parameters
c Jet radius:  1
c Jet pressure:  1
c Jet density:  5/3
c Jet velocity:  30
c Ambient density:  5/3
c Ambient pressure: 1
c ambient velocity: 0

C end gmurphy added


c      ambient_p=jet_p
c      ambient_den=jet_den
c      ambient_vy=0.0


      ambient_p=AMBIENT_P
      ambient_den=AMBIENT_DEN
      ambient_vx=AMBIENT_VX
      ambient_vy=AMBIENT_VY


      jet_p=JET_P
      jet_den=JET_DEN
      jet_vx=JET_VX
      jet_vy=JET_VY
      jet_vz=JET_VZ




      correction=(1.0-jet_beta**4.0)/(4.0*(1.0-jet_beta))

        do k=kl_bnd+ng0*k3d,ku_bnd-ng0*k3d
          do j=jl_bnd+ng0*k2d,ju_bnd-ng0*k2d
            do i=il_bnd+ng0,iu_bnd-ng0

               unk(:,i,j,k,ib) = 0.0


               xpos = xmin + delx*(i-0.5)
               ypos = ymin + dely*(j-0.5)*k2d
#if N_DIM == 3
               zpos = zmin + delz*(k-0.5)
#else
               zpos = 0.
#endif
c               rpos2 = (xpos**2) + (ypos**2) + (zpos**2)
               rpos2 = ((xpos)**2) + ((zpos)**2)
               rpos = sqrt(rpos2)
               distance=rpos
               
               rpos2b = ((xpos+30)**2) + ((zpos+30)**2)
               rposb = sqrt(rpos2b)
               distanceb=rposb
! The vanilla state
               unk(1,i,j,k,ib) = ambient_den
               rr=ambient_den
               unk(2,i,j,k,ib) = ambient_den*ambient_vx
               unk(3,i,j,k,ib) = ambient_den*ambient_vy
               ek=0.5*rr*(ambient_vx*ambient_vx+ ambient_vy*ambient_vy)
               eth=ambient_p/(gamma-1.)
               unk(5,i,j,k,ib) =ek+eth ! no fields so far

#ifdef MHD
               unk(7, i , j , k ,ib)=10
#endif /* MHD */
!     This puts in the jet
              if(xpos .gt. (JET_XMAX /24 + ypos/1.73205) )then

                  rr=jet_den
                  unk(1,i,j,k,ib) = rr
                  unk(2,i,j,k,ib) = rr * jet_vx
                  unk(3,i,j,k,ib) = rr *jet_vy
                  unk(4,i,j,k,ib) = 0
                  p=jet_p
                  ek=0.5*rr*(jet_vx * jet_vx + jet_vy*jet_vy)
                  eth=p/(gamma-1.)
                  unk(5,i,j,k,ib) =ek+eth ! no fields so far

              endif


              
            enddo
         enddo
      enddo


!     END OF JET INTIALIZATION
!
!*******


      
      if(nfacevar.gt.0) then
! set values for facevarx
        do k=kl_bnd+ng0*k3d,ku_bnd-ng0*k3d
          do j=jl_bnd+ng0*k2d,ju_bnd-ng0*k2d
            do i=il_bnd+ng0,iu_bnd-ng0+1
!             facevarx(1,i,j,k,ib) = 0.0
              facevarx(1,i,j,k,ib) = bx
            enddo
          enddo
        enddo

! set values for facevary
        do k=kl_bnd+ng0*k3d,ku_bnd-ng0*k3d
          do j=jl_bnd+ng0*k2d,ju_bnd+(-ng0+1)*k2d
            do i=il_bnd+ng0,iu_bnd-ng0
!             facevary(1,i,j,k,ib) = 0.0
              facevary(1,i,j,k,ib) = by
            enddo
          enddo
        enddo

! set values for facevarz
        do k=kl_bnd+ng0*k3d,ku_bnd+(-ng0+1)*k3d
          do j=jl_bnd+ng0*k2d,ju_bnd-ng0*k2d
            do i=il_bnd+ng0,iu_bnd-ng0
!             facevarz(1,i,j,k,ib) = 0.0
              facevarz(1,i,j,k,ib) = bz
            enddo
          enddo
        enddo

        endif


#ifndef ADVANCE_ALL_LEVELS
      endif
#endif
      enddo ! end loop over grid blocks
      endif

      return
      end
