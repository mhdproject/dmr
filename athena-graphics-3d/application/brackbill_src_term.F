#include "paramesh_preprocessor.fh"
      subroutine brackbill(v2,src,dt)
      use paramesh_dimensions
      use physcons
      use physicaldata
      use tree
      use strongfields
      use grid
      use athena_specific_parameters

      use paramesh_interfaces
#ifdef MPI_USED
      use paramesh_mpi_interfaces
#endif /* MPI_USED */
      implicit none
      real                                 ,intent(in)    ::dt
      real,dimension(ne,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,
     .                  kl_bnd1:ku_bnd1    ),intent(in)   ::v2
      real,dimension(ne,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,
     .                  kl_bnd1:ku_bnd1    ),intent(inout)::src
      real    :: bx,by,bz,bx1,by1,bz1
      real    :: vx,vy,vz
      real    :: divb,ub
      integer :: ix,iy,iz
!-----------------------------------------------------------------------
!
      do iz=1+nguard*k3d,nzb+nguard*k3d
       do iy=1+nguard*k2d,nyb+nguard*k2d
        do ix=1+nguard  ,nxb+nguard
 
         vx   = v2(2,ix,iy,iz)
         vy   = v2(3,ix,iy,iz)
         vz   = v2(4,ix,iy,iz)
         bx1  = v2(6,ix,iy,iz)
         by1  = v2(7,ix,iy,iz)
         bz1  = v2(8,ix,iy,iz)
         bx   = bx1+b_sf(1,ix,iy,iz)
         by   = by1+b_sf(2,ix,iy,iz)
         bz   = bz1+b_sf(3,ix,iy,iz)
         divb = 0.
         divb = divb + 0.5*(v2(6,ix+1,iy    ,iz    ) -
     &                      v2(6,ix-1,iy    ,iz    ))*rdx
         divb = divb + 0.5*(v2(7,ix  ,iy+k2d,iz    ) -
     &                      v2(7,ix  ,iy-k2d,iz    ))*rdy
         divb = divb + 0.5*(v2(8,ix  ,iy    ,iz+k3d) -
     &                      v2(8,ix  ,iy    ,iz-k3d))*rdz
         divb = divb*dt

         ub   = vx*bx1 + vy*by1 + vz*bz1 ! <== is this correct ?

         src(1,ix,iy,iz) =  0.0
         src(2,ix,iy,iz) = -divb*bx*pie4i
         src(3,ix,iy,iz) = -divb*by*pie4i
         src(4,ix,iy,iz) = -divb*bz*pie4i
c        src(5,ix,iy,iz) = -divb*ub*pie4i
         src(5,ix,iy,iz) =  0.0
         src(6,ix,iy,iz) = -divb*vx
         src(7,ix,iy,iz) = -divb*vy
         src(8,ix,iy,iz) = -divb*vz

        enddo
       enddo
      enddo
      return
      end subroutine brackbill
