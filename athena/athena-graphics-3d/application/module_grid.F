#include "paramesh_preprocessor.fh"
      module grid
      implicit none
      private
      save
        public grid_variables
        real,public:: delx,dely ,delz 
        real,public:: rdx ,rdy  ,rdz
        real,public:: rvx ,rvy  ,rvz
        real,public:: rv
        real,public:: dtdxi,dtdyi,dtdzi
        contains
        subroutine grid_variables(lb,dt)
        use paramesh_dimensions
        use tree
        implicit none
        real    :: dt
	integer :: lb
        delx = bsize(1,lb)/real(nxb-gc_off_x*k1d)
        dely = bsize(2,lb)/real(nyb-gc_off_y*k2d)
        delz = bsize(3,lb)/real(nzb-gc_off_z*k3d)
c        write(*,*)'lb=',lb
c        write(*,*)'delx,dely,delz=',delx,dely,delz
c        write(*,*)'gc_off_x,gc_off_y,nzb-gc_off_z=',
c     &             gc_off_x,gc_off_y,nzb-gc_off_z
c        write(*,*)'k1d,k2d,k3d=',k1d,k2d,k3d
#if N_DIM == 1
        dely = 1.0
        delz = 1.0
#endif
#if N_DIM == 2
        delz = 1.0
#endif
        rdx     = 1./delx
        rdy     = 1./dely
        rdz     = 1./delz
        rv      = rdx*rdy*rdz*dt
        rvx     = dely*delz
        rvy     = delx*delz
        rvz     = delx*dely
        rvx     = dt*rdx
        rvy     = dt*rdy
        rvz     = dt*rdz
        dtdxi   = dt*rdx
        dtdyi   = dt*rdy
        dtdzi   = dt*rdz
#if N_DIM == 2
        rvz     = 0.0
        dtdzi   = 0.0
#endif

        end subroutine grid_variables
        end module grid
