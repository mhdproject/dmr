#include "athena_preprocessor.fh"
      subroutine div_V(vv,divv,idir)
      use paramesh_dimensions
      use physcons
      use grid
      use athena_specific_parameters

      use paramesh_interfaces
#ifdef MPI_USED
      use paramesh_mpi_interfaces
#endif /* MPI_USED */
      implicit none
#define cell_s real,dimension(   iu_bnd1 ,ju_bnd1 ,ku_bnd1 )
#define cell_v real,dimension(ne,iu_bnd1 ,ju_bnd1 ,ku_bnd1 )
      cell_v,intent(in   )::vv
      cell_s,intent(inout)::divv
      integer :: i,j,k,idir

      rdx = 1./delx
      rdy = 1./dely
      rdz = 1./delz

      if(idir==1)then
        do k= 1+k3d, ku_bnd1-k3d
         do j= 1+k2d, ju_bnd1-k2d
          do i= 1+k1d, iu_bnd1-k1d
           divv(i,j,k)= 
     &      ( vv(2,i  ,j  ,k  ) - vv(2,i-1,j  ,k  ))*rdx
     & +0.5*((vv(3,i  ,j+1,k  ) + vv(3,i-1,j+1,k  ))  
     &     - (vv(3,i  ,j-1,k  ) + vv(3,i-1,j-1,k  )))*rdy 
#if N_DIM == 3
     & +0.5*((vv(4,i  ,j  ,k+1) + vv(4,i-1,j  ,k+1))  
     &     - (vv(4,i  ,j  ,k-1) + vv(4,i-1,j  ,k-1)))*rdz 
#endif
             divv(i,j,k)=divv(i,j,k)*delx 
          enddo
         enddo
        enddo
      elseif(idir==2)then
        do k= 1+k3d, ku_bnd1-k3d
         do j= 1+k2d, ju_bnd1-k2d
          do i= 1+k1d, iu_bnd1-k1d
           divv(i,j,k)= 
     &      ( vv(3,i  ,j  ,k  ) - vv(3,i  ,j-1,k  ) )*rdy
     &+ 0.5*((vv(2,i+1,j  ,k  ) + vv(2,i+1,j-1,k  ))  
     &     - (vv(2,i-1,j  ,k  ) + vv(2,i-1,j-1,k  )))*rdx 
#if N_DIM == 3
     &+ 0.5*((vv(4,i  ,j  ,k+1) + vv(4,i  ,j-1,k+1))  
     &     - (vv(4,i  ,j  ,k-1) + vv(4,i  ,j-1,k-1)))*rdz 
#endif
           divv(i,j,k)=divv(i,j,k)*dely
          end do
         end do
        end do
      elseif(idir==3)then
c       divv=0.0
        do k= 1+k3d, ku_bnd1-k3d
         do j= 1+k2d, ju_bnd1-k2d
          do i= 1+k1d, iu_bnd1-k1d
          divv(i,j,k)= 
     &       (vv(4,i  ,j  ,k  ) - vv(4,i  ,j  ,k-1) )*rdz
     &+ 0.5*((vv(2,i+1,j  ,k  ) + vv(2,i+1,j  ,k-1))  
     &     - (vv(2,i-1,j  ,k  ) + vv(2,i-1,j  ,k-1)))*rdx 
     &+ 0.5*((vv(3,i  ,j+1,k  ) + vv(3,i  ,j+1,k-1))  
     &     - (vv(3,i  ,j-1,k  ) + vv(3,i  ,j-1,k-1)))*rdy 
          divv(i,j,k)=divv(i,j,k)*delz 
         end do
        end do
       end do
      endif 
      return
      end subroutine div_V
