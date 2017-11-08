#include "athena_preprocessor.fh"
      subroutine plm_1d(vc,flat,bx,by,bz,vl,vr,n1,n2,dt,idir,lflag)
      use paramesh_dimensions
      use physcons
      use grid
      use strongfields
      use athena_specific_parameters
      implicit none
      integer :: i,j,k,kn,l,n1,n2,idir
      logical                         ,intent(in   )::lflag
      real,dimension(ne  ,     nmax  ),intent(in   )::vc
      real,dimension(ne  ,     nmax+1),intent(  out)::vl,vr
      real,dimension(          nmax  ),intent(in   )::flat
      real,dimension(nem1,     nmax  )              ::v,u,vd,dvd,dvm
      real,dimension(          nmax  )              ::bn
      real,dimension(          nmax  ),intent(in   )::bx,by,bz
      real,dimension(          nmax  )              ::bnd,cfast
      real,dimension(nem1,     nmax  )              ::eigv,dwdum
      real,dimension(   4,     nmax  )              ::speeds
      real,dimension(nem1,nem1,nmax  )              ::rhe ,lhe
      real,dimension(nem1,nem1,nmax  )              ::rhec,lhec
      real,dimension(nem1,     nmax  )              ::wp ,w, wm 
      real,dimension(nem1,     nmax  )              ::dw ,dwd,dvv
      real,dimension(          nmax  )              ::evmin,evmax
      real,dimension(          nmax  )              ::dvy,dvz
      real,dimension(nem1            )              ::chi
      real,dimension(nem1,     nmax  )              ::slv
      real :: dvl,dvr,dvc,dt,dtdx,signm,sig,sum1,sgnx,sgx,sc4
      real, parameter                               ::fact1=2.0/3.0
!
         chi(1)=2.
         chi(2)=1.
         chi(3)=2.
         chi(4)=1.
         chi(5)=2.
         chi(6)=1.
         chi(7)=2.

!
!-----first rotate conserved variables 
!
      select case (idir) 

      case (1)     ! x-direction
       dtdx = 0.5*dt*rdx
       do i=n1,n2
         v(1,i)   = vc(1,i)
         v(2,i)   = vc(2,i)
         v(3,i)   = vc(3,i)
         v(4,i)   = vc(4,i)
         v(5,i)   = vc(5,i)
         v(6,i)   = vc(7,i)
         v(7,i)   = vc(8,i)
         bn(i)    = vc(6,i)
       enddo
      case (2)     ! y-direction
       dtdx = 0.5*dt*rdy
       do i=n1,n2
         v(1,i)   = vc(1,i)
         v(2,i)   = vc(3,i)
         v(3,i)   = vc(4,i)
         v(4,i)   = vc(2,i)
         v(5,i)   = vc(5,i)
         v(6,i)   = vc(8,i)
         v(7,i)   = vc(6,i)
         bn(i)    = vc(7,i)
       enddo
      case (3)     ! z-direction
       dtdx = 0.5*dt*rdz
       do i=n1,n2
         v(1,i)   = vc(1,i)
         v(2,i)   = vc(4,i)
         v(3,i)   = vc(2,i)
         v(4,i)   = vc(3,i)
         v(5,i)   = vc(5,i)
         v(6,i)   = vc(6,i)
         v(7,i)   = vc(7,i)
         bn(i)    = vc(8,i)
       enddo
      end select
!
!-----compute eigenvectors and eigenvalues to construct reference states
!
       do i=n1,n2
         vd(1,i) = v(1,i)
         vd(2,i) = v(2,i)
         vd(3,i) = v(3,i)
         vd(4,i) = v(4,i)
         vd(5,i) = v(5,i)
         vd(6,i) = v(6,i)+by(i)
         vd(7,i) = v(7,i)+bz(i)
         bnd(i)  = bn(i) +bx(i)
       enddo
       dvy(:)=0.
       dvz(:)=0.
       do i=n1+1,n2-1
        dvy(i)=vd(3,i+1)-vd(3,i-1)
        dvz(i)=vd(4,i+1)-vd(4,i-1)
       enddo

      call eigenvalues(1,n2,vd,bnd,eigv,cfast)
      dvv=0.
      call mhd_eigenvectors ( rhe ,lhe ,rhec,lhec,eigv,vd   ,bnd   ,
     &                        dw  ,dvv ,1   ,1   ,n2  ,cfast,speeds,
     &                        dvy, dvz ,idir )
      do i=1,n2  
       evmax(i) = eigv(7,i)
       evmin(i) = eigv(1,i)
      end do

#ifdef SPATIAL_SECOND_ORDER
#ifdef MONOTONIZE_PRIMITIVES
      do i= 2,n2-1
       do j=1,ne-1
!     Compute left differences
        dvl = v(j,i  ) - v(j,i-1)
!     Compute right differences
        dvr = v(j,i+1) - v(j,i  )
!     Compute center differences
        dvc = (dvr+dvl)
!     Compute monotonized slope
        sgx       = sign(1.0,dvc    )
        sgnx     = sign(0.5,dvr*dvl)
        slv(j,i) = min(0.5*abs(dvc),2.0*abs(dvr),2.0*abs(dvl))*sgx
        slv(j,i) = (0.5+sgnx)*slv(j,i)

       enddo
      enddo
      do i=3,n2-2 !<======
       do j=1,ne-1
!     Compute left differences
        dvl   = v(j,i  ) - v(j,i-1)
!     Compute right differences
        dvr   = v(j,i+1) - v(j,i  )
        dvc   = (dvr+dvl)
!     Compute sign of central difference
        sc4   = fact1*((v(j,i+1)-0.25*slv(j,i+1))-
     &                 (v(j,i-1)+0.25*slv(j,i-1)))
        sgx   = sign(1.0,sc4    )
        sgnx  = sign(0.5,dvr*dvl)
        dvm(j,i)=min(abs(sc4),2.0*abs(dvr),2.0*abs(dvl))*sgx
        dvm(j,i) = (0.5+sgnx)*dvm(j,i)
        enddo
      enddo

      do i=3,n2-2
       do j=1,nem1
        dw(j,i) = dot_product(lhe(j,:,i),dvm(:,i))
       enddo
      end do

#else
!
!-----monotonize characteristics
!
      do i=1,n2
       do j=1,nem1
        w(j,i) = dot_product(lhe(j,:,i),v(:,i))
       enddo
      end do
!
      do i=2,n2
       do j=1,nem1
        wm(j,i) = dot_product(lhe(j,:,i),v(:,i-1))
       enddo
      end do
!
      do i=1,n2-1
       do j=1,nem1
        wp(j,i) = dot_product(lhe(j,:,i),v(:,i+1))
       enddo
      end do
!
#ifdef LIMITOR1
      do i=2,n2-1
       do j=1,nem1
         sgx(j) = max(0.0,sign(1.0,(wp(j,i)-w(j,i))*(w(j,i)-wm(j,i))))
!
         dw(j,i) = min(0.5*abs( wp(j,i)-wm(j,i)),
     &              chi(j)*abs( wp(j,i)-w (j,i)),
     &              chi(j)*abs( w (j,i)-wm(j,i)))*
     &                 sign(1.0,wp(j,i)-w (j,i))
         dw(j,i) = dw(j,i)*sgx(j)
       end do
      end do
#endif /*LIMITOR1*/
!
!     the changes in the primitives are:
!
      do i=2,n2-1
       dvm(:,i) = matmul(rhe(:,:,i),dw(:,i))
      end do

#endif
#else
      do i=2,n2-1
       do j=1,nem1
!     Compute left differences
          dvl      = v(j,i  ) - v(j,i-1)
!     Compute right differences
          dvr      = v(j,i+1) - v(j,i  )
!     Compute center differences
          dvm(j,i) = 0.5*(dvr+dvl)
       enddo
      end do

      dvm(:,:) = 0.0

#endif /*SPATIAL_SECOND_ORDER*/

!
!-----apply flattening here
!
      do i = 4,n2-2
       do j=1,nem1
        dvm(j,i) = flat(i)*dvm(j,i)
       enddo
      enddo
!
!-----left and right reference states
!
      do i = 4,n2-2
         vl(1,i)=v(1,i-1)+(0.5-max(evmax(i-1),0.0)*dtdx)*dvm(1,i-1)
         vl(2,i)=v(2,i-1)+(0.5-max(evmax(i-1),0.0)*dtdx)*dvm(2,i-1)
         vl(3,i)=v(3,i-1)+(0.5-max(evmax(i-1),0.0)*dtdx)*dvm(3,i-1)
         vl(4,i)=v(4,i-1)+(0.5-max(evmax(i-1),0.0)*dtdx)*dvm(4,i-1)
         vl(5,i)=v(5,i-1)+(0.5-max(evmax(i-1),0.0)*dtdx)*dvm(5,i-1)
         vl(6,i)=bn(i-1)
         vl(7,i)=v(6,i-1)+(0.5-max(evmax(i-1),0.0)*dtdx)*dvm(6,i-1)
         vl(8,i)=v(7,i-1)+(0.5-max(evmax(i-1),0.0)*dtdx)*dvm(7,i-1)

         vr(1,i)=v(1,i  )-(0.5+min(evmin(i  ),0.0)*dtdx)*dvm(1,i  )
         vr(2,i)=v(2,i  )-(0.5+min(evmin(i  ),0.0)*dtdx)*dvm(2,i  )
         vr(3,i)=v(3,i  )-(0.5+min(evmin(i  ),0.0)*dtdx)*dvm(3,i  )
         vr(4,i)=v(4,i  )-(0.5+min(evmin(i  ),0.0)*dtdx)*dvm(4,i  )
         vr(5,i)=v(5,i  )-(0.5+min(evmin(i  ),0.0)*dtdx)*dvm(5,i  )
         vr(6,i)=bn(i  )
         vr(7,i)=v(6,i  )-(0.5+min(evmin(i  ),0.0)*dtdx)*dvm(6,i  )
         vr(8,i)=v(7,i  )-(0.5+min(evmin(i  ),0.0)*dtdx)*dvm(7,i  )
      enddo

#ifdef CHARACTERISTIC_TRACING
!
!-----construct time-centered left state for riemann solver
!
      dvd(:,:)=0.0
      do i = 3,n2-2
       k = 1
       do kn=1,nem1; if(eigv(k,i-1)<0.0 )exit
        do l=1,nem1
         dvd(l,i-1)=dvd(l,i-1)+
     &      (evmax(i-1)-eigv(k,i-1))*dw(k,i-1)*rhe(l,k,i-1)
        enddo
        k = k+1
       enddo
      enddo


      do i = 4,n2-2
         vl(1,i)= vl(1,i)+dtdx*dvd(1,i-1) 
         vl(2,i)= vl(2,i)+dtdx*dvd(2,i-1) 
         vl(3,i)= vl(3,i)+dtdx*dvd(3,i-1) 
         vl(4,i)= vl(4,i)+dtdx*dvd(4,i-1) 
         vl(5,i)= vl(5,i)+dtdx*dvd(5,i-1) 
         vl(6,i)= bn(i-1) 
         vl(7,i)= vl(7,i)+dtdx*dvd(6,i-1) 
         vl(8,i)= vl(8,i)+dtdx*dvd(7,i-1) 
      enddo
!
!-----construct time-centered right state for riemann solver
!
      dvd(:,:)=0.0
      do i = 3,n2-2
       k = nem1
       do kn=1,nem1; if(eigv(k,i)>0.0 )exit
        do l=1,nem1
         dvd(l,i)=dvd(l,i)+
     &      (evmin(i)-eigv(k,i))*dw(k,i)*rhe(l,k,i)
        enddo
        k = k-1
       enddo
      enddo

      do i = 4,n2-2
         vr(1,i)= vr(1,i)+dtdx*dvd(1,i) 
         vr(2,i)= vr(2,i)+dtdx*dvd(2,i) 
         vr(3,i)= vr(3,i)+dtdx*dvd(3,i) 
         vr(4,i)= vr(4,i)+dtdx*dvd(4,i) 
         vr(5,i)= vr(5,i)+dtdx*dvd(5,i) 
         vr(6,i)= bn(i) 
         vr(7,i)= vr(7,i)+dtdx*dvd(6,i) 
         vr(8,i)= vr(8,i)+dtdx*dvd(7,i) 
      enddo

#endif /*CHARACTERISTIC_TRACING*/


      return
      end subroutine plm_1d
