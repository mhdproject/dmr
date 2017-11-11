#include "athena_preprocessor.fh"
      subroutine flaten_slopes(vx,vy,vz,chi)
      use paramesh_dimensions
      use physcons
      use athena_specific_parameters

      use paramesh_interfaces
#ifdef MPI_USED
      use paramesh_mpi_interfaces
#endif /* MPI_USED */


      implicit none
      integer :: i,j,k
      real,dimension(ne,iu_bnd1,ju_bnd1,ku_bnd1),intent(in)::vx
      real,dimension(ne,iu_bnd1,ju_bnd1,ku_bnd1),intent(in)::vy
      real,dimension(ne,iu_bnd1,ju_bnd1,ku_bnd1),intent(in)::vz
      real,dimension(iu_bnd1,ju_bnd1,ku_bnd1),intent(out)::chi
      real,dimension(iu_bnd1,ju_bnd1,ku_bnd1)::chix
      real,dimension(iu_bnd1,ju_bnd1,ku_bnd1)::chiy
      real,dimension(iu_bnd1,ju_bnd1,ku_bnd1)::chiz
      real,parameter        :: small = 1.e-5
      real,parameter        :: smlrho= 1.e-3
      real,parameter        :: smallp= 1.e-3
      real,parameter        :: smallu= 1.e-6
      real,parameter        :: smalle= 1.e-6
      real,parameter        :: z0    = 0.25
      real,parameter        :: z1    = 0.75
      real,parameter        :: z2    = 0.75
      real,parameter        :: z3    = 0.85
      real :: slim,dsig,dsig2,kbar,z,fx,fy,fz,cs
!
!     flaten zone structure in regions where shocks are too thin
!
!     Code based on that of Greg Miller, LBNL
!
      chix(:,:,:) =1.0 ; chiy(:,:,:) =1.0 ; chiz(:,:,:) =1.0
      do k=1+2*k3d,ku_bnd1-2*k3d
       do j=3,ju_bnd1-2
        do i=3,iu_bnd1-2
          if( vx(2,i-1,j,k).lt.vx(2,i+1,j,k)) then
            chix(i,j,k) = 1.0
          else
            dsig = abs(vx(5,i+1,j,k)-vx(5,i-1,j,k))
            dsig2= abs(vx(5,i+2,j,k)-vx(5,i-2,j,k))
            if( dsig .gt. z3*dsig2 ) then
              slim = 0.0
            elseif(dsig.le.z2*dsig2.or.dsig2.le.small) then
              slim = 1.0
            else
              slim = max(0.0,min(1.0,
     .		(z3-dsig/dsig2)/(z3-z2) ))
            endif
            cs          = sqrt(gamma*vx(5,i,j,k)/vx(1,i,j,k))
            kbar        = vx(1,i,j,k)*cs**2
            z           = dsig/kbar
            chix(i,j,k) = (z1-z)/(z1-z0)
#ifdef DEBUG
          if(i.eq.5.and.j.eq.6)then
          write(*,*)'chix(i,j,k)=',chix(i,j,k)
	  write(*,*)'slim,z,kbar,dsig,dsig2=',slim,z,kbar,dsig,dsig2
	  write(*,*)'z1,z2,z3=',z1,z2,z3
          endif
#endif
            chix(i,j,k) = max(slim,min(1.0,chix(i,j,k)))
#ifdef DEBUG
          if(i.eq.5.and.j.eq.6)then
          write(*,*)'chix(i,j,k)=',chix(i,j,k)
          write(*,*)'vx(2,i-1,j,k),vx(2,i+1,j,k)=',
     &               vx(2,i-1,j,k),vx(2,i+1,j,k)
          endif
#endif
          endif

          if( vy(3,i,j-1,k).lt.vy(3,i,j+1,k)) then
            chiy(i,j,k) = 1.0
          else
            dsig = abs(vy(5,i,j+1,k)-vy(5,i,j-1,k))
            dsig2= abs(vy(5,i,j+2,k)-vy(5,i,j-2,k))
            if( dsig .gt. z3*dsig2 ) then
              slim = 0.0
            elseif(dsig.le.z2*dsig2.or.dsig2.le.small) then
              slim = 1.0
            else
              slim = max(0.0,min(1.0,
     .		(z3-dsig/dsig2)/(z3-z2) ))
            endif
            cs          = sqrt(gamma*vy(5,i,j,k)/vy(1,i,j,k))
            kbar        = vy(1,i,j,k)*cs**2
            z           = dsig/kbar
            chiy(i,j,k) = (z1-z)/(z1-z0)
            chiy(i,j,k) = max(slim,min(1.0,chiy(i,j,k)))
          endif
#if N_DIM == 3
          if( vz(4,i,j,k-1).lt.vz(4,i,j,k+1)) then
            chiz(i,j,k) = 1.0
          else
            dsig = abs(vz(5,i,j,k+1)-vz(5,i,j,k-1))
            dsig2= abs(vz(5,i,j,k+2)-vz(5,i,j,k-2))
            if( dsig .gt. z3*dsig2 ) then
              slim = 0.0
            elseif(dsig.le.z2*dsig2.or.dsig2.le.small) then
              slim = 1.0
            else
              slim = max(0.0,min(1.0,
     .		(z3-dsig/dsig2)/(z3-z2) ))
            endif
            cs          = sqrt(gamma*vz(5,i,j,k)/vz(1,i,j,k))
            kbar        = vz(1,i,j,k)*cs**2
            z           = dsig/kbar
            chiz(i,j,k) = (z1-z)/(z1-z0)
            chiz(i,j,k) = max(slim,min(1.0,chiz(i,j,k)))
          endif
!
!-----Determine final flatening coefficient
!
        enddo
       enddo
      enddo

!
#else
        enddo
       enddo
      enddo
#endif

#if N_DIM == 3
      chi(:,:,:)=1.0
      do k=3,ku_bnd1-2
       do j=3,ju_bnd1-2
        do i=3,iu_bnd1-2
          fx = min(chix(i-1,j  ,k  ),chix(i,j,k),chix(i+1,j  ,k  ))
          fy = min(chiy(i  ,j-1,k  ),chiy(i,j,k),chiy(i  ,j+1,k  ))
          fz = min(chiz(i  ,j  ,k-1),chiz(i,j,k),chiz(i  ,j  ,k+1))
          chi(i,j,k) = min(fx,fy,fz)
        enddo
       enddo
      enddo
#else
      chi(:,:,:)=1.0
      do k=1,ku_bnd1
       do j=3,ju_bnd1-2
        do i=3,iu_bnd1-2
          fx = min(chix(i-1,j  ,k),chix(i,j,k),chix(i+1,j  ,k))
          fy = min(chiy(i  ,j-1,k),chiy(i,j,k),chiy(i  ,j+1,k))
          chi(i,j,k) = min(fx,fy)
#ifdef DEBUG
          if(i.eq.6.and.j.eq.6)then
          write(*,*)'chix(i-1,j  ,k),chix(i,j,k),chix(i+1,j  ,k)=',
     &               chix(i-1,j  ,k),chix(i,j,k),chix(i+1,j  ,k)
          write(*,*)'chiy(i  ,j-1,k),chiy(i,j,k),chiy(i  ,j+1,k)=',
     &               chiy(i  ,j-1,k),chiy(i,j,k),chiy(i  ,j+1,k)
          endif
#endif
        enddo
       enddo
      enddo
#endif


      return
      end subroutine flaten_slopes
!========================================================================================
!
!========================================================================================
      function cvmgm(aa,bb,cc)
      real :: aa,bb,cc,cvmgm
      cvmgm=bb
      if(cc.lt.0.) cvmgm=aa
      return
      end function cvmgm
!========================================================================================
!
!========================================================================================
      function cvmgp(aa,bb,cc)
      real :: aa,bb,cc,cvmgp
      cvmgp=bb
      if(cc.ge.0.) cvmgp=aa
      return
      end function cvmgp
!========================================================================================
!
!========================================================================================
      function cvmgz(aa,bb,cc)
      real :: aa,bb,cc,cvmgz
      cvmgz=bb
      if(cc.eq.0.) cvmgz=aa
      return
      end function cvmgz
