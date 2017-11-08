#include "athena_preprocessor.fh"
#define NEVER7 
      subroutine ppm_1d(vc,flat,bx,by,bz,wl,wr,n1,n2,dt,idir,lflag)
      use paramesh_dimensions
      use physcons
      use grid
      use strongfields
      use athena_specific_parameters
      implicit none
      integer :: i,j,k,n1,n2,l,kn,idir
      real    :: dt,dtdx,sl
      logical                         ,intent(in   )::lflag                        
      real,dimension(ne  ,     nmax  ),intent(in   )::vc
      real,dimension(ne  ,     nmax+1),intent(  out)::wl,wr
      real,dimension(nem1,     nmax  )              ::vn,vp
      real,dimension(          nmax  ),intent(in   )::flat
      real,dimension(          nmax  ),intent(in   )::bx,by,bz
      real,dimension(nem1,     nmax  )              ::vd,dvl,dvr,vl,vr
      real,dimension(nem1,     nmax  )              ::v,u
      real,dimension(          nmax  )              ::bn,cfast,bnd
      real,dimension(nem1,     nmax  )              ::vave
      real,dimension(nem1,     nmax  )              ::a,b,c
      real,dimension(nem1,     nmax  )              ::eigv,dw
      real,dimension(   4,     nmax  )              ::speeds
      real,dimension(nem1,nem1,nmax  )              ::rhe ,lhe
      real,dimension(nem1,nem1,nmax  )              ::vlt,vrt,dvt
      real,dimension(nem1,nem1,nmax  )              ::rhec,lhec
      real,dimension(nem1,     nmax  )              ::alpha_e ,alpha_lim 
      real,dimension(nem1,     nmax  )              ::eta_l,eta_r 
      real,dimension(nem1,     nmax  )              ::alph ,beta 
      real,dimension(nem1,     nmax  )              ::bsum,bdif 
      real,dimension(nem1,     nmax  )              ::dwe,dvd
      real,dimension(          nmax  )              ::evmin,evmax
      real,dimension(          nmax  )              ::dvy,dvz
      real,parameter                                ::third = 1.0/3.0
      real,parameter                                ::sixth = 1.0/6.0
      real,parameter                                ::twth  = 1.0/12.0
      real,parameter                                ::svth  = 7.0/12.0
      real,dimension(ne          )                  ::varmax,varmin
      real :: sum1,sum2,c1,c2,sgx,dum1,dum2,abslim
!
!-----first rotate conserved variables 
!
      if(idir.eq.1)then     ! x-direction
       dtdx = dt*rdx
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
      elseif(idir.eq.2)then ! y-direction
       dtdx = dt*rdy
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
      else                  ! z-direction
       dtdx = dt*rdz
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
      endif
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

#ifdef HIGH_ORDER
!
!-----compute eigenvectors and eigenvalues to construct reference states
!
      dvd=0.
      call mhd_eigenvectors(rhe ,lhe ,rhec,lhec,eigv,vd   ,bnd    ,
     &                      dw  ,dvd ,1   ,1   ,n2  ,cfast,speeds ,
     &                      dvy, dvz ,idir)
!
!-----first define the 4th order zone average at i+1/2: eq(2)
!
      do i=2,n2-2
       do j=1,nem1
        vave(j,i) = svth*(v(j,i  )+v(j,i+1))
     &            - twth*(v(j,i+2)+v(j,i-1))
       enddo
      enddo
      do i=3,n2-1
       do j=1,nem1
        vn(j,i) = vave(j,i-1)
       enddo
      enddo
      do i=2,n2-2
       do j=1,nem1
        vp(j,i) = vave(j,i  )
       enddo
      enddo
c
c-----constrain face values to be between adjacent cell values
c
      do i=3,n2-1
       do j=1,nem1
        vn(j,i) = max(min(v(j,i),v(j,i-1)),vn(j,i))
        vn(j,i) = min(max(v(j,i),v(j,i-1)),vn(j,i))
       enddo
      enddo
      do i=2,n2-2
       do j=1,nem1
        vp(j,i) = max(min(v(j,i),v(j,i+1)),vp(j,i))
        vp(j,i) = min(max(v(j,i),v(j,i+1)),vp(j,i))
       enddo
      enddo
c
c-----Insure face values allow a monotone parabolic profile
c
      do i=3,n2-2
       do j=1,nem1
        if( (vp(j,i)-v(j,i))*(v(j,i)-vn(j,i)) <= 0.0  )then
         vn(j,i) = v(j,i)
         vp(j,i) = v(j,i)
        endif
       enddo
      enddo
c
c-----Construct parabolic profile coefficients
c
      do i=3,n2-2
       do j=1,nem1
        a(j,i) = 3.0*(vn(j,i)+vp(j,i)-2.0*v(j,i))
        b(j,i) = vp(j,i)-vn(j,i)
        c(j,i) = v(j,i)-twth*a(j,i)
       enddo
      enddo

      do i=3,n2-2
       do j=1,nem1
        if((a(j,i)+b(j,i))*b(j,i)<0.0)then
         vn(j,i) = v(j,i)-2.0*(vp(j,i)-v(j,i))
        endif
        if((b(j,i)-a(j,i))*b(j,i)<0.0)then
         vp(j,i) = v(j,i)+2.0*(v(j,i)-vn(j,i))
        endif
       enddo
      enddo
#ifdef NEVER7
c
c-----monotonize face quantities vl and vr using characteristics
c
      do i=3,n2-2
        dwe(:,i)=(vp(:,i)-0.5*(v(:,i)+v(:,i+1)))
      enddo       
      do i=3,n2-2
       do k=1,nem1
        alpha_e(k,i) = dot_product(lhe(k,:,i),dwe(:,i))
       enddo
      enddo
      do i=3,n2-2
        dwe(:,i)=(v(:,i+1)-v(:,i))
      enddo       
      do i=3,n2-2
       do k=1, nem1
        alpha_lim(k,i) = dot_product(lhe(k,:,i),dwe(:,i))
       enddo
      enddo

      do i=3,n2-2
       do k=1,nem1
        abslim      = abs(alpha_lim(k,i))
        sgx         =sign(1.,alpha_lim(k,i)) 
        alpha_e(k,i)=min(0.5*abslim,max(-0.5*abslim,
     &               alpha_e(k,i)))
       enddo
      enddo

      do i=3,n2-2
       do j=1,nem1
        sum1 = 0.
        do k=1,nem1
         sum1=sum1+alpha_e(k,i)*rhe(j,k,i)
        enddo
         vp(j,i)=0.5*(v(j,i)+v(j,i+1))+sum1 ! first place vr(4).ne.0
       enddo
      enddo

      do i=3,n2-2
       do j=1,nem1
        vn(j,i)=vp(j,i-1)
       enddo
      enddo

c      return
c
c-----construct coefficients for parabolic interpolant
c       
      do i=3,n2-2
       do j=1,nem1
        a(j,i)=3.0*(vp(j,i)+vn(j,i)-2.0*v(j,i))
        b(j,i)=(vp(j,i)-vn(j,i))
        c(j,i)=v(j,i)-a(j,i)*twth
       enddo
      enddo       
     
      do i=3,n2-2
       do k=1, nem1
        alph(k,i) = dot_product(lhe(k,:,i),a(:,i))
        beta(k,i) = dot_product(lhe(k,:,i),b(:,i))
       enddo
      enddo

      do i=3,n2-2
       do k=1, nem1
        eta_l(k,i) = dot_product(lhe(k,:,i),( v(:,i)-vn(:,i)))
        eta_r(k,i) = dot_product(lhe(k,:,i),(vp(:,i)- v(:,i)))
        bsum(k,i)  = beta(k,i)+alph(k,i)
        bdif(k,i)  = beta(k,i)-alph(k,i)
       enddo
      enddo

      do i=3,n2-2
       do k=1, nem1
        if(eta_l(k,i)*eta_r(k,i)<0) then
         eta_l(k,i)=0.0
         eta_r(k,i)=0.0
        endif
        if(bsum(k,i)*bdif(k,i)<0) then
         if(beta(k,i)*bsum(k,i)<0) eta_l(k,i)=-2.0*eta_r(k,i)
         if(beta(k,i)*bdif(k,i)<0) eta_r(k,i)= 2.0*eta_l(k,i)
        endif
       enddo
      enddo

      do i=3,n2-2
       do j=1,nem1
        sum1 = 0.0
        sum2 = 0.0
        do k=1,nem1
         sum1=sum1+eta_l(k,i)*rhe(j,k,i)
         sum2=sum2+eta_r(k,i)*rhe(j,k,i)
        enddo
         vn(j,i)=v(j,i)-sum1
         vp(j,i)=v(j,i)+sum2
       enddo
      enddo
!      return
c      write(*,*)'passing 1'
#endif /*NEVER7*/
!
!-----apply flattening
!
      do i=3,n2-2
       do j=1,nem1
         vn(j,i)=(1.0-flat(i))*v(j,i)+flat(i)*vn(j,i)
         vp(j,i)=(1.0-flat(i))*v(j,i)+flat(i)*vp(j,i)
       enddo
      enddo
c
c-----construct new coefficients for parabolic interpolant
c  
      do i=3,n2-2
       do j=1,nem1
        a(j,i)=3.0*(vp(j,i)+vn(j,i)-2.0*v(j,i))
        b(j,i)=(vp(j,i)-vn(j,i))
        c(j,i)=v(j,i)-a(j,i)*twth
       enddo
      enddo
c      write(*,*)'passing 2'


      do i=3,n2-2
       do k=1,nem1
         if(eigv(k,i)>0)then
          c1     = eigv(k,i)*dtdx
          c2     = c1*c1
          do j=1,nem1
           vlt(k,j,i)=a(j,i)*(0.25-0.5*c1+third*c2)
     &               +0.5*b(j,i)*(1-c1)+c(j,i)
          enddo
         else
          do j=1,nem1
           vlt(k,j,i)=0.25*a(j,i)+0.5*b(j,i)+c(j,i)
          enddo
         endif
        enddo
      enddo
       
c      write(*,*)'passing 3'

      do i=3,n2-3
       do k=1,nem1
        if(eigv(k,i+1)<0)then
         c1     = eigv(k,i+1)*dtdx
         c2     = c1*c1
         do j=1,nem1
          vrt(k,j,i)=a(j,i+1)*(0.25+0.5*c1+third*c2)-0.5*b(j,i+1)*(1+c1)+c(j,i+1)
         enddo
        else
         do j=1,nem1
          vrt(k,j,i)=0.25*a(j,i+1)-0.5*b(j,i+1)+c(j,i+1)
         enddo
        endif
       enddo
      enddo
c
c-----perform characteristic tracing
c 
      do i=3,n2-3
       do k=1,nem1
        do j=1,nem1
         dvt(k,j,i)=vlt(k,j,i)-vlt(nem1,j,i)
        enddo
       enddo
      enddo

c      write(*,*)'passing 5'
      

      do i=3,n2-3
       do k=1, nem1
        dvl(k,i) = dot_product(lhe(k,:,i),dvt(k,:,i))
       enddo
      enddo


      do i=3,n2-3
       do j=1,nem1
        sum1 = 0.0
        do k=1,nem1
         if(eigv(k,i)>0)then
          sum1= sum1+dvl(k,i)*rhe(j,k,i)
         endif
        enddo
         vl(j,i)= vlt(nem1,j,i)+sum1
       enddo
      enddo
c      write(*,*)'passing 7'

      do i=3,n2-3
       do k=1,nem1
        do j=1,nem1
         dvt(k,j,i)=vrt(k,j,i)-vrt(1,j,i)
        enddo
       enddo
      enddo 
c      write(*,*)'passing 8'

      do i=3,n2-3
       do k=1, nem1
        dvr(k,i) = dot_product(lhe(k,:,i+1),dvt(k,:,i))
       enddo
      enddo

     
      do i=3,n2-3
       do j=1,nem1
        sum2 = 0.0
        do k=1,nem1
         if(eigv(k,i+1)<0)then
          sum2= sum2+dvr(k,i)*rhe(j,k,i+1)
         endif
        enddo
        vr(j,i)=vrt(1,j,i)+sum2
       enddo
      enddo

      do i = 3,n2-3
         wl(1,i)= vl(1,i-1)
         wl(2,i)= vl(2,i-1)
         wl(3,i)= vl(3,i-1) 
         wl(4,i)= vl(4,i-1)
         wl(5,i)= vl(5,i-1) 
         wl(6,i)= bn(i-1) 
         wl(7,i)= vl(6,i-1) 
         wl(8,i)= vl(7,i-1)

         wr(1,i)= vr(1,i-1)
         wr(2,i)= vr(2,i-1)
         wr(3,i)= vr(3,i-1) 
         wr(4,i)= vr(4,i-1)
         wr(5,i)= vr(5,i-1) 
         wr(6,i)= bn(i) 
         wr(7,i)= vr(6,i-1) 
         wr(8,i)= vr(7,i-1)
      enddo
#else 
      do i = 3,n2-3
         wl(1,i)= v(1,i-1)
         wl(2,i)= v(2,i-1)
         wl(3,i)= v(3,i-1) 
         wl(4,i)= v(4,i-1)
         wl(5,i)= v(5,i-1) 
         wl(6,i)= bn(i-1) 
         wl(7,i)= v(6,i-1) 
         wl(8,i)= v(7,i-1)

         wr(1,i)= v(1,i)
         wr(2,i)= v(2,i)
         wr(3,i)= v(3,i) 
         wr(4,i)= v(4,i)
         wr(5,i)= v(5,i) 
         wr(6,i)= bn(i) 
         wr(7,i)= v(6,i) 
         wr(8,i)= v(7,i)
      enddo

#endif  /*HIGH_ORDER*/      
        
      
      return
      end subroutine ppm_1d
