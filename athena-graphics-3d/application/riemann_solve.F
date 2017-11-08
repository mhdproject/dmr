#include "paramesh_preprocessor.fh"
#include "athena_preprocessor.fh"
#include "mydebug.fh"
#undef BFILTER1
#undef DEBUG
!=======================================================================
!
!=======================================================================
      subroutine riemann_solve(ud,divv,vxfl,vxfr,fx,vf,idir,is,ir,lflag,lb)      
      use paramesh_dimensions
      use physcons
      use physicaldata
      use tree
      use strongfields
      use interior_gbc
      use grid
      use athena_specific_parameters
      use interfaces, only :  roe_mhd,hlle_mhd

      use paramesh_interfaces
#ifdef MPI_USED
      use paramesh_mpi_interfaces
#endif /* MPI_USED */


      implicit none
#define cell_v real,dimension(ne,iu_bnd1 ,ju_bnd1 ,ku_bnd1 )
#define cell_s real,dimension(   iu_bnd1 ,ju_bnd1 ,ku_bnd1 )
#define face_v real,dimension(ne,iu_bnd1p,ju_bnd1p,ku_bnd1p)
#define face_f real,dimension(nf,iu_bnd1p,ju_bnd1p,ku_bnd1p)

      logical,intent(in   )::lflag
      integer,intent(in   )::idir,is,lb,ir
      cell_s ,intent(in   )::divv
      cell_v ,intent(in   )::ud
      face_v ,intent(inout)::vxfl,vxfr
      face_f ,intent(  out)::fx  ! fluid fluxes
      face_v ,intent(  out)::vf
      real,dimension(ne,nmax+1)::vl,vr 
      real,dimension(ne,nmax  )::uc,ur,vfd 
      real,dimension( 3,nmax+1)::bsf
      real,dimension(   nmax  )::div_x,bn
      real,dimension(11,nmax  )::flo
      integer :: n1,n2,n3,i,j,k,nl1,nh1,nl2,nh2,nl3,nh3,jj
      integer :: l3d,l2d,l1d
#ifdef DEBUG
      integer :: il,iu,jl,ju,kl,ku,l
      real, dimension(nf)::vmax   ,vmin
      real, dimension(nf)::varmax ,varmin
      real, dimension(nf)::varmaxr,varminr
      real, dimension(nf)::varmaxl,varminl
#endif /*DEBUG*/

c       write(*,*)'nx,ny,nz=',nx,ny,nz
!
!-----x-direction
!
      select case (idir)
       case (1)
        n1 =nx   ;n2 =ny ;n3=nz
        l1d=k1d  ;l2d=k2d   ;l3d=k3d
        if(is.eq.1)then
         nl1=1+norder*l1d;nh1=n1-(norder-1)*l1d
         nl2=1+l2d       ;nh2=n2-l2d
         nl3=1+l3d       ;nh3=n3-l3d
         nl2=1           ;nh2=n2  
         nl3=1           ;nh3=n3  
        elseif(is.eq.2)then
         nl1=nlx  ;nh1=nhx 
         nl2=nly  ;nh2=nhym1
         nl3=nlz  ;nh3=nhzm1
        elseif(is.eq.3)then
         nl1=nlxp1;nh1=nhxm1
         nl2=nlyp1;nh2=nhym1
         nl3=nlzp1;nh3=nhzm1
        endif
       case (2)
        n1 =ny ;n2 =nz ;n3=nx
        l1d=k2d;l2d=k3d;l3d=k1d
        if(is.eq.1)then
         nl1=1+norder*l1d;nh1=n1-(norder-1)*l1d
         nl2=1+l2d       ;nh2=n2-l2d
         nl3=1+l3d       ;nh3=n3-l3d
         nl2=1           ;nh2=n2    
         nl3=1           ;nh3=n3    
        elseif(is.eq.2)then
         nl1=nly  ;nh1=nhy 
         nl2=nlz  ;nh2=nhzm1
         nl3=nlx  ;nh3=nhxm1
        elseif(is.eq.3)then
         nl1=nlyp1;nh1=nhym1
         nl2=nlzp1;nh2=nhzm1
         nl3=nlxp1;nh3=nhxm1
        endif
       case (3)
        n1 =nz ;n2 =nx ;n3=ny
        l1d=k3d;l2d=k1d;l3d=k2d
        if(is.eq.1)then
         nl1=1+norder*l1d;nh1=n1-(norder-1)*l1d
         nl2=1+l2d       ;nh2=n2-l2d
         nl3=1+l3d       ;nh3=n3-l3d
         nl2=1           ;nh2=n2  
         nl3=1           ;nh3=n3  
        elseif(is.eq.2)then
         nl1=nlz  ;nh1=nhz 
         nl2=nlx  ;nh2=nhxm1
         nl3=nly  ;nh3=nhym1
        elseif(is.eq.3)then
         nl1=nlzp1;nh1=nhzm1
         nl2=nlxp1;nh2=nhxm1
         nl3=nlyp1;nh3=nhym1
        endif
      end select
c       write(*,*)'nl1,nh1,nl2,nh2,nl3,nh3,idir = ',nl1,nh1,nl2,nh2,nl3,nh3,idir
       do k = nl3,nh3   ! z-loop
        do j = nl2,nh2  ! y-loop

         select case (idir)
         case (1)
         do i=nl1,nh1
          bsf(1,i) = bxf_sf(1,i,j,k)
          bsf(2,i) = bxf_sf(2,i,j,k)
          bsf(3,i) = bxf_sf(3,i,j,k)
         enddo
         do i=2,nx-1
          div_x(i) = divv(i,j,k)
         enddo
         do i=1,nx
          uc(:ne,i) = ud(:ne,i,j,k)
         enddo
         call rotate_u(uc,ur,1,nx,idir)
         do i=nl1,nh1
          vl(:ne,i) = vxfl(:ne,i,j,k)
          vr(:ne,i) = vxfr(:ne,i,j,k)
         enddo
#ifdef BFILTER1
         do i=nl1,nh1
          bn(i) = facevarx1(1,i,j,k,1)
         enddo 
#else
         do i=nl1,nh1
          bn(i) = 0.5*(vl(6,i)+vr(6,i))
         enddo
#endif 
 
#ifdef INTERIOR_BOUNDARY
        if(lflag)call internal_hardwall_bc(i,j,k,1)
#endif /*INTERIOR_BOUNDARY	*/ 
         case (2)
         do i=nl1,nh1
          bsf(1,i) = byf_sf(1,k,i,j)
          bsf(2,i) = byf_sf(2,k,i,j)
          bsf(3,i) = byf_sf(3,k,i,j)
         enddo
         do i=2,ny-1
          div_x(i) = divv(k,i,j)
         enddo
         do i=1,ny
          uc(:ne,i) = ud(:ne,k,i,j)
         enddo
         call rotate_u(uc,ur,1,ny,idir)
         do i=nl1,nh1
          vl(:ne,i) = vxfl(:ne,k,i,j)
          vr(:ne,i) = vxfr(:ne,k,i,j)
         enddo 
#ifdef BFILTER1
         do i=nl1,nh1
          bn(i) = facevary1(1,k,i,j,1)
         enddo 
#else
         do i=nl1,nh1
          bn(i) = 0.5*(vl(6,i)+vr(6,i))
         enddo
#endif 
#ifdef INTERIOR_BOUNDARY
         if(lflag)call internal_hardwall_bc(k,i,j,2)
#endif /*INTERIOR_BOUNDARY	*/ 
         case (3)
         do i=nl1,nh1
          bsf(1,i) = bzf_sf(1,j,k,i)
          bsf(2,i) = bzf_sf(2,j,k,i)
          bsf(3,i) = bzf_sf(3,j,k,i)
         enddo
         do i=2,nz-1
          div_x(i) = divv(j,k,i)
         enddo
         do i=1,nz
          uc(:ne,i) = ud(:ne,j,k,i)
         enddo
         call rotate_u(uc,ur,1,nz,idir)
         do i=nl1,nh1
          vl(:ne,i) = vxfl(:ne,j,k,i)
          vr(:ne,i) = vxfr(:ne,j,k,i)
         enddo 
#ifdef BFILTER1
         do i=nl1,nh1
          bn(i) = facevarz1(1,i,j,k,1)
         enddo 
#else
         do i=nl1,nh1
          bn(i) = 0.5*(vl(6,i)+vr(6,i))
         enddo
#endif 
c         if(lb.eq.73)write(*,*)'entering internal_hardwall_bc:3'
#ifdef INTERIOR_BOUNDARY
         if(lflag)call internal_hardwall_bc(j,k,i,3)
#endif /*INTERIOR_BOUNDARY	*/ 
         end select
!
!-----call riemann solver
!
#ifdef DEBUG
         if(lb.eq.73.and.idir.eq.3.and.k.eq.10.and.j.eq.11)then
          do i=10,12
	   do jj=1,ne
	    write(*,*)'vl(',jj,',',i,'),vr(',jj,',',i,')=',vl(jj,i),vr(jj,i)
	   enddo
	   do jj=1,ne
	    write(*,*)'ur(',jj,',',i,')=',ur(jj,i)
	   enddo
	   do jj=1,3
	    write(*,*)'bsf(',jj,',',i,')=',bsf(jj,i)
	   enddo

	   write(*,*)'bn(',i,'),div_x(',i,')=',bn(i),div_x(i)
          enddo

	 endif 
#endif 
  
c#ifdef HLLE
cc         if(lb.eq.73.and.idir.eq.3.and.k.eq.10.and.j.eq.11)write(*,*)'entering hlle solver:3,k,j=',k,j
c         call  hlle_mhd(nl1,nh1,vl,vr,ur,bn,bsf,div_x,flo,vfd,idir,lb)
cc         if(lb.eq.73.and.idir.eq.3.and.k.eq.10.and.j.eq.11)write(*,*)'leaving hlle solver:3'
c#endif
c#ifdef ROE
c         call   roe_mhd(nl1,nh1,vl,vr,ur,bn,bsf,div_x,flo,vfd,idir)
c#endif
c#ifdef NON_LINEAR_CFD
c         call   nls_cfd(nl1,nh1,vl,vr,ur,bn,bsf,div_x,flo,vfd,idir)
c#endif
c#ifdef EXACT
c         call exact_mhd(nl1,nh1,vl,vr,ur,bn,bsf,div_x,flo,vfd,idir)
c#endif
c       if(lb.eq.73.and.idir.eq.3.and.k.eq.10.and.j.eq.11)write(*,*)'leaving riemann solver:3'

        select case (ir)
	 case (1)
          call  hlle_mhd(nl1,nh1,vl,vr,ur,bn,bsf,div_x,flo,vfd,idir,lb)
	 case (2)
          call   roe_mhd(nl1,nh1,vl,vr,ur,bn,bsf,div_x,flo,vfd,idir)
	end select
 
        select case (idir)
         case (1)
         do i=nl1,nh1
          fx( :ne,i,j,k) = flo( :ne,i)
#ifdef LORENTZ_FORCE
          fx(9:11,i,j,k) = flo(9:11,i)
#else
          fx(9:11,i,j,k) = 0.0
#endif
          vf( :ne,i,j,k) = vfd( :ne,i)
         enddo
         case (2)
         do i=nl1,nh1
          fx( :ne,k,i,j) = flo( :ne,i)
#ifdef LORENTZ_FORCE
          fx(9:11,k,i,j) = flo(9:11,i)
#else
          fx(9:11,k,i,j) = 0.0
#endif
          vf( :ne,k,i,j) = vfd( :ne,i)
         enddo
         case (3)
          do i=nl1,nh1
          fx( :ne,j,k,i) = flo( :ne,i)
#ifdef LORENTZ_FORCE
          fx(9:11,j,k,i) = flo(9:11,i)
#else
          fx(9:11,j,k,i) = 0.0
#endif
          vf( :ne,j,k,i) = vfd( :ne,i)
         enddo
         end select

       enddo  ! end of y-loop over j
      enddo  ! end of z-loop over k
#ifdef DEBUG
      if(lb.eq.73.and.idir.eq.3)then

      write(*,*)'nl1,nh1,nl2,nh2,nl3,nh3=',nl1,nh1,nl2,nh2,nl3,nh3
      select case (idir)
      case(1)
      do l=1,nf
       varmax(l) = maxval(fx(l,nl1:nh1,nl2:nh2,nl3:nh3) )
       varmin(l) = minval(fx(l,nl1:nh1,nl2:nh2,nl3:nh3) )
      enddo
      do l=1,nf
      write(*,*)'fx(',l,'),fx(',l,') = ',
     $               varmax(l),varmin(l)
      enddo
       do l=1,nf
       write(*,*)'fx(',l,'),i,j,k=',l,maxloc(fx(l,nl1:nh1,nl2:nh2,nl3:nh3))
       write(*,*)'fx(',l,'),i,j,k=',l,minloc(fx(l,nl1:nh1,nl2:nh2,nl3:nh3))
       enddo
      case(2)
c      write(*,*)'entering case 2'

      do l=1,nf
       varmax(l) = maxval(fx(l,nl2:nh2,nl1:nh1,nl1:nh3) )
       varmin(l) = minval(fx(l,nl2:nh2,nl1:nh1,nl1:nh3) )
      enddo
      do l=1,nf
      write(*,*)'fy(',l,'),fy(',l,') = ',
     $               varmax(l),varmin(l)
      enddo
      write(*,*)'leaving case 2'

      case(3)

      do l=1,nf
       varmax(l) = maxval(fx(l,nl2:nh2,nl3:nh3,nl1:nh1) )
       varmin(l) = minval(fx(l,nl2:nh2,nl3:nh3,nl1:nh1) )
      enddo
      do l=1,nf
      write(*,*)'fz(',l,'),fz(',l,') = ',
     $               varmax(l),varmin(l)
      enddo
      end select
      endif
#endif /*DEBUG*/

      return
      end subroutine riemann_solve
