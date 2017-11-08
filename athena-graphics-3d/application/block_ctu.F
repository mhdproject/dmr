#undef DEBUG
#include "paramesh_preprocessor.fh"
#include "athena_preprocessor.fh"
#include "mydebug.fh"
#undef STONE_FIX
#define CTU_CORRECTION

       subroutine block_ctu(v1    ,u1    ,
     &                      fx    ,fy    ,fz    ,
     &                      ex    ,ey    ,ez    ,
     &                      du    ,dl    ,
     &                      dt    ,time  ,lb    ,
     &                      lflag ,istep ,mype  )
!=======================================================================
      use paramesh_dimensions
      use physcons
      use physicaldata
      use tree
      use strongfields
      use interior_gbc
      use grid
      use athena_specific_parameters
      use interfaces, only :riemann_solve                 ,
     &                      reconstruction                ,
     &                      correct_rl_states_2d          ,
     &                      correct_rl_states_3d_xa       ,      
     &                      correct_rl_states_3d_ya       ,      
     &                      correct_rl_states_3d_za       ,      
     &                      correct_rl_states_3d_xb       ,      
     &                      correct_rl_states_3d_yb       ,      
     &                      correct_rl_states_3d_zb       ,
     &                      brackbill                     ,
     &                      stonefix                       
c     &                      interior_bc_fl                
c     &                      interior_bc_states_nplus_1half,
c     &                      interior_bc_states_n        

      use paramesh_interfaces
#ifdef MPI_USED
      use paramesh_mpi_interfaces
#endif /* MPI_USED */


      implicit none
!--------------------------------------------------------------
! include files for amr


#include "amr_shmem.fh"
#include "shmem_reduce.fh"
#define face_x  real,dimension(ne,iu_bnd1p,ju_bnd1p,ku_bnd1p)
#define face_v  real,dimension(ne,iu_bnd1p,ju_bnd1p,ku_bnd1p)
#define face_f  real,dimension(nf,iu_bnd1p,ju_bnd1p,ku_bnd1p)
#define cell_s  real,dimension(   iu_bnd1 ,ju_bnd1 ,ku_bnd1 )
#define cell_sp real,dimension(   iu_bnd1p,ju_bnd1p,ku_bnd1p)
#define cell_v  real,dimension(ne,iu_bnd1 ,ju_bnd1 ,ku_bnd1 )
#define cell_l  real,dimension( 3,iu_bnd1 ,ju_bnd1 ,ku_bnd1 )

      cell_v,intent(   in) :: u1
      cell_v,intent(   in) :: v1
      face_f,intent(  out) :: fx
      face_f,intent(  out) :: fy
      face_f,intent(  out) :: fz
      cell_v,intent(  out) :: du
      cell_l,intent(  out) :: dl
      real, dimension(il_bnd1:iu_bnd1p,jl_bnd1:ju_bnd1p,
     &                kl_bnd1:ku_bnd1p),intent(out) ::ex,ey,ez
      logical, intent(in   )::lflag 
      integer, intent(in   )::lb,istep ,mype 
      real   , intent(in   )::dt,time  

!-------------------------------------------------------------
!mk9
      real,dimension(ne,iu_bnd1p,ju_bnd1p ,ku_bnd1p ) :: xmy_vl
      real,dimension(ne,iu_bnd1p,ju_bnd1p ,ku_bnd1p ) :: xmy_vr
      real,dimension(ne,iu_bnd1p,ju_bnd1p ,ku_bnd1p ) :: ymy_vl 
      real,dimension(ne,iu_bnd1p,ju_bnd1p ,ku_bnd1p ) :: ymy_vr
      real,dimension( 3,iu_bnd1p,ju_bnd1p ,ku_bnd1p ) :: bxf_my
      real,dimension( 3,iu_bnd1p,ju_bnd1p ,ku_bnd1p ) :: byf_my
!-------------------------------------------------------------

      face_x::vfx
      face_x::vfy
      face_x::vfz

      integer ::  i,j,k,l,m,n
      integer :: il,iu,jl,ju,kl,ku,ii
!-----------------------------------------------------------------------
      cell_v               :: vx
      cell_v               :: vy
      cell_v               :: vz
      cell_v               :: ud
      cell_sp              :: xmask
      cell_sp              :: ymask
      cell_sp              :: zmask

!-----------------------------------------------------------------------
      cell_v::src,deu,dux,duy,duz
      cell_s::divv,chi
      face_x::vxfl ,vxfr
      face_x::vyfl ,vyfr
      face_x::vzfl ,vzfr
      face_x::vxyfl,vxyfr
      face_x::vxzfl,vxzfr
      face_x::vyxfl,vyxfr
      face_x::vyzfl,vyzfr
      face_x::vzxfl,vzxfr
      face_x::vzyfl,vzyfr
      face_f::fxz
      face_f::fxy
      face_f::fyx
      face_f::fyz
      face_f::fzx
      face_f::fzy

!
!-----------------------------------------------------------------------
!
!-----1d scratch arrays
!
      real, dimension(   nmax  )::div_x,div_y,div_z,flat
      real, dimension(   nmax  )::mas_x,mas_y,mas_z
      real, dimension(   nmax  )::bxs,bys,bzs,bn
      real, dimension(nf,nmax  )::flo
      real, dimension(ne,nmax  )::vfd,vc,uc,ur,vd,vdr
      real, dimension(ne,nmax+1)::vl,vr
!-----------------------------------------------------------------------
      real :: bx,by,bz,bxt,byt,bzt,bs
c      integer ::imin,imax ,jmin,jmax ,kmin,kmax	
      integer ::ilow,ihigh,jlow,jhigh,klow,khigh	
!
!
!
!------rotation indexes
!
       integer,dimension(ne+3)::index_x ,index_y ,index_z
       parameter (index_x =(/1,2,3,4,5,6,7,8, 9,10,11/))
       parameter (index_y =(/1,3,4,2,5,7,8,6,10,11, 9/))
       parameter (index_z =(/1,4,2,3,5,8,6,7,11, 9,10/))
!
!-----------------------------------------------------------------------
      real :: fmaxx,fminx
      real, dimension(nf)          :: varmax,varmin
      integer, parameter :: mx = max(iu_bnd1+1,ju_bnd1+k2d,ku_bnd1+k3d)
      real, dimension( 3,nmax     ) :: bsf
      integer :: ix,iy,iz,f,in
      integer ::iup1,jup1,kup1
      integer :: mylb,lb1,lb2
      real, parameter :: onethird = 1.0/3.0
      real            :: dtdx,dtdy,dtdz
!---------------------------------------------------------------------------
!
!-----x-direction from left supersonic inflow 
!
      real :: sl, rh, ph, eh
      real :: b2,bv
      real :: rlb,vxlb,vylb,vzlb,plb,elb,bxlb,bylb,bzlb
      common/left_state/rlb,vxlb,vylb,vzlb,plb,bxlb,bylb,bzlb
      real :: rr,bx1,by1,bz1,bx0,by0,bz0,ek,bk,et,px,py,pz
      integer, dimension(ne) :: imax,jmax,kmax
      integer, dimension(ne) :: imin,jmin,kmin
      real   , dimension(ne) :: fmax,fmin
      real :: dth,dbx,dby,dbz,vybx,vzbx,vxby,vzby,vxbz,vybz
      real :: bxby,byby,bzby,bzbz,vyby,bzbx,bybx,ri,bxbz,bybz,vzbz

!---------------------------------------------------------------------------
c        write(*,*)'within  block_ctu'
c
c-------initialize
c
        vc    = 0.0
        uc    = 0.0
        ur    = 0.0
        vl    = 0.0
        vr    = 0.0
        bsf   = 0.0
        div_x = 0.0
        mas_x = 0.0
        div_y = 0.0
        mas_y = 0.0
        div_z = 0.0
        mas_z = 0.0
        xmask = 1.0
        ymask = 1.0
        zmask = 1.0

#ifdef DEBUG
      if(lb.eq.73)then
      do j=1,ne
       varmax(j) = maxval(v1(j,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,
     &                          kl_bnd1:ku_bnd1) )
       varmin(j) = minval(v1(j,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,
     &                         kl_bnd1:ku_bnd1) )
      enddo
      do j=1,ne
c      write(*,*)'v1(',j,'),v1(',j,') = ',
     $               varmax(j),varmin(j)
      enddo
      endif
#endif /*DEBUG*/

      vx = v1
      vy = v1
      vz = v1
#ifdef INTERIOR_BOUNDARY
      if(lflag)then
       call interior_bc_states_n(vx,1,lb)
       call interior_bc_states_n(vy,2,lb)
#if N_DIM == 3
       call interior_bc_states_n(vz,3,lb)
#endif
      endif 
#endif /*INTERIOR_BOUNDARY*/

#ifdef FLATEN_SLOPES
!
!----Flatten slopes
!
      call flaten_slopes(vx,vy,vz,chi)
#else
      chi = 1.0
#endif /*FLATEN_SLOPES*/

       vxfl(:,:,:,:)=0.0
       vxfr(:,:,:,:)=0.0
       vyfl(:,:,:,:)=0.0
       vyfr(:,:,:,:)=0.0
#if N_DIM == 3
       vzfl(:,:,:,:)=0.0
       vzfr(:,:,:,:)=0.0
#endif /*N_DIM == 3*/
#if N_DIM == 3
      dtdx        =dt*rdx
      dtdy        =dt*rdy
      dtdz        =dt*rdz
#else
      dtdx        =dt*rdx
      dtdy        =dt*rdy
      dtdz        =0.0
#endif
#ifdef STONE_FIX
       call stonefix(vx,vy,vz,vxfl,vxfr,vyfl,vyfr,vzfl,vzfr,dt)
#endif /*STONE_FIX*/

!
!=============Compute First Approximation to Reconstructed left and Right States
!
!
!-----x-direction
!
c      if(lb.eq.73)write(*,*)'entrying reconstruction-x'
      call reconstruction(vx,chi,vxfl,vxfr,dt,1,lflag)
!
!-----y-direction
!
c      if(lb.eq.73)write(*,*)'entrying reconstruction-y'
      call reconstruction(vy,chi,vyfl,vyfr,dt,2,lflag)
#if N_DIM == 3
!
!-----z-direction
!
c      if(lb.eq.73)write(*,*)'entrying reconstruction-z'
      call reconstruction(vz,chi,vzfl,vzfr,dt,3,lflag)
#endif /*N_DIM == 3*/

!=======================================================================
      ud = u1
#ifdef INTERIOR_BOUNDARY
      if(lflag)then
       call interior_bc_states_n(ud,1,lb)
      endif 
#endif /*INTERIOR_BOUNDARY*/
      divv(:,:,:)= 0.0
      call div_V(vx,divv,1)
c     if(lb.eq.73)write(*,*)'entrying riemann_solve-x'
      call riemann_solve(ud,divv,vxfl,vxfr,fx,vfx,1,1,1,lflag,lb)
!
!----y-direction
!
      ud = u1
#ifdef INTERIOR_BOUNDARY
      call interior_bc_states_n(ud,2,lb)
#endif /*INTERIOR_BOUNDARY*/
      call div_V(vy,divv,2)
c      if(lb.eq.73)write(*,*)'entrying riemann_solve-y'
      call riemann_solve(ud,divv,vyfl,vyfr,fy,vfy,2,1,1,lflag,lb)
#if N_DIM == 3
!
!----z-direction
!
      ud = u1
#ifdef INTERIOR_BOUNDARY
      call interior_bc_states_n(ud,3,lb)
#endif /*INTERIOR_BOUNDARY*/
      call div_V(vz,divv,3)
c      if(lb.eq.73)write(*,*)'entrying riemann_solve-z'
      call riemann_solve(ud,divv,vzfl,vzfr,fz,vfz,3,1,1,lflag,lb)
c      write(*,*)'leaving riemann_solve-z'
#endif /*N_DIM == 3*/
      if(.not.lflag)then
#ifdef CTU_CORRECTION
!
!-----apply CTU correction
!
!
!-----correct x-face right and left states for transverse corrections
!
#if N_DIM == 3
      dtdx        =onethird*dt*rdx
      dtdy        =onethird*dt*rdy
      dtdz        =onethird*dt*rdz
#else
      dtdx        =0.5*dt*rdx
      dtdy        =0.5*dt*rdy
      dtdz        =0.0
#endif

#ifdef INTERIOR_BOUNDARY
      if(lflag)then
       call interior_bc_fl(fx,fy,fz,vfx,vfy,vfz,v1,lb)
      endif
#endif /*INTERIOR_BOUNDARY*/
   
#if N_DIM == 2
c      write(*,*)'entrying correct_rl_states_2d'
      call correct_rl_states_2d(vxfl   ,vxfr   ,vyfl   ,vyfr   ,
     &                          fx     ,fy     ,
     &                          index_x,index_y,
     &                          dtdx   ,dtdy   )
c      write(*,*)'leaving correct_rl_states_2d'
#endif /*N_DIM == 2*/

#if N_DIM == 3
!
!-----correct right and left states for transverse corrections
!
c      if(lb.eq.73)write(*,*)'entrying correct_rl_states_3d_xa'
      call correct_rl_states_3d_xa(vxfl   ,vxfr   ,
     &                             vxyfl  ,vxyfr  ,vxzfl  ,vxzfr  ,
     &                             fy     ,fz     ,index_x,
     &                             dtdy   ,dtdz   )
c      write(*,*)'leaving correct_rl_states_3d_xa'
      ud = u1
#ifdef INTERIOR_BOUNDARY
      if(lflag)then
       call interior_bc_states_n(ud,1,lb)
      endif
#endif /*INTERIOR_BOUNDARY*/

      call div_V(vx,divv,1)
c      if(lb.eq.73)write(*,*)'entrying riemann_solve ax1:'
      call riemann_solve(ud,divv,vxyfl,vxyfr,fxy,vfx,1,2,1,lflag,lb)
c      if(lb.eq.73)write(*,*)'entrying riemann_solve ax2:'
      call riemann_solve(ud,divv,vxzfl,vxzfr,fxz,vfx,1,2,1,lflag,lb)
!
!-----y-direction
!
c      if(lb.eq.73)write(*,*)'entrying correct_rl_states_3d_ya:'
      call correct_rl_states_3d_ya(vyfl   ,vyfr   ,
     &                             vyxfl  ,vyxfr  ,vyzfl  ,vyzfr  ,
     &                             fx     ,fz     ,index_y,
     &                             dtdx   ,dtdz   )
!
!-----y-direction
!
      ud = u1
#ifdef INTERIOR_BOUNDARY
      if(lflag)then
       call interior_bc_states_n(ud,2,lb)
      endif
#endif /*INTERIOR_BOUNDARY*/

      call div_V(vy,divv,2)
c     if(lb.eq.73)write(*,*)'entrying riemann_solve ay1:'
      call riemann_solve(ud,divv,vyxfl,vyxfr,fyx,vfy,2,2,1,lflag,lb)
c     if(lb.eq.73)write(*,*)'entrying riemann_solve ay2:'
      call riemann_solve(ud,divv,vyzfl,vyzfr,fyz,vfy,2,2,1,lflag,lb)
!
!-----z-direction
!
c      if(lb.eq.73)write(*,*)'entrying correct_rl_states_3d_za:'
      call correct_rl_states_3d_za(vzfl   ,vzfr   ,
     &                             vzxfl  ,vzxfr  ,vzyfl  ,vzyfr  ,
     &                             fx     ,fy     ,index_z,
     &                             dtdx   ,dtdy   )
!
!-----z-direction
!
      ud = u1
#ifdef INTERIOR_BOUNDARY
      if(lflag)then
       call interior_bc_states_n(ud,3,lb)
      endif
#endif /*INTERIOR_BOUNDARY*/

      call div_V(vz,divv,3)
c      if(lb.eq.73)write(*,*)'entrying riemann_solve az1:'
      call riemann_solve(ud,divv,vzxfl,vzxfr,fzx,vfz,3,2,1,lflag,lb)
c      if(lb.eq.73)write(*,*)'entrying riemann_solve az2:'
      call riemann_solve(ud,divv,vzyfl,vzyfr,fzy,vfz,3,2,1,lflag,lb)

      dtdx        =0.5*dt*rdx
      dtdy        =0.5*dt*rdy
      dtdz        =0.5*dt*rdz

!
!-----final transverse corrections to right and left states.
!
#ifdef INTERIOR_BOUNDARY
      if(lflag)then
       call interior_bc_fl(fxy,fyx,fzx,vfx,vfy,vfz,v1,lb)
       call interior_bc_fl(fxz,fyz,fzy,vfx,vfy,vfz,v1,lb)
      endif
#endif /*INTERIOR_BOUNDARY*/
c      if(lb.eq.73)write(*,*)'entrying correct_rl_states_3d_xb:'
      call correct_rl_states_3d_xb(vxfl   ,vxfr   ,
     &                             fyz    ,fzy    ,index_x,
     &                             dtdy   ,dtdz   )
      call interior_bc_states_nplus_1half(vxfl,vxfr,1,lb)
!
!-----y-direction
!
c      if(lb.eq.73)write(*,*)'entrying correct_rl_states_3d_yb:'
      call correct_rl_states_3d_yb(vyfl   ,vyfr   ,
     &                             fxz    ,fzx    ,index_y,
     &                             dtdx   ,dtdz   )
      call interior_bc_states_nplus_1half(vyfl,vyfr,2,lb)
!
!-----z-direction
!
c      if(lb.eq.73)write(*,*)'entrying correct_rl_states_3d_zb:'
      call correct_rl_states_3d_zb(vzfl   ,vzfr   ,
     &                             fxy    ,fyx    ,index_z,
     &                             dtdx   ,dtdy   )     
      call interior_bc_states_nplus_1half(vzfl,vzfr,3,lb)
#endif /*N_DIM == 3*/
!=======================================================================
!             Compute the Final Fluxes
!=======================================================================
      ud = u1
#ifdef INTERIOR_BOUNDARY
      if(lflag)then      
       call interior_bc_states_n(ud,1,lb)
      endif
#endif /*INTERIOR_BOUNDARY*/
       call div_V(vx,divv,1)
c      if(lb.eq.73)write(*,*)'entrying riemann_solve bx:'
       call riemann_solve(ud,divv,vxfl,vxfr,fx,vfx,1,3,2,lflag,lb)

      ud = u1
#ifdef INTERIOR_BOUNDARY
      if(lflag)then      
       call interior_bc_states_n(ud,2,lb)
      endif
#endif /*INTERIOR_BOUNDARY*/
      call div_V(vy,divv,2)
c      if(lb.eq.73)write(*,*)'entrying riemann_solve by:'
      call riemann_solve(ud,divv,vyfl,vyfr,fy,vfy,2,3,2,lflag,lb)

#if N_DIM == 3
!
!-----call riemann solver
!
      ud = u1
#ifdef INTERIOR_BOUNDARY
      if(lflag)then      
       call interior_bc_states_n(ud,3,lb)
      endif
#endif /*INTERIOR_BOUNDARY*/

      call div_V(vz,divv,3)
c      if(lb.eq.73)write(*,*)'entrying riemann_solve bz:lflag=',lflag
      call riemann_solve(ud,divv,vzfl,vzfr,fz,vfz,3,3,2,lflag,lb)
c      if(lb.eq.73)write(*,*)'leaving riemann_solve bz:'
#endif /*N_DIM == 3*/
#endif /*CTU_CORRECTION*/
       endif ! if(.not.lflag)
!==============================================================================================
!                   Set Boundary Conditions at Interior Boundary
!==============================================================================================
#ifdef INTERIOR_BOUNDARY
!
!-----set interior boundary conditions on the fluxes
!
      if(lflag)then
c        if(lb.eq.73)write(*,*)'entering interior_bc_fl'

       call interior_bc_fl(fx,fy,fz,vfx,vfy,vfz,v1,lb)
        if(lb.eq.73)write(*,*)'leaving interior_bc_fl'
      endif
#endif /*INTERIOR_BOUNDARY*/
!==============================================================================================
!                   Compute Changes to Conserved Varables
!==============================================================================================
      dux=0.
      duy=0.
      duz=0.
      du =0.

#ifdef DEBUG
      if(lb.eq.73)write(*,*)'entering DEBUG after CTU_CORRECTION'
            if(lb.eq.73)then
      write(*,*)'ne,nlxp1,nhxm2,nlyp1,nhym2,nlzp1,nhzm2=',
     &           ne,nlxp1,nhxm2,nlyp1,nhym2,nlzp1,nhzm2
      do k=nlzp1,nhzm2
       do j=nlyp1,nhym1
        do i=nlxp1,nhxm1
         do l=1,ne
          write(*,*)'fx(',l,',',i,',',j,',',k,')=',fx(l,i,j,k)
         enddo
        enddo
       enddo
      enddo

      do l=1,ne
       varmax(l) = maxval(fx(l,nlxp1:nhxm2,nlyp1:nhym2,nlzp1:nhzm2) )
       varmin(l) = minval(fx(l,nlxp1:nhxm1,nlyp1:nhym1,nlzp1:nhzm1) )
      enddo
      do l=1,ne
      write(*,*)'fx(',l,'),fx(',l,') = ',
     $               varmax(l),varmin(l)
      enddo
      do l=1,ne
       varmax(l) = maxval(fy(l,nlxp1:nhxm2,nlyp1:nhym2,nlzp1:nhzm2) )
       varmin(l) = minval(fy(l,nlxp1:nhxm1,nlyp1:nhym1,nlzp1:nhzm1) )
      enddo
      do l=1,ne
      write(*,*)'fy(',l,'),fy(',l,') = ',
     $               varmax(l),varmin(l)
      enddo
      do k=nlzp1,nhzm2
       do j=nlyp1,nhym1
        do i=nlxp1,nhxm1
          if(i.eq.5.and.j.eq.5)then
             write(*,*)'i,j,k        =',i,j,k
             write(*,*)'duu(1,i  ,j  ,k)=',du(1,i  ,j  ,k)
             write(*,*)'dux(1,i  ,j  ,k)=',dux(1,i  ,j  ,k)
             write(*,*)'duy(1,i  ,j  ,k)=',duy(1,i  ,j  ,k)
             write(*,*)' v1(1,i  ,j  ,k)=',v1(1,i  ,j  ,k)
             write(*,*)' v1(2,i  ,j  ,k)=',v1(2,i  ,j  ,k)
             write(*,*)' v1(3,i  ,j  ,k)=',v1(3,i  ,j  ,k)
             write(*,*)' v1(4,i  ,j  ,k)=',v1(4,i  ,j  ,k)
             write(*,*)' v1(5,i  ,j  ,k)=',v1(5,i  ,j  ,k)
             write(*,*)' v1(1,i+1,j  ,k)=',v1(1,i+1,j  ,k)
             write(*,*)' v1(2,i+1,j  ,k)=',v1(2,i+1,j  ,k)
             write(*,*)' v1(3,i+1,j  ,k)=',v1(3,i+1,j  ,k)
             write(*,*)' v1(4,i+1,j  ,k)=',v1(4,i+1,j  ,k)
             write(*,*)' v1(5,i+1,j  ,k)=',v1(5,i+1,j  ,k)
             write(*,*)' fx(1,i  ,j  ,k)=',fx(1,i  ,j  ,k)
             write(*,*)' fx(1,i+1,j  ,k)=',fx(1,i+1,j  ,k)
             write(*,*)' fy(1,i  ,j  ,k)=',fy(1,i  ,j  ,k)
             write(*,*)' fy(1,i  ,j+1,k)=',fy(1,i  ,j+1,k)
          endif
        enddo
       enddo
      enddo
#if N_DIM == 3

      do j=1,ne
       varmax(j) = maxval(fz(j,nlxp1:nhxm2,nlyp1:nhym2,nlzp1:nhzm2) )
       varmin(j) = minval(fz(j,nlxp1:nhxm1,nlyp1:nhym1,nlzp1:nhzm1) )
      enddo
      do j=1,ne
      write(*,*)'fz(',j,'),fz(',j,') = ',
     $               varmax(j),varmin(j)
      enddo
#endif
         endif
#endif /*DEBUG*/

      do k=nlzp1,nhzm2
       do j=nlyp1,nhym2
        do i=nlxp1,nhxm2
         do l=1,ne
          dux(l,i,j,k)= (fx(l,i+1,j  ,k  )-fx(l,i,j,k))*rvx
          duy(l,i,j,k)= (fy(l,i  ,j+1,k  )-fy(l,i,j,k))*rvy
#if N_DIM == 3         
          duz(l,i,j,k)= (fz(l,i  ,j  ,k+1)-fz(l,i,j,k))*rvz
           du(l,i,j,k)=-(dux(l,i,j,k)+duy(l,i,j,k)+duz(l,i,j,k))
#else
          duz(l,i,j,k)= 0.0
           du(l,i,j,k)=-(dux(l,i,j,k)+duy(l,i,j,k))
#endif
         enddo
        enddo
       enddo
      enddo


!mylb=1 for 1blk and >1 for multi-block
!for multi-block
!#ifdef multidebug
!       mylb=3
!if(lb.eq.3.and.istep.ge.5900) 
!       if(lb.eq.3) 
!    & call ch7_data(fx,fy,fz,xmy_vl,xmy_vr,ymy_vl,ymy_vr,
!    &           bxf_my, byf_my, mylb) 
!#endif
!#ifdef onedebug
!for 1block
!      mylb=1
!if(istep.ge.5900) 
!      call ch7_data(fx,fy,fz,xmy_vl,xmy_vr,ymy_vl,ymy_vr,
!    &           bxf_my, byf_my, mylb) 
!#endif



#ifdef BRACKBILL_SRC_TERM
       call brackbill(v1,src,dt)
#ifdef INTERIOR_BOUNDARY
      if(lflag)then
        do iz=1,izones
         i = izone(1,iz)
         j = izone(2,iz)
         k = izone(3,iz)
         scr(:ne,i,j,k) = 0.0
        enddo
       endif
#endif /*INTERIOR_BOUNDARY*/
      do k=1+nguard*k3d,nzb+nguard*k3d
         do j=1+nguard*k2d,nyb+nguard*k2d
            do i=1+nguard  ,nxb+nguard
             du(:ne,i,j,k) = du(:ne,i,j,k)+src(:ne,i,j,k) 
            enddo
         enddo
      enddo
#endif /*BRACKBILL_SRC_TERM*/

!=======================================================================
#ifdef INTERIOR_BOUNDARY
      if(lflag)then
       do iz=1,izones
         i = izone(1,iz)
         j = izone(2,iz)
         k = izone(3,iz)
         du(:8,i,j,k) = 0.0
       enddo
      endif
#endif
!=======================================================================
#ifdef LORENTZ_FORCE
c           if(lb.eq.73)write(*,*)'beginning lorentz force'

!
!-----now project out any components of Lorentz force produced
!     by discretization parallel to B. First compute deltas produced by
!     the Lorentz force
!
      dl(:,:,:,:)=0.0
      do k=nlzp1,nhzm2
       do j=nlyp1,nhym2
        do i=nlxp1,nhxm2
         do l=1,3
          m=ne+l
          dux(l,i,j,k)=-(fx(m,i+1,j  ,k  )-fx(m,i,j,k))*rvx
          duy(l,i,j,k)=-(fy(m,i  ,j+1,k  )-fy(m,i,j,k))*rvy
#if N_DIM == 3
          duz(l,i,j,k)=-(fz(m,i  ,j  ,k+1)-fz(m,i,j,k))*rvz
           dl(l,i,j,k)= (dux(l,i,j,k)+duy(l,i,j,k)+duz(l,i,j,k))
#else
          duz(l,i,j,k)= 0.0
           dl(l,i,j,k)= (dux(l,i,j,k)+duy(l,i,j,k))          
#endif /*N_DIM == 3*/
!=======================================================================
         enddo
        enddo   
       enddo   
      enddo
!=======================================================================
#ifdef INTERIOR_BOUNDARY
      if(lflag)then
        do iz=1,izones
         i = izone(1,iz)
         j = izone(2,iz)
         k = izone(3,iz)
         dl(1:3,i,j,k) = 0.0
        enddo
      endif

#endif
!=======================================================================

#endif /*LORENTZ_FORCE*/

!=======================================================================
!
!-----D i v e r g e n c e  F r e e  M a g n e t i c  F i e l d  T r a n s p o r t
!                          t o  R o u n d - O f f
!
!=======================================================================
#ifdef BFILTER
c           if(lb.eq.73)write(*,*)'beginning bfilter'
#ifdef BALSARA_SPICER
!
!----B A L S A R A  A N D  S P I C E R  F I E L D  T R A N S P O R T
!
#if N_DIM == 3
          do k= 1+nguard*k3d, nzb+(1+nguard)*k3d
            do j= 1+nguard*k2d, nyb+(1+nguard)*k2d
              do i= 1+nguard   , nxb+nguard   
                ex(i,j,k)=0.25*(fz(7,i,j,k)+fz(7,i,j-k2d,k    )
     $                        - fy(8,i,j,k)-fy(8,i,j    ,k-k3d))
              enddo
            enddo
          enddo
          do k= 1+nguard*k3d, nzb+(1+nguard)*k3d
            do j= 1+nguard*k2d, nyb+nguard*k2d
              do i= 1+nguard   , nxb+(1+nguard)   
                ey(i,j,k)=0.25*(fx(8,i,j,k  )+fx(8,i  ,j,k-k3d)
     $                        - fz(6,i,j,k  )-fz(6,i-1,j,k    ))
              enddo
            enddo
          enddo
          do k= 1+nguard*k3d, nzb+nguard*k3d
            do j= 1+nguard*k2d, nyb+(1+nguard)*k2d
              do i= 1+nguard   , nxb+(1+nguard)   
                ez(i,j,k)=0.25*(fy(6,i,j,k)+fy(6,i-1,j    ,k)
     $                        - fx(7,i,j,k)-fx(7,i  ,j-k2d,k))
              end do
            end do
          end do
#else   
          do k= 1+nguard*k3d, nzb+(1+nguard)*k3d
            do j= 1+nguard*k2d, nyb+(1+nguard)*k2d
              do i= 1+nguard   , nxb+nguard   
                ex(i,j,k)=-0.5*(fy(8,i,j,k)+fy(8,i,j,k-k3d))
              enddo
            enddo
          enddo
          do k= 1+nguard*k3d, nzb+(1+nguard)*k3d
            do j= 1+nguard*k2d, nyb+nguard*k2d
              do i= 1+nguard   , nxb+(1+nguard)   
                ey(i,j,k)=0.5*(fx(8,i,j,k  )+fx(8,i,j,k-k3d))
              enddo
            enddo
          enddo
          do k= 1+nguard*k3d, nzb+nguard*k3d
            do j= 1+nguard*k2d, nyb+(1+nguard)*k2d
              do i= 1+nguard   , nxb+(1+nguard)   
                ez(i,j,k)=0.25*(fy(6,i,j,k)+fy(6,i-1,j    ,k)
     $                        - fx(7,i,j,k)-fx(7,i  ,j-k2d,k))
              end do
            end do
          end do
#ifdef DEBUGN
          if(lb.eq.73)then
           i = 9+nguard
           j = 1+nguard
           k = 1
           write(*,*)'lb,ez(i,j,k)                =',lb,ez(i,j,k)
           write(*,*)'lb,fy(6,i,j,k),fy(6,i-1,j,k)=',lb,
     &                   fy(6,i,j,k),fy(6,i-1,j,k)
           write(*,*)'lb,fx(7,i,j,k),fx(7,i,j-1,k)=',lb,
     &                   fx(7,i,j,k),fx(7,i,j-1,k)
          endif
         
#endif /*DEBUGN*/ 
#endif
#endif /*BALSARA_SPICER*/
#ifdef INTERIOR_BOUNDARY
#ifndef IONOSPHERE
      if(lflag)then
      do in =1, izones
        i    = izone(1,in)
        j    = izone(2,in)
#if N_DIM == 3
        k    = izone(3,in)
        ex(i    ,j:j+1,k:k+1) = 0.0
        ey(i:i+1,j    ,k:k+1) = 0.0
#else
        k    = 1
        ex(i    ,j:j+1,k    ) = 0.0
        ey(i:i+1,j    ,k    ) = 0.0
#endif
        ez(i:i+1,j:j+1,k    ) = 0.0
      enddo
      endif
#endif /*IONOSPHERE*/
#endif /*INTERIOR_BOUNDARY*/

#ifdef EXTERNAL_BOUNDARY_INFLOW
        if(neigh(1,1,lb).eq.-30)then
#if N_DIM == 3

c          do k= 1+nguard*k3d, nzb+(1+nguard)*k3d
c            do j= 1+nguard*k2d, nyb+(1+nguard)*k2d
c              do i= 1,nguard 
c                bya      = 0.5(byf_sf(1,i,j,k)+byf_sf(1,i,j,k-1))  
c                bza      = 0.5(byf_sf(2,i,j,k)+byf_sf(2,i,j,k-1))  
c                ex(i,j,k)=-(vylb*(Bzlb-bza)-vzlb*(Bylb-bya))
c              enddo
c            enddo
c          enddo
          do k= 1+nguard*k3d, nzb+(1+nguard)*k3d
            do j= 1+nguard*k2d, nyb+nguard*k2d
              do i= 1,nguard+1 
                bxa      = 0.5(bxf_sf(1,i,j,k)+bxf_sf(1,i,j,k-1  ))  
                bza      = 0.5(bxf_sf(3,i,j,k)+bxf_sf(3,i,j,k-1  ))  
                ey(i,j,k)=-(vzlb*(Bxlb-bxa)-vxlb*(Bzlb-bza))
              enddo
            enddo
          enddo
          do k= 1+nguard*k3d, nzb+nguard*k3d
            do j= 1+nguard*k2d, nyb+(1+nguard)*k2d
              do i= 1,nguard+1  
                bxa      = 0.5(bxf_sf(1,i,j,k)+bxf_sf(1,i,j-1,k  ))  
                bya      = 0.5(bxf_sf(2,i,j,k)+bxf_sf(2,i,j-1,k  ))  
                ez(i,j,k)=-(vxlb*(Bylb-bya)-vylb*(Bxlb-bxa))
              end do
            end do
          end do
#else   
          do k= 1+nguard*k3d, nzb+(1+nguard)*k3d
            do j= 1+nguard*k2d, nyb+(1+nguard)*k2d
              do i= 1,nguard  
                ex(i,j,k)=-(vylb*Bzlb-vzlb*Bylb)
              enddo
            enddo
          enddo
          do k= 1+nguard*k3d, nzb+(1+nguard)*k3d
            do j= 1+nguard*k2d, nyb+nguard*k2d
              do i= 1,nguard   
                ey(i,j,k)=-(vzlb*Bxlb-vxlb*Bzlb)
              enddo
            enddo
          enddo
          do k= 1+nguard*k3d, nzb+nguard*k3d
            do j= 1+nguard*k2d, nyb+(1+nguard)*k2d
              do i= 1,nguard   
                ez(i,j,k)=-(vxlb*Bylb-vylb*Bxlb)
              end do
            end do
          end do
#endif /*EXTERNAL_BOUNDARY_INFLOW*/

        endif
#endif
#if N_DIM == 2
        ex(:,:,:)=0.0
        ey(:,:,:)=0.0
#endif
#else
        ex(:,:,:)=0.0
        ey(:,:,:)=0.0
        ez(:,:,:)=0.0

#endif /*BFILTER*/

      do k=nlzp1,nhzm2
       do j=nlyp1,nhym2
        do i=nlxp1,nhxm2
         do l=1,ne
          ud(l,i,j,k)=u1(l,i,j,k)+du(l,i,j,k)
!=======================================================================
         enddo
        enddo   
       enddo   
      enddo


c      if(lb.eq.73)write(*,*)'leaving block_ctu'
      return
      end subroutine block_ctu
!=======================================================================
!
!=======================================================================
      subroutine stonefix(vx,vy,vz,vxfl,vxfr,vyfl,vyfr,vzfl,vzfr,dt)      
      use paramesh_dimensions
      use physcons
      use physicaldata
      use tree
      use strongfields
      use interior_gbc
      use grid
      use athena_specific_parameters
      use paramesh_interfaces
#ifdef MPI_USED
      use paramesh_mpi_interfaces
#endif /* MPI_USED */

      implicit none
      real  ,intent(in   )::dt
      cell_v,intent(in   )::vx
      cell_v,intent(in   )::vy
      cell_v,intent(in   )::vz
      face_x,intent(inout)::vxfl ,vxfr
      face_x,intent(inout)::vyfl ,vyfr
      face_x,intent(inout)::vzfl ,vzfr
      integer :: n1,n2,n3,i,j,k,nl1,nh1
      real    :: dth,dbx,dby,dbz,vybx,vzbx,vxby,vzby,vxbz,vybz
      real    :: bxby,byby,bzby,vyby,bzbx,bybx,vzbz,bybz,bzbz,bxbz
      real    :: ri,bxt,byt,bzt,dtdx,dtdy,dtdz

#if N_DIM == 3
      dtdx        =dt*rdx
      dtdy        =dt*rdy
      dtdz        =dt*rdz
#else
      dtdx        =dt*rdx
      dtdy        =dt*rdy
      dtdz        =0.0
#endif
       dth = 0.25*dtdx
       do k = 1, nz   ! z-loop
        do j = 1, ny  ! y-loop
         do i = nlx,nhx ! x-loop
          ri            = 1.0/max(vx(1,i,j,k),rmin)
          byt             = b_sf(2,i,j,k)+vy(7,i,j,k)
          bzt             = b_sf(3,i,j,k)+vz(8,i,j,k)
          dbx             = (vx(6,i+1,j,k)-vx(6,i-1,j,k))
c          dbx             = facevarx1(1,i+1,j,k,1)-facevarx1(1,i,j,k,1)
          vybx            = dth*vy(3,i,j,k)*dbx
          vzbx            = dth*vz(4,i,j,k)*dbx
          bybx            = dth*byt*dbx*pie4i*ri
          bzbx            = dth*bzt*dbx*pie4i*ri

          vxfr(3,i  ,j,k) = vxfr(3,i  ,j,k)+bybx
          vxfr(4,i  ,j,k) = vxfr(4,i  ,j,k)+bzbx
          vxfl(3,i+1,j,k) = vxfl(3,i+1,j,k)+bybx
          vxfl(4,i+1,j,k) = vxfl(4,i+1,j,k)+bzbx

          vxfr(7,i  ,j,k) = vxfr(7,i  ,j,k)+vybx
          vxfr(8,i  ,j,k) = vxfr(8,i  ,j,k)+vzbx
          vxfl(7,i+1,j,k) = vxfl(7,i+1,j,k)+vybx
          vxfl(8,i+1,j,k) = vxfl(8,i+1,j,k)+vzbx
         enddo 
        enddo  ! end of y-loop over j
       enddo  ! end of z-loop over k

       dth = 0.25*dtdy
       do k = 1, nz   ! z-loop
        do i = 1, nx    ! x-loop
         do j = nly,nhy  ! y-loop
          ri              = 1.0/max(vy(1,i,j,k),rmin)

          bxt             = b_sf(1,i,j,k)+vx(6,i,j,k)
          bzt             = b_sf(3,i,j,k)+vz(8,i,j,k)
          dby             = (vy(7,i,j+1,k)-vy(7,i,j-1,k))
c          dby             = facevary1(1,i,j+1,k,1)-facevarx1(1,i,j,k,1)
          vxby            = dth*vx(2,i,j,k)*dby
          vyby            = 0.0
          vzby            = dth*vz(4,i,j,k)*dby
          bxby            = dth*bxt*dby*pie4i*ri
          byby            = 0.0
          bzby            = dth*bzt*dby*pie4i*ri

c         vyfr(2,i,j  ,k) = vyfr(2,i,j  ,k)+byby
          vyfr(3,i,j  ,k) = vyfr(3,i,j  ,k)+bzby
          vyfr(4,i,j  ,k) = vyfr(4,i,j  ,k)+bxby
c         vyfl(2,i,j+1,k) = vyfl(2,i,j+1,k)+byby
          vyfl(3,i,j+1,k) = vyfl(3,i,j+1,k)+bzby
          vyfl(4,i,j+1,k) = vyfl(4,i,j+1,k)+bxby

c         vyfr(6,i,j  ,k) = vyfr(6,i,j  ,k)+vyby
          vyfr(7,i,j  ,k) = vyfr(7,i,j  ,k)+vzby
          vyfr(8,i,j  ,k) = vyfr(8,i,j  ,k)+vxby
c         vyfl(6,i,j+1,k) = vyfl(6,i,j+1,k)+vyby
          vyfl(7,i,j+1,k) = vyfl(7,i,j+1,k)+vzby
          vyfl(8,i,j+1,k) = vyfl(8,i,j+1,k)+vxby
         enddo 
        enddo  ! end of y-loop over j
       enddo  ! end of z-loop over k       
#if N_DIM == 3
       dth = 0.25*dtdz
       do i = 1, nx   ! x-loop
        do j = 1, ny   ! y-loop
         do k = nlz,nhz ! z-loop
          ri              = 1.0/max(vz(1,i,j,k),rmin)

          bxt             = b_sf(1,i,j,k)+vx(6,i,j,k)
          byt             = b_sf(2,i,j,k)+vy(7,i,j,k)
          dbz             = (vz(8,i,j,k+1)-vz(8,i,j,k-1))
c          dbz             = facevarz1(1,i,j,k+1,1)-facevarz1(1,i,j,k,1)
          vxbz            = dth*vx(2,i,j,k)*dbz
          vybz            = dth*vy(3,i,j,k)*dbz
          vzbz            = 0.0
          bxbz            = dth*bxt*dbz*pie4i*ri
          bybz            = dth*byt*dbz*pie4i*ri
          bzbz            = 0.0

          vzfr(2,i,j,k  ) = vzfr(2,i,j,k  )+bxbz
          vzfr(3,i,j,k  ) = vzfr(3,i,j,k  )+bybz
c         vzfr(4,i,j,k  ) = vzfr(4,i,j,k  )+bzbz
          vzfl(2,i,j,k+1) = vzfl(2,i,j,k+1)+bxbz
          vzfl(3,i,j,k+1) = vzfl(3,i,j,k+1)+bybz
c         vzfl(4,i,j,k+1) = vzfl(4,i,j,k+1)+bzbz

          vzfr(6,i,j,k  ) = vzfr(6,i,j,k  )+vzbz
          vzfr(7,i,j,k  ) = vzfr(7,i,j,k  )+vxbz
c         vzfr(8,i,j,k  ) = vzfr(8,i,j,k  )+vybz
          vzfl(6,i,j,k+1) = vzfl(6,i,j,k+1)+vzbz
          vzfl(7,i,j,k+1) = vzfl(7,i,j,k+1)+vxbz
c         vzfl(8,i,j,k+1) = vzfl(8,i,j,k+1)+vybz

         enddo 
        enddo  ! end of y-loop over j
       enddo  ! end of z-loop over k       
#endif /*N_DIM == 3*/
      return
      end subroutine stonefix
