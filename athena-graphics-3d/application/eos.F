#include "athena_preprocessor.fh"
!#define DEBUG
!#define DEBUGG
c----------------------------------------------------------------------
      subroutine eos_u_to_v_pmn(u,v,cf,lb,mype)
      use physcons
      use paramesh_dimensions
      use physicaldata
      use tree
      use strongfields
      use athena_specific_parameters

      use paramesh_interfaces
#ifdef MPI_USED
      use paramesh_mpi_interfaces
#endif /* MPI_USED */

      integer                              ,intent(in):: lb,mype
      real, dimension(nvar,il_bnd1:iu_bnd1,
     .                   jl_bnd1:ju_bnd1,
     .                   kl_bnd1:ku_bnd1,2),intent(in)::u
      real, dimension(ne,il_bnd1:iu_bnd1,
     .                   jl_bnd1:ju_bnd1,
     .                   kl_bnd1:ku_bnd1),intent(out)::v
      real, dimension(   il_bnd1:iu_bnd1,
     .                   jl_bnd1:ju_bnd1,
     .                   kl_bnd1:ku_bnd1),intent(out)::cf
      logical :: lprint
c

       do k = kl_bnd1,ku_bnd1
        do j = jl_bnd1,ju_bnd1
         do i = il_bnd1,iu_bnd1
          rr           = u(1,i,j,k,1)
          ri           = 1.0/rr
          px           = u(2,i,j,k,1)
          py           = u(3,i,j,k,1)
          pz           = u(4,i,j,k,1)
          et           = u(5,i,j,k,1)
          bx           = u(6,i,j,k,1)
          by           = u(7,i,j,k,1)
          bz           = u(8,i,j,k,1)
          b2           = bx*bx+by*by+bz*bz
          bxt          = bx+b_sf(1,i,j,k)
          byt          = by+b_sf(2,i,j,k)
          bzt          = bz+b_sf(3,i,j,k)
          bsq          = bxt*bxt+byt*byt+bzt*bzt
          rbi          = 1.0/(rr + afact*bsq)
          bv           = (px*bxt+py*byt+pz*bzt)*rbi

          vx           = (px + afact*bv*bxt)*rbi
          vy           = (py + afact*bv*byt)*rbi
          vz           = (pz + afact*bv*bzt)*rbi
          vsq          = vx*vx+vy*vy+vz*vz 

          ek           = 0.5*rr*vsq
          b2           = pie8i*b2
          pt           = et - ek - b2 
          p            = max(gammam1*pt,pmin)
          csq          = gamma*p*ri
          vasq         = 2.0*bsq*rbi
          csq          = sqrt(csq+vasq)
          v(1,i,j,k)   = rr
          v(2,i,j,k)   = vx
          v(3,i,j,k)   = vy
          v(4,i,j,k)   = vz
          v(5,i,j,k)   = p
          v(6,i,j,k)   = bx
          v(7,i,j,k)   = by
          v(8,i,j,k)   = bz
          cf(i,j,k)    = csq
        enddo
       enddo
      enddo
      return
      end subroutine eos_u_to_v_pmn
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
      subroutine eos_u_to_v(u,v,lb)
      use physcons
      use paramesh_dimensions
      use physicaldata
      use strongfields
      use athena_specific_parameters
      integer :: lb
      real, dimension(ne,il_bnd1:iu_bnd1,
     .                   jl_bnd1:ju_bnd1,
     .                   kl_bnd1:ku_bnd1)::u,v

      real, dimension(ne)::varmax,varmin
      real :: bxt,byt,bzt,rbi,bv,vsq,bsq,b2
      real :: bx,by,bz,vx,vy,vz,p,pt,et,rr,px,py,pz
#ifdef DEBUG
      do j = 1, 3
        varmax(j) = maxval( b_sf(j,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,
     &                                             kl_bnd1:ku_bnd1))
        varmin(j) = minval( b_sf(j,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,
     &                                             kl_bnd1:ku_bnd1))
      enddo
      do j=1,3
       write(*,*)'b_sf(',j,'),b_sf(',j,') = ',
     $               varmax(j),varmin(j)
      enddo

       do j=1,3
       write(*,*)'b_sf(',j,'),i,j,k=',j,maxloc(b_sf(j,il_bnd1:iu_bnd1,
     &                                                jl_bnd1:ju_bnd1,
     &                                                kl_bnd1:ku_bnd1))
       write(*,*)'b_sf(',j,'),i,j,k=',j,minloc(b_sf(j,il_bnd1:iu_bnd1,
     &                                                jl_bnd1:ju_bnd1,
     &                                                kl_bnd1:ku_bnd1))
       enddo
#endif /*DEBUG*/
c
       do k = kl_bnd1,ku_bnd1
        do j = jl_bnd1,ju_bnd1
         do i = il_bnd1,iu_bnd1
          rr           = u(1,i,j,k)
          ri           = 1.0/rr
          px           = u(2,i,j,k)
          py           = u(3,i,j,k)
          pz           = u(4,i,j,k)
          et           = u(5,i,j,k)
          bx           = u(6,i,j,k)
          by           = u(7,i,j,k)
          bz           = u(8,i,j,k)
          b2           = bx*bx+by*by+bz*bz
          bxt          = bx+b_sf(1,i,j,k)
          byt          = by+b_sf(2,i,j,k)
          bzt          = bz+b_sf(3,i,j,k)
          bsq          = bxt*bxt+byt*byt+bzt*bzt
#ifdef BORIS_CORRECTION
          rbi          = 1.0/(rr + afact*bsq)
          bv           = (px*bxt+py*byt+pz*bzt)*rbi
          vx           = (px + afact*bv*bxt)*rbi
          vy           = (py + afact*bv*byt)*rbi
          vz           = (pz + afact*bv*bzt)*rbi
#else
          vx           = px*ri
          vy           = py*ri
          vz           = pz*ri
          rbi          = ri
#endif /*BORIS_CORRECTION*/
          vsq          = vx*vx+vy*vy+vz*vz 
          ek           = 0.5*rr*vsq
          b2           = pie8i*b2
          pt           = et - ek - b2
c          pt           = et - ek 
          p            = max(gammam1*pt,pmin)
          v(1,i,j,k)   = rr
          v(2,i,j,k)   = vx
          v(3,i,j,k)   = vy
          v(4,i,j,k)   = vz
          v(5,i,j,k)   = p
#ifdef CFD
          v(6,i,j,k)   = 0.0
          v(7,i,j,k)   = 0.0
          v(8,i,j,k)   = 0.0
#else
          v(6,i,j,k)   = u(6,i,j,k)
          v(7,i,j,k)   = u(7,i,j,k)
          v(8,i,j,k)   = u(8,i,j,k)
#endif /*CFD*/
#ifdef DEBUGG
          if(i.eq.3.and.j.eq.45)then
          write(*,*)'cell i,j,k   =',i,j,k
          write(*,*)'afact        =',afact
          write(*,*)'et,ek,b2     = ',et,ek,b2
          write(*,*)'pt           = ',pt
          write(*,*)'rr,vsq       = ',rr,vsq
          write(*,*)'u(1,i,j,k)   =',u(1,i,j,k)
          write(*,*)'u(2,i,j,k)   =',u(2,i,j,k)
          write(*,*)'u(3,i,j,k)   =',u(3,i,j,k)
          write(*,*)'u(4,i,j,k)   =',u(4,i,j,k)
          write(*,*)'u(5,i,j,k)   = ',u(5,i,j,k)
          write(*,*)'u(6,i,j,k)   =',u(6,i,j,k)
          write(*,*)'u(7,i,j,k)   =',u(7,i,j,k)
          write(*,*)'u(8,i,j,k)   =',u(8,i,j,k)
          write(*,*)'b_sf(1,i,j,k)=',b_sf(1,i,j,k)
          write(*,*)'b_sf(2,i,j,k)=',b_sf(2,i,j,k)
          write(*,*)'b_sf(3,i,j,k)=',b_sf(3,i,j,k)
          write(*,*)'bv,rbi,bsq   =',bv,rbi,bsq
          write(*,*)'px,py,pz     =',px,py,pz
          write(*,*)'vx,vy,vz     =',vx,vy,vz
          write(*,*)'bx,by,bz     =',bx,by,bz
          write(*,*)'bxt,byt,bzt  =',bxt,byt,bzt
          pause
          endif
#endif
        enddo
       enddo
      enddo
      return
      end subroutine eos_u_to_v
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
      subroutine eos_v_to_u (u,v)
      use physcons
      use paramesh_dimensions
      use physicaldata
      use strongfields
      use athena_specific_parameters
      real, dimension(ne,il_bnd1:iu_bnd1,
     .                   jl_bnd1:ju_bnd1,
     .                   kl_bnd1:ku_bnd1)::u,v
c
       do k = kl_bnd1,ku_bnd1
        do j = jl_bnd1,ju_bnd1
         do i = il_bnd1,iu_bnd1
          rr          = v(1,i,j,k)
          ri          = 1.0/rr
          vx          = v(2,i,j,k)
          vy          = v(3,i,j,k)
          vz          = v(4,i,j,k)
          bx          = v(6,i,j,k)
          by          = v(7,i,j,k)
          bz          = v(8,i,j,k)
          bsq         = bx*bx+by*by+bz*bz
          bxt         = bx+b_sf(1,i,j,k)
          byt         = by+b_sf(2,i,j,k)
          bzt         = bz+b_sf(3,i,j,k)
          p           = v(5,i,j,k)
          v2          = vx**2+vy**2+vz**2 
          b2          = bxt**2+byt**2+bzt**2 
          bv          = (bxt*vx+byt*vy+bzt*vz)
          ek          = 0.5*rr*v2
          et          = ek + p*gammam1i+bsq*pie8i
          u(1,i,j,k)  = rr
          u(2,i,j,k)  = rr*vx + afact*(b2*vx - bv*bxt)
          u(3,i,j,k)  = rr*vy + afact*(b2*vy - bv*byt)
          u(4,i,j,k)  = rr*vz + afact*(b2*vz - bv*bzt)
          u(5,i,j,k)  = et
          u(6,i,j,k)  = bx 
          u(7,i,j,k)  = by 
          u(8,i,j,k)  = bz 
        enddo
       enddo
      enddo
      return
      end subroutine eos_v_to_u
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
      subroutine eos_u_to_v_loc( u,v,lb)
      use physcons
      use paramesh_dimensions
      use physicaldata
      use tree
      use strongfields
      use athena_specific_parameters
      integer :: ii,jj,kk,lb,i,j,k
      real, dimension(ne,il_bnd:iu_bnd,jl_bnd:ju_bnd,
     .                   kl_bnd:ku_bnd    ) :: u,v
       do k = 1,nzb
          kk = nguard*k3d+k
        do j = 1,nyb
            jj = nguard*k2d+j
         do i = 1,nxb
              ii = nguard+i
          rr           = u(1,i,j,k)
          px           = u(2,i,j,k)
          py           = u(3,i,j,k)
          pz           = u(4,i,j,k)
          et           = u(5,i,j,k)
          bx           = u(6,i,j,k)
          by           = u(7,i,j,k)
          bz           = u(8,i,j,k)
          b2           = bx**2+by**2+bz**2
#ifdef BORIS_CORRECTION
          bxt          = bx+b_sf(1,ii,jj,kk)
          byt          = by+b_sf(2,ii,jj,kk)
          bzt          = bz+b_sf(3,ii,jj,kk)
          bsq          = bxt**2+byt**2+bzt**2
          rbi          = 1.0/(rr + afact*bsq)
          bv           = (px*bxt+py*byt+pz*bzt)*rbi
          vx           = (px + afact*bv*bxt)*rbi
          vy           = (py + afact*bv*byt)*rbi
          vz           = (pz + afact*bv*bzt)*rbi
#else
          ri           = 1.0/rr
          vx           = px*ri
          vy           = py*ri
          vz           = pz*ri
#endif
          vsq          = vx**2+vy**2+vz**2 
          ek           = 0.5*rr*vsq
          pt           = et - ek - pie8i*b2
c          pt           = et - ek 
          p            = max(gammam1*pt,pmin)
          v(1,i,j,k)   = rr
          v(2,i,j,k)   = vx
          v(3,i,j,k)   = vy
          v(4,i,j,k)   = vz
          v(5,i,j,k)   = p
          v(6,i,j,k)   = bx
          v(7,i,j,k)   = by
          v(8,i,j,k)   = bz

c          write(*,*)'eos:neg den:',v(1,21,31,1),u(1,21,31,1)
c          stop

        enddo
       enddo
      enddo
      return
      end subroutine eos_u_to_v_loc
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
      subroutine eos_vf_to_uf(vl,vr,nx1,nx2,ny1,ny2,nz1,nz2,idir)
      use physcons
      use paramesh_dimensions
      use strongfields
      use athena_specific_parameters
      integer :: nx1,nx2,ny1,ny2,nz1,nz2,idir
      real,dimension(ne,iu_bnd1p,ju_bnd1p,ku_bnd1p),
     &                                  intent(inout)::vl,vr
      real :: bxt,byt,bzt,bxs,bys,bzs
      real :: ux,uy,uz
       do k = nz1,nz2
        do j = ny1,ny2
         do i = nx1,nx2
           rr          = max(vl(1,i,j,k),rmin)
           ri          = 1.0/rr
           vx          = vl(2,i,j,k)
           vy          = vl(3,i,j,k)
           vz          = vl(4,i,j,k)
#ifdef CFD
           bx          = 0.0
           by          = 0.0
           bz          = 0.0
#else
           bx          = vl(6,i,j,k)
           by          = vl(7,i,j,k)
           bz          = vl(8,i,j,k)
#endif /*CFD*/
          bsq         = bx*bx+by*by+bz*bz
          select case (idir)
          case (1)
#ifdef CFD
           bxt         = 0.0
           byt         = 0.0
           bzt         = 0.0
#else
           bxt         = bx+bxf_sf(1,i,j,k)
           byt         = by+bxf_sf(2,i,j,k)
           bzt         = bz+bxf_sf(3,i,j,k)
#endif /*CFD*/
!---------------------------------------------------
          case (2)
#ifdef CFD
           bxt         = 0.0
           byt         = 0.0
           bzt         = 0.0
#else
           bxt         = bx+byf_sf(1,i,j,k)
           byt         = by+byf_sf(2,i,j,k)
           bzt         = bz+byf_sf(3,i,j,k)
#endif /*CFD*/
!---------------------------------------------------
          case (3)
#ifdef CFD
           bxt         = 0.0
           byt         = 0.0
           bzt         = 0.0
#else
           bxt         = bx+bzf_sf(1,i,j,k)
           byt         = by+bzf_sf(2,i,j,k)
           bzt         = bz+bzf_sf(3,i,j,k)
#endif /*CFD*/
          end select
#ifdef BORIS_CORRECTION
          b2           = bxt**2+byt**2+bzt**2 
          bv           = (bxt*vx+byt*vy+bzt*vz)
          ux           = rr*vx + afact*(b2*vx - bv*bxt)
          uy           = rr*vy + afact*(b2*vy - bv*byt)
          uz           = rr*vz + afact*(b2*vz - bv*bzt)
#else
          ux           = rr*vx 
          uy           = rr*vy 
          uz           = rr*vz 
#endif /*BORIS_CORRECTION*/
           vl(2,i,j,k)  = ux 
           vl(3,i,j,k)  = uy 
           vl(4,i,j,k)  = uz 
           vl(6,i,j,k)  = bx 
           vl(7,i,j,k)  = by 
           vl(8,i,j,k)  = bz 
           p            = vl(5,i,j,k)
           v2           = vx**2+vy**2+vz**2 
           ek           = 0.5*rr*v2
           et           = ek+p*gammam1i+bsq*pie8i
           vl(5,i,j,k)  = et
        enddo
       enddo
      enddo
      do k = nz1,nz2
       do j = ny1,ny2
        do i = nx1,nx2
           rr          = max(vr(1,i,j,k),rmin)
           ri          = 1.0/rr
           vx          = vr(2,i,j,k)
           vy          = vr(3,i,j,k)
           vz          = vr(4,i,j,k)
#ifdef CFD
           bx          = 0.0
           by          = 0.0
           bz          = 0.0
#else
           bx          = vr(6,i,j,k)
           by          = vr(7,i,j,k)
           bz          = vr(8,i,j,k)
#endif /*CFD*/
          bsq         = bx*bx+by*by+bz*bz
          select case (idir)
          case (1)
#ifdef CFD
           bxt         = 0.0
           byt         = 0.0
           bzt         = 0.0
#else
           bxt         = bx+bxf_sf(1,i,j,k)
           byt         = by+bxf_sf(2,i,j,k)
           bzt         = bz+bxf_sf(3,i,j,k)
#endif /*CFD*/
!---------------------------------------------------
          case (2)
#ifdef CFD
           bxt         = 0.0
           byt         = 0.0
           bzt         = 0.0
#else
           bxt         = bx+byf_sf(1,i,j,k)
           byt         = by+byf_sf(2,i,j,k)
           bzt         = bz+byf_sf(3,i,j,k)
#endif /*CFD*/
!---------------------------------------------------
          case (3)
#ifdef CFD
           bxt         = 0.0
           byt         = 0.0
           bzt         = 0.0
#else
           bxt         = bx+bzf_sf(1,i,j,k)
           byt         = by+bzf_sf(2,i,j,k)
           bzt         = bz+bzf_sf(3,i,j,k)
#endif /*CFD*/
          end select
#ifdef BORIS_CORRECTION
          b2           = bxt**2+byt**2+bzt**2 
          bv           = (bxt*vx+byt*vy+bzt*vz)
          ux           = rr*vx + afact*(b2*vx - bv*bxt)
          uy           = rr*vy + afact*(b2*vy - bv*byt)
          uz           = rr*vz + afact*(b2*vz - bv*bzt)
#else
          ux           = rr*vx 
          uy           = rr*vy 
          uz           = rr*vz 
#endif /*BORIS_CORRECTION*/
           vr(2,i,j,k)  = ux 
           vr(3,i,j,k)  = uy 
           vr(4,i,j,k)  = uz 
           vr(6,i,j,k)  = bx 
           vr(7,i,j,k)  = by 
           vr(8,i,j,k)  = bz 
           p            = vr(5,i,j,k)
           v2           = vx**2+vy**2+vz**2 
           ek           = 0.5*rr*v2
           et           = ek+p*gammam1i+bsq*pie8i
           vr(5,i,j,k)  = et
        enddo
       enddo
      enddo
      return
      end subroutine eos_vf_to_uf
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
      subroutine eos_uf_to_vf(vl,vr,nx1,nx2,ny1,ny2,nz1,nz2,idir)
      use physcons
      use paramesh_dimensions
      use strongfields
      use athena_specific_parameters
      integer :: nx1,nx2,ny1,ny2,nz1,nz2,idir
      real :: bsqw,bsq
      real,dimension(ne,iu_bnd1p,ju_bnd1p,ku_bnd1p),
     &                                  intent(inout)::vl,vr
       do k = nz1,nz2
        do j = ny1,ny2
         do i = nx1,nx2
          rr           = max(vl(1,i,j,k),rmin)
          ri           = 1.0/rr
          px           = vl(2,i,j,k)
          py           = vl(3,i,j,k)
          pz           = vl(4,i,j,k)
          bx           = vl(6,i,j,k)
          by           = vl(7,i,j,k)
          bz           = vl(8,i,j,k)
          select case (idir)
          case (1)
           bxt         = bx+bxf_sf(1,i,j,k)
           byt         = by+bxf_sf(2,i,j,k)
           bzt         = bz+bxf_sf(3,i,j,k)
          case (2)
           bxt         = bx+byf_sf(1,i,j,k)
           byt         = by+byf_sf(2,i,j,k)
           bzt         = bz+byf_sf(3,i,j,k)
          case (3)
           bxt         = bx+bzf_sf(1,i,j,k)
           byt         = by+bzf_sf(2,i,j,k)
           bzt         = bz+bzf_sf(3,i,j,k)
          end select
          bsq          = bx*bx+by*by+bz*bz
#ifdef BORIS_CORRECTION
          b2           = bxt*bxt+byt*byt+bzt*bzt
          rbi          = 1.0/(rr + afact*b2)
          bv           = (px*bxt+py*byt+pz*bzt)*rbi

          vx           = (px + afact*bv*bxt)*rbi
          vy           = (py + afact*bv*byt)*rbi
          vz           = (pz + afact*bv*bzt)*rbi
#else
          vx           = px*ri
          vy           = py*ri
          vz           = pz*ri
#endif /*BORIS_CORRECTION*/
          et           = vl(5,i,j,k)
          vsq          = vx*vx+vy*vy+vz*vz
          ek           = 0.5*rr*vsq
! Start GM Added
!          write (*,*)  rr, ek
! End GM Added
          bsq          = pie8i*bsq
          pt           = et - ek - bsq
          p            = max(gammam1*pt,pmin)
          vl(2,i,j,k)  = vx
          vl(3,i,j,k)  = vy
          vl(4,i,j,k)  = vz
          vl(5,i,j,k)  = p
          vl(6,i,j,k)  = bx
          vl(7,i,j,k)  = by
          vl(8,i,j,k)  = bz
        enddo
       enddo
      enddo
       do k = nz1,nz2
        do j = ny1,ny2
         do i = nx1,nx2
          rr           = max(vr(1,i,j,k),rmin)
          ri           = 1.0/rr
          px           = vr(2,i,j,k)
          py           = vr(3,i,j,k)
          pz           = vr(4,i,j,k)
          bx           = vr(6,i,j,k)
          by           = vr(7,i,j,k)
          bz           = vr(8,i,j,k)
          select case (idir)
          case (1)
           bxt         = bx+bxf_sf(1,i,j,k)
           byt         = by+bxf_sf(2,i,j,k)
           bzt         = bz+bxf_sf(3,i,j,k)
          case (2)
           bxt         = bx+byf_sf(1,i,j,k)
           byt         = by+byf_sf(2,i,j,k)
           bzt         = bz+byf_sf(3,i,j,k)
          case (3)
           bxt         = bx+bzf_sf(1,i,j,k)
           byt         = by+bzf_sf(2,i,j,k)
           bzt         = bz+bzf_sf(3,i,j,k)
          end select
          bsq          = bx*bx+by*by+bz*bz
#ifdef BORIS_CORRECTION
          b2           = bxt*bxt+byt*byt+bzt*bzt
          rbi          = 1.0/(rr + afact*b2)
          bv           = (px*bxt+py*byt+pz*bzt)*rbi

          vx           = (px + afact*bv*bxt)*rbi
          vy           = (py + afact*bv*byt)*rbi
          vz           = (pz + afact*bv*bzt)*rbi
#else
          vx           = px*ri
          vy           = py*ri
          vz           = pz*ri
#endif /*BORIS_CORRECTION*/
          et           = vr(5,i,j,k)
          vsq          = vx*vx+vy*vy+vz*vz
          ek           = 0.5*rr*vsq
          bsq          = pie8i*bsq
          pt           = et - ek - bsq
          p            = max(gammam1*pt,pmin)
          vr(2,i,j,k)  = vx
          vr(3,i,j,k)  = vy
          vr(4,i,j,k)  = vz
          vr(5,i,j,k)  = p
          vr(6,i,j,k)  = bx
          vr(7,i,j,k)  = by
          vr(8,i,j,k)  = bz
       enddo
       enddo
      enddo
      return
      end subroutine eos_uf_to_vf
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
