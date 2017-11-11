      module interfaces

      interface
      subroutine riemann_solve(ud,divv,vxfl,vxfr,fx,vf,idir,is,ir,lflag,lb)      
      use paramesh_dimensions
      use physcons
      use physicaldata
      use tree
      use strongfields
      use interior_gbc
      use grid
      use athena_specific_parameters
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

      end subroutine riemann_solve
      end interface

      interface
      subroutine reconstruction(v1,chi,vxfl,vxfr,dt,idir,lflag)      
      use paramesh_dimensions
      use physcons
      use physicaldata
      use tree
      use strongfields
      use interior_gbc
      use grid
      use athena_specific_parameters
#define cell_v real,dimension(ne,iu_bnd1 ,ju_bnd1 ,ku_bnd1 )
#define cell_s real,dimension(   iu_bnd1 ,ju_bnd1 ,ku_bnd1 )
#define face_v real,dimension(ne,iu_bnd1p,ju_bnd1p,ku_bnd1p)
#define face_f real,dimension(nf,iu_bnd1p,ju_bnd1p,ku_bnd1p)

      implicit none
      logical    ,intent(in   )::lflag
      integer    ,intent(in   )::idir
      real       ,intent(in   )::dt
      cell_s     ,intent(in   )::chi
      cell_v     ,intent(in   )::v1
      face_v     ,intent(inout)::vxfl,vxfr
      end subroutine reconstruction
      end interface

      interface
      subroutine correct_rl_states_2d(vxfl   ,vxfr   ,vyfl   ,vyfr   ,
     &                                fx     ,fy     ,
     &                                index_x,index_y,
     &                                dtdx   ,dtdy   )
      use paramesh_dimensions
      use physcons
      use physicaldata
      use tree
      use strongfields
      use interior_gbc
      use grid
      use athena_specific_parameters
      implicit none
#define cell_v real,dimension(ne,iu_bnd1 ,ju_bnd1 ,ku_bnd1 )
#define cell_s real,dimension(   iu_bnd1 ,ju_bnd1 ,ku_bnd1 )
#define face_v real,dimension(ne,iu_bnd1p,ju_bnd1p,ku_bnd1p)
#define face_f real,dimension(nf,iu_bnd1p,ju_bnd1p,ku_bnd1p)
      integer,dimension(ne+3),intent(in   )::index_x,index_y 
      face_v                 ,intent(inout)::vxfl,vxfr,vyfl,vyfr
      face_f                 ,intent(in   )::fx,fy
      real                   ,intent(in   )::dtdx,dtdy
      end subroutine correct_rl_states_2d
      end interface

      interface
      subroutine correct_rl_states_3d_xa(vxfl   ,vxfr   ,
     &                                   vxyfl  ,vxyfr  ,vxzfl  ,vxzfr  ,
     &                                   fy     ,fz     ,
     &                                   index_x,
     &                                   dtdy   ,dtdz   )
      use paramesh_dimensions
      use physcons
      use physicaldata
      use tree
      use strongfields
      use interior_gbc
      use grid
      use athena_specific_parameters
      implicit none
#define cell_v real,dimension(ne,iu_bnd1 ,ju_bnd1 ,ku_bnd1 )
#define cell_s real,dimension(   iu_bnd1 ,ju_bnd1 ,ku_bnd1 )
#define face_v real,dimension(ne,iu_bnd1p,ju_bnd1p,ku_bnd1p)
#define face_f real,dimension(nf,iu_bnd1p,ju_bnd1p,ku_bnd1p)
      integer,dimension(ne+3),intent(in   )::index_x 
      face_v                 ,intent(in   )::vxfl,vxfr
      face_v                 ,intent(  out)::vxyfl,vxyfr
      face_v                 ,intent(  out)::vxzfl,vxzfr
      face_f                 ,intent(in   )::fy,fz
      real                   ,intent(in   )::dtdy,dtdz
      end subroutine correct_rl_states_3d_xa
      end interface

      interface
      subroutine correct_rl_states_3d_ya(vyfl   ,vyfr ,
     &                                   vyxfl  ,vyxfr,
     &                                   vyzfl  ,vyzfr,
     &                                   fx     ,fz   ,
     &                                   index_y,
     &                                   dtdx   ,dtdz  )
      use paramesh_dimensions
      use physcons
      use physicaldata
      use tree
      use strongfields
      use interior_gbc
      use grid
      use athena_specific_parameters
      implicit none
#define cell_v real,dimension(ne,iu_bnd1 ,ju_bnd1 ,ku_bnd1 )
#define cell_s real,dimension(   iu_bnd1 ,ju_bnd1 ,ku_bnd1 )
#define face_v real,dimension(ne,iu_bnd1p,ju_bnd1p,ku_bnd1p)
#define face_f real,dimension(nf,iu_bnd1p,ju_bnd1p,ku_bnd1p)
      integer,dimension(ne+3),intent(in   )::index_y 
      face_v                 ,intent(in   )::vyfl ,vyfr
      face_v                 ,intent(  out)::vyxfl,vyxfr
      face_v                 ,intent(  out)::vyzfl,vyzfr
      face_f                 ,intent(in   )::fx,fz
      real                   ,intent(in   )::dtdx,dtdz
      end subroutine correct_rl_states_3d_ya
      end interface

      interface
      subroutine correct_rl_states_3d_za(vzfl  ,vzfr   ,
     &                                  vzxfl  ,vzxfr  ,vzyfl  ,vzyfr  ,
     &                                  fx     ,fy     ,
     &                                  index_z,
     &                                  dtdx   ,dtdy   )
      use paramesh_dimensions
      use physcons
      use physicaldata
      use tree
      use strongfields
      use interior_gbc
      use grid
      use athena_specific_parameters
      implicit none
#define cell_v real,dimension(ne,iu_bnd1 ,ju_bnd1 ,ku_bnd1 )
#define cell_s real,dimension(   iu_bnd1 ,ju_bnd1 ,ku_bnd1 )
#define face_v real,dimension(ne,iu_bnd1p,ju_bnd1p,ku_bnd1p)
#define face_f real,dimension(nf,iu_bnd1p,ju_bnd1p,ku_bnd1p)
      integer,dimension(ne+3),intent(in   )::index_z 
      face_v                 ,intent(in   )::vzfl ,vzfr
      face_v                 ,intent(  out)::vzxfl,vzxfr
      face_v                 ,intent(  out)::vzyfl,vzyfr
      face_f                 ,intent(in   )::fx,fy
      real                   ,intent(in   )::dtdx,dtdy
      end subroutine correct_rl_states_3d_za
      end interface

      interface
      subroutine correct_rl_states_3d_xb(vxfl   ,vxfr   ,
     &                                   fyz    ,fzy    ,
     &                                   index_x,
     &                                   dtdy   ,dtdz   )
      use paramesh_dimensions
      use physcons
      use physicaldata
      use tree
      use strongfields
      use interior_gbc
      use grid
      use athena_specific_parameters
      implicit none
#define cell_v real,dimension(ne,iu_bnd1 ,ju_bnd1 ,ku_bnd1 )
#define cell_s real,dimension(   iu_bnd1 ,ju_bnd1 ,ku_bnd1 )
#define face_v real,dimension(ne,iu_bnd1p,ju_bnd1p,ku_bnd1p)
#define face_f real,dimension(nf,iu_bnd1p,ju_bnd1p,ku_bnd1p)
      integer,dimension(ne+3),intent(in   )::index_x 
      face_v                 ,intent(inout)::vxfl,vxfr
      face_f                 ,intent(in   )::fzy,fyz
      real                   ,intent(in   )::dtdy,dtdz
      end subroutine correct_rl_states_3d_xb
      end interface

      interface
      subroutine correct_rl_states_3d_yb(vyfl   ,vyfr   ,
     &                                   fxz    ,fzx    ,
     &                                   index_y,
     &                                   dtdx   ,dtdz   )
      use paramesh_dimensions
      use physcons
      use physicaldata
      use tree
      use strongfields
      use interior_gbc
      use grid
      use athena_specific_parameters
      implicit none
#define cell_v real,dimension(ne,iu_bnd1 ,ju_bnd1 ,ku_bnd1 )
#define cell_s real,dimension(   iu_bnd1 ,ju_bnd1 ,ku_bnd1 )
#define face_v real,dimension(ne,iu_bnd1p,ju_bnd1p,ku_bnd1p)
#define face_f real,dimension(nf,iu_bnd1p,ju_bnd1p,ku_bnd1p)
      integer,dimension(ne+3),intent(in   )::index_y 
      face_v                 ,intent(inout)::vyfl,vyfr
      face_f                 ,intent(in   )::fxz,fzx
      real                   ,intent(in   )::dtdx,dtdz
      end subroutine correct_rl_states_3d_yb
      end interface


      interface
      subroutine correct_rl_states_3d_zb(vzfl   ,vzfr   ,
     &                                   fxy    ,fyx    ,
     &                                   index_z,
     &                                   dtdx   ,dtdy   )
      use paramesh_dimensions
      use physcons
      use physicaldata
      use tree
      use strongfields
      use interior_gbc
      use grid
      use athena_specific_parameters
      implicit none
#define cell_v real,dimension(ne,iu_bnd1 ,ju_bnd1 ,ku_bnd1 )
#define cell_s real,dimension(   iu_bnd1 ,ju_bnd1 ,ku_bnd1 )
#define face_v real,dimension(ne,iu_bnd1p,ju_bnd1p,ku_bnd1p)
#define face_f real,dimension(nf,iu_bnd1p,ju_bnd1p,ku_bnd1p)
      integer,dimension(ne+3),intent(in   )::index_z 
      face_v                 ,intent(inout)::vzfl,vzfr
      face_f                 ,intent(in   )::fxy,fyx
      real                   ,intent(in   )::dtdx,dtdy
      end subroutine correct_rl_states_3d_zb
      end interface


      interface
      subroutine block_ctu(v1    ,u1    ,
     &                     fx    ,fy    ,fz    ,
     &                     ex    ,ey    ,ez    ,
     &                     du    ,dl    ,
     &                     dt    ,time  ,lb    ,
     &                     lflag ,istep ,mype  )
!=======================================================================
      use paramesh_dimensions
      use physcons
      use physicaldata
      use tree
      use strongfields
      use interior_gbc
      use grid
      use athena_specific_parameters
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

      end subroutine block_ctu
      end interface

      interface
      subroutine roe_mhd(n1,n2,vml,vmr,uc,bn,bf,divv,flux,
     &                      vfout,idir)
      use paramesh_dimensions
      use physcons
      use grid
      use athena_specific_parameters
      integer                    ,intent(in   )::n1,n2,idir
      real, dimension(  8,nmax+1),intent(in   )::vml,vmr
      real, dimension(  8,nmax  ),intent(in   )::uc
      real, dimension(    nmax  ),intent(in   )::bn,divv 
      real, dimension( 11,nmax  ),intent(  out)::flux
      real, dimension(  3,nmax+1),intent(in   )::bf
      real, dimension(  8,nmax  ),intent(  out)::vfout
      end subroutine roe_mhd
      end interface

      interface
       subroutine hlle_mhd(n1,n2,vml,vmr,uc,bn,bf,divv,flux,
     &                     vfout,idir,lb)
       use physcons
       use paramesh_dimensions
       use athena_specific_parameters
       implicit none
       integer                    ,intent(in   )::n1,n2,idir,lb
       real, dimension(  8,nmax+1),intent(in   )::vml,vmr
       real, dimension(  8,nmax  ),intent(in   )::uc
       real, dimension(    nmax  ),intent(in   )::bn,divv 
       real, dimension( 11,nmax  ),intent(  out)::flux
       real, dimension(  3,nmax+1),intent(in   )::bf
       real, dimension(  8,nmax  ),intent(  out)::vfout
      end subroutine hlle_mhd
      end interface

      interface
      subroutine muscl_1d(vc,flat,bx,by,bz,vl,vr,n1,n2,dt,idir,lflag)
      use paramesh_dimensions
      use physcons
      use grid
      use strongfields
      use interior_gbc
      use athena_specific_parameters
      implicit none
      logical                         ,intent(in   )::lflag
      integer                         ,intent(in   )::n1,n2,idir
      real                            ,intent(in   )::dt
      real,dimension(ne  ,     nmax  ),intent(in   )::vc
      real,dimension(          nmax  ),intent(in   )::flat
      real,dimension(          nmax  ),intent(in   )::bx,by,bz
      real,dimension(ne  ,     nmax+1),intent(  out)::vl,vr
      end subroutine muscl_1d
      end interface

      interface
      subroutine stonefix(vx,vy,vz,vxfl,vxfr,vyfl,vyfr,vzfl,vzfr,dt)      
      use paramesh_dimensions
      use physcons
      use physicaldata
      use tree
      use strongfields
      use interior_gbc
      use grid
      use athena_specific_parameters
#define face_x  real,dimension(ne,iu_bnd1p,ju_bnd1p,ku_bnd1p)
#define cell_v  real,dimension(ne,iu_bnd1 ,ju_bnd1 ,ku_bnd1 )
      implicit none
      real  ,intent(in   )::dt
      cell_v,intent(in   )::vx
      cell_v,intent(in   )::vy
      cell_v,intent(in   )::vz
      face_x,intent(inout)::vxfl ,vxfr
      face_x,intent(inout)::vyfl ,vyfr
      face_x,intent(inout)::vzfl ,vzfr
      end subroutine stonefix
      end interface

      interface
      subroutine eos_u_to_v_pmn(u,v,cf,lb,mype)
      use physcons
      use paramesh_dimensions
      use physicaldata
      use tree
      use strongfields
      use athena_specific_parameters

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
      end subroutine eos_u_to_v_pmn
      end interface
 
      interface
      subroutine brackbill(v2,src,dt)
      use paramesh_dimensions
      use physcons
      use physicaldata
      use tree
      use strongfields
      use grid
      use athena_specific_parameters
      implicit none
      real                                 ,intent(in)    ::dt
      real,dimension(ne,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,
     .                  kl_bnd1:ku_bnd1    ),intent(in)   ::v2
      real,dimension(ne,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,
     .                  kl_bnd1:ku_bnd1    ),intent(inout)::src

      end subroutine brackbill
      end interface

      interface
      subroutine interior_bc_fl(fx,fy,fz,vfx,vfy,vfz,v2,lb)
      use physcons
      use paramesh_dimensions
      use physicaldata
      use tree
      use strongfields
      use athena_specific_parameters
      implicit none
c#define face_x  real,dimension(ne,iu_bnd1p,ju_bnd1p,ku_bnd1p)
c#define face_f  real,dimension(nf,iu_bnd1p,ju_bnd1p,ku_bnd1p)
c#define cell_v  real,dimension(ne,iu_bnd1 ,ju_bnd1 ,ku_bnd1 )
      integer,intent(in  ) ::  lb
      face_f ,intent(inout) :: fx
      face_f ,intent(inout) :: fy
      face_f ,intent(inout) :: fz
      face_x ,intent(inout) ::vfx
      face_x ,intent(inout) ::vfy
      face_x ,intent(inout) ::vfz
      cell_v ,intent(in   ) :: v2
      end subroutine interior_bc_fl
      end interface

      interface
      subroutine interior_bc_states_nplus_1half(vlx,vrx,idir,lb)
      use paramesh_dimensions
      use physicaldata
      use tree
      use athena_specific_parameters
      implicit none
      integer,intent(in   )::idir,lb
      face_x ,intent(inout)::vlx
      face_x ,intent(inout)::vrx
      end subroutine interior_bc_states_nplus_1half
      end interface

c      interface
c      subroutine interior_bc_states_n(v,idir,lb)
c      use paramesh_dimensions
c      use physicaldata
c      use tree
c      use athena_specific_parameters
c      implicit none
c      integer,intent(in   )::idir,lb
c      cell_v ,intent(inout)::v
c      end subroutine interior_bc_states_n
c      end interface


      end module interfaces
