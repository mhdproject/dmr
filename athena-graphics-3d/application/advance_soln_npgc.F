#include "paramesh_preprocessor.fh"
#include "athena_preprocessor.fh"
!#define DEBUGX
!#define DEBUGG
!#define DEBUGGXX
#undef OVERWRITE
      subroutine advance_soln(istep,mype,nprocs,time,dt)
      use physcons
      use paramesh_dimensions
      use physicaldata
      use interior_gbc
      use strongfields
      use grid
      use athena_specific_parameters

! include file defining the tree
      use tree
      use workspace
      use paramesh_interfaces, only : amr_1blk_copy_soln,
     &                                amr_1blk_guardcell,
     &                                amr_flux_conserve,
     &                                amr_edge_average,
     &                                amr_restrict      ,
     &                                amr_1blk_to_perm

c      use athena_paramesh_interfaces, only : amr_timestep
     

#ifdef MPI_USED
      use paramesh_mpi_interfaces, only : 
     &                                mpi_amr_comm_setup,
     &                                mpi_amr_tree_setup,
!     &                                mpi_amr_local_surr_blks

     &                                mpi_morton_bnd
!    &                                mpi_amr_1blk_guardcell,
!    &                                mpi_morton_bnd_prolong,
#endif /* MPI_USED */
      use interfaces, only :  block_ctu          

c      implicit none
#ifdef MPI_USED
      include "mpif.h"
#endif /* MPI_USED */

      integer :: tag_offset,max_blks_sent,iproc
      integer :: istep,mype,nprocs
      real    :: time,dt

!
!
!--------------------------------------------------------------

#include "amr_shmem.fh"
#include "shmem_reduce.fh"

!--------------------------------------------------------------

! Store the assignment between block boundary fluxes and the physical
! variables.
      common/flux_assign/ iflux_target(nfluxvar)
      integer iflux_target
!-----------------------------------------------------------------------
#ifdef IONOSPHERE
! include file defining the ionospheric parameters
#include "ionosphere.fh"
#endif /*IONOSPHERE*/ 

!-----------------------------------------------------------------------
      logical lcc,lfc,ldiag,l_srl_only,low_order_space,lec,lnc
      logical lfullblock

      real    dx,dy,dz,dxi,dyi,dzi
      real, dimension(nf,il_bnd1:iu_bnd1p,jl_bnd1:ju_bnd1p,
     .                   kl_bnd1:ku_bnd1p)::fx,fy,fz

      integer i,j,k,idest,nvarp,iv,iflx,ivar,ndel
!
      real, dimension(ne,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,
     .                   kl_bnd1:ku_bnd1    ) :: v1,u1
      real, dimension(ne,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,
     .                   kl_bnd1:ku_bnd1    ) :: delu
      real, dimension( 3,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,
     .                   kl_bnd1:ku_bnd1    ) :: dl
!
      real, dimension(il_bnd1:iu_bnd1p,jl_bnd1:ju_bnd1p,
     .                   kl_bnd1:ku_bnd1p    ) :: phi3d1
!
!-----1d scratch arrays
!
     &                                             lwall ,lcorner
!-----------------------------------------------------------------------

!
!-----Electric Field at edges
!
      real, dimension(il_bnd1:iu_bnd1p,jl_bnd1:ju_bnd1p,
     &                kl_bnd1:ku_bnd1p) ::ex,ey,ez
!-----------------------------------------------------------------------
!
!-----Array needed for projection of Lorentz force
!
      real, dimension(3,il_bnd:iu_bnd,jl_bnd:ju_bnd,
     &                  kl_bnd:ku_bnd,maxblocks) :: dmw

!
!-----------------------------------------------------------
!
      real, dimension(ne          ):: varmin,varmax
      real, dimension(ne          ):: uarmin,uarmax
      real, dimension(ne,maxblocks):: vmin  ,vmax  ,umin,umax
      save varmax,varmin
      integer, dimension(ne      ) :: ivmax,ivmin,iumax,iumin

      integer :: ic,lb,ii,jj,kk,ng0,l,nsub,lf,icoord,iopt,nlayers
      integer :: iempty,f,lb1,lb2,lb3,lb4
      integer :: ix,iy,iz,jx,jy,jz
      integer :: nxl,nxh,nxhm1,nyl,nyh,nyhm1,nzl,nzh,nzhm1
      real    :: dtold
      common/old_timestep/dtold
      integer :: nxc,nyc,nzc
      integer :: mxc,myc,mzc
      common/ncentral/nxc,nyc,nzc
      common/mcentral/mxc,myc,mzc
      integer :: ncells
      logical :: lflag
      logical :: lguard,lprolong,lflux,leedge,lrestrict,lfulltree
      integer :: errcode

      real :: bxmin,bxmax,bymin,bymax
      real :: bt1,bt2,bt3,bsqi,bdotd,bd
      real :: tf,dtmin,dtmax
      real,dimension(il_bnd:iu_bnd,jl_bnd:ju_bnd,
     .               kl_bnd:ku_bnd) :: b_sq_o,b_sq_n
      real :: bx,by,bz,bx1,by1,bz1,be
      real :: vx,vy,vz
      real :: px,py,pz
      real :: divb,ub
      real :: rr,vsq1,vsq2,bv,rbi,bsq

      integer :: mylb


!     open(5,FILE='pointers.dat')
!
!-----pointer declarations
!
!

      ng0  = nguard*npgs    
      nxl  = il_bnd1+nguard
      nxh  = iu_bnd1-(nguard-1)
      nxhm1= nxh-1
      nyl  = jl_bnd1+nguard*k2d
      nyh  = ju_bnd1-(nguard-1)*k2d
      nyhm1= nyh-1
      nzl  = kl_bnd1+nguard*k3d
      nzh  = ku_bnd1-(nguard-1)*k3d
      nzhm1= nzh-1
#if N_DIM < 3
      nzl  =1
      nzhm1=1
      nzh  =1
#endif
#ifdef LORENTZ_FORCE
      dmw(:,:,:,:,:)  = 0.0
#endif /*LORENTZ_FORCE*/

!--------------------------------------------------------------
! loop over leaf grid blocks
 
#ifndef ADVANCE_ALL_LEVELS
! A valid solution will be required on the parents of leaf blocks
! when refinement testing is done. See the comment before the call
! to amr_test_refinement.
        iempty = 0
        iopt   = 1
        call amr_restrict(mype,iopt,iempty)
#endif


c       if(mype.eq.0)write(*,*)'mk40: timestep: istep, dt=',istep, dt
        call amr_timestep(dt,dtmin,dtmax,mype)
c       if(mype.eq.0)write(*,*)'mk40: timestep: istep, dt=',istep, dt

#ifdef NO_PERMANENT_GUARDCELLS
! Store a copy of the current solution in gt_unk

        call amr_1blk_copy_soln(-1)

#endif

#ifdef MPI_USED
        call mpi_barrier (MPI_COMM_WORLD, errcode)
#else
        call shmem_barrier_all()
#endif /*MPI_USED*/

      iopt       = 1
      lcc        = .true.
      lfc        = .true.
#ifdef MPI_USED
      lguard    = .true.
      lprolong  = .false.
      lflux     = .false.
      leedge    = .false.
      lrestrict = .false.
      lec        =.false.
      lnc        =.false.

      lfulltree=.false.

      tag_offset= 100
      call mpi_amr_comm_setup(mype,nprocs,
     .                        lguard,lprolong,lflux,leedge,lrestrict,
     .                        lfulltree,
     .                        iopt,lcc,lfc,lec,lnc,tag_offset)

#endif /*MPI_USED*/

!-----------------------------------------------------------------------
! loop over leaf grid blocks


c      do lb=1,lnblocks
c         if(istep.eq.39.and.nodetype(lb).eq.1) then
c            write(*,*) mype, lb,time, 'adddv coord', coord(1:2,lb)
c         endif
c      enddo


      if(lnblocks.gt.0) then
       do lb=1,lnblocks
#ifndef ADVANCE_ALL_LEVELS
        if(nodetype(lb).eq.1) then
#endif


#ifdef NO_PERMANENT_GUARDCELLS
! Copy data from current block into working block and fill its guardcells

        iopt      = 1
        nlayers    = nguard
        lcc        = .true.
        lfc        = .false.
        l_srl_only = .false.
        icoord     = 0
        ldiag      = .true.
        lec        = .false.
        lnc        = .false.

!
!-----call amr_1blk_guardcell: selects a block including guard cells
!

      call grid_variables(lb,dt)
      call strongfields_wb(mype,lb)
      call amr_1blk_guardcell(mype,iopt,nlayers,lb,mype,lcc,lfc,
     .                        lec,lnc,
     .                        l_srl_only,icoord,ldiag)


!mk0
!     if(lb.eq.1) call ch_data(mype,lb)

!
!-----unk1 now has the unknown and the necessary guard cell information
!     for block lb
!
#endif
!
        u1(1:ne,:,:,:)=unk1(1:ne  ,:,:,:,1) ! time level n


c        write(*,*) lb,'IN ADV: SHOWING UNK'
c        do k=nzl,nzhm1
c           do j=nyl,nyhm1
c              do i=nxl,nxhm1
c                 do l=1,ne


!-----------------------------------------------------------------------
!
      call grid_variables(lb,dt)
      call strongfields_wb(mype,lb)
!
!-----convert conserved variables to primitives
!
        call eos_u_to_v(u1,v1,lb)
#ifdef INTERIOR_BOUNDARY
          lflag             = earth_blocks(lb)
          izones            = idzones(lb) 
          ifaces            = idfaces(lb)
          iedges            = idefields(lb)
          izone(:,:izones)  = p_zone(:,:izones,lb) 
          iface(:,:ifaces)  = p_face(:,:ifaces,lb)
          iedge(:,:iedges)  = p_efields(:,:iedges,lb)
#else
          lflag             = .false.
#endif /*INTERIOR_BOUNDARY*/
c       if(mype.eq.0)write(*,*)'entering block_ctu,lb=',lb
       call block_ctu(v1    ,u1    ,
     &                fx    ,fy    ,fz    ,
     &                ex    ,ey    ,ez    ,
     &                delu  ,dl    ,
     &                dt    ,time  ,lb    ,
     &                lflag, istep ,mype )



c       if(mype.eq.2.and.lb.eq.1) then
c          write(*,*) mype, time,'adv1',unk1(1,11,5,1,1),delu(1,11,5,1)
c       endif
c       if(mype.eq.0.and.lb.eq.5) then
c          write(*,*) mype, time,'adv1',unk1(1,6,5,1,1),delu(1,6,5,1)
c       endif



#ifdef PULSAR
c      call pulsar_injection(u1,dt,istep,lb,mype)
#endif

!-----------------------------------------------------------------------
c       if(mype.eq.0)write(*,*)'leaving block_ctu,lb=',lb

!
!-----apply boundary conditions to delu for internal boundaries
!
      do k=nzl,nzhm1
       do j=nyl,nyhm1
        do i=nxl,nxhm1
         do l=1,ne
          unk1(l,i,j,k,1) = unk1(l,i,j,k,1)+delu(l,i,j,k)
c          if(unk1(6,i,j,k,1).ne.0.0.or.unk1(7,i,j,k,1).ne.0.0)then
c           write(*,*)'i,j,k=',i,j,k
c           write(*,*)'unk1(6,i,j,k,1)=',unk1(6,i,j,k,1)
c           write(*,*)'unk1(7,i,j,k,1)=',unk1(7,i,j,k,1)
c          endif
         enddo
        enddo   
       enddo   
      enddo

c
c
c



!
#ifdef LORENTZ_FORCE
      do k=nzl,nzhm1
       do j=nyl,nyhm1
        do i=nxl,nxhm1
            ii = i-nguard 
            jj = j-nguard*k2d 
            kk = k-nguard*k3d
            dmw(1:3,ii,jj,kk,lb)  = dl(1:3,i,j,k)
        enddo   
       enddo   
      enddo
#else
      dmw(1:3,:,:,:,lb)  = 0.0
#endif /*LORENTZ_FORCE*/
#ifdef BFILTER
!
#ifdef INTERIOR_BOUNDARY
!
!-----set boundary conditions on ionospheric electric fields
!
#ifdef IONOSPHERE



! get phi3d1 for this block.
#ifdef NO_PERMANENT_GUARDCELLS
        phi3d1(1+nguard:nxb+nguard+1,1+nguard:nyb+nguard+k2d,
     &         1+nguard:nzb+nguard+k3d) = 
     &                     phi3d(1:nxb+1,1:nyb+k2d,1:nzb+k3d,lb)
#else
      phi3d1(:,:,:,) = phi3d(:,:,:,lb)
#endif /* NO_PERMANENT_GUARDCELLS */

      call bndry_electric_field(ex,ey,ez,learth, phi3d1(1,1,1),lb)
#endif /*IONOSPHERE*/
#endif /*INTERIOR_BOUNDARY*/

!
!  F = Bx dy dz + By dz dx + Bz dx dy + Ex dx dt + Ey dy dt + Ez dz dt
!  d F = (div(B)) dx dy dz + (Bx,t + Ez,y - Ey,z) dy dz dt
!                          + (By,t + Ex,z - Ez,x) dz dx dt
!                          + (Bz,t + Ey,x - Ex,y) dx dy dt
!  0 = Bx,t + Ez,y - Ey,z = Bx,t + (Vy Bx - Vx By),y + (Vz Bx - Vx Bz),z


            do iz= 1+nguard*k3d, nzb+nguard*k3d
               do iy= 1+nguard*k2d, nyb+nguard*k2d
                  do ix= 1+nguard  , nxb+(1+nguard) 
                     facevarx1(1,ix,iy,iz,1)=facevarx1(1,ix,iy,iz,1)
     &                    -(ez(ix,iy+k2d,iz)
     &                    - ez(ix,iy    ,iz))*dtdyi
#if N_DIM == 3
     &                    +(ey(ix,iy,iz+k3d)
     &                    - ey(ix,iy,iz    ))*dtdzi
#endif /*NDIM3*/
                  end do
               end do
            end do


!  0 = By,t + Ex,z - Ez,x = By,t + (Vz By - Vy Bz),z + (Vx By - Vy Bx),x
            do iz= 1+nguard*k3d, nzb+nguard*k3d
               do iy= 1+nguard*k2d, nyb+(1+nguard)*k2d
                  do ix= 1+nguard  , nxb+nguard

!------------------------------------------------------------------
                     facevary1(1,ix,iy,iz,1)=facevary1(1,ix,iy,iz,1)
     &                    +(ez(ix+ 1,iy,iz )
     &                    - ez(ix   ,iy,iz ))*dtdxi
#if N_DIM == 3
                     facevary1(1,ix,iy,iz,1)=facevary1(1,ix,iy,iz,1)
     &                    -(ex(ix,iy,iz+k3d)
     &                    - ex(ix,iy,iz    ))*dtdzi
#endif /*NDIM3*/

                  end do
               end do
            end do


!  0 = Bz,t + Ey,x - Ex,y = Bz,t + (Vx Bz - Vz Bx),x + (Vy Bz - Vz By),y
          do iz= 1+nguard*k3d, nzb+(1+nguard)*k3d
             do iy= 1+nguard*k2d, nyb+nguard*k2d
                do ix= 1+nguard  , nxb+nguard
                   facevarz1(1,ix,iy,iz,1)=facevarz1(1,ix,iy,iz,1)
     &                 -(ey(ix+1,iy    ,iz)
     &                 - ey(ix  ,iy    ,iz))*dtdxi
     &                 +(ex(ix  ,iy+k2d,iz)
     &                 - ex(ix  ,iy    ,iz))*dtdyi
                  end do
               end do
            end do
#endif /*BFILTER*/


c        if(mype.eq.0)write(*,*)'Capture fluxes at the block boundaries'


! Capture fluxes at the block boundaries
        do iflx = 1,nf

        flux_x(iflx,1,1+ng0*k2d:nyb+ng0*k2d,
     .                1+ng0*k3d:nzb+ng0*k3d,lb)= 
     .                fx(iflx,  1+nguard,
     .                1+nguard*k2d:nyb+nguard*k2d,
     .                1+nguard*k3d:nzb+nguard*k3d) 
        flux_x(iflx,2,1+ng0*k2d:nyb+ng0*k2d,
     .                1+ng0*k3d:nzb+ng0*k3d,lb)= 
     .                fx(iflx,  nxb+1+nguard,
     .                1+nguard*k2d:nyb+nguard*k2d,
     .                1+nguard*k3d:nzb+nguard*k3d) 
        flux_y(iflx,  1+ng0:nxb+ng0,
     .                1,
     .                1+ng0*k3d:nzb+ng0*k3d,lb)=
     .                fy(iflx,  1+nguard:nxb+nguard,
     .                1+nguard*k2d,
     .                1+nguard*k3d:nzb+nguard*k3d) 
        flux_y(iflx,  1+ng0:nxb+ng0,
     .                2,
     .                1+ng0*k3d:nzb+ng0*k3d,lb)=
     .                fy(iflx,  1+nguard:nxb+nguard,
     .                nyb+(1+nguard)*k2d,
     .                1+nguard*k3d:nzb+nguard*k3d) 
#if N_DIM == 3
        flux_z(iflx,  1+ng0:nxb+ng0,
     .                1+ng0*k2d:nyb+ng0*k2d,
     .                1,lb) = 
     .                fz(iflx,  1+nguard:nxb+nguard,
     .                1+nguard*k2d:nyb+nguard*k2d,
     .                1+nguard*k3d)
        flux_z(iflx,  1+ng0:nxb+ng0,
     .                1+ng0*k2d:nyb+ng0*k2d,
     .                2,lb) = 
     .                fz(iflx,  1+nguard:nxb+nguard,
     .                1+nguard*k2d:nyb+nguard*k2d,
     .                nzb+(1+nguard)*k3d)
#endif /*N_DIM == 3*/       
        enddo

#ifdef BFILTER
c       if(mype.eq.0)write(*,*)'entering a:BFILTER'
! Capture fluxes at the block boundaries
        bedge_facex_y(1,1,1+ng0*k2d:nyb+(ng0  )*k2d,
     .                    1+ng0*k3d:nzb+(ng0+1)*k3d,lb)=
     .    ey(1+nguard,1+nguard*k2d:nyb+(nguard  )*k2d,
     .                1+nguard*k3d:nzb+(nguard+1)*k3d)

        bedge_facex_y(1,2,1+ng0*k2d:nyb+(ng0  )*k2d,
     .                    1+ng0*k3d:nzb+(ng0+1)*k3d,lb)=
     .    ey(nxb+1+nguard,1+nguard*k2d:nyb+(nguard  )*k2d,
     .                    1+nguard*k3d:nzb+(nguard+1)*k3d)

        bedge_facex_z(1,1,1+ng0*k2d:nyb+(ng0+1)*k2d,
     .                    1+ng0*k3d:nzb+(ng0  )*k3d,lb)=
     .    ez(1+nguard,1+nguard*k2d:nyb+(nguard+1)*k2d,
     .                1+nguard*k3d:nzb+(nguard  )*k3d)
        bedge_facex_z(1,2,1+ng0*k2d:nyb+(ng0+1)*k2d,
     .                    1+ng0*k3d:nzb+(ng0  )*k3d,lb)=
     .    ez(nxb+1+nguard,1+nguard*k2d:nyb+(nguard+1)*k2d,
     .                    1+nguard*k3d:nzb+(nguard  )*k3d)

        bedge_facey_x(1,1+ng0:nxb+(ng0  ),1,
     .                  1+ng0*k3d:nzb+(ng0+1)*k3d,lb)=
     .    ex(1+nguard:nxb+(nguard ),1+nguard,
     .                1+nguard*k3d:nzb+(nguard+1)*k3d)
        bedge_facey_x(1,1+ng0:nxb+(ng0  ),2,
     .                  1+ng0*k3d:nzb+(ng0+1)*k3d,lb)=
     .    ex(1+nguard:nxb+(nguard ),
     .                nyb+1+nguard,
     .                1+nguard*k3d:nzb+(nguard+1)*k3d)
        bedge_facey_z(1,1+ng0:nxb+ng0+1,1,
     .                1+ng0*k3d:nzb+(ng0  )*k3d,lb)=
     .    ez(1+nguard:nxb+nguard+1,
     .                1+nguard,
     .                1+nguard*k3d:nzb+(nguard  )*k3d)
        bedge_facey_z(1,1+(ng0  ):nxb+ng0+1,2,
     .                1+ng0*k3d:nzb+(ng0  )*k3d,lb)=
     .    ez(1+nguard:nxb+nguard+1,nyb+1+nguard,
     .                1+nguard*k3d:nzb+(nguard  )*k3d)

        bedge_facez_x(1,1+ng0:nxb+(ng0  ),
     .                1+ng0*k2d:nyb+(ng0+1)*k2d,1,lb)=
     .    ex(1+nguard:nxb+(nguard  ),
     .                1+nguard*k2d:nyb+(nguard+1)*k2d,
     .                1+nguard*k3d)
        bedge_facez_x(1,1+ng0:nxb+(ng0  ),
     .                1+ng0*k2d:nyb+(ng0+1)*k2d,2,lb)=
     .    ex(1+nguard:nxb+(nguard  ),
     .                1+nguard*k2d:nyb+(nguard+1)*k2d,
     .                nzb+(nguard+1)*k3d)

        bedge_facez_y(1,1+ng0:nxb+ng0+1,
     .                1+ng0*k2d:nyb+(ng0  )*k2d,1,lb)=
     .    ey(1+nguard:nxb+nguard+1,
     .                1+nguard*k2d:nyb+(nguard  )*k2d,
     .                1+nguard*k3d)
        bedge_facez_y(1,1+ng0:nxb+ng0+1,
     .                1+ng0*k2d:nyb+(ng0  )*k2d,2,lb)=
     .    ey(1+nguard:nxb+nguard+1,
     .                1+nguard*k2d:nyb+(nguard  )*k2d,
     .                nzb+(nguard+1)*k3d)

#endif /*BFILTER*/



#ifdef NO_PERMANENT_GUARDCELLS
! Store new solution for this block in secondary copy
c       if(mype.eq.0)write(*,*)'entering NO_PERMANENT_GUARDCELLS'
        idest=1
        call amr_1blk_to_perm( lcc,lfc,lec,lnc,lb,iopt,idest)


#endif

c       if(mype.eq.2.and.lb.eq.1) then
c          write(*,*) mype, time,'adv2',unk(1,7,1,1,1),unk1(1,11,5,1,1)
c       endif

c        if(mype.eq.0.and.lb.eq.5) then
c           write(*,*) mype, time,'adv2',unk(1,2,1,1,5),unk1(1,6,5,1,1)
c        endif




#ifndef ADVANCE_ALL_LEVELS
      endif
#endif
      enddo             ! end loop over grid blocks
      endif

!     
c      if(mype.eq.0)write(*,*)'mk41: advance_soln: ',istep,mype

#ifdef MPI_USED
        call mpi_barrier (MPI_COMM_WORLD, errcode)
#else
        call shmem_barrier_all()
#endif /*MPI_USED*/

!--------------------
! Modify block boundary fluxes to ensure conservation
      nsub = 1
c       if(mype.eq.0)write(*,*)'entrying  amr_flux_conserve'
       call amr_flux_conserve(mype,nsub)
c       if(mype.eq.0)write(*,*)'leaving  amr_flux_conserve'

#ifdef BFILTER
      lfullblock = .false.
c      if(mype.eq.0)write(*,*)'entrying  amr_edge_average'

      call amr_edge_average(mype,lfullblock,nsub)
c      if(mype.eq.0)write(*,*)'leaving  amr_edge_average'
#endif /*BFILTER*/


#ifdef MPI_USED
        call mpi_barrier (MPI_COMM_WORLD, errcode)
#else
        call shmem_barrier_all()
#endif /*MPI_USED*/


!--------------------
! Apply any changes to the block boundary fluxes which may have 
! been made by AMR_FLUX_CONSERVE to enforce conservation.
! loop over leaf grid blocks

      if(lnblocks.gt.0) then
       do lb=1,lnblocks
#ifndef ADVANCE_ALL_LEVELS
        if(nodetype(lb).eq.1) then
#endif

! Adjust the solution at the block boundaries using 
! the corrected conservative fluxes.
          call grid_variables(lb,dt)

          ng0 = nguard*npgs



! Loop over the stored fluxes, applying them to correct the appropriate
! solution variables.
c          write(*,*)'rvx,rvy,rvz=',rvx,rvy,rvz
          do iflx = 1,ne
c          write(*,*)'iflx,ng0,nxb,lb=',iflx,ng0,nxb,lb

! modify flux/variable association for the current timestep 
c            ivar = iflx + ne
            ivar = iflx

            unk(ivar,1+ng0,:,:,lb) = 
     .                        unk(ivar,1+ng0,:,:,lb)
     .                    - ( tflux_x(iflx,1,:,:,lb)
     .                       - flux_x(iflx,1,:,:,lb) )*rvx
            unk(ivar,nxb+ng0,:,:,lb) = 
     .                        unk(ivar,nxb+ng0,:,:,lb)
     .                    + ( tflux_x(iflx,2,:,:,lb)
     .                       - flux_x(iflx,2,:,:,lb) )*rvx
            if(ndim.ge.2) then
            unk(ivar,:,1+ng0*k2d,:,lb) = 
     .                        unk(ivar,:,1+ng0*k2d,:,lb)
     .                    - ( tflux_y(iflx,:,1,:,lb)
     .                       - flux_y(iflx,:,1,:,lb) )*rvy
            unk(ivar,:,nyb+ng0*k2d,:,lb) = 
     .                        unk(ivar,:,nyb+ng0*k2d,:,lb)
     .                    + ( tflux_y(iflx,:,2,:,lb)
     .                       - flux_y(iflx,:,2,:,lb) )*rvy
            endif
#if N_DIM == 3
            unk(ivar,:,:,1+ng0*k3d,lb) = 
     .                        unk(ivar,:,:,1+ng0*k3d,lb)
     .                    - ( tflux_z(iflx,:,:,1,lb)
     .                       - flux_z(iflx,:,:,1,lb) )*rvz
            unk(ivar,:,:,nzb+ng0*k3d,lb) = 
     .                        unk(ivar,:,:,nzb+ng0*k3d,lb)
     .                    + ( tflux_z(iflx,:,:,2,lb)
     .                       - flux_z(iflx,:,:,2,lb) )*rvz
#endif /*N_DIM == 3*/
          enddo

c           if(mype.eq.0)write(*,*)'entering LORENTZ_FORCE'
 
#ifdef LORENTZ_FORCE
       do ivar=1,3
          lf = ne+ivar
            dmw(ivar,1  +ng0,:,:,lb)=dmw(ivar,1  +ng0,:,:,lb)
     &           -(tflux_x(lf,1,:,:,lb) - flux_x (lf,1,:,:,lb))*rvx
            dmw(ivar,nxb+ng0,:,:,lb)=dmw(ivar,nxb+ng0,:,:,lb)
     &           +(tflux_x(lf,2,:,:,lb) - flux_x (lf,2,:,:,lb))*rvx

            dmw(ivar,:,1  +ng0*k2d,:,lb)=dmw(ivar,:,1  +ng0*k2d,:,lb)
     &           -(tflux_y(lf,:,1,:,lb) - flux_y (lf,:,1,:,lb))*rvy
            dmw(ivar,:,nyb+ng0*k2d,:,lb)=dmw(ivar,:,nyb+ng0*k2d,:,lb)
     &           +(tflux_y(lf,:,2,:,lb) - flux_y (lf,:,2,:,lb))*rvy
            
#if N_DIM == 3
            dmw(ivar,:,:,1  +ng0*k3d,lb)=dmw(ivar,:,:,1  +ng0*k3d,lb)
     &           -(tflux_z(lf,:,:,1,lb) - flux_z (lf,:,:,1,lb))*rvz
            dmw(ivar,:,:,nzb+ng0*k3d,lb)=dmw(ivar,:,:,nzb+ng0*k3d,lb)
     &           +(tflux_z(lf,:,:,2,lb) - flux_z (lf,:,:,2,lb))*rvz
#else
            dmw(3,:,:,1  +ng0*k3d,lb)=0.0
            dmw(3,:,:,nzb+ng0*k3d,lb)=0.0
#endif /*N_DIM == 3*/
      enddo

      call strongfields_wb(mype,lb)

#ifdef PROJECTION
      do iz = 1+ng0*k3d,nzb+ng0*k3d
         jz = nguard*k3d+iz
       do iy = 1+ng0*k2d,nyb+ng0*k2d
          jy = nguard*k2d+iy
         do ix = 1+ng0,nxb+ng0
            jx   = nguard+ix

             bx   = unk(6,ix,iy,iz,lb)+b_sf(1,jx,jy,jz)
             by   = unk(7,ix,iy,iz,lb)+b_sf(2,jx,jy,jz)
             bz   = unk(8,ix,iy,iz,lb)+b_sf(3,jx,jy,jz)
             bsq  = bx*bx + by*by + bz*bz
             bsqi = 1.0/(bsq+1.0e-28)
             bdotd= 0.5*(bx*dmw(1,ix,iy,iz,lb)+dmw(1,ix,iy,iz,lb)*bx)+
     &              0.5*(by*dmw(2,ix,iy,iz,lb)+dmw(2,ix,iy,iz,lb)*by)+
     &              0.5*(bz*dmw(3,ix,iy,iz,lb)+dmw(3,ix,iy,iz,lb)*bz)
             bd   = 0.5*(bdotd*bsqi+bsqi*bdotd) 
!----------------------------------------------------------------------------         
             dmw(1,ix,iy,iz,lb)=dmw(1,ix,iy,iz,lb)-0.5*(bd*bx+bx*bd)
             dmw(2,ix,iy,iz,lb)=dmw(2,ix,iy,iz,lb)-0.5*(bd*by+by*bd)
             dmw(3,ix,iy,iz,lb)=dmw(3,ix,iy,iz,lb)-0.5*(bd*bz+bz*bd)
         enddo
        enddo
       enddo
#endif /*PROJECTION*/

      do iz = 1+ng0*k3d,nzb+ng0*k3d
         jz = nguard*k3d+iz
       do iy = 1+ng0*k2d,nyb+ng0*k2d
          jy = nguard*k2d+iy
         do ix = 1+ng0,nxb+ng0
            jx   = nguard+ix
             rr   = unk(1,ix,iy,iz,lb)
             px   = unk(2,ix,iy,iz,lb)
             py   = unk(3,ix,iy,iz,lb)
             pz   = unk(4,ix,iy,iz,lb)
!
!-----correct momentum for Lorentz force
!
             px   = px+dmw(1,ix,iy,iz,lb)
             py   = py+dmw(2,ix,iy,iz,lb)
             pz   = pz+dmw(3,ix,iy,iz,lb)
!
!     now correct momentum 
!
             unk(2,ix,iy,iz,lb)=px
             unk(3,ix,iy,iz,lb)=py
             unk(4,ix,iy,iz,lb)=pz
!-----------------------------------------------------------------------
         enddo
        enddo
       enddo
#endif /*LORENTZ_FORCE*/


#ifdef BFILTER
c       if(mype.eq.0)write(*,*)'entering c:BFILTER'
!------------------------------
!
! Final step, 
! apply patched up E-fields to advance solution at block boundaries.



!  The changes on the x faces from flux_conserve:

!  F = Bx dy dz + By dz dx + Bz dx dy + Ex dx dt + Ey dy dt + Ez dz dt
!  d F = (div(B)) dx dy dz + (Bx,t + Ez,y - Ey,z) dy dz dt
!                          + (By,t + Ex,z - Ez,x) dz dx dt
!                          + (Bz,t + Ey,x - Ex,y) dx dy dt

!  0 = Bx,t + Ez,y - Ey,z = Bx,t + (Vy Bx - Vx By),y + (Vz Bx - Vx Bz),z
!  The changes on the x faces from flux_conserve:
!  Correct the x faces



! face 1

          do k=kl_bnd+ng0*k3d,ku_bnd-ng0*k3d
            do j=jl_bnd+ng0*k2d,ju_bnd-ng0*k2d
              facevarx(1,1+ng0,j,k,lb) =
     .        facevarx(1,1+ng0,j,k,lb) + (
     .         ( ( bedge_facex_z(1,1,j    ,k    ,lb)
     .            -bedge_facex_z(1,1,j+k2d,k    ,lb) )*dtdyi
     .          -( bedge_facex_y(1,1,j    ,k    ,lb)
     .            -bedge_facex_y(1,1,j    ,k+k3d,lb) )*dtdzi )
     .       - ( (tbedge_facex_z(1,1,j    ,k    ,lb)
     .           -tbedge_facex_z(1,1,j+k2d,k    ,lb) )*dtdyi
     .          -(tbedge_facex_y(1,1,j    ,k    ,lb)
     .           -tbedge_facex_y(1,1,j    ,k+k3d,lb) )*dtdzi ) )
            enddo
          enddo


!-----------------------------------------------------------------------
          do k=kl_bnd+ng0*k3d,ku_bnd-ng0*k3d
            do j=jl_bnd+ng0*k2d+k2d,ju_bnd-ng0*k2d
              facevary(1,1+ng0,j,k,lb) = facevary(1,1+ng0,j,k,lb)-
     &      (bedge_facex_z(1,1,j,k,lb)-tbedge_facex_z(1,1,j,k,lb))*dtdxi
            enddo
          enddo
          do k=kl_bnd+ng0*k3d+k3d,ku_bnd-ng0*k3d
            do j=jl_bnd+ng0*k2d,ju_bnd-ng0*k2d
              facevarz(1,1+ng0,j,k,lb) = facevarz(1,1+ng0,j,k,lb)+
     &      (bedge_facex_y(1,1,j,k,lb)-tbedge_facex_y(1,1,j,k,lb))*dtdxi
            enddo
          enddo
!-----------------------------------------------------------------------
! face 2
          do k=kl_bnd+ng0*k3d,ku_bnd-ng0*k3d
            do j=jl_bnd+ng0*k2d,ju_bnd-ng0*k2d
              facevarx(1,nxb+1+ng0,j,k,lb) =
     &        facevarx(1,nxb+1+ng0,j,k,lb) + (
     &         ( ( bedge_facex_z(1,2,j    ,k    ,lb)
     &            -bedge_facex_z(1,2,j+k2d,k    ,lb) )*dtdyi
     &          -( bedge_facex_y(1,2,j    ,k    ,lb)
     &            -bedge_facex_y(1,2,j    ,k+k3d,lb) )*dtdzi)
     &       - ( (tbedge_facex_z(1,2,j    ,k    ,lb)
     &           -tbedge_facex_z(1,2,j+k2d,k    ,lb) )*dtdyi
     &          -(tbedge_facex_y(1,2,j    ,k    ,lb)
     &           -tbedge_facex_y(1,2,j    ,k+k3d,lb) )*dtdzi)
     &                                                  )
            enddo
          enddo
!-----------------------------------------------------------------------
          do k=kl_bnd+ng0*k3d,ku_bnd-ng0*k3d
           do j=jl_bnd+ng0*k2d+k2d,ju_bnd-ng0*k2d
            facevary(1,nxb+ng0,j,k,lb)=facevary(1,nxb+ng0,j,k,lb)+
     &      (bedge_facex_z(1,2,j,k,lb)-tbedge_facex_z(1,2,j,k,lb))*dtdxi
           enddo
          enddo
          do k=kl_bnd+ng0*k3d+k3d,ku_bnd-ng0*k3d
           do j=jl_bnd+ng0*k2d,ju_bnd-ng0*k2d
            facevarz(1,nxb+ng0,j,k,lb)=facevarz(1,nxb+ng0,j,k,lb)-
     &      (bedge_facex_y(1,2,j,k,lb)-tbedge_facex_y(1,2,j,k,lb))*dtdxi
           enddo
          enddo
!-----------------------------------------------------------------------
! face 3 
          do k=kl_bnd+ng0*k3d,ku_bnd-ng0*k3d
              do i=il_bnd+ng0,iu_bnd-ng0
              facevary(1,i,1+ng0,k,lb) =facevary(1,i,1+ng0,k,lb)+(
     &        ( ( bedge_facey_x(1,i  ,1,k    ,lb)
     &           -bedge_facey_x(1,i  ,1,k+k3d,lb) )*dtdzi
     &         -( bedge_facey_z(1,i  ,1,k    ,lb)
     &           -bedge_facey_z(1,i+1,1,k    ,lb) )*dtdxi )
     &      - ( (tbedge_facey_x(1,i  ,1,k    ,lb)
     &          -tbedge_facey_x(1,i  ,1,k+k3d,lb) )*dtdzi
     &         -(tbedge_facey_z(1,i  ,1,k    ,lb)
     &          -tbedge_facey_z(1,i+1,1,k    ,lb) )*dtdxi )
     &                                                    )
              enddo
          enddo
!-----------------------------------------------------------------------
          do k=kl_bnd+ng0*k3d,ku_bnd-ng0*k3d
           do i=il_bnd+ng0+1,iu_bnd-ng0
            facevarx(1,i,1+ng0,k,lb)=facevarx(1,i,1+ng0,k,lb) + 
     &      (bedge_facey_z(1,i,1,k,lb)-tbedge_facey_z(1,i,1,k,lb))*dtdxi
           enddo
          enddo
          do k=kl_bnd+ng0*k3d+k3d,ku_bnd-ng0*k3d
           do i=il_bnd+ng0,iu_bnd-ng0
            facevarz(1,i,1+ng0,k,lb) =facevarz(1,i,1+ng0,k,lb) - 
     &      (bedge_facey_x(1,i,1,k,lb)-tbedge_facey_x(1,i,1,k,lb))*dtdzi
           enddo
          enddo
!-----------------------------------------------------------------------
! face 4 
          do k=kl_bnd+ng0*k3d,ku_bnd-ng0*k3d
            do i=il_bnd+ng0,iu_bnd-ng0
             facevary(1,i,nyb+1+ng0,k,lb)=
     &       facevary(1,i,nyb+1+ng0,k,lb)+(
     &       ( ( bedge_facey_x(1,i  ,2,k    ,lb)
     &          -bedge_facey_x(1,i  ,2,k+k3d,lb) )*dtdzi
     &        -( bedge_facey_z(1,i  ,2,k    ,lb)
     &          -bedge_facey_z(1,i+1,2,k    ,lb) )*dtdxi )
     &      - ((tbedge_facey_x(1,i  ,2,k    ,lb)
     &         -tbedge_facey_x(1,i  ,2,k+k3d,lb) )*dtdzi
     &       -( tbedge_facey_z(1,i  ,2,k    ,lb)
     &         -tbedge_facey_z(1,i+1,2,k    ,lb) )*dtdxi )
     &                                                       )
              enddo
          enddo
!-----------------------------------------------------------------------
          do k=kl_bnd+ng0*k3d,ku_bnd-ng0*k3d
           do i=il_bnd+ng0+1,iu_bnd-ng0
            facevarx(1,i,nyb+ng0,k,lb) =facevarx(1,i,nyb+ng0,k,lb) - 
     &      (bedge_facey_z(1,i,2,k,lb)-tbedge_facey_z(1,i,2,k,lb))*dtdxi
           enddo
          enddo
          do k=kl_bnd+ng0*k3d+k3d,ku_bnd-ng0*k3d
           do i=il_bnd+ng0,iu_bnd-ng0
            facevarz(1,i,nyb+ng0,k,lb) =facevarz(1,i,nyb+ng0,k,lb) + 
     &      (bedge_facey_x(1,i,2,k,lb)-tbedge_facey_x(1,i,2,k,lb))*dtdzi
           enddo
          enddo


          if(ndim.eq.3.or.l2p5d.eq.1) then

! face 5 
          do j=jl_bnd+ng0*k2d,ju_bnd-ng0*k2d
           do i=il_bnd+ng0,iu_bnd-ng0
            facevarz(1,i,j,1+ng0,lb) =
     .                facevarz(1,i,j,1+ng0,lb) + (
     .        ( ( bedge_facez_y(1,i,j,1,lb)
     .           -bedge_facez_y(1,i+1,j,1,lb) )*dtdxi
     .         -( bedge_facez_x(1,i,j,1,lb)
     .           -bedge_facez_x(1,i,j+k2d,1,lb) )*dtdyi )
     .      - ( ( tbedge_facez_y(1,i,j,1,lb)
     .           -tbedge_facez_y(1,i+1,j,1,lb) )*dtdxi
     .         -( tbedge_facez_x(1,i,j,1,lb)
     .           -tbedge_facez_x(1,i,j+k2d,1,lb) )*dtdyi )
     .                                                       )
            enddo
          enddo
!-----------------------------------------------------------------------
          do j=jl_bnd+ng0*k2d,ju_bnd-ng0*k2d
           do i=il_bnd+ng0+1,iu_bnd-ng0
            facevarx(1,i,j,1+ng0,lb)=facevarx(1,i,j,1+ng0,lb) - 
     &     (bedge_facez_y(1,i,j,1,lb)-tbedge_facez_y(1,i,j,1,lb))*dtdxi
           enddo
          enddo
          do j=jl_bnd+ng0*k2d+k2d,ju_bnd-ng0*k2d
           do i=il_bnd+ng0,iu_bnd-ng0
            facevary(1,i,j,1+ng0,lb)=facevary(1,i,j,1+ng0,lb) + 
     &      (bedge_facez_x(1,i,j,1,lb)-tbedge_facez_x(1,i,j,1,lb))*dtdxi
           enddo
          enddo

          endif


! face 6 
          if(ndim.eq.3) then
          do j=jl_bnd+ng0*k2d,ju_bnd-ng0*k2d
           do i=il_bnd+ng0,iu_bnd-ng0
            facevarz(1,i,j,nzb+1+ng0,lb) =
     .      facevarz(1,i,j,nzb+1+ng0,lb) + (
     .        ((bedge_facez_y(1,i  ,j    ,2,lb)-
     .          bedge_facez_y(1,i+1,j    ,2,lb))*dtdxi
     .        -(bedge_facez_x(1,i  ,j    ,2,lb)
     .         -bedge_facez_x(1,i  ,j+k2d,2,lb))*dtdyi )
     .      -((tbedge_facez_y(1,i  ,j    ,2,lb)
     .        -tbedge_facez_y(1,i+1,j    ,2,lb))*dtdxi
     .       -(tbedge_facez_x(1,i  ,j    ,2,lb)
     .        -tbedge_facez_x(1,i  ,j+k2d,2,lb))*dtdyi )
     .                                                 )

            enddo
          enddo
!-----------------------------------------------------------------------
          do j=jl_bnd+ng0*k2d,ju_bnd-ng0*k2d
           do i=il_bnd+ng0+1,iu_bnd-ng0
            facevarx(1,i,j,nzb+ng0,lb)=facevarx(1,i,j,nzb+ng0,lb) + 
     .      (bedge_facez_y(1,i,j,2,lb)-tbedge_facez_y(1,i,j,2,lb))*dtdxi
           enddo
          enddo
         do j=jl_bnd+ng0*k2d+k2d,ju_bnd-ng0*k2d
           do i=il_bnd+ng0,iu_bnd-ng0
            facevary(1,i,j,nzb+ng0,lb)=facevary(1,i,j,nzb+ng0,lb) - 
     &      (bedge_facez_x(1,i,j,2,lb)-tbedge_facez_x(1,i,j,2,lb))*dtdxi
           end do
          end do

         endif


#ifdef OVERWRITE
          lflag             = earth_blocks(lb)
c          write(*,*)'entering overwrite'
! Overwrite magnetic fields at cell centers from those on faces 
         do iz= 1,ku_bnd
          do iy= 1,ju_bnd
           do ix= 1,iu_bnd
            b_sq_o(ix,iy,iz)=unk(6,ix,iy,iz,lb)**2+unk(7,ix,iy,iz,lb)**2 
     $                      +unk(8,ix,iy,iz,lb)**2
            unk(6,ix,iy,iz,lb)= 
     $         0.5*( facevarx(1,ix+1,iy    ,iz    ,lb)
     $             + facevarx(1,ix  ,iy    ,iz    ,lb))
            unk(7,ix,iy,iz,lb)= 
     $         0.5*( facevary(1,ix  ,iy+k2d,iz    ,lb)
     $             + facevary(1,ix  ,iy    ,iz    ,lb))
#if N_DIM == 3
            unk(8,ix,iy,iz,lb)= 
     $         0.5*( facevarz(1,ix  ,iy    ,iz+k3d,lb)
     $             + facevarz(1,ix  ,iy    ,iz    ,lb))
#else
!mk0
!           unk(4,ix,iy,iz,lb)= 0.0
!           unk(8,ix,iy,iz,lb)= 0.0

#endif /*N_DIM == 3*/
            b_sq_n(ix,iy,iz)=unk(6,ix,iy,iz,lb)**2+unk(7,ix,iy,iz,lb)**2 
     $                      +unk(8,ix,iy,iz,lb)**2
           end do
          end do
         end do
c       write(*,*)'u3(7,1,1,1),byf(1,1,1,1),byf(1,2,1,1)=',
c     &  unk(7,1,1,1,lb),facevary(1,1,1,1,lb),facevary(1,1,2,1,lb)
#ifdef INTERIOR_BOUNDARY1	 
         if(lflag)then
          irzones            = irdzones(lb) 
          irzone(:,:irzones)  = p_rzone(:,:irzones,lb) 
c           write(*,*)'lb,lflag,izones=',lb,lflag,izones
          do iz=1,irzones
           i = irzone(1,iz)
           j = irzone(2,iz)
           k = irzone(3,iz)
           b_sq_o(i,j,k)     =0.0
           unk(2:4,i,j,k,lb) =0.0
           unk(6:8,i,j,k,lb) =0.0
           b_sq_n(i,j,k)     =0.0
          enddo
         endif
#endif /*INTERIOR_BOUNDARY*/

         do iz= 1,ku_bnd
          do iy= 1,ju_bnd
           do ix= 1,iu_bnd
            unk(5,ix,iy,iz,lb)=unk(5,ix,iy,iz,lb)
     $        +(b_sq_n(ix,iy,iz) - b_sq_o(ix,iy,iz))*pie8i
           end do
          end do
         end do

#endif /*OVERWRITE*/
#endif /*BFILTER*/

#ifndef ADVANCE_ALL_LEVELS
      endif
#endif
      enddo             ! end loop over grid blocks
      endif

      call shmem_barrier_all()

!--------------------
!
!-----solution at time level n+1 is now stored at ne+1=>2*ne
!

       if(lnblocks.gt.0) then
        do lb=1,lnblocks
#ifndef ADVANCE_ALL_LEVELS
         if(nodetype(lb).eq.1) then
#endif
#ifdef CAVITATION_FLOOR
         do k=1,nzb
          do j=1,nyb
           do i=1,nxb
            if(unk(1,i,j,k,lb).le.rmin)then
             unk(1,i,j,k,lb)=rmin 
             unk(2,i,j,k,lb)=0.0 
             unk(3,i,j,k,lb)=0.0 
             unk(4,i,j,k,lb)=0.0 
             be             =(unk(6,i,j,k,lb)**2+unk(7,i,j,k,lb)**2+
     &                        unk(8,i,j,k,lb)**2)
             be             =be*pie8i
             unk(5,i,j,k,lb)=pmin*gammam1i+be
            endif
           enddo
          enddo
         enddo
#endif /*CAVITATION_FLOOR*/
#ifdef DEBUGGXX
         do j = 1, ne
          umax(j,lb) = maxval( unk(j,1:nxb,1:nyb,1:nzb,lb) )
          umin(j,lb) = minval( unk(j,1:nxb,1:nyb,1:nzb,lb) )
         enddo
#endif /*DEBUGGXX*/
#ifndef ADVANCE_ALL_LEVELS
         endif
#endif
       enddo
      endif  /* lnblocks>0 */
#ifdef DEBUGGXX

      if(lnblocks.gt.0) then
       do lb=1,lnblocks
#ifndef ADVANCE_ALL_LEVELS
        if(nodetype(lb).eq.1) then
#endif
         do j=1,ne
          uarmax(j) = maxval( umax(j,1:lnblocks) )
          uarmin(j) = minval( umin(j,1:lnblocks) )
         enddo

#ifndef ADVANCE_ALL_LEVELS
        endif
#endif
       enddo
      endif /* lnblocks>0 */
      call shmem_barrier_all()
      do j=1,ne
        call comm_real_min_to_all(uarmin(j),uarmin(j))
        call comm_real_max_to_all(uarmax(j),uarmax(j))
      enddo 
        call shmem_barrier_all()
       do j=1,ne
c       write(*,*)'uu(',j,'),uu(',j,') = ',
c     &               uarmax(j),uarmin(j)
       enddo
#endif /*DEBUGG*/

#ifdef MPI_USED
        call mpi_barrier (MPI_COMM_WORLD, errcode)
#else
        call shmem_barrier_all()
#endif /*MPI_USED*/



! Increment time
      time  = time + dt
      dtold =  dt


c      write(*,*)'mk42: istep,mype',istep,mype
      return
      end subroutine advance_soln

      subroutine ch_data(mype,ib)
      use paramesh_dimensions
      use physicaldata
      use strongfields
      

      implicit none
      integer :: mype,ib
      integer iunit1, iunit2,iunit3, i, j, n
      integer is1, ie1, is2, ie2
      if(mype.eq.0.and.ib.eq.1) then
         iunit1=70; iunit2=71; iunit3=72;
         open(unit=iunit1,file='out_unk1_wo_bfilter')
         open(unit=iunit2,file='out_facvary_wo_bfilter')
         open(unit=iunit3,file='out_sf_wo_bfilter')

         is1=1;ie1=2*nguard;  is2=nxb+1; ie2=nxb+2*nguard;
         do n = 7, 7
            if(n.eq.1)write(iunit1,*)'-------- rho --------'
            if(n.eq.2)write(iunit1,*)'-------- u   --------'
            if(n.eq.3)write(iunit1,*)'-------- v   --------'
            if(n.eq.4)write(iunit1,*)'-------- w   --------'
            if(n.eq.5)write(iunit1,*)'-------- p   --------'
            if(n.eq.6)write(iunit1,*)'-------- bx  --------'
            if(n.eq.7)write(iunit1,*)'-------- by  --------'
            if(n.eq.8)write(iunit1,*)'-------- bz  --------'

            if(n.eq.1)write(iunit2,*)'-------- rho --------'
            if(n.eq.2)write(iunit2,*)'-------- u   --------'
            if(n.eq.3)write(iunit2,*)'-------- v   --------'
            if(n.eq.4)write(iunit2,*)'-------- w   --------'
            if(n.eq.5)write(iunit2,*)'-------- p   --------'
            if(n.eq.6)write(iunit2,*)'-------- bx  --------'
            if(n.eq.7)write(iunit2,*)'-------- by  --------'
            if(n.eq.8)write(iunit2,*)'-------- bz  --------'
            do j = 1,nguard 
              write(iunit1,100) (unk1(n,i,j,1,1),i=is1,ie1),
     &                          (unk1(n,i,j,1,1),i=is2,ie2)
              write(iunit2,100) (facevary1(n,i,j,1,1),i=is1,ie1),
     &                          (facevary1(n,i,j,1,1),i=is2,ie2)
              write(iunit3,100) (byf_sf(1,i,j,1),i=is1,ie1),
     &                          (byf_sf(1,i,j,1),i=is2,ie2)
            end do
            write(iunit1,*)' '
            write(iunit2,*)' '
            write(iunit3,*)' '
            do j = nguard+1,2*nguard 
              write(iunit1,100) (unk1(n,i,j,1,1),i=is1,ie1),
     &                          (unk1(n,i,j,1,1),i=is2,ie2)
              write(iunit2,100) (facevary1(n,i,j,1,1),i=is1,ie1),
     &                          (facevary1(n,i,j,1,1),i=is2,ie2)
              write(iunit3,100) (byf_sf(1,i,j,1),i=is1,ie1),
     &                          (byf_sf(1,i,j,1),i=is2,ie2)
            end do
            write(iunit1,*)' '
            write(iunit2,*)' '
            write(iunit3,*)' '
            do j = nyb+1,nyb+nguard 
              write(iunit1,100) (unk1(n,i,j,1,1),i=is1,ie1),
     &                          (unk1(n,i,j,1,1),i=is2,ie2)
              write(iunit2,100) (facevary1(n,i,j,1,1),i=is1,ie1),
     &                          (facevary1(n,i,j,1,1),i=is2,ie2)
              write(iunit3,100) (byf_sf(1,i,j,1),i=is1,ie1),
     &                          (byf_sf(1,i,j,1),i=is2,ie2)
            end do
            write(iunit1,*)' '
            write(iunit2,*)' '
            write(iunit3,*)' '
            do j = nyb+nguard+1,nyb+2*nguard 
              write(iunit1,100) (unk1(n,i,j,1,1),i=is1,ie1),
     &                          (unk1(n,i,j,1,1),i=is2,ie2)
              write(iunit2,100) (facevary1(n,i,j,1,1),i=is1,ie1),
     &                          (facevary1(n,i,j,1,1),i=is2,ie2)
              write(iunit3,100) (byf_sf(1,i,j,1),i=is1,ie1),
     &                          (byf_sf(1,i,j,1),i=is2,ie2)
            end do
            write(iunit1,*)' '
            write(iunit2,*)' '
            write(iunit3,*)' '
             write(iunit2,100) 
     &          (facevary1(n,i,nyb+2*nguard+1,1,1),i=is1,ie1),
     &          (facevary1(n,i,nyb+2*nguard+1,1,1),i=is2,ie2)
            write(iunit1,*)' ' 
            write(iunit2,*)' '
            write(iunit3,*)' '
         end do
 100     format(4f8.4,"|",4f8.4,"|",4f8.4,"|",4f8.4)

         close(iunit1)
         close(iunit2)
         close(iunit3)
      endif

      return
      end

      subroutine ch2_data(mype,istep,time,dt)
      use paramesh_dimensions
      use physicaldata
      use strongfields
      

      implicit none
      integer :: istep
      real    :: time,dt
      integer iunit30, iunit31, iunit32, iunit33, iunit34, iunit35
      integer iunit40, iunit41, iunit42, iunit43, iunit44, iunit45
      integer i, j, n
      integer is1, ie1, is2, ie2

      integer :: mype,iopt,nlayers,lb,icoord
      logical :: lcc,lfc,lec,lnc,l_srl_only,ldiag


!for NO_PERMANENT_GUARDCELLS
        iopt       = 1
        nlayers    = nguard
        lcc        = .true.
        lfc        = .false.
        l_srl_only = .false.
        icoord     = 0
        ldiag      = .true.
        lec        =.false.
        lnc        =.false.

!lb=1 and 3
      call amr_1blk_copy_soln(-1)
      lb=3
      call amr_1blk_guardcell(mype,iopt,nlayers,lb,mype,lcc,lfc,
     &                        lec,lnc,
     &                        l_srl_only,icoord,ldiag)

!---------------------------------------------------------
      iunit30=30; iunit31=31; iunit32=32; iunit33=33; 
      iunit34=34; iunit35=35;

      write(iunit30,*)'----------------------------------------------'
      write(iunit31,*)'----------------------------------------------'
      write(iunit32,*)'----------------------------------------------'
      write(iunit33,*)'----------------------------------------------'
      write(iunit34,*)'----------------------------------------------'
      write(iunit35,*)'----------------------------------------------'
      write(iunit30,*)'istep,time,dt = ',istep,time,dt
      write(iunit31,*)'istep,time,dt = ',istep,time,dt
      write(iunit32,*)'istep,time,dt = ',istep,time,dt
      write(iunit33,*)'istep,time,dt = ',istep,time,dt
      write(iunit34,*)'istep,time,dt = ',istep,time,dt
      write(iunit35,*)'istep,time,dt = ',istep,time,dt

      is1=nxb/2+1; ie1=nxb;
      is2=1; ie2=1;
!is2=1; ie2=nxb/2;
      do j = 1,nyb 
         write(iunit30,100) (unk(1,i,j,1,1),i=is1,ie1),
     &                      (unk(1,i,j,1,1),i=is2,ie2)
         write(iunit31,100) (unk(2,i,j,1,1),i=is1,ie1),
     &                      (unk(2,i,j,1,1),i=is2,ie2)
         write(iunit32,100) (unk(3,i,j,1,1),i=is1,ie1),
     &                      (unk(3,i,j,1,1),i=is2,ie2)
         write(iunit33,100) (unk(5,i,j,1,1),i=is1,ie1),
     &                      (unk(5,i,j,1,1),i=is2,ie2)
         write(iunit34,100) (unk(6,i,j,1,1),i=is1,ie1),
     &                      (unk(6,i,j,1,1),i=is2,ie2)
         write(iunit35,100) (unk(7,i,j,1,1),i=is1,ie1),
     &                      (unk(7,i,j,1,1),i=is2,ie2)
      end do

!---------------------------------------------------------
!     iunit40=40; iunit41=41; iunit42=42; iunit43=43;
!     iunit44=44; iunit45=45;

!     write(iunit40,*)'----------------------------------------------'
!     write(iunit41,*)'----------------------------------------------'
!     write(iunit42,*)'----------------------------------------------'
!     write(iunit43,*)'----------------------------------------------'
!     write(iunit44,*)'----------------------------------------------'
!     write(iunit45,*)'----------------------------------------------'
!     write(iunit40,*)'istep,time,dt = ',istep,time,dt
!     write(iunit41,*)'istep,time,dt = ',istep,time,dt
!     write(iunit42,*)'istep,time,dt = ',istep,time,dt
!     write(iunit43,*)'istep,time,dt = ',istep,time,dt
!     write(iunit44,*)'istep,time,dt = ',istep,time,dt
!     write(iunit45,*)'istep,time,dt = ',istep,time,dt

!     is1=nxb/2+1; ie1=nxb;
!     is2=1; ie2=nguard;
!     do j = 1,nyb 
!        write(iunit40,100) (unk(1,i,j,1,1),i=is1,ie1),
!    &                      (unk1(1,i,j+nguard,1,1),i=is2,ie2)
!        write(iunit41,100) (unk(2,i,j,1,1),i=is1,ie1),
!    &                      (unk1(2,i,j+nguard,1,1),i=is2,ie2)
!        write(iunit42,100) (unk(3,i,j,1,1),i=is1,ie1),
!    &                      (unk1(3,i,j+nguard,1,1),i=is2,ie2)
!        write(iunit43,100) (unk(5,i,j,1,1),i=is1,ie1),
!    &                      (unk1(5,i,j+nguard,1,1),i=is2,ie2)
!        write(iunit44,100) (unk(6,i,j,1,1),i=is1,ie1),
!    &                      (unk1(6,i,j+nguard,1,1),i=is2,ie2)
!        write(iunit45,100) (unk(7,i,j,1,1),i=is1,ie1),
!    &                      (unk1(7,i,j+nguard,1,1),i=is2,ie2)
!     end do
 100  format(5(1x,f15.11))
 200  format(8(1x,f9.5))

!is1=9; ie1=16
!do j = 1,8 
      return
      end


      subroutine ch3_data(istep,time,dt)
!for 1block
      use paramesh_dimensions
      use physicaldata
      use strongfields
      

      implicit none
      integer :: istep
      real    :: time,dt
      integer iunit30, iunit31, iunit32, iunit33, iunit34, iunit35
      integer i, j
      integer is1, ie1


!---------------------------------------------------------
      iunit30=30; iunit31=31; iunit32=32; iunit33=33; 
      iunit34=34; iunit35=35;

      write(iunit30,*)'----------------------------------------------'
      write(iunit31,*)'----------------------------------------------'
      write(iunit32,*)'----------------------------------------------'
      write(iunit33,*)'----------------------------------------------'
      write(iunit34,*)'----------------------------------------------'
      write(iunit35,*)'----------------------------------------------'
      write(iunit30,*)'1_blk: istep,time,dt = ',istep,time,dt
      write(iunit31,*)'1_blk: istep,time,dt = ',istep,time,dt
      write(iunit32,*)'1_blk: istep,time,dt = ',istep,time,dt
      write(iunit33,*)'1_blk: istep,time,dt = ',istep,time,dt
      write(iunit34,*)'1_blk: istep,time,dt = ',istep,time,dt
      write(iunit35,*)'1_blk: istep,time,dt = ',istep,time,dt

      is1=5; ie1=9;
      do j = 1,8 
         write(iunit30,100) (unk(1,i,j,1,1),i=is1,ie1)
         write(iunit31,100) (unk(2,i,j,1,1),i=is1,ie1)
         write(iunit32,100) (unk(3,i,j,1,1),i=is1,ie1)
         write(iunit33,100) (unk(5,i,j,1,1),i=is1,ie1)
         write(iunit34,100) (unk(6,i,j,1,1),i=is1,ie1)
         write(iunit35,100) (unk(7,i,j,1,1),i=is1,ie1)
      end do

 100  format(5(1x,f15.11))

      return
      end

      subroutine ch4_data(istep,time,dt,ib)
      implicit none
      integer :: istep,ib
      real    :: time,dt

      integer iunit30, iunit31, iunit32, iunit33, iunit34, iunit35
!---------------------------------------------------------
      iunit30=30; iunit31=31; iunit32=32; iunit33=33;
      iunit34=34; iunit35=35;

      write(iunit30,*)'----------------------------------------------'
      write(iunit30,*)'istep,time,dt = ',istep,time,dt

      return
      end

      subroutine ch5_data(fx,fy,fz,ib)

      use paramesh_dimensions
      use physicaldata
      use athena_specific_parameters
      implicit none
      real,dimension(nf,iu_bnd1p,ju_bnd1 ,ku_bnd1 ) :: fx
      real,dimension(nf,iu_bnd1 ,ju_bnd1p,ku_bnd1 ) :: fy
      real,dimension(nf,iu_bnd1 ,ju_bnd1 ,ku_bnd1p) :: fz

      integer :: ib

      integer iunit30, iunit31, iunit32, iunit33, iunit34, iunit35
      integer is1, ie1;
      integer i, j
!---------------------------------------------------------
      iunit30=30; iunit31=31; iunit32=32; iunit33=33;
      iunit34=34; iunit35=35;

      if(ib.ne.1) then
!for lb=3
         do j= 8, 4, -1
            write(iunit30,100) (fx(1,i,j,1),i=5,6),
     &                      fy(1,5,j,1),fy(1,5,j+1,1)
            write(iunit31,100) (fx(2,i,j,1),i=5,6),
     &                      fy(2,5,j,1),fy(2,5,j+1,1)
            write(iunit32,100) (fx(3,i,j,1),i=5,6),
     &                      fy(3,5,j,1),fy(3,5,j+1,1)
            write(iunit33,100) (fx(5,i,j,1),i=5,6),
     &                      fy(5,5,j,1),fy(5,5,j+1,1)
            write(iunit34,100) (fx(6,i,j,1),i=5,6),
     &                      fy(6,5,j,1),fy(6,5,j+1,1)
            write(iunit35,100) (fx(7,i,j,1),i=5,6),
     &                      fy(7,5,j,1),fy(7,5,j+1,1)
         end do
      else
!1blk
         do j= 8, 4, -1
            write(iunit30,100) (fx(1,i,j,1),i=13,14),
     &                      fy(1,13,j,1),fy(1,13,j+1,1)
            write(iunit31,100) (fx(2,i,j,1),i=13,14),
     &                      fy(2,13,j,1),fy(2,13,j+1,1)
            write(iunit32,100) (fx(3,i,j,1),i=13,14),
     &                      fy(3,13,j,1),fy(3,13,j+1,1)
            write(iunit33,100) (fx(5,i,j,1),i=13,14),
     &                      fy(5,13,j,1),fy(5,13,j+1,1)
            write(iunit34,100) (fx(6,i,j,1),i=13,14),
     &                      fy(6,13,j,1),fy(6,13,j+1,1)
            write(iunit35,100) (fx(7,i,j,1),i=13,14),
     &                      fy(7,13,j,1),fy(7,13,j+1,1)
         end do
      endif

 100  format(5(1x,e18.11))

      return
      end


      subroutine ch6_data(fx,fy,fz,xmy_vl,xmy_vr,ymy_vl,ymy_vr,
     &              bxf_my, byf_my, ib)

      use paramesh_dimensions
      use physicaldata
      use athena_specific_parameters
      implicit none
      real,dimension(nf,iu_bnd1p,ju_bnd1 ,ku_bnd1 ) :: fx
      real,dimension(nf,iu_bnd1 ,ju_bnd1p,ku_bnd1 ) :: fy
      real,dimension(nf,iu_bnd1 ,ju_bnd1 ,ku_bnd1p) :: fz

      real,dimension(ne,iu_bnd1p,ju_bnd1p ,ku_bnd1p ) :: xmy_vl
      real,dimension(ne,iu_bnd1p,ju_bnd1p ,ku_bnd1p ) :: xmy_vr
      real,dimension(ne,iu_bnd1p,ju_bnd1p ,ku_bnd1p ) :: ymy_vl
      real,dimension(ne,iu_bnd1p,ju_bnd1p ,ku_bnd1p ) :: ymy_vr

      real,dimension(3,iu_bnd1p,ju_bnd1p ,ku_bnd1p ) :: bxf_my
      real,dimension(3,iu_bnd1p,ju_bnd1p ,ku_bnd1p ) :: byf_my

      integer :: ib

      integer iunit30, iunit31, iunit32, iunit33, iunit34, iunit35
      integer is, js
      integer i, j, m, n
!---------------------------------------------------------
      iunit30=30; iunit31=31; iunit32=32; iunit33=33;
      iunit34=34; iunit35=35;
      
      if(ib.ne.1) then
!for lb=3
         write(iunit30,100)
         is=5
         do m = 1,ne
          write(iunit30,101) m
!do j= 8, 4, -1
          do j= 6, 6, -1
            write(iunit30,102) j
            write(iunit30,103) xmy_vl(m,is,j,1),xmy_vr(m,is,j,1),
     &                xmy_vl(m,is+1,j,1),xmy_vr(m,is+1,j,1) 
            write(iunit30,103) ymy_vl(m,is,j,1),ymy_vr(m,is,j,1),
     &                ymy_vl(m,is,j+1,1),ymy_vr(m,is,j+1,1) 
            write(iunit30,103) (fx(m,i,j,1),i=is,is+1),
     &                      fy(m,is,j,1),fy(m,is,j+1,1)
!           write(iunit30,103) (bxf_my(n,is,j,1),n=1,3)
!           write(iunit30,103) (byf_my(n,is,j,1),n=1,3)
            write(iunit30,104)
          end do 
         end do
      else
!1blk
         write(iunit30,105)
         is=13
         do m = 1,ne
          write(iunit30,101) m
!do j= 8, 4, -1
          do j= 6, 6, -1
            write(iunit30,102) j
            write(iunit30,103) xmy_vl(m,is,j,1),xmy_vr(m,is,j,1),
     &                xmy_vl(m,is+1,j,1),xmy_vr(m,is+1,j,1) 
            write(iunit30,103) ymy_vl(m,is,j,1),ymy_vr(m,is,j,1),
     &                ymy_vl(m,is,j+1,1),ymy_vr(m,is,j+1,1) 
            write(iunit30,103) (fx(m,i,j,1),i=is,is+1),
     &                      fy(m,is,j,1),fy(m,is,j+1,1)
!           write(iunit30,103) (bxf_my(n,is,j,1),n=1,3)
!           write(iunit30,103) (byf_my(n,is,j,1),n=1,3)
            write(iunit30,104)
          end do 
         end do
      endif

 100  format(2x,"multi-block")
 101  format(2x,"variable # = ", i2)
 102  format(2x,"j = ", i2)
 103  format(5(1x,e18.11))
 104  format(" ")
 105  format(2x,"1-block")

      return
      end

      subroutine ch7_data(fx,fy,fz,xmy_vl,xmy_vr,ymy_vl,ymy_vr,
     &              bxf_my, byf_my, ib)

#include "mydebug.fh"
      use paramesh_dimensions
      use physicaldata
      use athena_specific_parameters
      implicit none
      real,dimension(nf,iu_bnd1p,ju_bnd1 ,ku_bnd1 ) :: fx
      real,dimension(nf,iu_bnd1 ,ju_bnd1p,ku_bnd1 ) :: fy
      real,dimension(nf,iu_bnd1 ,ju_bnd1 ,ku_bnd1p) :: fz

      real,dimension(ne,iu_bnd1p,ju_bnd1p ,ku_bnd1p ) :: xmy_vl
      real,dimension(ne,iu_bnd1p,ju_bnd1p ,ku_bnd1p ) :: xmy_vr
      real,dimension(ne,iu_bnd1p,ju_bnd1p ,ku_bnd1p ) :: ymy_vl
      real,dimension(ne,iu_bnd1p,ju_bnd1p ,ku_bnd1p ) :: ymy_vr

      real,dimension(3,iu_bnd1p,ju_bnd1p ,ku_bnd1p ) :: bxf_my
      real,dimension(3,iu_bnd1p,ju_bnd1p ,ku_bnd1p ) :: byf_my

      integer :: ib, lb

      integer iunit30, iunit31, iunit32, iunit33, iunit34, iunit35
      integer is, js
      integer i, j, m, n
!---------------------------------------------------------
      iunit30=30; iunit31=31; iunit32=32; iunit33=33;
      iunit34=34; iunit35=35;

!-------------------------------------------------------------
#ifdef onedebug
!1blk
      write(iunit30,100)
      is=13;
 100   format(2x,"1blk")
#endif

#ifdef multidebug
!multi block
      write(iunit30,200)
      is=5;
 200  format(2x,"multi-block")
#endif
!-------------------------------------------------------------

      do j=6,4,-1
         write(iunit30,201)is,j
!x-component
         write(iunit30,202)
         write(iunit30,203)
         do n = 1, ne
            write(iunit30,204) xmy_vl(n,is,j,1),xmy_vr(n,is,j,1),
     &                            fx(n,is,j,1)
         end do
         write(iunit30,205)
         do n = 1, ne
            write(iunit30,204) xmy_vl(n,is+1,j,1),xmy_vr(n,is+1,j,1),
     &                            fx(n,is+1,j,1)
         end do
         write(iunit30,206)
!y-component
         write(iunit30,207)
         write(iunit30,208)
         do n = 1, ne
            write(iunit30,204) ymy_vl(n,is,j,1),ymy_vr(n,is,j,1),
     &                            fy(n,is,j,1)
         end do
         write(iunit30,206)
         write(iunit30,209)
         do n = 1, ne
            write(iunit30,204) ymy_vl(n,is,j+1,1),ymy_vr(n,is,j+1,1),
     &                            fy(n,is,j+1,1)
         end do
      end do

 201  format(2x,"------------------ I=",i2," J=",i2," ---------------")
 202  format(2x,"  x-component left and right states          flux")
 203  format(2x,"I, J ")
 204  format(5(1x,e18.11))
 205  format(2x,"I + 1, J ")
 206  format(2x,"   ")
 207  format(2x,"  y-component left and right states          flux")
 208  format(2x,"I, J ")
 209  format(2x,"I, J + 1 ")

      return
      end
