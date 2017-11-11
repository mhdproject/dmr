#include "paramesh_preprocessor.fh"
#include "athena_preprocessor.fh"
#undef DEBUG
       subroutine amr_timestep(dt,dtmin,dtmax,mype)
       use physcons
       use paramesh_dimensions
       use physicaldata
       use grid
       use interior_gbc
       use athena_specific_parameters
       use strongfields
! include file defining the tree
      use tree


      use paramesh_interfaces
#ifdef MPI_USED
      use paramesh_mpi_interfaces
#endif /* MPI_USED */


#ifdef MPI_USED
      include "mpif.h"
      integer :: max_blks_sent
#endif /* MPI_USED */

!-----------------------------------------------------------------------
! This routine computes the hydro timestep, and distributes
! it to all the processors.
!
! Written:      Peter MacNeice February 1997
!
!
#include "amr_shmem.fh"
      include "shmem_reduce.fh"

!------------------------------------------
      real    :: dt,dtmin,dtmax
      integer :: mype
! local variables

      real dtl,dtmaxl
      save dtl,dtmaxl
      integer :: i,j,k
      integer npes,shmem_n_pes
      integer, parameter :: ng0 = nguard*npgs
      integer, parameter :: large_integer = 10000

      logical :: lcc,lfc,ldiag,l_srl_only,lec,lnc

      real :: cfl
      common/courant/cfl
!-----------------------------------------------------------------------

      real, dimension(ne,iu_bnd1,ju_bnd1,ku_bnd1) :: vkn
      real, dimension( 3,iu_bnd1,ju_bnd1,ku_bnd1) :: dum
      real, dimension(iu_bnd1,ju_bnd1,ku_bnd1)    :: cf
      real, dimension(ne,nxb,nyb,nzb)             :: u2
      real :: frx,fry,frz,fr
      real, parameter::  eps = 1.0e-30
      real, dimension(ne) :: varmax,varmin
      logical :: lguard,lprolong,lflux,ledgel,lrestrict,lflag,lfulltree
      integer :: tag_offset,errcode,nprocs
!     integer :: errcode,nprocs,iproc
#ifdef MPI_USED
        call mpi_barrier (MPI_COMM_WORLD, errcode)
#else
        call shmem_barrier_all()
#endif /*MPI_USED*/
       npes = shmem_n_pes()
       nprocs = npes

      nxm2 = nxb+nguard
      nym2 = nyb+nguard*k2d
      nzm2 = nzb+nguard*k3d

      dtmax      = 0.
      dtmin      = 1.e10
      dtl        = 1.e10
      dtlevel(:) = 1.e10


#ifdef NO_PERMANENT_GUARDCELLS
c       write(*,*)'entering amr_1blk_copy_soln:mype=',mype
       call amr_1blk_copy_soln(-1)
c       write(*,*)'leaving amr_1blk_copy_soln:mype=',mype
  
#endif
!        write(*,*)'mk10: amr_timestep: 79 ',mype

#ifdef MPI_USED
        call mpi_barrier (MPI_COMM_WORLD, errcode)
#else
        call shmem_barrier_all()
#endif /*MPI_USED*/
      iopt      = 1
      lcc       = .true.
      lfc       = .false.
#ifdef MPI_USED
      lguard    = .true.
      lprolong  = .false.
      lflux     = .false.
      ledgel    = .false.
      lrestrict = .false.
      lec        =.false.
      lnc        =.false.
      tag_offset= 100
      lfulltree=.false.

c      write(*,*)'mk10: amr_timestep: 80 ',mype
c      write(*,*)'mk10: mype,nprocs,lguard =',mype,nprocs,lguard
c      write(*,*)'mk10: lprolong,lflux,ledgel =',lprolong,lflux,ledgel
c      write(*,*)'mk10: lrestrict,iopt,lcc,lfc =',lrestrict,iopt,lcc,lfc
c      write(*,*)'mk10: lec,lnc,tag_offset =',lec,lnc,tag_offset
c      write(*,*)'entering mpi_amr_comm_setup:mype=',mype
      call mpi_amr_comm_setup(mype,nprocs,
     .                        lguard,lprolong,lflux,ledgel,lrestrict,
     .                        lfulltree,
     .                        iopt,lcc,lfc,lec,lnc,tag_offset)
c       write(*,*)'leaving mpi_amr_comm_setup:mype=',mype

!        write(*,*)'mk10: amr_timestep: 81 ',mype
         call mpi_barrier (MPI_COMM_WORLD, errcode)
!       write(*,*)'mpi_barrier:errcode,mype: ',errcode,mype
#endif /*MPI_USED*/

!       write(*,*)'mk10: amr_timestep: 100:errcode: ',mype


!--------------------------------------------------------------
       do iproc = 0,nprocs-1
       if(mype.eq.iproc)then

! loop over leaf grid blocks
       if(lnblocks.gt.0) then
       do lb=1,lnblocks
        if(nodetype(lb).eq.1) then
#ifdef NO_PERMANENT_GUARDCELLS

#ifdef DEBUG
        if(mype.eq.1.and.lb.eq.2)then
        do k=1,nzb
         do j=1,nyb
          do i=1,nxb
           u2(1:ne,i,j,k)=unk(1:ne,i,j,k,lb) ! time level: n 
          enddo
         enddo
        enddo
        do l = 1, ne
         varmax(l) = maxval( u2(l,:nxb,:nyb,:nzb) )
         varmin(l) = minval( u2(l,:nxb,:nyb,:nzb) )
        enddo
       do l=1,ne
        write(*,*)'a:mype,lb,varmax(',l,'),varmin(',l,') = ',
     &               mype,lb,varmax(l),varmin(l)
       enddo

        endif
#endif 
        iopt       = 1
        nlayers    = nguard
        lcc        =.true.
        lfc        =.false.
        lec        =.false.
        lnc        =.false.
        l_srl_only =.false.
        icoord     = 0
        ldiag      =.true.


c      write(*,*)'mk10: entrying amr_1blk_guardcell:lb,mype=',lb,mype
      call amr_1blk_guardcell(mype,iopt,nlayers,lb,mype,lcc,lfc,
     &                        lec,lnc,l_srl_only,icoord,ldiag)
c      write(*,*)'mk10: leaving amr_1blk_guardcell:lb,mype=',lb,mype

#ifdef DEBUG
      if(mype.eq.1.and.lb.eq.2)then
      do k = kl_bnd1,ku_bnd1
        do j = jl_bnd1,ju_bnd1
         do i = il_bnd1,iu_bnd1
           vkn(1:ne,i,j,k)=unk1(1:ne,i,j,k,1) ! time level: n 
          enddo
         enddo
        enddo
        do l = 1, ne
         varmax(l) = maxval( vkn(l,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1) )
         varmin(l) = minval( vkn(l,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1) )
        enddo
       do l=1,ne
        write(*,*)'b:mype,lb,varmax(',l,'),varmin(',l,') = ',
     &             mype,lb,varmax(l),varmin(l)
       enddo
c      do k = kl_bnd1,ku_bnd1
c        do j = jl_bnd1,ju_bnd1
c         do i = il_bnd1,iu_bnd1
c         if(unk1(1,i,j,k,1).eq.0.0)then
c          write(*,*)'i,j,k=',i,j,k
c          write(*,*)'unk1(1,',i,',',j,',',k,',1)=',unk1(1,i,j,k,1)
c          endif
c  
c         enddo
c       enddo
c      enddo

      endif
#endif
#endif /* NO_PERMANENT_GUARDCELLS */

!----
c      write(*,*)'entering grid_variables:lb,mype=',lb,mype
      call grid_variables(lb,dt)
c      write(*,*)'entering strongfields_wb:lb,mype=',lb,mype
      call strongfields_wb(mype,lb)
c      write(*,*)'leaving strongfields_wb:lb,mype=',lb,mype
   
c      write(*,*)'entering eos_u_to_v_pmn:lb,mype=',lb,mype
      call eos_u_to_v_pmn (unk1(1,1,1,1,1),vkn,cf,lb,mype)
c      write(*,*)'leaving eos_u_to_v_pmn:lb,mype=',lb,mype

!-----new
!
#ifdef INTERIOR_BOUNDARY
!     As the fields are strong within the interior boundary but are not being
!     used to update solutions the sonic speed is set to zero
!     within those cells making up the interior boundary
      lflag             = earth_blocks(lb)
      izones            = idzones(lb) 
      izone(:,:izones)  = p_zone(:,:izones,lb) 
      if(lflag)then
        do iz=1,izones
          i = izone(1,iz)
          j = izone(2,iz)
          k = izone(3,iz)
           vkn(2:4,i,j,k) =0.0
           cf(i,j,k)      =0.0
       enddo
      endif
#endif /*INTERIOR_BOUNDARY*/
!-----new
      frx  = eps
      fry  = eps
      frz  = eps
      do k=1+nguard*k3d,nzm2
       do j=1+nguard*k2d,nym2
        do i=1+nguard    ,nxm2
         frx = max( frx, (abs(vkn(2,i,j,k))+cf(i,j,k))*rdx  )
         fry = max( fry, (abs(vkn(3,i,j,k))+cf(i,j,k))*rdy  )
#if N_DIM < 3
         frz = 0.0
#else
         frz = max( frz, (abs(vkn(4,i,j,k))+cf(i,j,k))*rdz  )
#endif
        enddo
       enddo
      enddo
!--------------------------------------------------------------
      dtold = dt
      fr    = max(frx,fry,frz)
      dtmax = 1.0/fr
      dtl   = cfl*dtmax
      dtl   = min(dtl,1.25*dtold)
ceec
c      dtl   = min (dtl, 4.e-2)
!----


!------------------------------------------
! NAG f95 has a problem with "min(dtlevel(lrefine(lb)), dtl)", so
! I broke it up with xleft and xright -- John.
      xleft                = dtlevel(lrefine(lb))
      xright               = dtl
      dtlevel(lrefine(lb)) = min(xleft, xright)

      endif
      enddo
      endif


      endif  ! loop over processors
c       call mpi_barrier (MPI_COMM_WORLD, errcode)
      enddo  ! loop over processors

!       write(*,*)'mk10: amr_timestep: 101: ',mype

#ifdef MPI_USED
        call mpi_barrier (MPI_COMM_WORLD, errcode)
#else
        call shmem_barrier_all()
#endif /*MPI_USED*/



!------------------------------------------


! find smallest timesteps for each refinement level  across all processors
      do i=1,maxlevels
        call comm_real_min_to_all(dtlevel(i),dtlevel(i))
      enddo

#ifdef MPI_USED
        call mpi_barrier (MPI_COMM_WORLD, errcode)
#else
        call shmem_barrier_all()
#endif /*MPI_USED*/


! ensure that the timestep does not increase as refinement level increases.
       do i=2,maxlevels
        dtlevel(i) = min(dtlevel(i),dtlevel(i-1))
      enddo

! ensure that each timestep is either equal to or a factor 2 larger
! than the timestep at the next higher refinement level.
      dtmin = dtlevel(maxlevels)
      dtmaxl = real(large_integer)*dtmin
      do i=1,maxlevels-1
         xleft = dtlevel(i)
         xright = dtmaxl
        dtlevel(i)=min(xleft, xright)
      enddo
      do i = maxlevels-1,1,-1
        ratio      = dtlevel(i)/dtmin
        iratio     = max((int(ratio)/2)*2,1)
        iratio     = min(iratio,2)
        dtlevel(i) = real(iratio)*dtmin
      enddo



      dtmaxl = 0.
      if(lnblocks.gt.0) then
      do l=1,lnblocks
      if(nodetype(l).eq.1) then
         xleft = dtmaxl
         xright = dtlevel(lrefine(l))
         dtmaxl = max(xleft, xright)
      endif
      enddo
      endif
#ifdef MPI_USED
        call mpi_barrier (MPI_COMM_WORLD, errcode)
#else
        call shmem_barrier_all()
#endif /*MPI_USED*/


! find largest timestep for any leaf node across all processors
      call comm_real_max_to_all(dtmaxl,dtmaxl)
#ifdef MPI_USED
        call mpi_barrier (MPI_COMM_WORLD, errcode)
#else
        call shmem_barrier_all()
#endif /*MPI_USED*/


      dtmax = dtmaxl


#ifdef VAR_DT
      dt    = dtmax
#else
      dt    = dtmin
      dtmax = dtmin
#endif
c      if(mype.eq.0) write(*,*) 'proc ',mype,' dt ',dt
c      write(*,*)'in dtset'
c      write(*,*)'dt,cfl,rdx,rdy,rdz=',dt,cfl,rdx,rdy,rdz
c      pause
#ifdef MPI_USED
        call mpi_barrier (MPI_COMM_WORLD, errcode)
#else
        call shmem_barrier_all()
#endif /*MPI_USED*/


      return
      end
