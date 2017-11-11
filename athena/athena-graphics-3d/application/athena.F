#include "paramesh_preprocessor.fh"
#include "athena_preprocessor.fh"
#include "jet_conditions.fh"
#undef NEVER200

#define CHOMBO_VIS


      program athena

      use physcons
      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace
      use grid
      use strongfields
      use interior_gbc
      use athena_specific_parameters

      use paramesh_interfaces
#ifdef MPI_USED
      use paramesh_mpi_interfaces
#endif /* MPI_USED */


! include file required for shmem library.
#include "amr_shmem.fh"

#ifdef MPI_USED
      include 'mpif.h'
      integer :: tag_offset,max_blks_sent,ierror
      logical :: lguard,lprolong,lflux,leedge,lrestrict
#endif /* MPI_USED */
      integer :: ierr
      save ierr
      real br(4)
      real x1, x2, y1, y2, z1, z2
      real xbr(8), ybr(8), zbr(8)
      real br3(8)

! local amr variables
      integer nprocs,errcode
      integer shmem_my_pe,shmem_n_pes
      integer :: i,j,k,l,lb,ib
      logical lcc,lfc,ldiag,l_srl_only
      logical lec,lnc, lfulltree

      common/evolution/time
      real :: dt,timstp,deltim
      save dt,timstp,deltim
      integer, parameter :: ng0  = nguard*npgs
      real :: cfl
      common/courant/cfl
      real, dimension(ne,il_bnd:iu_bnd,jl_bnd:ju_bnd,
     .                   kl_bnd:ku_bnd    ) :: u2,v2

      real*8  text(20)
      character*80  atext
      integer :: ifilter 
      integer, parameter :: muref = 20
      integer :: nuref,luref
      real    :: ctore,ctode,epsil
      common/refine_parameters/ctore,ctode,epsil,luref(muref),nuref
      integer :: nurefd,lurefd(muref)
      integer :: icheck,icheckg
      save icheck,icheckg
      integer :: iunit1,iunit2
      integer :: iunit3,iunit4,iunit5,iunit6,iunit7
      save nstep_tecplot , nstep_times , iunit4, iunit7

#ifdef IONOSPHERE
! include file defining the ionospheric parameters
#include "ionosphere.fh"
#endif /*IONOSPHERE*/ 
!
!------------------------------------------------------------------------
!     G r a p h i c s 
!------------------------------------------------------------------------
c      real :: xmin,ymin,zmin,xmax,ymax,zmax
      real :: xminb,yminb,zminb,xmaxb,ymaxb,zmaxb
      common/box/xminb,yminb,zminb,xmaxb,ymaxb,zmaxb
      real, dimension(iu_bnd,ju_bnd,ku_bnd)::temp_mask
      real :: vmax(ne,maxblocks)=-1.e30,vmin(ne,maxblocks)=+1.e30
      real :: umax(ne,maxblocks)=-1.e30,umin(ne,maxblocks)=+1.e30
      real :: bxmax(maxblocks),bxmin(maxblocks)
      real :: bymax(maxblocks),bymin(maxblocks)
      real :: bzmax(maxblocks),bzmin(maxblocks)
      real :: gbxmax,gbxmin,gbymax,gbymin,gbzmax,gbzmin
      save gbxmin1,gbxmax1,gbymin1,gbymax1,gbzmin1,gbzmax1
      save gbxmin,gbxmax,gbymin,gbymax,gbzmin,gbzmax
      real :: varmax(ne)=+1.0e+30, varmin(ne)=-1.0e+30
      real :: uarmax(ne)=+1.0e+30, uarmin(ne)=-1.0e+30
 !MG  save varmax,varmin
      integer, dimension(ne) :: imax,jmax,kmax
      integer, dimension(ne) :: imin,jmin,kmin
      real   , dimension(ne) :: fmax,fmin

      logical :: lmark,lmark1,lmark2,lmark3,lmark4
      logical :: lheavy

c#include "interior_boundary.fh"


      real ebox(3,2)
      real xpos1,ypos1,zpos1
      real xpos2,ypos2,zpos2
      save ebox
!
!-----------------------------------------------------------------------
      logical :: lmovie
      integer :: nmovie, movie_si
      integer :: maxstp,minstep
      integer :: nstep_dump,nstep_amr,nstep_graphics,nstep_times,
     &           nstep_tecplot
      integer :: ltecplot
!     logical :: restart
#ifdef INTERACTIVE_GRAPHICS
      integer, parameter :: iunitc=10

#ifdef GRAPHICS
      real*4 wxmn, wymn, wxmx, wymx
      real,dimension(ne,il_bnd:iu_bnd,
     &                  jl_bnd:ju_bnd,
     &                  kl_bnd:ku_bnd,maxblocks)::vnk
      save vnk
      real, dimension(ne,il_bnd:iu_bnd,jl_bnd:ju_bnd,
     .                   kl_bnd:ku_bnd    ) :: ud
#endif /*GRAPHICS*/
!
!-----------------------------------------------------------------------
      namelist /DRIVERPARMS/
     &     nstep_dump, nstep_amr, nstep_graphics,nstep_print,
     &     ltecplot,nstep_tecplot,nstep_times,cfl,
     &     ctore,ctode,epsil,nuref,luref,
     &     iunit4, iunit7
!    &     iunit1, iunit2, iunit3, iunit5, iunit6 
#else
#endif /*INTERACTIVE*/
#ifdef SGI_SHMEM
      save maxstp,minstp
      save nstep_dump,nstep_amr,nstep_graphics,nstep_print,ltecplot,
     & nstep_tecplot
!     save restart
#endif /*SGI_SHMEM*/
!
!------------------------------------------------------------------------
!	L o g i c a l s
!------------------------------------------------------------------------
!
      logical :: lgraph, restart,lrefine_again
      save lgraph, restart,lrefine_again
!      
!------------------------------------------------------------------------
!	I n t e g e r s
!------------------------------------------------------------------------
!
      integer :: istep,istep_max,loop_start_index,no_of_steps,nchombo



      real :: jet_den, jet_p, jet_vx,jet_vy,jet_vz, jet_radius, jet_beta,
     .     ambient_den, ambient_vy, ambient_p
      common/jet_state/jet_den, jet_p, jet_vx,jet_vy,jet_vz,jet_radius, jet_beta,
     .     ambient_den, ambient_vy, ambient_p

      

      real :: delx0, delxn, xmn1, xmn2, xmx1, xmx2
      integer :: ijet_radius


c
c     cooling data structures
c

      real, dimension(50) :: Ttab, Ltab
      common/coolstuff/Ttab,Ltab

C
C     wrong cooling function: should use zero field
C



      Ltab(1:50)=(/-22.97,-22.87,-22.77,-22.64,-22.52,-22.42,-22.31,-22.22,-22.14,-22.07,-22.02,-21.97,-21.91,-21.86,-21.81,
     .     -21.75,-21.69,-21.63,-21.58,-21.53,-21.49,-21.45,-21.41,-21.38,-21.37,-21.36,-21.37,-21.4,-21.45,-21.52,-21.61,
     .     -21.69,-21.75,-21.78,-21.81,-21.85,-21.91,-21.96,-21.98,-21.98,-21.98,-21.97,-21.95,-21.95,-21.99,-22.1,-22.23,
     .     -22.34,-22.41,-22.44/)

      do j=1,50
         Ltab(j)=Ltab(j)+22.0 ! removing 10^-22 factor
      enddo

      Ttab(1:50)=(/4.04,4.08,4.12,4.17,4.21,4.25,4.3,4.34,4.38,4.44,4.48,4.53,4.58,4.62,4.67,4.71,4.77,4.82,4.88,4.93,4.98,
     .     5.04,5.09,5.14,5.2,5.25,5.3,5.36,5.41,5.46,5.51,5.57,5.62,5.67,5.73,5.78,5.83,5.88,5.94,5.99,6.04,6.09,6.14,6.2,
     .     6.25,6.3,6.35,6.4,6.45,6.5/)



!     http://www.journals.uchicago.edu/ApJ/journal/issues/ApJ/v560n2/53692/53692.html
      jet_beta=0.6
      jet_radius=JET_XMAX/24.!  will be adjusted to fit cell spacings later!



c      jet_p=1.0
c      jet_den=1.0
c      jet_vy=10.0


      jet_p=JET_P
      jet_den=JET_DEN
      jet_vy=JET_VY

      


!
!------------------------------------------------------------------------
!	F o r m a t s
!------------------------------------------------------------------------
!
  111 format(40a4)
    1 format(40a4)
!------------------------------------------------------------------------
!

! make initialization call for amr package
#ifdef SGI_SHMEM
            write(*,*)'calling comm_start: 95'
        call comm_start(maxprocs,nprocs,mype)
         write(*,*)'leaving comm_start'
         write(*,*)'nprocs =',nprocs
#endif


        call amr_initialize



czzc
c        interp_mask_unk(:)=2
c        interp_mask_facex(:)=2
c        interp_mask_facey(:)=2
c        interp_mask_facez(:)=2


        mype   = shmem_my_pe()
        nprocs = shmem_n_pes()
        if(mype.eq.0) write(*,*)'exiting amr_ini: 945',mype

        call MPI_COMM_SIZE(MPI_COMM_WORLD, numb_processors, ierrcode)
        call MPI_COMM_RANK(MPI_COMM_WORLD, id, ierrcode)
        if(mype.eq.0)
     &  write(*,*)'mpi-call: numb_processors, id ',
     &          numb_processors, id

!       call amr_close
!       stop

#ifdef IONOSPHERE
       phi3d(:,:,:,:) = 0.
#endif /*IONOSPHERE*/
!---------------------------------------------------------------
!
!-----Default settings
!

      time         = 0.
      nstep_total  = 7000
      nstep_dump   = 100
!
! I/O Unit numbers
      iunit1 = 10
      iunit2 = 11
      iunit3 = 49
      iunit4 = 50
      iunit5 = 51
      iunit6 = 52
      iunit7 = 53
      loop_start_index = 0
!      nstep_amr = 10
czzc
      nstep_amr = 1

      nchombo   = 0
!
! BLOCK 1
!       print*,'! BLOCK 1'
        restart = .false.


      if(mype==0) then
        minstp = 1
        maxstp = 20
        call read_control('athena.control', maxstp,
     &       lgraph,nstep_graphics,restart,
     &       nstep_dump,nstep_print,
     &       ltecplot,nstep_tecplot,
     &       rmin,pmin,smin,
     &       cfl,dt,deltim,
     &       ctore,ctode,epsil,nuref,
     &       luref, muref,
     &       nmovie,
     &       movie_si,
     &       lmovie,
     &       arfact2)
       timstp = time+deltim
       write(*,*)'lgraph,nstep_graphics,restart,
     &              nstep_dump,nstep_print= ',
     &             lgraph,nstep_graphics,
     &             restart,nstep_dump,nstep_print
        write(*,*)'ltecplot,nstep_print=',ltecplot,nstep_print
        write(*,*)'rmin,pmin,smin= ',rmin,pmin,smin
        write(*,*)'cfl = ',cfl
        write(*,*)'ctore,ctode,nuref=',ctore,ctode,epsil,nuref
        if(nuref.gt.muref) then
          write(*,*)' please increase muref, needed:',nuref
          stop
        endif
        do ir=1,nuref        
         write(*,*)'luref(',ir,')=',luref(ir)
        enddo
        write(*,*)'nmovie       =',nmovie
        write(*,*)'movie_si     =',movie_si
        write(*,*)'lmovie       =',lmovie
        write(*,*)'earth_radius =',arfact2


      endif           ! mype==0



#ifndef GRAPHICS
CXXC      lgraph = .false.
#endif /*GRAPHICS*/
#ifdef SHMEM_USED

        call shmem_barrier_all()
        write(*,*) '-1:Running on ',nprocs,' processors,mype=',mype
!       call flush(6)
        call shmem_barrier_all()

      if (mype.ne.0) then
         call shmem_integer_get(nstep_print   ,nstep_print   ,1,0)
         call shmem_integer_get(maxstp        ,maxstp        ,1,0)
         call shmem_integer_get(minstp        ,minstp        ,1,0)
c         call shmem_integer_get(nprocs        ,nprocs        ,1,0)
         call shmem_integer_get(nstep_amr     ,nstep_amr     ,1,0)
         call shmem_integer_get(nstep_graphics,nstep_graphics,1,0)
         call shmem_integer_get(nstep_dump    ,nstep_dump    ,1,0)
         call shmem_integer_get(nstep_tecplot ,nstep_tecplot ,1,0)
         call shmem_logical_get(restart       ,restart       ,1,0)
         call shmem_logical_get(lgraph        ,lgraph       ,1,0)
         call shmem_integer_get(ltecplot      ,ltecplot      ,1,0)
         call shmem_real_get(arfact2      ,arfact2     ,1,0)
         call shmem_real_get(rmin     ,rmin     ,1,0)
         call shmem_real_get(pmin     ,pmin     ,1,0)
         call shmem_real_get(smin     ,smin     ,1,0)
         call shmem_real_get(cfl      ,cfl      ,1,0)
         call shmem_real_get(dt       ,dt       ,1,0)
         call shmem_real_get(deltim   ,deltim    ,1,0)
         call shmem_real_get(timstp   ,timstp    ,1,0)
          call shmem_real_get(ctore    ,ctore      ,1,0)
         call shmem_real_get(ctode    ,ctode      ,1,0)
         call shmem_real_get(epsil    ,epsil      ,1,0)
         call shmem_integer_get(nuref ,nuref      ,1,0)
         call shmem_integer_get(luref ,luref      ,20,0)
      end if
        write(*,*) '0:Running on ',nprocs,' processors,mype=',mype
!       call flush(6)
        call shmem_barrier_all()
#endif /*SHMEM_USED*/
#ifdef MPI_USED
      call MPI_Barrier   (MPI_COMM_WORLD,ierror)
        call mpi_bcast(maxstp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call mpi_bcast(minstp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call mpi_bcast(lgraph,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
        call mpi_bcast(nstep_graphics,1,
     &                 MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call mpi_bcast(restart,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
        call mpi_bcast(nstep_dump,1,
     &                 MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call mpi_bcast(nstep_print,1,
     &                 MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call mpi_bcast(nstep_amr,1,
     &                 MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call mpi_bcast(ltecplot,1,
     &                 MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call mpi_bcast(nstep_tecplot,1,
     &                 MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call mpi_bcast(rmin,1,MPI_DOUBLE_PRECISION,0,
     &                     MPI_COMM_WORLD,ierr)
        call mpi_bcast(pmin,1,MPI_DOUBLE_PRECISION,0,     
     &                     MPI_COMM_WORLD,ierr)
        call mpi_bcast(smin,1,MPI_DOUBLE_PRECISION,0,
     &                     MPI_COMM_WORLD,ierr)
        call mpi_bcast(cfl,1,MPI_DOUBLE_PRECISION,0,
     &                     MPI_COMM_WORLD,ierr)
        call mpi_bcast(dt,1,MPI_DOUBLE_PRECISION,0,
     &                     MPI_COMM_WORLD,ierr)
         call mpi_bcast(deltim,1,MPI_DOUBLE_PRECISION,0,
     &                     MPI_COMM_WORLD,ierr)
        call mpi_bcast(timstp,1,MPI_DOUBLE_PRECISION,0,
     &                     MPI_COMM_WORLD,ierr)
        call mpi_bcast(arfact2,1,MPI_DOUBLE_PRECISION,0,
     &                     MPI_COMM_WORLD,ierr)
        call mpi_bcast(ctore,1,MPI_DOUBLE_PRECISION,0,
     &                     MPI_COMM_WORLD,ierr)
        call mpi_bcast(ctode,1,MPI_DOUBLE_PRECISION,0,
     &                     MPI_COMM_WORLD,ierr)
        call mpi_bcast(epsil,1,MPI_DOUBLE_PRECISION,0,
     &                     MPI_COMM_WORLD,ierr)
        call mpi_bcast(nuref,1,MPI_INTEGER,0,
     &                     MPI_COMM_WORLD,ierr)
        call mpi_bcast(luref,nuref,
     &                 MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call mpi_bcast(nmovie,1,
     &                 MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call mpi_bcast(movie_si,1,
     &                 MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call mpi_bcast(lmovie,1,
     &                 MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)

      call MPI_Barrier   (MPI_COMM_WORLD,ierror)
 !MG  write(*,*) 'pe ',mype,' minstp ',minstp
 !MG  write(*,*) 'pe ',mype,' maxstp ',maxstp
 !MG  write(*,*) 'pe ',mype,' rmin ',rmin
 !MG  write(*,*) 'pe ',mype,' nuref ',nuref
 !MG  write(*,*) 'pe ',mype,' luref ',luref
#endif /* MPI_USED */
!
! BLOCK 2
!       print*,'! BLOCK 2'

c        if(mype==0) write(*,*) 'Running on ',nprocs,' processors'
c        write(*,*) '1:Running on ',nprocs,' processors,mype=',mype
!       call flush(6)
        call shmem_barrier_all()

!---------------------------------------------------------------
!
! BLOCK 3  INITIAL GRID REFINEMENT
!       print*,'! BLOCK 3'


! set a limit on the refinement level
#if N_DIM < 3
      lrefine_initial = 7
      lrefine_initial = 4        ! 4^3 cells for the Earth
      lrefine_initial = 4        ! 4^3 cells for the Earth
#else /* N_DIM == 3*/
      lrefine_initial = 4        ! 4^3 cells for the Earth
      lrefine_initial = 4        ! 4^3 cells for the Earth
c      lrefine_initial = 4        ! 4^3 cells for the Earth
#endif       
        lrefine_max = 7
        lrefine_min = 2
        
c        lrefine_initial = lrefine_min


#if N_DIM == 3
        roa =  3.50*rearth
#else
        roa = 15.0*rearth
#endif

!
! set up initial physical state or read in checkpoint file
!---------------------------------------------------------


!
! Start a new calculation or restart from a previous checkpoint file ?
!------
      if (restart) then
!------

!
! If a restart read the existing checkpoint file
#ifdef MAGNETOSPHERE
      call initialmesh(mype,xminb,yminb,zminb,xmaxb,ymaxb,zmaxb)
c      call strongfield_initial(mype)
      bflags(1,:) = 0
#endif

      call amr_checkpoint_re_a(dt,time,loop_start_index,nchombo,iunit4)
      write(*,601) ebox_xmin, ebox_xmax, mype
      write(*,602) ebox_ymin, ebox_ymax, mype
      write(*,603) ebox_zmin, ebox_zmax, mype
 601  format(2x,"check: ebox-x ",f10.5,3x,f10.5, 5x,i2)
 602  format(2x,"check: ebox-y ",f10.5,3x,f10.5, 5x,i2)
 603  format(2x,"check: ebox-z ",f10.5,3x,f10.5, 5x,i2)
      call shmem_barrier_all()

        go to 1000
!------
      else  ! not restart
!------
      
! set up initial grid state.
c        write(*,*) '2:Running on ',nprocs,' processors,mype=',mype
!       call flush(6)
        call shmem_barrier_all()

! set up step function on a single block covering the whole square domain
        lnblocks = 0
        if(mype==0) then
                lnblocks = 1
                coord(:,1) = 0.
		bnd_box(1,:,1) =-4.
		bnd_box(2,:,1) = 4.
                bsize(:,1) = bnd_box(2,:,1) - bnd_box(1,:,1)
                nodetype(1) = 1
                lrefine(1) = 1

                neigh(1,:,1) = 1                  ! initial block is its own
                neigh(2,:,1) = 0                  ! neighbor ie periodic bc's

                refine(1)=.true.
                refine(1)=.false.
        endif

!===========================>>>> PETER CUT HERE <<<<<<<==========


#ifdef EXPLOSION
! Establish solution on the initial grid.
	 call initialmesh(mype,xminb,yminb,zminb,xmaxb,ymaxb,zmaxb)

! Now cycle over blocks, refining all existing leaf blocks.
#ifdef INITIAL_REFINE_MESH

        do loop_count=1,4
c        do loop_count=1,1

         refine(1:lnblocks) = .true.

         call shmem_barrier_all()

! refine grid and apply morton reordering to grid blocks if necessary
         write(*,*)'explosion:entering amr_refine_derefine I'
         call amr_refine_derefine
         write(*,*)'explosion:leaving amr_refine_derefine I'

! prolong solution to any new leaf blocks if necessary
         iopt    = 1
         nlayers = nguard
c         call amr_prolong(mype,iopt,nlayers)

        enddo

          write(*,*) '*****************************************'
!          do i = 1, lnblocks
!             write(*,*) i, nodetype(i),(coord(j,i),j=1,2),
!     &   (bnd_box(j,1,i),j=1,2),(bnd_box(j,2,i),j=1,2), mype
!          end do
          write(*,*) '*****************************************'
         write(*,*)'explosion:entering amr_initial_bomb I'
	call amr_initial_bomb(mype,xminb,yminb,zminb,xmaxb,ymaxb,zmaxb)
         write(*,*)'explosion:leaving amr_initial_bomb I'

        call shmem_barrier_all()
c        ldiag = .true.
c        call amr_gsurrounding_blks(mype,ldiag)
#else
	call amr_initial_bomb(mype,xminb,yminb,zminb,xmaxb,ymaxb,zmaxb)

        call shmem_barrier_all()
c        ldiag = .true.
c        call amr_gsurrounding_blks(mype,ldiag)

#endif /*INITIAL_REFINE_MESH*/

#endif /*EXPLOSION*/
!===========================>>>> PETER CUT HERE <<<<<<<==========

!===========================>>>> PETER CUT HERE <<<<<<<==========
#ifdef JET



        write(*,*)'SET UP JET REFINEMENT' 
        write(*,*)'------------------- 707 ------------',mype
      

!
! Establish solution on the initial grid.
!-----------------------------------------------------------------------
!
!-------Create initial mesh
!
      write(*,*)mype,'entering initialmesh'
      call initialmesh(mype,xminb,yminb,zminb,xmaxb,ymaxb,zmaxb)
!-----------------------------------------------------------------------
cyyc      write(*,*)'1:entering initialstate'
cyyc      call initialstate(mype,xmin,ymin,zmin,xmax,ymax,zmax)
cyyc      write(*,*)'leaving initialstate'

      call mpi_amr_global_domain_limits

      boundary_index(1:6 ) = -21
!      boundary_index(1   ) = -22
      boundary_index(1   ) = -38
      boundary_index(3   ) = -40
      boundary_index(4   ) = -50

! x boundaries 
      boundary_box(1,2:3,1:2) = -1.e30
      boundary_box(2,2:3,1:2) =  1.e30
      boundary_box(1,1,1    ) = -1.e30
      boundary_box(2,1,1    ) = grid_xmin
      boundary_box(1,1,2    ) = grid_xmax
      boundary_box(2,1,2    ) =  1.e30
! y boundaries 
      if(ndim.ge.2) then
      boundary_box(1,1,3:4) = -1.e30
      boundary_box(2,1,3:4) =  1.e30
      boundary_box(1,1,3:4) =  grid_xmin
      boundary_box(2,1,3:4) =  grid_xmax
      boundary_box(1,3,3:4) = -1.e30
      boundary_box(2,3,3:4) =  1.e30
      boundary_box(1,2,3  ) = -1.e30
      boundary_box(2,2,3  ) =  grid_ymin
      boundary_box(1,2,4  ) =  grid_ymax
      boundary_box(2,2,4  ) =  1.e30
      endif
! z boundaries 
      if(ndim.eq.3) then
c      boundary_box(1,1:2,5:6) = -1.e30
c      boundary_box(2,1:2,5:6) =  1.e30
      boundary_box(1,1,5:6) =  grid_xmin
      boundary_box(2,1,5:6) =  grid_xmax
      boundary_box(1,2,5:6) =  grid_ymin
      boundary_box(2,2,5:6) =  grid_ymax
       
      boundary_box(1,3,5  ) = -1.e30
      boundary_box(2,3,5  ) = grid_zmin
      boundary_box(1,3,6  ) = grid_zmax
      boundary_box(2,3,6  ) =  1.e30
      endif

      call shmem_barrier_all()

! Now cycle over blocks, refining all existing leaf blocks.
!
!-----store information concerning refinement
        lurefd(:nuref)= luref(:nuref)
        nurefd        = nuref
        luref(1)      = 1
        nuref         = 1

! Now cycle over blocks adjusting refinement of initial setup as required

         lrefine_again = .true.
         refine(:)=.true. ! refine master block
         istep         = 0
         istep_max     = lrefine_max
         
         write(*,*)'entering do while loop,lnblocks=',lnblocks
         

         do while(lrefine_again.and.(istep.lt.istep_max))
            
        write(*,*) 'mype=',mype,'isteppp=',istep,istep_max,'lnblocks=',lnblocks

! refine grid and apply morton reordering to grid blocks if necessary

c                              refine(:)=.true.  ! temp test mod - uniform
c       if(istep.lt.istep_max-1) refine(:)=.true.  ! temp test mod - nonuniform
c       if(istep.lt.istep_max-2) refine(:)=.true.  ! temp test mod - nonuniform

            write(*,*)'entering amr_refine_derefine'
            call amr_refine_derefine
            write(*,*)'leaving amr_refine_derefine'

            if(lnblocks.gt.0) then
               refine  (1:lnblocks) = .false.
               derefine(1:lnblocks) = .false.
            endif
      
            call shmem_barrier_all()

cyyc      write(*,*)'2:entering initialstate'
cyyc      call initialstate(mype,xmin,ymin,zmin,xmax,ymax,zmax)

cyyc      where (lrefine.eq.lrefine_max) refine = .FALSE.
cyyc      lrefine_again=.false.



            call shmem_barrier_all()


#if N_DIM == 3

      do ib = 1,lnblocks

         bflags(1,ib) = 0   

         if (nodetype(ib).eq.1) then

            xmx1=coord(1,ib) -JET_X_0 + bsize(1,ib)/2.
            xmn1=coord(1,ib) -JET_X_0 - bsize(1,ib)/2.
            if(xmx1*xmn1.le.0.) then
               dx1=0.
            elseif (coord(1,ib)-JET_X_0.gt.0.) then
               dx1 = xmn1
            else
               dx1 = xmx1
            end if

            xmx2=coord(1,ib) +JET_X_0 + bsize(1,ib)/2.
            xmn2=coord(1,ib) +JET_X_0 - bsize(1,ib)/2.
            if(xmx2*xmn2.le.0.) then
               dx2=0.
            elseif (coord(1,ib)+JET_X_0.gt.0.) then
               dx2 = xmn2
            else
               dx2 = xmx2
            end if

            ymx=coord(2,ib) + bsize(2,ib)/2.
            ymn=coord(2,ib) - bsize(2,ib)/2.
            if(ymx*ymn.le.0.) then
               dy=0.
            elseif (coord(2,ib).gt.0.) then
               dy = ymn
            else
               dy = ymx
            end if

            zmx=coord(3,ib) + bsize(3,ib)/2.
            zmn=coord(3,ib) - bsize(3,ib)/2.
            if(zmx*zmn.le.0.) then
               dz=0.
            elseif (coord(3,ib).gt.0.) then
               dz = zmn
            else
               dz = zmx
            end if
            
            dx = abs(dx)
            dy = abs(dy)
            dz = abs(dz)

CYYC              force refinement around jet inflow aperture
            r=4.0*jet_radius


                  if ((dx1.le.r.or.dx2.le.r).and.dy.le.r.and.dz.le.r) then
                     refine  (ib) = .TRUE.
                     derefine(ib) = .FALSE.
! Prevent blocks in immediate neighborhood of the earth
! from derefining at any time. To do this we set bflags(1,.)=1
! for these blocks, are test this flag in test_refinement.
                     bflags(1,ib) = 1
c                     write(*,*)'bflags(1,ib) = ',bflags(1,ib) 
                  end if

               end if
               
            end do

#else /*NOT NDIM3*/
            do ib = 1,lnblocks
               
               bflags(1,ib) = 0 ! initialise to non-earth block

               if (nodetype(ib).eq.1) then
                  
                  xmx=coord(1,ib) + bsize(1,ib)/2.
                  xmn=coord(1,ib) - bsize(1,ib)/2.
                  if(xmx*xmn.le.0.) then
                     dx=0.
                  elseif (coord(1,ib).gt.0.) then
                     dx = xmn
                  else
                     dx = xmx
                  end if

                  ymx=coord(2,ib) + bsize(2,ib)/2.
                  ymn=coord(2,ib) - bsize(2,ib)/2.
                  if(ymx*ymn.le.0.) then
                     dy=0.
                  elseif (coord(2,ib).gt.0.) then
                     dy = ymn
                  else
                     dy = ymx
                  end if
                  
                  dx = abs(dx)
                  dy = abs(dy)
                  

                  r  = 2.1*jet_radius
                  
c                  write(*,*) ib, 'bflags', r, dx,dy
                  if (dx.le.r.and.dy.le.r) then
                     refine(ib)   = .TRUE.
                     derefine(ib) = .FALSE.
! Prevent blocks in immediate neighborhood of the earth
! from derefining at any time. To do this we set bflags(1,.)=1
! for these blocks, are test this flag in test_refinement.
                     bflags(1,ib) = 1
c                     write(*,*)ib,'bflags(1,ib) = ',bflags(1,ib) 
                   
                  end if

c                  r = 10.0*rearth
c                  r = 4.*r
                  
               end if
            
            end do
#endif /*NDIM3*/

! Final safety check to prevent refinement below lrefine_max
            
            do ib=1,lnblocks
               if (lrefine(ib).ge.lrefine_max) refine(ib) = .FALSE.
            enddo
            
            icheckg=0
            icheck=0
            do ib=1,lnblocks
               if( (refine(ib)).or.(derefine(ib)) ) icheck=1
            enddo

            call COMM_INT_MAX_TO_ALL(icheckg,icheck)

            if(icheckg.eq.1) lrefine_again=.true.
            istep=istep+1
            if (mype==0) print *,'main:  ISTEP = ',ISTEP

! prolong solution to any new leaf blocks if necessary
         iopt    = 1
         nlayers = nguard
c         call amr_prolong(mype,iopt,nlayers)

         if(istep.lt.istep_max-1) lrefine_again = .true. ! temp - nonuniform
c      lrefine_again = .true. ! temp - uniform

      enddo !do while


c      stop
      

CYYC
C      have set up refinement over inflow, now specify solution on leaf nodes

c      do ib = 1,lnblocks
c         if (nodetype(ib).eq.1) then
            write(*,*) mype,lnblocks,'...3:entering initialstate'
            call initialstate(mype,xmin,ymin,zmin,xmax,ymax,zmax)
            write(*,*) mype,lnblocks,'...3:entering initialstate'
c         endif
c      enddo




        luref(:nurefd)= lurefd(:nurefd)
        nuref         = nurefd
! prolong solution to any new leaf blocks if necessary
        iopt    = 1
        nlayers = nguard
!pmn        call amr_prolong(mype,iopt,nlayers)

! exchange guardcell information
      lcc = .true.
      lfc = .true.
      tag_offset = 100

      lguard    = .true.
      lprolong  = .false.
      lflux     = .false.
      leedge    = .false.
      lrestrict = .false.
      lec        =.false.
      lnc        =.false.
      lfulltree=.false.

      call mpi_amr_comm_setup(mype,nprocs,
     .                        lguard,lprolong,lflux,leedge,lrestrict,
     .                        lfulltree,
     .                        iopt,lcc,lfc,lec,lnc,tag_offset)

      write(*,*) mype,lnblocks,'...3:entering guardcell'
!       call mpi_amr_guardcell(mype,iopt,nlayers,nprocs)
        call amr_guardcell(mype,iopt,nlayers)

      write(*,*) mype,lnblocks,'...3:leaving guardcell'

#ifdef OLD


! Establish solution on the initial grid.
	 call initialmesh(mype,xminb,yminb,zminb,xmaxb,ymaxb,zmaxb)


      delx0=bsize(1,1)/real(nxb-gc_off_x)
      delxn=delx0/2.0**real(lrefine_max-1)  
      ijet_radius=int(jet_radius/delxn+0.5)
      jet_radius=(ijet_radius+0.01)*delxn !  round jet radius to integer number of finest grid cell widths, plus offset to avoid rounding discrepencies
      write(*,*) 'jet_radius', jet_radius, 'xmax,xmin', xmax,xmin



! Now cycle over blocks, refining all existing leaf blocks.

       do loop_count=1,4
c        do loop_count=1,lrefine_initial

         refine(1:lnblocks) = .true.

         call shmem_barrier_all()

! refine grid and apply morton reordering to grid blocks if necessary
         write(*,*)'jet:entering amr_refine_derefine I'
         call amr_refine_derefine
         write(*,*)'jet:leaving amr_refine_derefine I'

! prolong solution to any new leaf blocks if necessary
         iopt    = 1
         nlayers = nguard
c         call amr_prolong(mype,iopt,nlayers)

        enddo

          write(*,*) '*****************************************'
!          do i = 1, lnblocks
!             write(*,*) i, nodetype(i),(coord(j,i),j=1,2),
!     &   (bnd_box(j,1,i),j=1,2),(bnd_box(j,2,i),j=1,2), mype
!          end do
          write(*,*) '*****************************************'
         write(*,*)'jet:entering amr_initial_jet I'
	call amr_initial_jet(mype,xminb,yminb,zminb,xmaxb,ymaxb,zmaxb)
         write(*,*)'jet:leaving amr_initial_jet I'

        call shmem_barrier_all()
c        ldiag = .true.
c        call amr_gsurrounding_blks(mype,ldiag)

#endif /* OLD */



#endif /*JET*/

!===========================>>>> PETER CUT HERE <<<<<<<==========

!===========================>>>> PETER CUT HERE <<<<<<<==========
#ifdef CIRCULATION
! Establish solution on the initial grid.
	 call initialmesh(mype,xminb,yminb,zminb,xmaxb,ymaxb,zmaxb)

! Now cycle over blocks, refining all existing leaf blocks.
#ifdef INITIAL_REFINE_MESH

c       do loop_count=1,4
        do loop_count=1,4

         refine(1:lnblocks) = .true.

         call shmem_barrier_all()

! refine grid and apply morton reordering to grid blocks if necessary
         write(*,*)'circulation:entering amr_refine_derefine I'
         call amr_refine_derefine
         write(*,*)'circulation:leaving amr_refine_derefine I'

! prolong solution to any new leaf blocks if necessary
         iopt    = 1
         nlayers = nguard
c         call amr_prolong(mype,iopt,nlayers)

        enddo

          write(*,*) '*****************************************'
!          do i = 1, lnblocks
!             write(*,*) i, nodetype(i),(coord(j,i),j=1,2),
!     &   (bnd_box(j,1,i),j=1,2),(bnd_box(j,2,i),j=1,2), mype
!          end do
          write(*,*) '*****************************************'
         write(*,*)'circulation:entering amr_initial_circulation I'
	call amr_initial_circulation(mype,xminb,yminb,zminb,xmaxb,ymaxb,zmaxb)
         write(*,*)'circulation:leaving amr_initial_circulation I'

        call shmem_barrier_all()
c        ldiag = .true.
c        call amr_gsurrounding_blks(mype,ldiag)
#else
	call amr_initial_circulation(mype,xminb,yminb,zminb,xmaxb,ymaxb,zmaxb)

        call shmem_barrier_all()
c        ldiag = .true.
c        call amr_gsurrounding_blks(mype,ldiag)

#endif /*INITIAL_REFINE_MESH*/

#endif /*CIRCULATION*/
!===========================>>>> PETER CUT HERE <<<<<<<==========

!===========================>>>> PETER CUT HERE <<<<<<<==========
#ifdef PULSAR
! Establish solution on the initial grid.
	 call initialmesh(mype,xminb,yminb,zminb,xmaxb,ymaxb,zmaxb)

! Now cycle over blocks, refining all existing leaf blocks.
#ifdef INITIAL_REFINE_MESH

c       do loop_count=1,4
        do loop_count=1,4

         refine(1:lnblocks) = .true.

         call shmem_barrier_all()

! refine grid and apply morton reordering to grid blocks if necessary
         write(*,*)'pulsar:entering amr_refine_derefine I'
         call amr_refine_derefine
         write(*,*)'pulsar:leaving amr_refine_derefine I'

! prolong solution to any new leaf blocks if necessary
         iopt    = 1
         nlayers = nguard
c         call amr_prolong(mype,iopt,nlayers)

        enddo

          write(*,*) '*****************************************'
!          do i = 1, lnblocks
!             write(*,*) i, nodetype(i),(coord(j,i),j=1,2),
!     &   (bnd_box(j,1,i),j=1,2),(bnd_box(j,2,i),j=1,2), mype
!          end do
          write(*,*) '*****************************************'
         write(*,*)'pulsar:entering amr_initial_pulsar I'
	call amr_initial_pulsar(mype,xminb,yminb,zminb,xmaxb,ymaxb,zmaxb)
         write(*,*)'pulsar:leaving amr_initial_pulsar I'

        call shmem_barrier_all()
c        ldiag = .true.
c        call amr_gsurrounding_blks(mype,ldiag)
#else
	call amr_initial_pulsar(mype,xminb,yminb,zminb,xmaxb,ymaxb,zmaxb)

        call shmem_barrier_all()
c        ldiag = .true.
c        call amr_gsurrounding_blks(mype,ldiag)

#endif /*INITIAL_REFINE_MESH*/

#endif /*PULSAR*/
!===========================>>>> PETER CUT HERE <<<<<<<==========



c        write(*,*) '3:Running on ',nprocs,' processors,mype=',mype
!       call flush(6)
        call shmem_barrier_all()

#ifdef MAGNETOSPHERE
#ifdef REFINE_UNIFORMILY
        if(mype.eq.0)
     &  write(*,*)'------------------- 704 ------------',mype
        if(mype.eq.0)
     &  write(*,*)'REFINE_UNIFORMILY:mype=',mype  
	call initialmesh(mype,xminb,yminb,zminb,xmaxb,ymaxb,zmaxb)

c        call strongfield_initial(mype)

#ifdef INITIAL_REFINE_MESH
        if(mype.eq.0) 
     &  write(*,*)'1:REFINE_UNIFORMILY:mype=',mype 
        
        do loop_count=1,3
!print*,'loop_count=',loop_count
            refine(1:lnblocks) = .true.
            call shmem_barrier_all()
!print*,'called shmem_barrier_all()'
            call amr_refine_derefine
!print*,'called amr_refine_derefine'
! prolong solution to any new leaf blocks if necessary
            iopt    = 1
            nlayers = nguard
!call amr_prolong(mype,iopt,nlayers)

        enddo
!-----------------------------------------------------
!           write(*,*)'# of blocks is = ', lnblocks
!        do i =1, lnblocks
!           if(nodetype(i).eq.1) then
!        write(*,52) i, ((bnd_box(j,m,i),j=1,2),m=1,3) 
!           endif
!        end do
!50      format(2x,i3,5x,2(f13.5,2x),5x,2(f13.5,2x),2(f13.5,2x))
52      format(2x,i3,5x,2(f13.5,2x),5x,2(f13.5,2x),2(f13.5,2x))

!-----------------------------------------------------
         do j = 1, 0
            x1=-0.5*300.*rearth
            x1=-0.5*100.*rearth
            x1= 0.5*xminb
            write(*,*)'xminb=',xminb

            refine(1:lnblocks) = .false.
            derefine(1:lnblocks) = .false.
            do i = 1, lnblocks
             if(nodetype(i).eq.1) then
              if(coord(1,i).gt.x1) refine(i)=.true.
             endif
            end do

            call shmem_barrier_all()
            call amr_refine_derefine
         end do

            call shmem_barrier_all()
         eps=1.0e-05
         do j = 1,1
!write(*,*)'105: j, lnblocks is = ',j,  lnblocks
            refine(1:lnblocks) = .false.
            derefine(1:lnblocks) = .false.
            do i = 1, lnblocks
            if(nodetype(i).eq.1) then
#if N_DIM == 2
              x1=coord(1,i)-bsize(1,i)/2; x2=coord(1,i)+bsize(1,i)/2;
              y1=coord(2,i)-bsize(2,i)/2; y2=coord(2,i)+bsize(2,i)/2;
              br(1)=sqrt(x1*x1+y1*y1); br(2)=sqrt(x2*x2+y1*y1);
              br(3)=sqrt(x2*x2+y2*y2); br(4)=sqrt(x1*x1+y2*y2);
              brmin=minval(br(1:4)); brmax=maxval(br(1:4))
#endif

#if N_DIM == 3
              x1=coord(1,i)-bsize(1,i)/2; x2=coord(1,i)+bsize(1,i)/2;
              y1=coord(2,i)-bsize(2,i)/2; y2=coord(2,i)+bsize(2,i)/2;
              z1=coord(3,i)-bsize(3,i)/2; z2=coord(3,i)+bsize(3,i)/2;

              xbr(1)=x1;  ybr(1)=y1;  zbr(1)=z1;
              xbr(2)=x2;  ybr(2)=y1;  zbr(2)=z1;
              xbr(3)=x2;  ybr(3)=y2;  zbr(3)=z1;
              xbr(4)=x1;  ybr(4)=y2;  zbr(4)=z1;
              xbr(5)=x1;  ybr(5)=y1;  zbr(5)=z2;
              xbr(6)=x2;  ybr(6)=y1;  zbr(6)=z2;
              xbr(7)=x2;  ybr(7)=y2;  zbr(7)=z2;
              xbr(8)=x1;  ybr(8)=y2;  zbr(8)=z2;
              do ibr = 1, 8
                 br3(ibr)=sqrt(xbr(ibr)*xbr(ibr)+
     $                    ybr(ibr)*ybr(ibr)+zbr(ibr)*zbr(ibr))
              end do
              brmin=minval(br3(1:8)); brmax=maxval(br3(1:8))
#endif

!             if(brmin.le.roa.and.brmax.ge.roa) refine(i)=.true.
              
                    ! write(*,*)'brmin,block,j,mype',brmin,i,j,mype
              if(brmin.le.roa) then
                     refine(i)=.true.
                    ! write(*,*)'refining,block,j,mype',brmin,i,j,mype
              endif

              if(lrefine(i).gt.lrefine_max) then
                   write(*,*)'refinement level exceeded for ',i
                   stop
              endif

            endif
            end do
            call shmem_barrier_all()
            call amr_refine_derefine
        end do
        x1=coord(1,5)-bsize(1,5)/2; x2=coord(1,5)+bsize(1,5)/2;
!write(*,*)'block 5:------x1, x2 ', x1,x2,mype


!write(*,*)'main:surr_blks ',surr_blks(1,:,:,1,8)
        !  pause5
        call amr_initial_magnetosphere(mype,xminb,yminb,zminb,
     &                                      xmaxb,ymaxb,zmaxb)

         call amr_restrict(mype,iopt,iempty)
!       write(*,*)'leaving amr_initial_magnetosphere'
         bflags(1,:lnblocks) = 1

        call shmem_barrier_all()
c        ldiag = .true.
c        call amr_gsurrounding_blks(mype,ldiag)
#else
        write(*,*)'2:REFINE_UNIFORMILY' 
	call amr_initial_magnetosphere(mype,xminb,yminb,zminb,
     &                                      xmaxb,ymaxb,zmaxb)
c         bflags(1,:lnblocks) = 1
#if N_DIM == 3
            do ib = 1,lnblocks
               
               if (nodetype(ib).eq.1) then
                  rmax = maxval(unk(1,:,:,:,ib) )
 !MG              write(*,*)'1:rmax,ib=',rmax,ib
                  
                  if (coord(1,ib).gt.0.) then
                     dx = coord(1,ib) - bsize(1,ib)/2.
                  else
                     dx = coord(1,ib) + bsize(1,ib)/2.
                  end if
                  if (coord(2,ib).gt.0.) then
                     dy = coord(2,ib) - bsize(2,ib)/2.
                  else
                     dy = coord(2,ib) + bsize(2,ib)/2.
                  end if
                  if (coord(3,ib).gt.0.) then
                     dz = coord(3,ib) - bsize(3,ib)/2.
                  else
                     dz = coord(3,ib) + bsize(3,ib)/2.
                  end if
                  
                  dx = abs(dx)
                  dy = abs(dy)
                  dz = abs(dz)
                  
                  r = 3.5*rearth
                  r = 2.*r
                  
                  if (dx.le.r.and.dy.le.r.and.dz.le.r) then
                     refine  (ib) = .TRUE.
                     derefine(ib) = .FALSE.
! Prevent blocks in immediate neighborhood of the earth
! from derefining at any time. To do this we set bflags(1,.)=1
! for these blocks, are test this flag in test_refinement.
                     bflags(1,ib) = 1
 !MG                 write(*,*)'bflags(1,ib) = ',bflags(1,ib) 
                  end if

                  r = 3.5*rearth
                  r = 3.*r
                  
                  if (dx.le.r.and.dy.le.r.and.dz.le.r) then
! To make blocks in immediate neighborhood of the earth
! use a first order spatial scheme, set bflags(2,ib) = 1
                     bflags(2,ib) = 1
                  end if

                  
               end if
               
            end do
#else /*NOT NDIM3*/
            do ib = 1,lnblocks
               
               if (nodetype(ib).eq.1) then
                  
                  if (coord(1,ib).gt.0.) then
                     dx = coord(1,ib) - bsize(1,ib)/2.
                  else
                     dx = coord(1,ib) + bsize(1,ib)/2.
                  end if
                  if (coord(2,ib).gt.0.) then
                     dy = coord(2,ib) - bsize(2,ib)/2.
                  else
                     dy = coord(2,ib) + bsize(2,ib)/2.
                  end if
                  
                  dx = abs(dx)
                  dy = abs(dy)
                  
                  r  = 10.*rearth
                  r  = 2.*r
                  
                  r  = 30.*rearth
                  
                  if (dx.le.r.and.dy.le.r) then
                     refine(ib)   = .TRUE.
                     derefine(ib) = .FALSE.
! Prevent blocks in immediate neighborhood of the earth
! from derefining at any time. To do this we set bflags(1,.)=1
! for these blocks, are test this flag in test_refinement.
                     bflags(1,ib) = 1
 !MG                 write(*,*)'bflags(1,ib) = ',bflags(1,ib) 
                   
                  end if

                  r = 10.0*rearth
                  r = 4.*r
                  
                  if (dx.le.r.and.dy.le.r) then
! To make blocks in immediate neighborhood of the earth
! use a first order spatial scheme set bflags(2,ib) = 1
                     bflags(2,ib) = 1
                  end if
                  
               end if
            
            end do
#endif /*NDIM3*/

        call shmem_barrier_all()
c        ldiag = .true.
c        call amr_gsurrounding_blks(mype,ldiag)

#endif /*INITIAL_REFINE_MESH*/
!       write(*,*)'LEAVING REFINE_UNIFORMILY' 

#else
        write(*,*)'3:REFINE_UNIFORMILY' 
!
! Establish solution on the initial grid.
!
!-------Create initial mesh
!
      write(*,*)'entering initialmesh'
      call initialmesh(mype,xminb,yminb,zminb,xmaxb,ymaxb,zmaxb)
c      write(*,*)'entering strongfield_initial'
c      call strongfield_initial(mype)
c      write(*,*)'leaving strongfield_initial'
!-----------------------------------------------------------------------
      write(*,*)'1:entering initialstate'
      call initialstate(mype,xmin,ymin,zmin,xmax,ymax,zmax)
      write(*,*)'leaving initialstate'

#ifdef MPI_USED
      call mpi_amr_global_domain_limits

      boundary_index(1:6) = -21
!      boundary_index(1) = -22
      boundary_index(1) = -38
      boundary_index(3) = -40
      boundary_index(4) = -50

! x boundaries 
      boundary_box(1,2:3,1:2) = -1.e30
      boundary_box(2,2:3,1:2) =  1.e30
      boundary_box(1,1,1) = -1.e30
      boundary_box(2,1,1) = grid_xmin
      boundary_box(1,1,2) = grid_xmax
      boundary_box(2,1,2) =  1.e30
! y boundaries 
      if(ndim.ge.2) then
      boundary_box(1,1,3:4) = -1.e30
      boundary_box(2,1,3:4) =  1.e30
      boundary_box(1,3,3:4) = -1.e30
      boundary_box(2,3,3:4) =  1.e30
      boundary_box(1,2,3) = -1.e30
      boundary_box(2,2,3) = grid_ymin
      boundary_box(1,2,4) = grid_ymax
      boundary_box(2,2,4) =  1.e30
      endif
! z boundaries 
      if(ndim.eq.3) then
      boundary_box(1,1:2,5:6) = -1.e30
      boundary_box(2,1:2,5:6) =  1.e30
      boundary_box(1,3,5) = -1.e30
      boundary_box(2,3,5) = grid_zmin
      boundary_box(1,3,6) = grid_zmax
      boundary_box(2,3,6) =  1.e30
      endif
#endif /* MPI_USED */

      call shmem_barrier_all()

! Now cycle over blocks, refining all existing leaf blocks.
!
!-----store information concerning refinement
        lurefd(:nuref)= luref(:nuref)
        nurefd        = nuref
        luref(1)      = 1
        nuref         = 1
#ifndef INITIAL_REFINE_MESH
         bflags(1,:lnblocks) = 1
#else
! Now cycle over blocks adjusting refinement of initial setup as required

         lrefine_again = .true.
         istep         = 0
         istep_max     = lrefine_initial
         
      write(*,*)'entering do while loop'
         
      do while(lrefine_again.and.(istep.lt.istep_max))

! refine grid and apply morton reordering to grid blocks if necessary

c                              refine(:)=.true.  ! temp test mod - uniform
 !PMN if(istep.lt.istep_max-1) refine(:)=.true.  ! temp test mod - nonuniform
      if(istep.lt.istep_max-2) refine(:)=.true.  ! temp test mod - nonuniform

c      write(*,*)'entering amr_refine_derefine'
      call amr_refine_derefine
c      write(*,*)'leaving amr_refine_derefine'

      if(lnblocks.gt.0) then
        refine  (1:lnblocks) = .false.
        derefine(1:lnblocks) = .false.
      endif

      call shmem_barrier_all()

c      write(*,*)'entering initialstate:2'
      call initialstate(mype,xmin,ymin,zmin,xmax,ymax,zmax)

      where (lrefine.eq.lrefine_max) refine = .FALSE.
      lrefine_again=.false.

      call amr_test_refinement_dan(mype,nprocs,time)

      call shmem_barrier_all()

#ifdef REFINE_EARTH_INITIAL
#if N_DIM == 3
            do ib = 1,lnblocks
               
               if (nodetype(ib).eq.1) then
                  rmax = maxval(unk(1,:,:,:,ib) )
 !MG              write(*,*)'1:rmax,ib=',rmax,ib
                  
                  if (coord(1,ib).gt.0.) then
                     dx = coord(1,ib) - bsize(1,ib)/2.
                  else
                     dx = coord(1,ib) + bsize(1,ib)/2.
                  end if
                  if (coord(2,ib).gt.0.) then
                     dy = coord(2,ib) - bsize(2,ib)/2.
                  else
                     dy = coord(2,ib) + bsize(2,ib)/2.
                  end if
                  if (coord(3,ib).gt.0.) then
                     dz = coord(3,ib) - bsize(3,ib)/2.
                  else
                     dz = coord(3,ib) + bsize(3,ib)/2.
                  end if
                  
                  dx = abs(dx)
                  dy = abs(dy)
                  dz = abs(dz)
                  
                  r = 3.5*rearth
                  r = 2.*r
                  
                  if (dx.le.r.and.dy.le.r.and.dz.le.r) then
                     refine  (ib) = .TRUE.
                     derefine(ib) = .FALSE.
! Prevent blocks in immediate neighborhood of the earth
! from derefining at any time. To do this we set bflags(1,.)=1
! for these blocks, are test this flag in test_refinement.
                     bflags(1,ib) = 1
 !MG                 write(*,*)'bflags(1,ib) = ',bflags(1,ib) 
                  end if

                  r = 3.5*rearth
                  r = 3.*r
                  
                  if (dx.le.r.and.dy.le.r.and.dz.le.r) then
! To make blocks in immediate neighborhood of the earth
! use a first order spatial scheme, set bflags(2,ib) = 1
                     bflags(2,ib) = 1
                  end if

                  
               end if
               
            end do
#else /*NOT NDIM3*/
            do ib = 1,lnblocks
               
               if (nodetype(ib).eq.1) then
                  
                  if (coord(1,ib).gt.0.) then
                     dx = coord(1,ib) - bsize(1,ib)/2.
                  else
                     dx = coord(1,ib) + bsize(1,ib)/2.
                  end if
                  if (coord(2,ib).gt.0.) then
                     dy = coord(2,ib) - bsize(2,ib)/2.
                  else
                     dy = coord(2,ib) + bsize(2,ib)/2.
                  end if
                  
                  dx = abs(dx)
                  dy = abs(dy)
                  
                  r  = 10.*rearth
                  r  = 2.*r
                  
                  r  = 30.*rearth
                  
                  if (dx.le.r.and.dy.le.r) then
                     refine(ib)   = .TRUE.
                     derefine(ib) = .FALSE.
! Prevent blocks in immediate neighborhood of the earth
! from derefining at any time. To do this we set bflags(1,.)=1
! for these blocks, are test this flag in test_refinement.
                     bflags(1,ib) = 1
 !MG                 write(*,*)'bflags(1,ib) = ',bflags(1,ib) 
                   
                  end if

                  r = 10.0*rearth
                  r = 4.*r
                  
                  if (dx.le.r.and.dy.le.r) then
! To make blocks in immediate neighborhood of the earth
! use a first order spatial scheme set bflags(2,ib) = 1
                     bflags(2,ib) = 1
                  end if
                  
               end if
            
            end do
#endif /*NDIM3*/
#endif /* REFINE_EARTH_INITIAL */

      write(*,*)'leaving REFINE_EARTH_INITIAL'

! Final safety check to prevent refinement below lrefine_max
            
            do ib=1,lnblocks
               if (lrefine(ib).ge.lrefine_max) refine(ib) = .FALSE.
            enddo
            
            icheckg=0
            icheck=0
            do ib=1,lnblocks
               if( (refine(ib)).or.(derefine(ib)) ) icheck=1
            enddo

            call COMM_INT_MAX_TO_ALL(icheckg,icheck)

            if(icheckg.eq.1) lrefine_again=.true.
            istep=istep+1
            if (mype==0) print *,'main:  ISTEP = ',ISTEP

! prolong solution to any new leaf blocks if necessary
         iopt    = 1
         nlayers = nguard
c         call amr_prolong(mype,iopt,nlayers)

         if(istep.lt.istep_max-1) lrefine_again = .true. ! temp - nonuniform
c      lrefine_again = .true. ! temp - uniform

      enddo !do while



#endif /*INITIAL_REFINE_MESH*/
        luref(:nurefd)= lurefd(:nurefd)
        nuref         = nurefd
! prolong solution to any new leaf blocks if necessary
        iopt    = 1
        nlayers = nguard
!pmn        call amr_prolong(mype,iopt,nlayers)

! exchange guardcell information
#ifdef MPI_USED
      lcc = .true.
      lfc = .true.
      tag_offset = 100

      lguard    = .true.
      lprolong  = .false.
      lflux     = .false.
      leedge    = .false.
      lrestrict = .false.
      lec        =.false.
      lnc        =.false.
      lfulltree=.false.

      call mpi_amr_comm_setup(mype,nprocs,
     .                        lguard,lprolong,lflux,leedge,lrestrict,
     .                        lfulltree,
     .                        iopt,lcc,lfc,lec,lnc,tag_offset)

!       call mpi_amr_guardcell(mype,iopt,nlayers,nprocs)
        call amr_guardcell(mype,iopt,nlayers)
#else
        write(*,*)'entering amr_guardcell'
        call     amr_guardcell(mype,iopt,nlayers)
        write(*,*)'leaving amr_guardcell'
#endif /*MPI_USED*/
#endif /*REFINE_UNIFORMILY*/
#endif /*MAGNETOSPHERE*/
!write(*,*)'xminb=',xminb
c         pause
#ifdef DIVERGENCE_TEST
#ifndef REFINE_DIVB_INITIAL
! Establish solution on the initial grid.
	 call initialmesh(mype,xminb,yminb,zminb,xmaxb,ymaxb,zmaxb)

! Now cycle over blocks, refining all existing leaf blocks.
#ifdef INITIAL_REFINE_MESH

        do loop_count=1,2

         refine(1:lnblocks) = .true.

         call shmem_barrier_all()

! refine grid and apply morton reordering to grid blocks if necessary
         call amr_refine_derefine

! prolong solution to any new leaf blocks if necessary
         iopt    = 1
         nlayers = nguard
c         call amr_prolong(mype,iopt,nlayers)

        enddo
	call amr_initial_divb(mype,xminb,yminb,zminb,xmaxb,ymaxb,zmaxb)

        call shmem_barrier_all()

#else
	call amr_initial_divb(mype,xminb,yminb,zminb,xmaxb,ymaxb,zmaxb)

        call shmem_barrier_all()
#endif

#else
!
c      call strongfield_initial(mype)
! Establish solution on the initial grid.
!
!-------Create initial mesh
!
      call initialmesh(mype,xminb,yminb,zminb,xmaxb,ymaxb,zmaxb)
!-----------------------------------------------------------------------
      call initialstate(mype,xmin,ymin,zmin,xmax,ymax,zmax)

! Now cycle over blocks, refining all existing leaf blocks.
!
!-----store information concerning refinement
        lurefd(:nuref)= luref(:nuref)
        nurefd        = nuref
        luref(1)      = 1
        nuref         = 1

! Now cycle over blocks adjusting refinement of initial setup as required

         lrefine_again = .true.
         istep         = 0
         istep_max     = lrefine_initial
         
      do while(lrefine_again.and.(istep.lt.istep_max))

! refine grid and apply morton reordering to grid blocks if necessary

      call amr_refine_derefine

      if(lnblocks.gt.0) then
        refine  (1:lnblocks) = .false.
        derefine(1:lnblocks) = .false.
      endif

      call shmem_barrier_all()

      call initialstate(mype,xmin,ymin,zmin,xmax,ymax,zmax)

      where (lrefine.eq.lrefine_max) refine = .FALSE.
        lrefine_again=.false.

        call amr_test_refinement_dan(mype,nprocs,time)

#if N_DIM == 3
            do ib = 1,lnblocks
               
               if (nodetype(ib).eq.1) then
                  if (coord(1,ib)<0.0) then
                     refine  (ib) = .TRUE.
                     derefine(ib) = .FALSE.
! Prevent blocks in immediate neighborhood of the earth
! from derefining at any time. To do this we set bflags(1,.)=1
! for these blocks, are test this flag in test_refinement.
                     bflags(1,ib) = 1
                  end if
                  
               end if
               
            end do
#else /*NOT NDIM3*/
            do ib = 1,lnblocks
               
               if (nodetype(ib).eq.1) then
                  if (coord(1,ib)<0.0) then
                     refine(ib) = .TRUE.
                     derefine(ib) = .FALSE.
! Prevent blocks in immediate neighborhood of the earth
! from derefining at any time. To do this we set bflags(1,.)=1
! for these blocks, are test this flag in test_refinement.
                     bflags(1,ib) = 1
                  end if
                  
               end if
            
            end do
#endif /*NDIM3*/


! Final safety check to prevent refinement below lrefine_max
            
            do ib=1,lnblocks
               if (lrefine(ib).ge.lrefine_max) refine(ib) = .FALSE.
            enddo
            
            icheckg=0
            icheck=0
            do ib=1,lnblocks
               if( (refine(ib)).or.(derefine(ib)) ) icheck=1
            enddo

            call COMM_INT_MAX_TO_ALL(icheckg,icheck)

            if(icheckg.eq.1) lrefine_again=.true.
            istep=istep+1
            if (mype==0) print *,' ISTEP = ',ISTEP

! prolong solution to any new leaf blocks if necessary
         iopt    = 1
         nlayers = nguard
c         call amr_prolong(mype,iopt,nlayers)

      enddo !do while



        luref(:nurefd)= lurefd(:nurefd)
        nuref         = nurefd
! prolong solution to any new leaf blocks if necessary
        iopt    = 1
        nlayers = nguard
        call amr_prolong(mype,iopt,nlayers)

! exchange guardcell information
#ifdef MPI_USED
      lcc = .true.
      lfc = .true.
      tag_offset = 100

      lguard    = .true.
      lprolong  = .false.
      lflux     = .false.
      leedge    = .false.
      lrestrict = .false.
      lec        =.false.
      lnc        =.false.
      lfulltree=.false.
      call mpi_amr_comm_setup(mype,nprocs,
     .                        lguard,lprolong,lflux,leedge,lrestrict,
     .                        lfulltree,
     .                        iopt,lcc,lfc,lec,lnc,tag_offset)
!       call mpi_amr_guardcell(mype,iopt,nlayers,nprocs)
        call amr_guardcell(mype,iopt,nlayers)
#else
        write(*,*)'entering amr_guardcell'
        call     amr_guardcell(mype,iopt,nlayers)
#endif /*MPI_USED*/
#endif /*REFINE_DIVB_INITIAL*/
#endif /*DIVERGENCE_TEST*/


! end of initialization
!------------------------------------------------------------------

!------
      end if                                           ! restart
!------
#ifdef INTERIOR_BOUNDARY1
!     write(*,*)'entering INTERIOR_BOUNDARY'

!     locate and flag those blocks that are include portions
!     of the interior boundary
!
! set marker flag to denote blocks to be acquired later
        call shmem_barrier_all()

        if(lnblocks.gt.0) then
          do lb = 1,lnblocks
           if(nodetype(lb).eq.1) then
           do k=kl_bnd,ku_bnd
            do j=jl_bnd,ju_bnd
             do i=il_bnd,iu_bnd
              if(unk(1,i,j,k,lb)>=0.01*rwall)unk(1,i,j,k,lb)=2.0
             enddo
            enddo
           enddo

          endif
         enddo
        endif
        call shmem_barrier_all()

!     write(*,*)'leaving INTERIOR_BOUNDARY'
#endif /* INTERIOR_BOUNDARY */
!
!---------------------------------------------------------------
!
! BLOCK 4
!       print*,'! BLOCK 4'


! exchange guardcell information
        iopt = 1
        nlayers = nguard
        call     amr_guardcell(mype,iopt,nlayers)

#ifdef NO_PERMANENT_GUARDCELLS
! Store a global copy of the current solution in gt_unk
!     write(*,*)'entering amr_1blk_copy_soln(-1)'
       call amr_1blk_copy_soln(-1)
!     write(*,*)'leaving amr_1blk_copy_soln(-1)'
#endif

       call shmem_barrier_all()

#ifdef MG_COMMENTS
      vmin(:,:) =  1.e30
      vmax(:,:) = -1.e30
      umin(:,:) =  1.e30
      umax(:,:) = -1.e30
#endif

      if(lnblocks.gt.0) then
      do lb=1,lnblocks
#ifndef ADVANCE_ALL_LEVELS
      if(nodetype(lb).eq.1) then
#endif
        idest = 1
        do k=1,nzb
         do j=1,nyb
          do i=1,nxb
           u2(1:ne,i,j,k)=unk(1:ne,i,j,k,lb) ! time level: n 
          enddo
         enddo
        enddo

c      write(*,*)'entering strongfields_wb'
        call strongfields_wb(mype,lb)
c      write(*,*)'leaving strongfields_wb'
        nxmax = nxb+2*nguard
        nymax = nyb+2*nguard*k2d
        nzmax = nzb+2*nguard*k3d
        bxmax(lb)=maxval(b_sf(1,:nxmax,:nymax,:nzmax))
        bxmin(lb)=minval(b_sf(1,:nxmax,:nymax,:nzmax))
        bymax(lb)=maxval(b_sf(2,:nxmax,:nymax,:nzmax))
        bymin(lb)=minval(b_sf(2,:nxmax,:nymax,:nzmax))
        bzmax(lb)=maxval(b_sf(3,:nxmax,:nymax,:nzmax))
        bzmin(lb)=minval(b_sf(3,:nxmax,:nymax,:nzmax))

        call eos_u_to_v_loc(u2,v2,lb)

        do l = 1, ne
         vmax(l,lb) = maxval( u2(l,:nxb,:nyb,:nzb) )
         vmin(l,lb) = minval( u2(l,:nxb,:nyb,:nzb) )
!-----------------------------------------------------------------------
        enddo
#ifndef ADVANCE_ALL_LEVELS
      endif
#endif
      enddo ! end loop over grid blocks
      endif
c      write(*,*)'-1:shmem_my_pe()=',shmem_my_pe()

      call shmem_barrier_all()
#ifdef NEVER20
      gbxmax1=-1.0e+30
      gbxmin1= 1.0e+30
      gbymax1=-1.0e+30
      gbymin1= 1.0e+30
      gbzmax1=-1.0e+30
      gbzmin1= 1.0e+30
#ifdef MG_COMMENTS
      varmax(:)=-1.0e+30
      varmin(:)= 1.0e+30
#endif

      if(lnblocks.gt.0) then
      do lb=1,lnblocks
#ifndef ADVANCE_ALL_LEVELS
      if(nodetype(lb).eq.1) then
#endif
        gbxmax1 = max( bxmax(lb),gbxmax1 )
        gbxmin1 = min( bxmin(lb),gbxmin1 )
        gbymax1 = max( bymax(lb),gbymax1 )
        gbymin1 = min( bymin(lb),gbymin1 )
        gbzmax1 = max( bzmax(lb),gbzmax1 )
        gbzmin1 = min( bzmin(lb),gbzmin1 )
        do j=1,ne
         varmax(j) = max( vmax(j,lb),varmax(j) )
         varmin(j) = min( vmin(j,lb),varmin(j) )
        enddo
#ifndef ADVANCE_ALL_LEVELS
      endif
#endif
      enddo ! end loop over grid blocks
      endif /* lnblocks>0 */
c      write(*,*)'1:shmem_my_pe()=',shmem_my_pe()
      call shmem_barrier_all()
      call comm_real_max_to_all(gbxmax1,gbxmax)
      call shmem_barrier_all()

      call comm_real_min_to_all(gbxmin1,gbxmin)
      call shmem_barrier_all()

      call comm_real_max_to_all(gbymax1,gbymax)
      call shmem_barrier_all()

      call comm_real_min_to_all(gbymin1,gbymin)
      call shmem_barrier_all()

      call comm_real_max_to_all(gbzmax1,gbzmax)
      call shmem_barrier_all()

      call comm_real_min_to_all(gbzmin1,gbzmin)
      call shmem_barrier_all()

      if(mype==0)then
       write(*,*)'gbxmin,gbxmax=',gbxmin,gbxmax
       write(*,*)'gbymin,gbymax=',gbymin,gbymax
       write(*,*)'gbzmin,gbzmax=',gbzmin,gbzmax
      endif
#endif /* NEVER20 */

c      write(*,*)'2:shmem_my_pe()=',shmem_my_pe()
      do j=1,ne
        call comm_real_min_to_all(varmin(j),varmin(j))
        call shmem_barrier_all()
        call comm_real_max_to_all(varmax(j),varmax(j))
        call shmem_barrier_all()
      enddo

c      if(mype==0)then
c        write(*,*)'a istep,time,dt:',istep,time,dt
c       do j=1,ne
c        write(*,'(" a varmax(",i2,"),varmin(",i2,") = ",1p2e16.7)')
c     &               j,j,varmax(j),varmin(j)
c       enddo
c      endif


! Output details of the current grid
 !MG    if(mype.eq.0) write(*,*) 'pe / blk / blk-coords / blk-sizes'
 !MG    call shmem_barrier_all()

 !MG    do lb=1,lnblocks
 !MG    write(*,*) mype,lb,(coord(i,lb),i=1,ndim),
 !MG .             (bsize(j,lb),j=1,ndim)
 !MG    enddo

c        write(*,*)'3:shmem_my_pe()=',shmem_my_pe()
        call shmem_barrier_all()

!----------------------------------------------------------
!
! BLOCK 6
!       print*,'! BLOCK 6'
!
#ifdef INTERIOR_BOUNDARY
#ifdef IONOSPHERE
      call ionosphere(mype,nprocs,0)
#endif /*IONOSPHERE*/
#endif /*INTERIOR_BOUNDARY*/

 1000 continue  ! restart

!
!----------------------------------------------------------
!
! BLOCK 7
!       print*,'! BLOCK 7'

! Initialize Graphics
!--------------------
#ifdef INTERACTIVE_GRAPHICS
      if(lgraph)then
c      write(*,*)'entering open:graphics_parms 1'
           if(mype.eq.0) then
                open(unit=iunitc,form='formatted',
     &          action='read',file='graphics_parms')
                ierr=0
                read(iunitc,nml=DRIVERPARMS,ERR=131)
                close(iunitc)
                go to 132
 131            ierr=1
                write(*,*)'ERROR: in reading graphics_parms'
                close(iunitc)
 132            continue
           endif
#ifdef SGI_SHMEM
      call shm_bcast_g_parms(ierr,mype,
     &     ltecplot, nstep_tecplot, maxstp, nstep_dump, nstep_amr,
     &     nstep_graphics, nstep_times, nstep_print,
     &     iunit4, iunit7, cfl,ctore,ctode,epsil,nuref,luref)
#endif /*SGI_SHMEM*/
#ifdef MPI_USED
           call bcast_g_parms(ierr,mype,
     &     ltecplot, nstep_tecplot, maxstp, nstep_dump, nstep_amr,
     &     nstep_graphics, nstep_times, nstep_print,
     &     iunit4, iunit7,cfl,ctore,ctode,epsil,nuref,luref)
#endif /* MPI_USED */
!----------------------------------------------------------

           print*,'.........1.... lgraph,lmovie',lgraph,lmovie
#ifdef GRAPHICS
           print*,'.........2.... lgraph,lmovie',lgraph,lmovie
#ifdef QUK_GRAPHICS
           print*,'.........3.... lgraph,lmovie',lgraph,lmovie
           if(mype==0) call init_window(lgraph, lmovie)
#endif /*QUK_GRAPHICS*/

         ib   = 1
         delx = bsize(1,ib)/nxb
         dely = bsize(2,ib)/nyb
         
         wxmn = coord(1,ib) + delx*(1 - (nguard + nxb/2) - 0.5)
         wxmx = coord(1,ib) + delx*(1+nxb+2*nguard 
     &        - (nguard + nxb/2) - 0.5)
         wymn = coord(2,ib) + dely*(1 - (nguard + nyb/2) - 0.5)
         wymx = coord(2,ib) + dely*(1+nyb+2*nguard 
     &        - (nguard + nyb/2) - 0.5)
         
         wxmn = xminb
         wxmx = xmaxb
         wymn = yminb
         wymx = ymaxb

!           write(*,*)'completed 901', mype
         if(mype==0) print *,
     &           ' wxmn, wymn, wxmx, wymx=',wxmn, wymn, wxmx, wymx

!           write(*,*)'completed 902', mype

        if(ltecplot>0)then
c         call output_tecplot_fe(0,time,dt)
        endif

#ifdef MPI_USED
        call mpi_barrier (MPI_COMM_WORLD, errcode)
#endif MPI_USED


#if N_DIM == 2
!       call amr_fieldline_trace(0)
#endif

#endif /*GRAPHICS*/
       endif ! lgraph
#endif /*INTERACTIVE*/

#ifdef INTERIOR_BOUNDARY
      if(mype.eq.0)
     & write(*,*)'entering create_interior_bc'
c      call create_interior_bc
      call earth_boundary(mype)
      call shmem_barrier_all()
      if(mype.eq.0)
     &   write(*,*)'leaving create_interior_bc'
#else
      idomain = lnblocks
      iearth  = 0
      if(lnblocks.gt.0) then
      do lb=1,lnblocks
#ifndef ADVANCE_ALL_LEVELS
      if(nodetype(lb).eq.1) then
#endif
!       earth_blocks(lb)=.false.
#ifndef ADVANCE_ALL_LEVELS
      endif
#endif
      enddo ! end loop over grid blocks
      endif /* lnblocks>0 */
      call shmem_barrier_all()

#endif /*INTERIOR_BOUNDARY*/

      no_of_steps = maxstp - minstp
      minstp      = loop_start_index+1 
      maxstp      = minstp + no_of_steps

!---------------------------------------------------------

      if(mype.eq.0) then
      write(*,501) ebox_xmin, ebox_xmax, mype
      write(*,502) ebox_ymin, ebox_ymax, mype
      write(*,503) ebox_zmin, ebox_zmax, mype
      endif
 501  format(2x,"ebox-x ",f10.5,3x,f10.5, 5x,i2)
 502  format(2x,"ebox-y ",f10.5,3x,f10.5, 5x,i2)
 503  format(2x,"ebox-z ",f10.5,3x,f10.5, 5x,i2)
! loop over time steps.
      if(mype.eq.0)
     & write(*,*)'----------minstp, maxstp,lnblocks,mype=',
     &           minstp, maxstp,lnblocks,mype

!---------------------------------------------------------
      iunit30=30; iunit31=31; iunit32=32; iunit33=33; 
      iunit34=34; iunit35=35;
!---------------------------------------------------------

#ifdef onedebug
!for 1blk
!outfiles for fx, fy in block_ctu
!     open(unit=iunit30,file='unk1_1blk_f2')
!open(unit=iunit31,file='unk2_1blk_f2')
!open(unit=iunit32,file='unk3_1blk_f2')
!open(unit=iunit33,file='unk5_1blk_f2')
!open(unit=iunit34,file='unk6_1blk_f2')
!open(unit=iunit35,file='unk7_1blk_f2')
#endif

#ifdef multidebug
!for multiple blk
!     open(unit=iunit30,file='unk1_f2')
!open(unit=iunit31,file='unk2_f2')
!open(unit=iunit32,file='unk3_f2')
!open(unit=iunit33,file='unk5_f2')
!open(unit=iunit34,file='unk6_f2')
!open(unit=iunit35,file='unk7_f2')
#endif
!---------------------------------------------------------
#ifdef CHECK_DIVB

c        write(*,*)'main:in gtest_neigh_data',mype
        call gtest_neigh_data(mype,0,0.0)

#endif /*CHECK_DIVB*/
!---------------------------------------------------------
      do istep = minstp, maxstp
!---------------------------------------------------------
#ifdef onedebug
!for output of fx, fy in block_ctu
      mylb = 1
!if(istep.ge.5900)
!     call ch4_data(istep,time,dt,mylb)
#endif
#ifdef multidebug
      mylb = 1
!if(istep.ge.5900)
!     call ch4_data(istep,time,dt,mylb)
#endif
!--------------------------------------------
     
!
      vmin(:,:) = +1.e30
      vmax(:,:) = -1.e30
      call shmem_barrier_all()
!-----timeadvance

      if(lnblocks.gt.0) then
      do lb=1,lnblocks
#ifndef ADVANCE_ALL_LEVELS
      if(nodetype(lb).eq.1) then
#endif
        do k=1,nzb
         do j=1,nyb
          do i=1,nxb
           u2(1:ne,i,j,k)=unk(1:ne,i,j,k,lb) ! time level: n 
          enddo
         enddo
        enddo

        call strongfields_wb(mype,lb)

        call eos_u_to_v_loc(u2,v2,lb)
        
c        write(*,*)'athena:neg den:',v2(1,21,31,1),u2(1,21,31,1)
c        stop
        
        do l = 1, ne
         vmax(l,lb) = maxval( v2(l,:nxb,:nyb,:nzb) )
         vmin(l,lb) = minval( v2(l,:nxb,:nyb,:nzb) )
        enddo
        do l = 1, ne
         umax(l,lb) = maxval( u2(l,:nxb,:nyb,:nzb) )
         umin(l,lb) = minval( u2(l,:nxb,:nyb,:nzb) )
        enddo
#ifdef INTERIOR_BOUNDARY
#define MASKIT
#endif 
#ifdef MASKIT
        do k=1,nzb
          kk = nguard*k3d+k
         do j=1,nyb
            jj = nguard*k2d+j
          do i=1,nxb
             ii = nguard*k1d+i
           v2(:,i,j,k)=v2(:,i,j,k)*bmask(ii,jj,kk,lb) ! time level: n 
          enddo
         enddo
        enddo
#endif

 
#ifdef GRAPHICS
        do k=1,nzb
         do j=1,nyb
          do i=1,nxb
           vnk(1:ne,i,j,k,lb)=v2(1:ne,i,j,k)
          enddo
         enddo
        enddo

#endif /*GRAPHICS*/

#ifndef ADVANCE_ALL_LEVELS
      endif
#endif
      enddo ! end loop over grid blocks
      endif /* lnblocks>0 */
      call shmem_barrier_all()

c      write(*,*)'@ mod: shmem_my_pe()=',
c     &     shmem_my_pe(),istep,nstep_print

      if(mod(istep,nstep_print) .eq. 0 ) then
c        write(*,*)'istep=',istep, shmem_my_pe()


       varmin(:) =  1.0e+30
       varmax(:) = -1.0e+30

      if(lnblocks.gt.0) then
       do lb=1,lnblocks
#ifndef ADVANCE_ALL_LEVELS
        if(nodetype(lb).eq.1) then
#endif
         do j=1,ne
          varmax(j) = maxval( vmax(j,1:lnblocks) )
          varmin(j) = minval( vmin(j,1:lnblocks) )
         enddo
#ifdef CONSERVED_DIAGNOSTICS
         do j=1,ne
          uarmax(j) = maxval( umax(j,1:lnblocks) )
          uarmin(j) = minval( umin(j,1:lnblocks) )
         enddo
#endif /*CONSERVED_DIAGNOSTICS*/

#ifndef ADVANCE_ALL_LEVELS
        endif
#endif
       enddo
      endif /* lnblocks>0 */
      call shmem_barrier_all()

!mk10
      do j=1,ne
         call shmem_barrier_all()
        call comm_real_min_to_all(varmin(j),varmin(j))
         call shmem_barrier_all()
        call comm_real_max_to_all(varmax(j),varmax(j))
         call shmem_barrier_all()
#ifdef CONSERVED_DIAGNOSTICS
        call comm_real_min_to_all(uarmin(j),uarmin(j))
         call shmem_barrier_all()
        call comm_real_max_to_all(uarmax(j),uarmax(j))
         call shmem_barrier_all()
#endif /*CONSERVED_DIAGNOSTICS*/
      enddo 
         call shmem_barrier_all()
c        write(*,*)'60:mype,nprocs=',mype,nprocs

         do i=1,nxb
            do j=1,nyb
               if(v2(1,i,j,1).lt.0.) then
c                  write(*,*) nxb,nyb,i,j,'athena:neg den:',v2(1,i,j,1),u2(1,i,j,1)
c                  stop
               endif
            enddo
         enddo

      if (mype==0) then
       write(*,*)'istep,time,dt = ',istep,time,dt
       do j=1,ne
       write(*,*)'v1(',j,') = ',
     &               varmax(j),varmin(j)
       enddo
#ifdef CONSERVED_DIAGNOSTICS
       do j=1,ne
       write(*,*)'u1(',j,'),u1(',j,') = ',
     &               uarmax(j),uarmin(j)
       enddo
#endif /*CONSERVED_DIAGNOSTICS*/

      endif

      endif  ! nstep_print

c      write(*,*)'leaving lnblocks loop'

      if(lgraph)then
#ifdef INTERACTIVE_GRAPHICS
!-------------------------------------------------------
           if(mype.eq.0) then
                open(unit=iunitc,form='formatted',
     &          action='read',file='graphics_parms')
                ierr=0
                read(iunitc,nml=DRIVERPARMS,ERR=133)
                close(iunitc)
                go to 134
 133            ierr=1
                write(*,*)'ERROR: in reading graphics_parms'
                close(iunitc)
 134            continue
           endif
#ifdef SGI_SHMEM
      call shm_bcast_g_parms(ierr,mype,
     &     ltecplot, nstep_tecplot, maxstp, nstep_dump, nstep_amr,
     &     nstep_graphics, nstep_times, nstep_print,
     &     iunit4, iunit7, cfl,ctore,ctode,epsil,nuref,luref)
#endif /*SGI_SHMEM*/
#ifdef MPI_USED
           call bcast_g_parms(ierr,mype,
     &     ltecplot, nstep_tecplot, maxstp, nstep_dump, nstep_amr,
     &     nstep_graphics, nstep_times, nstep_print,
     &     iunit4, iunit7, cfl,ctore,ctode,epsil,nuref,luref)
#endif /* MPI_USED */
!-------------------------------------------------------

#ifdef GRAPHICS
         if(lgraph.and. mod(istep,nstep_graphics) .eq. 0 ) then
         ebox(1,1) = ebox_xmin
         ebox(2,1) = ebox_ymin
         ebox(3,1) = ebox_zmin
         ebox(1,2) = ebox_xmax
         ebox(2,2) = ebox_ymax
         ebox(3,2) = ebox_zmax
#if N_DIM == 2
         write(*,*)'68:mype,nprocs=',mype,nprocs
!        call amr_fieldline_trace(istep)
#endif
c        write(*,*)'Entering grdrive'
c         call grdrive(nprocs,mype,vnk,varmin,varmax,ebox,
c     $          istep,lgraph, nstep_graphics, lmovie,nmovie,movie_si)
c        write(*,*)'Leaving grdrive ... the graphics have been drawn.'
         endif
#endif /*GRAPHICS*/
c        write(*,*)'70:mype,nprocs=',mype,nprocs

#endif /*INTERACTIVE*/
      endif ! lgraph

        if(ltecplot>0.and.mod(istep,nstep_tecplot) .eq. 0 ) then
c         call output_tecplot_fe(istep,time,dt)
        endif
!       print*,'JHM 1'
#ifdef INTERIOR_BOUNDARY
#ifdef IONOSPHERE
        call ionosphere(mype,nprocs,istep)
#endif /*IONOSPHERE*/
#endif /*INTERIOR_BOUNDARY*/

! advance the solution using a User provided routine
      vmin(:,:) =  1.e30
      vmax(:,:) = -1.e30

      if(lnblocks.gt.0) then
      do lb=1,lnblocks
#ifndef ADVANCE_ALL_LEVELS
      if(nodetype(lb).eq.1) then
#endif
        do k=1,nzb
         do j=1,nyb
          do i=1,nxb
           u2(1:ne,i,j,k)=unk(1:ne,i,j,k,lb) ! time level: n 
          enddo
         enddo
        enddo
        do l = 1, ne
         vmax(l,lb) = maxval( u2(l,:nxb,:nyb,:nzb) )
         vmin(l,lb) = minval( u2(l,:nxb,:nyb,:nzb) )
        enddo
#ifndef ADVANCE_ALL_LEVELS
      endif
#endif
      enddo ! end loop over grid blocks
      endif /* lnblocks>0 */
c        write(*,*)'80:mype,nprocs=',mype,nprocs
!       print*,'JHM 2'

#ifdef NEVER200
      do j=1,ne
       varmax(j) = maxval( vmax(j,1:lnblocks) )
       varmin(j) = minval( vmin(j,1:lnblocks) )
      enddo

      do j=1,ne
        call comm_real_min_to_all(varmin(j),varmin(j))
        call shmem_barrier_all()
        call comm_real_max_to_all(varmax(j),varmax(j))
        call shmem_barrier_all()
      enddo

      if(mype.eq.0)then
        write(*,*)'g istep,time,dt:',istep,time,dt
       do j=1,ne
        write(*,*)'g varmax(',j,'),varmin(',j,') = ',
     &               varmax(j),varmin(j)
       enddo
      endif
#endif /*NEVER200*/



c      write(*,*) bflags(1,:20)
c      write(*,*)
c      write(*,*) lrefine(:20)


!       write(*,*)'mk01: main: mype, istep= ',mype, istep

        call advance_soln(istep,mype,nprocs,time,dt)

!       write(*,*)'mk02: main: mype, istep= ',mype, istep


#ifdef NEVER200

      vmin(:,:) =  1.e30
      vmax(:,:) = -1.e30

      if(lnblocks.gt.0) then
      do lb=1,lnblocks
#ifndef ADVANCE_ALL_LEVELS
      if(nodetype(lb).eq.1) then
#endif
!       print*,'JHM 3'
        do k=1,nzb
         do j=1,nyb
          do i=1,nxb
           u2(1:ne,i,j,k)=unk(1:ne,i,j,k,lb) ! time level: n 
          enddo
         enddo
        enddo
        do l = 1, ne
         vmax(l,lb) = maxval( u2(l,:nxb,:nyb,:nzb) )
         vmin(l,lb) = minval( u2(l,:nxb,:nyb,:nzb) )
        enddo
#ifndef ADVANCE_ALL_LEVELS
      endif
#endif
      enddo ! end loop over grid blocks
      endif /* lnblocks>0 */

c        write(*,*)'90:mype,nprocs=',mype,nprocs

      do j=1,ne
       write(*,*)'f:max:j,lb=',j,maxloc( vmax(j,1:lnblocks) )
       write(*,*)'f:min:j,lb=',j,minloc( vmin(j,1:lnblocks) )
      enddo

      do j=1,ne
       varmax(j) = maxval( vmax(j,1:lnblocks) )
       varmin(j) = minval( vmin(j,1:lnblocks) )
      enddo
      do j=1,ne
        call shmem_barrier_all()
        call comm_real_min_to_all(varmin(j),varmin(j))
        call shmem_barrier_all()
        call comm_real_max_to_all(varmax(j),varmax(j))
      enddo
      call shmem_barrier_all()

      if(mype==0)then
        write(*,*)'f istep,time,dt:',istep,time,dt
       do j=1,ne
        write(*,*)'f varmax(',j,'),varmin(',j,') = ',
     &               varmax(j),varmin(j)
       enddo
      endif
#endif /*NEVER200*/

#ifdef CHECK_DIVB
        write(*,*)'main:in gtest_neigh_data',mype
        call gtest_neigh_data(mype,istep,0.0)
#endif /*CHECK_DIVB*/
      
c        if(lgraph .and. mod(istep,nmovie).eq.0)then
        if(lgraph .and. (time>=timstp))then
#if N_DIM == 2
!       if(mod(istep,50000) .eq. 0 ) 
!    &    call amr_fieldline_trace(istep)
#endif

#ifdef CHOMBO_VIS

         write(*,*)'writing ChomboVis file'
         nchombo=nchombo+1
         call plotfile_chombo(nprocs,mype,nchombo,time)
         timstp = time+deltim
#endif /*CHOMBO_VIS*/
        endif
!       print*,'JHM 4'
!
! test to see if additional refinement or derefinement is necessary
! note - a call to guardcell must come before this call to ensure that
! the refinement test can be done on parents of leafs also. This avoids
! a potential refinement/derefinement flip-flop happenning on successive
! timesteps.
c        write(*,*)'istep,nstep_amr=',istep,nstep_amr
c        write(*,*)'100:mype,nprocs=',mype,nprocs
        if( istep.ne.1.and.mod(istep, nstep_amr) .eq. 0 ) then
c        write(*,*)'inside refinement loop:istep,nstep_amr=',istep,nstep_amr
 
#ifdef RUNTIME_REFINE_MESH

#ifndef ADVANCE_ALL_LEVELS
! A valid solution will be required on the parents of leaf blocks
! when refinement testing is done. See the comment before the call
! to amr_test_refinement.
        iempty = 0
        iopt   = 1
#ifdef MPI_USED
!call amr_restrict(mype,iopt,iempty,nprocs)
        call amr_restrict(mype,iopt,iempty)
#else
c        write(*,*)'entering amr_restrict:mype,nprocs=',mype,nprocs
        call amr_restrict(mype,iopt,iempty)
#endif /* MPI_USED */
#endif /* ADVANCE_ALL_LEVELS */

        call amr_test_refinement_dan(mype,nprocs,time)

! refine grid and apply morton reordering to grid blocks if necessary
c        write(*,*)'a:entering amr_refine_derefine'
        call amr_refine_derefine
c        write(*,*)'a:leaving amr_refine_derefine'

! Output details of the current grid
 !MG    if(mype==0) write(*,*) 'pe / blk / blk-coords / blk-sizes'
         call shmem_barrier_all()

 !MG    do lb=1,lnblocks
 !MG    write(*,*) mype,lb,(coord(i,lb),i=1,ndim),
 !MG .             (bsize(j,lb),j=1,ndim)
 !MG    enddo

! prolong solution to any new leaf blocks if necessary
        iopt    = 1
        nlayers = nguard
c        write(*,*)'amr_prolong'
        call amr_prolong(mype,iopt,nlayers)


#ifdef INTERIOR_BOUNDARY
      write(*,*)'entering earth_boundary'
c      call create_interior_bc
      call earth_boundary(mype)
      call shmem_barrier_all()
!     print*,'JHM 5'

      write(*,*)'leaving earth_boundary'
#endif /*INTERIOR_BOUNDARY*/

#endif /* RUNTIME_REFINE_MESH */

! exchange guardcell information
!#ifdef MPI_USED
!        call mpi_amr_guardcell(mype,iopt,nlayers,nprocs)
!#else


        call     amr_guardcell(mype,iopt,nlayers)
!#endif /*MPI_USED*/
        endif /* istep loop */ ! refinement loop

!-------------------------------------------------------------
!


#ifdef NEVER20
      vmin(:,:) =  1.e30
      vmax(:,:) = -1.e30

      if(lnblocks.gt.0) then
      do lb=1,lnblocks
#ifndef ADVANCE_ALL_LEVELS
      if(nodetype(lb).eq.1) then
#endif
        do k=1,nzb
         do j=1,nyb
          do i=1,nxb
           u2(1:ne,i,j,k)=unk(1:ne,i,j,k,lb) ! time level: n 
          enddo
         enddo
        enddo
        do l = 1, ne
         vmax(l,lb) = maxval( u2(l,:nxb,:nyb,:nzb) )
         vmin(l,lb) = minval( u2(l,:nxb,:nyb,:nzb) )
        enddo
#ifndef ADVANCE_ALL_LEVELS
      endif
#endif
      enddo ! end loop over grid blocks
      endif /* lnblocks>0 */

      do j=1,ne
       varmax(j) = maxval( vmax(j,1:lnblocks) )
       varmin(j) = minval( vmin(j,1:lnblocks) )
      enddo
      do j=1,ne
        call shmem_barrier_all()
        call comm_real_min_to_all(varmin(j),varmin(j))
        call shmem_barrier_all()
        call comm_real_max_to_all(varmax(j),varmax(j))
      enddo
        call shmem_barrier_all()

!     print*,'JHM 6'
      gbxmax1=-1.0e+30
      gbxmin1= 1.0e+30
      gbymax1=-1.0e+30
      gbymin1= 1.0e+30
      gbzmax1=-1.0e+30
      gbzmin1= 1.0e+30
#ifdef MG_COMMENTS
      varmax(:)=-1.0e+30
      varmin(:)= 1.0e+30
#endif

      if(lnblocks.gt.0) then
      do lb=1,lnblocks
#ifndef ADVANCE_ALL_LEVELS
      if(nodetype(lb).eq.1) then
#endif

        call strongfields_wb(mype,lb)
        nxmax = nxb+2*nguard
        nymax = nyb+2*nguard*k2d
        nzmax = nzb+2*nguard*k3d
        bxmax(lb)=maxval(b_sf(1,:nxmax,:nymax,:nzmax))
        bxmin(lb)=minval(b_sf(1,:nxmax,:nymax,:nzmax))
        bymax(lb)=maxval(b_sf(2,:nxmax,:nymax,:nzmax))
        bymin(lb)=minval(b_sf(2,:nxmax,:nymax,:nzmax))
        bzmax(lb)=maxval(b_sf(3,:nxmax,:nymax,:nzmax))
        bzmin(lb)=minval(b_sf(3,:nxmax,:nymax,:nzmax))

#ifndef ADVANCE_ALL_LEVELS
      endif
#endif
      enddo ! end loop over grid blocks
      endif /* lnblocks>0 */

       gbxmax = maxval( bxmax(1:lnblocks) )
       gbxmin = minval( bxmin(1:lnblocks) )
       gbymax = maxval( bymax(1:lnblocks) )
       gbymin = minval( bymin(1:lnblocks) )
       gbzmax = maxval( bzmax(1:lnblocks) )
       gbzmin = minval( bzmin(1:lnblocks) )

#endif /* NEVER20 */
!-------------------------------------------------------------
!
!----- Is it time to dump a checkpoint file?
!
c        write(*,*)'Is it time to dump a checkpoint file?'
        if( mod(istep,2*nstep_dump) .eq. 0) then
         if(mype==0) then
          print *,' Checkpointing: istep,time,dt= ',
     &              istep,time,dt,' on unit',iunit5
         endif
         call amr_checkpoint_wr_a(dt,time,istep,nchombo,iunit5)
        endif
        if( mod(istep,2*nstep_dump) .eq. nstep_dump) then
         if(mype==0) then
          print *,' Checkpointing: istep,time,dt= ',
     &              istep,time,dt,' on unit',iunit6
         endif
        call amr_checkpoint_wr_a(dt,time,istep,nchombo,iunit6)
        endif


        end do  ! time loop
!       print*,'JHM 7'
!
!---------------------------------------------------------------
        if(mype==0) then
          print *,' Checkpointing: istep,time,dt= ',
     &              istep,time,dt,' on unit',iunit7
        endif
        call amr_checkpoint_wr_a(dt,time,istep-1,nchombo,iunit7)

c        call output_tecplot_fe(istep,time,dt)
        call amr_close

        stop
        end

#ifdef MPI_USED
      subroutine bcast_g_parms(ierr,mype,
     &     ltecplot, nstep_tecplot, maxstp, nstep_dump, nstep_amr,
     &     nstep_graphics, nstep_times, nstep_print,
     &     iunit4, iunit7, cfl,ctore,ctode,epsil,nuref,luref)

      implicit none
      include "mpif.h"

      integer :: ltecplot, nstep_tecplot, maxstp, nstep_dump
      integer :: nstep_amr, nstep_graphics, nstep_times
      integer :: nstep_print, iunit4, iunit7, ierr, mype
      integer :: nuref,luref(20)
      real :: cfl,ctore,ctode,epsil 
      integer iroot, ierror, num, ib(12),nmovie

      call shmem_barrier_all()
      iroot = 0
      ib(1)=ierr; ib(2)=ltecplot; ib(3)=nstep_tecplot; 
      ib(4)=nstep_dump; ib(5)=nstep_amr;
      ib(6)=nstep_graphics; ib(7)=nstep_times; ib(8)=nstep_print;
      ib(9)=iunit4; ib(10)=iunit7;
      ib(11)=nuref; 

      num=11
      call mpi_bcast(ib(1),num,MPI_INTEGER,iroot,MPI_COMM_WORLD,
     &         ierror)
      call mpi_bcast(luref(1),20,MPI_INTEGER,iroot,MPI_COMM_WORLD,
     &         ierror)
      num=1
      call mpi_bcast(cfl,num,MPI_DOUBLE_PRECISION,
     &         iroot,MPI_COMM_WORLD,ierror)
      call mpi_bcast(ctore,num,MPI_DOUBLE_PRECISION,
     &         iroot,MPI_COMM_WORLD,ierror)
      call mpi_bcast(ctode,num,MPI_DOUBLE_PRECISION,
     &         iroot,MPI_COMM_WORLD,ierror)
      call mpi_bcast(epsil,num,MPI_DOUBLE_PRECISION,
     &         iroot,MPI_COMM_WORLD,ierror)
      
      ierr=ib(1)
      if(ierr.eq.1) then
           write(*,*)'ERROR in reading graphics.parms'
           call amr_close
           stop
      endif

      ltecplot      =ib( 2);   nstep_tecplot=ib( 3); 
      nstep_dump    =ib( 4);   nstep_amr    =ib( 5);
      nstep_graphics=ib( 6);   nstep_times  =ib( 7);   nstep_print=ib( 8);
      iunit4        =ib( 9);   iunit7       =ib(10);   nuref      =ib(11);
      nmovie        =ib(12);


!      write(*,*)'-----------------------------------------'
!      write(*,*)'           (mpi) mype = ',mype
!      write(*,*)ltecplot,nstep_tecplot,maxstp,nstep_dump,nstep_amr,
!    & mype
!      write(*,*)nstep_graphics,nstep_times,nstep_print,iunit4,iunit7,
!    & mype
!      write(*,*)cfl,mype

      return
      end
#endif

#ifdef SGI_SHMEM
      subroutine shm_bcast_g_parms(ierr,mype,
     &     ltecplot, nstep_tecplot, maxstp, nstep_dump, nstep_amr,
     &     nstep_graphics, nstep_times, nstep_print,
     &     iunit4, iunit7, cfl,ctore,ctode,epsil,nuref,luref)

      implicit none

      integer :: ltecplot, nstep_tecplot, maxstp, nstep_dump
      integer :: nstep_amr, nstep_graphics, nstep_times
      integer :: nstep_print, iunit4, iunit7, ierr, mype
      integer :: nuref,luref(20)
      real    :: cfl,ctore,ctode,epsil 
      integer iroot, ierror, num


      iroot = 0; num=1
      call shmem_barrier_all()
      if(mype.gt.0)  then
         call shmem_integer_get(ierr ,ierr  ,num,iroot)
         call shmem_integer_get(ltecplot ,ltecplot  ,num,iroot)
         call shmem_integer_get(nstep_tecplot,nstep_tecplot,num,iroot)
         call shmem_integer_get(maxstp ,maxstp  ,num,iroot)
         call shmem_integer_get(nstep_dump ,nstep_dump  ,num,iroot)
         call shmem_integer_get(nstep_amr ,nstep_amr  ,num,iroot)
         call shmem_integer_get(nstep_graphics,nstep_graphics,num,iroot)
         call shmem_integer_get(nstep_times ,nstep_times  ,num,iroot)
         call shmem_integer_get(nstep_print ,nstep_print  ,num,iroot)
         call shmem_integer_get(iunit4 ,iunit4  ,num,iroot)
         call shmem_integer_get(iunit7 ,iunit7  ,num,iroot)
         call shmem_integer_get(nuref  ,nuref   ,num,iroot)
         call shmem_integer_get(luref  ,luref   ,20,iroot)
         call shmem_real_get(cfl     ,cfl     ,num,iroot)
         call shmem_real_get(ctore   ,ctore   ,num,iroot)
         call shmem_real_get(ctode   ,ctode   ,num,iroot)
         call shmem_real_get(epsil   ,epsil   ,num,iroot)
      endif

      if(ierr.eq.1) then
           write(*,*)'ERROR in reading graphics.parms'
           call amr_close
           stop
      endif

!      write(*,*)'-----------------------------------------'
!      write(*,*)'             mype = ',mype
!      write(*,*)ltecplot,nstep_tecplot,maxstp,nstep_dump,nstep_amr,
!    & mype
!      write(*,*)nstep_graphics,nstep_times,nstep_print,iunit4,iunit7,
!    & mype
!      write(*,*)cfl,mype

      return
      end
#endif

