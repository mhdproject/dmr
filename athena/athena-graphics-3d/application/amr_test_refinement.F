#include "paramesh_preprocessor.fh"
#define DONT_REFINE_EARTH
      subroutine amr_test_refinement_dan(mype,nprocs,time)
      use physcons
      use paramesh_dimensions
      use physicaldata
      use strongfields
      use grid
      use interior_gbc

! include file defining the tree
      use tree
      use workspace
      use athena_specific_parameters




      use paramesh_interfaces
#ifdef MPI_USED
      use paramesh_mpi_interfaces
      include "mpif.h"
      integer :: tag_offset,max_blks_sent
#endif /* MPI_USED */

! $RCSfile: amr_test_refinement.F,v $
! $Revision: 1.3 $
! $Date: 2003/11/13 17:13:43 $


!--------------------------------------------------------------------------
! 
! This is a template to assist in constructing the routine AMR_TEST_REFINEMENT
! for use in your application. In this illustration we use the workspace
! array WORK to store the data which is used in computing the error measure
! at each grid point. This gives us the freedom to extend the testing
! beyond the normal bounds of individual blocksi, since WORK is declared
! with NGUARD_WORK guard cells at each boundary, which can be set to a
! larger number than NGUARD.

! Arguments:
!      mype integer local processor number
!      lrefine_min integer minimum refinement level to be permitted
!      lrefine_max integer maximum refinement level to be permitted

!--------------------------------------------------------------------------

      integer :: mype,nprocs

#include "amr_shmem.fh"
      integer :: iref
      integer, parameter :: ng0 = nguard*npgs
      
      real,save,dimension(   il_bnd:iu_bnd,
     &                  jl_bnd:ju_bnd,
     &                  kl_bnd:ku_bnd)::Omegax,Omegay,Omegaz
      real,save,dimension( 3,il_bnd:iu_bnd,
     &                  jl_bnd:ju_bnd,
     &                  kl_bnd:ku_bnd)::vel
      real,save,dimension(ne,il_bnd:iu_bnd,
     &                  jl_bnd:ju_bnd,
     &                  kl_bnd:ku_bnd)::u1,v1
      real, dimension(maxblocks) :: eerror
      real :: errort
      save errort, eerror
!
!-----refinement parameters
!
      integer, parameter :: muref = 20
      integer :: nuref,luref
      real    :: ctore,ctode,epsil
      common/refine_parameters/ctore,ctode,epsil,luref(muref),nuref
c      real :: delx,dely,delz

      common/work_bc_control/ ibc_work
      integer :: ibc_work

!-----MPI parameter
      logical :: lguard,lprolong,lflux,ledge,lrestrict,lcc,lfc
      logical :: lec,lnc, lfulltree
      integer :: iopt, ii,jj,kk

      integer, parameter :: ndel = 0

!      write (*,*) 'where athena 2 d fails on 27th file'


!-----------------------------------------------------------------------
!
! Error56

!
! Re-initialize the refinement and derefinement flag arrays
      refine(:)   = .false.
      derefine(:) = .false.

c      ctode = 0.50
c      ctore = 0.35
c      ctode = 0.4
c      ctore = .8
c       if(mype.eq.0)then
c                write(*,*)'nuref,mype,ir=',nuref,mype,luref
c       endif
c       if(mype.eq.1)then
c                write(*,*)'nuref,mype,ir=',nuref,mype,luref
c       endif
      eerror(:)=0.0
      do iuref=1,nuref
       ir = luref(iuref)  ! specifies what needs to be refined on

!
! Set up the workspace array WORK to store the variable we wish to examine 
! in order to test the refinement level.

      if(ir.eq.5)then
!
!-----refine on current density
!
c      write(*,*)'entering block loop'
      if(lnblocks.gt.0) then
      do lb=1,lnblocks
      if(nodetype(lb).eq.1.or.nodetype(lb).eq.2) then
        Omegax = 0.
        Omegay = 0.
        Omegaz = 0.
        delx   = bsize(1,lb)/float(nxb)
        dely   = bsize(2,lb)/float(nyb)
#if N_DIM == 3
        delz   = bsize(3,lb)/float(nzb)
#endif
        delx   = 1./(2.*delx)
        dely   = 1./(2.*dely)
#if N_DIM == 3 
        delz   = 1./(2.*delz)
#else
        delz   = 0.0
#endif

        call strongfields_wb(mype,lb)
        do k=kl_bnd+ndel*k3d,ku_bnd+ndel*k3d
         do j=jl_bnd+ndel*k2d,ju_bnd+ndel*k2d
          do i=il_bnd+ndel*k1d,iu_bnd+ndel*k1d
           vel(1,i,j,k) = unk(6,i-ndel,j-ndel*k2d,k-ndel*k3d,lb)
           vel(2,i,j,k) = unk(7,i-ndel,j-ndel*k2d,k-ndel*k3d,lb)
           vel(3,i,j,k) = unk(8,i-ndel,j-ndel*k2d,k-ndel*k3d,lb) 
          end do
         end do
        end do
        do k=kl_bnd+ndel*k3d+k3d,ku_bnd+ndel*k3d-k3d
         do j=jl_bnd+ndel*k2d+k2d,ju_bnd+ndel*k2d-k2d
          do i=il_bnd+ndel*k1d+k1d,iu_bnd+ndel*k1d-k1d

! Curl V = Omega

                 Omegax(i,j,k) = 
     &           (vel(3,i,j+1,k) - vel(3,i,j-1,k))*dely
#if N_DIM == 3
                 Omegax(i,j,k) = Omegax(i,j,k) -
     &           (vel(2,i,j,k+1) - vel(2,i,j,k-1))*delz
#endif

#if N_DIM == 3
                 Omegay(i,j,k) = 
     &           (vel(1,i,j,k+1) - vel(1,i,j,k-1))*delz
#endif
                 Omegay(i,j,k) = Omegay(i,j,k) -
     &           (vel(3,i+1,j,k) - vel(3,i-1,j,k))*delx

                 Omegaz(i,j,k) = 
     &           (vel(2,i+1,j,k) - vel(2,i-1,j,k))*delx
                 Omegaz(i,j,k) = Omegaz(i,j,k) -
     &           (vel(1,i,j+1,k) - vel(1,i,j-1,k))*dely

                 Omegax(i,j,k) = Omegax(i,j,k)*cboris*pie4i
                 Omegay(i,j,k) = Omegay(i,j,k)*cboris*pie4i
                 Omegaz(i,j,k) = Omegaz(i,j,k)*cboris*pie4i

              end do
           end do
        end do

! Compute |Omega|
        do k=kl_bnd+ndel*k3d+k3d,ku_bnd+ndel*k3d-k3d
         kk = nguard*k3d+k
         do j=jl_bnd+ndel*k2d+k2d,ju_bnd+ndel*k2d-k2d
          jj = nguard*k2d+j
          do i=il_bnd+ndel*k1d+k1d,iu_bnd+ndel*k1d-k1d
            ii = nguard*k1d+i
            bx = unk(6,i-ndel,j-ndel*k2d,k-ndel*k3d,lb)
            by = unk(7,i-ndel,j-ndel*k2d,k-ndel*k3d,lb)
            bz = unk(8,i-ndel,j-ndel*k2d,k-ndel*k3d,lb) 
            bxt= bx+b_sf(1,ii,jj,kk)
            byt= by+b_sf(2,ii,jj,kk)
            bzt= bz+b_sf(3,ii,jj,kk)
            bsq= bxt**2+byt**2+bzt**2+1.0e-28
            Omega_sq   = Omegax(i,j,k)**2 + Omegay(i,j,k)**2 + 
     $                   Omegaz(i,j,k)**2
            eerror(lb) = max(eerror(lb),sqrt(Omega_sq/bsq))
          end do
         end do
        end do
      endif
      end do
      endif
!       write(*,*)'leaving block loop'
       go to 56
      endif

      if(ir.eq.1)then  
!
! refine on density
!
! Set up the workspace array to store the current solution.
      if(lnblocks.gt.0) then
      do lb=1,lnblocks
      if(nodetype(lb).eq.1.or.nodetype(lb).eq.2) then
         do k = 1,nzb
            do j = 1,nyb
               do i = 1,nxb
                  work(i,j,k,lb,1) = unk(ir,i,j,k,lb)
!-----------------------------------------------------------------------
               end do
            end do
         end do
      endif
      
      end do
      endif

      elseif(ir.eq.2)then
!
! refine on pressure
!
      write(*,*)'refining on pressure:indices wrong: see ir.eq.1';stop
      if(lnblocks.gt.0) then
      do lb=1,lnblocks
      if(nodetype(lb).eq.1.or.nodetype(lb).eq.2) then
#ifdef BORIS_CORRECTION
        call strongfields_wb(mype,lb)
#endif
        do k=kl_bnd+ndel*k3d+k3d,ku_bnd+ndel*k3d-k3d
         kk = nguard*k3d+k
         do j=jl_bnd+ndel*k2d+k2d,ju_bnd+ndel*k2d-k2d
          jj = nguard*k2d+j
          do i=il_bnd+ndel*k1d+k1d,iu_bnd+ndel*k1d-k1d
           ii = nguard*k1d+i
           rr = unk(1,i-ndel,j-ndel*k2d,k-ndel*k3d,lb)
           px = unk(2,i-ndel,j-ndel*k2d,k-ndel*k3d,lb)
           py = unk(3,i-ndel,j-ndel*k2d,k-ndel*k3d,lb)
           pz = unk(4,i-ndel,j-ndel*k2d,k-ndel*k3d,lb)
           rr = unk(1,i-ndel,j-ndel*k2d,k-ndel*k3d,lb)
#ifndef MHD
           b2 = 0.0
           bx = 0.0
           by = 0.0
           bz = 0.0
#else
           bx = unk(6,i-ndel,j-ndel*k2d,k-ndel*k3d,lb)
           by = unk(7,i-ndel,j-ndel*k2d,k-ndel*k3d,lb)
           bz = unk(8,i-ndel,j-ndel*k2d,k-ndel*k3d,lb)
           b2 = bx**2+by**2+bz**2
#endif 
#ifdef BORIS_CORRECTION
           bxt= bx+b_sf(1,ii,jj,kk)
           byt= by+b_sf(2,ii,jj,kk)
           bzt= bz+b_sf(3,ii,jj,kk)
           bsq= bxt**2+byt**2+bzt**2
           rbi= 1.0/(rr + afact*bsq)
           bv = (px*bxt+py*byt+pz*bzt)*rbi
           vx = (px + afact*bv*bxt)*rbi
           vy = (py + afact*bv*byt)*rbi
           vz = (pz + afact*bv*bzt)*rbi
#else
           ri = 1.0/rr
           vx = px*ri
           vy = py*ri
           vz = pz*ri
#endif           
           et = unk(5,i-ndel,j-ndel*k2d,k-ndel*k3d,lb)
           vs = vx*vx+vy*vy+vz*vz
           ek = 0.5*rr*vs
           p  = max(gammam1*(et-ek-b2*pie8i),pmin)
           work(i,j,k,lb,1) = p
          end do
         end do
        end do
      endif
      end do
      endif
      write(*,*)'completed work on pressure'
 

      elseif(ir.eq.3)then
!
!-----refine on velocity magnitude
!
      if(lnblocks.gt.0) then
      do lb=1,lnblocks
      if(nodetype(lb).eq.1.or.nodetype(lb).eq.2) then
#ifdef BORIS_CORRECTION
        call strongfields_wb(mype,lb)
#endif
        do k=kl_bnd+ndel*k3d+k3d,ku_bnd+ndel*k3d-k3d
         kk = nguard*k3d+k
         do j=jl_bnd+ndel*k2d+k2d,ju_bnd+ndel*k2d-k2d
          jj = nguard*k2d+j
          do i=il_bnd+ndel*k1d+k1d,iu_bnd+ndel*k1d-k1d
            ii = nguard*k1d+i
           rr = unk(1,i-ndel,j-ndel*k2d,k-ndel*k3d,lb)
           px = unk(2,i-ndel,j-ndel*k2d,k-ndel*k3d,lb)
           py = unk(3,i-ndel,j-ndel*k2d,k-ndel*k3d,lb)
           pz = unk(4,i-ndel,j-ndel*k2d,k-ndel*k3d,lb)
           rr = unk(1,i-ndel,j-ndel*k2d,k-ndel*k3d,lb)
           bx = unk(6,i-ndel,j-ndel*k2d,k-ndel*k3d,lb)
           by = unk(7,i-ndel,j-ndel*k2d,k-ndel*k3d,lb)
           bz = unk(8,i-ndel,j-ndel*k2d,k-ndel*k3d,lb) 
#ifdef BORIS_CORRECTION
           bxt= bx+b_sf(1,ii,jj,kk)
           byt= by+b_sf(2,ii,jj,kk)
           bzt= bz+b_sf(3,ii,jj,kk)
           bsq= bxt**2+byt**2+bzt**2
           rbi= 1.0/(rr + afact*bsq)
           bv = (px*bxt+py*byt+pz*bzt)*rbi
           vx = (px + afact*bv*bxt)*rbi
           vy = (py + afact*bv*byt)*rbi
           vz = (pz + afact*bv*bzt)*rbi
#else
           ri = 1.0/rr
           vx = px*ri
           vy = py*ri
           vz = pz*ri
#endif           
           work(i,j,k,lb,1) = sqrt(vx*vx+vy*vy+vz*vz)
!-----------------------------------------------------------------------

          end do
         end do
        end do
      endif
      end do
      endif

      elseif(ir.eq.4)then
!
!-----refine on vorticity
!
      if(lnblocks.gt.0) then
      do lb=1,lnblocks
      if(nodetype(lb).eq.1.or.nodetype(lb).eq.2) then
        Omegax = 0.
        Omegay = 0.
        Omegaz = 0.
        delx   = bsize(1,lb)/float(nxb)
        dely   = bsize(2,lb)/float(nyb)
#if N_DIM == 3
        delz   = bsize(3,lb)/float(nzb)
#endif
        delx   = 1./(2.*delx)
        dely   = 1./(2.*dely)
#if N_DIM == 3 
        delz   = 1./(2.*delz)
#endif


#ifdef BORIS_CORRECTION
        call strongfields_wb(mype,lb)
#endif
        do k=kl_bnd+ndel*k3d+k3d,ku_bnd+ndel*k3d-k3d
         kk = nguard*k3d+k
         do j=jl_bnd+ndel*k2d+k2d,ju_bnd+ndel*k2d-k2d
          jj = nguard*k2d+j
          do i=il_bnd+ndel*k1d+k1d,iu_bnd+ndel*k1d-k1d
           ii = nguard*k1d+i
           rr = unk(1,i-ndel,j-ndel*k2d,k-ndel*k3d,lb)
           px = unk(2,i-ndel,j-ndel*k2d,k-ndel*k3d,lb)
           py = unk(3,i-ndel,j-ndel*k2d,k-ndel*k3d,lb)
           pz = unk(4,i-ndel,j-ndel*k2d,k-ndel*k3d,lb)
#ifdef BORIS_CORRECTION
           bx = unk(6,i-ndel,j-ndel*k2d,k-ndel*k3d,lb)
           by = unk(7,i-ndel,j-ndel*k2d,k-ndel*k3d,lb)
           bz = unk(8,i-ndel,j-ndel*k2d,k-ndel*k3d,lb) 
           bxt= bx+b_sf(1,ii,jj,kk)
           byt= by+b_sf(2,ii,jj,kk)
           bzt= bz+b_sf(3,ii,jj,kk)
           bsq= bxt**2+byt**2+bzt**2
           rbi= 1.0/(rr + afact*bsq)
           bv = (px*bxt+py*byt+pz*bzt)*rbi
           vx = (px + afact*bv*bxt)*rbi
           vy = (py + afact*bv*byt)*rbi
           vz = (pz + afact*bv*bzt)*rbi
#else
           ri = 1.0/rr
           vx = px*ri
           vy = py*ri
           vz = pz*ri
#endif            
           vel(1,i,j,k) = vx
           vel(2,i,j,k) = vy
           vel(3,i,j,k) = vz
          end do
         end do
        end do
         do k = 1+k3d,ku_bnd-k3d
            do j = 2,ju_bnd-1
               do i = 2,iu_bnd-1

! Curl V = Omega

                 Omegax(i,j,k) = 
     &           (vel(3,i,j+1,k) - vel(3,i,j-1,k))*dely
#if N_DIM == 3
                 Omegax(i,j,k) = Omegax(i,j,k) -
     &           (vel(2,i,j,k+1) - vel(2,i,j,k-1))*delz
                 Omegay(i,j,k) = 
     &           (vel(1,i,j,k+1) - vel(1,i,j,k-1))*delz
#endif
                 Omegay(i,j,k) = Omegay(i,j,k) -
     &           (vel(3,i+1,j,k) - vel(3,i-1,j,k))*delx

                 Omegaz(i,j,k) = 
     &           (vel(2,i+1,j,k) - vel(2,i-1,j,k))*delx
                 Omegaz(i,j,k) = Omegaz(i,j,k) -
     &           (vel(1,i,j+1,k) - vel(1,i,j-1,k))*dely


              end do
           end do
        end do

! Compute |Omega|

        do k=kl_bnd+ndel*k3d+k3d,ku_bnd+ndel*k3d-k3d
         do j=jl_bnd+ndel*k2d+k2d,ju_bnd+ndel*k2d-k2d
          do i=il_bnd+ndel*k1d+k1d,iu_bnd+ndel*k1d-k1d
            Omega_sq = Omegax(i,j,k)**2 + Omegay(i,j,k)**2 + 
     $                 Omegaz(i,j,k)**2
            work(i,j,k,lb,1) = sqrt(Omega_sq)
            if(Omega_sq.eq.0.0)then
c      write(*,*)'lb,delx,dely,delz=',lb,delx,dely,delz
            endif
             end do
          end do
        end do

      endif
      end do
      endif

      
!
      endif



!
! Fill the guard cell layers of the workspace array.
      iopt=2
      nlayers=nguard_work

! set the type of boundary condition to be applied to work
      ibc_work = 1


      if(ir.ne.ibc_work) then
         write(*,*) 'AMR_TEST_REFINEMENT IR .NE. IBC_WORK', ir, ibc_work
      endif


#ifndef ADVANCE_ALL_LEVELS
! A valid solution will be required on the parents of leaf blocks
! when refinement testing is done. See the comment before the call
! to amr_test_refinement.
        iempty = 0
!       iopt   = 1
c      write(*,*)'entering amr_restrict'

        call amr_restrict(mype,iopt,iempty)
c      write(*,*)'leaving amr_restrict'

#endif



#ifndef NO_PERMANENT_GUARDCELLS

c      write(*,*)'entering amr_guardcell'
        call amr_guardcell(mype,iopt,nlayers)
c      write(*,*)'leaving amr_guardcell'


#else /* NO_PERMANENT_GUARDCELLS */

      iopt       = 2
      lcc        = .true.
      lfc        = .false.
#ifdef MPI_USED
      lguard    = .true.
      lprolong  = .false.
      lflux     = .false.
      ledge     = .false.
      lrestrict = .false.
      lec        =.false.
      lnc        =.false.
      lfulltree=.false.
      tag_offset= 100
      call mpi_amr_comm_setup(mype,nprocs,
     .                        lguard,lprolong,lflux,ledge,lrestrict,
     .                        lfulltree,
     .                        iopt,lcc,lfc,lec,lnc,tag_offset)
#endif /*MPI_USED*/

#endif /* NO_PERMANENT_GUARDCELLS */


!
! Loop over all leaf blocks and all parents of leaf blocks

      if(lnblocks.gt.0) then
      do lb=1,lnblocks
      if(nodetype(lb).eq.1.or.nodetype(lb).eq.2) then
!
! User provided routine which returns an array error, which has some error
! measure computed for each grid cell, based on some computation on the 
! input array WORK.
       error = 0.0
c      write(*,*)'entering error_measure'
       call error_measure(mype,nprocs,iopt,lcc,lfc,error,epsil,lb)
c      write(*,*)'leaving error_measure'
c      write(*,*)'error=',error

       error_max = error 

! Does the error measure on this block anywhere exceed the limit which 
! should trigger refinement?

      if( lrefine(lb).lt.lrefine_max ) then
         if ( error_max . ge. ctore) then
            refine(lb) = .true.
c            write(*,*) mype,lb,'amr_test_refine: REFINING'
         endif
      endif


! Can we derefine this block?

      if( lrefine(lb).gt.lrefine_min .and. (.not.refine(lb)) ) then
        if ( error_max . lt. ctode) then
           derefine(lb) = .true.
c           write(*,*) mype,lb,'amr_test_refine: DEREFINING'
        endif
      endif



      endif
      end do                                   ! end of loop over blocks
      endif
      call shmem_barrier_all()


      enddo   ! end of loop over refinement parameters
c      stop

! Check if blocks are 'near' the earth

      if(lb.eq.21) write(*,*) 'end of test_ref : blk ',lb,
     .                 refine(lb),derefine(lb)


      return
 56   continue
c      write(*,*)'entering ir.eq.4.or.ir.eq.5 loop'
      if (ir.eq.4.or.ir.eq.5.or.ir.eq.6) then
         
         error_max = -1.e20
         do lb = 1,lnblocks
            
            if (nodetype(lb).eq.1.or.nodetype(lb).eq.2) then
             error_max = max(eerror(lb),error_max)
            else
             eerror(lb) = 0.
            end if
            
         end do
c      write(*,*)'a:leaving ir.eq.4.or.ir.eq.5 loop'
         
         call shmem_barrier_all()
         call COMM_REAL_MAX_TO_ALL(error_max,error_max)
         call shmem_barrier_all()
        
         do lb = 1,lnblocks
          eerror(lb) = eerror(lb)/error_max
         end do

      end if
      call shmem_barrier_all()
c      write(*,*)'b:leaving ir.eq.4.or.ir.eq.5 loop'

      do lb = 1,lnblocks
            
       if (nodetype(lb).eq.1) then
               
! fetch error of parent
        if (parent(1,lb).ge.1) then
      write(*,*)'before SHMEM_REAL_GET'
      write(*,*)'errort              =',errort
      write(*,*)'parent(1,lb)        =',parent(1,lb)
      write(*,*)'eerror(parent(1,lb))=',eerror(parent(1,lb))
      write(*,*)'parent(2,lb)        =',parent(2,lb)
      call SHMEM_REAL_GET (errort,eerror(parent(1,lb)),1,
     &      parent(2,lb))
      write(*,*)'after SHMEM_REAL_GET'
                  
          if (.not.refine(lb).and.eerror(lb).lt.ctode.and.
     &                 errort.lt.ctode) then
	  print *,' MARKING FOR DE-REFINEMENT;eerror(lb)= ',eerror(lb)
	  print *,' MARKING FOR DE-REFINEMENT;ctode= ',ctode
          derefine(lb) = .TRUE.
          end if
         end if
               
! test for refinement

        if (eerror(lb).gt.ctore) then
!	   print *,' MARKING FOR REFINEMENT ',eerror(lb)
           derefine(lb) = .FALSE.
           refine(lb)   = .TRUE.
           nref         = nref + 1
         end if
               
         if (lrefine(lb).ge.lrefine_max) 
     &      refine(lb) = .FALSE.
         end if

c         if(bflags(1,lb).eq.1)then
c           derefine(lb) = .FALSE.
c           refine(lb) = .FALSE.
c         endif    

      end do
c      write(*,*)'leaving ir.eq.4.or.ir.eq.5 loop'


      return
      end subroutine amr_test_refinement_dan
