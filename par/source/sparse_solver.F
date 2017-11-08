#include "paramesh_preprocessor.fh"
#ifdef CLEAN_DIVB
!     Michael L. Rilee, December 2002, *clean_divb*
!        Support for projecting field onto divergenceless field

! Reads file 'sparskit.ini' to initialize solver.

       module sparse_solver

       use constants
       use io, only : iu_log

       implicit none

       integer, save :: nmax, nzmax, maxits,lwk
       integer, save, allocatable, dimension(:) :: ia,ja,jau,ju,iw
       integer, save ::  ipar(16),nx,ny,nz,i,lfil,nwk,nrow,ierr
       real(kind=kind(1.0d0)), save, allocatable, dimension(:) ::  
     &   a,sol,rhs,au,wk, xran, al
       real(kind=kind(1.0d0)), save ::   fpar(16)
       real(kind=kind(1.0d0)), save ::   tol

       character(len=32) :: solver_name = 'gmres'

       logical, save :: sparse_solver_verbosep
       private :: ldbg
       logical, save :: ldbg = .true.

       external cg,bcg,dbcg,bcgstab,tfqmr,gmres,fgmres,dqgmres
       external cgnr, fom, ilut

       contains
  
       subroutine init_sparse_solver0(nmax_,nzmax_,maxits_)
        implicit none
        integer :: nmax_,nzmax_,maxits_
        nmax=nmax_; nzmax=nzmax_; maxits=maxits_
        lwk=nmax*40
       end subroutine init_sparse_solver0

       subroutine init_sparse_solver1
        call init_sparse_solver0(5000,100000,60)
       end subroutine init_sparse_solver1

       subroutine init_sparse_solver
        implicit none
        integer :: itmp
        real(kind=kind(1.0d0)) :: ftmp
    	character(len=32) :: ctmp
    	character(len=80) :: l80tmp
    	logical :: ltmp
    	open(unit=83,file='sparskit.ini',status='old')
    	read(83,*)l80tmp ! header
    	read(83,*)l80tmp; read(83,*)itmp ! preconditioning
    	ipar(2)=itmp
    	read(83,*)l80tmp; read(83,*)itmp ! comment and data
    	ipar(3)=itmp
    	read(83,*)l80tmp; read(83,*)itmp ! comment and data
    	ipar(4)=itmp
    	read(83,*)l80tmp; read(83,*)itmp ! comment and data
    	ipar(5)=itmp
    	read(83,*)l80tmp; read(83,*)itmp ! comment and data
    	ipar(6)=itmp
    	read(83,*)l80tmp; read(83,*)ftmp ! comment and data
    	fpar(1)=ftmp
    	read(83,*)l80tmp; read(83,*)ftmp ! comment and data
    	fpar(2)=ftmp
    	read(83,*)l80tmp; read(83,*)ctmp ! solver to use
    	solver_name=ctmp
    	read(83,*)l80tmp ! header
    	read(83,*)l80tmp ! header
    	read(83,*)l80tmp; read(83,*)ftmp ! comment and data
    	tol=ftmp
    	read(83,*)l80tmp; read(83,*)itmp ! comment and data
    	lfil=itmp
    	read(83,*)l80tmp ! header
    	read(83,*)l80tmp; read(83,*)itmp ! comment and data
    	nmax=itmp
    	lwk=nmax*40
    	read(83,*)l80tmp; read(83,*)itmp ! comment and data
    	nzmax=itmp
    	nwk=nzmax
    	read(83,*)l80tmp; read(83,*)itmp ! comment and data
    	maxits=itmp
    	read(83,*)l80tmp; read(83,*)ltmp  ! read verbosity control flag
    	sparse_solver_verbosep=ltmp
    	close(83)
    	if(ipar(4).lt.0)ipar(4) = lwk
    	if(ipar(6).lt.0)ipar(6) = maxits
    	!print *,'sn: ',solver_name,ipar,fpar
       end subroutine init_sparse_solver
  
       subroutine allocate_sparse_solver
        implicit none
    	allocate( 
     &    ia(nmax),ja(nzmax),jau(nzmax),ju(nzmax),iw(nmax*3))
        allocate( 
     &    a(nzmax),sol(nmax),rhs(nmax),au(nzmax),wk(nmax*40), &
     &    xran(nmax), al(nmax))
       end subroutine allocate_sparse_solver

       subroutine deallocate_sparse_solver
        deallocate(ia,ja,jau,ju,iw)
        deallocate(a,sol,rhs,au,wk,xran,al)
       end subroutine deallocate_sparse_solver

       subroutine init_run_solver( 
     &   n1,n2,n3,rhs_array,lfil_,tol_,ipar_,fpar_,solver_name_)
        implicit none
        integer :: n1, n2, n3, lfil_, ipar_(:)
        real(kind=kind(1.0d0)) :: rhs_array(:,:,:),tol_, fpar_(:)
        integer :: mode
        character(len=*) :: solver_name_
        ! The grid 
        nx = n1; ny = n2; nz = n3;
        ! Derived vals.
        nrow=nx * ny * nz
        ! Preconditioner (ILUT)
        lfil = lfil_
        tol = tol_
        nwk = nzmax
        ! Solver
        ipar(:) = ipar_(:)
        fpar(:) = fpar_(:)
        if(ipar(4).lt.0)ipar(4) = lwk
        if(ipar(6).lt.0)ipar(6) = maxits
        solver_name = solver_name_
       end subroutine init_run_solver

       subroutine init_run_solver_brief( 
     &   n1,n2,n3,rhs_array,firstp)
        implicit none
        integer :: n1, n2, n3
        real(kind=kind(1.0d0)) :: rhs_array(:,:,:)
        logical :: firstp
        ! The grid 
        nx = n1; ny = n2; nz = n3;
        ! Derived vals.
        nrow=nx * ny * nz
        ! Solver
        ipar(1)=0
       end subroutine init_run_solver_brief

       subroutine runrc( 
     &   n,rhs,sol,ipar,fpar,wk,guess,a,ja,ia, 
     &   au,jau,ju,solver,verbosep,iou)
        ! cf. runrc in SPARSPAK's "itaux.f".
        implicit none
        integer n,ipar(16),ia(n+1),ja(:),ju(:),jau(:),iou
        real(kind=kind(1.0d0)) :: 
     &    fpar(16),rhs(n),sol(n),guess(n),wk(:),a(:),au(:)
        logical :: verbosep
    
        ! local variables
        integer i,its
        real(kind=kind(1.0d0)) :: res, dnrm2
        external :: dnrm2
        save :: its,res

        ipar(1) = 0
    
        if(ipar(2)>2)then
         write (iou,*) 
     &   'sparse_solver::runrc:: 
     & Cannot perform both left and right precond.'
         stop
        end if
    
        its=0; res=0.0d0
        do i = 1,n; sol(i)=guess(i)
        end do

        loop_solver: do 
        call solver(n,rhs,sol,ipar,fpar,wk)
        if (ipar(7).ne.its) then
          if (verbosep)
     & write (iou, *) 'runrc:i7,f5: ',ipar(7), real(fpar(6))
           its = ipar(7)
          endif
          res = fpar(5)

          select case (ipar(1))
          case (1)
           call amux(n, wk(ipar(8):), wk(ipar(9):), a, ja, ia)
          case (2)
           call amux(n, wk(ipar(8):), wk(ipar(9):), a, ja, ia)
          case (3,5)
           call lusol(n,wk(ipar(8):),wk(ipar(9):),au,jau,ju)
          case (4,6)
           call lutsol(n,wk(ipar(8):),wk(ipar(9):),au,jau,ju)
          case (:0)
           if(verbosep)then
              select case (ipar(1))
              case (0)
                 write(iou,*)
     & 'runrc: Iterative solver satisfied convergence test.'
              case (-1)  
                 write(iou,*)
     & 'runrc: Iterative solver exceeded iteration limit.'
              case (-2)  
                 write(iou,*)'runrc: Iterative solver workspace too small.'
                 write(iou,*)
     & 'runrc:    Workspace should have at least ',ipar(4),
     & ' elements.'
              case (-3)
                 write(iou,*)'runrc: Iterarive solver facing breakdown.'
              case default
                write(iou,*)
     & 'runrc: Iterative solver terminated, code= ',ipar(1)
              end select
           end if ! verbosep
           exit loop_solver
          end select
         end do loop_solver

         if(verbosep)then
          write (iou, *) 'runrc:i7,f6: ',ipar(7), real(fpar(6))
          write (iou, *) 'runrc:# return code =', ipar(1), 
     &                   ' convergence rate =', fpar(7)
         end if

         if(ldbg)print *,'runrc:1000'
         call amux(n,sol(:),wk,a,ja,ia)
         if(ldbg)print *,'runrc:1001'
         wk(1:n) = sol(1:n) - 1.0d0
         if(ldbg)print *,'runrc:1002'
         wk(1:n) = wk(1:n) - rhs(1:n)
         if(ldbg)print *,'runrc:1003'
    
         if(verbosep)then
          write(iou,*)
     & 'runrc:# the actual residual norm is',dnrm2(n,wk(1:),1)
          write(iou,*)
     & 'runrc:# the error norm is',dnrm2(n,wk(1+n:),1)
         end if

        end subroutine runrc

        subroutine run_solver(ret_sol,verbose)
         implicit none
         real(kind=kind(1.0d0)) :: ret_sol(:,:,:)
         integer :: ix, iy, iz, i, jx, jy, jz
         logical :: verbose

         if(ldbg)print *,'rs:1000'

         ! Preconditioning step.
         call ilut (nrow,a,ja,ia,lfil,tol,au,jau,ju,nwk, 
     & wk,iw,ierr)

         if(ldbg)print *,'rs:2000'
         if(ldbg)print *,'rs:2001:solver_name: ',solver_name

         ! Generate a linear solution
         select case (solver_name)
         !! GMRES and variants seem to be preferred for elliptic PDE work.
         case ('gmres')
          call runrc(nrow,rhs,sol,ipar,fpar,wk,xran,a,ja,ia,au,jau,ju, 
     & gmres,verbose,iu_log)
         case ('fgmres')
          call runrc(nrow,rhs,sol,ipar,fpar,wk,xran,a,ja,ia,au,jau,ju, 
     & fgmres,verbose,iu_log) 
         case ('dqgmres')
          call runrc(nrow,rhs,sol,ipar,fpar,wk,xran,a,ja,ia,au,jau,ju, 
     & dqgmres,verbose,iu_log) 
         !! CG-based iterations for completeness and experimentation.
         case ('cg')
          print *,'cg-100'
          call runrc(nrow,rhs,sol,ipar,fpar,wk,xran,a,ja,ia,au,jau,ju, 
     & cg,verbose,iu_log)
          print *,'cg-101'
         case ('bcg')
          call runrc(nrow,rhs,sol,ipar,fpar,wk,xran,a,ja,ia,au,jau,ju, 
     & bcg,verbose,iu_log)
         case ('dbcg')
          call runrc(nrow,rhs,sol,ipar,fpar,wk,xran,a,ja,ia,au,jau,ju, 
     & dbcg,verbose,iu_log)
         case ('bcgstab')
          call runrc(nrow,rhs,sol,ipar,fpar,wk,xran,a,ja,ia,au,jau,ju, 
     & bcgstab,verbose,iu_log)
         case ('tfqmr')
          call runrc(nrow,rhs,sol,ipar,fpar,wk,xran,a,ja,ia,au,jau,ju, 
     & tfqmr,verbose,iu_log)
         case default
          print *,'sparse_solver::run_solver--error: solver_name not found'
          stop
         end select

         if(ldbg)print *,'rs:3000'
    
         i=0
         do iz=1,nz
          do iy=1,ny
           do ix=1,nx
             i=i+1
             ret_sol(ix,iy,iz) = sol(i) ! is this correct? 
             ! can we recode without the copy?
           end do
          end do
         end do

        end subroutine run_solver

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  subroutine sparse_test1(ret_sol,n1,n2,n3,rhs_array,maxit_,genmatrix)
!!    ! An example of how to use these routines to find a solution.
!!    ! Not tested recently...
!!    real(kind=kind(1.0d0)) :: ret_sol(:,:,:),rhs_array(:,:,:)
!!    integer :: n1,n2,n3,maxit_
!!
!!    integer :: lfil_st, ipar_st(16)
!!    real(kind=kind(1.0d0)) :: tol_st, fpar_st(16)
!!    logical :: verbosep
!!    
!!    external genmatrix
!!
!!
!!    call init_spars(5000,100000,maxit_)
!!    call allocate_spars
!!
!!    lfil_st=8           ! Number of nonzero row elements to keep
!!    tol_st=0.1          ! Tolerance for dropping elements from calc.
!!    ipar_st(1)=0        ! first run, set to zero (redundant)
!!    ipar_st(2)=2        ! right preconditioning
!!    ipar_st(3)=1        
!!    ipar_st(4)=-1       ! automatically set workspace size
!!    ipar_st(5)=8        ! Krylov subspace dim: higher is better
!!    ipar_st(6)=-1       ! automatically set maximum iterations
!!    fpar_st(1)=zero     ! relative tolerance (stopping criteria)
!!    fpar_st(2)=1.0d-15  ! absolute tolerance
!!    call init_run_solver( &
!!         & n1,n2,n3,rhs_array,lfil_st,tol_st,ipar_st,fpar_st,'gmres')
!!
!!    call genmatrix(nx,ny,nz,al,2,nrow,a,ja,ia,ju,rhs,rhs_array) ! set A
!!    call genmatrix(nx,ny,nz,al,3,nrow,a,ja,ia,ju,rhs,rhs_array) ! set RHS
!!
!!    verbosep = sparse_solver_verbosep ! print iteration status at each iter.
!!    call run_solver(ret_sol,verbosep)
!!    call deallocate_spars
!!  end subroutine sparse_test1
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       end module sparse_solver
#else
       subroutine sparse_solver
       end subroutine sparse_solver
#endif
