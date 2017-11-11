#include "paramesh_preprocessor.fh"
#ifdef CLEAN_DIVB
!     Michael L. Rilee, December 2002, *clean_divb*
!        Support for projecting field onto divergenceless field

       module clean_divb

       use constants

       use sparse_solver, only : 
     & init_sparse_solver, allocate_sparse_solver, 
     & init_run_solver, init_run_solver_brief,  
     & run_solver, deallocate_sparse_solver, 
     & al,a,ja,ia,rhs,ju,nrow, 
     & sparse_solver_verbosep, 
     & ipar,fpar
       implicit none

       real(kind=kind(1.0d0)), save, allocatable :: 
     & divb(:,:,:),phi(:,:,:)

! Set acceptable tolerance on original nett flux (ie B x Area)
! through sub-grid block surface
       real(kind=kind(1.0d0)), save :: b_da_tolerance = 1.0e-15

  
       integer, save :: n1,n2,n3

       ! use to toggle calls to these routines
       logical, save :: clean_divb_testflag=.false.

! Set the factor by which the cleanup will improve the value of divB
! beyond its value before the cleanup step.
       real(kind=kind(1.0d0)), save :: 
     & divb_improvement_tolerance = 1.0d-6
       integer, save :: clean_divb_test_failures
  
       contains
  
       subroutine start_clean_field(n1_,n2_,n3_)
        implicit none
        integer :: n1_,n2_,n3_
        n1=n1_; n2=n2_; n3=n3_
        allocate(divb(n1,n2,n3),phi(n1,n2,n3))
        call init_sparse_solver
        call allocate_sparse_solver
       end subroutine start_clean_field

       subroutine stop_clean_field
        implicit none
        call deallocate_sparse_solver
        deallocate(divb,phi)
       end subroutine stop_clean_field

       subroutine clean_divb_test_report(nerrors)
        implicit none
        integer, intent(inout) :: nerrors
        if(clean_divb_testflag)then
         print *,'clean_divb failed ',clean_divb_test_failures, 'times.'
         nerrors=nerrors+clean_divb_test_failures
        end if
       end subroutine clean_divb_test_report

       subroutine clean_field(f1,f2,f3, 
     &                 x1,x2,x3, 
     &                 xc1,xc2,xc3, 
     &                 calc_matrixp,status)
       implcit none
       real(kind=kind(1.0d0)), dimension(:,:,:) :: f1,f2,f3
       real(kind=kind(1.0d0)), dimension(:) :: x1,x2,x3
       real(kind=kind(1.0d0)), dimension(:) :: xc1,xc2,xc3
       real(kind=kind(1.0d0)) :: dx1i,dx2i,dx3i,da,b_da, 
     & f1save0,f1save1,b_da_save,max_divb0_save,max_divb0,max_divb1
       integer :: i1,i2,i3
       logical :: calc_matrixp, status
       logical :: dbg = .true.
       logical :: dbg1 = .true.

       ! PM: look at firstp...
       ! MLR: Changed name to calc_matrixp.

       ! PM: #ifdef curvilinear
       !     print warning if this code is compiled & run with curvilinear on
       ! MLR:  Put 'include clean_divb.fh' in ``initialize.F''.

       if(dbg1)print *,'cf-1000'

       ! test execution only !
       if(clean_divb_testflag)then
! compute initial divergence for later testing of how much cleaning will
! have been accomplished.
       divb=0.0d0
       do i3=1,n3 
          dx3i=one/(x3(i3+1)-x3(i3))    ! coord difference in z direction
          do i2=1,n2 
             dx2i=one/(x2(i2+1)-x2(i2))
             do i1=1,n1
                dx1i=one/(x1(i1+1)-x1(i1))
                divb(i1,i2,i3) = &
                     &   (f1(i1+1,i2,i3) - f1(i1,i2,i3))*dx1i &
                     &+  (f2(i1,i2+1,i3) - f2(i1,i2,i3))*dx2i &
                     &+  (f3(i1,i2,i3+1) - f3(i1,i2,i3))*dx3i
             end do
          end do
       end do
       max_divb0_save=maxval(abs(divb))

       f1save0=f1(n1/2+1,n2/2+1,n3/2+1)
! introduce deliberate error for testing cleanup
       f1(n1/2+1,n2/2+1,n3/2+1)=f1save0+5.0d-11
       f1save1=f1(n1/2+1,n2/2+1,n3/2+1)
       end if

! verify that initial nett flux is zero before starting cleanup
       b_da=0.0d0
       do i3=1,n3
        do i2=1,n2
          da=(x2(i2+1)-x2(i2))*(x3(i3+1)-x3(i3))      ! face area
          do i1=1,1
             b_da=b_da-f1(i1,i2,i3)*da
          end do
         end do
       end do
       do i3=1,n3
        do i2=1,n2
          da=(x2(i2+1)-x2(i2))*(x3(i3+1)-x3(i3))
          do i1=n1+1,n1+1
             b_da=b_da+f1(i1,i2,i3)*da
          end do
         end do
       end do

       do i3=1,n3
        do i1=1,n1
          da=(x1(i1+1)-x1(i1))*(x3(i3+1)-x3(i3))
          do i2=1,1
             b_da=b_da-f2(i1,i2,i3)*da
          end do
        end do
       end do
       do i3=1,n3
        do i1=1,n1
          da=(x1(i1+1)-x1(i1))*(x3(i3+1)-x3(i3))
          do i2=n2+1,n2+1
             b_da=b_da+f2(i1,i2,i3)*da
          end do
        end do
       end do

       do i1=1,n1
        do i2=1,n2
          da=(x1(i1+1)-x1(i1))*(x2(i2+1)-x2(i2))
          do i3=1,1
             b_da=b_da-f3(i1,i2,i3)*da
          end do
        end do
       end do
       do i1=1,n1
        do i2=1,n2
          da=(x1(i1+1)-x1(i1))*(x2(i2+1)-x2(i2))
          do i3=n3+1,n3+1
             b_da=b_da+f3(i1,i2,i3)*da
          end do
        end do
       end do
       b_da_save=b_da
       if(dbg)print *,'cf-b_da,1: ',b_da_save
       if(abs(b_da)>b_da_tolerance)
     & print *,'cf-warning:divergence out of block detected! abs ' 
     & ,b_da,'>',b_da_tolerance

       ! compute divb. This array will be the source for the potential
       ! equation.
       divb=0.0d0
       do i3=1,n3 
        dx3i=one/(x3(i3+1)-x3(i3))
        do i2=1,n2 
          dx2i=one/(x2(i2+1)-x2(i2))
          do i1=1,n1
             dx1i=one/(x1(i1+1)-x1(i1))
             divb(i1,i2,i3) = &
                  &   (f1(i1+1,i2,i3) - f1(i1,i2,i3))*dx1i &
                  &+  (f2(i1,i2+1,i3) - f2(i1,i2,i3))*dx2i &
                  &+  (f3(i1,i2,i3+1) - f3(i1,i2,i3))*dx3i
          end do
        end do
       end do
       max_divb0=maxval(abs(divb))

! initialize the potential
       phi=0.0d0

       if(dbg1)print *,'cf-2000'

       call init_run_solver_brief(n1,n2,n3,divb,calc_matrixp)

       if(dbg1)print *,'cf-3000'

       al(1:6) = 0.0           !This array is obsolete
       if(calc_matrixp)then
        call genmatrix1(n1,n2,n3,al,2,nrow,a,ja,ia,ju,rhs,
     & xc1,xc2,xc3,divb,f1,f2,f3,.false.) ! set A, no n.f source.
       end if

       if(dbg1)print *,'cf-4000'

       !call dump(1,10,.true.,a,ja,ia,6)
       !call dump(60,70,.true.,a,ja,ia,6)
       call genmatrix1(n1,n2,n3,al,3,nrow,a,ja,ia,ju,rhs,
     & xc1,xc2,xc3,divb,f1,f2,f3,.false.) ! set RHS, no n.f source.

       if(dbg1)print *,'cf-5000'

       call run_solver(phi,sparse_solver_verbosep)

       if(dbg1)print *,'cf-6000'

       !  Are we approaching a solution? ( The vectors ipar and fpar
       !  are exchanged with SPARSKIT and are associated with monitoring
       ! the quality of the elliptic solve. Some are input from
       ! sparskit.ini)
       status=(ipar(1).eq.0).or. 
     & ((ipar(1).eq.-1).and.(fpar(7).gt.zero))
       ! If so...
       if(status)then
        do i3=1,n3
          do i2=1,n2
             do i1=2,n1
                dx1i=one/(xc1(i1)-xc1(i1-1))
                f1(i1,i2,i3)=f1(i1,i2,i3) &
                     & -(phi(i1,i2,i3)-phi(i1-1,i2,i3))*dx1i
             end do
          end do
        end do

        do i3=1,n3
          do i2=2,n2
             dx2i=one/(xc2(i2)-xc2(i2-1))
             do i1=1,n1
                f2(i1,i2,i3)=f2(i1,i2,i3) &
                     & -(phi(i1,i2,i3)-phi(i1,i2-1,i3))*dx2i
             end do
          end do
        end do

        do i3=2,n3
          dx3i=one/(xc3(i3)-xc3(i3-1))
          do i2=1,n2
             do i1=1,n1
                f3(i1,i2,i3)=f3(i1,i2,i3) &
                     & -(phi(i1,i2,i3)-phi(i1,i2,i3-1))*dx3i
             end do
          end do
        end do
        if(sparse_solver_verbosep)print *,'cf-success'
       else
        print *,'cf-failure:ipar: ',ipar
       end if
       ! We could actually see if the divb was reduced and if so, pass that back.

       ! recompute divb to verify that it has been improved
       divb=0.0d0
       do i3=1,n3 
        dx3i=one/(x3(i3+1)-x3(i3))
        do i2=1,n2 
          dx2i=one/(x2(i2+1)-x2(i2))
          do i1=1,n1
             dx1i=one/(x1(i1+1)-x1(i1))
             divb(i1,i2,i3) = &
                  &   (f1(i1+1,i2,i3) - f1(i1,i2,i3))*dx1i &
                  &+  (f2(i1,i2+1,i3) - f2(i1,i2,i3))*dx2i &
                  &+  (f3(i1,i2,i3+1) - f3(i1,i2,i3))*dx3i
          end do
        end do
       end do
       max_divb1=maxval(abs(divb))
       if(dbg)print *,'cf-max-divb,1:             ',max_divb0
       if(dbg)print *,'cf-max-divb,2:             ',max_divb1
       ! test !
!      if(clean_divb_testflag)then
         if(divb_improvement_tolerance*max_divb1>max_divb0)then
          clean_divb_test_failures=clean_divb_test_failures+1
          print *,'input divb,maxval(abs(divb)):             ',max_divb0
          print *,'correctd divb,maxval(abs(divb)):          ',max_divb1
         if(clean_divb_testflag)then
          print *,'original divb,maxval(abs(divb)):          ',
     & max_divb0_save
          print *,'original,f1(n1/2+1,n2/2+1,n3/2+1):        ',f1save0
          print *,'corrupted,f1(n1/2+1,n2/2+1,n3/2+1),delta: ',
     & f1save1,f1save1-f1save0
          print *,'corrected,f1(n1/2+1,n2/2+1,n3/2+1),delta: ', 
     & f1(n1/2+1,n2/2+1,n3/2+1),f1(n1/2+1,n2/2+1,n3/2+1)-f1save0
         end if
       end if
!    end if


       if(dbg1)print *,'cf-8000'

       end subroutine clean_field

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! FDIF stencil adapted to divergenceless finite volume situation.
  ! Equivalent to second-order finite difference Laplace problem with
  ! constant coefficients and nonzero right hand side.  The solution
  ! volume is a 3-rectangular region discretized in the 3 Cartesian
  ! dimensions as x1(1:nx), x2(1:ny), x3(1:nz).  FD-nodes correspond
  ! to cell centers of the finite volumes.  Boundaries correspond to
  ! the cell faces. 

  ! The normal derivatives on these boundaries are set to zero because
  ! the MHD codes we work with already guarantee zero net flux through
  ! block boundaries. This modifies the structure of the matrix as
  ! follows.

  ! The net result of the solve is the return of the scalar potential
  ! at cell centers within the block.

       function west1(i1,i2,i3,x1,x2,x3) result(sten)
       real(kind=kind(1.0d0)) :: x1(:),x2(:),x3(:),sten,tmpsten,hi
       integer i1,i2,i3
       if(i1.eq.1)then
          ! na
          print *,'west-stop'; stop
       else if(i1.eq.n1)then
          hi=one/(x1(i1)-x1(i1-1))
          tmpsten=hi/(x1(i1)-x1(i1-1))
       else
          hi=two/(x1(i1+1)-x1(i1-1))
          tmpsten=hi/(x1(i1)-x1(i1-1))
       end if
       sten=tmpsten
       end function west1

       function south1(i1,i2,i3,x1,x2,x3) result(sten)
       real(kind=kind(1.0d0)) :: x1(:),x2(:),x3(:),sten,tmpsten,hi
       integer :: i1,i2,i3
       if(i2.eq.1)then
          ! na
          print *,'south-stop'; stop
       else if(i2.eq.n2)then
          hi=one/(x2(i2)-x2(i2-1))
          tmpsten=hi/(x2(i2)-x2(i2-1))
       else
          hi=two/(x2(i2+1)-x2(i2-1))
          tmpsten=hi/(x2(i2)-x2(i2-1))
       end if
       sten=tmpsten
       end function south1

       function front1(i1,i2,i3,x1,x2,x3) result(sten)
       real(kind=kind(1.0d0)) :: x1(:),x2(:),x3(:),sten,tmpsten,hi
       integer :: i1,i2,i3
       if(i3.eq.1)then
          ! na
          print *,'front-stop'; stop
       else if(i3.eq.n3)then
          hi=one/(x3(i3)-x3(i3-1))
          tmpsten=hi/(x3(i3)-x3(i3-1))
       else
          hi=two/(x3(i3+1)-x3(i3-1))
          tmpsten=hi/(x3(i3)-x3(i3-1))
       end if
       sten=tmpsten
       end function front1

       function center1(i1,i2,i3,x1,x2,x3) result(sten)
       real(kind=kind(1.0d0)) :: x1(:),x2(:),x3(:),
     &  sten,tmpsten,hi,h2pi,h2mi
       integer :: i1,i2,i3

       tmpsten=zero

       if(n1>1)then
       if(i1.eq.1)then
          hi=one/(x1(i1+1)-x1(i1))
          h2pi=hi/(x1(i1+1)-x1(i1))
          tmpsten=tmpsten-h2pi
       else if(i1.eq.n1)then
          hi=one/(x1(i1)-x1(i1-1))
          h2mi=hi/(x1(i1)-x1(i1-1))
          tmpsten=tmpsten-h2mi
       else
          hi=two/(x1(i1+1)-x1(i1-1))
          h2pi=hi/(x1(i1+1)-x1(i1))
          h2mi=hi/(x1(i1)-x1(i1-1))
          tmpsten=tmpsten-(h2pi+h2mi)
       end if
       end if

       if(n2>1)then
       if(i2.eq.1)then
          hi=one/(x2(i2+1)-x2(i2))
          h2pi=hi/(x2(i2+1)-x2(i2))
          tmpsten=tmpsten-h2pi
       else if(i2.eq.n2)then
          hi=one/(x2(i2)-x2(i2-1))
          h2mi=hi/(x2(i2)-x2(i2-1))
          tmpsten=tmpsten-h2mi
       else
          hi=two/(x2(i2+1)-x2(i2-1))
          h2pi=hi/(x2(i2+1)-x2(i2))
          h2mi=hi/(x2(i2)-x2(i2-1))
          tmpsten=tmpsten-(h2pi+h2mi)
       end if
       end if

       if(n3>1)then
       if(i3.eq.1)then
          hi=one/(x3(i3+1)-x3(i3))
          h2pi=hi/(x3(i3+1)-x3(i3))
          tmpsten=tmpsten-h2pi
       else if(i3.eq.n3)then
          hi=one/(x3(i3)-x3(i3-1))
          h2mi=hi/(x3(i3)-x3(i3-1))
          tmpsten=tmpsten-h2mi
       else
          hi=two/(x3(i3+1)-x3(i3-1))
          h2pi=hi/(x3(i3+1)-x3(i3))
          h2mi=hi/(x3(i3)-x3(i3-1))
          tmpsten=tmpsten-(h2pi+h2mi)
       end if
       end if

       sten=tmpsten
       end function center1

       function east1(i1,i2,i3,x1,x2,x3) result(sten)
       real(kind=kind(1.0d0)) :: x1(:),x2(:),x3(:),sten,tmpsten,hi
       integer :: i1,i2,i3
       if(i1.eq.1)then
          hi=one/(x1(i1+1)-x1(i1))
          tmpsten=hi/(x1(i1+1)-x1(i1))
       else if(i1.eq.n1)then
          ! na
          print *,'east-stop'; stop
       else
          hi=two/(x1(i1+1)-x1(i1-1))
          tmpsten=hi/(x1(i1+1)-x1(i1))
       end if
       sten=tmpsten
       end function east1

       function north1(i1,i2,i3,x1,x2,x3) result(sten)
       real(kind=kind(1.0d0)) :: x1(:),x2(:),x3(:),sten,tmpsten,hi
       integer :: i1,i2,i3
       if(i2.eq.1)then
          hi=one/(x2(i2+1)-x2(i2))
          tmpsten=hi/(x2(i2+1)-x2(i2))
       else if(i2.eq.n2)then
          ! na
          print *,'north-stop'; stop
       else
          hi=two/(x2(i2+1)-x2(i2-1))
          tmpsten=hi/(x2(i2+1)-x2(i2))
       end if
       sten=tmpsten
       end function north1

       function back1(i1,i2,i3,x1,x2,x3) result(sten)
       real(kind=kind(1.0d0)) :: x1(:),x2(:),x3(:),sten,tmpsten,hi
       integer :: i1,i2,i3
       if(i3.eq.1)then
          hi=one/(x3(i3+1)-x3(i3))
          tmpsten=hi/(x3(i3+1)-x3(i3))
       else if(i3.eq.n3)then
          ! na
          print *,'back-stop'; stop
       else
          hi=two/(x3(i3+1)-x3(i3-1))
          tmpsten=hi/(x3(i3+1)-x3(i3))
       end if
       sten=tmpsten
       end function back1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! cf gen57pt in genmat.f from SPARSKIT2
       subroutine genmatrix1(nx,ny,nz,al,mode,n,a,ja,ia,iau,rhs, 
     & x1,x2,x3,rhs_array,f1,f2,f3,ndotf_bc_rhs)
       ! Generates a sparse matrix based on finite difference stencil.
       ! Based on "genmat.f" from SPARSKIT2.
       !
       implicit none
       integer :: ja(:),ia(:),iau(:), nx, ny, nz, mode, n
       real(kind=kind(1.0d0)) ::  a(:), rhs(:), al(6),rhs_array(:,:,:)
       real(kind=kind(1.0d0)) :: x1(:), x2(:), x3(:)
       real(kind=kind(1.0d0)) :: f1(:,:,:), f2(:,:,:), f3(:,:,:)

       integer :: ix,iy,iz,kx,ky,kz
       integer :: iedge, node

       real(kind=kind(1.0d0)) :: h, x, y, z, hscale, hscale2, ftmp

       logical :: ndotf_bc_rhs ! =.false. ! .true.

       hscale=x1(2)-x1(1); hscale2=hscale*hscale

!       if (nx.le.1) then
       if (nx.le.0) then
        print *,'genmatrix1-stop: nx',nx; stop
       end if

       kx=1; ky=nx; kz=nx*ny  ! strides through the array for index-deltas
       iedge=1                ! essentially the row
       node=1                 ! The element of the state-vector to consider (row)

       if(mode==1)then ! both A and rhs
       iz_loop1: do iz=1,nz
        iy_loop1: do iy=1,ny
          ix_loop1: do ix=1,nx
             ia(node) = iedge ! the "i" or row index
             ! begin add/enter array elements and increment iedge after each
             ! west
             if(ix>1)then
                ja(iedge)=node-kx
                a(iedge) = hscale2*west1(ix,iy,iz,x1,x2,x3)
                iedge=iedge+1
             end if
             ! south
             if(iy>1)then
                ja(iedge)=node-ky
                a(iedge) = hscale2*south1(ix,iy,iz,x1,x2,x3)
                iedge=iedge+1
             end if
             ! front
             if(iz>1)then
                ja(iedge)=node-kz
                a(iedge) = hscale2*front1(ix,iy,iz,x1,x2,x3)
                iedge=iedge+1
             end if
             ! center
             ja(iedge)=node
             iau(node)=iedge
             a(iedge) = hscale2*center1(ix,iy,iz,x1,x2,x3)
             iedge=iedge+1
             ! east
             if(ix<nx)then
                ja(iedge)=node+kx
                a(iedge) = hscale2*east1(ix,iy,iz,x1,x2,x3)
                iedge=iedge+1
             end if
             ! north
             if(iy<ny)then
                ja(iedge)=node+ky
                a(iedge) = hscale2*north1(ix,iy,iz,x1,x2,x3)
                iedge=iedge+1
             end if
             ! back
             if(iz<nz)then
                ja(iedge)=node+kz
                a(iedge) = hscale2*back1(ix,iy,iz,x1,x2,x3)
                iedge=iedge+1
             end if
             ! end add/enter array elements
             ! set rhs(node) if necessary
             ! Note that the following lines setting rhs implicitly define
             ! the connection between {x,y,z} and {ix,iy,iz}.  It
             ! should be consistent with the definition of h above.
! Was...
!             rhs(node)=h*h*hfun(x,y,z)
! And then...
!             rhs(node)=h*h*rhs_array(ix,iy,iz) 
! But now...
             rhs(node)=hscale2*rhs_array(ix,iy,iz) 
! Add in boundary terms to get d/dn phi = n dot B.
! Q:  For what problems is this a valid boundary condition?
! A:  Tokomaks.  Situations with fieldless conductors during
!     the short term.
             if(ndotf_bc_rhs)then
             if(ix==1)then
                rhs(node)=rhs(node)+hscale*f1(ix,iy,iz)
             else if(ix==nx)then
                rhs(node)=rhs(node)-hscale*f1(ix+1,iy,iz)
             end if
             if(iy==1)then
                rhs(node)=rhs(node)+hscale*f2(ix,iy,iz)
             else if(iy==ny)then
                rhs(node)=rhs(node)-hscale*f2(ix,iy+1,iz)
             end if
             if(iz==1)then
                rhs(node)=rhs(node)+hscale*f3(ix,iy,iz)
             else if(iz==nz)then
                rhs(node)=rhs(node)-hscale*f3(ix,iy,iz+1)
             end if
             end if
             ! done !
             node=node+1 ! increment 
          end do ix_loop1
         end do iy_loop1
        end do iz_loop1
        ia(node)=iedge
       else if(mode==2)then ! both A only
        iz_loop2: do iz=1,nz
         iy_loop2: do iy=1,ny
          ix_loop2: do ix=1,nx
             ia(node) = iedge ! the "i" or row index
             ! begin add/enter array elements and increment iedge after each
             ! west
             if(ix>1)then
                ja(iedge)=node-kx
                a(iedge) = hscale2*west1(ix,iy,iz,x1,x2,x3)
                iedge=iedge+1
             end if
             ! south
             if(iy>1)then
                ja(iedge)=node-ky
                a(iedge) = hscale2*south1(ix,iy,iz,x1,x2,x3)
                iedge=iedge+1
             end if
             ! front
             if(iz>1)then
                ja(iedge)=node-kz
                a(iedge) = hscale2*front1(ix,iy,iz,x1,x2,x3)
                iedge=iedge+1
             end if
             ! center
             ja(iedge)=node
             iau(node)=iedge
             a(iedge) = hscale2*center1(ix,iy,iz,x1,x2,x3)
             iedge=iedge+1
             ! east
             if(ix<nx)then
                ja(iedge)=node+kx
                a(iedge) = hscale2*east1(ix,iy,iz,x1,x2,x3)
                iedge=iedge+1
             end if
             ! north
             if(iy<ny)then
                ja(iedge)=node+ky
                a(iedge) = hscale2*north1(ix,iy,iz,x1,x2,x3)
                iedge=iedge+1
             end if
             ! back
             if(iz<nz)then
                ja(iedge)=node+kz
                a(iedge) = hscale2*back1(ix,iy,iz,x1,x2,x3)
                iedge=iedge+1
             end if
             ! end add/enter array elements
             ! set rhs(node) if necessary
             ! Note that the following lines setting rhs implicitly define
             ! the connection between {x,y,z} and {ix,iy,iz}.  It
             ! should be consistent with the definition of h above.
! Do not set...
!             rhs(node)=h*h*rhs_array(ix,iy,iz) 
             ! done !
             node=node+1 ! increment 
          end do ix_loop2
         end do iy_loop2
        end do iz_loop2
        ia(node)=iedge
       else if(mode==3)then ! rhs only
       iz_loop3: do iz=1,nz
        iy_loop3: do iy=1,ny
          ix_loop3: do ix=1,nx
!             rhs(node)=h*h*rhs_array(ix,iy,iz) 
             rhs(node)=hscale2*rhs_array(ix,iy,iz) 
             if(ndotf_bc_rhs)then
! Add in boundary terms to get d/dn phi = n dot B.
! Q:  For what problems is this a valid boundary condition?
! A:  Tokomaks.  Situations with fieldless conductors during
!     the short term.
             if(ix==1)then
                rhs(node)=rhs(node)+hscale*f1(ix,iy,iz)
             else if(ix==nx)then
                rhs(node)=rhs(node)-hscale*f1(ix+1,iy,iz)
             end if
             if(iy==1)then
                rhs(node)=rhs(node)+hscale*f2(ix,iy,iz)
             else if(iy==ny)then
                rhs(node)=rhs(node)-hscale*f2(ix,iy+1,iz)
             end if
             if(iz==1)then
                rhs(node)=rhs(node)+hscale*f3(ix,iy,iz)
             else if(iz==nz)then
                rhs(node)=rhs(node)-hscale*f3(ix,iy,iz+1)
             end if
             end if
             node=node+1 ! increment 
          end do ix_loop3
         end do iy_loop3
        end do iz_loop3
!        ia(node)=iedge
       else
        print *,'clean_divb::genmatrix-no-such-mode! stopping...'
        stop
       end if
      end subroutine genmatrix1

      end module clean_divb
#else
      subroutine clean_divb
      end subroutine clean_divb
#endif
