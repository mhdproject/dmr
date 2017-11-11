	subroutine pim_1d(v1,vl,vr,nxl,nxh)
        use paramesh_dimensions
        use physicaldata
	use athena_specific_parameters
	implicit none
        integer :: i,nxl,nxh,l      
        real, dimension(nem1,nmax), intent(in   ) :: v1
        real, dimension(nem1,nmax), intent(  out) :: vr,vl
!----------------------------------------------------------------
        real, dimension(nem1,nmax)                :: vf
        real, dimension(nem1,nmax)                :: f0
        real, dimension(nem1,nmax)                :: f
        real, dimension(nem1,nmax)                :: a2,a3
        real, dimension(nem1,nmax)                :: xl,xr
        real :: f1,f2,f3,f4,a1
        real, parameter :: alpha  = 0.7
        real, parameter :: pdmb   = (1.0-alpha)/alpha
        real,parameter  :: cf1    = 37.0/60.0
        real,parameter  :: cf2    =  2.0/15.0
        real,parameter  :: cf3    =  1.0/60.0
!
!------compute high order face quantities vf at face i from the cell averaged quantities v1
!
!
!-----compute the high-order primitives at the faces vf ala Zalesak(IMACS 1981)
!
        do i=nxl,nxh
         vf(:,i) = cf1*(v1(:,i-1)+v1(:,i  ))
     &            -cf2*(v1(:,i-2)+v1(:,i+1))
     &            +cf3*(v1(:,i-3)+v1(:,i+2))
        enddo
!
!-----next clip vf with respect to their old values at time level n-1
!

        do i = nxl,nxh
         vf(:,i)=max(min( v1(:,i),v1(:,i-1) ),
     &               min( vf(:,i),
     &               max( v1(:,i),v1(:,i-1) )))
        enddo
!
!-----apply pdm limitor and set up left and right states for riemann solver
!
        a1   = 0.

        do i = nxl-1,iu_bnd1
         a2(:,i)=pdmb*( v1(:,i-1)-v1(:,i-2)) 
        enddo

        do i = nxl,nxh
         do l=1,nem1
          a3(l,i) = vf(l,i) - v1(l,i-1)
          f1      = max ( a1, a2(l,i) )
          f2      = min ( a1, a2(l,i) )
          f3      = min ( f1, a3(l,i) )
          f4      = max ( f2, f3      )
          xl(l,i) = f4
          a3(l,i) = vf(l,i) - v1(l,i)
          f1      = max ( a1, - a2(l,i+2) )
          f2      = min ( a1, - a2(l,i+2) )
          f3      = min ( f1,   a3(l,i  ) )
          f4      = max ( f2,   f3        )
          xr(l,i) = f4
         enddo
        enddo
!
!-----finally set up left and right states for riemann solver
!
        do i = nxl,nxh
         do l=1,ne
          vl(:nem1,i)=v1(:nem1,i-1)+xl(:nem1,i)
          vr(:nem1,i)=v1(:nem1,i  )+xr(:nem1,i)
         enddo
        enddo


	return
	end subroutine pim_1d
