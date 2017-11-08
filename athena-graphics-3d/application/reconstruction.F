#include "paramesh_preprocessor.fh"
#include "athena_preprocessor.fh"
#include "mydebug.fh"
!=======================================================================
!
!=======================================================================
      subroutine reconstruction(v1,chi,vxfl,vxfr,dt,idir,lflag)      
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
#define cell_v  real,dimension(ne,iu_bnd1 ,ju_bnd1 ,ku_bnd1 )
#define cell_s  real,dimension(   iu_bnd1 ,ju_bnd1 ,ku_bnd1 )
#define face_v  real,dimension(ne,iu_bnd1p,ju_bnd1p,ku_bnd1p)
#define face_f  real,dimension(nf,iu_bnd1p,ju_bnd1p,ku_bnd1p)
      logical                            ,intent(in   )::lflag
      integer                            ,intent(in   )::idir
      real                               ,intent(in   )::dt
      cell_s     ,intent(in   )::chi
      cell_v     ,intent(in   )::v1
      face_v     ,intent(inout)::vxfl,vxfr
      real,dimension(ne,nmax+1)::vl,vr 
      real,dimension(ne,nmax  )::vc 
      real,dimension(   nmax  )::bxs,bys,bzs
      real,dimension(   nmax  )::flat
      integer :: n1,n2,n3,i,j,k,nl1,nh1
!
!-----x-direction
!
      select case (idir)
       case (1)
        n1 =nx ;n2 =ny ;n3=nz
        nl1=nlx;nh1=nhx
       case (2)
        n1 =ny ;n2 =nz ;n3=nx
        nl1=nly;nh1=nhy
       case (3)
        n1 =nz ;n2 =nx ;n3=ny
        nl1=nlz;nh1=nhz
      end select

      do k = 1, n3		! z-loop
       do j = 1, n2		! y-loop
	      
        select case (idir)
         case (1)
	 do i = 1, n1	! x-loop
	  vc(:ne,i)= v1(:ne,i,j,k)
	  bxs(i)   = b_sf(1,i,j,k)
	  bys(i)   = b_sf(2,i,j,k)
	  bzs(i)   = b_sf(3,i,j,k)
	 enddo
	 do i=2,n1-1
	  flat(i) = chi(i,j,k)
	 enddo
#ifdef INTERIOR_BOUNDARY
        if(lflag)call internal_hardwall_bc(i,j,k,1)
#endif /*INTERIOR_BOUNDARY	*/    
        case (2)
	 do i = 1, n1	! y-loop
	  vc(:ne,i)= v1(:ne,k,i,j)
	  bxs(i)   = b_sf(2,k,i,j)
	  bys(i)   = b_sf(3,k,i,j)
	  bzs(i)   = b_sf(1,k,i,j)
	 enddo
	 do i=2,n1-1
	  flat(i) = chi(k,i,j)
	 enddo
#ifdef INTERIOR_BOUNDARY
        if(lflag)call internal_hardwall_bc(k,i,j,2)
#endif /*INTERIOR_BOUNDARY	*/    
        case (3)
	 do i = 1, n1	! z-loop
	  vc(:ne,i)= v1(:ne,j,k,i)
	  bxs(i)   = b_sf(3,j,k,i)
	  bys(i)   = b_sf(1,j,k,i)
	  bzs(i)   = b_sf(2,j,k,i)
	 enddo
	 do i=2,n1-1
	  flat(i) = chi(j,k,i)
	 enddo
#ifdef INTERIOR_BOUNDARY
        if(lflag)call internal_hardwall_bc(j,k,i,3)
#endif /*INTERIOR_BOUNDARY	*/    
       end select


        if(.not.lflag)then 
#ifdef MUSCL
	call muscl_1d(vc,flat,bxs,bys,bzs,vl,vr,1,n1,dt,idir,lflag)
#endif
#ifdef PPM
        call   ppm_1d(vc,flat,bxs,bys,bzs,vl,vr,1,n1,dt,idir,lflag)
#endif /*PPM*/
#ifdef PLM
        call   plm_1d(vc,flat,bxs,bys,bzs,vl,vr,1,n1,dt,idir,lflag)
#endif /*PLM*/
        else
        select case (idir)
         case (1)
          do i = 1,n1 ! x-loop
           vc(:ne,i)= v1(:ne,i,j,k)
          enddo
          call rotate_v(vc,vl,vr,1,n1,idir)
         case (2)
          do i = 1,n1 ! y-loop
           vc(:ne,i)= v1(:ne,k,i,j)
          enddo
          call rotate_v(vc,vl,vr,1,n1,idir)
         case (3)
          do i = 1,n1 ! z-loop
           vc(:ne,i)= v1(:ne,j,k,i)
          enddo
          call rotate_v(vc,vl,vr,1,n1,idir)
        end select

        endif  
!
!-----save left and right states for CTU correction
!         
        select case (idir)
         case (1)
         do i=nl1,nh1
          vxfl(:ne,i,j,k) = vxfl(:ne,i,j,k)+vl(:ne,i)
          vxfr(:ne,i,j,k) = vxfr(:ne,i,j,k)+vr(:ne,i)
         enddo 
         case (2)
         do i=nl1,nh1
          vxfl(:ne,k,i,j) = vxfl(:ne,k,i,j)+vl(:ne,i)
          vxfr(:ne,k,i,j) = vxfr(:ne,k,i,j)+vr(:ne,i)
         enddo 
         case (3)
         do i=nl1,nh1
          vxfl(:ne,j,k,i) = vxfl(:ne,j,k,i)+vl(:ne,i)
          vxfr(:ne,j,k,i) = vxfr(:ne,j,k,i)+vr(:ne,i)
         enddo 
        end select

       enddo  ! end of y-loop over j
      enddo  ! end of z-loop over k
      return
      end subroutine reconstruction
!=======================================================================
!
!=======================================================================
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
      integer,dimension(ne+3),intent(in   )::index_x,index_y 
      face_v                 ,intent(inout)::vxfl,vxfr,vyfl,vyfr
      face_f                 ,intent(in   )::fx,fy
      real                   ,intent(in   )::dtdx,dtdy
      integer :: i,j,k,l,n
!
!-----x-direction
!
      call eos_vf_to_uf(vxfl,vxfr,nlxp1,nhxm1,nlyp1,nhym1,nlzp1,nhzm1,1)
#ifdef SECOND_ORDER_CORRECTION
      do k = nlzp1,nhzm1! z-loop
       do j = nlyp1,nhym1! y-loop
        do i = nlxp1,nhxm1
         do l = 1,ne
          n            =index_x(l)
          vxfl(l,i,j,k)=vxfl(l,i,j,k)-dtdy*(fy(n,i-1,j+1,k)-fy(n,i-1,j,k)) 
          vxfr(l,i,j,k)=vxfr(l,i,j,k)-dtdy*(fy(n,i  ,j+1,k)-fy(n,i  ,j,k))
         enddo
#ifdef LORENTZ_FORCE
!
!------correct for Lorentz force fluxes
!
         do l=2,4		
	  n            =index_x(ne+l-1)
	  vxfl(l,i,j,k)=vxfl(l,i,j,k)-dtdy*(fy(n,i-1,j+1,k)-fy(n,i-1,j,k)) 
	  vxfr(l,i,j,k)=vxfr(l,i,j,k)-dtdy*(fy(n,i  ,j+1,k)-fy(n,i  ,j,k))
	 enddo
#endif
        enddo
       enddo  ! end of y-loop over j
      enddo  ! end of z-loop over k

#endif /*SECOND_ORDER_CORRECTION*/
      call eos_uf_to_vf(vxfl,vxfr,nlxp1,nhxm1,nlyp1,nhym1,nlzp1,nhzm1,1)
!
!-----y-direction
!
      call eos_vf_to_uf(vyfl,vyfr,nlxp1,nhxm1,nlyp1,nhym1,nlzp1,nhzm1,2)
#ifdef SECOND_ORDER_CORRECTION

      do k = nlzp1,nhzm1    ! z-loop
       do i = nlxp1,nhxm1   ! x-loop
        do j = nlyp1,nhym1  ! y-loop
         do l = 1,ne
          n = index_y(l)
        vyfl(l,i,j,k)=vyfl(l,i,j,k)-dtdx*(fx(n,i+1,j-1,k)-fx(n,i,j-1,k)) 
        vyfr(l,i,j,k)=vyfr(l,i,j,k)-dtdx*(fx(n,i+1,j  ,k)-fx(n,i,j  ,k))
         enddo
#ifdef LORENTZ_FORCE
!
!------correct for Lorentz force fluxes
!
         do l=2,4
        n            = index_y(ne+l-1)
        vyfl(l,i,j,k)=vyfl(l,i,j,k)-dtdx*(fx(n,i+1,j-1,k)-fx(n,i,j-1,k)) 
        vyfr(l,i,j,k)=vyfr(l,i,j,k)-dtdx*(fx(n,i+1,j  ,k)-fx(n,i,j  ,k))
         enddo
#endif
        enddo
       enddo  ! end of x-loop over i
      enddo  ! end of z-loop over k
#endif /*SECOND_ORDER_CORRECTION*/
      call eos_uf_to_vf(vyfl,vyfr,nlxp1,nhxm1,nlyp1,nhym1,nlzp1,nhzm1,2)

      return
      end subroutine correct_rl_states_2d
!=======================================================================
!
!=======================================================================
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
      integer,dimension(ne+3),intent(in   )::index_x 
      face_v                 ,intent(in   )::vxfl,vxfr
      face_v                 ,intent(  out)::vxyfl,vxyfr
      face_v                 ,intent(  out)::vxzfl,vxzfr
      face_f                 ,intent(in   )::fy,fz
      real                   ,intent(in   )::dtdy,dtdz
      integer :: i,j,k,l,n
      call eos_vf_to_uf(vxfl,vxfr,nlx,nhx,nly,nhym1,nlz,nhzm1,1)
#ifdef SECOND_ORDER_CORRECTION
c      write(*,*)'x-nlx,nhx,nly,nhym1,nlz,nhzm1 =',nlx,nhx,nly,nhym1,nlz,nhzm1
      do k = nlz,nhzm1  
       do j = nly,nhym1  
        do i = nlx,nhx
         do l = 1,ne
          vxyfl(l,i,j,k)= vxfl(l,i,j,k)-dtdy*(fy(l,i-1,j+1,k)-fy(l,i-1,j,k)) 
          vxyfr(l,i,j,k)= vxfr(l,i,j,k)-dtdy*(fy(l,i  ,j+1,k)-fy(l,i  ,j,k))
         enddo
#ifdef LORENTZ_FORCE
!
!------correct for Lorentz force fluxes
!
         do l = 2,4
          n             =index_x(ne+l-1)
          vxyfl(l,i,j,k)=vxyfl(l,i,j,k)-dtdy*(fy(n,i-1,j+1,k)-
     &                                        fy(n,i-1,j  ,k)) 
          vxyfr(l,i,j,k)=vxyfr(l,i,j,k)-dtdy*(fy(n,i  ,j+1,k)-
     &                                        fy(n,i  ,j  ,k))
         enddo
#endif
        enddo
       enddo  ! end of y-loop over j
      enddo  ! end of z-loop over k
	
#else
      vxyfl(:ne,nlx:nhx,nly:nhym1,nlz:nhzm1)=vxfl(:ne,nlx:nhx,nly:nhym1,nlz:nhzm1)
      vxyfr(:ne,nlx:nhx,nly:nhym1,nlz:nhzm1)=vxfr(:ne,nlx:nhx,nly:nhym1,nlz:nhzm1)
#endif /*SECOND_ORDER_CORRECTION*/
!
!-----convert conserved variables vxyfl & vxyfr back to primitives
!
      call eos_uf_to_vf(vxyfl,vxyfr,nlx,nhx,nly,nhym1,nlz,nhzm1,1)

#ifdef SECOND_ORDER_CORRECTION

      do k = nlz,nhzm1  
       do j = nly,nhym1  
        do i = nlx,nhx
         do l = 1,ne
          vxzfl(l,i,j,k)= vxfl(l,i,j,k)-dtdz*(fz(l,i-1,j,k+1)-fz(l,i-1,j,k)) 
          vxzfr(l,i,j,k)= vxfr(l,i,j,k)-dtdz*(fz(l,i  ,j,k+1)-fz(l,i  ,j,k))
         enddo
#ifdef LORENTZ_FORCE
!
!------correct for Lorentz force fluxes
!
         do l = 2,4
          n             =index_x(ne+l-1)
          vxzfl(l,i,j,k)=vxzfl(l,i,j,k)-dtdz*(fz(n,i-1,j,k+1)-
     &                                        fz(n,i-1,j,k  )) 
          vxzfr(l,i,j,k)=vxzfr(l,i,j,k)-dtdz*(fz(n,i  ,j,k+1)-
     &                                        fz(n,i  ,j,k  ))
         enddo
#endif
        enddo
       enddo  ! end of y-loop over j
      enddo  ! end of z-loop over k


#else
      vxzfl(:ne,nlx:nhx,nly:nhym1,nlz:nhzm1)=vxfl(:ne,nlx:nhx,nly:nhym1,nlz:nhzm1)
      vxzfr(:ne,nlx:nhx,nly:nhym1,nlz:nhzm1)=vxfr(:ne,nlx:nhx,nly:nhym1,nlz:nhzm1)
#endif /*SECOND_ORDER_CORRECTION*/
!
!-----convert conserved variables vxzfl & vxzfr back to primitives
!
      call eos_uf_to_vf(vxzfl,vxzfr,nlx,nhx,nly,nhym1,nlz,nhzm1,1)

      return
      end subroutine correct_rl_states_3d_xa
!=======================================================================
!
!=======================================================================
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
      integer,dimension(ne+3),intent(in   )::index_y 
      face_v                 ,intent(in   )::vyfl ,vyfr
      face_v                 ,intent(  out)::vyxfl,vyxfr
      face_v                 ,intent(  out)::vyzfl,vyzfr
      face_f                 ,intent(in   )::fx,fz
      real                   ,intent(in   )::dtdx,dtdz
      integer :: i,j,k,l,n
!
!-----convert primitives to conserved variables and rotate if needed
!
      call eos_vf_to_uf(vyfl,vyfr,nlx,nhxm1,nly,nhy ,nlz,nhzm1,2)
#ifdef SECOND_ORDER_CORRECTION
c      write(*,*)'y-nlx,nhxm1,nly,nhy,nlz,nhzm1=',nlx,nhxm1,nly,nhy,nlz,nhzm1
      do k = nlz,nhzm1   ! z-loop
       do j = nly,nhy  ! y-loop
        do i = nlx,nhxm1
         do l = 1,ne
          n             =index_y(l)
          vyxfl(l,i,j,k)=vyfl(l,i,j,k)-dtdx*(fx(n,i+1,j-1,k)-
     &                                       fx(n,i  ,j-1,k)) 
          vyxfr(l,i,j,k)=vyfr(l,i,j,k)-dtdx*(fx(n,i+1,j  ,k)-
     &                                       fx(n,i  ,j  ,k))
         enddo
#ifdef LORENTZ_FORCE
!
!------correct for Lorentz force fluxes
!
         do l = 2,4
          n             =index_y(ne+l-1)
          vyxfl(l,i,j,k)=vyxfl(l,i,j,k)-dtdx*(fx(n,i+1,j-1,k)-
     &                                        fx(n,i  ,j-1,k)) 
          vyxfr(l,i,j,k)=vyxfr(l,i,j,k)-dtdx*(fx(n,i+1,j  ,k)-
     &                                        fx(n,i  ,j  ,k))
         enddo
#endif
        enddo
       enddo   
      enddo   


#else
      vyxfl(:ne,nlx:nhxm1,nly:nhy,nlz:nhzm1)=vyfl(:ne,nlx:nhxm1,nly:nhy,nlz:nhzm1)
      vyxfr(:ne,nlx:nhxm1,nly:nhy,nlz:nhzm1)=vyfr(:ne,nlx:nhxm1,nly:nhy,nlz:nhzm1)
#endif /*SECOND_ORDER_CORRECTION*/
!
!-----convert conserved variables vyxfl & vyxfr back to primitives
!
      call eos_uf_to_vf(vyxfl,vyxfr,nlx,nhxm1,nly,nhy ,nlz,nhzm1,2)

#ifdef SECOND_ORDER_CORRECTION
      do k = nlz,nhzm1   ! z-loop
       do j = nly,nhy  ! y-loop
        do i = nlx,nhxm1
         do l = 1,ne
          n             =index_y(l)
          vyzfl(l,i,j,k)=vyfl(l,i,j,k)-dtdz*(fz(n,i,j-1,k+1)-
     &                                       fz(n,i,j-1,k  )) 
          vyzfr(l,i,j,k)=vyfr(l,i,j,k)-dtdz*(fz(n,i,j  ,k+1)-
     &                                       fz(n,i,j  ,k  ))
         enddo
#ifdef LORENTZ_FORCE
!
!------correct for Lorentz force fluxes
!
         do l = 2,4
          n             =index_y(ne+l-1)
          vyzfl(l,i,j,k)=vyzfl(l,i,j,k)-dtdz*(fz(n,i,j-1,k+1)-
     &                                        fz(n,i,j-1,k  )) 
          vyzfr(l,i,j,k)=vyzfr(l,i,j,k)-dtdz*(fz(n,i,j  ,k+1)-
     &                                        fz(n,i,j  ,k  ))
         enddo
#endif
        enddo
       enddo  ! end of x-loop over i
      enddo  ! end of z-loop over k

#else
      vyzfl(:ne,nlx:nhxm1,nly:nhy,nlz:nhzm1)=vyfl(:ne,nlx:nhxm1,nly:nhy,nlz:nhzm1)
      vyzfr(:ne,nlx:nhxm1,nly:nhy,nlz:nhzm1)=vyfr(:ne,nlx:nhxm1,nly:nhy,nlz:nhzm1)
#endif /*SECOND_ORDER_CORRECTION*/
!
!-----convert conserved variables vyzfl & vyzfr back to primitives
!
      call eos_uf_to_vf(vyzfl,vyzfr,nlx,nhxm1,nly,nhy,nlz,nhzm1,2)

      return
      end subroutine correct_rl_states_3d_ya
!=======================================================================
!
!=======================================================================
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
      integer,dimension(ne+3),intent(in   )::index_z 
      face_v                 ,intent(in   )::vzfl ,vzfr
      face_v                 ,intent(  out)::vzxfl,vzxfr
      face_v                 ,intent(  out)::vzyfl,vzyfr
      face_f                 ,intent(in   )::fx,fy
      real                   ,intent(in   )::dtdx,dtdy
      integer :: i,j,k,l,n
!
!-----convert primitives to conserved variables and rotate if needed
!
!
!-----z-direction
!
!
!-----convert primitives to conserved variables and rotate if needed
!
      do k = nlzp1,nhzm1   
       do j = nlyp1,nhym1  
        do i = nlxp1,nhxm1
      if(vzfl(1,i,j,k).le.0.0.or.vzfl(1,i,j,k).le.0.0)then
      write(*,*)'vzfl(1,',i,',',j,',',k,'),vzfl(1,',i,',',j,',',k,'),aa=',
     &vzfl(1,i,j,k),vzfl(1,i,j,k)
          stop
	  endif
	  enddo
       enddo  ! end of x-loop over i
      enddo  ! end of y-loop over j
      call eos_vf_to_uf(vzfl,vzfr,nlx,nhxm1,nly,nhym1,nlz,nhz,3)
#ifdef SECOND_ORDER_CORRECTION
c      write(*,*)'x-nlx,nhxm1,nly,nhym1,nlz,nhz =',nlx,nhxm1,nly,nhym1,nlz,nhz

      do k = nlz,nhz   
       do j = nly,nhym1  
        do i = nlx,nhxm1
         do l = 1,ne
          n             =index_z(l)
      if(vzfl(1,i,j,k).le.0.0.or.vzfl(1,i,j,k).le.0.0)then
      write(*,*)'vzfl(1,',i,',',j,',',k,'),vzfl(1,',i,',',j,',',k,'),a=',
     &vzfl(1,i,j,k),vzfl(1,i,j,k)	     
          stop
	  endif
          vzxfl(l,i,j,k)=vzfl(l,i,j,k)-dtdx*(fx(n,i+1,j,k-1)-
     &                                       fx(n,i  ,j,k-1)) 
          vzxfr(l,i,j,k)=vzfr(l,i,j,k)-dtdx*(fx(n,i+1,j,k  )-
     &                                       fx(n,i  ,j,k  ))
      if(vzfl(1,i,j,k).le.0.0.or.vzfl(1,i,j,k).le.0.0)then
      write(*,*)'vzxfl(1,',i,',',j,',',k,'),vzxfl(1,',i,',',j,',',k,'),b=',
     &vzxfl(1,i,j,k),vzxfl(1,i,j,k)	     
	  endif
         enddo
#ifdef LORENTZ_FORCE
!
!------correct for Lorentz force fluxes
!
         do l = 2,4
          n             =index_z(ne+l-1)
          vzxfl(l,i,j,k)=vzxfl(l,i,j,k)-dtdx*(fx(n,i+1,j,k-1)-
     &                                        fx(n,i  ,j,k-1)) 
          vzxfr(l,i,j,k)=vzxfr(l,i,j,k)-dtdx*(fx(n,i+1,j,k  )-
     &                                        fx(n,i  ,j,k  ))
         enddo
#endif
        enddo
       enddo   
      enddo   
#else
      vzxfl(:ne,nlx:nhxm1,nly:nhym1,nlz:nhz)=vzfl(:ne,nlx:nhxm1,nly:nhym1,nlz:nhz)
      vzxfr(:ne,nlx:nhxm1,nly:nhym1,nlz:nhz)=vzfr(:ne,nlx:nhxm1,nly:nhym1,nlz:nhz)
#endif /*SECOND_ORDER_CORRECTION*/
!
!-----convert conserved variables vzxfl & vzxfr back to primitives
!
      call eos_uf_to_vf(vzxfl,vzxfr,nlx,nhxm1,nly,nhym1,nlz,nhz,3)

!-----------------------------------------------------------------------------
#ifdef SECOND_ORDER_CORRECTION
      do k = nlz,nhz   
       do j = nly,nhym1  
        do i = nlx,nhxm1
         do l = 1,ne
          n             =index_z(l)
      if(vzfl(1,i,j,k).le.0.0.or.vzfl(1,i,j,k).le.0.0)then
      write(*,*)'vzfl(1,',i,',',j,',',k,'),vzfl(1,',i,',',j,',',k,'),c=',
     &vzyfl(1,i,j,k),vzyfl(1,i,j,k)	     
          stop
	  endif
          vzyfl(l,i,j,k)=vzfl(l,i,j,k)-dtdy*(fy(n,i,j+1,k-1)-
     &                                       fy(n,i,j  ,k-1)) 
          vzyfr(l,i,j,k)=vzfr(l,i,j,k)-dtdy*(fy(n,i,j+1,k  )-
     &                                       fy(n,i,j  ,k  ))
       if(vzyfl(1,i,j,k).le.0.0.or.vzyfl(1,i,j,k).le.0.0)then
      write(*,*)'vzyfl(1,',i,',',j,',',k,'),vzyfl(1,',i,',',j,',',k,'),d=',
     &vzyfl(1,i,j,k),vzyfl(1,i,j,k)	     
          stop
	  endif
        enddo
#ifdef LORENTZ_FORCE
!
!------correct for Lorentz force fluxes
!
         do l = 2,4
          n             =index_z(ne+l-1)
          vzyfl(l,i,j,k)=vzyfl(l,i,j,k)-dtdy*(fy(n,i,j+1,k-1)-
     &                                        fy(n,i,j  ,k-1)) 
          vzyfr(l,i,j,k)=vzyfr(l,i,j,k)-dtdy*(fy(n,i,j+1,k  )-
     &                                        fy(n,i,j  ,k  ))
         enddo
#endif
        enddo
       enddo  ! end of x-loop over i
      enddo  ! end of y-loop over j
#else
      vzyfl(:ne,nlx:nhxm1,nly:nhym1,nlz:nhz)=vzfl(:ne,nlx:nhxm1,nly:nhym1,nlz:nhz)
      vzyfr(:ne,nlx:nhxm1,nly:nhym1,nlz:nhz)=vzfr(:ne,nlx:nhxm1,nly:nhym1,nlz:nhz)
#endif /*SECOND_ORDER_CORRECTION*/
!
!-----convert conserved variables vzyfl & vzyfr back to primitives
!
       call eos_uf_to_vf(vzyfl,vzyfr,nlx,nhxm1,nly,nhym1,nlz,nhz,3)
      return
      end subroutine correct_rl_states_3d_za
!=======================================================================
!
!=======================================================================
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
      integer,dimension(ne+3),intent(in   )::index_x 
      face_v                 ,intent(inout)::vxfl,vxfr
      face_f                 ,intent(in   )::fzy,fyz
      real                             ,intent(in   )::dtdy,dtdz
      integer :: i,j,k,l,n
#ifdef SECOND_ORDER_CORRECTION
      do k = nlzp1,nhzm1   
       do j = nlyp1,nhym1  
        do i = nlxp1,nhxm1
         do l=1,ne
          vxfl(l,i,j,k)= vxfl(l,i,j,k)-dtdy*(fyz(l,i-1,j+1,k)-fyz(l,i-1,j,k)) 
          vxfr(l,i,j,k)= vxfr(l,i,j,k)-dtdy*(fyz(l,i  ,j+1,k)-fyz(l,i  ,j,k))
         enddo
#ifdef LORENTZ_FORCE
!
!------correct for Lorentz force fluxes
!
         do l=2,4
        n             =index_x(ne+l-1)
        vxfl(l,i,j,k)=vxfl(l,i,j,k)-dtdy*(fyz(n,i-1,j+1,k)-
     &                                    fyz(n,i-1,j  ,k)) 
        vxfr(l,i,j,k)=vxfr(l,i,j,k)-dtdy*(fyz(n,i  ,j+1,k)-
     &                                    fyz(n,i  ,j  ,k))
         enddo
#endif
!------------------------------------------------------------------------------
         do l=1,ne
          vxfl(l,i,j,k)= vxfl(l,i,j,k)-dtdz*(fzy(l,i-1,j,k+1)-fzy(l,i-1,j,k)) 
          vxfr(l,i,j,k)= vxfr(l,i,j,k)-dtdz*(fzy(l,i  ,j,k+1)-fzy(l,i  ,j,k))
         enddo
#ifdef LORENTZ_FORCE
!
!------correct for Lorentz force fluxes
!
         do l=2,4
        n            =index_x(ne+l-1)
        vxfl(l,i,j,k)=vxfl(l,i,j,k)-dtdz*(fzy(n,i-1,j,k+1)-
     &                                    fzy(n,i-1,j,k  )) 
        vxfr(l,i,j,k)=vxfr(l,i,j,k)-dtdz*(fzy(n,i  ,j,k+1)-
     &                                    fzy(n,i  ,j,k  ))
         enddo
#endif
        enddo
       enddo  ! end of y-loop over j
      enddo  ! end of z-loop over k	


#endif /*SECOND_ORDER_CORRECTION*/
      call eos_uf_to_vf(vxfl,vxfr,nlxp1,nhxm1,nlyp1,nhym1,nlzp1,nhzm1,1)

      return
      end subroutine correct_rl_states_3d_xb
!=======================================================================
!
!=======================================================================
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
      integer,dimension(ne+3),intent(in   )::index_y 
      face_v                 ,intent(inout)::vyfl,vyfr
      face_f                 ,intent(in   )::fxz,fzx
      real                   ,intent(in   )::dtdx,dtdz
      integer :: i,j,k,l,n
#ifdef SECOND_ORDER_CORRECTION
      do k = nlzp1,nhzm1   
       do j = nlyp1,nhym1  
        do i = nlxp1,nhxm1
         do l=1,ne
        n            =index_y(l)
        vyfl(l,i,j,k)=vyfl(l,i,j,k)-dtdx*(fxz(n,i+1,j-1,k)-
     &                                    fxz(n,i  ,j-1,k)) 
        vyfr(l,i,j,k)=vyfr(l,i,j,k)-dtdx*(fxz(n,i+1,j  ,k)-
     &                                    fxz(n,i  ,j  ,k))
         enddo
#ifdef LORENTZ_FORCE
!
!------correct for Lorentz force fluxes
!
         do l=2,4
          n            =index_y(ne+l-1)
          vyfl(l,i,j,k)=vyfl(l,i,j,k)-dtdx*(fxz(n,i+1,j-1,k)-
     &                                      fxz(n,i  ,j-1,k)) 
          vyfr(l,i,j,k)=vyfr(l,i,j,k)-dtdx*(fxz(n,i+1,j  ,k)-
     &                                      fxz(n,i  ,j  ,k))
         enddo
#endif
!------------------------------------------------------------------------------
         do l=1,ne
          n            =index_y(l)
          vyfl(l,i,j,k)=vyfl(l,i,j,k)-dtdz*(fzx(n,i,j-1,k+1)-
     &                                      fzx(n,i,j-1,k  )) 
          vyfr(l,i,j,k)=vyfr(l,i,j,k)-dtdz*(fzx(n,i,j  ,k+1)-
     &                                      fzx(n,i,j  ,k  ))
         enddo
#ifdef LORENTZ_FORCE
!
!------correct for Lorentz force fluxes
!
         do l=2,4
          n            =index_y(ne+l-1)
          vyfl(l,i,j,k)=vyfl(l,i,j,k)-dtdz*(fzx(n,i,j-1,k+1)-
     &                                      fzx(n,i,j-1,k  )) 
          vyfr(l,i,j,k)=vyfr(l,i,j,k)-dtdz*(fzx(n,i,j  ,k+1)-
     &                                      fzx(n,i,j  ,k  ))
         enddo
#endif

        enddo
       enddo  ! end of x-loop over i
      enddo  ! end of z-loop over k

#endif /*SECOND_ORDER_CORRECTION*/
      call eos_uf_to_vf(vyfl,vyfr,nlxp1,nhxm1,nlyp1,nhym1,nlzp1,nhzm1,2)


      return
      end subroutine correct_rl_states_3d_yb
!=======================================================================
!
!=======================================================================
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
      integer,dimension(ne+3),intent(in   )::index_z 
      face_v                 ,intent(inout)::vzfl,vzfr
      face_f                 ,intent(in   )::fxy,fyx
      real                   ,intent(in   )::dtdx,dtdy
      integer :: i,j,k,l,n
!#undef SECOND_ORDER_CORRECTION
#ifdef SECOND_ORDER_CORRECTION
      do k = nlzp1,nhzm1   
       do j = nlyp1,nhym1  
        do i = nlxp1,nhxm1
         do l = 1,ne
          n            =index_z(l)
          if(vzfl(1,i,j,k).le.0.0.or.vzfl(1,i,j,k).le.0.0)then
	  write(*,*)'a:vzfl(1,',i,',',j,',',k,'),vzfl(1,',i,',',j,',',k,'),a=',
     &vzfl(1,i,j,k),vzfl(1,i,j,k)	     
          stop
	  endif
          vzfl(l,i,j,k)=vzfl(l,i,j,k)-dtdx*(fxy(n,i+1,j,k-1)-
     &                                      fxy(n,i  ,j,k-1)) 
          vzfr(l,i,j,k)=vzfr(l,i,j,k)-dtdx*(fxy(n,i+1,j,k  )-
     &                                      fxy(n,i  ,j,k  ))
          if(vzfl(1,i,j,k).le.0.0.or.vzfl(1,i,j,k).le.0.0)then
       write(*,*)'b:vzfl(1,',i,',',j,',',k,'),vzfl(1,',i,',',j,',',k,'),b=',
     & vzfl(1,i,j,k),vzfl(1,i,j,k)
       write(*,*)'dtdx   ,dtdy=',dtdx   ,dtdy	     
       write(*,*)'fxy(n,i+1,j,k-1),fxy(n,i  ,j,k-1)=',fxy(n,i+1,j,k-1),fxy(n,i  ,j,k-1)	     
       write(*,*)'fxy(n,i+1,j,k  ),fxy(n,i  ,j,k  )=',fxy(n,i+1,j,k  ),fxy(n,i  ,j,k  )	     
       write(*,*)'dtdx*(fxy(n,i+1,j,k-1)-fxy(n,i  ,j,k-1))=',dtdx*(fxy(n,i+1,j,k-1)-fxy(n,i  ,j,k-1))	     
       write(*,*)'dtdx*(fxy(n,i+1,j,k  )-fxy(n,i  ,j,k  ))=',dtdx*(fxy(n,i+1,j,k  )-fxy(n,i  ,j,k  ))	     
          stop
	  endif
         enddo
#ifdef LORENTZ_FORCE
!
!------correct for Lorentz force fluxes
!
         do l=2,4
          n            =index_z(ne+l-1)
          vzfl(l,i,j,k)=vzfl(l,i,j,k)-dtdx*(fxy(n,i+1,j,k-1)-
     &                                      fxy(n,i  ,j,k-1)) 
          vzfr(l,i,j,k)=vzfr(l,i,j,k)-dtdx*(fxy(n,i+1,j,k  )-
     &                                      fxy(n,i  ,j,k  ))
         enddo
#endif
!-----------------------------------------------------------------------------
         do l = 1,ne
          n            =index_z(l)
          if(vzfl(1,i,j,k).le.0.0.or.vzfl(1,i,j,k).le.0.0)then
      write(*,*)'vzfl(1,',i,',',j,',',k,'),vzfl(1,',i,',',j,',',k,'),c=',
     &vzfl(1,i,j,k),vzfl(1,i,j,k)	     
	  endif
          vzfl(l,i,j,k)=vzfl(l,i,j,k)-dtdy*(fyx(n,i,j+1,k-1)-
     &                                      fyx(n,i,j  ,k-1)) 
          vzfr(l,i,j,k)=vzfr(l,i,j,k)-dtdy*(fyx(n,i,j+1,k  )-
     &                                      fyx(n,i,j  ,k  ))
          if(vzfl(1,i,j,k).le.0.0.or.vzfl(1,i,j,k).le.0.0)then
      write(*,*)'vzfl(1,',i,',',j,',',k,'),vzfl(1,',i,',',j,',',k,'),d=',
     &vzfl(1,i,j,k),vzfl(1,i,j,k)	     
          stop
	  endif
         enddo
#ifdef LORENTZ_FORCE
!
!------correct for Lorentz force fluxes
!
         do l = 2,4
          n            =index_z(ne+l-1)
          vzfl(l,i,j,k)=vzfl(l,i,j,k)-dtdy*(fyx(n,i,j+1,k-1)-
     &                                      fyx(n,i,j  ,k-1)) 
          vzfr(l,i,j,k)=vzfr(l,i,j,k)-dtdy*(fyx(n,i,j+1,k  )-
     &                                      fyx(n,i,j  ,k  ))
         enddo
#endif
        enddo
       enddo  ! end of x-loop over i
      enddo  ! end of y-loop over j
#endif /*SECOND_ORDER_CORRECTION*/
      call eos_uf_to_vf(vzfl,vzfr,nlxp1,nhxm1,nlyp1,nhym1,nlzp1,nhzm1,3)
      do k = nlzp1,nhzm1   
       do j = nlyp1,nhym1  
        do i = nlxp1,nhxm1
      if(vzfl(1,i,j,k).le.0.0.or.vzfl(1,i,j,k).le.0.0)then
      write(*,*)'vzfl(1,',i,',',j,',',k,'),vzfl(1,',i,',',j,',',k,'),e=',
     &vzfl(1,i,j,k),vzfl(1,i,j,k)	     
          stop
	  endif
	  enddo
       enddo  ! end of x-loop over i
      enddo  ! end of y-loop over j
      return
      end subroutine correct_rl_states_3d_zb

