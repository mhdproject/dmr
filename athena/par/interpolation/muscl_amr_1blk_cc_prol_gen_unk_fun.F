!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"

      subroutine amr_1blk_cc_prol_gen_unk_fun(recv,
     .         ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff,mype)





!------------------------------------------------------------------------
!
! This routine takes data from the array recv, originally extracted 
! from the solution array unk, and performs a prolongation
! operation on it. The data in recv is from a parent block and the
! result of the prolongation operation is written into layer idest
! of the working block array unk1.
! The position of the child, block isg,  within the 
! parent block is specified by the ioff, joff and koff arguments.
! The argument jface allows the call to limit its effect to a specific
! face of the block if required. If jface is set to a value between 1 and
! 6, then guard cells for that face are set. If jface is not between 1 to
! 6 then the prolongation operation is applied to the whole block.
!
! This particular prolongation implements a Muscl style monotonization
! which guarantees a conservative second order interpolation. It can
! only be used for blocks with an even number of grid cells.
!
! It will only work if nguard > 1.
!
!
! The 1D definition is as follows:
!
! if child grid cell ic1 and ic2 are the two child cells corresponding
! to cell ip of the parent, then
!
!                 U(ic1) = U(ip) - gradm * .5 * dxc
!                 U(ic2) = U(ip) + gradm * .5 * dxc
!
! where dxc is the cell size of the children, and
!
!        gradm = s * max( 0. , min( abs(gradc) , 2*s*gradl, 2*s*gradr ) )
!
! with
!        gradc = (U(ip+1) - U(ip-1) )/ (2.*dxp)
!        gradl = (U(ip)   - U(ip-1) )/ dxp
!        gradr = (U(ip+1) - U(ip)   )/ dxp
!
! where dxp (=2*dxc) is the cell size of the children.
! The multidimensional implementation applies this 1D operation first
! in the x direction, then y and finally z if required.
!
!
!
! Written :     Peter MacNeice          September 1999
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use prolong_arrays

!------------------------------------
! local arrays
      real :: recv1(nvar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,
     .                   kl_bnd1:ku_bnd1)
      real :: temp(nvar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,
     .                  kl_bnd1:ku_bnd1)


      integer, intent(in) :: ia,ib,ja,jb,ka,kb,idest
      integer, intent(in) :: ioff,joff,koff,mype
      real,    intent(in) :: recv(:,:,:,:)

!------------------------------------
! local scalars

! declare large even integer
      integer,parameter :: largei=100

      real	:: dx,dy,dz,gradl,gradr,gradc,gradm,ss
      integer	:: ii,jj,kk,ivar

      integer   :: ierrorcode,ierr

!------------------------------------
        if(nguard.le.1) then
          write(*,*) ' Error - the muscl interpolation version of',
     .               ' the routine amr_1blk_cc_prol_gen_unk_fun ',
     .               ' requires nguard > 1.'
#ifdef SHMEM_USED      
          stop
#endif
#ifdef MPI_USED      
          call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
#endif
        endif

! reciprocal of parent cell size
        dxpr = .5
        dypr = .5
        dzpr = .5


	ilow = (ia-nguard-1+largei)/2+nguard+ioff-largei/2
	ihi  = (ib-nguard-1+largei)/2+nguard+ioff-largei/2+2
	jlow = ((ja-nguard-1+largei)/2+nguard+joff-largei/2-1)*k2d+1
	jhi  = ((jb-nguard-1+largei)/2+nguard+joff-largei/2+1)*k2d+1
	klow = ((ka-nguard-1+largei)/2+nguard+koff-largei/2-1)*k3d+1
	khi  = ((kb-nguard-1+largei)/2+nguard+koff-largei/2+1)*k3d+1

!
! Perform sweep in x direction
        do k=klow,khi
        do j=jlow,jhi
        do i=ia,ib
           ii = (i-nguard-1+largei)/2+nguard+1+ioff-largei/2
           do ivar=1,nvar
              gradl = (recv(ivar,ii,j,k)-recv(ivar,ii-1,j,k))
     .                                         *dxpr
              gradr = (recv(ivar,ii+1,j,k)-recv(ivar,ii,j,k))
     .                                         *dxpr
              gradc = (recv(ivar,ii+1,j,k)-recv(ivar,ii-1,j,k))
     .                                         *dxpr*.5
              ss = sign(1.,gradc)
              gradm = ss*max(0.,min(abs(gradc),2.*gradr*ss,
     .                                         2.*gradl*ss))
              is = mod(i-nguard-1+largei,2)
              sdx = real(2*is-1)*.5
              temp(ivar,i,j,k) = recv(ivar,ii,j,k)+gradm*sdx
           enddo
        enddo
        enddo
        enddo

!
! Perform sweep in y direction
	if( ndim.ge.2) then

        do k=klow,khi
        do j=jlow,jhi
        do i=ia,ib
        do ivar=1,nvar
	recv1(ivar,i,j,k) = temp(ivar,i,j,k)
        enddo
        enddo
        enddo
        enddo

        do k=klow,khi
        do i=ia,ib
        do j=ja,jb
           jj = (j-nguard-1+largei)/2+nguard+1+joff-largei/2
           do ivar=1,nvar
              gradl = (recv1(ivar,i,jj,k)-recv1(ivar,i,jj-k2d,k))
     .                                         *dypr
              gradr = (recv1(ivar,i,jj+k2d,k)-recv1(ivar,i,jj,k))
     .                                         *dypr
              gradc = (recv1(ivar,i,jj+k2d,k)-recv1(ivar,i,jj-k2d,k))
     .                                         *dypr*.5
              ss = sign(1.,gradc)
              gradm = ss*max(0.,min(abs(gradc),2.*gradr*ss,
     .                                         2.*gradl*ss))
              is = mod(j-nguard-1+largei,2)
              sdy = real(2*is-1)*.5
              temp(ivar,i,j,k) = recv1(ivar,i,jj,k)+gradm*sdy
           enddo
        enddo
        enddo
        enddo

	endif

!
! Perform sweep in z direction
	if(ndim.eq.3) then

        do k=klow,khi
        do j=ja,jb
        do i=ia,ib
        do ivar=1,nvar
        recv1(ivar,i,j,k) = temp(ivar,i,j,k)
        enddo
        enddo
        enddo
        enddo

        do i=ia,ib
        do j=ja,jb
        do k=ka,kb
           kk = (k-nguard-1+largei)/2+nguard+1+koff-largei/2
           do ivar=1,nvar
              gradl = (recv1(ivar,i,j,kk)-recv1(ivar,i,j,kk-k3d))
     .                                         *dzpr
              gradr = (recv1(ivar,i,j,kk+k3d)-recv1(ivar,i,j,kk))
     .                                         *dzpr
              gradc = (recv1(ivar,i,j,kk+k3d)-recv1(ivar,i,j,kk-k3d))
     .                                         *dzpr*.5
              ss = sign(1.,gradc)
              gradm = ss*max(0.,min(abs(gradc),2.*gradr*ss,
     .                                         2.*gradl*ss))
              is = mod(k-nguard-1+largei,2)
              sdz = real(2*is-1)*.5
              temp(ivar,i,j,k) = recv1(ivar,i,j,kk)+gradm*sdz
           enddo
        enddo
        enddo
        enddo

	endif

! save new interpolated solution

        do k=ka,kb
        do j=ja,jb
        do i=ia,ib
        do ivar=1,nvar
           unk1(ivar,i,j,k,idest)=temp(ivar,i,j,k)
        enddo
        enddo
        enddo
        enddo


      return
      end subroutine amr_1blk_cc_prol_gen_unk_fun
