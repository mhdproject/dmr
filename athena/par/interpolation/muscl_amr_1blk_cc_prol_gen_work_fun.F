!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"


      subroutine amr_1blk_cc_prol_gen_work_fun(ia,ib,ja,jb,ka,kb,
     .       idest,ioff,joff,koff,mype)




!------------------------------------------------------------------------
!
! This routine takes data from the array recvw1, originally extracted 
! from the solution array work1, and performs a prolongation
! operation on it. The data in recvw1 is from a parent block and the
! result of the prolongation operation is written directly into layer
! idest of the working block work array work1.
! The position of the child block, isg, within the 
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
! Written :     Peter MacNeice          September 1999
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace
      use prolong_arrays

!------------------------------------
! local arrays


      integer, intent(in) :: ia,ib,ja,jb,ka,kb,idest
      integer, intent(in) :: ioff,joff,koff,mype

      real :: recvwloc(ilw1:iuw1,jlw1:juw1,klw1:kuw1)

!------------------------------------
! local scalars

! declare large even integer
      integer largei
      parameter(largei=100)

      real	dx,dy,dz,gradl,gradr,gradc,gradm,ss
      integer	ii,jj,kk

!------------------------------------

        if(nguard_work.le.1) then
          write(*,*) ' Error - the muscl interpolation version of',
     .               ' the routine amr_1blk_cc_prol_gen_work_fun ',
     .               ' requires nguard_work > 1.'
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


	ilow = (ia-nguard_work-1+largei)/2+nguard_work+ioff
     .                              -largei/2
	ihi  = (ib-nguard_work-1+largei)/2+nguard_work+ioff
     .                               -largei/2+2
	jlow = ((ja-nguard_work-1+largei)/2+nguard_work+joff
     .                               -largei/2-1)*k2d+1
	jhi  = ((jb-nguard_work-1+largei)/2+nguard_work+joff
     .                               -largei/2+1)*k2d+1
	klow = ((ka-nguard_work-1+largei)/2+nguard_work+koff
     .                               -largei/2-1)*k3d+1
	khi  = ((kb-nguard_work-1+largei)/2+nguard_work+koff
     .                               -largei/2+1)*k3d+1

        recvwloc(ilw1:iuw1,jlw1:juw1,klw1:kuw1) =
     .        recvw1(ilw1:iuw1,jlw1:juw1,klw1:kuw1,1)
!
! Perform sweep in x direction
        do k=klow,khi
        do j=jlow,jhi
        do i=ia,ib
           ii = (i-nguard_work-1+largei)/2+nguard_work+1
     .                                 +ioff-largei/2
              gradl = (recvwloc(ii,j,k)-recvwloc(ii-1,j,k))*dxpr
              gradr = (recvwloc(ii+1,j,k)-recvwloc(ii,j,k))*dxpr
              gradc = (recvwloc(ii+1,j,k)-recvwloc(ii-1,j,k))
     .                                                  *dxpr*.5
              ss = sign(1.,gradc)
              gradm = ss*max(0.,min(abs(gradc),2.*gradr*ss,
     .                                         2.*gradl*ss))
              is = mod(i-nguard_work-1+largei,2)
              sdx = real(2*is-1)*.5
              tempw1(i,j,k) = recvwloc(ii,j,k)+gradm*sdx
        enddo
        enddo
        enddo

!
! Perform sweep in y direction
	if( ndim.ge.2) then

        do k=klow,khi
        do j=jlow,jhi
        do i=ia,ib
	recvwloc(i,j,k) = tempw1(i,j,k)
        enddo
        enddo
        enddo

        do k=klow,khi
        do i=ia,ib
        do j=ja,jb
           jj = (j-nguard_work-1+largei)/2+nguard_work+1
     .                          +joff-largei/2
              gradl = (recvwloc(i,jj,k)-recvwloc(i,jj-k2d,k))
     .                                                  *dypr
              gradr = (recvwloc(i,jj+k2d,k)-recvwloc(i,jj,k))
     .                                                  *dypr
              gradc = (recvwloc(i,jj+k2d,k)-recvwloc(i,jj-k2d,k))
     .                                                  *dypr*.5
              ss = sign(1.,gradc)
              gradm = ss*max(0.,min(abs(gradc),2.*gradr*ss,
     .                                         2.*gradl*ss))
              is = mod(j-nguard_work-1+largei,2)
              sdy = real(2*is-1)*.5
              tempw1(i,j,k) = recvwloc(i,jj,k)+gradm*sdy
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
        recvwloc(i,j,k) = tempw1(i,j,k)
        enddo
        enddo
        enddo

        do i=ia,ib
        do j=ja,jb
        do k=ka,kb
           kk = (k-nguard_work-1+largei)/2+nguard_work
     .                           +1+koff-largei/2
              gradl = (recvwloc(i,j,kk)-recvwloc(i,j,kk-k3d))*dzpr
              gradr = (recvwloc(i,j,kk+k3d)-recvwloc(i,j,kk))*dzpr
              gradc = (recvwloc(i,j,kk+k3d)-recvwloc(i,j,kk-k3d))
     .                                                   *dzpr*.5
              ss = sign(1.,gradc)
              gradm = ss*max(0.,min(abs(gradc),2.*gradr*ss,
     .                                         2.*gradl*ss))
              is = mod(k-nguard_work-1+largei,2)
              sdz = real(2*is-1)*.5
              tempw1(i,j,k) = recvwloc(i,j,kk)+gradm*sdz
        enddo
        enddo
        enddo

	endif

! save new interpolated solution

        do k=ka,kb
        do j=ja,jb
        do i=ia,ib
           work1(i,j,k,idest)=tempw1(i,j,k)
        enddo
        enddo
        enddo


      return
      end subroutine amr_1blk_cc_prol_gen_work_fun
