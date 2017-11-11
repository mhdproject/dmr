#include "paramesh_preprocessor.fh"
      subroutine error_measure(mype,nprocs,iopt,lcc,lfc,error,epsil,lb)
!---------------------------------------------------------------------
      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace
      use grid

      use paramesh_interfaces
#ifdef MPI_USED
      use paramesh_mpi_interfaces
#endif /* MPI_USED */
      implicit none
#ifdef MPI_USED
      include "mpif.h"
      integer :: tag_offset,max_blks_sent
#endif /* MPI_USED */
     
      integer :: mype,nprocs,iopt,idest,icoord,nlayers

#include "amr_shmem.fh"
      integer, parameter :: ndim2=ndim*ndim
      integer :: AMR_SOURCE_FILES,errcode
      real,save,dimension(ndim, il_bnd1:iu_bnd1,
     &                      jl_bnd1:ju_bnd1,
     &                      kl_bnd1:ku_bnd1) ::delu,delua
      real,save,dimension(ndim2,il_bnd1:iu_bnd1,
     &                      jl_bnd1:ju_bnd1,
     &                      kl_bnd1:ku_bnd1) ::delu2,delu3,delu4


      logical :: lguard,lprolong,lflux,ledge,lrestrict

      integer :: lb,i,j,k,ii,jj,kk,ndel
      logical :: lcc,lfc,l_srl_only,ldiag,lec,lnc
c      real :: delx,dely,delz
      real,save :: errort,num,denom
      real  :: error,epsil
c      delx = bsize(1,lb)/float(nxb)
c      dely = bsize(2,lb)/float(nyb)
c#if N_DIM == 3
c      delz = bsize(3,lb)/float(nzb)
c#endif
c      delx = 1./(2.*delx)
c      dely = 1./(2.*dely)
c#if N_DIM == 3 
c      delz = 1./(2.*delz)
c#endif
       call grid_variables(lb,0.0)
#ifdef NO_PERMANENT_GUARDCELLS
! Copy data from current block into working block and fill its guardcells
        idest      = 1
c        iopt       = 2
        nlayers    = nguard
c        lcc        = .true.
c        lfc        = .false.
        l_srl_only = .false.
        icoord     = 0
        ldiag      = .true.
        lec        = .false.
        lnc        = .false.
!
!-----call amr_1blk_guardcell: selects a block including guard cells
!
c      write(*,*)'entrying amr_1blk_guardcell',mype,lb

!#ifdef MPI_USED
!      call mpi_amr_1blk_guardcell(mype,iopt,nlayers,lb,mype,lcc,lfc,
!     .                              l_srl_only,icoord,ldiag,
!     .                              nprocs)
!#else  /* MPI_USED */

      call amr_1blk_guardcell(mype,iopt,nlayers,lb,mype,lcc,lfc,
     .                        lec,lnc,
     .                        l_srl_only,icoord,ldiag)
!#endif /* MPI_USED */
c      write(*,*)'leaving amr_1blk_guardcell',mype,lb
!
!-----unk1 now has the unknown and the necessary guard cell information
!     for block lb
!
#endif
c       if(mype.eq.0.and.lb.eq.30)then
c                write(*,*)'lb,delx,dely,delz=',lb,delx,dely,delz
c       endif
c       if(mype.eq.1.and.lb.eq.42)then
c                write(*,*)'lb,delx,dely,delz=',lb,delx,dely,delz
c       endif

! Compute first derivatives

      ndel=nguard-nguard_work
      do k = 1+k3d*(nguard-1),nzb+k3d*(nguard+1)
         kk=k-k3d*ndel
         do j = nguard,nyb+nguard+1
            jj=j-ndel
            do i = nguard,nxb+nguard+1
               ii=i-ndel

c            write(*,*)'i,j,k,lb=',i,j,k,lb
c            pause
c       if(mype.eq.0.and.lb.eq.30)then
c                write(*,*)'mype,lb,work1=',mype,lb,work1
c       endif
c       if(mype.eq.1.and.lb.eq.42)then
c                write(*,*)'mype,lb,work1=',mype,lb,work1
c       endif

!-----------------------------------------------------------------------
! d/dx                
            delu( 1,i,j,k) = work1(ii+1,jj,kk,1)-work1(ii-1,jj,kk,1)
            delu( 1,i,j,k) = delu(1,i,j,k)*delx
            delua(1,i,j,k) = abs(work1(ii+1,jj,kk,1))+
     &                       abs(work1(ii-1,jj,kk,1))
            delua(1,i,j,k) = delua(1,i,j,k)*delx
! d/dy
            delu( 2,i,j,k) = work1(ii,jj+1,kk,1)-work1(ii,jj-1,kk,1)
            delu( 2,i,j,k) = delu(2,i,j,k)*dely
            delua(2,i,j,k) = abs(work1(ii,jj+1,kk,1))+
     &                       abs(work1(ii,jj-1,kk,1))
            delua(2,i,j,k) = delua(2,i,j,k)*dely
! d/dz

#if N_DIM == 3
            delu( 3,i,j,k) = work1(ii,jj,kk+1,1)-work1(ii,jj,kk-1,1)
            delu( 3,i,j,k) = delu(3,i,j,k)*delz
            delua(3,i,j,k) = abs(work1(ii,jj,kk+1,1))+
     &                       abs(work1(ii,jj,kk-1,1))
            delua(3,i,j,k) = delua(3,i,j,k)*delz
!-----------------------------------------------------------------------
#endif

           end do
          end do
         end do

! Compute second derivatives

         do k = 1+(nguard*k3d),nzb+(nguard*k3d)
          do j = 1+nguard,nyb+nguard
           do i = 1+nguard,nxb+nguard

            do kk = 1,ndim2
             delu2(kk,i,j,k) = 0.
             delu3(kk,i,j,k) = 0.
             delu4(kk,i,j,k) = 0.
            end do
! d/dxdx
            delu2(1,i,j,k) = delu(1,i+1,j,k)-delu(1,i-1,j,k)
            delu2(1,i,j,k) = delu2(1,i,j,k)*delx
            delu3(1,i,j,k) = abs(delu(1,i+1,j,k))+abs(delu(1,i-1,j,k))
            delu3(1,i,j,k) = delu3(1,i,j,k)*delx
            delu4(1,i,j,k) = delua(1,i+1,j,k)+delua(1,i-1,j,k)
            delu4(1,i,j,k) = delu4(1,i,j,k)*delx
! d/dydx
            delu2(2,i,j,k) = delu(1,i,j+1,k)-delu(1,i,j-1,k)
            delu2(2,i,j,k) = delu2(2,i,j,k)*dely
            delu3(2,i,j,k) = abs(delu(1,i,j+1,k))+abs(delu(1,i,j-1,k))
            delu3(2,i,j,k) = delu3(2,i,j,k)*dely
            delu4(2,i,j,k) = delua(1,i,j+1,k)+delua(1,i,j-1,k)
            delu4(2,i,j,k) = delu4(2,i,j,k)*dely
! d/dxdy
            delu2(3,i,j,k) = delu(2,i+1,j,k)-delu(2,i-1,j,k)
            delu2(3,i,j,k) = delu2(3,i,j,k)*delx
            delu3(3,i,j,k) = abs(delu(2,i+1,j,k))+abs(delu(2,i-1,j,k))
            delu3(3,i,j,k) = delu3(3,i,j,k)*delx
            delu4(3,i,j,k) = delua(2,i+1,j,k)+delua(2,i-1,j,k)
            delu4(3,i,j,k) = delu4(3,i,j,k)*delx
! d/dydy
            delu2(4,i,j,k) = delu(2,i,j+1,k)-delu(2,i,j-1,k)
            delu2(4,i,j,k) = delu2(4,i,j,k)*dely
            delu3(4,i,j,k) = abs(delu(2,i,j+1,k))+abs(delu(2,i,j-1,k))
            delu3(4,i,j,k) = delu3(4,i,j,k)*dely
            delu4(4,i,j,k) = delua(2,i,j+1,k)+delua(2,i,j-1,k)
            delu4(4,i,j,k) = delu4(4,i,j,k)*dely
! d/dzdx

#if N_DIM == 3
            delu2(5,i,j,k) = delu(1,i,j,k+1)-delu(1,i,j,k-1)
            delu2(5,i,j,k) = delu2(5,i,j,k)*delz
            delu3(5,i,j,k) = abs(delu(1,i,j,k+1))+abs(delu(1,i,j,k-1))
            delu3(5,i,j,k) = delu3(5,i,j,k)*delz
            delu4(5,i,j,k) = delua(1,i,j,k+1)+delua(1,i,j,k-1)
            delu4(5,i,j,k) = delu4(5,i,j,k)*delz
! d/dzdy
            delu2(6,i,j,k) = delu(2,i,j,k+1)-delu(2,i,j,k-1)
            delu2(6,i,j,k) = delu2(6,i,j,k)*delz
            delu3(6,i,j,k) = abs(delu(2,i,j,k+1))+abs(delu(2,i,j,k-1))
            delu3(6,i,j,k) = delu3(6,i,j,k)*delz
            delu4(6,i,j,k) = delua(2,i,j,k+1)+delua(2,i,j,k-1)
            delu4(6,i,j,k) = delu4(6,i,j,k)*delz
! d/dxdz
            delu2(7,i,j,k) = delu(3,i+1,j,k)-delu(3,i-1,j,k)
            delu2(7,i,j,k) = delu2(7,i,j,k)*delx
            delu3(7,i,j,k) = abs(delu(3,i+1,j,k))+abs(delu(3,i-1,j,k))
            delu3(7,i,j,k) = delu3(7,i,j,k)*delx
            delu4(7,i,j,k) = delua(3,i+1,j,k)+delua(3,i-1,j,k)
            delu4(7,i,j,k) = delu4(7,i,j,k)*delx
! d/dydz
            delu2(8,i,j,k) = delu(3,i,j+1,k)-delu(3,i,j-1,k)
            delu2(8,i,j,k) = delu2(8,i,j,k)*dely
            delu3(8,i,j,k) = abs(delu(3,i,j+1,k))+abs(delu(3,i,j-1,k))
            delu3(8,i,j,k) = delu3(8,i,j,k)*dely
            delu4(8,i,j,k) = delua(3,i,j+1,k)+delua(3,i,j-1,k)
            delu4(8,i,j,k) = delu4(8,i,j,k)*dely
! d/dzdz
            delu2(9,i,j,k) = delu(3,i,j,k+1)-delu(3,i,j,k-1)
            delu2(9,i,j,k) = delu2(9,i,j,k)*delz
            delu3(9,i,j,k) = abs(delu(3,i,j,k+1))+abs(delu(3,i,j,k-1))
            delu3(9,i,j,k) = delu3(9,i,j,k)*delz
            delu4(9,i,j,k) = delua(3,i,j,k+1)+delua(3,i,j,k-1)
            delu4(9,i,j,k) = delu4(9,i,j,k)*delz
!-----------------------------------------------------------------------

#endif

               end do
            end do
         end do

! COMPUTE ERROR 
c         error = 0.
c         write(*,*)'error loop:error,epsil,lb=',error,epsil,lb
         do k = 1+(nguard*k3d),nzb+(nguard*k3d)
            do j = 1+nguard,nyb+nguard
               do i = 1+nguard,nxb+nguard

                  errort = 0.
                  
                  do kk = 1,ndim2
                     num    = delu2(kk,i,j,k)**2
                     denom  = delu3(kk,i,j,k)+(epsil*delu4(kk,i,j,k))
                     denom  = denom*denom
                     errort = errort+(num/(denom+1.e-20))
                  end do
c                  if(mype.eq.0.and.lb.eq.30)then
c                   write(*,*)'i,j,lb,mype,errort =',i,j,lb,mype,errort
c                  endif
c                  if(mype.eq.1.and.lb.eq.42)then
c                   write(*,*)'i,j,lb,mype,errort =',i,j,lb,mype,errort
c                  endif

                  error = max(error,sqrt(abs(errort)))

               end do
            end do
         end do
c         if(lb.eq.1)stop

 56      continue



      return
      end subroutine error_measure
