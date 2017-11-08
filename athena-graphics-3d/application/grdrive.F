#include "paramesh_preprocessor.fh"
#include "athena_preprocessor.fh"
      subroutine grdrive(nprocs,mype,vnk,varmin,varmax,ebox,
     $               istep,lgraph,nstep_graphics,lmovie,nmovie,movie_si)
#ifdef INTERACTIVE_GRAPHICS
!
! Subroutine to parallel graphics using AMR package.
! Writes out tree data structure and data stored in blocks
! Currently graphing are done serially by processor 0 collecting data from
! other processors and then writing it out.
!
! Written: K. Olson and C. Mobarry 9/97
      use physcons
      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace
      use strongfields
      use interior_gbc
      use athena_specific_parameters

      use paramesh_interfaces
#ifdef MPI_USED
      use paramesh_mpi_interfaces
#endif /* MPI_USED */


      implicit none
#ifdef MPI_USED
      include "mpif.h"
      integer status(MPI_STATUS_SIZE)
#endif /* MPI_USED */

      real,dimension(ne,il_bnd:iu_bnd,
     &                  jl_bnd:ju_bnd,
     &                  kl_bnd:ku_bnd,maxblocks)::vnk
      real,dimension(ne,il_bnd:iu_bnd,
     &                  jl_bnd:ju_bnd,
     &                  kl_bnd:ku_bnd,maxblocks)::vnk_temp

      real varmax(ne),varmin(ne)
      real*4 bbox(6,maxblocks_tr)
      real*4 bboxc(6)
      integer jproc,nprocs,mype
      integer ix,iy,iz,ixo,iyo,izo
      integer nodetypet(maxblocks_tr)
      integer lnblockst
      integer :: i,j,k,ib, iflagblock, iblocks
      integer icomm, idest, itag, isource, num, ierr 

      integer mbblocks
      parameter(mbblocks=108)
      integer ibmask(iu_bnd,ju_bnd,ku_bnd,mbblocks)
      integer ibflags(mbblocks+1), n, ibb

      real ebox(3,2)
      
      real,dimension(mdim,maxblocks_tr)             ::coordt,sizet
      real,dimension(ne,il_bnd:iu_bnd,jl_bnd:ju_bnd,
     &                                kl_bnd:ku_bnd)::unkf

      save coordt,sizet,nodetypet,
     .     lnblockst

      save unkf

      real*4,dimension(il_bnd:iu_bnd,jl_bnd:ju_bnd,
     &                               kl_bnd:ku_bnd)::wrk
      real*4 wrkmin,wrkmax

      integer iv,iunitc,igblock,i_d_processor, i_continue
      integer mv1,mx1,my1,nx1,ny1,ng0

      character*41 my_frame
      integer istep,nmovie,movie_si, nstep_graphics
      logical lmovie, lgraph
      logical draw_to_screen, save_to_file
      character*2 name1
      character*3 name2
      character*4 name3
      character*5 name4
      integer ig_draw
      real*4 bbox_ig(6)

      parameter(iunitc=10)
      logical logarithm,primitives,blockoutline,totalb,
     $     xcutflag,ycutflag,zcutflag
      logical xview, yview, zview 
      real*4
     $     cutplanex,cutplaney,cutplanez,fac_n,
     $     xmin,ymin,zmin,
     $     xmax,ymax,zmax,
     $     xtrans,ytrans,ztrans,
     $     xscale,yscale,zscale,
     $     anglex, angley, anglez,
     $     xaxis,yaxis,zaxis
      save
     $     iv,
     $     logarithm, primitives,blockoutline,totalb,
     $     wrkmin, wrkmax,
     $     cutplanex,cutplaney,cutplanez,fac_n,
     $     xcutflag,ycutflag,zcutflag,
     $     xmin,ymin,zmin,
     $     xmax,ymax,zmax,
     $     xtrans,ytrans,ztrans,
     $     xscale,yscale,zscale,
     $     anglex, angley, anglez,
     $     xview, yview, zview,
     $     xaxis,yaxis,zaxis

      real
     $     w_xmin,w_ymin,w_zmin,
     $     w_xmax,w_ymax,w_zmax
#if N_DIM == 2
!-----------------------
! pmn^s fieldline tracer
#include "fieldline_tracer.fh"
        real*4 :: xxx(npt_max)
        real*4 :: yyy(npt_max)
        integer :: jline,isign_f
        integer :: my_isign(nflines),my_f_npts(nflines)
! pmn^s fieldline tracer
!-----------------------
#endif
      integer :: ii,jj,kk
      integer :: errcode
      integer :: idothis
      integer :: iplane, jplane, kplane, ijkplane
      integer :: it, nv1, nv2, nv3
      logical :: plot_block
      real*4 xp_loc,yp_loc,zp_loc,dx1,dy1,dz1
      integer :: mb, i0, j0, k0
!-----------------------------------------------------------------------


      logical flagblock,blines

      DATA iv,igblock,i_d_processor, i_continue /1,0,0,0/
      DATA logarithm,primitives,blockoutline,totalb,blines,
     $     xcutflag,ycutflag,zcutflag 
     $     / .false.,.true.,.true.,.true.,
     $     .false.,.false.,.false.,.false. /
      DATA wrkmin,wrkmax /0.,11./
      DATA xmin,ymin,zmin,xmax,ymax,zmax /0.,1.,0.,1.,0.,1./
      DATA xtrans,ytrans,ztrans /0.,0.,0./
      DATA xscale,yscale,zscale /1.,1.,1./
      DATA anglex, angley, anglez /0.,0.,0./
      DATA xaxis,yaxis,zaxis /0.,0.,0./

      namelist /QUKPARMS/ 
     $     iv,igblock,i_d_processor, i_continue,
     $     logarithm, primitives,blockoutline,totalb,blines,
     $     wrkmin, wrkmax,
     $     xcutflag,ycutflag,zcutflag,
     $     cutplanex,cutplaney,cutplanez,fac_n,
     $     xmin,ymin,zmin,
     $     xmax,ymax,zmax,
     $     xtrans,ytrans,ztrans,
     $     xscale,yscale,zscale,
     $     anglex, angley, anglez,
     $     xview, yview, zview,
     $     xaxis,yaxis,zaxis

! compute n_to_left for each processor
      ng0 = nguard*npgs
#ifdef MPI_USED
        call mpi_barrier (MPI_COMM_WORLD, errcode)
#else
        call shmem_barrier_all()
#endif /*MPI_USED*/
!save solution
       do mb=1,maxblocks
          do k0=kl_bnd,ku_bnd
             do j0=jl_bnd,ju_bnd
                do i0=il_bnd,iu_bnd
                   do n=1,ne
                      vnk_temp(n,i0,j0,k0,mb)=vnk(n,i0,j0,k0,mb)
                   end do
                end do
             end do
          end do
       end do
 101   continue

!copy solution
       do mb=1,maxblocks
          do k0=kl_bnd,ku_bnd
             do j0=jl_bnd,ju_bnd
                do i0=il_bnd,iu_bnd
                   do n=1,ne
                      vnk(n,i0,j0,k0,mb)=vnk_temp(n,i0,j0,k0,mb)
                   end do
                end do
             end do
          end do
       end do

      ig_draw=0
      if (mype .eq. 0) then
           open(unit=iunitc,form='formatted',
     $        action='read',file='qukparms')
           read(iunitc,nml=QUKPARMS,ERR=135) 
           close(iunitc)
           !write(6,nml=QUKPARMS,ERR=135)
 135     continue
         w_xmin=xmin; w_ymin=ymin; w_zmin=zmin;
         w_xmax=xmax; w_ymax=ymax; w_zmax=zmax;

!Clear the graphics window here.
         if(xview) then
            anglex=0.; angley=90.; anglez=0.;
            xcutflag=.true.;ycutflag=.false.;zcutflag=.false.;
         endif
         if(yview) then
            anglex=-90.; angley=0.; anglez=0.;
            xcutflag=.false.;ycutflag=.true.;zcutflag=.false.;
         endif
         if(zview) then
            anglex=0.; angley=0.; anglez=0.;
            xcutflag=.false.;ycutflag=.false.;zcutflag=.true.;
         endif

         call qukframe1(
     $        xmin,ymin,zmin,
     $        xmax,ymax,zmax,
     $        xtrans,ytrans,ztrans,
     $        xscale,yscale,zscale,
     $        anglex,angley,anglez,
     $        xaxis,yaxis,zaxis
     $        )
      endif
!mk5do jproc = 0,nprocs-1
    
!        call qukframe8(
!    $        xmin,ymin,zmin,
!    $        xmax,ymax,zmax)
!     call qukframe4
            


!--------- save data to be displayed into vnk array --------------
!--note: original vnk original data is lost

       do ib = 1,lnblocks
       if(.not.primitives) then
          if(totalb)then
                call strongfields_wb(mype,ib)
                do k=1,nzb
                   do j=1,nyb
                      do i=1,nxb
                         unkf(1:5,i,j,k)=unk(1:5,i,j,k,ib)
                      enddo
                   enddo
                enddo
                do k=1,nzb
                   do j=1,nyb
                      do i=1,nxb
                      unkf(6:8,i,j,k)=unk(6:8,i,j,k,ib)+b_sf(1:3,i,j,k)  
                      enddo
                   enddo
                enddo
          else
                do k=1,nzb
                   do j=1,nyb
                      do i=1,nxb
                         unkf(1:ne,i,j,k)=unk(1:ne,i,j,k,ib)  
                      enddo
                   enddo
                enddo
          endif
       else
          if(totalb)then
                call strongfields_wb(mype,ib)
                do k=1,nzb
                   do j=1,nyb
                      do i=1,nxb
                         unkf(1:5,i,j,k)=vnk(1:5,i,j,k,ib)  
                      enddo
                   enddo
                enddo
                do k=1,nzb
                   kk = nguard*k3d+k
                   do j=1,nyb
                      jj = nguard*k2d+j
                      do i=1,nxb
                         ii = nguard+i
                   unkf(6:8,i,j,k)=vnk(6:8,i,j,k,ib)+b_sf(1:3,ii,jj,kk)
                      enddo
                   enddo
                enddo
          else
                do k=1,nzb
                   do j=1,nyb
                      do i=1,nxb
                         unkf(1:ne,i,j,k)=vnk(1:ne,i,j,k,ib)  
                      enddo
                   enddo
                enddo
          endif
       endif
!save unkf for block ib in vnk
       do k=1,nzb
          do j=1,nyb
             do i=1,nxb
                   vnk(1:ne,i,j,k,ib) = unkf(1:ne,i,j,k)
             end do
          end do
       end do
       end do  ! ib loop
!--------- end save data into vnk array ---------------------


!------------------------------------------------!
!MPI version, send data from each processor,     |
!receive and display in root                     |
!   send:                                        |
!           (lnblockst,lnblocks,1)               |
!           (coordt,coord,mdim*lnblockst)        |
!           (sizet,bsize,mdim*lnblockst)         |
!           (nodetypet,nodetype,lnblockst)       |
!           (vnk,  vnk, nxb*nyb*nzb*lnblockst)   |
!------------------------------------------------!
       !number of blocks that have interior boundary
       ! ibmask(:,:,:,:)=1
       ibflags(:)=0; n=0;
       do ib = 1,lnblocks
          if(bflags(1,ib).eq.1) then 
              n=n+1
              ibflags(n+1)=ib
          endif
       end do
       ibflags(1)=n
       if(n.gt.mbblocks) then
           write(*,*)'ERROR: grdrive: pl. increase mbblocks to ',n
           call amr_close
           stop
       endif

       !form ibmask for block interior
       do ibb = 1,n
          ib=ibflags(ibb+1)
          do k=1,nzb
             do j=1,nyb
                do i=1,nxb
                   ibmask(i,j,k,ibb)
     &               =bmask(i+nguard,j+nguard,k+nguard*k3d,ib)
                end do
             end do
          end do
       end do
#ifdef MPI_USED
       icomm = MPI_COMM_WORLD
       if(mype.ge.1) then
!..send data
!set for a max of 200 messages
            idest=0
            itag=5000+200*(mype-1)+1; num=1;
            call MPI_SEND(lnblocks,num,MPI_INTEGER,idest,
     .                 itag,icomm,ierr)
            itag=5000+200*(mype-1)+2; num=mdim*lnblocks;
            call MPI_SEND(coord(1,1),num,MPI_DOUBLE_PRECISION,idest,
     .                 itag,icomm,ierr)
            itag=5000+200*(mype-1)+3; num=mdim*lnblocks;
            call MPI_SEND(bsize(1,1),num,MPI_DOUBLE_PRECISION,idest,
     .                 itag,icomm,ierr)
            itag=5000+200*(mype-1)+4; num=lnblocks;
            call MPI_SEND(nodetype(1),num,MPI_INTEGER,idest,
     .                 itag,icomm,ierr)

            itag=5000+200*(mype-1)+5; 

            num = (iu_bnd-il_bnd+1)*(ju_bnd-jl_bnd+1)*(ku_bnd-kl_bnd+1)
            num=num*ne*lnblocks;
            call MPI_SEND(vnk(1,1,1,1,1),num,MPI_DOUBLE_PRECISION,idest,
     .                 itag,icomm,ierr)

            itag=5000+200*(mype-1)+6; num=mbblocks+1;
            call MPI_SEND(ibflags(1),num,MPI_INTEGER,idest,
     .                 itag,icomm,ierr)
           
            if(ibflags(1).gt.0) then
            !send ibmask
            itag=5000+200*(mype-1)+7; num=nxb*nyb*nzb*ibflags(1);
            call MPI_SEND(ibmask(1,1,1,1),num,MPI_INTEGER,idest,
     .                 itag,icomm,ierr)
            endif

       endif
#endif /*MPI_USED*/
!-----------------------------------------------------------------
       if(mype.eq.0) then
          iblocks=0
          do jproc = 0,nprocs-1
             if(jproc.eq.0) then
               lnblockst = lnblocks
               coordt(:,:)=coord(:,:)
               sizet(:,:)=bsize(:,:)
               nodetypet(:)=nodetype(:)
             else 
#ifdef MPI_USED
               isource=jproc
               itag=5000+200*(isource-1)+1; num=1;
               call MPI_RECV(lnblockst,num,MPI_INTEGER,
     .                isource,itag,icomm,status,ierr)
               itag=5000+200*(isource-1)+2; num=mdim*lnblockst;
               call MPI_RECV(coordt(1,1),num,MPI_DOUBLE_PRECISION,
     .                isource,itag,icomm,status,ierr)
               itag=5000+200*(isource-1)+3; num=mdim*lnblockst;
               call MPI_RECV(sizet(1,1),num,MPI_DOUBLE_PRECISION,
     .                isource,itag,icomm,status,ierr)
               itag=5000+200*(isource-1)+4; num=lnblockst;
               call MPI_RECV(nodetypet(1),num,MPI_INTEGER,
     .                isource,itag,icomm,status,ierr)

               itag=5000+200*(isource-1)+5;

               num=(iu_bnd-il_bnd+1)*(ju_bnd-jl_bnd+1)*(ku_bnd-kl_bnd+1)
               num=num*ne*lnblockst;
               call MPI_RECV(vnk(1,1,1,1,1),num,MPI_DOUBLE_PRECISION,
     .                isource,itag,icomm,status,ierr)

               itag=5000+200*(isource-1)+6; num=mbblocks+1;
               call MPI_RECV(ibflags(1),num,MPI_INTEGER,
     .                isource,itag,icomm,status,ierr)
            
               if(ibflags(1).gt.0) then
                  !receive ibmask
                itag=5000+200*(isource-1)+7; num=nxb*nyb*nzb*ibflags(1);
                  call MPI_RECV(ibmask(1,1,1,1),num,MPI_INTEGER,
     .                isource,itag,icomm,status,ierr)
               endif
#endif /*MPI_USED*/

#ifdef SHMEM_USED
            call SHMEM_INTEGER_GET (lnblockst,lnblocks,1,jproc)
            call SHMEM_REAL_GET(coordt,coord,mdim*lnblockst,jproc)
            call SHMEM_REAL_GET(sizet,bsize,mdim*lnblockst,jproc)
            call SHMEM_INTEGER_GET(nodetypet,nodetype,lnblockst,jproc)
            num=(iu_bnd-il_bnd+1)*(ju_bnd-jl_bnd+1)*(ku_bnd-kl_bnd+1)
            num=num*ne*lnblockst
            call SHMEM_REAL_GET(vnk(1,1,1,1,1),vnk(1,1,1,1,1),
     &              num,jproc)
#endif /*SHMEM_USED*/

             endif   ! jproc.ne.0

             idothis=1
             if(idothis.eq.1) then
!        ----------------------------------------------------
            do ib = 1,lnblockst
               bbox(1,ib) = coordt(1,ib) - 0.5*sizet(1,ib)
               bbox(2,ib) = coordt(1,ib) + 0.5*sizet(1,ib)
               bbox(3,ib) = coordt(2,ib) - 0.5*sizet(2,ib)
               bbox(4,ib) = coordt(2,ib) + 0.5*sizet(2,ib)
#if N_DIM == 3
               bbox(5,ib) = coordt(3,ib) - 0.5*sizet(3,ib)
               bbox(6,ib) = coordt(3,ib) + 0.5*sizet(3,ib)
#endif
            enddo

            do ib = 1,lnblockst
               do k=1,nzb
                  do j=1,nyb
                     do i=1,nxb
                        ! copy block data
                        unkf(1:ne,i,j,k)= vnk(1:ne,i,j,k,ib)
                     end do
                  end do
               end do

               flagblock = nodetypet(ib).eq.1
               iflagblock=0; if(flagblock)iflagblock=1;

        if(flagblock) then
            iblocks=iblocks+1

! 2D plots:
#if N_DIM == 2
            ixo=1+ng0; iyo=1+ng0; izo=(1+nzb+ng0*2*k3d)/2;
            bboxc(1)  = bbox(1,ib); bboxc(2)  = bbox(2,ib);
            bboxc(3)  = bbox(3,ib); bboxc(4)  = bbox(4,ib);
            mv1= 1; mx1= (nxb+ng0*2); my1= (nyb+ng0*2);
            nv1= 1; nx1= nxb; ny1= nyb;

            if(iv .le. ne .and. iv .ge.1) then
                    wrk   = unkf(iv,:,:,:)
            endif         ! if(iv .le. nvar .and. iv .ge.1)

               if(logarithm) then
                  wrk = alog10(abs(wrk)+1.e-20)
               endif
               wrk = max(wrkmin,wrk)
               wrk = min(wrkmax,wrk)

               call qukframe2(
     $              1,
     $              iflagblock, 
     $              bboxc,
     $              mv1,mx1,my1,
     $              nv1,nx1,ny1,
     $              wrk(ixo,iyo,izo),
     $              wrkmin,wrkmax)
               

#endif
! 3D plots:
#if N_DIM == 3
            ixo=1; iyo=1; izo=1; nv1=nxb; nv2=nyb; nv3=nzb;
            do it = 1, 6
               bboxc(it) = bbox(it,ib)
            end do
            if(cutplanex.eq.bboxc(1).or.
     $         cutplanex.eq.bboxc(2))cutplanex=cutplanex+0.0001
            if(cutplaney.eq.bboxc(3).or.
     $         cutplaney.eq.bboxc(4))cutplaney=cutplaney+0.0001
            if(cutplanez.eq.bboxc(5).or.
     $         cutplanez.eq.bboxc(6))cutplanez=cutplanez+0.0001

            do ijkplane=1,3
               iplane=0; jplane=0; kplane=0; plot_block=.false.;

               if(ijkplane.eq.1.and.xcutflag) then
                  ! z-y plane, x-plane
                  xp_loc = cutplanex
                  if(xp_loc.ge.bboxc(1).and.xp_loc.le.bboxc(2)) then
                      plot_block=.true.     
                      !location of x-plane 
                      dx1=(bboxc(2)-bboxc(1))/nxb;
                      n = (xp_loc-bboxc(1))/dx1;
                      iplane=n+1; if(iplane.gt.nxb)iplane=nxb;
            
                     if(iv .le. ne .and. iv .ge.1) then
                       do iz=1,nzb
                          do iy=1,nyb
                             wrk(iy,iz,1)=unkf(iv,iplane,iy,iz)
                        wrk(iy,iz,1) = max(wrkmin,wrk(iy,iz,1))
                        wrk(iy,iz,1) = min(wrkmax,wrk(iy,iz,1))
                          end do
                       end do
                     endif         ! if(iv .le. nvar .and. iv .ge.1)

                  endif    ! xp_loc
               endif   ! xcutflag

               if( ijkplane.eq.2.and.ycutflag) then
                  ! x-z plane, y-plane
                  yp_loc = cutplaney
                  if(yp_loc.ge.bboxc(3).and.yp_loc.le.bboxc(4)) then
                     plot_block=.true.     
                     !location of y-plane 
                     dy1=(bboxc(4)-bboxc(3))/nyb;
                     n = (yp_loc-bboxc(3))/dy1;
                     jplane=n+1; if(jplane.gt.nyb)jplane=nyb;
            
                     if(iv .le. ne .and. iv .ge.1) then
                       do iz=1,nzb
                          do ix=1,nxb
                             wrk(ix,iz,1)=unkf(iv,ix,jplane,iz)
                        wrk(ix,iz,1) = max(wrkmin,wrk(ix,iz,1))
                        wrk(ix,iz,1) = min(wrkmax,wrk(ix,iz,1))
                          end do
                       end do
                     endif         ! if(iv .le. nvar .and. iv .ge.1)

                  endif    ! yp_loc
               endif   ! ycutflag

               if( ijkplane.eq.3.and.zcutflag) then
                  !x-y plane, z-plane
                  zp_loc = cutplanez
                  if(zp_loc.ge.bboxc(5).and.zp_loc.le.bboxc(6)) then
                     plot_block=.true.     
                     !location of z-plane 
                     dz1=(bboxc(6)-bboxc(5))/nzb;
                     n = (zp_loc-bboxc(5))/dz1;
                     kplane=n+1; if(kplane.gt.nzb)kplane=nzb;
            
                     if(iv .le. ne .and. iv .ge.1) then
                       do iy=1,nyb
                          do ix=1,nxb
                             wrk(ix,iy,1)=unkf(iv,ix,iy,kplane)
                        wrk(ix,iy,1) = max(wrkmin,wrk(ix,iy,1))
                        wrk(ix,iy,1) = min(wrkmax,wrk(ix,iy,1))
                          end do
                       end do
                     endif         ! if(iv .le. nvar .and. iv .ge.1)

                  endif    ! zp_loc
               endif   ! zcutflag


           !   if(logarithm) then
           !      wrk = alog10(abs(wrk)+1.e-20)
           !   endif

               if(plot_block)
     $         call qukframe2_3d(
     $              ijkplane,iplane,jplane,kplane,
     $              cutplanex,cutplaney,cutplanez,
     $              bboxc,
     $              nv1,nv2,nv3,
     $              wrk,
     $              wrkmin,wrkmax)

               if(plot_block.and.blockoutline)
     $         call qukframe3_3d(
     $              ijkplane,iplane,jplane,kplane,
     $              nv1,nv2,nv3,
     $              cutplanex,cutplaney,cutplanez,fac_n,
     $              bboxc)

            end do   ! ijkplane
#endif
            endif  ! flagblock

         end do                 ! ib

c      Plots the block boundaries         
! 2D plots:
#if N_DIM == 2
         if(blockoutline) then

            call qukframe3(
     $           lnblockst,
     $           nodetypet, 
     $           bbox)

         endif

! draw a black box around the cells representing the earth

         if(ibflags(1).gt.0) then
             do ibb = 1,ibflags(1)
                ib=ibflags(ibb+1) 
               call qukframe6(nxb,nyb,nzb,bbox(1,ib),ibmask(1,1,1,ibb))
             end do
         endif

         if( (i_d_processor-1).eq.jproc.and.igblock.gt.0) then
             ! save data
              if(igblock.le.lnblockst) then
                   bbox_ig(1)=bbox(1,igblock) 
                   bbox_ig(2)=bbox(2,igblock)
                   bbox_ig(3)=bbox(3,igblock)
                   bbox_ig(4)=bbox(4,igblock)
                   bbox_ig(5)=bbox(5,igblock)
                   bbox_ig(6)=bbox(6,igblock)
                   ig_draw=1
              endif
         endif
#endif  !N_DIM == 2

! 3D plots:
#if N_DIM == 3
! draw a black box around the cells representing the earth

         if(ibflags(1).gt.0) then
             do ibb = 1,ibflags(1)
                ib=ibflags(ibb+1) 
               call qukframe6_3d(nxb,nyb,nzb,bbox(1,ib),
     $              ibmask(1,1,1,ibb))
             end do
         endif

!        if( (i_d_processor-1).eq.jproc.and.igblock.gt.0) then
             ! save data
!             if(igblock.le.lnblockst) then
!                  bbox_ig(1)=bbox(1,igblock) 
!                  bbox_ig(2)=bbox(2,igblock)
!                  bbox_ig(3)=bbox(3,igblock)
!                  bbox_ig(4)=bbox(4,igblock)
!                  bbox_ig(5)=bbox(5,igblock)
!                  bbox_ig(6)=bbox(6,igblock)
!                  ig_draw=1
!             endif
!        endif
#endif  !N_DIM == 3

!fieldlines
!#if N_DIM == 2
!      if( blines) then
!        do jline = 1,nflines
!           xxx(1:fline_npts(jline)) = xf(1:fline_npts(jline),jline)
!           yyy(1:fline_npts(jline)) = yf(1:fline_npts(jline),jline)
!           if(mod(jline,2).eq.1) then
!                isign_f = 1
!           else
!                isign_f = -1
!           endif
!           call quk_fieldlines_2d(
!     .        fline_npts(jline),xxx,yyy,isign_f )
!        enddo
!      endif
!#endif /*N_DIM*/

             endif    !idothis
           end do   ! jproc
 

      if(ig_draw.eq.1) call qukframe5(bbox_ig(1))

         call qukframe8(
     $        xmin,ymin,zmin,
     $        xmax,ymax,zmax)
! Windows are numbered starting at 1.
      if (lmovie .or. lgraph) then
         if (mod(istep,nmovie).eq.0 .or.
     +       mod(istep,nstep_graphics) .eq. 0) then
            draw_to_screen = mod(istep,nstep_graphics) .eq. 0
            save_to_file   = mod(istep,nmovie).eq.0
            call qukframe4(1, draw_to_screen, save_to_file)
         endif
      endif
    

!         if(lmovie) then
!            if(mod(istep,nmovie).eq.0) then
!                 my_frame(1:8)="scrsave "
!                 my_frame(9:21)="file_0000.rgb"
!                 if(movie_si.ge.0.and.movie_si.lt.10) then
!                     write(name1,"(i1)") movie_si
!                     my_frame(17:17)=name1
!                 endif
!                 if(movie_si.gt.9.and.movie_si.lt.100) then
!                     write(name2,"(i2)") movie_si
!                     my_frame(16:17)=name2
!                 endif
!                 if(movie_si.gt.99.and.movie_si.lt.1000) then
!                     write(name3,"(i3)") movie_si
!                     my_frame(15:17)=name3
!                 endif
!                 if(movie_si.gt.999.and.movie_si.lt.10000) then
!                     write(name4,"(i4)") movie_si
!                     my_frame(14:17)=name4
!                 endif

!note: total length should be 39
!mk                 my_frame(22:39)=" 174 576 387 796 " // char(0)
!mk                my_frame(22:39)="  50 700 100 900 " // char(0)
!mk                my_frame(22:39)="  70 576 330 855 " // char(0)
                  my_frame(22:39)="  70 650 300 885 " // char(0)

              !     call system(my_frame)
!             movie_si = movie_si + 1
!            endif
!         endif

        endif   ! mype = 0

        call mpi_bcast(i_continue,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         if(i_continue.ge.1) go to 101
      return
#endif
      end
