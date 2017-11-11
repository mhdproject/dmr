!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#define RICK_DIVB

#include "paramesh_preprocessor.fh"


      subroutine amr_prolong_fc_divbconsist(mype,level,nfield)




!------------------------------------------------------------------------
!
! Like the routine amr_prolong_fc_consist this routine checks for
! existing neighbor to newly created child blocks, and where found,
! uses facevar data from the existing neighbor at the common block
! boundary, in place of interpolation from the new childs parent.
! In addition this routine makes this change while adjusting values
! immediately interior to the new face to guarantee that div B is
! kept at zero.
!
! Only 1 field constructed from facevar components can be modified 
! during this call.
! To select it set the values for nfield and i_divf_fc_vars.
! Note, the algorithm which modifies the field is not optimal. In particular
! it has a tendency to drive oscillations in the neighborhood of strong
! field gradients. It is our intenetion to improve this in the future.
!
! Written :     Peter MacNeice          April 1998
!------------------------------------------------------------------------


      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton

      use paramesh_mpi_interfaces, only : 
     .                       mpi_amr_get_remote_block_fvar

#include "amr_shmem.fh"

      integer, intent(in) ::  mype
      integer, intent(in) ::  level
      integer, intent(in) ::  nfield

!------------------------------------
! local variables

      integer remote_pe,remote_block
      logical cnewchild
      save    cnewchild
      integer :: ierrorcode,ierr

! local arrays
      real recvx(nbndvar,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,
     .       kl_bnd:ku_bnd)
      real recvy(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,
     .       kl_bnd:ku_bnd)
      real recvz(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,
     .       kl_bnd:ku_bnd+k3d)

      real tempx(jl_bnd:ju_bnd,kl_bnd:ku_bnd)
      real tempy(il_bnd:iu_bnd,kl_bnd:ku_bnd)
      real tempz(il_bnd:iu_bnd,jl_bnd:ju_bnd)

      real e1,e2,ea,eb,efact
      real b11,b12,b21,b22

      real :: dx,dy,dz
      real :: area11,area12,area21,area22
      real :: areax1,areax2,areay1,areay2,areaz1,areaz2

      real :: bsum,divbmax,divb
      integer :: isw,idvx,idvy,idvz,isg,jf,i,j,k
      integer :: idest,i_dest,j_dest,i_source,j_source
      logical :: anewchild(1)

      integer nguard0
      parameter(nguard0 = nguard*npgs)
      integer nd
      parameter(nd = nguard - nguard0)

!------------------------------------
!!!#if N_DIM > 1
      if (ndim > 1) then

      call shmem_barrier_all()

!------------------------------------
!
! set components of facevar to which div correction is to be applied.
      idvx = i_divf_fc_vars(1,nfield)
      idvy = i_divf_fc_vars(2,nfield)
      idvz = i_divf_fc_vars(3,nfield)
!
!------------------------------------



      isw = 0


      efact = 1.
      if(ndim.eq.2) efact = .5

! cycle through the grid blocks on this processor
      if(lnblocks.gt.0) then
      do isg = 1,lnblocks

! Is this a newly created leaf block ?
!      if(nodetype(isg).eq.1.and.newchild(isg)) then
      if(nodetype(isg).eq.1.and.newchild(isg)
     .   .and. lrefine(isg).eq.level) then

!
! get block geometry information
#ifdef CURVILINEAR
       call amr_block_geometry(isg,mype)
#else
       dx = (bnd_box(2,1,isg)-bnd_box(1,1,isg))/real(nxb)
       dy = 1.
       dz = 1.
       if(ndim.ge.2)
     .         dy = (bnd_box(2,2,isg)-bnd_box(1,2,isg))/real(nyb)
       if(ndim.eq.3)
     .         dz = (bnd_box(2,3,isg)-bnd_box(1,3,isg))/real(nzb)
       
#endif /* CURVILINEAR */

! Cycle over the blocks faces
       do jf = 1,nfaces

          remote_pe = neigh(2,jf,isg)
          remote_block  = neigh(1,jf,isg)

#ifdef MPI_USED
          if(remote_block.gt.0.and.remote_pe.ne.mype) then
            do iblk = strt_buffer,last_buffer
              if( (remote_pe.eq.laddress(2,iblk)) .and.
     .            (remote_block.eq.laddress(1,iblk)) ) then
                remote_pe = mype
                remote_block = iblk
              endif
            enddo
          endif
          if(remote_block.gt.0.and.remote_pe.ne.mype) then
            write(*,*) 'Error : amr_prolong_fc_divbconsist : ',
     .       ' pe ',mype,' current blk ',isg,
     .       ' Remote block ',remote_block,remote_pe,' not ',
     .       'found locally.'
            call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
          endif
#endif /* MPI_USED */

! Is the neighbor to this face a pre-existing block?
          cnewchild = .true.
          if(remote_block.gt.0) then
              call shmem_logical_get(anewchild,
     .                      newchild(remote_block),1,remote_pe)
              cnewchild =anewchild(1)
          endif
          
          if(.not.cnewchild) then

! If the neighbor block is pre-existing then get its facevar data on the
! shared block boundary.

            idest = isg

            if(jf.eq.1) then
               i_dest   = nguard0 + 1 + iface_off
               i_source = nxb+nguard0 + 1 -gc_off_x + iface_off
               tempx(:,:) = facevarx(idvx,i_dest,:,:,idest)

               call mpi_amr_get_remote_block_fvar(mype,
     .                     remote_pe,remote_block,1,
     .                     recvx,recvy,recvz,idest)
               facevarx(idvx,i_dest,:,:,idest) = 
     .                             recvx(idvx,i_source,:,:)

               do k=nguard0*k3d+1,nguard0*k3d + nzb-k3d,2
               do j=nguard0+1,nguard0 + nyb-1,2
#ifdef CURVILINEAR
                area11= cell_area1(1+nguard,j+nd*k2d,k+nd*k3d)
                area21= cell_area1(1+nguard,j+nd*k2d,k+(nd+1)*k3d)
                area12= cell_area1(1+nguard,j+(nd+1)*k2d,k+nd*k3d)
                area22= cell_area1(1+nguard,j+(nd+1)*k2d,k+(nd+1)*k3d)
#else
                area11 = dy*dz
                area21 = dy*dz
                area12 = dy*dz
                area22 = dy*dz
#endif /* CURVILINEAR */
                b11 = (facevarx(idvx,nguard0+1,j,k,idest)
     .                       -tempx(j,k))*area11
                b21 = (facevarx(idvx,nguard0+1,j,k+k3d,idest)
     .                       - tempx(j,k+k3d))*area21
                b12 = (facevarx(idvx,nguard0+1,j+1,k,idest)
     .                       - tempx(j+1,k))*area12
                b22 = (facevarx(idvx,nguard0+1,j+1,k+k3d,idest)
     .                       - tempx(j+1,k+k3d))*area22
                bsum = b11+b12+b21+b22

                call compute_evalues(b11,b12,b21,b22,ea,eb,e1,e2,isw)

#ifdef CURVILINEAR
                areay1= cell_area2(1+nguard,j+(nd+1)*k2d,k+nd*k3d) 
                areay2= cell_area2(1+nguard,j+(nd+1)*k2d,k+(nd+1)*k3d) 
                areaz1= cell_area3(1+nguard,j+k2d,k+(nd+1)*k3d)
                areaz2= cell_area3(1+nguard,j+(nd+1)*k2d,k+(nd+1)*k3d)
#else
                areay1= dx*dz
                areay2= dx*dz
                areaz1= dx*dy
                areaz2= dx*dy
#endif /* CURVILINEAR */
                 facevary(idvy,nguard0+1,j+1,k,idest) =
     .             facevary(idvy,nguard0+1,j+1,k,idest)
     .             - e1*efact/areay1
                 facevary(idvy,nguard0+1,j+1,k+k3d,idest) =
     .             facevary(idvy,nguard0+1,j+1,k+k3d,idest)
     .             - e2*efact/areay2
                 facevarz(idvz,nguard0+1,j,k+k3d,idest) =
     .             facevarz(idvz,nguard0+1,j,k+k3d,idest)
     .             + ea/areaz1
                 facevarz(idvz,nguard0+1,j+1,k+k3d,idest) =
     .             facevarz(idvz,nguard0+1,j+1,k+k3d,idest)
     .             + eb/areaz2

               enddo
               enddo


            elseif(jf.eq.2) then

               i_dest   = nxb+1+nguard0 + iface_off
               i_source = 1+nguard0+gc_off_x + iface_off
               tempx(:,:) = facevarx(idvx,i_dest,:,:,idest)
               call mpi_amr_get_remote_block_fvar(mype,
     .                     remote_pe,remote_block,1,
     .                     recvx,recvy,recvz,idest)
               facevarx(idvx,i_dest,:,:,idest) = 
     .                             recvx(idvx,i_source,:,:)

               do k=nguard0*k3d+1,nguard0*k3d + nzb-k3d,2
               do j=nguard0+1,nguard0 + nyb-1,2
#ifdef CURVILINEAR
                area11= cell_area1(nxb+1+nguard,j+nd*k2d,k+nd*k3d)
                area21= cell_area1(nxb+1+nguard,j+nd*k2d,k+(nd+1)*k3d)
                area12= cell_area1(nxb+1+nguard,j+(nd+1)*k2d,k+nd*k3d)
                area22= 
     .              cell_area1(nxb+1+nguard,j+(nd+1)*k2d,k+(nd+1)*k3d)
#else
                area11 = dy*dz
                area21 = dy*dz
                area12 = dy*dz
                area22 = dy*dz
#endif /* CURVILINEAR */

                b11 = (facevarx(idvx,nguard0+nxb+1,j,k,idest)
     .                         - tempx(j,k))*area11
                b21 = (facevarx(idvx,nguard0+nxb+1,j,k+k3d,idest)
     .                         - tempx(j,k+k3d))*area21
                b12 = (facevarx(idvx,nguard0+nxb+1,j+1,k,idest)
     .                         - tempx(j+1,k))*area12
                b22 = (facevarx(idvx,nguard0+nxb+1,j+1,k+k3d,idest)
     .                         - tempx(j+1,k+k3d))*area22
                bsum = b11+b12+b21+b22
                call compute_evalues(b11,b12,b21,b22,ea,eb,e1,e2,isw)

#ifdef CURVILINEAR
                areay1= cell_area2(nxb+nguard,j+(nd+1)*k2d,k+nd*k3d)
                areay2= cell_area2(nxb+nguard,j+(nd+1)*k2d,k+nd*k3d)
                areaz1= cell_area3(nxb+nguard,j+k2d,k+(nd+1)*k3d)
                areaz2= cell_area3(nxb+nguard,j+nd*k2d,k+(nd+1)*k3d)
#else
                areay1= dx*dz
                areay2= dx*dz
                areaz1= dx*dy
                areaz2= dx*dy
#endif /* CURVILINEAR */

                facevary(idvy,nguard0+nxb,j+1,k,idest) =
     .             facevary(idvy,nguard0+nxb,j+1,k,idest)
     .             + e1*efact/areay1
                facevary(idvy,nguard0+nxb,j+1,k+k3d,idest) =
     .             facevary(idvy,nguard0+nxb,j+1,k+k3d,idest)
     .             + e2*efact/areay2
                facevarz(idvz,nguard0+nxb,j,k+k3d,idest) =
     .             facevarz(idvz,nguard0+nxb,j,k+k3d,idest)
     .             - ea/areaz1
                facevarz(idvz,nguard0+nxb,j+1,k+k3d,idest) =
     .             facevarz(idvz,nguard0+nxb,j+1,k+k3d,idest)
     .             - eb/areaz2
               enddo
               enddo

            elseif(jf.eq.3) then

               j_dest   = nguard0*k2d + 1 + iface_off*k2d
               j_source = nyb+nguard0 + 1 -gc_off_y + iface_off
               tempy(:,:) = facevary(idvy,:,j_dest,:,idest)
               call mpi_amr_get_remote_block_fvar(mype,
     .                     remote_pe,remote_block,2,
     .                     recvx,recvy,recvz,idest)
               facevary(idvy,:,j_dest,:,idest) =
     .                                recvy(idvy,:,j_source,:)

               do k=nguard0*k3d+1,nguard0*k3d + nzb-k3d,2
               do i=nguard0+1,nguard0 + nxb-1,2
#ifdef CURVILINEAR
                area11= cell_area2(i+nd,1+nguard,k+nd*k3d)
                area21= cell_area2(i+nd,1+nguard,k+(nd+1)*k3d)
                area12= cell_area2(i+1+nd,1+nguard,k+nd*k3d)
                area22= cell_area2(i+1+nd,1+nguard,k+(nd+1)*k3d)
#else
                area11 = dx*dz
                area21 = dx*dz
                area12 = dx*dz
                area22 = dx*dz
#endif /* CURVILINEAR */

                b11 = (facevary(idvy,i,nguard0*k2d+1,k,idest) 
     .                 - tempy(i,k) )*area11
                b21 = (facevary(idvy,i,nguard0*k2d+1,k+k3d,idest)
     .                 - tempy(i,k+k3d) )*area21
                b12 = (facevary(idvy,i+1,nguard0*k2d+1,k,idest)
     .                 - tempy(i+1,k) )*area12
                b22 = (facevary(idvy,i+1,nguard0*k2d+1,k+k3d,idest)
     .                 - tempy(i+1,k+k3d) )*area22
                bsum = b11+b12+b21+b22
                call compute_evalues(b11,b12,b21,b22,ea,eb,e1,e2,isw)

#ifdef CURVILINEAR
                areax1= cell_area1(i+1+nd,1+nguard,k+nd*k3d)
                areax2= cell_area1(i+1+nd,1+nguard,k+nd*k3d)
                areaz1= cell_area3(i+nd,1+nguard,k+(nd+1)*k3d)
                areaz2= cell_area3(i+nd,1+nguard,k+(nd+1)*k3d)
#else
                areax1= dy*dz
                areax2= dy*dz
                areaz1= dx*dy
                areaz2= dx*dy
#endif /* CURVILINEAR */

                facevarx(idvx,i+1,nguard0*k2d+1,k,idest) =
     .             facevarx(idvx,i+1,nguard0*k2d+1,k,idest)
     .             - e1*efact/areax1
                facevarx(idvx,i+1,nguard0*k2d+1,k+k3d,idest) =
     .             facevarx(idvx,i+1,nguard0*k2d+1,k+k3d,idest)
     .             - e2*efact/areax2
                facevarz(idvz,i,nguard0*k2d+1,k+k3d,idest) =
     .             facevarz(idvz,i,nguard0*k2d+1,k+k3d,idest)
     .             + ea/areaz1
                facevarz(idvz,i+1,nguard0*k2d+1,k+k3d,idest) =
     .             facevarz(idvz,i+1,nguard0*k2d+1,k+k3d,idest)
     .             + eb/areaz2
               enddo
               enddo

            elseif(jf.eq.4) then

               j_dest   = nyb*k2d + 1 + nguard0*k2d + iface_off*k2d
               j_source = 1+nguard0+gc_off_y + iface_off
               tempy(:,:) = facevary(idvy,:,j_dest,:,idest)
               call mpi_amr_get_remote_block_fvar(mype,
     .                     remote_pe,remote_block,2,
     .                     recvx,recvy,recvz,idest)
               facevary(idvy,:,j_dest,:,idest) = 
     .                           recvy(idvy,:,j_source,:)

               do k=nguard0*k3d+1,nguard0*k3d + nzb-k3d,2
               do i=nguard0+1,nguard0 + nxb-1,2
#ifdef CURVILINEAR
                area11= cell_area2(i+nd,nyb+1+nguard,k+nd*k3d)
                area21= cell_area2(i+nd,nyb+1+nguard,k+(nd+1)*k3d)
                area12= cell_area2(i+1+nd,nyb+1+nguard,k+nd*k3d)
                area22= cell_area2(i+1+nd,nyb+1+nguard,k+(nd+1)*k3d)
#else
                area11 = dx*dz
                area21 = dx*dz
                area12 = dx*dz
                area22 = dx*dz
#endif /* CURVILINEAR */

                b11 = ( facevary(idvy,i,nguard0*k2d+nyb+k2d,k,idest)
     .                      - tempy(i,k) )*area11
                b21 = ( facevary(idvy,i,nguard0*k2d+nyb+k2d,k+k3d,idest)
     .                      - tempy(i,k+k3d) )*area21
                b12 = ( facevary(idvy,i+1,nguard0*k2d+nyb+k2d,k,idest)
     .                      - tempy(i+1,k) )*area12
                b22 = ( facevary(idvy,i+1,nguard0*k2d+nyb+k2d,
     .                        k+k3d,idest)
     .                      - tempy(i+1,k+k3d) )*area22
                bsum = b11+b12+b21+b22
                call compute_evalues(b11,b12,b21,b22,ea,eb,e1,e2,isw)

#ifdef CURVILINEAR
                areax1= cell_area1(i+1+nd,nyb+1+nguard,k+nd*k3d)
                areax2= cell_area1(i+1+nd,nyb+1+nguard,k+nd*k3d)
                areaz1= cell_area3(i+nd,nyb+1+nguard,k+(nd+1)*k3d)
                areaz2= cell_area3(i+nd,nyb+1+nguard,k+(nd+1)*k3d)
#else
                areax1= dy*dz
                areax2= dy*dz
                areaz1= dx*dy
                areaz2= dx*dy
#endif /* CURVILINEAR */

                facevarx(idvx,i+1,nguard0*k2d+nyb,k,idest) =
     .             facevarx(idvx,i+1,nguard0*k2d+nyb,k,idest)
     .                  + e1*efact/areax1
                facevarx(idvx,i+1,nguard0*k2d+nyb,k+k3d,idest) =
     .             facevarx(idvx,i+1,nguard0*k2d+nyb,k+k3d,idest)
     .                  + e2*efact/areax2
                facevarz(idvz,i,nguard0*k2d+nyb,k+k3d,idest) =
     .             facevarz(idvz,i,nguard0*k2d+nyb,k+k3d,idest)
     .                  - ea/areaz1
                facevarz(idvz,i+1,nguard0*k2d+nyb,k+k3d,idest) =
     .             facevarz(idvz,i+1,nguard0*k2d+nyb,k+k3d,idest)
     .                  - eb/areaz2
               enddo
               enddo

!!!#if N_DIM == 3
            elseif(jf.eq.5 .and. ndim == 3) then

               k_dest   = nguard0*k3d + 1 + iface_off*k3d
               k_source = nzb+nguard0 + 1 -gc_off_z + iface_off
               tempz(:,:) = facevarz(idvz,:,:,k_dest,idest)
               call mpi_amr_get_remote_block_fvar(mype,
     .                     remote_pe,remote_block,3,
     .                     recvx,recvy,recvz,idest)
               facevarz(idvz,:,:,k_dest,idest) =
     .                            recvz(idvz,:,:,k_source)

               do j=nguard0+1,nguard0 + nyb-1,2
               do i=nguard0+1,nguard0 + nxb-1,2
#ifdef CURVILINEAR
                area11= cell_area3(i+nd  ,j+nd  ,1+nguard*k3d)
                area21= cell_area3(i+nd  ,j+1+nd,1+nguard*k3d)
                area12= cell_area3(i+1+nd,j+nd  ,1+nguard*k3d)
                area22= cell_area3(i+1+nd,j+1+nd,1+nguard*k3d)
#else
                area11 = dx*dy
                area21 = dx*dy
                area12 = dx*dy
                area22 = dx*dy
#endif /* CURVILINEAR */

                b11 = (facevarz(idvz,i,j,nguard0*k3d+1,idest) 
     .                  - tempz(i,j) )*area11
                b21 = (facevarz(idvz,i,j+1,nguard0*k3d+1,idest)
     .                  - tempz(i,j+1) )*area21
                b12 = (facevarz(idvz,i+1,j,nguard0*k3d+1,idest)
     .                  - tempz(i+1,j) )*area12
                b22 = (facevarz(idvz,i+1,j+1,nguard0*k3d+1,idest)
     .                  - tempz(i+1,j+1) )*area22
                bsum = b11+b12+b21+b22
                call compute_evalues(b11,b12,b21,b22,ea,eb,e1,e2,isw)

#ifdef CURVILINEAR
                areax1= cell_area1(i+1+nd,j+nd  ,1+nguard*k3d)
                areax2= cell_area1(i+1+nd,j+1+nd,1+nguard*k3d)
                areay1= cell_area2(i+nd  ,j+1+nd,1+nguard*k3d)
                areay2= cell_area2(i+nd+1,j+1+nd,1+nguard*k3d)
#else
                areax1= dy*dz
                areax2= dy*dz
                areay1= dx*dz
                areay2= dx*dz
#endif /* CURVILINEAR */

                facevarx(idvx,i+1,j,nguard0*k3d+1,idest) =
     .             facevarx(idvx,i+1,j,nguard0*k3d+1,idest)
     .             - e1/areax1
                facevarx(idvx,i+1,j+1,nguard0*k3d+1,idest) =
     .             facevarx(idvx,i+1,j+1,nguard0*k3d+1,idest)
     .             - e2/areax2
                facevary(idvy,i,j+1,nguard0*k3d+1,idest) =
     .             facevary(idvy,i,j+1,nguard0*k3d+1,idest)
     .             + ea/areay1
                facevary(idvy,i+1,j+1,nguard0*k3d+1,idest) =
     .             facevary(idvy,i+1,j+1,nguard0*k3d+1,idest)
     .             + eb/areay2
               enddo
               enddo

            elseif(jf.eq.6 .and. ndim == 3) then

               k_dest   = nzb*k3d+1+nguard0*k3d + iface_off*k3d
               k_source = 1+nguard0+gc_off_z + iface_off
               tempz(:,:) = facevarz(idvz,:,:,k_dest,idest)
               call mpi_amr_get_remote_block_fvar(mype,
     .                     remote_pe,remote_block,3,
     .                     recvx,recvy,recvz,idest)
               facevarz(idvz,:,:,k_dest,idest) =
     .                             recvz(idvz,:,:,k_source)

               do j=nguard0+1,nguard0 + nyb -1,2
               do i=nguard0+1,nguard0 + nxb -1,2

#ifdef CURVILINEAR
                area11= cell_area3(i+nd  ,j+nd  ,nzb+(1+nguard)*k3d)
                area21= cell_area3(i+nd  ,j+1+nd,nzb+(1+nguard)*k3d)
                area12= cell_area3(i+1+nd,j+nd  ,nzb+(1+nguard)*k3d)
                area22= cell_area3(i+1+nd,j+1+nd,nzb+(1+nguard)*k3d)
#else
                area11 = dx*dy
                area21 = dx*dy
                area12 = dx*dy
                area22 = dx*dy
#endif /* CURVILINEAR */

                b11 = (facevarz(idvz,i,j,nzb+(nguard0+1)*k3d,idest)
     .                       - tempz(i,j) )*area11
                b21 = (facevarz(idvz,i,j+1,nzb+(nguard0+1)*k3d,idest)
     .                       - tempz(i,j+1) )*area21
                b12 = (facevarz(idvz,i+1,j,nzb+(nguard0+1)*k3d,idest)
     .                       - tempz(i+1,j) )*area12
                b22 = (facevarz(idvz,i+1,j+1,nzb+(nguard0+1)*k3d,idest)
     .                       - tempz(i+1,j+1) )*area22
                bsum = b11+b12+b21+b22
                call compute_evalues(b11,b12,b21,b22,ea,eb,e1,e2,isw)

#ifdef CURVILINEAR
                areax1= cell_area1(i+1+nd,j+nd  ,nzb+(1+nguard)*k3d)
                areax2= cell_area1(i+1+nd,j+1+nd,nzb+(1+nguard)*k3d)
                areay1= cell_area2(i+nd  ,j+1+nd,nzb+(1+nguard)*k3d)
                areay2= cell_area2(i+nd+1,j+1+nd,nzb+(1+nguard)*k3d)
#else
                areax1= dy*dz
                areax2= dy*dz
                areay1= dx*dz
                areay2= dx*dz
#endif /* CURVILINEAR */

                facevarx(idvx,i+1,j,nguard0*k3d+nzb,idest) =
     .             facevarx(idvx,i+1,j,nguard0*k3d+nzb,idest)
     .             + e1/areax1
                facevarx(idvx,i+1,j+1,nguard0*k3d+nzb,idest) =
     .             facevarx(idvx,i+1,j+1,nguard0*k3d+nzb,idest)
     .             + e2/areax2
                facevary(idvy,i,j+1,nguard0*k3d+nzb,idest) =
     .             facevary(idvy,i,j+1,nguard0*k3d+nzb,idest)
     .             - ea/areay1
                facevary(idvy,i+1,j+1,nguard0*k3d+nzb,idest) =
     .             facevary(idvy,i+1,j+1,nguard0*k3d+nzb,idest)
     .             - eb/areay2
               enddo
               enddo


!!!#endif                      /* N_DIM == 3 */
            endif

          endif

       enddo

      endif
      enddo
      endif

!
! div B check

      if(lnblocks.gt.0) then
      do isg = 1,lnblocks

! Is this a newly created leaf block ?
      if(nodetype(isg).eq.1) then

               divbmax = 0.
               do k=nguard0*k3d+1,nguard0*k3d + nzb
               do j=nguard0+1,nguard0 + nyb
               do i=nguard0+1,nguard0 + nxb
                 divb = (
     .            facevarx(1,i+1,j,k,isg) - facevarx(1,i,j,k,isg)
     .            + facevary(1,i,j+1,k,isg) - facevary(1,i,j,k,isg)
     .            + facevarz(1,i,j,k+k3d,isg) - facevarz(1,i,j,k,isg))
                divbmax = max(divbmax,abs(divb))
               enddo
               enddo
               enddo

!               write(*,*) 'Max div B on block ',isg,' is ',divbmax
      endif
      enddo
      endif

!!!#endif   /* N_DIM > 1 */
      end if

      return
      end subroutine amr_prolong_fc_divbconsist



      subroutine compute_evalues(b11,b12,b21,b22,ea,eb,e1,e2,isw)

!
! Computes virtual electric field values required to produce
! changes in facevar's in a manner which will preserve div = 0
! constraints. The b11, b12 etc input arguments specify the
! required change in the vector component on the chosen face.
! The ea,eb,e1,e2 arguments return the virtual electric field
! values required to achieve this adjustment.
!
!
! The relationship between b's and e's is
!
!               -------------------------------------
!              |                  |                  |
!              |                  |                  |
!              |                  ^                  |
!              |       b21        e2     b22         |
!              |                  |                  |
!              |                  |                  |
!              |------ ea ->------|----- eb ->-------|
!              |                  |                  |
!              |                  ^                  |
!              |       b11        e1     b12         |
!              |                  |                  |
!              |                  |                  |
!              |                  |                  |
!               -------------------------------------
!
! These electric fields must be applied to adjust the appropriate
! components of the vector field in the planes perpendicular to
! the chosen face, immediately inside this face.


      use paramesh_dimensions
      use physicaldata

      integer :: isw

! local variables

      real e1,e2,ea,eb
      real b11,b12,b21,b22

#ifdef ORIGINAL_DIVB
      ea = ( 2.*b11 +    b12 -    b21          )*.25
      eb = (    b11 + 2.*b12             - b22 )*.25
      e1 = (-2.*b11 +    b12 -    b21          )*.25
      e2 = (   -b11          - 2.*b21 + b22    )*.25
#endif /* ORIGINAL_DIVB */
#ifdef RICK_DIVB
      ea = (+ 5. * b11 + 3. * b12 - 1. * b21 + 1. * b22)*.125
      eb = (- 1. * b11 + 1. * b12 - 3. * b21 - 5. * b22)*.125
      e1 = (- 1. * b11 + 5. * b12 + 1. * b21 + 3. * b22)*.125
      e2 = (- 3. * b11 - 1. * b12 - 5. * b21 + 1. * b22)*.125
#endif /* RICK_DIVB */

      if(isw.eq.1) then
       write(*,*) 'compute : b ',b11,b12,b21,b22
       write(*,*) 'compute : e ',ea,eb,e1,e2
       isw = 0
      endif

      return
      end subroutine compute_evalues
