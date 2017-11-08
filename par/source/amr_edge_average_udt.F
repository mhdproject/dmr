!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"


      subroutine amr_edge_average_udt(mype)




!------------------------------------------------------------------------
!
! This routine gets cell edge-based data at block boundaries from 
! neighbors who are parents of leaf blocks. 
!
! The data structure used to store and pass this data is defined
! in the include file 'block_boundary_data.fh' which can be included
! in 'physicaldata.fh'.
!
! This version is called when uniform timesteps are being used across
! the blocks in the computation.
!
!
! Written :     Peter MacNeice          August 1997
!------------------------------------------------------------------------


      use paramesh_dimensions
      use physicaldata
      use tree

      use paramesh_interfaces, only : amr_restrict_edge_data


#include "amr_shmem.fh"

      integer, intent(in)  ::  mype

!------------------------------------
! local variables

      integer nguard0
      parameter(nguard0 = nguard*npgs)

      integer ng_off
      parameter(ng_off = nguard0+iface_off)

      integer kup,klo
      parameter(klo  = 1+k3d*nguard0)
      parameter(kup  = 1+k3d*(nzb+nguard0-1))
      parameter(kup1 = k3d+nzb+k3d*nguard0)


      integer remote_pe,remote_block
      integer cnodetype

      integer :: lb,jf

!------------------------------------


!!!#if N_DIM >= 2
      if (ndim >= 2) then


! all leaf blocks provide reduced boundary edge data to their parents
      call amr_restrict_edge_data(mype)

      call shmem_barrier_all()


! cycle through the grid blocks on this processor
      if(lnblocks.gt.0) then
      do lb = 1,lnblocks

! Is this a leaf block and not at the original refinement level ?
!      if(nodetype(lb).eq.1.and.lrefine(lb).gt.1) then
      if(nodetype(lb).eq.1) then


! Cycle over the blocks faces
        do jf = 1,nfaces

          remote_pe = neigh(2,jf,lb)
          remote_block  = neigh(1,jf,lb)

! Is the neighbor to this face a parent of a leaf block?
          cnodetype = 0
          if(remote_block.gt.0) call shmem_integer_get(cnodetype,
     .                         nodetype(remote_block),1,remote_pe)
          if(cnodetype.eq.2) then

! If yes then copy the appropriate layer from its boundary variable data 

            if(jf.eq.1) then

              call shmem_real_get(recvarx1e(1,1,1,1),
     .                       bedge_facex_y(1,1,1,1,remote_block),
     .                       len_block_ex*nedges,remote_pe)
              bedge_facex_y(:,1,:,:,lb) = recvarx1e(:,2,:,:)

              if((ndim.eq.3).or.(l2p5d.eq.1)) then
                call shmem_real_get(recvarx1e(1,1,1,1),
     .                        bedge_facex_z(1,1,1,1,remote_block),
     .                        len_block_ex*nedges,remote_pe)
                bedge_facex_z(:,1,:,:,lb) = recvarx1e(:,2,:,:)
              endif

! make common variables on an edge consistent
              bedge_facey_z(:,1+nguard0,1,klo:kup,lb) =
     .               bedge_facex_z(:,1,1+nguard0,klo:kup,lb)
       
              bedge_facey_z(:,1+nguard0,2,klo:kup,lb) =
     .               bedge_facex_z(:,1,1+nyb+nguard0,klo:kup,lb)

              if(ndim.eq.3) then
                bedge_facez_y(:,1+nguard0,1+nguard0:nyb+nguard0,1,lb)
     .             = bedge_facex_y(:,1,1+nguard0:nyb+nguard0,klo,lb)

                bedge_facez_y(:,1+nguard0,1+nguard0:nyb+nguard0,2,lb)
     .             = bedge_facex_y(:,1,1+nguard0:nyb+nguard0,kup1,lb)
              endif


            elseif(jf.eq.2) then

              call shmem_real_get(recvarx1e(1,1,1,1),
     .               bedge_facex_y(1,1,1,1,remote_block),
     .               len_block_ex*nedges,remote_pe)
              bedge_facex_y(:,2,:,:,lb)=recvarx1e(:,1,:,:)

              if((ndim.eq.3).or.(l2p5d.eq.1)) then
                call shmem_real_get(recvarx1e(1,1,1,1),
     .               bedge_facex_z(1,1,1,1,remote_block),
     .               len_block_ex*nedges,remote_pe)
                bedge_facex_z(:,2,:,:,lb)=recvarx1e(:,1,:,:)
              endif

! make common variables on an edge consistent
              bedge_facey_z(:,1+nxb+nguard0,1,klo:kup,lb) =
     .            bedge_facex_z(:,2,1+nguard0,klo:kup,lb)
                
              bedge_facey_z(:,1+nxb+nguard0,2,klo:kup,lb) =
     .            bedge_facex_z(:,2,1+nyb+nguard0,klo:kup,lb)

              if(ndim.eq.3) then
                bedge_facez_y(:,1+nxb+nguard0,1+nguard0:nyb+nguard0,
     .                                                       1,lb)=
     .            bedge_facex_y(:,2,1+nguard0:nyb+nguard0,klo,lb)

                bedge_facez_y(:,1+nxb+nguard0,1+nguard0:nyb+nguard0,
     .                                                       2,lb)=
     .            bedge_facex_y(:,2,1+nguard0:nyb+nguard0,kup1,lb)


              endif

            elseif(jf.eq.3) then

              if((ndim.eq.3).or.(l2p5d.eq.1)) then
                call shmem_real_get(recvary1e(1,1,1,1),
     .               bedge_facey_z(1,1,1,1,remote_block),
     .               len_block_ey*nedges,remote_pe)
                bedge_facey_z(:,:,1,:,lb) = recvary1e(:,:,2,:)
              endif

              call shmem_real_get(recvary1e(1,1,1,1),
     .               bedge_facey_x(1,1,1,1,remote_block),
     .               len_block_ey*nedges,remote_pe)
              bedge_facey_x(:,:,1,:,lb) = recvary1e(:,:,2,:)


! make common variables on an edge consistent
              bedge_facex_z(:,1,1+nguard0,klo:kup,lb) =
     .          bedge_facey_z(:,1+nguard0,1,klo:kup,lb)

              bedge_facex_z(:,2,1+nguard0,klo:kup,lb) =
     .          bedge_facey_z(:,1+nxb+nguard0,1,klo:kup,lb)

              if(ndim.eq.3) then
                bedge_facez_x(:,1+nguard0:nxb+nguard0,1+nguard0,1,lb)=
     .          bedge_facey_x(:,1+nguard0:nxb+nguard0,1,klo,lb)

                bedge_facez_x(:,1+nguard0:nxb+nguard0,1+nguard0,2,lb)=
     .          bedge_facey_x(:,1+nguard0:nxb+nguard0,1,kup1,lb)
                endif

              elseif(jf.eq.4) then

                if((ndim.eq.3).or.(l2p5d.eq.1)) then
                  call shmem_real_get(recvary1e(1,1,1,1),
     .               bedge_facey_z(1,1,1,1,remote_block),
     .               len_block_ey*nedges,remote_pe)
                  bedge_facey_z(:,:,2,:,lb) = recvary1e(:,:,1,:)
                endif

                call shmem_real_get(recvary1e(1,1,1,1),
     .               bedge_facey_x(1,1,1,1,remote_block),
     .               len_block_ey*nedges,remote_pe)
                bedge_facey_x(:,:,2,:,lb) = recvary1e(:,:,1,:)


! make common variables on an edge consistent
                bedge_facex_z(:,1,1+nyb+nguard0,klo:kup,lb) =
     .              bedge_facey_z(:,1+nguard0,2,klo:kup,lb)

                bedge_facex_z(:,2,1+nyb+nguard0,klo:kup,lb) =
     .              bedge_facey_z(:,1+nxb+nguard0,2,klo:kup,lb)
                if(ndim.eq.3) then
                bedge_facez_x(:,1+nguard0:nxb+nguard0,1+nyb+nguard0,
     .                                                       1,lb)=
     .              bedge_facey_x(:,1+nguard0:nxb+nguard0,2,klo,lb)

                bedge_facez_x(:,1+nguard0:nxb+nguard0,1+nyb+nguard0,
     .                                                       2,lb)=
     .              bedge_facey_x(:,1+nguard0:nxb+nguard0,2,kup1,lb)
                endif

              elseif(jf.eq.5) then

                call shmem_real_get(recvarz1e(1,1,1,1),
     .               bedge_facez_x(1,1,1,1,remote_block),
     .               len_block_ez*nedges,remote_pe)
                bedge_facez_x(:,:,:,1,lb) = recvarz1e(:,:,:,2)
                call shmem_real_get(recvarz1e(1,1,1,1),
     .               bedge_facez_y(1,1,1,1,remote_block),
     .               len_block_ez*nedges,remote_pe)
                bedge_facez_y(:,:,:,1,lb) = recvarz1e(:,:,:,2) 

! make common variables on an edge consistent
                bedge_facey_x(:,1+nguard0:nxb+nguard0,1,klo,lb)=
     .          bedge_facez_x(:,1+nguard0:nxb+nguard0,1+nguard0,
     .                                                     1,lb)

                bedge_facey_x(:,1+nguard0:nxb+nguard0,2,klo,lb)=
     .          bedge_facez_x(:,1+nguard0:nxb+nguard0,1+nyb+nguard0
     .                                                     ,1,lb)

                bedge_facex_y(:,1,1+nguard0:nyb+nguard0,klo,lb)=
     .          bedge_facez_y(:,1+nguard0,1+nguard0:nyb+nguard0,
     .                                                     1,lb)

                bedge_facex_y(:,2,1+nguard0:nyb+nguard0,klo,lb)=
     .          bedge_facez_y(:,1+nxb+nguard0,1+nguard0:nyb+nguard0
     .                                                    ,1,lb)


              elseif(jf.eq.6) then 

                call shmem_real_get(recvarz1e(1,1,1,1),
     .               bedge_facez_x(1,1,1,1,remote_block),
     .               len_block_ez*nedges,remote_pe)
                bedge_facez_x(:,:,:,2,lb) = recvarz1e(:,:,:,1)
                call shmem_real_get(recvarz1e(1,1,1,1),
     .               bedge_facez_y(1,1,1,1,remote_block),
     .               len_block_ez*nedges,remote_pe)
                bedge_facez_y(:,:,:,2,lb) = recvarz1e(:,:,:,1)

! make common variables on an edge consistent
                bedge_facey_x(:,1+nguard0:nxb+nguard0,1,kup1,lb)=
     .          bedge_facez_x(:,1+nguard0:nxb+nguard0,1+nguard0,2,lb)

                bedge_facey_x(:,1+nguard0:nxb+nguard0,2,kup1,lb)=
     .          bedge_facez_x(:,1+nguard0:nxb+nguard0,1+nyb+nguard0,
     .                                                        2,lb)

                bedge_facex_y(:,1,1+nguard0:nyb+nguard0,kup1,lb)=
     .          bedge_facez_y(:,1+nguard0,1+nguard0:nyb+nguard0,2,lb)

                bedge_facex_y(:,2,1+nguard0:nyb+nguard0,kup1,lb)=
     .          bedge_facez_y(:,1+nxb+nguard0,1+nguard0:nyb+nguard0,
     .                                                          2,lb)



              endif


          endif

        enddo

      endif
      enddo
      endif

      call shmem_barrier_all()

!!!#endif
      end if

      return
      end subroutine amr_edge_average_udt
