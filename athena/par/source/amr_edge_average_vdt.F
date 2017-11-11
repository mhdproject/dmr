!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"


      subroutine amr_edge_average_vdt(mype,nsub)




!------------------------------------------------------------------------
!
! This routine gets cell edge-based data at block boundaries from 
! neighbors who are parents of leaf blocks. 
!
! The data structure used to store and pass this data is defined
! in the include file 'block_boundary_data.fh' which can be included
! in 'physicaldata.fh'.
!
! This version is used when variable timesteps are allowed across the
! blocks in the computation.
!
! Arguments:
!       mype            processor number
!       nsub            current time subcycle. If this is 1 then this
!                       info is used to reset the temporary boundary edge
!                       data arrays to 0.
!
!
! Written :     Peter MacNeice          August 1997
!------------------------------------------------------------------------


      use paramesh_dimensions
      use physicaldata
      use tree

      use paramesh_interfaces, only : amr_restrict_edge_data_vdt


#include "amr_shmem.fh"

      integer, intent(in)  ::  mype,nsub

!------------------------------------
! local variables

        integer remote_pe,remote_block
        integer cnodetype,cchild(2),ich
        logical lnodetime
        save    lnodetime,cchild,cnodetype

      integer :: lb,jf

      integer nguard0
      parameter(nguard0 = nguard*npgs)

      integer ng_off
      parameter(ng_off = nguard0+iface_off)

      integer kup,klo
      parameter(klo  = 1+k3d*nguard0)
      parameter(kup  = 1+k3d*(nzb+nguard0-1))
      parameter(kup1 = k3d+nzb+k3d*nguard0)

!------------------------------------

#ifdef VAR_DT


        if(lnblocks.gt.0) then
        do lb = 1,lnblocks

! Is this a parent of at least one leaf block ?
      if(nodetype(lb).eq.2) then


! Set timestep phases for the current block, and for the next finer level.
        lcycle = loc_cycle(lrefine(lb))
        phase0 = phase_dt(lrefine(lb))
        phase1 = phase_dt(lrefine(lb)+1)

! At start of the current blocks timestep zero out the arrays used to
! accumulate boundary edge data from its children.
        if(lcycle.eq.1) then
           ttbedge_facex_y(:,:,:,:,lb) = 0.
           ttbedge_facey_x(:,:,:,:,lb) = 0.
           if(ndim.eq.3) then
             ttbedge_facex_z(:,:,:,:,lb) = 0.
             ttbedge_facey_z(:,:,:,:,lb) = 0.
             ttbedge_facez_x(:,:,:,:,lb) = 0.
             ttbedge_facez_y(:,:,:,:,lb) = 0.
           endif
        endif

      endif
      enddo
      endif
      call shmem_barrier_all()
!------------------------------------

! Leaf blocks which have completed their timestep provide reduced
! boundary edge data to their parents.
! Edge values are accumulated in the ttbedge_face arrays.
      call amr_restrict_edge_data_vdt(mype)


      call shmem_barrier_all()
!------------------------------------

! Parents who have completed their timestep and border a leaf block
! update their edges.
      do lb = 1,lnblocks


! Is this a parent block of at least one leaf node?
      if((nodetype(lb).eq.2).and.ldtcomplete(lb)) then

! If yes then cycle through its neighbors.
        do iface=1,nfaces

! If this neighbor is a leaf block or an external boundary then
! replace fluxes with restricted fluxes.
          cnodetype = 1
          if(neigh(1,iface,lb).ge.1) then
            remote_pe    = neigh(2,iface,lb)
            remote_block = neigh(1,iface,lb)

            call shmem_integer_get(cnodetype,nodetype(remote_block),
     .                             1,remote_pe)
          endif
          if(cnodetype.eq.1) then

            if(iface.eq.1) then

              bedge_facex_y(:,1,:,:,lb)=ttbedge_facex_y(:,1,:,:,lb)
              if((ndim.eq.3).or.(l2p5d.eq.1)) 
     .        bedge_facex_z(:,1,:,:,lb)=ttbedge_facex_z(:,1,:,:,lb)

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

            elseif(iface.eq.2) then

              bedge_facex_y(:,2,:,:,lb)=ttbedge_facex_y(:,2,:,:,lb)
              if((ndim.eq.3).or.(l2p5d.eq.1)) 
     .        bedge_facex_z(:,2,:,:,lb)=ttbedge_facex_z(:,2,:,:,lb)
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

            elseif(iface.eq.3) then

              bedge_facey_x(:,:,1,:,lb)=ttbedge_facey_x(:,:,1,:,lb)
              if((ndim.eq.3).or.(l2p5d.eq.1)) 
     .        bedge_facey_z(:,:,1,:,lb)=ttbedge_facey_z(:,:,1,:,lb)
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

            elseif(iface.eq.4) then
              bedge_facey_x(:,:,2,:,lb)=ttbedge_facey_x(:,:,2,:,lb)
              if((ndim.eq.3).or.(l2p5d.eq.1)) 
     .        bedge_facey_z(:,:,2,:,lb)=ttbedge_facey_z(:,:,2,:,lb)
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

            elseif(iface.eq.5) then
              bedge_facez_x(:,:,:,1,lb)=ttbedge_facez_x(:,:,:,1,lb)
              bedge_facez_y(:,:,:,1,lb)=ttbedge_facez_y(:,:,:,1,lb)
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

            elseif(iface.eq.6) then
              bedge_facez_x(:,:,:,2,lb)=ttbedge_facez_x(:,:,:,2,lb)
              bedge_facez_y(:,:,:,2,lb)=ttbedge_facez_y(:,:,:,2,lb)
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


      call shmem_barrier_all()
!------------------------------------

! cycle through the grid blocks on this processor
      if(lnblocks.gt.0) then
      do lb = 1,lnblocks

! Is this a leaf block and not at the original refinement level ?
!      if(nodetype(lb).eq.1.and.lrefine(lb).gt.1) then
      if(nodetype(lb).eq.1) then

! Has this block completed its timestep?
      if(ldtcomplete(lb)) then

! Cycle over the blocks faces
       do jf = 1,nfaces

          remote_pe = neigh(2,jf,lb)
          remote_block  = neigh(1,jf,lb)

! Is the neighbor to this face a parent of a leaf block?
          cnodetype = 0
          if(remote_block.gt.0) then
             call shmem_integer_get(cnodetype,
     .                     nodetype(remote_block),1,remote_pe)
          endif


       if(cnodetype.eq.2) then


          if(jf.eq.1) then

            call shmem_real_get(recvarx1e(1,1,1,1),
     .                          bedge_facex_y(1,1,1,1,remote_block),
     .                          len_block_ex*nedges,remote_pe)
            bedge_facex_y(:,1,:,:,lb) = recvarx1e(:,2,:,:)

            if(ndim.eq.3) then
              call shmem_real_get(recvarx1e(1,1,1,1),
     .                            bedge_facex_z(1,1,1,1,remote_block),
     .                            len_block_ex*nedges,remote_pe)
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
     .                          bedge_facex_y(1,1,1,1,remote_block),
     .                          len_block_ex*nedges,remote_pe)
            bedge_facex_y(:,2,:,:,lb) = recvarx1e(:,1,:,:)

            if(ndim.eq.3) then
              call shmem_real_get(recvarx1e(1,1,1,1),
     .                            bedge_facex_z(1,1,1,1,remote_block),
     .                            len_block_ex*nedges,remote_pe)
              bedge_facex_z(:,2,:,:,lb) = recvarx1e(:,1,:,:)
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

            if(ndim.eq.3) then
              call shmem_real_get(recvary1e(1,1,1,1),
     .                            bedge_facey_z(1,1,1,1,remote_block),
     .                            len_block_ey*nedges,remote_pe)
              bedge_facey_z(:,:,1,:,lb) = recvary1e(:,:,2,:)
            endif
            call shmem_real_get(recvary1e(1,1,1,1),
     .                          bedge_facey_x(1,1,1,1,remote_block),
     .                          len_block_ey*nedges,remote_pe)
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

!            call shmem_real_get(recvary1e(1,1,1,1),
!     .                          bedge_facey_x(1,1,1,1,remote_block),
!     .                          len_block_ey*nedges,remote_pe)
!            bedge_facey_x(:,:,1,:,lb) = recvary1e(:,:,2,:)


          elseif(jf.eq.4) then

            if(ndim.eq.3) then
              call shmem_real_get(recvary1e(1,1,1,1),
     .                            bedge_facey_z(1,1,1,1,remote_block),
     .                            len_block_ey*nedges,remote_pe)
              bedge_facey_z(:,:,2,:,lb) = recvary1e(:,:,1,:)
            endif

            call shmem_real_get(recvary1e(1,1,1,1),
     .                          bedge_facey_x(1,1,1,1,remote_block),
     .                          len_block_ey*nedges,remote_pe)
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
     .                          bedge_facez_x(1,1,1,1,remote_block),
     .                          len_block_ez*nedges,remote_pe)
            bedge_facez_x(:,:,:,1,lb) = recvarz1e(:,:,:,2)


            call shmem_real_get(recvarz1e(1,1,1,1),
     .                          bedge_facez_y(1,1,1,1,remote_block),
     .                          len_block_ez*nedges,remote_pe)
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
     .                          bedge_facez_x(1,1,1,1,remote_block),
     .                          len_block_ez*nedges,remote_pe)
            bedge_facez_x(:,:,:,2,lb) = recvarz1e(:,:,:,1)


            call shmem_real_get(recvarz1e(1,1,1,1),
     .                          bedge_facez_y(1,1,1,1,remote_block),
     .                          len_block_ez*nedges,remote_pe)
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

      endif                      ! end of ldtcomplete if test

      endif
      enddo
      endif

      call shmem_barrier_all()
!------------------------------------

#endif /* VAR_DT */

      return
      end subroutine amr_edge_average_vdt
