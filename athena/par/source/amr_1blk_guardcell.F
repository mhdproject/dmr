!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"

!#define DEBUG

      subroutine amr_1blk_guardcell(mype,iopt,nlayers,lb,pe,
     .                              lcc,lfc,lec,lnc,
     .                              l_srl_only,icoord,ldiag,
     .                              nlayersx,nlayersy,nlayersz)


      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace

      use paramesh_interfaces, only : amr_gsurrounding_blks,
     &                                amr_perm_to_1blk,
     &                                amr_1blk_guardcell_srl,
     &                                amr_1blk_guardcell_c_to_f,
     &                                amr_1blk_guardcell_f_to_c

!------------------------------------------------------------------------

!
! This routine manages the transfer of guard cell data for a
! specific single block.
!
!
! IMPORTANT NOTE :
! This routine must be used with care!
! This routine was written to be used in a code as illustrated
! in the following snippet of pseudo-code
!
!              .
!              .
!              .
!        synchronization pt
!
!        loop over grid blocks
!          (copy current block into working block and fill its guardcells)
!          (perform some set of operations on block)
!          (store result from working block)
!        end loop
!
!        synchronization pt
!              .
!              .
!              .
!
! Caveat 1:
! If you are using this routine, you must remember to save the solution
! at the first synchronization point (ie call amr_1blk_copy_soln), so 
! that each block uses the same time synchronized solution during its 
! guardcell filling.
!
! Caveat 2:
! It is implicitly assumed that the parent blocks on all leaf nodes
! have valid data. (This is necessary to ensure that a general restriction
! operator can be supported in the neighborhood of a jump in refinement.)
! If ADVANCE_ALL_LEVELS is defined, then this will generally be true. 
! However, if the solution is being time-advanced on leaf blocks only 
! this may not be true. In this case you should call amr_restrict  
! before the first synchronization pt in the example above.
! If you are using blocks with an even number of grid points and the
! default restriction operators this is not necessary.
!
!
! Written :     Peter MacNeice          July 1998
!------------------------------------------------------------------------
!
! Arguments:
!      mype           local processor number
!      iopt           a switch to control which data source is to be used
!                      iopt=1 will use 'unk'
!                      iopt>=2 will use 'work'
!      nlayers        the number of guard cell layers at each boundary
!      lb             the selected block
!      pe             processor storing the selected block
!      lcc            a logical switch controlling whether unk or work data
!                      is filled
!      lfc            a logical switch controlling whether facevar data
!                      is filled
!      lec            a logical switch controlling whether unk_e_x(y)(z) data
!                      is filled
!      lnc            a logical switch controlling whether unk_n data
!                      is filled
!      l_srl_only     a logical switch which, if true, switches off the
!                      filling from coarse neighbors. This is used during
!                      restriction when odd  block sizes are in use.
!      icoord         an integer switch used to select which faces of
!                      the block are to be considered. If icoord=0 all
!                      faces are considered. If icoord=1 only faces perp.
!                      to the x-axis are considered, if icoord=2 only faces
!                      perp. to the y-axis are considered, and if icoord=3
!                      only faces perp. to the z-axis are considered.
!      ldiag          a logical switch which controls whether guardcells
!                      corresponding to neighbor blocks diagonally opposite
!                      block edges and corners are filled.
!      nlayersx       Optional arguments !!!!! 
!      nlayersy        Specify the number of guardcells to fill in each
!      nlayersz        direction.
!
!------------------------------------

#include "amr_shmem.fh"

      integer, intent(in) :: mype,iopt,nlayers,lb,pe,icoord
      logical, intent(in) :: lcc,lfc,lec,lnc,l_srl_only,ldiag
      integer, intent(in), optional :: nlayersx,nlayersy,nlayersz

!------------------------------------


      integer,save :: parent_lb,parent_pe
      integer remote_blk,remote_pe,cneigh(2,mfaces)
      save    remote_blk,remote_pe,cneigh

      logical :: lcoarse
      logical :: loc_lcc,loc_lfc,loc_lec,loc_lnc


! local arrays
      integer :: surrblks(3,3,3,3)
      integer :: psurrblks(3,3,3,3)
      integer :: csurr_blks(3,3,1+2*k2d,1+2*k3d)
      save csurr_blks
      integer pcache_pe,pcache_blk

      integer :: icoord_loc
      logical :: ldiag_loc

      real recv(nvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd)
      real recvx(nbndvar,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,
     .       kl_bnd:ku_bnd)
      real recvy(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,
     .       kl_bnd:ku_bnd)
      real recvz(nbndvar,il_bnd:iu_bnd,jl_bnd:ju_bnd,
     .       kl_bnd:ku_bnd+k3d)
      real recvex(nbndvare,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,
     .                    kl_bnd:ku_bnd+k3d)
      real recvey(nbndvare,il_bnd:iu_bnd+1,jl_bnd:ju_bnd,
     .                    kl_bnd:ku_bnd+k3d)
      real recvez(nbndvare,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+k2d,
     .                    kl_bnd:ku_bnd)
      real recvn(nbndvarc,il_bnd:iu_bnd+1,jl_bnd:ju_bnd+k2d,
     .                    kl_bnd:ku_bnd+k3d)
      save recv
      save recvx,recvy,recvz,recvex,recvey,recvez,recvn

      integer :: nsurr_blks,i,j,k,idest,iopt0,iblock

      integer nguard0,nguard_work0
      parameter (nguard0 = nguard*npgs )
      parameter (nguard_work0 = nguard_work*npgs )


!------------------------------------
!
! The sequence of operations required to fill the guardcells of an
! individual block is :
!
! Step 1:
! Construct a list of blocks surrounding block lb
!
! Step 2:
! Check for coarse neighbors
!
! Step 3:
! Put leaf block data into the data_1blk.fh datastructures, with the
! appropriate guardcell padding.
!
! Step 4:
! Put data from leaf blocks parent into the data_1blk.fh datastructures, 
! with the appropriate guardcell padding. Check to see if this data is 
! currently cached.
!
! Step 5:
! Construct a list of blocks surrounding block lb's parent
!
! Step 6:
! Do guardcell filling for lb's parent from any surrounding blocks at 
! the same refinement level as this parent.
!
! Step 7:
! Do guardcell filling from coarse neigbors into the current block
!
! Step 8:
! Do guardcell filling from any surrounding blocks at the same refinement
! level as the leaf block.
!
! Step 9:
! Apply boundary conditions.
!
!
!------------------------------------

        l_f_to_c = .false.

      if (iopt == 1) then
         if (.not.present(nlayersx)) then
            nlayers0x = nguard
         else
            nlayers0x = nlayersx
         end if
         if (.not.present(nlayersy)) then
            nlayers0y = nguard
         else
            nlayers0y = nlayersy
         end if
         if (.not.present(nlayersz)) then
            nlayers0z = nguard
         else
            nlayers0z = nlayersz
         end if
      else
         if (.not.present(nlayersx)) then
            nlayers0x = nguard_work
         else
            nlayers0x = nlayersx
         end if
         if (.not.present(nlayersy)) then
            nlayers0y = nguard_work
         else
            nlayers0y = nlayersy
         end if
         if (.not.present(nlayersz)) then
            nlayers0z = nguard_work
         else
            nlayers0z = nlayersz
         end if
      end if


#ifdef AMR_ERROR_CHECKING
        loc_lcc = .false.
        loc_lfc = .false.
        loc_lnc = .false.
        loc_lec = .false.
        if(nvar.gt.0)     loc_lcc = .true.
        if(nfacevar.gt.0) loc_lfc = .true.
        if(nvaredge.gt.0) loc_lec = .true.
        if(nvarcorn.gt.0) loc_lnc = .true.
        if( (lcc.and.(.not.loc_lcc)) .or. (lfc.and.(.not.loc_lfc)) 
     .  .or.(lec.and.(.not.loc_lec)) .or. (lnc.and.(.not.loc_lnc)) 
     .    ) then
          if(mype.eq.0) then
            write(*,*)
     .         'Paramesh Error: for a call to amr_1blk_guardcell',
     .         ' one of more of the arguments lcc/lfc/lec/lnc are not',
     .         ' consistent with nvar/nfacevar/nvaredge/nvarcorn.'
          endif
          call shmem_barrier_all()
          call amr_abort
        endif
#endif /* AMR_ERROR_CHECKING */


!
! Construct a list of blocks surrounding block lb
        if(gsurrblks_set.ne.1) then
          if(mype.eq.0) then
            write(*,*) 'Paramesh Error: You are probably missing ',
     .                 'a call to amr_gsurrounding_blks.'
          endif
          call shmem_barrier_all()
          call amr_abort
        endif

        nsurr_blks = 3*3*(1+2*k2d)*(1+2*k3d)
        call shmem_integer_get(csurr_blks,surr_blks(1,1,1,1,lb),
     .                         nsurr_blks,pe)
        surrblks(:,:,2-k2d:2+k2d,2-k3d:2+k3d) = 
     .                csurr_blks(:,:,1:1+2*k2d,1:1+2*k3d)


        if(.not.l_srl_only) then

!
! Are there any coarse neighbors?
        lcoarse = .false.

        if(parent(1,lb).gt.0) then
        do k = 2-k3d,2+k3d
        do j = 2-k2d,2+k2d
        do i = 1,3
        if(surrblks(1,i,j,k).gt.-20.and.surrblks(1,i,j,k).lt.0)
     .                   lcoarse=.true.
        enddo
        enddo
        enddo
        endif 

        endif                        ! end of l_srl_only if test

!------------------------------------
!
! Put leaf block lb's data into the data_1blk.fh datastructures, with the
! appropriate guardcell padding.
          idest = 1
          call amr_perm_to_1blk(lcc,lfc,lec,lnc,lb,pe,iopt,idest)


!------------------------------------
!
        if(iopt.eq.1) then
          pcache_pe = pcache_pe_u
          pcache_blk = pcache_blk_u
        elseif(iopt.ge.2) then
          pcache_pe = pcache_pe_w
          pcache_blk = pcache_blk_w
        endif


        if(.not.l_srl_only) then
        if(lcoarse) then
!
!
! Put data from lb's parent into the data_1blk.fh datastructures, with the
! appropriate guardcell padding. Check to see if data is currently cached.
        if (mype == pe) then
           parent_lb = parent(1,lb)
           parent_pe = parent(2,lb)
        else
           call shmem_integer_get(parent_lb,parent(1,lb),1,pe)
           call shmem_integer_get(parent_pe,parent(2,lb),1,pe)
        end if

#ifdef DEBUG
        write(*,*) 'pe ',mype,
     .      'parent blk ',parent_lb,parent_pe,' cached ',
     .      pcache_blk,pcache_pe
#endif
        if( (parent_lb.gt.0) .and.
     .      ((parent_lb.ne.pcache_blk).or.(parent_pe.ne.pcache_pe) )
     .      ) then

! record id of new parent block placed in cache
          lnew_parent = .true.
          pcache_blk = parent_lb
          pcache_pe = parent_pe

#ifdef DEBUG
        write(*,*) 'pe ',mype,
     .      'fetching parent blk ',parent_lb,parent_pe,
     .      'new cache ',pcache_blk,pcache_pe,' lnew_parent ',
     .      lnew_parent
#endif

          if(iopt.eq.1) then
          if(lcc) then
            unk1(:,:,:,:,2) = 0.
#ifdef NO_PERMANENT_GUARDCELLS
            call shmem_real_get(recv(1,1,1,1),
     .                          gt_unk(1,1,1,1,parent_lb),
     .                          len_block,parent_pe)
#else
            call shmem_real_get(recv(1,1,1,1),
     .                          unk(1,1,1,1,parent_lb),
     .                          len_block,parent_pe)
#endif
            unk1(:,1+nguard:nxb+nguard,1+nguard*k2d:nyb+nguard*k2d,
     .            1+nguard*k3d:nzb+nguard*k3d,2) =
     .  recv(:,1+nguard0:nxb+nguard0,1+nguard0*k2d:nyb+nguard0*k2d,
     .            1+nguard0*k3d:nzb+nguard0*k3d)
          endif

          if(lfc) then
            facevarx1(:,:,:,:,2) = 0.
#ifdef NO_PERMANENT_GUARDCELLS
            call shmem_real_get(recvx(1,1,1,1),
     .                          gt_facevarx(1,1,1,1,parent_lb),
     .                          len_blockfx*nfacevar,parent_pe)
#else
            call shmem_real_get(recvx(1,1,1,1),
     .                          facevarx(1,1,1,1,parent_lb),
     .                          len_blockfx*nfacevar,parent_pe)
#endif
            facevarx1(:,1+nguard:nxb+nguard+1,
     .                  1+nguard*k2d:nyb+nguard*k2d,
     .                  1+nguard*k3d:nzb+nguard*k3d,2) =
     .        recvx(:,1+nguard0:nxb+nguard0+1,
     .                  1+nguard0*k2d:nyb+nguard0*k2d,
     .                  1+nguard0*k3d:nzb+nguard0*k3d)


            if(ndim.ge.2) then
              facevary1(:,:,:,:,2) = 0.
#ifdef NO_PERMANENT_GUARDCELLS
              call shmem_real_get(recvy(1,1,1,1),
     .                          gt_facevary(1,1,1,1,parent_lb),
     .                          len_blockfx*nfacevar,parent_pe)
#else
              call shmem_real_get(recvy(1,1,1,1),
     .                          facevary(1,1,1,1,parent_lb),
     .                          len_blockfx*nfacevar,parent_pe)
#endif
              facevary1(:,1+nguard:nxb+nguard,
     .                   1+nguard*k2d:nyb+(nguard+1)*k2d,
     .                   1+nguard*k3d:nzb+nguard*k3d,2) =
     .         recvy(:,1+nguard0:nxb+nguard0,
     .                   1+nguard0*k2d:nyb+(nguard0+1)*k2d,
     .                   1+nguard0*k3d:nzb+nguard0*k3d)
            endif
            if(ndim.eq.3) then
              facevarz1(:,:,:,:,2) = 0.
#ifdef NO_PERMANENT_GUARDCELLS
              call shmem_real_get(recvz(1,1,1,1),
     .                          gt_facevarz(1,1,1,1,parent_lb),
     .                          len_blockfx*nfacevar,parent_pe)
#else
              call shmem_real_get(recvz(1,1,1,1),
     .                          facevarz(1,1,1,1,parent_lb),
     .                          len_blockfx*nfacevar,parent_pe)
#endif
              facevarz1(:,1+nguard:nxb+nguard,
     .                   1+nguard*k2d:nyb+nguard*k2d,
     .                   1+nguard*k3d:nzb+(nguard+1)*k3d,2) =
     .         recvz(:,1+nguard0:nxb+nguard0,
     .                   1+nguard0*k2d:nyb+nguard0*k2d,
     .                   1+nguard0*k3d:nzb+(nguard0+1)*k3d)
            endif
          endif


!!!#if N_DIM > 1
          if (ndim > 1) then
          if(lec) then
            unk_e_x1(:,:,:,:,2) = 0.
#ifdef NO_PERMANENT_GUARDCELLS
            call shmem_real_get(recvex(1,1,1,1),
     .                          gt_unk_e_x(1,1,1,1,parent_lb),
     .                          len_blockex*nvaredge,parent_pe)
#else
            call shmem_real_get(recvex(1,1,1,1),
     .                          unk_e_x(1,1,1,1,parent_lb),
     .                          len_blockex*nvaredge,parent_pe)
#endif
            unk_e_x1(:,1+nguard:nxb+nguard,
     .                  1+nguard*k2d:nyb+(nguard+1)*k2d,
     .                  1+nguard*k3d:nzb+(nguard+1)*k3d,2) =
     .        recvex(:,1+nguard0:nxb+nguard0,
     .                  1+nguard0*k2d:nyb+(nguard0+1)*k2d,
     .                  1+nguard0*k3d:nzb+(nguard0+1)*k3d)

            unk_e_y1(:,:,:,:,2) = 0.


#ifdef NO_PERMANENT_GUARDCELLS
            call shmem_real_get(recvey(1,1,1,1),
     .                          gt_unk_e_y(1,1,1,1,parent_lb),
     .                          len_blockey*nvaredge,parent_pe)
#else
            call shmem_real_get(recvey(1,1,1,1),
     .                          unk_e_y(1,1,1,1,parent_lb),
     .                          len_blockey*nvaredge,parent_pe)
#endif
            unk_e_y1(:,1+nguard:nxb+nguard+1,
     .                  1+nguard*k2d:nyb+nguard*k2d,
     .                  1+nguard*k3d:nzb+(nguard+1)*k3d,2) =
     .        recvey(:,1+nguard0:nxb+nguard0+1,
     .                  1+nguard0*k2d:nyb+nguard0*k2d,
     .                  1+nguard0*k3d:nzb+(nguard0+1)*k3d)


!!!#if N_DIM == 3
            if (ndim == 3) then
            unk_e_z1(:,:,:,:,2) = 0.
#ifdef NO_PERMANENT_GUARDCELLS
            call shmem_real_get(recvez(1,1,1,1),
     .                          gt_unk_e_z(1,1,1,1,parent_lb),
     .                          len_blockez*nvaredge,parent_pe)
#else
            call shmem_real_get(recvez(1,1,1,1),
     .                          unk_e_z(1,1,1,1,parent_lb),
     .                          len_blockez*nvaredge,parent_pe)
#endif
            unk_e_z1(:,1+nguard:nxb+nguard+1,
     .                  1+nguard*k2d:nyb+(nguard+1)*k2d,
     .                  1+nguard*k3d:nzb+nguard*k3d,2) =
     .        recvez(:,1+nguard0:nxb+nguard0+1,
     .                  1+nguard0*k2d:nyb+(nguard0+1)*k2d,
     .                  1+nguard0*k3d:nzb+nguard0*k3d)
!!!#endif /* N_DIM == 3 */
            end if

          endif
!!!#endif /* N_DIM > 1 */
          end if


          if(lnc) then
            unk_n1(:,:,:,:,2) = 0.
#ifdef NO_PERMANENT_GUARDCELLS
            call shmem_real_get(recvn(1,1,1,1),
     .                          gt_unk_n(1,1,1,1,parent_lb),
     .                          len_blockn*nvarcorn,parent_pe)
#else
            call shmem_real_get(recvn(1,1,1,1),
     .                          unk_n(1,1,1,1,parent_lb),
     .                          len_blockn*nvarcorn,parent_pe)
#endif
            unk_n1(:,1+nguard:nxb+nguard+1,
     .                  1+nguard*k2d:nyb+(nguard+1)*k2d,
     .                  1+nguard*k3d:nzb+(nguard+1)*k3d,2) =
     .        recvn(:,1+nguard0:nxb+nguard0+1,
     .                  1+nguard0*k2d:nyb+(nguard0+1)*k2d,
     .                  1+nguard0*k3d:nzb+(nguard0+1)*k3d)

          endif


          elseif(iopt.ge.2) then
            iopt0 = iopt-1
            work1(:,:,:,2) = 0.
            call shmem_real_get(recvw(1,1,1),
     .                          work(1,1,1,parent_lb,iopt0),
     .                          len_wblock,parent_pe)
            work1(1+nguard_work:nxb+nguard_work,
     .            1+nguard_work*k2d:nyb+nguard_work*k2d,
     .            1+nguard_work*k3d:nzb+nguard_work*k3d,2) =
     .         recvw(1+nguard_work0:nxb+nguard_work0,
     .               1+nguard_work0*k2d:nyb+nguard_work0*k2d,
     .               1+nguard_work0*k3d:nzb+nguard_work0*k3d)

          endif                 ! end of iopt if test

!------------------------------------
!
! Construct a list of blocks surrounding block lb's parent.

        nsurr_blks = 3*3*(1+2*k2d)*(1+2*k3d)
        call shmem_integer_get(csurr_blks,
     .                         surr_blks(1,1,1,1,parent_lb),
     .                         nsurr_blks,parent_pe)
        psurrblks(:,:,2-k2d:2+k2d,2-k3d:2+k3d) = 
     .                csurr_blks(:,:,1:1+2*k2d,1:1+2*k3d)

!------------------------------------
! Do guardcell filling for lb's parent from any surrounding blocks at 
! the same refinement level as this parent.
! Diagonal elements are required to ensure that all cells are filled
! correctly when icoord is non-zero.
          iblock=2
          icoord_loc = 0
          ldiag_loc = .true.
          call amr_1blk_guardcell_srl(mype,parent_pe,parent_lb,
     .                                iblock,iopt,nlayers,psurrblks,
     .                                lcc,lfc,lec,lnc,
     .                                icoord_loc,ldiag_loc,
     .                                nlayers0x,nlayers0y,nlayers0z)


        endif       ! end if parents data not previously cached


!------------------------------------
!
! Do guardcell filling from coarse neigbors into the current block
        call amr_1blk_guardcell_c_to_f( mype,lb,pe,iopt,nlayers,
     .                                  surrblks,lcc,lfc,lec,lnc,
     .                                  icoord,ldiag,
     .                             nlayers0x,nlayers0y,nlayers0z)


!------------------------------------

        endif                       ! end of lcoarse if test

        endif                       ! end of l_srl_only if test


!
!------------------------------------
!
! Do guardcell filling from any surrounding blocks at the same refinement
! level as block lb.
        iblock = 1
        call amr_1blk_guardcell_srl(mype,pe,lb,
     .                              iblock,iopt,nlayers,surrblks,
     .                              lcc,lfc,lec,lnc,icoord,ldiag,
     .                              nlayers0x,nlayers0y,nlayers0z)

!        call amr_1blk_guardcell_f_to_c(mype,pe,lb,
!     .                              iblock,iopt,nlayers,surrblks,
!     .                              lcc,lfc,lec,lnc,icoord,ldiag,
!     .                              nlayers0x,nlayers0y,nlayers0z)

        if(iopt.eq.1) then
          pcache_pe_u = pcache_pe
          pcache_blk_u = pcache_blk
        elseif(iopt.ge.2) then
          pcache_pe_w = pcache_pe
          pcache_blk_w = pcache_blk
        endif


      return
      end subroutine amr_1blk_guardcell
