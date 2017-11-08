!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------
!
! Modification history:
!     Michael L. Rilee, November 2002, *dbz*
!        Initial support for divergenceless prolongation
!     Michael L. Rilee, December 2002, *clean_divb*
!        Support for projecting field onto divergenceless field
!
!----------------------------------------------------------------------
#include "paramesh_preprocessor.fh"
!#define MULTIBLOCK_TEST_2D

#undef DEBUG
!#define DEBUG
#undef MLRDBG
!#define MLRDBG
#undef MLRDBG1
!#define MLRDBG1
      subroutine amr_prolong(mype,iopt,nlayers)


!------------------------------------------------------------------------
!
! This routine interpolates data from parent blocks to any
! newly created child blocks. 
!
! Note, this is controlled by the NEWCHILD array, which is 
! reset to .false. everywhere on exit. Thus if you wish to
! make another call to AMC_PROLONG using the same refinement
! change defined by NEWCHILD you must save a copy of NEWCHILD
! in the calling program.
! 
!
! This routine calls  routines prolong_cc and prolong_fc which
! actually perform the prolongation on the cell center and cell face
! center datastructures respectively.
!
! Arguments :
!      mype     integer        Current processor number
!      iopt     integer        Switch to select which datastructures
!                               are updated. If iopt=1 the prolong_cc
!                               acts on UNK, and if NBNDVAR is > 0 
!                               prolong_fc updates the FACEVARX(Y)(Z)
!                               arrays. If iopt=2 only WORK is updated.
!      nlayers  integer        Number of layers of guard cells at a
!                               block boundary.
!
!
! Written :     Peter MacNeice          July 1997
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace

      use prolong_arrays, only : 
     &     prol_fc_dbz,
     &     prol_fc_dbz_ivar,
     &     prol_fc_dbz_n, 
     &     prol_fc_clean_divb

      use paramesh_interfaces, only : 
     &                  amr_1blk_cc_prol_gen_work_fun,
     &                  amr_1blk_fc_prol_gen_fun,
     &                  amr_1blk_fc_prol_dbz,
     &                  amr_1blk_ec_prol_gen_fun,
     &                  amr_1blk_nc_prol_gen_fun,
     &                  amr_1blk_copy_soln,
     &                  amr_1blk_guardcell_reset,
     &                  amr_1blk_guardcell,
     &                  amr_1blk_cc_prol_gen_unk_fun,
     &                  comm_int_max_to_all,
     &                  comm_int_min_to_all



#include "amr_shmem.fh"

      integer, intent(in) ::  mype,iopt,nlayers

!------------------------------------
! local arrays
      real :: recv1(nvar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,
     .                                      kl_bnd1:ku_bnd1)
      real :: recvf(nbndvar,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,
     .                                      kl_bnd1:ku_bnd1+k3d)

!mlr begin
      real :: recvfx(nbndvar,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,
     .                                      kl_bnd1:ku_bnd1)
      real :: recvfy(nbndvar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,
     .                                      kl_bnd1:ku_bnd1)
      real :: recvfz(nbndvar,il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,
     .                                      kl_bnd1:ku_bnd1+k3d)
!mlr end

      real :: recve(nbndvare,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,
     .                                      kl_bnd1:ku_bnd1+k3d)
      real :: recvn(nbndvarc,il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,
     .                                      kl_bnd1:ku_bnd1+k3d)

      integer p_cache_addr(2),idest,nguard0
      integer parent_blk,parent_pe
      integer parent_list(2,maxblocks)
      integer indx(maxblocks),index(maxblocks)
      integer parent_id(maxblocks),icoord

      logical :: lcc,lfc,lec,lnc,l_srl_only,ldiag

      integer,save :: lref_min,lref_max
      integer,save :: lref_mint,lref_maxt
      integer nnewchildg,nnewchild
      logical lnewchild
      save nnewchildg,nnewchild

      integer :: ia,ib,ja,jb,ka,kb,iblock,lb,lreflevel,lbi
      integer :: ioff,joff,koff

      integer :: iib,ilb,nfield

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer :: imlr, jmlr, kmlr,idim_mlr
      integer :: iprol, iv1, iv2, iv3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef TIMINGS
#include "timer.fh"
#endif

!------------------------------------

#ifdef TIMINGS
      call shmem_barrier_all()
      itimer1 = irtc()
#endif

!------------------------------------
!
! Are there any new children?
       nnewchild = 0
       lnewchild = any(newchild)
       if(lnewchild) nnewchild = 1
       call shmem_barrier_all()
       call comm_int_max_to_all (nnewchildg,nnewchild)
       call shmem_barrier_all()
#ifdef DEBUG
       write(*,*) 'prolong: nnewchildg ',nnewchildg
#endif /* DEBUG */
       if(nnewchildg.eq.0) return


! reset cache addresses
       call amr_1blk_guardcell_reset


!#ifdef NO_PERMANENT_GUARDCELLS
! Store a copy of the current solution in gt_unk
!        call amr_1blk_copy_soln
!#endif



! Identify variables to be prolonged
       lcc = .false.
       lfc = .false.
       lec = .false.
       lnc = .false.
       if(iopt.eq.1) then
         if(nvar.gt.0)     lcc = .true.
         if(nfacevar.gt.0) lfc = .true.
         if(nvaredge.gt.0) lec = .true.
         if(nvarcorn.gt.0) lnc = .true.
         nguard0 = nguard
       elseif(iopt.ge.2) then
         lcc = .true.
         nguard0 = nguard_work
       endif

       ia = 1+nguard0       
       ib = nxb+nguard0       
       ja = 1+nguard0*k2d      
       jb = nyb+nguard0*k2d
       ka = 1+nguard0*k3d      
       kb = nzb+nguard0*k3d

! construct a list of parent blocks for new children
       parent_list(:,:) = -1
       iblock=0
       lref_min = 100
       lref_max = 1
       if(lnblocks.gt.0) then
       do lb = 1,lnblocks
       if(newchild(lb)) then
         iblock = iblock+1
         index(iblock) = lb
         indx(iblock) = iblock
         parent_list(:,iblock) = parent(:,lb)
         parent_id(iblock) = parent(1,lb)+maxblocks*parent(2,lb)
         lref_min = min(lrefine(lb),lref_min)
         lref_max = max(lrefine(lb),lref_max)
       endif

       enddo
       endif
       lref_min = min(lref_max,lref_min)

       lref_maxt = lref_max
       lref_mint = lref_min
       call shmem_barrier_all()
       call comm_int_max_to_all (lref_max,lref_maxt)
       call shmem_barrier_all()
       call comm_int_min_to_all (lref_min,lref_mint)
       call shmem_barrier_all()


       if(iblock.gt.nchild) then
! sort the list of newchildren according to their parents ids
! This will enable us to avoid costly extra guardcell filling operations
! on parent blocks
         call quicksort_index(iblock,parent_id,indx)
       endif

       p_cache_addr(:) = -1

       call shmem_barrier_all()

! prolongation must be applied in descending order of refinement level
! so that the case where neighboring blocks at different refinement
! level are required to be refined at the same time is handled correctly.
       do lreflevel = lref_min,lref_max

        call amr_1blk_guardcell_reset
        call shmem_barrier_all()
#ifdef NO_PERMANENT_GUARDCELLS
! Store a copy of the current solution in gt_unk
        level = -1                      ! copy all refinement levels
        call amr_1blk_copy_soln(level)
        call shmem_barrier_all()
#endif


#ifndef NO_PERMANENT_GUARDCELLS
#ifdef FORCE_CONSISTENCY_AT_SRL_INTERFACES
      if(lfc) then  
! If using facevars and permanent guardcells then we will need to ensure
! that gt_facevarx, etc are filled before prolongation begins
      if(nfacevar.gt.0) then
        do lb = 1,lnblocks
          gt_facevarx(:,1,:,:,lb) = facevarx(:,1+nguard,:,:,lb)
          gt_facevarx(:,2,:,:,lb) = facevarx(:,nxb+1+nguard,:,:,lb)
          if(ndim.ge.2) then
          gt_facevary(:,:,1,:,lb) = facevary(:,:,1+nguard*k2d,:,lb)
          gt_facevary(:,:,1+k2d,:,lb) =
     .                        facevary(:,:,nyb+(1+nguard)*k2d,:,lb)
          endif
          if(ndim.eq.3) then
          gt_facevarz(:,:,:,1,lb) = facevarz(:,:,:,1+nguard*k3d,lb)
          gt_facevarz(:,:,:,1+k3d,lb) =
     .                        facevarz(:,:,:,nzb+(1+nguard)*k3d,lb)
          endif
        enddo
      endif
      endif
      call shmem_barrier_all()
#endif  /* FORCE_CONSISTENCY_AT_SRL_INTERFACES */
#endif  /* NO_PERMANENT_GUARDCELLS */


       if(iblock.gt.0) then
       do lbi = 1,iblock
       lb = index(indx(lbi))

       if(lrefine(lb).eq.lreflevel) then


! compute offset for child cell inside parent
       ioff = mod(which_child(lb)-1,2)*nxb/2
       joff = mod((which_child(lb)-1)/2,2)*nyb/2
       koff = mod((which_child(lb)-1)/4,2)*nzb/2

! get address of parent block
       parent_blk = parent(1,lb)
       parent_pe  = parent(2,lb)

! Is parent data currently cached?
       if( parent_blk.ne.p_cache_addr(1) .or.
     .     parent_pe.ne.p_cache_addr(2) ) then


! Fetch parent data block into layer 1 of 1blk data structure and fill
! its guardcells 
         ldiag = .false.
#ifdef DIAGONALS
         ldiag = .true.
#endif
!        l_srl_only = .false.
         l_srl_only = .true.           ! seems to be all that is required ??
         icoord = 0

         call amr_1blk_guardcell(mype,iopt,nlayers,parent_blk,
     .                           parent_pe,lcc,lfc,lec,lnc,
     .                           l_srl_only,
     .                           icoord,ldiag)


! update address of cached parent
         p_cache_addr(1) = parent_blk
         p_cache_addr(2) = parent_pe

!      endif                       ! end of p_cache_addr if test
       else
       endif


! Prolong data from working block to the new child block
        idest = 2

! cell-centered data
       if(lcc) then

         if(iopt.eq.1) then
           recv1(:,:,:,:) = unk1(:,:,:,:,1)
 

           call amr_1blk_cc_prol_gen_unk_fun(recv1,
     .              ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff,mype,
     .              lb,parent_pe,parent_blk)


         elseif(iopt.ge.2) then
           recvw1(:,:,:,1) = work1(:,:,:,1)
           call amr_1blk_cc_prol_gen_work_fun(
     .              recvw1(:,:,:,1),
     .              ia,ib,ja,jb,ka,kb,idest,ioff,joff,koff,mype,
     .              lb,parent_pe,parent_blk,interp_mask_work(iopt-1))
         endif                           ! end of iopt if test
       endif                             ! end of lcc if test


#ifdef MLRDBG
       print *,'MLRDBG: 1000'
#endif

! cell face-centered data 
       if(lfc) then 
!mlr begin
          ! do ibndvar=1,nbndvar...
          if(
     .         prol_fc_dbz
     .         )then     
                                ! prolong divergenceless b on fc

             if(ndim.ne.3)then
                print *,'map:dbz:ndim.ne.3! halting...'
                stop
             endif
             
             recvfx(1:nfacevar, il_bnd1:iu_bnd1+1,
     .            jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1)
     .          = facevarx1(1:nfacevar, il_bnd1:iu_bnd1+1,
     .            jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1, 1)

             recvfy(1:nfacevar, il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,
     .            kl_bnd1:ku_bnd1)
     .          = facevary1(1:nfacevar, il_bnd1:iu_bnd1,
     .            jl_bnd1:ju_bnd1+k2d,kl_bnd1:ku_bnd1, 1)

             recvfz(1:nfacevar, il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,
     .            kl_bnd1:ku_bnd1+k3d)
     .          = facevarz1(1:nfacevar, il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,
     .            kl_bnd1:ku_bnd1+k3d, 1)


       do iprol = 1, prol_fc_dbz_n
          iv1 = prol_fc_dbz_ivar(1,iprol)
          iv2 = prol_fc_dbz_ivar(2,iprol)
          iv3 = prol_fc_dbz_ivar(3,iprol)
             
             call amr_1blk_fc_prol_dbz( 
     .            recvfx, recvfy, recvfz,       
     .            nfacevar, iv1, iv2, iv3,
     .            ia,ib,ja,jb,ka,kb,   
     .            idest,ioff,joff,koff,         
     .            mype,lb,parent_pe,parent_blk  
     .            )

      end do
          
       end if ! prol_fc_dbz

!          else
                                ! prolong using other methods
!mlr end
! x-face
         recvf(1:nfacevar, il_bnd1:iu_bnd1+1,
     .               jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1)
     .      = facevarx1(1:nfacevar, il_bnd1:iu_bnd1+1,
     .               jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1, 1)

         call amr_1blk_fc_prol_gen_fun(recvf,
     .              ia,ib+1,ja,jb,ka,kb,idest,ioff,joff,koff,
     .              mype,lb,parent_pe,parent_blk,1)

!!!#if N_DIM > 1
         if (ndim > 1) then
! y-face
         recvf(1:nfacevar, il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1+k2d,
     .              kl_bnd1:ku_bnd1)
     .      = facevary1(1:nfacevar, il_bnd1:iu_bnd1,
     .             jl_bnd1:ju_bnd1+k2d,kl_bnd1:ku_bnd1, 1)
         call amr_1blk_fc_prol_gen_fun(recvf,
     .              ia,ib,ja,jb+1,ka,kb,idest,ioff,joff,koff,
     .              mype,lb,parent_pe,parent_blk,2)
!!!#endif
         end if

!!!#if N_DIM == 3
         if (ndim == 3) then
! z-face
         recvf(1:nfacevar, il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,
     .              kl_bnd1:ku_bnd1+k3d)
     .      = facevarz1(1:nfacevar, il_bnd1:iu_bnd1,jl_bnd1:ju_bnd1,
     .              kl_bnd1:ku_bnd1+k3d, 1)
         call amr_1blk_fc_prol_gen_fun(recvf,
     .              ia,ib,ja,jb,ka,kb+1,idest,ioff,joff,koff,
     .              mype,lb,parent_pe,parent_blk,3)
!!!#endif

!+ CLEAN_DIVB
#ifdef CLEAN_DIVB
#ifdef MLRDBG
       print *,'MLRDBG: 1009'
#endif

       if(prol_fc_clean_divb)then

#ifdef MLRDBG
       print *,'MLRDBG: 1010'
#endif

          call amr_1blk_fc_clean_divb( 
     &         nfacevar,
     &         ia,ib,ja,jb,ka,kb,
     &         0, 0, 0, 0, 0, 0,
     &         idest,ioff,joff,koff,
     &         mype,lb,parent_pe,parent_blk
     &         )

#ifdef MLRDBG
       print *,'MLRDBG: 1020'
#endif

       end if
#endif
!- CLEAN_DIVB

         end if

        endif                     ! end of lfc if test


!!!#if N_DIM > 1
       if (ndim > 1) then
! cell edge-centered data
       if(lec) then

! x-edge
         recve(1:nvaredge, il_bnd1:iu_bnd1,
     .               jl_bnd1:ju_bnd1+k2d,kl_bnd1:ku_bnd1+k3d)
     .      = unk_e_x1(1:nvaredge, il_bnd1:iu_bnd1,
     .               jl_bnd1:ju_bnd1+k2d,kl_bnd1:ku_bnd1+k3d, 1)
        
         call amr_1blk_ec_prol_gen_fun(recve,
     .                                 ia,ib,ja,jb+k2d,ka,kb+k3d,
     .                                 idest,ioff,joff,koff,mype,1)

         
! y-edge
         recve(1:nvaredge, il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1,
     .              kl_bnd1:ku_bnd1+k3d)
     .      = unk_e_y1(1:nvaredge, il_bnd1:iu_bnd1+1,
     .             jl_bnd1:ju_bnd1,kl_bnd1:ku_bnd1+k3d, 1)
         call amr_1blk_ec_prol_gen_fun(recve,
     .                                 ia,ib+1,ja,jb,ka,kb+k3d,
     .                                 idest,ioff,joff,koff,mype,2)

!!!#if N_DIM == 3
         if (ndim == 3) then
! z-edge
         recve(1:nvaredge, il_bnd1:iu_bnd1+1,jl_bnd1:ju_bnd1+k2d,
     .              kl_bnd1:ku_bnd1)
     .      = unk_e_z1(1:nvaredge, il_bnd1:iu_bnd1+1,
     .             jl_bnd1:ju_bnd1+k2d,kl_bnd1:ku_bnd1, 1)
         call amr_1blk_ec_prol_gen_fun(recve,
     .                                 ia,ib+1,ja,jb+k2d,ka,kb,
     .                                 idest,ioff,joff,koff,mype,3)
!!!#endif
         end if

        endif                     ! end of lec if test
!!!#endif
        end if

! cell corner data
       if(lnc) then

           recvn(:,:,:,:) = unk_n1(:,:,:,:,1)


           call amr_1blk_nc_prol_gen_fun(recvn,
     .                                   ia,ib+1,ja,jb+k2d,ka,kb+k3d,
     .                                   idest,ioff,joff,koff,mype)


       endif                             ! end of lnc if test

! copy data back to permanent storage arrays
        call amr_1blk_to_perm( lcc,lfc,lec,lnc,lb,iopt,idest )

        

#ifdef MULTIBLOCK_TEST_2D
!-- i/o
        if(
     .      ( abs(coord(1,lb)-.625).lt..01)
     .    .and.
     .      abs(coord(2,lb)-.375).lt..01 .and.
     .      lrefine(lb).eq.lreflevel ) then
          do j=1+nguard*npgs,nyb+nguard*npgs
           write(*,*) 'before modif facevarx blk ',lb,mype,' j ',j,
     .                 facevarx(1,1+nguard*npgs,j,1,lb)
     .                 ,' instance ',instance
          enddo
          do j=1+nguard*npgs,nyb+nguard*npgs
           write(*,*) 'before modif facevary blk ',lb,mype,' j ',j,
     .                 facevary(1,1+nguard*npgs,j,1,lb)
     .                 ,' instance ',instance
          enddo
        endif

        if( 
     .       ( abs(coord(1,lb)-.375).lt..01) 
     .    .and.
     .      abs(coord(2,lb)-.375).lt..01 .and.
     .      lrefine(lb).eq.lreflevel ) then
          do j=1+nguard*npgs,nyb+nguard*npgs
           write(*,*) 'before modif facevarx blk ',lb,mype,' j ',j,
     .                 facevarx(1,nxb+1+nguard*npgs,j,1,lb)
          enddo
        endif
!-- i/o

! introduce an error in facevars on the face of a block to test
! whether the divergence consistency forcing and the
! srl forced consistency are working together properly.
        if(lfc .and. abs(coord(1,lb)-.625).lt..01 .and.
     .      abs(coord(2,lb)-.375).lt..01 .and.
     .      lrefine(lb).eq.lreflevel) then
          do j = 1+nguard*npgs,nyb+nguard*npgs,2
            facevarx(1,1+nguard*npgs,j,1,lb) =
     .       facevarx(1,1+nguard*npgs,j,1,lb) + 1.e-4
            facevarx(1,1+nguard*npgs,j+1,1,lb) =
     .       facevarx(1,1+nguard*npgs,j+1,1,lb) - 1.e-4
            facevary(1,1+nguard*npgs,j+1,1,lb) =
     .       facevary(1,1+nguard*npgs,j+1,1,lb) + 1.e-4
          enddo
        endif


!-- i/o
        if( abs(coord(1,lb)-.625).lt..01 .and.
     .      abs(coord(2,lb)-.375).lt..01 .and.
     .      lrefine(lb).eq.lreflevel) then
          do j=1+nguard*npgs,nyb+nguard*npgs
           write(*,*) 'before divbcon facevary blk ',lb,' j ',j,
     .                 facevary(1,1+nguard*npgs,j,1,lb)
          enddo
          do j=1+nguard*npgs,nyb+nguard*npgs
           write(*,*) 'before divbcon facevarx blk ',lb,' j ',j,
     .                 facevarx(1,1+nguard*npgs,j,1,lb),
     .                ' facevarx blk ',neigh(1,1,lb),
     .             facevarx(1,nxb+1+nguard*npgs,j,1,neigh(1,1,lb))
          enddo
          do j=1+nguard*npgs,nyb+nguard*npgs
           write(*,*) 'before divbcon facevarx blk ',lb,' i = 2 j ',j,
     .                 facevarx(1,2+nguard*npgs,j,1,lb)
          enddo
          do j=1+nguard*npgs,nyb+nguard*npgs
           write(*,*) 'before divbcon blk ',lb,' j ',j,' divb ',
     .    facevarx(1,2+nguard*npgs,j,1,lb)
     .    - facevarx(1,1+nguard*npgs,j,1,lb)
     .    + facevary(1,1+nguard*npgs,j+1,1,lb)
     .    - facevary(1,1+nguard*npgs,j,1,lb)
          enddo
        endif
!-- i/o
#endif /* MULTIBLOCK_TEST_2D */


      endif                       ! end of lrefine if test

      enddo                       ! end of loop over blocks
      endif                       ! end of lnblocks if test


! Ensure new blocks inherit data on block face shared with an old
! existing neighbor, instead of filling from parent by interpolation.
#ifdef DIVERGENCE_FREE
        if(lfc) then
        do nfield = 1,nfield_divf
          call amr_prolong_fc_divbconsist(mype,lreflevel,
     .                              nfield)
        enddo
        endif
#endif /* DIVERGENCE_FREE */



#ifdef MULTIBLOCK_TEST_2D
!-- i/o
        iib = -1
        do ilb = 1,lnblocks
        if( abs(coord(1,ilb)-.625).lt..01 .and.
     .      abs(coord(2,ilb)-.375).lt..01 .and.
     .      lrefine(ib).eq.lreflevel) iib = ilb
        enddo
        if(iib.gt.0) then
          write(*,*) 'lreflevel ',lreflevel
          do j=1+nguard*npgs,nyb+nguard*npgs
           write(*,*) 'after divbcon facevary blk ',iib,' j ',j,
     .                 facevary(1,1+nguard*npgs,j,1,iib),' divb ',
     .    facevarx(1,2+nguard*npgs,j,1,iib)
     .    - facevarx(1,1+nguard*npgs,j,1,iib)
     .    + facevary(1,1+nguard*npgs,j+1,1,iib)
     .    - facevary(1,1+nguard*npgs,j,1,iib)
          enddo
          do j=1+nguard*npgs,nyb+nguard*npgs
           write(*,*) 'after divbcon facevarx blk ',iib,' j ',j,
     .                 facevarx(1,1+nguard*npgs,j,1,iib),
     .                ' facevarx blk ',neigh(1,1,iib),
     .           facevarx(1,nxb+1+nguard*npgs,j,1,neigh(1,1,iib))
          enddo
        endif
!-- i/o
#endif /* MULTIBLOCK_TEST_2D */



      enddo                       ! end of loop over levels

! reset cache addresses
      call amr_1blk_guardcell_reset

!------------------------------------

      call shmem_barrier_all()
      newchild(:) = .false.
      call shmem_barrier_all()


#ifdef TIMINGS
      itimer2 = irtc()
      call shmem_barrier_all()
      irtc_prolong = itimer2-itimer1+irtc_prolong
#endif
      return
      end subroutine amr_prolong
