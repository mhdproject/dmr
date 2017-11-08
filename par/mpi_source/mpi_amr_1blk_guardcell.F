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
!#define DEBUGX

      subroutine amr_1blk_guardcell(
     .                              mype,iopt,nlayers,lb,pe,
     .                              lcc,lfc,lec,lnc,
     .                              l_srl_only,icoord,ldiag,
     .                              nlayersx,nlayersy,nlayersz)




!------------------------------------------------------------------------
!
! This routine manages the transfer of guard cell data for a
! specific single block. It uses the morton numbering scheme for the
! determination of the neighbor relations and differs in that from
! the older routine amr_1blk_guardcell.
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
! Written :     Michael Gehmeyr & Peter MacNeice         November 1999
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

      use paramesh_dimensions
      use physicaldata
      use tree
      use timings
      use workspace
      use mpi_morton

      use paramesh_mpi_interfaces, only : 
     &                                mpi_amr_local_surr_blks_lkup,
     &                                mpi_amr_1blk_guardcell_c_to_f,
     &                                mpi_amr_get_remote_block

      use paramesh_interfaces, only : amr_perm_to_1blk,
     &                                amr_1blk_guardcell_srl,
     &                                amr_1blk_guardcell_f_to_c

      implicit none

#include "amr_shmem.fh"
      include 'mpif.h'

      logical, intent(in) :: lcc,lfc,lec,lnc,l_srl_only,ldiag
      integer, intent(in) :: mype,lb,pe,iopt,nlayers,icoord
      integer, intent(in), optional :: nlayersx,nlayersy,nlayersz

!------------------------------------

! local arrays and variables

      integer ::  nprocs
      integer, external :: shmem_n_pes

      integer :: parent_lb,parent_pe
      integer :: i,j,k, ll, idest, iblock

      logical :: lcoarse,l_parent
      logical :: loc_lcc,loc_lfc,loc_lec,loc_lnc
      logical :: lfound

      integer :: surrblks(3,3,3,3), tsurrblks(3,3,3,3)
      integer :: psurrblks(3,3,3,3)
      integer :: pcache_pe,pcache_blk
      integer :: ierrorcode,ierr

      integer :: icoord_loc
      logical :: ldiag_loc

      integer :: nlayers0x, nlayers0y, nlayers0z
      integer :: nguard0,nguard_work0
      parameter (nguard0 = nguard*npgs )
      parameter (nguard_work0 = nguard_work*npgs )
      logical, save :: first_cc = .true.
      logical, save :: first_nc = .true.
      logical, save :: first_ec = .true.
      logical, save :: first_fc = .true.


#ifdef TIMING_MPI
      double precision :: time1
      double precision :: time2
#endif /* TIMING_MPI */

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
#ifdef TIMING_MPI
      time1 = mpi_wtime()
      time2 = mpi_wtime()
#endif /* TIMING_MPI */


      nprocs = shmem_n_pes()

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

      if(pe.ne.mype) then
          write(*,*) 'Error : trying to fill guardcells for a ',
     .               'remote block - not supported with the mpi ',
     .               'version. '
          call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif
      if( 
     .    mpi_pattern_id.eq.10  .or.
     .    (mpi_pattern_id.eq.20.and.lprolong_in_progress) .or.
     .    (mpi_pattern_id.eq.40.and.lrestrict_in_progress)
     .                            ) then
      elseif(nprocs.gt.1) then
        write(*,*) 'Paramesh error : amr_1blk_guardcell : ',
     . ' wrong pattern being',
     . ' used for pre-communication for guardcell fill : Fix',
     . ' - insert appropriate call to mpi_amr_comm_setup '
     . ,' mpi_pattern_id ',mpi_pattern_id,
     . ' lprolong_in_progress ',lprolong_in_progress,
     . ' lrestrict_in_progress ',lrestrict_in_progress
        call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
      endif

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
            write(*,*) 'Paramesh Error: for a call to',
     .       '  mpi_amr_1blk_guardcell',
     .       ' one of more of the arguments lcc/lfc/lec/lnc are not',
     .       ' consistent with nvar/nfacevar/nvaredge/nvarcorn.'
          endif
          call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
        endif
#endif /* AMR_ERROR_CHECKING */


      surrblks  = -1          ! catchy initialization
      tsurrblks = -1  

! construct a list of blocks surrounding local block lb
      l_parent = .true.
      if(l_srl_only) l_parent = .false.
      call mpi_amr_local_surr_blks_lkup(mype,lb,
     &                          surrblks,l_parent,psurrblks)

! relate surrblks with the guard block indices stored implicitly 
! in laddress and update tsurrblks

      tsurrblks = surrblks    ! set up temporary array
#ifdef DEBUG
       if(mype.eq.0)
     .        write(*,*) '1blk_guardcell : pe ',mype,
     .         ' working on lb ',lb,' tsurrblks ',
     .         tsurrblks(:,:,:,2),
     .         ' lfc ',lfc
#endif /*DEBUG */


#ifdef TIMING_MPI
      time2 = mpi_wtime()
#endif /* TIMING_MPI */
! guard block indeces start at strt_buffer after lnblocks, and end at
! last_buffer as determined in subroutine mpi_commatrix.


       do k = 2-k3d,2+k3d     ! loop over all its surrounding blocks
       do j = 2-k2d,2+k2d
       do i = 1,3
         if( surrblks(2,i,j,k).ne.mype .and.
     .       surrblks(2,i,j,k).ge.0 ) then

!-
         lfound = .false.
         ll = ladd_strt(surrblks(2,i,j,k))
         do while(.not.lfound.and.ll.le.ladd_end(surrblks(2,i,j,k)))
!         ll = strt_buffer
!         do while(.not.lfound.and.ll.le.last_buffer)
           if( (surrblks(2,i,j,k).eq.laddress(2,ll))  .and.
     .         (surrblks(1,i,j,k).eq.laddress(1,ll)) ) then
                              ! found the corresponding block id ll
             tsurrblks(1,i,j,k) = ll
             tsurrblks(2,i,j,k) = mype
             tsurrblks(3,i,j,k) = nodetype(ll)     !?????
             lfound = .true.
#ifdef DEBUG
       if(mype.eq.0)
     .       write(*,*) 'pe ',mype,' looking for ',surrblks(:,i,j,k)
     .             ,' in slot ',ll,' FOUND ',
     .             laddress(:,ll)
#endif /*DEBUG */
           else
#ifdef DEBUG
       if(mype.eq.0)
     .       write(*,*) 'pe ',mype,' looking for ',surrblks(:,i,j,k)
     .             ,' in slot ',ll,' found ',
     .             laddress(:,ll)
#endif /*DEBUG */
             ll = ll+1  
           endif
         enddo
!-

         endif
         if( (tsurrblks(2,i,j,k).ne.mype) .and.
     .       (tsurrblks(2,i,j,k).ne.-1) .and.
     .       (tsurrblks(1,i,j,k).gt.-20) ) then
             write(*,*) '1blk_guardcell : pe ',mype,
     .         ' working on lb ',lb,' neigh ',i,j,k,
     .         ' cannot find surrblk ',
     .         surrblks(:,i,j,k),' on this proc ',
     .         ' laddress ',laddress(:,strt_buffer:last_buffer),
     .         ' strt_buffer,last_buffer ',strt_buffer,last_buffer,
     .         ' tsurrblks ',tsurrblks(1:2,i,j,k)
     .       ,' ladd_strt ',ladd_strt,' ladd_end ',ladd_end
             call mpi_abort(MPI_COMM_WORLD,ierrorcode,ierr)
         endif
        enddo
        enddo
        enddo


! update surrblks with the local guard block info

      surrblks = tsurrblks

#ifdef DEBUG
       if(mype.eq.0) then
        print *, 'MPI SURRBLKS on proc=',mype,' and block=',lb
        do j = 2+k2d,2-k2d,-1
        write(*,'(i4,2x,i4,2x,i4)') (surrblks(1,i,j,2),i=1,3) 
        enddo
       endif
#endif

#ifdef TIMING_MPI
      timer_amr_1blk_guardcell(1) = timer_amr_1blk_guardcell(1)
     .                           +  mpi_wtime() - time2
      time2 = mpi_wtime()
#endif /* TIMING_MPI */
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if(.not.l_srl_only) then
!
! Are there any coarse neighbors?
        lcoarse = .false.


        if(parent(1,lb).gt.0) then
        do k = 2-k3d,2+k3d
        do j = 2-k2d,2+k2d
        do i = 1,3
          if(surrblks(1,i,j,k).gt.-20.and.surrblks(1,i,j,k).lt.0)
     .                   lcoarse = .true.
        enddo
        enddo
        enddo
        endif 

        endif                        ! end of l_srl_only if test

!------------------------------------
!
! Put leaf block lb's data into the data_1blk.fh datastructures, 
! with the appropriate guardcell padding.

          idest = 1

          call amr_perm_to_1blk(lcc,lfc,lec,lnc,lb,pe,iopt,idest)

#ifdef DEBUG
       if(mype.eq.0) then
          write(*,*) 'after perm : pe ',mype,' blk ',lb,pe,
     .    ' facevarx1 ',facevarx1(1,:,:,1,1)
       endif
#endif /* DEBUG */

!------------------------------------
!
        if(iopt.eq.1) then
          pcache_pe  = pcache_pe_u
          pcache_blk = pcache_blk_u
        elseif(iopt.ge.2) then
          pcache_pe  = pcache_pe_w
          pcache_blk = pcache_blk_w
        endif

        if(.not.l_srl_only) then
        if(lcoarse) then

!
! Put data from lb's parent into the data_1blk.fh datastructures, with the
! appropriate guardcell padding. Check to see if data is currently cached.
          parent_lb = parent(1,lb)
          parent_pe = parent(2,lb)


#ifdef DEBUG
       if(mype.eq.0) then
        write(*,*) 'pe ',mype,'blk ',lb,
     .      'parent blk ',parent_lb,parent_pe,' cached ',
     .      pcache_blk,pcache_pe
        endif
#endif
        if( (parent_lb.gt.0) .and.
     .      ((parent_lb.ne.pcache_blk).or.(parent_pe.ne.pcache_pe) )
     .      ) then

! record id of new parent block placed in cache
          lnew_parent = .true.
          pcache_blk = parent_lb
          pcache_pe  = parent_pe

#ifdef DEBUG
       if(mype.eq.0) then
        write(*,*) 'pe ',mype,'blk ',lb,
     .      'fetching parent blk ',parent_lb,parent_pe,
     .      'new cache ',pcache_blk,pcache_pe,' lnew_parent ',
     .      lnew_parent,' parent(:,lb) ',parent(:,lb)
        endif
#endif

        if(lcc) then
           if (first_cc) then
              unk1(:,:,:,:,2) = 0.
              work1(:,:,:,2) = 0.
              first_cc = .false.
           end if
        end if
        if(lnc) then
           if (first_nc) then
              unk_n1(:,:,:,:,2) = 0.
              first_nc = .false.
           end if
        end if
        if(lec) then
           if (first_ec) then
             unk_e_x1(:,:,:,:,2) = 0.
             unk_e_y1(:,:,:,:,2) = 0.
             unk_e_z1(:,:,:,:,2) = 0.
             first_ec = .false.
           endif
        endif
        if(lfc) then
           if (first_fc) then
             facevarx1(:,:,:,:,2) = 0.
             facevary1(:,:,:,:,2) = 0.
             facevarz1(:,:,:,:,2) = 0.
             first_fc = .false.
           endif
        end if

        idest = 2
        call mpi_amr_get_remote_block(mype,parent_pe,parent_lb,
     .                                idest,iopt,lcc,lfc,lec,lnc,
     .                            nlayers0x,nlayers0y,nlayers0z)

!
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
     .                              nlayers0x,nlayers0y,nlayers0z)

        endif       ! end if parents data not previously cached
!
!------------------------------------
!
! Do guardcell filling from coarse neigbors into the current block


        call mpi_amr_1blk_guardcell_c_to_f( mype,lb,pe,iopt,nlayers,
     .                                      surrblks,
     .                                      lcc,lfc,lec,lnc,
     .                                      icoord,ldiag,
     .                                      nlayers0x,
     .                                      nlayers0y,
     .                                      nlayers0z)

!------------------------------------

        endif                       ! end of lcoarse if test

        endif                       ! end of l_srl_only if test

#ifdef TIMING_MPI
      timer_amr_1blk_guardcell(2) = timer_amr_1blk_guardcell(2)
     .                           +  mpi_wtime() - time2
      time2 = mpi_wtime()
#endif /* TIMING_MPI */
!
!------------------------------------
!
! Do guardcell filling from any surrounding blocks at the same refinement
! level as block lb.
        iblock = 1
#ifdef DEBUG
          if(iopt.gt.1) 
     .     write(*,*) 'in guardcell_srl : pe ',mype,' blk ',lb,
     .        ' surrblks(:,:,:,2) ',surrblks(:,:,:,2)
#endif /*DEBUG */

        call amr_1blk_guardcell_srl(mype,mype,lb,
     .                              iblock,iopt,nlayers,surrblks,
     .                              lcc,lfc,lec,lnc,
     .                              icoord,ldiag,
     .                              nlayers0x,nlayers0y,nlayers0z)

#ifdef TIMING_MPI
      timer_amr_1blk_guardcell(3) = timer_amr_1blk_guardcell(3)
     .                           +  mpi_wtime() - time2
#endif /* TIMING_MPI */

#ifndef ADVANCE_ALL_LEVELS
!        call amr_1blk_guardcell_f_to_c(mype,pe,lb,
!     .                              iblock,iopt,nlayers,surrblks,
!     .                              lcc,lfc,lec,lnc,icoord,ldiag,
!     .                              nlayers0x,nlayers0y,nlayers0z)
#endif

#ifdef DEBUG
          write(*,*) 'in guardcell_srl : pe ',mype,' blk ',lb,
     .        ' surrblks(:,:,:,2) ',surrblks(:,:,:,2)
#endif /*DEBUG */
    
#ifdef DEBUG
          write(*,*) 'after guardcell_srl : pe ',mype,' blk ',lb,
     .    ' unk1 ',unk1(1,:,:,1,1)
#endif /*DEBUG */

        if(iopt.eq.1) then
          pcache_pe_u  = pcache_pe
          pcache_blk_u = pcache_blk
        elseif(iopt.ge.2) then
          pcache_pe_w  = pcache_pe
          pcache_blk_w = pcache_blk
        endif

#ifdef TIMING_MPI
      timer_amr_1blk_guardcell(0) = timer_amr_1blk_guardcell(0)
     .                           +  mpi_wtime() - time1
#endif /* TIMING_MPI */

      return
      end subroutine amr_1blk_guardcell



