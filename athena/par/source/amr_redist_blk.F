!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"


      subroutine amr_redist_blk(mype,new_loc,lnblocks_old)




!------------------------------------------------------------------------
!
! This routine redistributes blocks to their new locations
! stored in the array new_loc(2,:). The elements new_loc(1,:)
! store the new local address and new_loc(2,:) store the 
! new processor address of the block.
!
! It basically works by copying any blocks which are destined 
! to be sent off-processor, to storage blocks at the end of 
! the original list of blocks. Then any on-processor movement
! of blocks is done. Finally the off-processor destined blocks
! are put into their correct places. The routine cycles over
! these steps as many times as is necessary to cope with cases 
! in which a shortage of storage space for the blocks means 
! that not all the blocks get sent in the first sweep.
!
! Written :       Peter MacNeice December 1996
!------------------------------------------------------------------------
!
! Arguments:
!      mype                      local processor number
!      new_loc(2,maxblocks_tr)   new locations for grid blocks.
!      lnblocks_old              the value of lnblocks before new child
!                                 blocks were spawned during this refinement
!                                 step.
!------------------------------------
!#define      DIAG

      use paramesh_dimensions
      use physicaldata
      use tree

      use paramesh_interfaces, only :  comm_int_max_to_all

      include  'tracer.fh'

#include "amr_shmem.fh"
      include  'shmem_reduce.fh'


!      integer, intent(inout) ::  new_loc(2,maxblocks_tr)
      integer, intent(inout) ::  new_loc(:,:)
      integer, intent(in)    ::  lnblocks_old
      integer, intent(in)    ::  mype


!------------------------------------
! local arrays

! vector to mark blocks with good data
      integer ialive(maxblocks_alloc)
      integer kalive(maxblocks_alloc)


! new_loct is a static copy of new_loc required for use in shmem calls
      integer new_loct(2,maxblocks_alloc)
      integer src_loc(2,maxblocks_alloc),source_id(2)


! indirection index for accessing the storage blocks
      integer indx(maxblocks_alloc)
      integer icount,igcount

      integer shmem_n_pes

      save icount,igcount,ialive,jalive
      save new_loct,src_loc,source_id


      integer :: jalive,nprocs,i,isg,istart_store,nstore,istore
      integer :: iloop,i_transfer,ilim,i_curr_blk,i_store_blk

      integer,save :: idataout(1)

!------------------------------------


      nprocs = shmem_n_pes()


      if(lnblocks.gt.maxblocks) write(*,*) 'pe ',mype,
     .      ' lnblocks too large ',lnblocks,
     .      ' lnblocks_old ',lnblocks_old


! initialize the marker array
      ialive(1:maxblocks_alloc) = -1
      pesource(:,:) = -1
      src_loc(:,:) = -1
      new_loct(:,:) = -1
      do i=1,lnblocks_old
      new_loct(1,i) = new_loc(1,i)
      new_loct(2,i) = new_loc(2,i)
      enddo
      source_id(:)=0


! mark live blocks
      if(lnblocks.gt.0) then
        do isg = 1,lnblocks_old
          ialive(isg) = 0
        enddo
      endif


! where is the first extra block available for use as storage?
      istart_store = max(lnblocks_old,new_lnblocks) + 1

#ifdef DIAG
           write(*,*) 'pe ',mype,' lnblocks_old ',lnblocks_old,
     .       ' lnblocks ',lnblocks
           write(*,*) 'pe ',mype,' new_lnblocks ',new_lnblocks
           write(*,*) 'pe ',mype,' istart_store ',istart_store
      call shmem_barrier_all()
#endif

! How many blocks can be stored at the end of the list?
!      nstore = maxblocks_alloc - istart_store + 1
      nstore = maxblocks - istart_store + 1

! initialize indirection index for storage blocks
      do i= 1,nstore
        indx(i) = istart_store + i - 1
      enddo

      istore = 0

! set icount to a value greater than 0 to ensure that the processor
! executes at least one sweep through the list of blocks.
      icount = 1
      iloop  = 0

      call shmem_barrier_all()

!------------------------------------
! start of loop over steps
10    continue
      iloop = iloop+1
      if(iloop.gt.10.and.mype.eq.0) then
       write(*,*) 'PARAMESH ERROR !'
       write(*,*) 'redist_blk : Error '
       write(*,*) 'Movement of data blocks failed to'
       write(*,*) 'complete in 10 cycles. This is probably'
       write(*,*) 'because MAXBLOCKS is set too small.'
       call amr_abort
      endif
      if(icount.eq.0) goto 20


! How many blocks need to be transferred ?
      i_transfer = 0
!      do isg = 1,maxblocks_alloc
      do isg = 1,maxblocks
         if( ( (new_loct(1,isg).ne.isg).or.(new_loct(2,isg).ne.mype) )
     .          .and. (ialive(isg).eq.0) ) 
     .       i_transfer = i_transfer + 1
      enddo

#ifdef DIAG
           write(*,*) 'pe ',mype,' loop ',iloop,' i_transfer ',
     .       i_transfer,' istart_store ',istart_store
#endif

! If there are any blocks to be moved then
      if(i_transfer.gt.0) then

! initialize current storage block number
      istore = 1

! loop over original list of blocks, moving their data into the next
! temporary storage slot, marking that slot alive and the original
! location dead.
         ilim = max(1,istart_store-1)
         do i_curr_blk = 1,ilim


              if( ((new_loct(2,i_curr_blk).ne.mype).or. 
     .            (new_loct(1,i_curr_blk).ne.i_curr_blk)) .and.
     .            (ialive(i_curr_blk).eq.0)  ) then

                if(istore.le.nstore) then
                  i_store_blk = indx(istore)
                  if(ialive(i_store_blk).eq.-1) then

                    if(nvar.gt.0) unk(:,:,:,:,i_store_blk) = 
     .                            unk(:,:,:,:,i_curr_blk)
                    if(nfacevar.gt.0) then
                      facevarx(:,:,:,:,i_store_blk) =
     .                        facevarx(:,:,:,:,i_curr_blk)
                      if(ndim.ge.2) facevary(:,:,:,:,i_store_blk) =
     .                        facevary(:,:,:,:,i_curr_blk)
                      if(ndim.eq.3) facevarz(:,:,:,:,i_store_blk) =
     .                        facevarz(:,:,:,:,i_curr_blk)
                    endif
                    if(nvaredge.gt.0) then
                      unk_e_x(:,:,:,:,i_store_blk) =
     .                        unk_e_x(:,:,:,:,i_curr_blk)
                      unk_e_y(:,:,:,:,i_store_blk) =
     .                        unk_e_y(:,:,:,:,i_curr_blk)
                      unk_e_z(:,:,:,:,i_store_blk) =
     .                        unk_e_z(:,:,:,:,i_curr_blk)
                    endif
                    if(nvarcorn.gt.0) unk_n(:,:,:,:,i_store_blk) =
     .                                unk_n(:,:,:,:,i_curr_blk)
                    ialive(i_curr_blk) = -1
                    ialive(i_store_blk) = 0
                    new_loct(1,i_store_blk) = new_loct(1,i_curr_blk)
                    new_loct(2,i_store_blk) = new_loct(2,i_curr_blk)
                    new_loct(1,i_curr_blk) = -1
                    new_loct(2,i_curr_blk) = -1

                    istore = istore + 1
                  endif
                endif
              endif

         enddo

      endif

20    continue

! Initialize temporary used to record changes to ialive.
      kalive(:) = 0


! synchronize all processors.
      call shmem_barrier_all()


! Now loop over the entire list of blocks sending a notice to the
! destination block that it should fetch data and giving the location
! of the source of the data.
!      do isg = 1,maxblocks_alloc
      do isg = 1,maxblocks

! Does this block need to be sent?
      if(ialive(isg).eq.0) then
        if( (new_loct(1,isg).ne.isg).or.
     .               (new_loct(2,isg).ne.mype) ) then

#ifdef DIAG
      if(new_loct(1,isg).gt.maxblocks.or.new_loct(1,isg).le.0
     .      .or.new_loct(2,isg).ge.nprocs.or.new_loct(2,isg).lt.0) then
           write(*,*) 'PARAMESH ERROR !'
           write(*,*) 'pe ',mype,' blk ',isg,' getting ialive from ',
     .       'a bad address ',new_loct(1,isg),new_loct(2,isg)
           call amr_abort
      endif
#endif

        call shmem_integer_get(idataout(1),ialive(new_loct(1,isg)),1,
     .       new_loct(2,isg))
        jalive = idataout(1)

#ifdef DIAG
               write(*,*) 'pe ',mype,' blk ',isg,' querying pe/blk ',
     .       new_loct(2,isg),new_loct(1,isg),' alive? ',jalive
#endif
        if(jalive.eq.-1) then
          source_id(1) = isg
          source_id(2) = mype
          call shmem_integer_put(src_loc(1,new_loct(1,isg)),
     .       source_id,2,new_loct(2,isg))
          kalive(isg) = -1
!         ialive(isg) = -1
#ifdef DIAG
            write(*,*) 'pe ',mype,' blk ',isg,' marking pe/blk ',
     .       new_loct(2,isg),new_loct(1,isg),' to fetch data.'
#endif
        endif

        endif
      endif

      enddo
      call shmem_barrier_all()
      call shmem_udcflush()

! If a block to be moved has an empty destination into which
! it can be moved, then reset it^s ialive marker to record the
! fact that its current location will be avilable for storage use
! on the next cycle.
!      do isg=1,maxblocks_alloc
      do isg=1,maxblocks
        if(kalive(isg).eq.-1) ialive(isg) = -1
      enddo

      call shmem_barrier_all()

! Now loop over the entire list of blocks and for any block which
! has been marked to fetch data (ie has src_loc(1,*) not equal 
! to -1) do a data fetch.
      do isg = 1,maxblocks

! Does this block need to fetch data?
      if(src_loc(1,isg).ne.-1) then
#ifdef DIAG
       write(*,*) 'pe ',mype,' blk ',isg,' fetching from pe/blk ',
     .       src_loc(2,isg),src_loc(1,isg)
#endif
        if(nvar.gt.0) call shmem_real_get(unk(1,1,1,1,isg),
     .       unk(1,1,1,1,src_loc(1,isg)),
     .       len_block,src_loc(2,isg))
        if(nfacevar.gt.0) then
          call shmem_real_get(facevarx(1,1,1,1,isg),
     .        facevarx(1,1,1,1,src_loc(1,isg)),
     .        len_blockfx*nbndvar,src_loc(2,isg))
          if(ndim.ge.2) call shmem_real_get(facevary(1,1,1,1,isg),
     .        facevary(1,1,1,1,src_loc(1,isg)),
     .        len_blockfy*nbndvar,src_loc(2,isg))
          if(ndim.eq.3) call shmem_real_get(
     .        facevarz(1,1,1,1,isg),
     .        facevarz(1,1,1,1,src_loc(1,isg)),
     .        len_blockfz*nbndvar,src_loc(2,isg))
        endif
        if(nvaredge.gt.0) then
          call shmem_real_get(unk_e_x(1,1,1,1,isg),
     .        unk_e_x(1,1,1,1,src_loc(1,isg)),
     .        len_blockex*nbndvare,src_loc(2,isg))
          call shmem_real_get(unk_e_y(1,1,1,1,isg),
     .        unk_e_y(1,1,1,1,src_loc(1,isg)),
     .        len_blockey*nbndvare,src_loc(2,isg))
          call shmem_real_get(unk_e_z(1,1,1,1,isg),
     .        unk_e_z(1,1,1,1,src_loc(1,isg)),
     .        len_blockez*nbndvare,src_loc(2,isg))
        endif
        if(nvarcorn.gt.0) call shmem_real_get(unk_n(1,1,1,1,isg),
     .        unk_n(1,1,1,1,src_loc(1,isg)),
     .        len_blockn*nvarcorn,src_loc(2,isg))

        new_loct(1,isg) = isg
        new_loct(2,isg) = mype
        ialive(isg) = 0
        pesource(1,isg) = src_loc(1,isg)
        pesource(2,isg) = src_loc(2,isg)
        src_loc(1,isg) = -1
      endif
      enddo


! count any remaining blocks that need to be sent
      icount = 0
!      do isg = 1,maxblocks_alloc
      do isg = 1,maxblocks

! Does this block need to be sent?
      if(ialive(isg).eq.0) then
        if( (new_loct(1,isg).ne.isg).or.
     .       (new_loct(2,isg).ne.mype) ) then

          icount = icount + 1
        endif
      endif

      enddo

      call shmem_barrier_all()

! if there are any rebuild the storage index, recompute nstore,
! then cycle back.
      if(icount.gt.0) then
        istore = 1
!        do i= istart_store,maxblocks_alloc
        do i= istart_store,maxblocks
          if(ialive(i).eq.-1) then
            indx(istore) = i
            istore = istore + 1
          endif
        enddo
        nstore = istore - 1
      endif

! because there is a shmem_barrier_all above, if any processor repeats the 
! cycle, all the processors must. those that no longer have any work to do
! will jump straight to the shmem_barrier_all, and when the shmem_barrier_all 
! clears will jump to here.
30    continue
      call shmem_barrier_all()
      call comm_int_max_to_all(igcount,icount)
      call shmem_barrier_all()
      if(mype.eq.0) write(*,*) 'redist: proc ',mype,' igcount = ',
     .       igcount
      call shmem_barrier_all()

      if(igcount.gt.0) goto 10

      do i=1,maxblocks
        new_loc(1,i) = new_loct(1,i)
        new_loc(2,i) = new_loct(2,i)
      enddo
      call shmem_barrier_all()

      return
      end subroutine amr_redist_blk
