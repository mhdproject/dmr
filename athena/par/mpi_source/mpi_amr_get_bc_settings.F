#include "paramesh_preprocessor.fh"



      subroutine mpi_amr_get_bc_settings(blk,
     .                           ibxl,ibxr,ibyl,ibyr,ibzl,ibzr)

!----------------------------------------------------------------
!
! This routine searches through a list of block faces at external
! boundaries, so that boundary condition flags can be identified
! when needed, even for off-processor blocks. This needs the 
! routine mpi_amr_boundary_block_info to have been called whenever
! the grid is modified.
!
! Written : Peter MacNeice         November 2002
!
!----------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree

      implicit none

      integer,intent(in)    :: blk
      integer,intent(inout) :: ibxl,ibxr,ibyl,ibyr,ibzl,ibzr

      logical :: lfound
      integer :: ib,istart,iend

!----------------------------------------------------------------
#ifdef MPI_USED


      if (bc_block_neighs_status.ne.100) then
        write(*,*) 
     .   'Error : You must call mpi_amr_boundary_block_info!'
        call amr_abort()
      endif

! locate first entry for block (blk,pe) in bc_block_neighs list
      ib = 1
      lfound = .false.
      do
        if( (ib.ge.bc_block_neighs_length) .or. lfound ) exit

        if( (blk.eq.bc_block_neighs(1,ib)) .and.
     .      (blk.eq.bc_block_neighs(2,ib)) ) then
          istart = ib
          lfound = .true.
        endif
        ib = ib + 1
      enddo


! locate last entry for block (blk,pe) in bc_block_neighs list
      if(lfound) then
        ib = istart + 1
        do
          if(ib.ge.bc_block_neighs_length) exit
          if( (blk.ne.bc_block_neighs(1,ib)) .or.
     .        (blk.ne.bc_block_neighs(2,ib)) ) exit
          ib = ib + 1
        enddo
      iend = ib - 1
      endif

!
      if(lfound) then
        do ib = istart,iend
          if( (bc_block_neighs(3,ib).eq.1).and.
     .        (bc_block_neighs(4,ib).eq.1+k2d).and.
     .        (bc_block_neighs(5,ib).eq.1+k3d) ) 
     .                           ibxl = bc_block_neighs(6,ib)
          if( (bc_block_neighs(3,ib).eq.3).and.
     .        (bc_block_neighs(4,ib).eq.1+k2d).and.
     .        (bc_block_neighs(5,ib).eq.1+k3d) ) 
     .                           ibxr = bc_block_neighs(6,ib)
          if(ndim.ge.2) then
          if( (bc_block_neighs(3,ib).eq.2).and.
     .        (bc_block_neighs(4,ib).eq.1).and.
     .        (bc_block_neighs(5,ib).eq.1+k3d) ) 
     .                           ibyl = bc_block_neighs(6,ib)
          if( (bc_block_neighs(3,ib).eq.2).and.
     .        (bc_block_neighs(4,ib).eq.3).and.
     .        (bc_block_neighs(5,ib).eq.1+k3d) ) 
     .                           ibyr = bc_block_neighs(6,ib)
          endif
          if(ndim.eq.3) then
          if( (bc_block_neighs(3,ib).eq.2).and.
     .        (bc_block_neighs(4,ib).eq.2).and.
     .        (bc_block_neighs(5,ib).eq.1) ) 
     .                           ibzl = bc_block_neighs(6,ib)
          if( (bc_block_neighs(3,ib).eq.2).and.
     .        (bc_block_neighs(4,ib).eq.2).and.
     .        (bc_block_neighs(5,ib).eq.3) ) 
     .                           ibzr = bc_block_neighs(6,ib)
          endif
        enddo
      else
        ibxl = -1
        ibxr = -1
        ibyl = -1
        ibyr = -1
        ibzl = -1
        ibzr = -1
      endif

#endif /* MPI_USED */
      return
      end subroutine mpi_amr_get_bc_settings
