#include "paramesh_preprocessor.fh"
        subroutine field_at_cell_corners

!
! This routine computes the magnetic field at cell corners.
!
! Written:       Peter MacNeice
!                March 2000
!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! include file to define physical qualities of the model and mesh
      use paramesh_dimensions
      use physicaldata
      use tree
      use workspace

      use paramesh_interfaces
#ifdef MPI_USED
      use paramesh_mpi_interfaces
#endif /* MPI_USED */
      use strongfields

#ifdef MPI_USED
      include "mpif.h"
      integer :: errcode
#endif /* MPI_USED */

! include file required for shmem library.
#include "amr_shmem.fh"
#include "shmem_reduce.fh"


#include "b_corner.common"

      integer :: nprocs,mype
      integer :: shmem_my_pe,shmem_n_pes

      logical :: lcc,lfc,ldiag,l_srl_only,lec,lnc, lfulltree
      logical :: lguard,lprolong,lflux,ledge,lrestrict
      integer :: tag_offset


      integer, parameter :: ng0 = nguard*npgs


      real :: xmin,ymin
      save mype

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        mype   = shmem_my_pe()
        nprocs = shmem_n_pes()


!--------------------
! Perform a global guardcell filling, if permanent guardcell storage
! is allocated - otherwise make the necessary copy of the solution
! which the 1blk_guardcell routines use.
 
#ifdef NO_PERMANENT_GUARDCELLS
! Store a copy of the current solution in GT_UNK, GT_FACEVAR*
        call amr_1blk_copy_soln(-1)
#else
        iopt = 1
        nlayers = nguard
#ifdef MPI_USED
!        call mpi_amr_guardcell(mype,iopt,nlayers,nprocs)
        call amr_guardcell(mype,iopt,nlayers)
#else
        call     amr_guardcell(mype,iopt,nlayers)
#endif /*MPI_USED*/
#endif


!--------------------
      iopt       = 1
      lcc        = .true.
      lfc        = .true.
#ifdef MPI_USED
      lguard    = .true.
      lprolong  = .false.
      lflux     = .false.
      ledge     = .false.
      lrestrict = .false.
      lec        =.false.
      lnc        =.false.
      lec        =.false.
      lnc        =.false.
      lfulltree=.false.
      tag_offset= 100
      call mpi_amr_comm_setup(mype,nprocs,
     .                        lguard,lprolong,lflux,ledge,lrestrict,
     .                        lfulltree,
     .                        iopt,lcc,lfc,lec,lnc,tag_offset)
#endif /*MPI_USED*/


        if(lnblocks.gt.0) then
        do lb = 1,lnblocks
        if(nodetype(lb).eq.1) then


          idest = 1
          lec   = .false.
          lnc   = .false.
c          lcc = .true.
c          lfc = .true.
c          iopt = 1
#ifdef NO_PERMANENT_GUARDCELLS
! Copy data from current block into working block and fill its guardcells
          nlayers = nguard
          l_srl_only = .false.
          icoord = 0
          ldiag = .true.
!#ifdef MPI_USED
!      call mpi_amr_1blk_guardcell(mype,iopt,nlayers,lb,mype,lcc,lfc,
!     .                              l_srl_only,icoord,ldiag,
!     .                              nprocs)
!#else  /* MPI_USED */
      call amr_1blk_guardcell(mype,iopt,nlayers,lb,mype,lcc,lfc,
     .                        lec,lnc,
     .                        l_srl_only,icoord,ldiag)
!#endif /* MPI_USED */

#else

          call amr_perm_to_1blk( lcc,lfc,lec,lnc,lb,mype,iopt,idest)
#endif


! IMPORTANT : This call to strongfields_wb must come AFTER
! the call to amr_1blk_guardcell. If the order is reversed, 
! because amr_1blk_guardcell uses
! the working block arrays to set strongfield in the parent block
! of block lb, if this has a zero left BC for guardcells, the
! parents left BC inadvertantly overwrite the
! children's own BCs.
!#ifdef STRONGFIELD
! compute the strong field for this block
          call strongfields_wb(mype,lb)
!#endif /* STRONGFIELD */



          do k=1,nzb+k3d
          do j=1,nyb+1
          do i=1,nxb+1

          bbx(i,j,k,lb) = 
     .    + .25*(
     .      bxf_sf(1,i+nguard,j+nguard,k+nguard*k3d)
     .    + bxf_sf(1,i+nguard,j+nguard-1,k+nguard*k3d)
     .    + bxf_sf(1,i+nguard,j+nguard,k+(nguard-1)*k3d)
     .    + bxf_sf(1,i+nguard,j+nguard-1,k+(nguard-1)*k3d)
     .                   )
     .    + .25*(
     .      facevarx1(1,i+nguard,j+nguard,k+nguard*k3d,1)
     .    + facevarx1(1,i+nguard,j+nguard-1,k+nguard*k3d,1)
     .    + facevarx1(1,i+nguard,j+nguard,k+(nguard-1)*k3d,1)
     .    + facevarx1(1,i+nguard,j+nguard-1,k+(nguard-1)*k3d,1)
     .                   )
          bby(i,j,k,lb) =
     .    + .25*(
     .      byf_sf(1,i+nguard,j+nguard,k+nguard*k3d)
     .    + byf_sf(1,i+nguard-1,j+nguard,k+nguard*k3d)
     .    + byf_sf(1,i+nguard,j+nguard,k+(nguard-1)*k3d)
     .    + byf_sf(1,i+nguard-1,j+nguard,k+(nguard-1)*k3d)
     .                   )
     .    + .25*(
     .      facevary1(1,i+nguard,j+nguard,k+nguard*k3d,1)
     .    + facevary1(1,i+nguard-1,j+nguard,k+nguard*k3d,1)
     .    + facevary1(1,i+nguard,j+nguard,k+(nguard-1)*k3d,1)
     .    + facevary1(1,i+nguard-1,j+nguard,k+(nguard-1)*k3d,1)
     .                   )
          bbz(i,j,k,lb) =
     .    + .25*(
     .      bzf_sf(1,i+nguard,j+nguard,k+nguard*k3d)
     .    + bzf_sf(1,i+nguard-1,j+nguard,k+nguard*k3d)
     .    + bzf_sf(1,i+nguard,j+nguard-1,k+nguard*k3d)
     .    + bzf_sf(1,i+nguard-1,j+nguard-1,k+nguard*k3d)
     .                   )
     .    + .25*(
     .      facevarz1(1,i+nguard,j+nguard,k+nguard*k3d,1)
     .    + facevarz1(1,i+nguard-1,j+nguard,k+nguard*k3d,1)
     .    + facevarz1(1,i+nguard,j+nguard-1,k+nguard*k3d,1)
     .    + facevarz1(1,i+nguard-1,j+nguard-1,k+nguard*k3d,1)
     .                   )

          enddo
          enddo
          enddo


        endif
        enddo
        endif

#ifdef MPI_USED
        call mpi_barrier (MPI_COMM_WORLD, errcode)
#else
        call shmem_barrier_all()
#endif /*MPI_USED*/
      return
      end
