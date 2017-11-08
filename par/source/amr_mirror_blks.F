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

      subroutine amr_mirror_blks(mype,nprocs,mirror_blks)



!
! DESIGN ISSUES :
!  
!  Some neighbors of parents are requested unnecessarily because we
!  cannot verify that the corresponding child neighbors actually
!  exist until after we have received morton lists.
!  Can we improve list of requested blocks by identifying neighbors
!  of parents which may not be needed?
!
!
!------------------------------------------------------------------------
!
! This routine calculates the morton number for each block on mype.
! It stores the result along with the refinement level of each block into
! the array mortonbnd, and distributes this array among all processors.
!
!
! Written :     Peter MacNeice         December 2002
!------------------------------------------------------------------------
!
! Arguments:
!      mype           rank of local processor
!
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
#ifdef MPI_USED
      use mpi_morton
#endif /* MPI_USED */

      implicit none

#ifdef MPI_USED
      include 'mpif.h'
#endif /* MPI_USED */


      integer, intent(in)    ::  mype,nprocs
      integer, intent(out)   ::  mirror_blks(2,maxblocks)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local variables

      integer      :: lb,mb,iproc


      integer,save :: clnblocks

      real         :: coordm(3)
      real,save    :: ccoord(3,maxblocks)
      real,save    :: ccoordm(3,maxblocks)
      logical      :: ldone(maxblocks),lcontinue_search


!--------------------------------------------------
!
! mark all mirror blocks unfound
      ldone = .false.

      lcontinue_search = .true.


! account for any off processor blocks

      iproc = 0

      do while( lcontinue_search .and. (iproc.lt.nprocs))
        lcontinue_search = .false.

        call shmem_real_get(clnblocks,lnblocks,
     .                      1,iproc)
        call shmem_real_get(ccoord(1,1),coord(1,1),
     .                      3*maxblocks,iproc)

        do lb=1,lnblocks
        if(nodetype(lb).eq.1) then
        if(.not.ldone(lb)) then

!--
! Mirror coordinate location

! equatorial symmetry
          coordm(1) = coord(1,lb)
          coordm(2) = grid_ymax - (coord(2,lb) - grid_ymin)
          coordm(3) = coord(3,lb)
!--

          do mb = 1,clnblocks
            if(   (abs(coordm(1)-ccoord(1,mb)).lt.1.e-6)
     .      .and. (abs(coordm(2)-ccoord(2,mb)).lt.1.e-6)
     .      .and. (abs(coordm(3)-ccoord(3,mb)).lt.1.e-6) )
     .                                      then
              mirror_blks(1,lb) = mb
              mirror_blks(2,lb) = iproc
              ldone(lb) = .true.
            endif
          enddo
        endif
        if(.not.ldone(lb)) lcontinue_search = .true.
        endif
        enddo

        iproc = iproc+1
      enddo


      return
      end subroutine amr_mirror_blks
