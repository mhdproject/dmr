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

      subroutine mpi_amr_morton_limits(mype)



!------------------------------------------------------------------------
!
! This routine calculates the range of (morton no., ref. level) pairs on
! all processors and communicates that information so that each
! processor knows the ranges on every other processor.
!
!
! Written :     Peter MacNeice          June 2000
!------------------------------------------------------------------------
!
! Arguments:
!      mype            rank of local processor
!
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use mpi_morton
      use paramesh_mpi_interfaces, only : morton_number

      implicit none

      include 'mpif.h'

      integer, intent(in)    ::  mype

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! local variables

      real    :: dx,dy,dz, x0,y0,z0
      real    :: xmin, ymin, zmin
      real    :: xmax, ymax, zmax

      integer :: lb, ierror
      integer :: morton(6)
#ifdef DEBUG_MORTON
      integer :: inxt_bit,nbitshft,nbits
#endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!
!
! This routine assumes that the grid blocks are ordered by morton
! number and that any blocks with different refinement levels but
! the same morton number are ordered from coarse to fine.


! mark morton data out of date.
      morton_limits_set = .false.

!--------------------------------------------------
!
! Compute the list of (morton number, refinement level) pairs for all 
! local blocks.

!
! Compute xmin,ymin,zmin,xmax,ymax,zmax or get them from storage
      xmin = grid_xmin
      ymin = grid_ymin
      zmin = grid_zmin
      xmax = grid_xmax
      ymax = grid_ymax
      zmax = grid_zmax

      mortonbnd = -1

      do lb=1,lnblocks

        dx = bsize(1,lb)
        dy = bsize(2,lb)
        dz = bsize(3,lb)

!#ifndef ADVANCE_ALL_LEVELS
!        if(nodetype(lb).le.2) then
!#endif

! set the block center
          x0 = coord(1,lb) - xmin
          y0 = coord(2,lb) - ymin
          z0 = coord(3,lb) - zmin


! compute the morton number for block lb 
          call morton_number (x0,y0,z0,bsize(:,lb),ndim,
     .                        lrefine_max,lrefine(lb),morton)


! For each block the first entry in mortonbnd is the block^s morton
! number, the second entry is it^s refinement level, and the third
! is an integer flag indicating whether or not the block is a new child.

          mortonbnd(:,1,lb) = morton(:)
          mortonbnd(:,2,lb) = lrefine(lb)
          mortonbnd(:,3,lb) = nodetype(lb)

!#ifndef ADVANCE_ALL_LEVELS
!        endif
!#endif
      enddo


!--------------------------------------------------
!
! Step 2.
! Store the limits of the lists constructed in step 1 and share them
! with all processors.

!
! Store the local limits of the (morton number, refinement level) list.
      morton_limits(:,1:2,1,mype+1) = mortonbnd(:,1:2,1)
      if(lnblocks.gt.0) then
        morton_limits(:,1:2,2,mype+1) = mortonbnd(:,1:2,lnblocks)
      else
        morton_limits(:,1:2,2,mype+1) = -100
      endif

     
#ifdef DEBUG 
       write(*,*) 'proc ',mype,' mortonbnd ',
     .          mortonbnd(:,1:lnblocks),' lnblocks ',lnblocks,
     .       ' mortonbnd ',mortonbnd(:,lnblocks+1:)
       write(*,*) 'proc ',mype,' morton_limits ',
     .                        morton_limits(1:2,:,mype+1)
#endif /* DEBUG  */
!--------------------------------------------------

      call MPI_Allgather(morton_limits(1,1,1,mype+1),  6*4,MPI_INTEGER,
     .                   morton_limits,6*4,MPI_INTEGER,
     .                   MPI_COMM_WORLD,ierror)

!--------------------------------------------------

! Mark morton data up to date
       morton_limits_set = .true.


      return
      end subroutine mpi_amr_morton_limits
