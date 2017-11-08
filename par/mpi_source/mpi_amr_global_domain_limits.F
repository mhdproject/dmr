!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"

         subroutine mpi_amr_global_domain_limits


!
! This routine computes the coordinate ranges for the grid
! as it currently exists.
!

         use paramesh_dimensions
         use physicaldata
         use tree

         implicit none

         include 'mpif.h'

         real :: xmin,ymin,zmin,xmax,ymax,zmax
         real :: xmin1,ymin1,zmin1,xmax1,ymax1,zmax1
         integer :: ierr


! Find the coordinate ranges
         xmin1 = minval(bnd_box(1,1,1:lnblocks))
         ymin1 = minval(bnd_box(1,2,1:lnblocks))
         zmin1 = minval(bnd_box(1,3,1:lnblocks))
         xmax1 = maxval(bnd_box(2,1,1:lnblocks))
         ymax1 = maxval(bnd_box(2,2,1:lnblocks))
         zmax1 = maxval(bnd_box(2,3,1:lnblocks))
         xmin = min(1.e30, xmin1)
         ymin = min(1.e30, ymin1)
         zmin = min(1.e30, zmin1)
         xmax = max(-1.e30, xmax1)
         ymax = max(-1.e30, ymax1)
         zmax = max(-1.e30, zmax1)
      call MPI_ALLREDUCE (xmin,grid_xmin,1,
     &                    MPI_REAL,
     &                    MPI_MIN,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE (ymin,grid_ymin,1,
     &                    MPI_REAL,
     &                    MPI_MIN,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE (zmin,grid_zmin,1,
     &                    MPI_REAL,
     &                    MPI_MIN,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE (xmax,grid_xmax,1,
     &                    MPI_REAL,
     &                    MPI_MAX,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE (ymax,grid_ymax,1,
     &                    MPI_REAL,
     &                    MPI_MAX,MPI_COMM_WORLD,ierr)
      call MPI_ALLREDUCE (zmax,grid_zmax,1,
     &                    MPI_REAL,
     &                    MPI_MAX,MPI_COMM_WORLD,ierr)

      return
      end subroutine mpi_amr_global_domain_limits
