#include "paramesh_preprocessor.fh"
#include "athena_preprocessor.fh"
      subroutine initialstate(mype,xmin,ymin,zmin,xmax,ymax,zmax)
      implicit none
      integer :: mype
      real    :: xmin,ymin,zmin,xmax,ymax,zmax
#ifdef MPI_USED
      include "mpif.h"
      integer :: errcode
#endif /* MPI_USED */

#ifdef JET
         write(*,*)'jet:entering amr_initial_jet I'
	call amr_initial_jet(mype,xmin,ymin,zmin,xmax,ymax,zmax)
         write(*,*)'jet:leaving amr_initial_jet I'
#endif /*JET*/

#ifdef EXPLOSION
      call amr_initial_bomb(mype,xmin,ymin,zmin,xmax,ymax,zmax)
#endif
#ifdef PULSAR
      call amr_initial_jet(mype,xmin,ymin,zmin,xmax,ymax,zmax)
#endif
#ifdef MAGNETOSPHERE
c      write(*,*)'entering amr_initial_magnetosphere' 
      call amr_initial_magnetosphere
     &       (mype,xmin,ymin,zmin,xmax,ymax,zmax)
#endif
#ifdef DIVERGENCE_TEST
      call amr_initial_divb(mype,xmin,ymin,zmin,xmax,ymax,zmax)
#endif

#ifdef MPI_USED
        call mpi_barrier (MPI_COMM_WORLD, errcode)
#else
        call shmem_barrier_all()
#endif /*MPI_USED*/

      return
      end subroutine initialstate
