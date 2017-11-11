#include "athena_preprocessor.fh"
!     subroutine initialmesh(mype,xmin,xmax,ymin,ymax,zmin,zmax)
      subroutine initialmesh(mype,xmin,ymin,zmin,xmax,ymax,zmax)
      integer :: mype
      real    :: xmin,ymin,zmin,xmax,ymax,zmax

#ifdef MAGNETOSPHERE 
      call amr_initial_grid_magnetosphere(
     &         mype,xmin,ymin,zmin,xmax,ymax,zmax)
#endif
#ifdef EXPLOSION
!-----------------------------------------------------------------------
      call amr_initial_grid_bomb(
     &         mype,xmin,ymin,zmin,xmax,ymax,zmax)
#endif
#ifdef PULSAR
!-----------------------------------------------------------------------
      call amr_initial_grid_pulsar(
     &         mype,xmin,ymin,zmin,xmax,ymax,zmax)
#endif
#ifdef DIVERGENCE_TEST
!-----------------------------------------------------------------------
      call amr_initial_grid_divb(
     &         mype,xmin,ymin,zmin,xmax,ymax,zmax)
#endif
#ifdef JET
!-----------------------------------------------------------------------
      call amr_initial_grid_jet(
     &         mype,xmin,ymin,zmin,xmax,ymax,zmax)
#endif
#ifdef CIRCULATION
!-----------------------------------------------------------------------
      call amr_initial_grid_circulation(
     &         mype,xmin,ymin,zmin,xmax,ymax,zmax)
#endif
      return
      end

