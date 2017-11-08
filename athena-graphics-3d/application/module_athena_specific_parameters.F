#include "paramesh_preprocessor.fh"
      module athena_specific_parameters

      use paramesh_dimensions

      integer, parameter :: ne      = N_E
      integer, parameter :: nem1    = N_E-1
      integer, parameter :: nf      = ne+3 ! number of fluxes
      integer, parameter :: iu_bnd1p=iu_bnd1+k1d
      integer, parameter :: ju_bnd1p=ju_bnd1+k2d
      integer, parameter :: ku_bnd1p=ku_bnd1+k3d
      integer, parameter :: nmax    =max(iu_bnd1,ju_bnd1,ku_bnd1)
#ifdef MUSCL  
      integer, parameter :: norder = 2
#endif
#ifdef PLM  
      integer, parameter :: norder = 3
#endif
#ifdef PPM  
      integer, parameter :: norder = 4
#endif
      integer, parameter :: nx     = iu_bnd1
      integer, parameter :: ny     = ju_bnd1
      integer, parameter :: nz     = ku_bnd1
      integer, parameter :: nlx    = 1+norder*k1d
      integer, parameter :: nhx    = nx-(norder-1)*k1d

      integer, parameter :: nlxp1  = 1 +(norder+1)*k1d
      integer, parameter :: nlxp2  = 1 +(norder+2)*k1d
      integer, parameter :: nhxm1  = nx-(norder  )*k1d
      integer, parameter :: nhxm2  = nx-(norder+1)*k1d
      integer, parameter :: nly    = 1+norder*k2d
      integer, parameter :: nhy    = ny-(norder-1)*k2d
      integer, parameter :: nlyp1  = 1 +(norder+1)*k2d
      integer, parameter :: nlyp2  = 1 +(norder+2)*k2d
      integer, parameter :: nhym1  = ny-(norder  )*k2d
      integer, parameter :: nhym2  = ny-(norder+1)*k2d
      integer, parameter :: nlz    = 1+norder*k3d
      integer, parameter :: nhz    = nz-(norder-1)*k3d
      integer, parameter :: nlzp1  = 1 +(norder+1)*k3d
      integer, parameter :: nlzp2  = 1 +(norder+2)*k3d
      integer, parameter :: nhzm1  = nz-(norder  )*k3d
      integer, parameter :: nhzm2  = nz-(norder+1)*k3d

      end
