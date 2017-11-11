#include "athena_preprocessor.fh"
!#define BRIO_WU
      module physcons
      implicit none
c      save
      real pie,pie4,pie4i,pie8,pie8i
      real fmhd
      real gamma,gammai,gammam1,gammam1i,gamp
#ifdef BRIO_WU
      parameter(
     $     pie     = 3.14159265358979323846,
     $     pie4    = 1.0               ,
     $     pie4i   = 1.0              ,
     $     pie8    = 2.0               ,
     $     pie8i   = 0.5              
     $                                   )
      parameter (gamma = 2.0)
      parameter(
     $     gammai  = 1.0/gamma,
     $     gammam1 = gamma-1.0,
     $     gammam1i= 1.0/gammam1,
     $     gamp    = gamma*gammam1 )
#else
      parameter(
     $     pie     = 3.14159265358979323846,
     $     pie4    = pie*4.0               ,
     $     pie4i   = 1.0/pie4              ,
     $     pie8    = pie*8.0               ,
     $     pie8i   = 1.0/pie8              
     $                                   )
      parameter (gamma = 1.4/1.0)
      parameter(
     $     gammai  = 1.0/gamma,
     $     gammam1 = gamma-1.0,
     $     gammam1i= 1.0/gammam1 ,
     $     gamp    = gamma*gammam1 )
#endif
      real, parameter,public :: rfloor=0.4
      real, parameter,public :: rearth=6.371
      real, parameter,public :: clight=300.0 ! in megameters/2
      real, parameter,public :: rwall =1.0e+15
#ifdef BORIS_CORRECTION
      real, parameter,public :: cboris= 5.0      ! in megameters
      real, parameter,public :: afact  = pie4i/(cboris*cboris)
#else
      real, parameter,public :: cboris= clight   ! in megameters
      real, parameter,public :: afact  = 0.0
#endif

      real,public :: rmin,pmin,smin

      real, parameter,public :: pulsar_radius=0.
      real, parameter,public :: pulsar_beta=0.6


      end module physcons
