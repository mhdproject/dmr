#define GCELL_ON
#define GCELL_ON_FC
#define GCELL_ON_EC
#define GCELL_ON_NC
!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"


      subroutine amr_1blk_to_perm( lcc,lfc,lec,lnc,lb,iopt,idest)




!------------------------------------------------------------------------
!
! This routine copies data from the 1-block working arrays with guardcells
! to the permanent data arrays, which may or may not have permanent
! guardcells, depending on whether NO_PERMANENT_GUARDCELLS is defined 
! in physicaldata.fh.
!
!
! Arguments :
!      lcc          logical       copies cell centered data if true
!      lfc          logical       copies cell face-centered data if true
!      lec          logical       copies cell edge-centered data if true
!      lnc          logical       copies cell corner data if true
!      lb           integer       block into which data is to be copied
!      iopt         integer       data structure to be copied
!      idest        integer       sets value for last dimension index
!                                  in the 1-blk data arrays
!
!
! Written :     Peter MacNeice          February 1999
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree
      use timings
      use workspace

      implicit none

#include "amr_shmem.fh"

!------------------------------------

      integer, intent(in) :: lb,iopt,idest
      logical, intent(in) :: lcc,lfc,lec,lnc

      integer :: iopt0, ivar, ivar_next

#ifdef TIMING_MPI
      include 'mpif.h'
      double precision :: time1
#endif /* TIMING_MPI */
!------------------------------------


#ifdef TIMING_MPI
      time1 = mpi_wtime()
#endif /* TIMING_MPI */


! cell-centered data
       if(lcc) then

         if(iopt.eq.1) then

#ifdef NO_PERMANENT_GUARDCELLS
#ifdef GCELL_ON
           if(ngcell_on_cc < nvar) then
             do ivar=1,ngcell_on_cc
               ivar_next = gcell_on_cc_pointer(ivar)
               unk(ivar_next,il_bnd:iu_bnd,jl_bnd:ju_bnd,
     .                                     kl_bnd:ku_bnd,lb)
     .         = unk1(ivar_next,il_bnd+nguard:iu_bnd+nguard,
     .               jl_bnd+nguard*k2d:ju_bnd+nguard*k2d,
     .               kl_bnd+nguard*k3d:ku_bnd+nguard*k3d,idest)
             enddo
           else
             unk(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,lb)
     .         = unk1(:,il_bnd+nguard:iu_bnd+nguard,
     .               jl_bnd+nguard*k2d:ju_bnd+nguard*k2d,
     .               kl_bnd+nguard*k3d:ku_bnd+nguard*k3d,idest)
           endif
#else
           unk(:,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,lb)
     .       = unk1(:,il_bnd+nguard:iu_bnd+nguard,
     .               jl_bnd+nguard*k2d:ju_bnd+nguard*k2d,
     .               kl_bnd+nguard*k3d:ku_bnd+nguard*k3d,idest)
#endif /* GCELL_ON */
#else
#ifdef GCELL_ON
           if(ngcell_on_cc < nvar) then
             do ivar=1,ngcell_on_cc
               ivar_next = gcell_on_cc_pointer(ivar)
               unk(ivar_next,:,:,:,lb) = unk1(ivar_next,:,:,:,idest)
             enddo
           else
             unk(:,:,:,:,lb) = unk1(:,:,:,:,idest)
           endif
#else
           unk(:,:,:,:,lb) = unk1(:,:,:,:,idest)
#endif /* GCELL_ON */
#endif

         elseif(iopt.ge.2) then
           iopt0 = iopt-1

#ifdef NO_PERMANENT_GUARDCELLS
           work(ilw:iuw,jlw:juw,klw:kuw,lb,iopt0)
     .       = work1(ilw+nguard_work:iuw+nguard_work,
     .               jlw+nguard_work*k2d:juw+nguard_work*k2d,
     .               klw+nguard_work*k3d:kuw+nguard_work*k3d,
     .               idest)
#else
           work(:,:,:,lb,iopt0) = work1(:,:,:,idest)
#endif
         endif                           ! end of iopt if test
       endif                             ! end of lcc if test



! cell face-centered data
       if(lfc) then

! x-face
#ifdef NO_PERMANENT_GUARDCELLS
#ifdef GCELL_ON_FC
         if(ngcell_on_fc(1) < nfacevar) then
           do ivar=1,ngcell_on_fc(1)
             ivar_next = gcell_on_fc_pointer(1,ivar)
             facevarx(ivar_next,il_bnd:iu_bnd+1,
     .                jl_bnd:ju_bnd,kl_bnd:ku_bnd,lb)
     .         = facevarx1(ivar_next,il_bnd+nguard:iu_bnd+nguard+1,
     .                jl_bnd+nguard*k2d:ju_bnd+nguard*k2d,
     .                kl_bnd+nguard*k3d:ku_bnd+nguard*k3d,idest)
           enddo
         else
           facevarx(1:nfacevar,il_bnd:iu_bnd+1,
     .                       jl_bnd:ju_bnd,kl_bnd:ku_bnd,lb)
     .       = facevarx1(1:nfacevar,il_bnd+nguard:iu_bnd+nguard+1,
     .               jl_bnd+nguard*k2d:ju_bnd+nguard*k2d,
     .               kl_bnd+nguard*k3d:ku_bnd+nguard*k3d,idest)
         endif
#else
         facevarx(1:nfacevar,il_bnd:iu_bnd+1,
     .                       jl_bnd:ju_bnd,kl_bnd:ku_bnd,lb)
     .       = facevarx1(1:nfacevar,il_bnd+nguard:iu_bnd+nguard+1,
     .               jl_bnd+nguard*k2d:ju_bnd+nguard*k2d,
     .               kl_bnd+nguard*k3d:ku_bnd+nguard*k3d,idest)
#endif /* GCELL_ON_FC */
#else
#ifdef GCELL_ON_FC
         if(ngcell_on_fc(1) < nfacevar) then
           do ivar=1,ngcell_on_fc(1)
             ivar_next = gcell_on_fc_pointer(1,ivar)
             facevarx(ivar_next,:,:,:,lb) =
     .                        facevarx1(ivar_next,:,:,:,idest)
           enddo
         else
           facevarx(1:nfacevar,:,:,:,lb) =
     .                        facevarx1(1:nfacevar,:,:,:,idest)
         endif
#else
         facevarx(1:nfacevar,:,:,:,lb) =
     .                        facevarx1(1:nfacevar,:,:,:,idest)
#endif /* GCELL_ON_FC */
#endif

         if (ndim > 1) then
! y-face
#ifdef NO_PERMANENT_GUARDCELLS
#ifdef GCELL_ON_FC
         if(ngcell_on_fc(2) < nfacevar) then
           do ivar=1,ngcell_on_fc(2)
             ivar_next = gcell_on_fc_pointer(2,ivar)
             facevary(ivar_next,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,
     .                         kl_bnd:ku_bnd,lb)
     .       = facevary1(ivar_next,il_bnd+nguard:iu_bnd+nguard,
     .               jl_bnd+nguard*k2d:ju_bnd+nguard*k2d+k2d,
     .               kl_bnd+nguard*k3d:ku_bnd+nguard*k3d,idest)
           enddo
         else
           facevary(1:nfacevar,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,
     .                         kl_bnd:ku_bnd,lb)
     .       = facevary1(1:nfacevar,il_bnd+nguard:iu_bnd+nguard,
     .               jl_bnd+nguard*k2d:ju_bnd+nguard*k2d+k2d,
     .               kl_bnd+nguard*k3d:ku_bnd+nguard*k3d,idest)
         endif
#else
         facevary(1:nfacevar,il_bnd:iu_bnd,jl_bnd:ju_bnd+k2d,
     .                         kl_bnd:ku_bnd,lb)
     .       = facevary1(1:nfacevar,il_bnd+nguard:iu_bnd+nguard,
     .               jl_bnd+nguard*k2d:ju_bnd+nguard*k2d+k2d,
     .               kl_bnd+nguard*k3d:ku_bnd+nguard*k3d,idest)
#endif /* GCELL_ON_FC */
#else
#ifdef GCELL_ON_FC
         if(ngcell_on_fc(2) < nfacevar) then
           do ivar=1,ngcell_on_fc(2)
             ivar_next = gcell_on_fc_pointer(2,ivar)
             facevary(ivar_next,:,:,:,lb) =
     .                        facevary1(ivar_next,:,:,:,idest)
           enddo
         else
           facevary(1:nfacevar,:,:,:,lb) =
     .                        facevary1(1:nfacevar,:,:,:,idest)
         endif
#else
         facevary(1:nfacevar,:,:,:,lb) =
     .                        facevary1(1:nfacevar,:,:,:,idest)
#endif /* GCELL_ON_FC */
#endif

         end if

         if (ndim == 3) then
! z-face
#ifdef NO_PERMANENT_GUARDCELLS
#ifdef GCELL_ON_FC
         if(ngcell_on_fc(3) < nfacevar) then
           do ivar=1,ngcell_on_fc(3)
             ivar_next = gcell_on_fc_pointer(3,ivar)
             facevarz(ivar_next,il_bnd:iu_bnd,jl_bnd:ju_bnd,
     .                         kl_bnd:ku_bnd+k3d,lb)
     .       = facevarz1(ivar_next,il_bnd+nguard:iu_bnd+nguard,
     .               jl_bnd+nguard*k2d:ju_bnd+nguard*k2d,
     .               kl_bnd+nguard*k3d:ku_bnd+nguard*k3d+k3d,idest)
           enddo
         else
           facevarz(1:nfacevar,il_bnd:iu_bnd,jl_bnd:ju_bnd,
     .                         kl_bnd:ku_bnd+k3d,lb)
     .       = facevarz1(1:nfacevar,il_bnd+nguard:iu_bnd+nguard,
     .               jl_bnd+nguard*k2d:ju_bnd+nguard*k2d,
     .               kl_bnd+nguard*k3d:ku_bnd+nguard*k3d+k3d,idest)
         endif
#else
         facevarz(1:nfacevar,il_bnd:iu_bnd,jl_bnd:ju_bnd,
     .                         kl_bnd:ku_bnd+k3d,lb)
     .       = facevarz1(1:nfacevar,il_bnd+nguard:iu_bnd+nguard,
     .               jl_bnd+nguard*k2d:ju_bnd+nguard*k2d,
     .               kl_bnd+nguard*k3d:ku_bnd+nguard*k3d+k3d,idest)
#endif /* GCELL_ON_FC */
#else
#ifdef GCELL_ON_FC
         if(ngcell_on_fc(3) < nfacevar) then
           do ivar=1,ngcell_on_fc(3)
             ivar_next = gcell_on_fc_pointer(3,ivar)
             facevarz(ivar_next,:,:,:,lb) =
     .                        facevarz1(ivar_next,:,:,:,idest)
           enddo
         else
           facevarz(1:nfacevar,:,:,:,lb) =
     .                        facevarz1(1:nfacevar,:,:,:,idest)
         endif
#else
         facevarz(1:nfacevar,:,:,:,lb) =
     .                        facevarz1(1:nfacevar,:,:,:,idest)
#endif /* GCELL_ON_FC */
#endif

         end if

        endif                     ! end of lfc if test


! cell edge-centered data

       if (ndim > 1) then
       if(lec) then
! x-edge
#ifdef NO_PERMANENT_GUARDCELLS
#ifdef GCELL_ON_EC
         if(ngcell_on_ec(1) < nvaredge) then
           do ivar=1,ngcell_on_ec(1)
             ivar_next = gcell_on_ec_pointer(1,ivar)
             unk_e_x(ivar_next,il_bnd:iu_bnd,
     .               jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd+k3d,lb)
     .       = unk_e_x1(ivar_next,il_bnd+nguard:iu_bnd+nguard,
     .               jl_bnd+nguard*k2d:ju_bnd+(nguard+1)*k2d,
     .               kl_bnd+nguard*k3d:ku_bnd+(nguard+1)*k3d,idest)
           enddo
         else
           unk_e_x(1:nvaredge,il_bnd:iu_bnd,
     .               jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd+k3d,lb)
     .       = unk_e_x1(1:nvaredge,il_bnd+nguard:iu_bnd+nguard,
     .               jl_bnd+nguard*k2d:ju_bnd+(nguard+1)*k2d,
     .               kl_bnd+nguard*k3d:ku_bnd+(nguard+1)*k3d,idest)
         endif
#else
         unk_e_x(1:nvaredge,il_bnd:iu_bnd,
     .               jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd+k3d,lb)
     .       = unk_e_x1(1:nvaredge,il_bnd+nguard:iu_bnd+nguard,
     .               jl_bnd+nguard*k2d:ju_bnd+(nguard+1)*k2d,
     .               kl_bnd+nguard*k3d:ku_bnd+(nguard+1)*k3d,idest)
#endif /* GCELL_ON_EC */
#else
#ifdef GCELL_ON_EC
         if(ngcell_on_ec(1) < nvaredge) then
           do ivar=1,ngcell_on_ec(1)
             ivar_next = gcell_on_ec_pointer(1,ivar)
             unk_e_x(ivar_next,:,:,:,lb) =
     .                        unk_e_x1(ivar_next,:,:,:,idest)

           enddo
         else
           unk_e_x(1:nvaredge,:,:,:,lb) =
     .                        unk_e_x1(1:nvaredge,:,:,:,idest)
         endif
#else
         unk_e_x(1:nvaredge,:,:,:,lb) =
     .                        unk_e_x1(1:nvaredge,:,:,:,idest)
#endif /* GCELL_ON_EC */
#endif
! y-edge
#ifdef NO_PERMANENT_GUARDCELLS
#ifdef GCELL_ON_EC
         if(ngcell_on_ec(2) < nvaredge) then
           do ivar=1,ngcell_on_ec(2)
             ivar_next = gcell_on_ec_pointer(2,ivar)
             unk_e_y(ivar_next,il_bnd:iu_bnd+1,
     .               jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d,lb)
     .       = unk_e_y1(ivar_next,il_bnd+nguard:iu_bnd+1+nguard,
     .               jl_bnd+nguard*k2d:ju_bnd+nguard*k2d,
     .               kl_bnd+nguard*k3d:ku_bnd+(nguard+1)*k3d,idest)

           enddo
         else
           unk_e_y(1:nvaredge,il_bnd:iu_bnd+1,
     .               jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d,lb)
     .       = unk_e_y1(1:nvaredge,il_bnd+nguard:iu_bnd+1+nguard,
     .               jl_bnd+nguard*k2d:ju_bnd+nguard*k2d,
     .               kl_bnd+nguard*k3d:ku_bnd+(nguard+1)*k3d,idest)
         endif
#else
         unk_e_y(1:nvaredge,il_bnd:iu_bnd+1,
     .               jl_bnd:ju_bnd,kl_bnd:ku_bnd+k3d,lb)
     .       = unk_e_y1(1:nvaredge,il_bnd+nguard:iu_bnd+1+nguard,
     .               jl_bnd+nguard*k2d:ju_bnd+nguard*k2d,
     .               kl_bnd+nguard*k3d:ku_bnd+(nguard+1)*k3d,idest)
#endif /* GCELL_ON_EC */
#else
#ifdef GCELL_ON_EC
         if(ngcell_on_ec(2) < nedgevar) then
           do ivar=1,ngcell_on_ec(2)
             ivar_next = gcell_on_ec_pointer(2,ivar)
             unk_e_y(ivar_next,:,:,:,lb) =
     .                        unk_e_y1(ivar_next,:,:,:,idest)

           enddo
         else
           unk_e_y(1:nvaredge,:,:,:,lb) =
     .                        unk_e_y1(1:nvaredge,:,:,:,idest)
         endif
#else
         unk_e_y(1:nvaredge,:,:,:,lb) =
     .                        unk_e_y1(1:nvaredge,:,:,:,idest)
#endif /* GCELL_ON_EC */
#endif


         if (ndim == 3) then
! z-edge
#ifdef NO_PERMANENT_GUARDCELLS
#ifdef GCELL_ON_EC
         if(ngcell_on_ec(3) < nvaredge) then
           do ivar=1,ngcell_on_ec(3)
             ivar_next = gcell_on_ec_pointer(3,ivar)
             unk_e_z(ivar_next,il_bnd:iu_bnd+1,
     .               jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd,lb)
     .       = unk_e_z1(ivar_next,il_bnd+nguard:iu_bnd+1+nguard,
     .               jl_bnd+nguard*k2d:ju_bnd+(nguard+1)*k2d,
     .               kl_bnd+nguard*k3d:ku_bnd+nguard*k3d,idest)

           enddo
         else
           unk_e_z(1:nvaredge,il_bnd:iu_bnd+1,
     .               jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd,lb)
     .       = unk_e_z1(1:nvaredge,il_bnd+nguard:iu_bnd+1+nguard,
     .               jl_bnd+nguard*k2d:ju_bnd+(nguard+1)*k2d,
     .               kl_bnd+nguard*k3d:ku_bnd+nguard*k3d,idest)
         endif
#else
         unk_e_z(1:nvaredge,il_bnd:iu_bnd+1,
     .               jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd,lb)
     .       = unk_e_z1(1:nvaredge,il_bnd+nguard:iu_bnd+1+nguard,
     .               jl_bnd+nguard*k2d:ju_bnd+(nguard+1)*k2d,
     .               kl_bnd+nguard*k3d:ku_bnd+nguard*k3d,idest)
#endif /* GCELL_ON_EC */
#else
#ifdef GCELL_ON_EC
         if(ngcell_on_ec(3) < nvaredge) then
           do ivar=1,ngcell_on_ec(3)
             ivar_next = gcell_on_ec_pointer(3,ivar)
             unk_e_z(ivar_next,:,:,:,lb) =
     .                        unk_e_z1(ivar_next,:,:,:,idest)

           enddo
         else
           unk_e_z(1:nvaredge,:,:,:,lb) =
     .                        unk_e_z1(1:nvaredge,:,:,:,idest)
         endif
#else
         unk_e_z(1:nvaredge,:,:,:,lb) =
     .                        unk_e_z1(1:nvaredge,:,:,:,idest)
#endif /* GCELL_ON_EC */
#endif

         end if

        endif                     ! end of lec if test

        end if

! cell corner data
       if(lnc) then
#ifdef NO_PERMANENT_GUARDCELLS
#ifdef GCELL_ON_NC
         if(ngcell_on_nc < nvarcorn) then
           do ivar=1,ngcell_on_nc
             ivar_next = gcell_on_nc_pointer(ivar)
             unk_n(ivar_next,il_bnd:iu_bnd+1,
     .               jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd+k3d,lb)
     .       = unk_n1(ivar_next,il_bnd+nguard:iu_bnd+1+nguard,
     .               jl_bnd+nguard*k2d:ju_bnd+(nguard+1)*k2d,
     .               kl_bnd+nguard*k3d:ku_bnd+(nguard+1)*k3d,idest)
           enddo
         else
           unk_n(1:nvarcorn,il_bnd:iu_bnd+1,
     .               jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd+k3d,lb)
     .       = unk_n1(1:nvarcorn,il_bnd+nguard:iu_bnd+1+nguard,
     .               jl_bnd+nguard*k2d:ju_bnd+(nguard+1)*k2d,
     .               kl_bnd+nguard*k3d:ku_bnd+(nguard+1)*k3d,idest)
         endif
#else
         unk_n(1:nvarcorn,il_bnd:iu_bnd+1,
     .               jl_bnd:ju_bnd+k2d,kl_bnd:ku_bnd+k3d,lb)
     .       = unk_n1(1:nvarcorn,il_bnd+nguard:iu_bnd+1+nguard,
     .               jl_bnd+nguard*k2d:ju_bnd+(nguard+1)*k2d,
     .               kl_bnd+nguard*k3d:ku_bnd+(nguard+1)*k3d,idest)
#endif /* GCELL_ON_NC */
#else
#ifdef GCELL_ON_NC
         if(ngcell_on_nc < nvarcorn) then
           do ivar=1,ngcell_on_nc
             ivar_next = gcell_on_nc_pointer(ivar)
             unk_n(ivar_next,:,:,:,lb) =
     .                        unk_n1(ivar_next,:,:,:,idest)
           enddo
         else
           unk_n(1:nvarcorn,:,:,:,lb) =
     .                        unk_n1(1:nvarcorn,:,:,:,idest)
         endif
#else
         unk_n(1:nvarcorn,:,:,:,lb) =
     .                        unk_n1(1:nvarcorn,:,:,:,idest)
#endif /* GCELL_ON_NC */
#endif
        endif                     ! end of lnc if test


#ifdef TIMING_MPI
              timer_amr_1blk_to_perm(iopt) =  
     .                           timer_amr_1blk_to_perm(iopt)
     .                          + mpi_wtime() - time1
#endif /* TIMING_MPI */
      return
      end subroutine amr_1blk_to_perm
