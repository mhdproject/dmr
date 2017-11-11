!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"


      subroutine amr_restrict(mype,iopt,iempty,filling_guardcells)




!------------------------------------------------------------------------
!
! This routine does the data averaging required when a child block
! passes data back to its parent. The parent receives interior data
! only, not guard cell data. 
! This routine calls amr_1blk_restrict, which in turn calls the
! following user provided routines
!    amr_restrict_unk_fun
!    amr_restrict_unk_fun_recip
!    amr_restrict_fc_fun
!    amr_restrict_work_fun
!    amr_restrict_work_fun_recip.
! which defines the pattern of restriction.
! Default versions of these routines are provided. These can be easily
! modified by the user to suit their needs.
!
!
! Arguments :
!       mype    integer         Current processor number
!       iopt    integer         Switch to select which datastructures
!                                are updated. If iopt=1 amr_1blk_restrict
!                                acts on UNK and/or FACEVAR[X][Y][Z]
!                                and/or UNK_E_[X][Y][Z] and/or UNK_N.
!                                If iopt=2 only WORK is updated.
!       iempty  integer         NO FUNCTION !!!! Kept for compatibility with
!                                earlier versions
!
! Written :     Peter MacNeice          July 1997
!------------------------------------------------------------------------

      use paramesh_dimensions
      use physicaldata
      use tree

      use paramesh_mpi_interfaces, only : 
     &                  mpi_amr_1blk_restrict

      implicit none
#include "amr_shmem.fh"

#ifdef TIMINGS
#include "timer.fh"
#endif

      integer, intent(in)    :: mype,iopt,iempty
      logical, optional, intent(in)  :: filling_guardcells

      logical :: lcc,lfc,lec,lnc,lfulltree,fillingt

!------------------------------------------------------------------------

      lrestrict_in_progress = .true.

      if (present(filling_guardcells)) then
         fillingt = filling_guardcells
      else
         fillingt = .false.
      end if

      if( (iopt.gt.1) .and. (mod(nxb,2).ne.0) 
     .     .and. (mype.eq.0) )  write(*,*) 
     .      'Restriction Warning ! Applying restriction to ',
     .      'WORK may lead to errors at external boundaries! '


#ifdef TIMINGS
      call shmem_barrier_all()
      itimer1 = irtc()
#endif


      lcc = .false.
      lfc = .false.
      lec = .false.
      lnc = .false.
      if(nvar.gt.0.or.iopt.ge.2) lcc = .true.
      if(nfacevar.gt.0.and.iopt.eq.1) lfc = .true.
      if(nvaredge.gt.0.and.iopt.eq.1) lec = .true.
      if(nvarcorn.gt.0.and.iopt.eq.1) lnc = .true.

!
! Perform restriction
      lfulltree = .false.
      call mpi_amr_1blk_restrict(mype,iopt,lcc,lfc,lec,lnc,
     .                           lfulltree,fillingt)

#ifdef TIMINGS
      itimer2 = irtc()
      call shmem_barrier_all()
      irtc_restrict = itimer2-itimer1+irtc_restrict
#endif

      lrestrict_in_progress = .false.

      return
      end subroutine amr_restrict
