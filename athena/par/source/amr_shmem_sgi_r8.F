!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"

! this subroutine was added to give this file some fortran statements
! since otherwise it is ignored by the T3E is never compiled
      subroutine For_cray
      end subroutine For_cray

#if defined(SGI_SHMEM) | defined(SHMEM_DUMMYLIB)


      subroutine shmem_real_get (target,source,nwords,pe)
      real target,source
      integer nwords,pe
#include "amr_shmem.fh"

#ifdef REAL8
      call shmem_get64(target,source,nwords,pe)
#else
      call shmem_get32(target,source,nwords,pe)
#endif

      return
      end

      subroutine shmem_real_put (target,source,nwords,pe)
      real target,source
      integer nwords,pe
#include "amr_shmem.fh"

#ifdef REAL8
      call shmem_put64(target,source,nwords,pe)
#else
      call shmem_put32(target,source,nwords,pe)
#endif

      return
      end

#endif /* SGI_SHMEM */
