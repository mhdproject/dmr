!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

#include "paramesh_preprocessor.fh"



	subroutine comm_start(MaxProcs,nprocs,mype)

        integer, intent(out) :: nprocs,mype
        integer, intent(in)  :: MaxProcs

#include "amr_shmem.fh"
	include 'shmem_reduce.fh'
#ifdef MPI_USED
        include 'mpif.h'
#endif

	integer :: shmem_n_pes,shmem_my_pe,ierror

       
#ifdef SGI_SHMEM
	call start_pes(0)
#endif
#ifdef MPI_USED
        call MPI_Init(ierror)
        if(ierror.ne.MPI_SUCCESS) write(*,*) 'MPI_init error'
#endif
	do i = 1,shmem_reduce_sync_size
		ipsync(i) = shmem_sync_value
	end do
	nprocs = shmem_n_pes()
	mype = shmem_my_pe()
	pestart = 0
	pestride = 0
	pesize = nprocs
        write(*,*) 'comm_start : shmem_reduce_sync_size ',
     .                           shmem_reduce_sync_size
        write(*,*) 'comm_start : shmem_reduce_min_wrkdata_size ',
     .                           shmem_reduce_min_wrkdata_size
        write(*,*) 'comm_start : pestart pestride pesize ',
     .        pestart, pestride, pesize

        return
	end subroutine comm_start


	subroutine comm_finish()
        integer :: ierror
#ifdef MPI_USED
        include 'mpif.h'

        write(*,*) 'calling MPI_Finalize'
        call MPI_Finalize(ierror)
#endif
	return
	end subroutine comm_finish


	subroutine comm_logical_or_to_all(target,source)
        logical, intent(in)  :: source
        logical, intent(out) :: target
	integer,save         :: target_l,source_l
	integer nred

#include "amr_shmem.fh"
	include 'shmem_reduce.fh'

	nred = 1
        source_l = 0
        if(source) source_l = 1
	call shmem_barrier_all()
#if defined(CRAY)
	call shmem_int8_sum_to_all(target_l,source_l,nred,pestart,
     -				pestride,pesize,ipwrk,ipsync)
#else
	call shmem_int4_sum_to_all(target_l,source_l,nred,pestart,
     -				pestride,pesize,ipwrk,ipsync)
#endif

        target = .false.
        if(target_l.gt.0) target = .true.

	call shmem_barrier_all()
	return
	end subroutine comm_logical_or_to_all

	subroutine comm_real_sum_to_all(target,source)
        real, intent(in)  :: source
        real, intent(out) :: target
	real,save :: target_l,source_l
	integer nred

#include "amr_shmem.fh"
	include 'shmem_reduce.fh'

	nred = 1
        source_l = source
	call shmem_barrier_all()
#if defined(CRAY) || defined(REAL8)
	call shmem_real8_sum_to_all(target_l,source_l,nred,pestart,
     -				pestride,pesize,pwrk,ipsync)
#else
	call shmem_real4_sum_to_all(target_l,source_l,nred,pestart,
     -				pestride,pesize,pwrk,ipsync)
#endif

        target = target_l

	call shmem_barrier_all()
	return
	end subroutine comm_real_sum_to_all




	subroutine comm_int_sum_to_all(target,source)
        integer, intent(in)  :: source
        integer, intent(out) :: target
	integer,save :: target_l,source_l
	integer nred

#include "amr_shmem.fh"
	include 'shmem_reduce.fh'

	nred = 1
        source_l = source
	call shmem_barrier_all()
#if defined(CRAY)
	call shmem_int8_sum_to_all(target_l,source_l,nred,pestart,
     -				pestride,pesize,ipwrk,ipsync)
#else
	call shmem_int4_sum_to_all(target_l,source_l,nred,pestart,
     -				pestride,pesize,ipwrk,ipsync)
#endif

        target = target_l

	call shmem_barrier_all()
	return
	end subroutine comm_int_sum_to_all







	subroutine comm_real_min_to_all(target,source)
        real, intent(in)  :: source
        real, intent(out) :: target
	real,save :: target_l,source_l
	integer nred

#include "amr_shmem.fh"
	include 'shmem_reduce.fh'

	nred = 1
        source_l = source
	call shmem_barrier_all()
#if defined(CRAY) || defined(REAL8)
	call shmem_real8_min_to_all(target_l,source_l,nred,pestart,
     -				pestride,pesize,pwrk,ipsync)
#else
	call shmem_real4_min_to_all(target_l,source_l,nred,pestart,
     -				pestride,pesize,pwrk,ipsync)
#endif


        target = target_l

	call shmem_barrier_all()
	return
	end subroutine comm_real_min_to_all




	subroutine comm_int_min_to_all(target,source)
        integer, intent(in)  :: source
        integer, intent(out) :: target
	integer,save :: target_l,source_l
	integer nred

#include "amr_shmem.fh"
	include 'shmem_reduce.fh'

	nred = 1
        source_l = source
	call shmem_barrier_all()
#if defined(CRAY) 
	call shmem_int8_min_to_all(target_l,source_l,nred,pestart,
     -				pestride,pesize,ipwrk,ipsync)
#else
	call shmem_int4_min_to_all(target_l,source_l,nred,pestart,
     -				pestride,pesize,ipwrk,ipsync)
#endif

        target = target_l

	call shmem_barrier_all()
	return
	end subroutine comm_int_min_to_all











	subroutine comm_real_max_to_all(target,source)
        real, intent(in)  :: source
        real, intent(out) :: target
	real,save :: target_l,source_l
	integer nred

#include "amr_shmem.fh"
	include 'shmem_reduce.fh'

	nred = 1
        source_l = source
	call shmem_barrier_all()
#if defined(CRAY) || defined(REAL8)
	call shmem_real8_max_to_all(target_l,source_l,nred,pestart,
     -				pestride,pesize,pwrk,ipsync)
#else
	call shmem_real4_max_to_all(target_l,source_l,nred,pestart,
     -				pestride,pesize,pwrk,ipsync)
#endif
        target = target_l

	call shmem_barrier_all()
	return
	end subroutine comm_real_max_to_all




	subroutine comm_int_max_to_all(target,source)
        integer, intent(in)  :: source
        integer, intent(out) :: target
	integer,save :: target_l,source_l
	integer nred

#include "amr_shmem.fh"
	include 'shmem_reduce.fh'

	nred = 1
        source_l = source
	call shmem_barrier_all()
#if defined(CRAY)
	call shmem_int8_max_to_all(target_l,source_l,nred,pestart,
     -				pestride,pesize,ipwrk,ipsync)
#else
	call shmem_int4_max_to_all(target_l,source_l,nred,pestart,
     -				pestride,pesize,ipwrk,ipsync)
#endif
        target = target_l

	call shmem_barrier_all()
	return
	end subroutine comm_int_max_to_all



        subroutine comm_int_to_all(target,peroot)
        integer, intent(in)  :: peroot
        integer, intent(inout) :: target
        integer,save :: target_l,source_l
        integer :: nred

#include "amr_shmem.fh"
        include 'shmem_reduce.fh'

        nred = 1
        source_l = target
        call shmem_barrier_all()
#ifdef SHMEM_DUMMYLIB
        target_l = source_l
#else
#if defined(CRAY)
        call shmem_broadcast64(target_l,source_l,nred,peroot,
     -                          pestart,pestride,pesize,ibsync)
#else
        call shmem_broadcast32(target_l,source_l,nred,peroot,
     -                          pestart,pestride,pesize,ibsync)
#endif
#endif
        target = target_l

        call shmem_barrier_all()
        return
        end subroutine comm_int_to_all

        subroutine comm_real_to_all(target,peroot)
        real, intent(inout) :: target
        integer, intent(in)  :: peroot
        real,save :: target_l,source_l
        integer :: nred

#include "amr_shmem.fh"
        include 'shmem_reduce.fh'

        nred = 1
        source_l = target
        call shmem_barrier_all()
#ifdef SHMEM_DUMMYLIB
        target_l = source_l
#else
#if defined(REAL8)
        call shmem_broadcast64(target_l,source_l,nred,peroot,
     -                          pestart,pestride,pesize,ibsync)
#else
        call shmem_broadcast32(target_l,source_l,nred,peroot,
     -                          pestart,pestride,pesize,ibsync)
#endif
#endif
        target = target_l

        call shmem_barrier_all()
        return
        end subroutine comm_real_to_all


