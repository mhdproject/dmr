!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

!
! This file contains a replacement library for all the shmem routines
! called by the package, which will enable the user to compile and
! run their code on a single processor machine which does not have
! the true shmem library or mpi. This may be useful for debugging
! purposes.
!




#include "paramesh_preprocessor.fh"

#ifndef SHMEM_USED
#ifndef MPI_USED

	subroutine shmem_logical_get(target,source,len,pe)
	logical	target(len),source(len)
	integer	len,pe

	target(:) = source(:)

	return
	end

#ifndef REAL8
	subroutine shmem_real_get(target,source,len,pe)
	real	target(len),source(len)
	integer	len,pe

	target(:) = source(:)

	return
	end
#endif /* REAL8 */


	subroutine shmem_integer_get(target,source,len,pe)
	integer	target(len),source(len)
	integer	len,pe

	target(:) = source(:)

	return
	end

	subroutine shmem_logical_put(target,source,len,pe)
	logical	target(len),source(len)
	integer	len,pe

	target(:) = source(:)

	return
	end

#ifndef REAL8
	subroutine shmem_real_put(target,source,len,pe)
	real	target(len),source(len)
	integer	len,pe

	target(:) = source(:)

	return
	end
#endif /* REAL8 */

	subroutine shmem_integer_put(target,source,len,pe)
	integer	target(len),source(len)
	integer	len,pe

	target(:) = source(:)

	return
	end

	subroutine shmem_udcflush()
	return
	end

	subroutine barrier()
	return
	end

	subroutine shmem_barrier_all()
	return
	end

	subroutine shmem_real8_min_to_all(target,source,nred,pestart,
     .		pestride,pesize,pwrk,ipsync)
	real	target,source,pwrk
	integer	nred,pestart,pestride,pesize,ipsync

	target = source
	return
	end

	subroutine shmem_real8_max_to_all(target,source,nred,pestart,
     .		pestride,pesize,pwrk,ipsync)
	real	target,source,pwrk
	integer	nred,pestart,pestride,pesize,ipsync

	target = source
	return
	end

	subroutine shmem_real8_sum_to_all(target,source,nred,pestart,
     .		pestride,pesize,pwrk,ipsync)
	real	target,source,pwrk
	integer	nred,pestart,pestride,pesize,ipsync

	target = source
	return
	end

        subroutine shmem_real4_min_to_all(target,source,nred,pestart,
     .          pestride,pesize,pwrk,ipsync)
        real    target,source,pwrk
        integer nred,pestart,pestride,pesize,ipsync

        target = source
        return
        end

        subroutine shmem_real4_max_to_all(target,source,nred,pestart,
     .          pestride,pesize,pwrk,ipsync)
        real    target,source,pwrk
        integer nred,pestart,pestride,pesize,ipsync

        target = source
        return
        end

        subroutine shmem_real4_sum_to_all(target,source,nred,pestart,
     .          pestride,pesize,pwrk,ipsync)
        real    target,source,pwrk
        integer nred,pestart,pestride,pesize,ipsync

        target = source
        return
        end

	subroutine shmem_int8_min_to_all(target,source,nred,pestart,
     .		pestride,pesize,ipwrk,ipsync)
	integer		target,source,ipwrk
	integer	nred,pestart,pestride,pesize,ipsync

	target = source
	return
	end

	subroutine shmem_int8_max_to_all(target,source,nred,pestart,
     .		pestride,pesize,ipwrk,ipsync)
	integer		target,source,ipwrk
	integer	nred,pestart,pestride,pesize,ipsync

	target = source
	return
	end

	subroutine shmem_int8_sum_to_all(target,source,nred,pestart,
     .		pestride,pesize,ipwrk,ipsync)
	integer		target,source,ipwrk
	integer	nred,pestart,pestride,pesize,ipsync

	target = source
	return
	end

        subroutine shmem_int4_min_to_all(target,source,nred,pestart,
     .          pestride,pesize,ipwrk,ipsync)
        integer         target,source,ipwrk
        integer nred,pestart,pestride,pesize,ipsync

        target = source
        return
        end

        subroutine shmem_int4_max_to_all(target,source,nred,pestart,
     .          pestride,pesize,ipwrk,ipsync)
        integer         target,source,ipwrk
        integer nred,pestart,pestride,pesize,ipsync

        target = source
        return
        end

        subroutine shmem_int4_sum_to_all(target,source,nred,pestart,
     .          pestride,pesize,ipwrk,ipsync)
        integer         target,source,ipwrk
        integer nred,pestart,pestride,pesize,ipsync

        target = source
        return
        end

	subroutine shmem_put32(target,source,len,pe)
	real	target(len),source(len)
	integer	len,pe

	target(:) = source(:)

	return
	end

	subroutine shmem_get32(target,source,len,pe)
	real	target(len),source(len)
	integer	len,pe

	target(:) = source(:)

	return
	end

	subroutine shmem_put64(target,source,len,pe)
	real	target(len),source(len)
	integer	len,pe

	target(:) = source(:)

	return
	end

	subroutine shmem_get64(target,source,len,pe)
	real	target(len),source(len)
	integer	len,pe

	target(:) = source(:)

	return
	end


	integer function shmem_my_pe()
	shmem_my_pe = 0
	return
	end

	integer function shmem_n_pes()
	shmem_n_pes = 1
	return
	end


        subroutine start_pes(nprocs)
        integer nprocs
        return
        end


#endif /* SHMEM_USED */
#endif /* MPI_USED */

        subroutine comm_package_used
#ifdef SHMEM_USED
        write(*,*) 'Using SHMEM library.'
#endif /* SHMEM_USED */
#ifdef MPI_USED
        write(*,*) 'Using MPI library.'
#endif /* MPI_USED */
        return
        end subroutine comm_package_used
