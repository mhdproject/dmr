!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

! timings module
!------------------------------------------------------------------------------

#include "paramesh_preprocessor.fh"

      module timings

        use paramesh_dimensions

        double precision,public :: timer_amr_initialize
        double precision,public :: start_time, end_time


        double precision,public :: timer_amr_refine_derefine
        double precision,public :: timer_amr_check_derefine(0:4)
        double precision,public :: timer_amr_check_refine
        double precision,public :: timer_amr_refine_blocks
        double precision,public :: timer_amr_derefine_blocks
        double precision,public :: timer_amr_morton_order
        double precision,public :: timer_amr_morton_process
        double precision,public :: timer_amr_boundary_block_info
        double precision,public :: timer_amr_global_domain_limits
        double precision,public :: timer_mort_comm_for_surrblks
        double precision,public :: timer_mpi_setup
        double precision,public :: timer_amr_morton_limits
        double precision,public :: timer_amr_gsurr_blks
        double precision,public :: timer_mpi_morton_bnd(0:20)
        double precision,public :: timer_mpi_morton_bnd3(1:7)
        double precision,public :: timer_mpi_morton_bnd_prolong1
        double precision,public :: timer_mpi_morton_bnd_fluxcon
        double precision,public :: timer_mpi_morton_bnd_restrict

        double precision,public :: timer_amr_guardcell
        double precision,public :: timer_amr_1blk_guardcell(0:3)
        double precision,public :: timer_amr_1blk_guardcell_c_to_f
        double precision,public :: timer_amr_1blk_guardcell_srl

        double precision,public :: timer_amr_1blk_cc_cp_remote(0:3)

        double precision,public :: timer_amr_1blk_copy_soln
#ifndef LIBRARY
        double precision,public :: timer_amr_1blk_to_perm(0:1+nvar_work)
#else
        double precision,allocatable,public :: timer_amr_1blk_to_perm(:)
#endif
        double precision,public :: timer_amr_comm_setup(0:9)
        double precision,public :: timer_amr_1blk_cc_prol_gen_unk
        double precision,public :: timer_amr_1blk_cc_prol_gen_work
        double precision,public :: timer_amr_1blk_fc_prol_gen
        double precision,public :: timer_amr_1blk_ec_prol_gen
        double precision,public :: timer_amr_1blk_nc_prol_gen

        double precision,public :: timer_amr_prolong
        double precision,public :: timer_amr_restrict
        double precision,public :: timer_amr_1blk_restrict
        double precision,public :: timer_amr_test_refinement(0:4)
        double precision,public :: timer_advance_soln(0:6)
        double precision,public :: no_of_flops_advance

        integer,public,parameter  :: addflops = 1
        integer,public,parameter  :: mulflops = 1
        integer,public,parameter  :: divflops = 1

        integer,public          :: no_of_calls_check_derefine
        integer,public          :: mess_counter_chk_deref


      end module timings
!-----------------------------------------------------------------
