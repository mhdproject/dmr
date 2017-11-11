
#!
#! Modification history:
#!     Michael L. Rilee, November 2002, *dbz*
#!        Initial support for divergenceless prolongation
#!     Michael L. Rilee, December 2002, *clean_divb*
#!        Support for projecting field onto divergenceless field
#!

.SUFFIXES :  
.SUFFIXES : .o .f .F .F90 .f90 .fh .a

sources := \
  amr_1blk_cc_cp_remote.F  \
  amr_1blk_cc_prol_gen_unk_fun.F  \
  amr_1blk_cc_prol_inject.F  \
  amr_1blk_cc_prol_linear.F  \
  amr_1blk_cc_prol_muscl.F  \
  amr_1blk_cc_prol_genorder.F \
  amr_1blk_cc_prol_gen_work_fun.F \
  amr_1blk_cc_prol_work_inject.F \
  amr_1blk_cc_prol_work_linear.F \
  amr_1blk_cc_prol_work_muscl.F \
  amr_1blk_cc_prol_work_genorder.F \
  amr_1blk_copy_soln.F    \
  amr_1blk_ec_cp_remote.F \
  amr_1blk_ec_prol_gen_fun.F \
  amr_1blk_ec_prol_genorder.F \
  amr_1blk_ec_prol_linear.F \
  amr_1blk_fc_prol_dbz.F \
  sparse_solver.F \
  clean_divb.F \
  amr_1blk_fc_clean_divb.F \
  amr_1blk_fc_cp_remote.F \
  amr_1blk_fc_prol_gen_fun.F \
  amr_1blk_fc_prol_inject.F \
  amr_1blk_fc_prol_linear.F \
  amr_1blk_fc_prol_muscl.F \
  amr_1blk_fc_prol_genorder.F \
  amr_1blk_guardcell_reset.F \
  amr_1blk_guardcell_srl.F \
  amr_1blk_guardcell_f_to_c.F \
  amr_1blk_guardcell_f_to_c_fil.F \
  amr_1blk_guardcell_f_to_c_set.F \
  amr_1blk_nc_cp_remote.F \
  amr_1blk_nc_prol_gen_fun.F \
  amr_1blk_nc_prol_genorder.F \
  amr_1blk_nc_prol_linear.F \
  amr_1blk_save_soln.F \
  amr_1blk_t_to_perm.F \
  amr_1blk_to_perm.F \
  amr_bcset_init.F \
  amr_block_geometry.F \
  amr_close.F \
  amr_initialize.F \
  amr_set_runtime_parameters.F \
  amr_mpi_find_blk_in_buffer.F \
  amr_perm_to_1blk.F \
  amr_prolong_cc_fun_init.F \
  amr_prolong_face_fun_init.F \
  amr_prolong_fun_init.F \
  amr_reorder_grid.F \
  amr_restrict_ec_fun.F \
  amr_restrict_edge.F \
  amr_restrict_fc_fun.F \
  amr_restrict_fulltree.F \
  amr_restrict_red.F \
  amr_restrict_unk_fun.F \
  amr_restrict_unk_fun_recip.F \
  amr_restrict_work_fun.F \
  amr_restrict_work_fun_recip.F \
  amr_system_calls.F \
  gtest_neigh_data1.F \
  mesh_test.F \
  quicksort_index.F

ifdef SHMEM
  sources += \
    amr_1blk_guardcell.F \
    amr_1blk_guardcell_c_to_f.F \
    amr_1blk_restrict.F \
    amr_bc_block.F  \
    amr_bsort.F \
    amr_checkpoint.F \
    amr_derefine_blocks.F \
    amr_edge_average.F \
    amr_edge_average_udt.F \
    amr_edge_average_vdt.F \
    amr_edge_diagonal_check.F \
    amr_flux_conserve.F \
    amr_flux_conserve_udt.F \
    amr_flux_conserve_vdt.F \
    amr_gsurrounding_blks.F \
    amr_guardcell.F  \
    amr_mirror_blks.F \
    amr_morton.F \
    amr_prolong.F \
    amr_prolong_fc_divbconsist.F \
    amr_redist_blk.F \
    amr_refine_blocks.F \
    amr_refine_derefine.F \
    amr_restrict.F \
    amr_restrict_bnd_data.F \
    amr_restrict_bnd_data_vdt.F \
    amr_restrict_edge_data.F \
    amr_restrict_edge_data_vdt.F \
    amr_shmem_lib_1proc.F \
    amr_shmem_reduce.F \
    amr_shmem_sgi_r8.F \
    gtest_neigh_data1.F \
    guardcell_test.F \
    test_neigh_data.F

#    amr_surrounding_blks.F
 
ifndef SHMEM_DUMMYLIB

ifdef SGI
LDFLAGS += -lsma
endif
ifdef COMPAQ
LDFLAGS += -lshmem
endif

endif

endif

objects := \
	$(patsubst %.F90,%.o, \
	$(patsubst %.f90,%.o, \
	  $(patsubst %.F,%.o,$(sources))))

vpath %.fh ../headers

libparamesh.a: $(objects)
	$(AR) $(ARFLAGS) $@ $^

ifdef MY_CPP
GNUmakefile.include: $(sources)
	find . -name \*[.f90,.F90,.F] | xargs $(MY_CPP) > $@
include GNUmakefile.include
else
$(objects): $(wildcard *.fh)
endif

.PHONY: clean
clean:
	$(RM) libparamesh.a *.o *~ GNUmakefile.include

#--------------------------------------------------------------------------

.SUFFIXES : .F90 .f90 .o

.F.o:
	$(FC) $(FFLAGS) -c -o $*.o $*.F

.f90.o:
	$(FC) $(FFLAGS) -c -o $*.o $*.f90

.F90.o:
	$(FC) $(FFLAGS) -c -o $*.o $*.F90

# .fh:;
# .mod:;

#--------------------------------------------------------------------------
