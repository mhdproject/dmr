sources := \
mpi_amr_1blk_guardcell.F \
mpi_amr_1blk_guardcell_c_to_f.F \
mpi_amr_1blk_restrict.F \
mpi_amr_bsort.F \
mpi_amr_comm_setup.F \
mpi_amr_edge_average.F \
mpi_amr_edge_average_udt.F \
mpi_amr_edge_average_vdt.F \
mpi_amr_edge_diagonal_check.F \
mpi_amr_flux_conserve.F \
mpi_amr_flux_conserve_udt.F \
mpi_amr_flux_conserve_vdt.F \
mpi_amr_get_remote_block.F \
mpi_amr_get_remote_block_fvar.F \
mpi_amr_global_domain_limits.F \
mpi_amr_gsurr_blks.F \
mpi_amr_guardcell.F \
mpi_amr_local_surr_blks.F \
mpi_amr_morton_limits.F \
mpi_amr_prolong.F \
mpi_amr_prolong_fc_divbconsist.F \
mpi_amr_refine_derefine.F \
mpi_amr_restrict.F \
mpi_amr_restrict_fulltree.F \
mpi_amr_restrict_bnd_data_vdt.F \
mpi_amr_restrict_edge_data_vdt.F \
mpi_amr_shmem_lib_1proc.F \
mpi_amr_store_comm_info.F \
mpi_amr_timing_report.F \
mpi_amr_tree_setup.F \
mpi_get_buffer.F \
mpi_get_edge_buffer.F \
mpi_get_flux_buffer.F \
mpi_amr_mirror_blks.F \
mpi_mort_comm_for_surrblks.F \
mpi_morton_bnd.F \
mpi_morton_bnd_fluxcon.F \
mpi_morton_bnd_prolong1.F \
mpi_morton_bnd_restrict.F \
mpi_pack_blocks.F \
mpi_unpack_blocks.F \
mpi_put_buffer.F \
mpi_pack_edges.F \
mpi_pack_fluxes.F \
mpi_put_edge_buffer.F \
mpi_put_flux_buffer.F \
mpi_set_message_limits.F \
mpi_set_message_sizes.F \
mpi_unpack_edges.F \
mpi_unpack_fluxes.F \
rationalize_list.F \
mpi_amr_checkpoint.F \
mpi_amr_derefine_blocks.F \
mpi_amr_morton.F \
mpi_amr_redist_blk.F \
mpi_amr_refine_blocks.F \
mpi_amr_restrict_bnd_data.F \
mpi_amr_restrict_edge_data.F \
mpi_amr_boundary_block_info.F \
mpi_amr_get_bc_settings.F \
mpi_lib.F \
mpi_pack_tree_info.F \
mpi_unpack_tree_info.F \
mpi_wrapper_int.F \
mpi_wrapper_logical.F \
mpi_wrapper_real.F


objects := $(sources:.F=.o)

vpath %.fh ../headers

libmpi_paramesh.a: $(objects)
	$(AR) $(ARFLAGS) $@ $^

ifdef MY_CPP
GNUmakefile.include: $(sources)
	find . -name \*.F | xargs $(MY_CPP) > $@
include GNUmakefile.include
else
$(objects): $(wildcard *.fh)
endif

.PHONY: clean
clean:
	$(RM) libmpi_paramesh.a *.o *.i *~ GNUmakefile.include
