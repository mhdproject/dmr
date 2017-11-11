!----------------------------------------------------------------------
! PARAMESH - an adaptive mesh library.
! Copyright (C) 2003
!
! Use of the PARAMESH software is governed by the terms of the
! usage agreement which can be found in the file
! 'PARAMESH_USERS_AGREEMENT' in the main paramesh directory.
!----------------------------------------------------------------------

# include "paramesh_preprocessor.fh"
      subroutine amr_close

      use paramesh_dimensions
      use physicaldata
      use workspace
#ifdef MPI_USED
      use mpi_morton
#endif /* MPI_USED */
      use tree
      use prolong_arrays
      use timings

      use paramesh_interfaces, only : comm_finish

      implicit none


! This subroutine closes the amr package. 

!----------------------------------------------------------------

#ifdef LIBRARY
      deallocate(unk)
      deallocate(interp_mask_unk)
      deallocate(facevarx)
      deallocate(facevary)
      deallocate(facevarz)
      deallocate(interp_mask_facex)
      deallocate(interp_mask_facey)
      deallocate(interp_mask_facez)
      deallocate(unk_e_x)
      deallocate(unk_e_y)
      deallocate(unk_e_z)
      deallocate(interp_mask_ec)
      deallocate(unk_n)
      deallocate(interp_mask_nc)
      deallocate(unk1)
      deallocate(facevarx1)
      deallocate(facevary1)
      deallocate(facevarz1)
      deallocate(unk_e_x1)
      deallocate(unk_e_y1)
      deallocate(unk_e_z1)
      deallocate(unk_n1)
      deallocate(time_loc)
      deallocate(ldtcomplete)
#if defined(VAR_DT) || defined(PRED_CORR)
      deallocate(t_unk)
      deallocate(tfacevarx)
      deallocate(tfacevary)
      deallocate(tfacevarz)
      deallocate(t_unk_e_x)
      deallocate(t_unk_e_y)
      deallocate(t_unk_e_z)
      deallocate(t_unk_n)
#endif
      deallocate(gt_unk)
#ifdef NO_PERMANENT_GUARDCELLS
      deallocate(gt_facevarx)
      deallocate(gt_facevary)
      deallocate(gt_facevarz)
#else
      deallocate(gt_facevarx)
      deallocate(gt_facevary)
      deallocate(gt_facevarz)
#endif
      deallocate(gt_unk_e_x)
      deallocate(gt_unk_e_y)
      deallocate(gt_unk_e_z)
      deallocate(gt_unk_n)
      deallocate(flux_x)
      deallocate(flux_y)
      deallocate(flux_z)
      deallocate(tflux_x)
      deallocate(tflux_y)
      deallocate(tflux_z)
#ifdef VAR_DT
      deallocate(ttflux_x)
      deallocate(ttflux_y)
      deallocate(ttflux_z)
#endif
      deallocate(bedge_facex_y)
      deallocate(bedge_facex_z)
      deallocate(bedge_facey_x)
      deallocate(bedge_facey_z)
      deallocate(bedge_facez_x)
      deallocate(bedge_facez_y)
      deallocate(recvarx1e)
      deallocate(recvary1e)
      deallocate(recvarz1e)
      deallocate(recvarx2e)
      deallocate(recvary2e)
      deallocate(recvarz2e)
      deallocate(tbedge_facex_y)
      deallocate(tbedge_facex_z)
      deallocate(tbedge_facey_x)
      deallocate(tbedge_facey_z)
      deallocate(tbedge_facez_x)
      deallocate(tbedge_facez_y)
#ifdef VAR_DT
      deallocate(ttbedge_facex_y)
      deallocate(ttbedge_facex_z)
      deallocate(ttbedge_facey_x)
      deallocate(ttbedge_facey_z)
      deallocate(ttbedge_facez_x)
      deallocate(ttbedge_facez_y)
#endif      
#ifdef CURVILINEAR
      deallocate(cell_vol)
      deallocate(cell_area1)
      deallocate(cell_area2)
      deallocate(cell_area3)
      deallocate(cell_leng1)
      deallocate(cell_leng2)
      deallocate(cell_leng3)
#endif
      deallocate(recvarx1)
      deallocate(recvary1)
      deallocate(recvarz1)
      deallocate(recvarxf)
      deallocate(recvaryf)
      deallocate(recvarzf)
      deallocate(bndtempx1)
      deallocate(bndtempy1)
      deallocate(bndtempz1)

! tree data

      deallocate(neigh)
      deallocate(child)
      deallocate(which_child)
      deallocate(type_of_neigh)
      deallocate(type_of_child)
      deallocate(parent)
      deallocate(lrefine)
      deallocate(nodetype)
      deallocate(empty)
      deallocate(bflags)
      deallocate(newchild)
      deallocate(derefine)
      deallocate(refine)
      deallocate(stay)
      deallocate(work_block)
      deallocate(coord)
      deallocate(bsize)
      deallocate(bnd_box)
      deallocate(level_cell_sizes)
      deallocate(laddress)
      deallocate(surr_blks)
#ifdef SAVE_MORTS
      deallocate(surr_morts)
#endif
      deallocate(boundary_box)
      deallocate(boundary_index)

! workspace data

      deallocate(work)
      deallocate(interp_mask_work)
      deallocate(recvw)
      deallocate(sendw)
      deallocate(tempw)
      deallocate(work1)
      deallocate(recvw1)
      deallocate(tempw1)
#ifdef CURVILINEAR
      deallocate(cell_vol_w)
#endif

! morton data

      deallocate(mortonbnd)
      deallocate(laddress_guard)
      deallocate(laddress_prol)
      deallocate(laddress_flux)
      deallocate(laddress_restrict)

! prolong_arrays data

      deallocate(prol_dx)
      deallocate(prol_dy)
      deallocate(prol_dz)
      deallocate(prol_indexx)
      deallocate(prol_indexy)
      deallocate(prol_indexz)
      deallocate(prol_f_dx)
      deallocate(prol_f_dy)
      deallocate(prol_f_dz)
      deallocate(prol_f_indexx)
      deallocate(prol_f_indexy)
      deallocate(prol_f_indexz)
      deallocate(prolw_dx)
      deallocate(prolw_dy)
      deallocate(prolw_dz)
      deallocate(prolw_indexx)
      deallocate(prolw_indexy)
      deallocate(prolw_indexz)


      deallocate(ladd_strt)
      deallocate(ladd_end)

      deallocate(unk1_fl)
      deallocate(facevarx1_fl)
      deallocate(facevary1_fl)
      deallocate(facevarz1_fl)
      deallocate(unk_e_x1_fl)
      deallocate(unk_e_y1_fl)
      deallocate(unk_e_z1_fl)
      deallocate(unk_n1_fl)
      deallocate(work1_fl)

      deallocate(timer_amr_1blk_to_perm)


#endif

      deallocate(i_divf_fc_vars)


! Call the machine/software environment specific closure routine.
! Different versions of comm_finish are provided for use with machines
! which run shmem or mpi. Make sure to compile with the appropriate version
! for your environment.
      call comm_finish()

      return
      end subroutine amr_close



