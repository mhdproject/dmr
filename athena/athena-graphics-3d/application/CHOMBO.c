/* This file contains routines which write out a file which can be read by the  *
 * chomboVis package.  The routines here convert the PARAMESH data structures   *
 * into those needed by chombVis.                                               */

/* Written: Brian vanStrahlen (LBL) and K. Olson (NASA/GSFC and UMBC).          */

/* Use of this software should credit the chombo project                        */
/*
**   _______              __
**  / ___/ /  ___  __ _  / /  ___
** / /__/ _ \/ _ \/  ' \/ _ \/ _ \
** \___/_//_/\___/_/_/_/_.__/\___/ 
**
**
** This software is copyright (C) by the Lawrence Berkeley
** National Laboratory.  Permission is granted to reproduce
** this software for non-commercial purposes provided that
** this notice is left intact.
** 
** It is acknowledged that the U.S. Government has rights to
** this software under Contract DE-AC03-765F00098 between
** the U.S. Department of Energy and the University of
** California.
**
** This software is provided as a professional and academic
** contribution for joint exchange.  Thus it is experimental,
** is provided ``as is'', with no warranties of any kind
** whatsoever, no support, no promise of updates, or printed
** documentation.  By using this software, you acknowledge
** that the Lawrence Berkeley National Laboratory and
** Regents of the University of California shall have no
** liability with respect to the infringement of other
** copyrights by any part of this software.
**
*/

#include <stdio.h>
#include <stdlib.h>
#include "chomboHDF5.h"
#include "chombo_decs.h"


void chomboinitializefile (int* ndim,
			   char* fileout2,
			   int* numLevels,
			   int* nvar,
			   int* iteration,
			   double* simtime,
			   int* tot_blocks,
			   double* size,
			   int* lrefine,
			   int* nxb,
			   int* nyb,
			   int* nzb,
			   int* min_coordinate,
			   int* max_coordinate,
			   double* dt,
			   int* no_at_level,
			   double* dx,
			   char* compNames)

{
  static int initialized = 1;
  int error;
  int i, j, ivar;
  char fileout[80];
  char compTmp[80];
  int block, level;
  box domain;
  int* refinementRatios;
 
  strncpy(fileout, fileout2, 80);
  for(i=79;i>=0;i--) {
	if (fileout[i] != ' ') { fileout[i+1]=0; break; }
  }
  for(i=(80* *nvar)-1;i>=0;i--) {
	if (compNames[i] == ' ') { compNames[i]=0; }
  }

  refinementRatios = malloc(sizeof(int) * *numLevels);
  for(i = 0; i < *numLevels ; i++)
    {
      refinementRatios[i] = 2;
    }

  if(initialized == 1){
    initialized = 0;
    intvect2d_id = H5Tcreate (H5T_COMPOUND, sizeof(intvect2d));
    H5Tinsert (intvect2d_id, "intvecti", HOFFSET(intvect2d, i), H5T_NATIVE_INT);
    H5Tinsert (intvect2d_id, "intvectj", HOFFSET(intvect2d, j), H5T_NATIVE_INT);

    intvect3d_id = H5Tcreate (H5T_COMPOUND, sizeof(intvect3d));
    H5Tinsert (intvect3d_id, "intvecti", HOFFSET(intvect3d, i), H5T_NATIVE_INT);
    H5Tinsert (intvect3d_id, "intvectj", HOFFSET(intvect3d, j), H5T_NATIVE_INT);
    H5Tinsert (intvect3d_id, "intvectk", HOFFSET(intvect3d, k), H5T_NATIVE_INT);

    box2d_id = H5Tcreate (H5T_COMPOUND, sizeof(box));
    H5Tinsert (box2d_id, "lo_i", HOFFSET(box2d, lo.i), H5T_NATIVE_INT);
    H5Tinsert (box2d_id, "lo_j", HOFFSET(box2d, lo.j), H5T_NATIVE_INT);
    H5Tinsert (box2d_id, "hi_i", HOFFSET(box2d, hi.i), H5T_NATIVE_INT);
    H5Tinsert (box2d_id, "hi_j", HOFFSET(box2d, hi.j), H5T_NATIVE_INT);

    box3d_id = H5Tcreate (H5T_COMPOUND, sizeof(box));
    H5Tinsert (box3d_id, "lo_i", HOFFSET(box3d, lo.i), H5T_NATIVE_INT);
    H5Tinsert (box3d_id, "lo_j", HOFFSET(box3d, lo.j), H5T_NATIVE_INT);
    H5Tinsert (box3d_id, "lo_k", HOFFSET(box3d, lo.k), H5T_NATIVE_INT);
    H5Tinsert (box3d_id, "hi_i", HOFFSET(box3d, hi.i), H5T_NATIVE_INT);
    H5Tinsert (box3d_id, "hi_j", HOFFSET(box3d, hi.j), H5T_NATIVE_INT);
    H5Tinsert (box3d_id, "hi_k", HOFFSET(box3d, hi.k), H5T_NATIVE_INT);
  }
  attr.dim = *ndim;
  attr.precision = precision;
  handle.dim = *ndim;
  handle.precision = precision;
  initializeHDF5datatypes();
  error = Handleopen(&handle, fileout, H5F_ACC_TRUNC);

  ADD_ATTRIBUTE(attr, int, INTEGER, "num_levels", *numLevels);
  ADD_ATTRIBUTE(attr, int, INTEGER, "iteration", *iteration);
  ADD_ATTRIBUTE(attr, int, INTEGER, "num_components", *nvar);
  ADD_ATTRIBUTE(attr, double, DOUBLE, "time", *simtime);
  for(i=0; i < *nvar; ++i)
    {
      sprintf(compTmp, "component_%i",i);
      ADD_CHAR_ATTRIBUTE(attr, compTmp, &compNames[i*80]);
    }
  error = writeHDF5attributes(&attr, handle);
  freeHDF5attributes(&attr);


  if (*ndim == 3) {
    domain.b3.lo.i = min_coordinate[0];
    domain.b3.lo.j = min_coordinate[1];
    domain.b3.lo.k = min_coordinate[2];
    domain.b3.hi.i = max_coordinate[0]+ *nxb-1;
    domain.b3.hi.j = max_coordinate[1]+ *nyb-1;
    domain.b3.hi.k = max_coordinate[2]+ *nzb-1;
  }
  else {
    domain.b2.lo.i = min_coordinate[0];
    domain.b2.lo.j = min_coordinate[1];
    domain.b2.hi.i = max_coordinate[0]+ *nxb-1;
    domain.b2.hi.j = max_coordinate[1]+ *nyb-1;
  }

  /* open and write out the CHOMBO HDF5 formatted file*/
  
  for (i = 0; i < *numLevels; ++i) 
    {
      r =  refinementRatios[i];
      sprintf(compTmp, "/level_%i",i);
      HandleCreateGroup(&handle, compTmp);
      ADD_ATTRIBUTE(attr, double, DOUBLE, "dx", dx[i]);
      ADD_ATTRIBUTE(attr, double, DOUBLE, "dt", *dt);
      if(handle.dim == 2)
	{
	  ADD_ATTRIBUTE(attr, box2d, BOX2D, "prob_domain", domain.b2);
	  refine2(&(domain.b2), r);
	}
      else
	{
	  ADD_ATTRIBUTE(attr, box3d, BOX3D, "prob_domain", domain.b3);
	  refine3(&(domain.b3), r);
	}
      ADD_ATTRIBUTE(attr, int, INTEGER, "ref_ratio",  r);
      error = writeHDF5attributes(&attr, handle);
      freeHDF5attributes(&attr);
      
      /* Create box data space */

      flatdims[0] = no_at_level[i];
      boxdataspace[i] = H5Screate_simple(1, flatdims, NULL);
      
      /* Create doubleing point data space*/
      
      comps = *nvar;
      if(handle.precision == Double) datatype = H5T_NATIVE_DOUBLE;
      else                           datatype = H5T_NATIVE_FLOAT;
      attributeGroup = H5Gcreate(handle.group_ID, "data_attributes", 0);
      aid  = H5Screate(H5S_SCALAR);
      compattr = H5Acreate(attributeGroup, "comps", H5T_NATIVE_INT, aid, H5P_DEFAULT);
      err = H5Awrite(compattr, H5T_NATIVE_INT, &comps);
      H5Gclose(attributeGroup);
      H5Sclose(aid);
      H5Aclose(compattr);    
      flatdims[0] = no_at_level[i]* *nvar * *nxb * *nyb;
      if (handle.dim == 3) flatdims[0] = flatdims[0] * *nzb;
      dataspace[i] = H5Screate_simple(1, flatdims, NULL);
      dataset[i]   = H5Dcreate(handle.group_ID, "data:datatype=0",
			       datatype, dataspace[i], H5P_DEFAULT);
      if (handle.dim == 2) {
	boxdataset[i]   = H5Dcreate(handle.group_ID, "boxes",  box2d_id,
				    boxdataspace[i], H5P_DEFAULT);
      }
      else {
	boxdataset[i]   = H5Dcreate(handle.group_ID, "boxes",  box3d_id,
				    boxdataspace[i], H5P_DEFAULT);
      }
    }
  free(refinementRatios);
}


/* The routine write_block_to_chombo_file gets called from inside the PARAMESH  *
 * source code and writes out the data for a single PARAMESH block              */

void writeblocktochombofile_ (int* block_loc,
				  int* block_tot,
				  int* no_at_level,
				  char* fileout2,
                                  int* tot_blocks,
                                  int* numLevels,
				  int* ndim,
				  int* nvar,
				  int* nxb,
				  int* nyb,
				  int* nzb,
                                  int* lrefine,
                                  int* min_coordinate,
                                  int* max_coordinate,
                                  double* coord,
                                  double* size,
                                  double* unk,
                                  double* simtime,
                                  double* dt,
				  int* first_call,
				  int* no_at_level_tot,
				  double* dx,
				  char* compNames,
				  double* min_domain)
{
  double* unk2;
  int error;
  box  domain;
  int  iteration = 0;
  int ivar;
  int i, j, k, l, var, d, ii;
  int iii, jjj, kkk;
  int coordinate[3];
  int level;
  int offset1, offset2;
  box3d boxp3d;
  box2d boxp2d;

  unk2 = malloc(sizeof(double) * *nxb * *nyb * *nzb * *nvar);

  /**************************************************************************/
  /* Initialize the CHOMBO file                                             */

  if (*first_call == 1) {
    chomboinitializefile(ndim,
			fileout2,
			numLevels,
			nvar,
			&iteration,
			simtime,
			tot_blocks,
			size,
			lrefine,
			nxb,
			nyb,
			nzb,
			min_coordinate,
			max_coordinate,
			dt,
			no_at_level_tot,
			dx,
		        compNames); 
  }
  *first_call = 0;
  
  /**************************************************************************/


  i = lrefine[*block_loc]-1;
  
  if(handle.dim == 2)
    {
      /* Write out box geometry data */
      
      count[0] = 1;
      memdataspace = H5Screate_simple(1, count, NULL);
      for(d = 0; d < *ndim; d++)
	{
	  if (d == 0) {
	    coordinate[d] =
	      (coord[3* *block_loc+d]-(size[3* *block_loc+d]/2)-min_domain[d] +
	       size[3* *block_loc+d]/(2* *nxb))
		  / (size[3* *block_loc+d]/ *nxb);
	  }else if (d == 1) {
	    coordinate[d] =
	      (coord[3* *block_loc+d]-(size[3* *block_loc+d]/2)-min_domain[d] +
	       size[3* *block_loc+d]/(2* *nyb))
	      / (size[3* *block_loc+d]/ *nyb);
	  }
	}
      boxp2d.lo.i = coordinate[0];
      boxp2d.lo.j = coordinate[1];
      boxp2d.hi.i = coordinate[0]+ *nxb-1;
      boxp2d.hi.j = coordinate[1]+ *nyb-1;
      offset[0] = *no_at_level-1;
      ret = H5Sselect_hyperslab (boxdataspace[i], H5S_SELECT_SET, offset, NULL,
				 count, NULL);
      ret = H5Dwrite(boxdataset[i], box2d_id, memdataspace, boxdataspace[i],
		     H5P_DEFAULT, &boxp2d);
      H5Sclose(memdataspace);
	  
      /* write out doubleing data */
      
      offset[0] = (*no_at_level-1) * *nvar * *nxb * *nyb;
      count[0] = *nvar * *nxb * *nyb;
      if (handle.dim == 3) {
	offset[0] = offset[0] * *nzb;
	count[0] = count[0] * *nzb;
      }
      err =  H5Sselect_hyperslab(dataspace[i], H5S_SELECT_SET,
				 offset, NULL,
				 count, NULL);
      if(err != 0)
	{
	  printf("Error selecting %i to %lld size %lld",j, offset[0], count[0]);
	}
      memdataspace = H5Screate_simple(1, count, NULL);
      for(kkk = 0; kkk < *nzb; kkk++)
	{
	  for(jjj = 0; jjj < *nyb; jjj++)
	    {
	      for(iii = 0; iii < *nxb; iii++)
		{
		  for(var = 0; var < *nvar; var++)
		    {
		      offset1 = kkk * *nyb * *nxb * *nvar +
			        jjj * *nxb* *nvar +
			        iii * *nvar +
			        var;
		      offset2 = var * *nzb * *nyb * *nxb +
			        kkk * *nyb * *nxb +
			        jjj * *nxb +
			        iii;
		      unk2[offset2] = unk[offset1];
		    }
		}
	    }
	}
      err = H5Dwrite(dataset[i], datatype, memdataspace, dataspace[i],
		     H5P_DEFAULT, unk2);
      if(err != 0)
	{
	  printf("Error writing %i to %lld size %lld",i, offset[0], count[0]);
	}
      H5Sclose(memdataspace);
    }
  else
    {
      
      /* Write out box geometry data */
      
      count[0] = 1;
      memdataspace = H5Screate_simple(1, count, NULL);
      for(d=0; d<*ndim; d++)
	{
	  if (d == 0) {
	    coordinate[d] =
	      (coord[3* *block_loc+d]-(size[3* *block_loc+d]/2)-min_domain[d] +
	       size[3* *block_loc+d]/(2* *nxb))
	      / (size[3* *block_loc+d]/ *nxb);
	  }else if (d == 1) {
	    coordinate[d] =
	      (coord[3* *block_loc+d]-(size[3* *block_loc+d]/2)-min_domain[d] +
	       size[3* *block_loc+d]/(2* *nyb))
	      / (size[3* *block_loc+d]/ *nyb);
	  }else if (d == 2) {
	    coordinate[d] =
	      (coord[3* *block_loc+d]-(size[3* *block_loc+d]/2)-min_domain[d] +
	       size[3* *block_loc+d]/(2* *nzb))
	      / (size[3* *block_loc+d]/ *nzb);
	  }
	}
      boxp3d.lo.i = coordinate[0];
      boxp3d.lo.j = coordinate[1];
      boxp3d.lo.k = coordinate[2];
      boxp3d.hi.i = coordinate[0]+ *nxb-1;
      boxp3d.hi.j = coordinate[1]+ *nyb-1;
      boxp3d.hi.k = coordinate[2]+ *nzb-1;
      offset[0] = *no_at_level-1;
      ret = H5Sselect_hyperslab (boxdataspace[i], H5S_SELECT_SET, offset, NULL,
				 count, NULL);
      ret = H5Dwrite(boxdataset[i], box3d_id, memdataspace, boxdataspace[i],
		     H5P_DEFAULT, &boxp3d);
      H5Sclose(memdataspace);
      
      
      /* Write out doubleing point data */
      

      offset[0] = (*no_at_level-1) * *nvar * *nxb * *nyb;
      count[0] = *nvar * *nxb * *nyb;
      if (handle.dim == 3) {
	offset[0] = offset[0] * *nzb;
	count[0] = count[0] * *nzb;
      }
      err =  H5Sselect_hyperslab(dataspace[i], H5S_SELECT_SET,
				 offset, NULL,
				 count, NULL);
      if(err != 0)
	{
	  printf("Error selecting %i to %lld size %lld",j, offset[0], count[0]);
	}
      memdataspace = H5Screate_simple(1, count, NULL);
      for(kkk = 0; kkk < *nzb; kkk++)
	{
	  for(jjj = 0; jjj < *nyb; jjj++)
	    {
	      for(iii = 0; iii < *nxb; iii++)
		{
		  for(var = 0; var < *nvar; var++)
		    {
		      offset1 = kkk * *nyb * *nxb * *nvar +
			        jjj * *nxb* *nvar +
			        iii * *nvar +
			        var;
		      offset2 = var * *nzb * *nyb * *nxb +
			        kkk * *nyb * *nxb +
			        jjj * *nxb +
			        iii;
		      unk2[offset2] = unk[offset1];
		    }
		}
	    }
	}
      err = H5Dwrite(dataset[i], datatype, memdataspace, dataspace[i],
		     H5P_DEFAULT, unk2);
      if(err != 0)
	{
	  printf("Error writing %i to %lld size %lld",i, offset[0], count[0]);
	}
      H5Sclose(memdataspace);
    }
  
  if (*block_tot == *tot_blocks-1) {
    for (i = 0; i < *numLevels; i++)
      {
	H5Sclose(boxdataspace[i]);
	H5Dclose(boxdataset[i]);
	H5Sclose(dataspace[i]);
	H5Dclose(dataset[i]); 
      }
    error = Handleclose(handle);
  }
  
  free(unk2);

}




