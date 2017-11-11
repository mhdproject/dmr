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


/*  source file for chomboHDF5

Since this code is meant to be linked in with tcl, it
needs to be compiled 'shared'

hence, for Sun cc, you need

cc -KPIC

and for gcc your need

gcc -fpic

Then use your system dependent linker to 
combine this object and your other properly compiled
objects into a nice shared library.

*/

#include "chomboHDF5.h"
#include <malloc.h>
#include <string.h>
#include <stdio.h>


static const char* TYPE_NAMES[ChTYPES] = {"INTEGER", "FLOAT", 
                                "DOUBLE", "CHAR",  "INTVECT2D", "INTVECT3D", 
                                "BOX2D", "BOX3D" };
void initializeHDF5datatypes()
{
  static int initialized = 1;
  if(initialized == 1){
  initialized = 0;
  intvect2d_id = H5Tcreate (H5T_COMPOUND, sizeof(intvect2d));
  H5Tinsert (intvect2d_id, "intvecti", HOFFSET(intvect2d, i), H5T_NATIVE_INT);
  H5Tinsert (intvect2d_id, "intvectj", HOFFSET(intvect2d, j), H5T_NATIVE_INT);

  intvect3d_id = H5Tcreate (H5T_COMPOUND, sizeof(intvect3d));	
  H5Tinsert (intvect3d_id, "intvecti", HOFFSET(intvect3d, i), H5T_NATIVE_INT);
  H5Tinsert (intvect3d_id, "intvectj", HOFFSET(intvect3d, j), H5T_NATIVE_INT);
  H5Tinsert (intvect3d_id, "intvectk", HOFFSET(intvect3d, k), H5T_NATIVE_INT);

  /* old composite-of-composite style boxes
  box2d_id = H5Tcreate (H5T_COMPOUND, sizeof(box));
  H5Tinsert (box2d_id, "smallend", HOFFSET(box2d, lo), intvect2d_id);
  H5Tinsert (box2d_id, "bigend",   HOFFSET(box2d, hi), intvect2d_id);

  box3d_id = H5Tcreate (H5T_COMPOUND, sizeof(box));
  H5Tinsert (box3d_id, "smallend", HOFFSET(box3d, lo), intvect3d_id);
  H5Tinsert (box3d_id, "bigend",   HOFFSET(box3d, hi), intvect3d_id);
  */
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
}


void printAttributes(HDF5attributes* attr)
{
  Attribute* at;
  int i;
  for(i=0; i<ChTYPES; ++i)
    {
      at = attr->accessByType[i];
      while(at != NULL)
	{
	  printf("%s attribute %s \n",
		TYPE_NAMES[i], 
		 at->name);
	  at = at->next;
	}
    }
}


void initHDF5attributes(HDF5attributes* attr)
{
  int i;
  for(i=0; i<ChTYPES; ++i)
    {
       attr->accessByType[i] = NULL;
       attr->numByType[i] = 0;
    }
}

void freeDataArray(FAB* ArrayBoxCArray,
		   int length)
{

  /* exploit the fact that I allocate the entire array and hand
     the pointer off to the first entry in the array */
   if(length > 0)
     {
       free(ArrayBoxCArray->f2.data);
     }
}

int numPnts2(const box2d* box)
{
  return (box->hi.i-box->lo.i+1)*(box->hi.j-box->lo.j+1);
}

int numPnts3(const box3d* box)
{
  return (box->hi.i-box->lo.i+1)*(box->hi.j-box->lo.j+1)*(box->hi.k-box->lo.k+1);
}

void refine2(box2d* b2, int refinement)
{
  b2->lo.i*=refinement;
  b2->lo.j*=refinement;
  b2->hi.i = (b2->hi.i + 1)*refinement - 1;
  b2->hi.j = (b2->hi.j + 1)*refinement - 1;
}

void refine3(box3d* b3, int refinement)
{
  b3->lo.i*=refinement;
  b3->lo.j*=refinement;
  b3->lo.k*=refinement;
  b3->hi.i = (b3->hi.i + 1)*refinement - 1;
  b3->hi.j = (b3->hi.j + 1)*refinement - 1;
  b3->hi.k = (b3->hi.k + 1)*refinement - 1;
}

int Handleopen(HDF5Handle* handle, const char* filename, hid_t access)
{
  herr_t ret = 0;     
  hid_t attr, datatype, aid; 
  float fl = 0;
  double dl = 0;

  if(access == H5F_ACC_RDONLY)
    {
      handle->file_ID = H5Fopen(filename, access, H5P_DEFAULT);
      handle->group_ID  =  H5Gopen(handle->file_ID, "Chombo_global");
      if(handle->group_ID < 0) ret = 1;
      attr = H5Aopen_name(handle->group_ID, "SpaceDim");
      ret  = H5Aread(attr, H5T_NATIVE_INT, &(handle->dim));
      ret  = H5Aclose(attr);
      
      attr = H5Aopen_name(handle->group_ID, "testReal");
      datatype = H5Aget_type(attr);
      if(H5Tget_precision(datatype) == H5Tget_precision(H5T_NATIVE_FLOAT))
	handle->precision = Float;
      else
	handle->precision = Double;
    }
  else if(access == H5F_ACC_TRUNC)
    {
      handle->file_ID = H5Fcreate(filename,  H5F_ACC_TRUNC,  H5P_DEFAULT, H5P_DEFAULT);
      handle->group_ID  =  H5Gcreate(handle->file_ID,  "Chombo_global", 0);
      aid  = H5Screate(H5S_SCALAR);
      attr = H5Acreate( handle->group_ID, "SpaceDim", H5T_NATIVE_INT, aid, H5P_DEFAULT);
      ret = H5Awrite(attr, H5T_NATIVE_INT, &(handle->dim));
      H5Aclose(attr);
      if(handle->precision == Float)
	{
	  
	  attr =  H5Acreate(handle->group_ID, "testReal", 
			    H5T_NATIVE_FLOAT, aid, H5P_DEFAULT);
	  ret = H5Awrite(attr, H5T_NATIVE_FLOAT, &fl);
	}
      else
	{
	  attr =  H5Acreate(handle->group_ID, "testReal", 
			    H5T_NATIVE_DOUBLE, aid, H5P_DEFAULT);
	  ret = H5Awrite(attr, H5T_NATIVE_DOUBLE, &dl);
	} 
    }    
  if(handle->group_ID > 0) H5Gclose(handle->group_ID);
  handle->group_ID = H5Gopen(handle->file_ID, "/");
  return 0;
}

int HandlesetGroup(HDF5Handle* handle, const char* group)
{
  H5Gclose(handle->group_ID);
  handle->group_ID = H5Gopen(handle->file_ID, group);
  if(handle->group_ID > 0) return 0;
  return -1;
}

int HandleCreateGroup(HDF5Handle* handle, const char* group)
{
  H5Gclose(handle->group_ID);
  handle->group_ID = H5Gcreate(handle->file_ID, group, 0);
  if(handle->group_ID > 0) return 0;
  return -1;
}

int Handleclose(HDF5Handle handle)
{
  H5Gclose(handle.group_ID);
  return H5Fclose(handle.file_ID);
}


void freeHDF5attributes(HDF5attributes* attributes)
{
  int i;  Attribute *next, *tmp;
  for(i=0; i<ChTYPES; ++i)
    {
      next = attributes->accessByType[i];
      while(next != NULL)
	{
	  tmp = next;
	  next=tmp->next;
	  free(tmp->data);
	  free(tmp->name);
	  free(tmp);
	}
      attributes->accessByType[i] = NULL;
      attributes->numByType[i]    = 0;
    }

}

herr_t attributeScan(hid_t loc_id, const char *name, void *opdata)
{
  herr_t ret = 0;
  HDF5attributes* data = opdata;
  hid_t attr, atype, aclass;
  size_t size = 0;
  Attribute* attribute, *head; 

  size = strlen(name) + 1;
  attribute = malloc(sizeof(Attribute));
  attribute->name = malloc(size);
  memcpy(attribute->name, name, size);
 

  attr   = H5Aopen_name(loc_id, name);
  atype  = H5Aget_type(attr);
  aclass = H5Tget_class(atype);

  switch(aclass){
  case H5T_INTEGER :
    data->numByType[INTEGER]++;
    attribute->attributeType = INTEGER;
    attribute->data = malloc(sizeof(int));
    ret = H5Aread(attr, H5T_NATIVE_INT, attribute->data);
    break;
  case H5T_FLOAT:
    if(data->precision == Float)
      {
	data->numByType[FLOAT]++;
	attribute->attributeType = FLOAT;
	attribute->data = malloc(sizeof(float));
	ret = H5Aread(attr, H5T_NATIVE_FLOAT, attribute->data);
      }
    else
      {
	data->numByType[DOUBLE]++;
	attribute->attributeType = DOUBLE;
	attribute->data = malloc(sizeof(double));
	ret = H5Aread(attr, H5T_NATIVE_DOUBLE, attribute->data);
      }	
    break;
  case H5T_STRING:
    data->numByType[CHAR]++;
    size = H5Tget_size(atype);
    attribute->attributeType = CHAR;
    attribute->data = malloc(size+1);
    ret = H5Aread(attr, atype, attribute->data );
    if(ret < 0) break;
    ((char *)(attribute->data))[size] = '\0';
    break;
  case H5T_COMPOUND:
    if(strcmp(H5Tget_member_name(atype, 0), "lo_i") == 0)
      {
	if(data->dim == 2)
	  {
	    data->numByType[BOX2D]++;
	    attribute->attributeType = BOX2D;
	    attribute->data = malloc(sizeof(box2d));
	    ret = H5Aread(attr, box2d_id, attribute->data);
	  }
	else if(data->dim == 3)
	  {
	    data->numByType[BOX3D]++;
	    attribute->attributeType = BOX3D;
	    attribute->data = malloc(sizeof(box3d));
	    ret = H5Aread(attr, box3d_id, attribute->data);
	  }
	break;
      }
    else if(strcmp(H5Tget_member_name(atype, 0), "intvecti") == 0)
      {
	if(data->dim == 2)
	  {
	    data->numByType[INTVECT2D]++;
	    attribute->attributeType = INTVECT2D;
	    attribute->data = malloc(sizeof(intvect2d));
	    ret = H5Aread(attr, intvect2d_id, attribute->data);
	  }
	else if(data->dim == 3)
	  {
	    data->numByType[INTVECT3D]++;
	    attribute->attributeType = INTVECT3D;
	    attribute->data = malloc(sizeof(intvect3d));
	    ret = H5Aread(attr, intvect3d_id, attribute->data);
	  }
	break;
      }
  default:
    /* don't know what the hell this thing is */
    free(attribute->name);
    free(attribute);
    return -1;
  }
  
  /* OK, lets tack this attribute to the right linked-list */
  head = data->accessByType[attribute->attributeType];
  data->accessByType[attribute->attributeType] = attribute;
  attribute->next = head;
  return ret;
}

int readHDF5attributes(HDF5attributes* attr, HDF5Handle handle)
{
  int i;

  attr->dim = handle.dim;
  attr->precision = handle.precision;
  for(i=0; i<ChTYPES; ++i)
    {
      attr->numByType[i] = 0;
      attr->accessByType[i] = NULL;
    }
  return H5Aiterate(handle.group_ID, NULL, attributeScan, attr);
}

int readBoxes(box** boxes, int* length, HDF5Handle handle)
{
  hid_t boxdataset, boxdataspace, memdataspace;
  hsize_t dims[1], maxdims[1];
  herr_t error;

  boxdataset = H5Dopen(handle.group_ID, "boxes");
  if(boxdataset < 0) return boxdataset;
  boxdataspace =  H5Dget_space(boxdataset);
  if(boxdataspace < 0) return boxdataspace;

  H5Sget_simple_extent_dims(boxdataspace, dims, maxdims);
  
  memdataspace = H5Screate_simple(1, dims, NULL);

  *length = dims[0];
  *boxes = malloc(dims[0]*sizeof(box));
  if(handle.dim == 2)
    {
      error = H5Dread(boxdataset, box2d_id, memdataspace, boxdataspace,
		      H5P_DEFAULT, *boxes);
    }
  else if (handle.dim == 3)
    {
      error = H5Dread(boxdataset, box3d_id, memdataspace, boxdataspace,
		      H5P_DEFAULT, *boxes);
    }
  
  H5Dclose(boxdataset);
  H5Sclose(boxdataspace);
  H5Sclose(memdataspace);
  return 0;
}

int writeBoxes(box* boxes, int length, HDF5Handle handle)
{
  herr_t ret;
  hssize_t offset[1];
  hsize_t  flatdims[1], count[1];
  hid_t boxdataspace, boxdataset, memdataspace;
  int i;

  count[0] = 1;
  flatdims[0] = length;

  boxdataspace = H5Screate_simple(1, flatdims, NULL);
  memdataspace = H5Screate_simple(1, count, NULL);
  if(handle.dim == 2)
    {
      boxdataset   = H5Dcreate(handle.group_ID, "boxes",  box2d_id, 
				     boxdataspace, H5P_DEFAULT);
      if(boxdataset < 0) return boxdataset;
      for(i=0; i<length; ++i)
	{
	   offset[0] = i;
	   ret = H5Sselect_hyperslab (boxdataspace, H5S_SELECT_SET, offset, NULL, 
				      count, NULL);
	   if(ret < 0) return ret;
	   ret = H5Dwrite(boxdataset, box2d_id, memdataspace, boxdataspace,
			  H5P_DEFAULT, boxes + i);
	   if(ret < 0) return ret;
	}
    }      
  else
    {
      boxdataset   = H5Dcreate(handle.group_ID, "boxes",  box3d_id, 
				     boxdataspace, H5P_DEFAULT);
      if(boxdataset < 0) return boxdataset;
      for(i=0; i<length; ++i)
	{
	  offset[0] = i;
	  ret = H5Sselect_hyperslab (boxdataspace, H5S_SELECT_SET, offset, NULL, 
				     count, NULL);
	  if(ret < 0) return ret;
	  ret = H5Dwrite(boxdataset, box3d_id, memdataspace, boxdataspace,
			 H5P_DEFAULT, boxes + i);
	  if(ret < 0) return ret;
	}
    }

  H5Sclose(boxdataspace);
  H5Sclose(memdataspace);
  H5Dclose(boxdataset);
  return 0;
}

int writeFABBoxes(FAB* FABArray, int length, HDF5Handle handle)
{
  herr_t ret;
  hssize_t offset[1];
  hsize_t  flatdims[1], count[1];
  hid_t boxdataspace, boxdataset, memdataspace;
  int i;

  count[0] = 1;
  flatdims[0] = length;

  boxdataspace = H5Screate_simple(1, flatdims, NULL);
  memdataspace = H5Screate_simple(1, count, NULL);
  if(handle.dim == 2)
    {
      boxdataset   = H5Dcreate(handle.group_ID, "boxes",  box2d_id, 
				     boxdataspace, H5P_DEFAULT);
      if(boxdataset < 0) return boxdataset;
      for(i=0; i<length; ++i)
	{
	   offset[0] = i;
	   ret = H5Sselect_hyperslab (boxdataspace, H5S_SELECT_SET, offset, NULL, 
				      count, NULL);
	   if(ret < 0) return ret;
	   ret = H5Dwrite(boxdataset, box2d_id, memdataspace, boxdataspace,
			  H5P_DEFAULT, &(FABArray[i].f2.box));
	   if(ret < 0) return ret;
	}
    }      
  else
    {
      boxdataset   = H5Dcreate(handle.group_ID, "boxes",  box3d_id, 
				     boxdataspace, H5P_DEFAULT);
      if(boxdataset < 0) return boxdataset;
      for(i=0; i<length; ++i)
	{
	  offset[0] = i;
	  ret = H5Sselect_hyperslab (boxdataspace, H5S_SELECT_SET, offset, NULL, 
				     count, NULL);
	  if(ret < 0) return ret;
	  ret = H5Dwrite(boxdataset, box3d_id, memdataspace, boxdataspace,
			 H5P_DEFAULT, &(FABArray[i].f3.box));
	  if(ret < 0) return ret;
	}
    }

  H5Sclose(boxdataspace);
  H5Sclose(memdataspace);
  H5Dclose(boxdataset);
  return 0;
}


int readData(FAB** arrayboxes, box* boxes, int length, HDF5Handle handle)
{
  hid_t attributeGroup, dataset, dataspace, memdataspace, compattr;
  hsize_t dims[1], maxdims[1];
  herr_t error = 0; 
  int comps, i;
  FloatArrayBox2D* fa2d;
  DoubleArrayBox2D* da2d;
  FloatArrayBox3D* fa3d;
  DoubleArrayBox3D* da3d;
  
  attributeGroup = H5Gopen(handle.group_ID, "data_attributes");
  compattr = H5Aopen_name(attributeGroup, "comps");
  if(compattr < 0) return compattr;
  error = H5Aread(compattr, H5T_NATIVE_INT, &comps);
  if(error != 0) return error;

  dataset = H5Dopen(handle.group_ID, "data:datatype=0");
  if(dataset < 0) return dataset;
  dataspace =  H5Dget_space(dataset);
  if(dataspace < 0) return dataspace;
  
  H5Sget_simple_extent_dims(dataspace, dims, maxdims);
  
  memdataspace = H5Screate_simple(1, dims, NULL);
  *arrayboxes = malloc(length*sizeof(FAB));

  if(handle.dim == 2)
    {
      if(handle.precision == Float)
	{
	  fa2d = &((**arrayboxes).f2);
	  fa2d->data = malloc(dims[0]*sizeof(float));
	  fa2d->box = (boxes[0]).b2;
	  fa2d->ncomps = comps;
	  error = H5Dread(dataset, H5T_NATIVE_FLOAT, memdataspace, dataspace,
			  H5P_DEFAULT, fa2d->data);
	  for(i=1; i<length; ++i){
	    ((*arrayboxes)[i]).f2.ncomps = comps;
	    ((*arrayboxes)[i]).f2.box = (boxes[i]).b2;
	    ((*arrayboxes)[i]).f2.data = 
	      ((*arrayboxes)[i-1]).f2.data + 
	      numPnts2(&(boxes[i-1]).b2)*comps;
	  }

	}
      else
	{
	  da2d = &((**arrayboxes).d2);
	  da2d->data = malloc(dims[0]*sizeof(double));
	  da2d->box = (boxes[0]).b2;
	  da2d->ncomps = comps;
	  error = H5Dread(dataset, H5T_NATIVE_DOUBLE, memdataspace, dataspace,
			  H5P_DEFAULT, da2d->data);
	  for(i=1; i<length; ++i){
	    ((*arrayboxes)[i]).d2.ncomps = comps;
	    ((*arrayboxes)[i]).d2.box = (boxes[i]).b2;
	    ((*arrayboxes)[i]).d2.data = 
	      ((*arrayboxes)[i-1]).d2.data + 
	      numPnts2(&(boxes[i-1]).b2)*comps;
	  }

	}

    }
  else if(handle.dim == 3)
    {
     if(handle.precision == Float)
	{
	  fa3d = &((**arrayboxes).f3);
	  fa3d->data = malloc(dims[0]*sizeof(float));
	  fa3d->box = (boxes[0]).b3;
	  fa3d->ncomps = comps;
	  error = H5Dread(dataset, H5T_NATIVE_FLOAT, memdataspace, dataspace,
			  H5P_DEFAULT, fa3d->data);
	  for(i=1; i<length; ++i){
	    ((*arrayboxes)[i]).f3.ncomps = comps;
	    ((*arrayboxes)[i]).f3.box = (boxes[i]).b3;
	    ((*arrayboxes)[i]).f3.data = 
	      ((*arrayboxes)[i-1]).f3.data + 
	      numPnts3(&(boxes[i-1]).b3)*comps;
	  }

	}
      else
	{
	  da3d = &((**arrayboxes).d3);
	  da3d->data = malloc(dims[0]*sizeof(double));
	  da3d->box = (boxes[0]).b3;
	  da3d->ncomps = comps;
	  error = H5Dread(dataset, H5T_NATIVE_DOUBLE, memdataspace, dataspace,
			  H5P_DEFAULT, da3d->data);
	  for(i=1; i<length; ++i){
	    (((*arrayboxes))[i]).d3.ncomps = comps;
	    (((*arrayboxes))[i]).d3.box = (boxes[i]).b3;
	    (((*arrayboxes))[i]).d3.data = 
	      (((*arrayboxes))[i-1]).d3.data + 
	      numPnts3(&(boxes[i-1]).b3)*comps;
	  }

	}


    }
  else return 2;
  
  H5Gclose(attributeGroup);
  H5Aclose(compattr);
  H5Dclose(dataset);
  H5Sclose(dataspace);
  H5Sclose(memdataspace);

  return error;
}

int writeHDF5attributes(HDF5attributes* attrib, HDF5Handle handle)
{
  H5E_auto_t efunc; void* edata;
  herr_t  ret;
  hid_t aid, attr, s_type;
  Attribute* atr;
  Attribute* at;

  
  H5Eget_auto(&efunc, &edata);
  
#define INSERT(Ttype, attributePtr, H5Ttype)                              \
  at = attributePtr;                                                      \
  while(at != NULL)                                                       \
    {                                                                     \
      aid  = H5Screate(H5S_SCALAR);                                       \
      H5Eset_auto(NULL, NULL);                                            \
      attr = H5Acreate(handle.group_ID, at->name, H5Ttype,                \
			     aid, H5P_DEFAULT);                           \
      if(attr < 0)                                                        \
	{                                                                 \
	  H5Adelete(handle.group_ID, at->name);                           \
	  attr = H5Acreate(handle.group_ID, at->name, H5Ttype,            \
			   aid, H5P_DEFAULT);                             \
	  if(attr < 0)                                                    \
	    {                                                             \
	      return -1;	       			                  \
	    }                                                             \
	}                                                                 \
      H5Eset_auto(efunc, edata);                                          \
      ret = H5Awrite(attr, H5Ttype, at->data);                            \
      if(ret < 0) return ret;                                             \
      H5Sclose(aid);                                                      \
      H5Aclose(attr);                                                     \
      at = (Attribute *)at->next;                                         \
    }                                                                     \
                                                              
    atr = attrib->accessByType[INTEGER];
    INSERT(int,  atr, H5T_NATIVE_INT);

    atr = attrib->accessByType[FLOAT];
    INSERT(float,   atr,   H5T_NATIVE_FLOAT);

    atr = attrib->accessByType[DOUBLE];
    INSERT(double,  atr,  H5T_NATIVE_DOUBLE);

    atr = attrib->accessByType[INTVECT2D];
    INSERT(intvect2d, atr,  intvect2d_id);

    atr = attrib->accessByType[INTVECT3D];
    INSERT(intvect3d, atr,  intvect3d_id);

    atr = attrib->accessByType[BOX2D];
    INSERT(box2d, atr,  box2d_id);

    atr = attrib->accessByType[BOX3D];
    INSERT(box3d, atr,  box3d_id);
    
    /* string is different, of course */
    at = attrib->accessByType[CHAR];
     while(at != NULL)                                                  
    {                                             
      aid    = H5Screate(H5S_SCALAR); 
      s_type = H5Tcopy(H5T_C_S1);
      H5Tset_size(s_type, strlen((char*)at->data)); /*extra requirement for strings*/
      H5Eset_auto(NULL, NULL);                                        
      attr = H5Acreate(handle.group_ID, at->name, s_type,         
			     aid, H5P_DEFAULT);                      
      if(attr < 0)                                                      
	{                                                              
	  H5Adelete(handle.group_ID, at->name);                         
	  attr = H5Acreate(handle.group_ID, at->name, s_type,        
			   aid, H5P_DEFAULT);                        
	  if(attr < 0)                                                 
	    {                                                           
	      return -1;	       					 
	    }                                                          
	}                                                                
      H5Eset_auto(efunc, edata);                                      
      ret = H5Awrite(attr, s_type, at->data);                          
      if(ret < 0) return ret;                                          
      H5Sclose(aid);                                                   
      H5Aclose(attr); 
      H5Tclose(s_type);
      at = (Attribute *)at->next;
    }                                 
 

    return 0;

}

void calcOffsets(hsize_t* offsets, FAB* FABArray, int length, int dim, ChPrecision precision )
{
  int i;
  int comps = FABArray[0].f2.ncomps;
  offsets[0] = 0;
  
    
  if(dim == 2)
    {
      if(precision == Float)
	{
	  for(i=1; i<=length; ++i)
	    {
	      box2d* box = &(FABArray[i-1].f2.box);
	      offsets[i] = offsets[i-1] + numPnts2(box)*comps;
	    }
	}
      else
	{
	  for(i=1; i<=length; ++i)
	    {
	      box2d* box = &(FABArray[i-1].d2.box);
	      offsets[i] = offsets[i-1] + numPnts2(box)*comps;
	    }
	}
    }
  else /* dim == 3 */
    {
      if(precision == Float)
	{
	  for(i=1; i<=length; ++i)
	    {
	      box3d* box = &(FABArray[i-1].f3.box);
	      offsets[i] =  offsets[i-1] +  numPnts3(box)*comps;
	    }
	}
      else
	{
	  for(i=1; i<=length; ++i)
	    {
	      box3d* box = &(FABArray[i-1].d3.box);
	      offsets[i] = offsets[i-1] +  numPnts3(box)*comps;
	    }
	}
    }
}


int writeData(FAB* FABArray, int length, HDF5Handle handle)
{
  hsize_t* FABoffsets;
  hsize_t flatdims[1];
  herr_t err;
  hsize_t count[1];
  hssize_t offset[1];
  int i;
  hid_t attributeGroup, compattr, aid;

  hid_t dataspace, dataset, memdataspace, datatype;
  int comps = FABArray[0].f2.ncomps;

  if(handle.precision == Float) datatype = H5T_NATIVE_FLOAT;
  else                          datatype = H5T_NATIVE_DOUBLE;

  FABoffsets = (hsize_t*)malloc((length+1)*sizeof(hsize_t));
  calcOffsets(FABoffsets, FABArray, length, handle.dim, handle.precision);

  attributeGroup = H5Gcreate(handle.group_ID, "data_attributes", 0);
  aid  = H5Screate(H5S_SCALAR);
  compattr = H5Acreate(attributeGroup, "comps", H5T_NATIVE_INT, aid, H5P_DEFAULT);
  err = H5Awrite(compattr, H5T_NATIVE_INT, &comps);
  H5Gclose(attributeGroup);
  H5Sclose(aid);
  H5Aclose(compattr);

  flatdims[0] = FABoffsets[length];
  dataspace = H5Screate_simple(1, flatdims, NULL);
  dataset   = H5Dcreate(handle.group_ID, "data:datatype=0",  
			datatype, dataspace, H5P_DEFAULT);
  
  for(i=0 ; i<length; ++i)
    {
      offset[0] = FABoffsets[i];
      count[0] = FABoffsets[i+1] - offset[0];
      if(count[0] <= 0) return -4;
      err =  H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, 
				 offset, NULL, 
				 count, NULL);
	  if(err != 0)
		{
		  printf("Error selecting %i to %lld size %lld",i, offset[0], count[0]);
		  return -8;
		}
      memdataspace = H5Screate_simple(1, count, NULL);
      err = H5Dwrite(dataset, datatype, memdataspace, dataspace,
		     H5P_DEFAULT, FABArray[i].f2.data);
	  if(err != 0)
		{
		  printf("Error writing %i to %lld size %lld",i, offset[0], count[0]);
		  return -5;
		}
      H5Sclose(memdataspace);
    }


  free(FABoffsets);
  H5Sclose(dataspace);
  H5Dclose(dataset);

  

  return 0;
}

int writeAMRHierarchy(const char* filename,
		      int dim,
		      ChPrecision precision,
		      int numLevels,
		      FAB** levels,
		      int* length,
		      char** compNames,
		      box  domain,
		      int iteration,
		      double dx,
		      double dt,
		      double time,
		      int* refinementRatios)
{
  int error;
  int i, r, l;

  char compTmp[20];
  HDF5Handle handle;
  HDF5attributes attr;

  initHDF5attributes(&attr);

  attr.dim = dim;
  attr.precision = precision;
  handle.dim = dim;
  handle.precision = precision;

  initializeHDF5datatypes();

  error = Handleopen(&handle, filename, H5F_ACC_TRUNC);
  if (error != 0) {
    return(-1);
  }

  ADD_ATTRIBUTE(attr, int, INTEGER, "num_levels", numLevels);
  ADD_ATTRIBUTE(attr, int, INTEGER, "iteration", iteration);
  ADD_ATTRIBUTE(attr, int, INTEGER, "num_components", levels[0][0].f2.ncomps);
  ADD_ATTRIBUTE(attr, double, DOUBLE, "time", time);

  for(i=0; i< levels[0][0].f2.ncomps; ++i)
    {
      sprintf(compTmp, "component_%i",i);
      ADD_CHAR_ATTRIBUTE(attr, compTmp, compNames[i]);
    }

  error = writeHDF5attributes(&attr, handle);
  if(error != 0){
    return -2;
  }

  freeHDF5attributes(&attr);
  for(i=0; i<numLevels; ++i)
    {
      r =  refinementRatios[i];
      sprintf(compTmp, "/level_%i",i);
      HandleCreateGroup(&handle, compTmp);
      ADD_ATTRIBUTE(attr, double, DOUBLE, "dx", dx);
      ADD_ATTRIBUTE(attr, double, DOUBLE, "dt", dt);
      dx/=r;
      dt/=r;
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
      if(error != 0){
	return -2;
      }
      freeHDF5attributes(&attr);

      error = writeFABBoxes(levels[i], length[i], handle);
      if(error != 0)
	{
	  return -4;
	}
      error = writeData(levels[i], length[i], handle);
      if(error != 0)
	{
	  return -5;
	} 
    }
  error = Handleclose(handle);
  if(error != 0)
    return -20;

  return 0;
}

int readAMRHierarchy(const char* filename,
		     int* dim,
		     ChPrecision* precision,
		     int* num_levels,
		     FAB*** levels,
		     int** length,
		     char*** compNames,
		     box*  domain,
		     int* iteration,
		     double* dx,
		     double* dt,
		     double* time,
		     int** refinementRatios)
{
  int error;
  int i;

  box* boxes;
  HDF5Handle handle;
  int num_components = -1;
  HDF5attributes attr;
  Attribute *at;
  char tempo[10];
  initializeHDF5datatypes();

  error = Handleopen(&handle, filename, H5F_ACC_RDONLY);
  if (error != 0) {
    return(-1);
  }

  *dim = handle.dim;
  *precision = handle.precision;

  error = readHDF5attributes(&attr, handle);

  if (error != 0) {
    return -2;
  }

  at = attr.accessByType[INTEGER];
  while (at != NULL) {
    if (strcmp(at->name, "num_levels") == 0) {
      *num_levels = *((int*)(at->data));
    }
    if (strcmp(at->name, "iteration") == 0) {
      *iteration = *((int*)(at->data));
    }
    if (strcmp(at->name, "num_components") == 0) {
      num_components = *((int*)(at->data));
    }
    
    at = (Attribute *)at->next;
  }

  at = attr.accessByType[DOUBLE];
  while (at != NULL) {
    if (strcmp(at->name, "time") == 0) {
      *time = *((double*)(at->data));
    }
    at = (Attribute *)at->next;
  }

  if (num_levels == NULL) {
    return -3;
  }

  if (num_components == -1) {
    return -4;
  }

  if (*num_levels <= 0) {
    return -5;
  }

  if (num_components <= 0) {
    return -6;
  }

  *levels = (FAB**)malloc(*num_levels * sizeof(FAB*));
  *compNames = (char**)malloc(num_components * sizeof(char*));
  *length = (int*)malloc(sizeof(int) * *num_levels);
  *refinementRatios = (int*)malloc(sizeof(int) * *num_levels);

  for (i = 0; i < num_components; i++) {
    (*compNames)[i] = NULL;
  }
  
  at = attr.accessByType[CHAR];
  while (at != NULL) {
    if (strncmp(at->name, "component_", 10) == 0) {
      sscanf(at->name+10,"%d",&i);
      (*compNames)[i] = malloc(strlen((const char*)(at->data)));
      strcpy((*compNames)[i], (const char*)(at->data));
    }

    at = (Attribute *)at->next;
  }

  for (i = 0; i < num_components; i++) {
    if ( (*compNames)[i] == NULL) {
      return -7;
    }
  }


  freeHDF5attributes(&attr);



  for(i=0; i<*num_levels; ++i)
    {
      sprintf(tempo, "level_%i",i);
      HandlesetGroup(&handle, tempo);
      error = readHDF5attributes(&attr, handle);
      if( i == 0){
	at = attr.accessByType[DOUBLE];
	while (at != NULL) {
	  if (strcmp(at->name, "dx") == 0) {
	    *dx = *((double*)(at->data));
	  }
	  if (strcmp(at->name, "dt") == 0) {
	    *dt = *((double*)(at->data));
	  }
	  
	at = (Attribute *)at->next;
	}
	if(handle.dim == 2)
	  {
	    at = attr.accessByType[BOX2D];
	    while (at != NULL) {
	      if (strcmp(at->name, "prob_domain") == 0) {
		(*domain).b2 = *((box2d*)(at->data));
	      }
	      at = (Attribute *)at->next;
	    }
	  }
	else
	  {
	    at = attr.accessByType[BOX3D];
	    while (at != NULL) {
	      if (strcmp(at->name, "prob_domain") == 0) {
		(*domain).b3 = *((box3d*)(at->data));
	      }
	      at = (Attribute *)at->next;
	    }
	  }
      }
      at = attr.accessByType[INTEGER];
      while (at != NULL) {
	if (strcmp(at->name, "ref_ratio") == 0) {
	  (*refinementRatios)[i] = *((int*)(at->data));
	} 
	at = (Attribute *)at->next;	
      }
     
      error = readBoxes(&boxes, &((*length)[i]), handle);

      if(error != 0)
	{
	  return -8;
	}

      error = readData(&((*levels)[i]), boxes, (*length)[i], handle);

      if(error != 0)
	{
	  return -100+error;
	}

      free(boxes);
      freeHDF5attributes(&attr);
    }
  error = Handleclose(handle);
  if(error != 0)
    return -20;

  return 0;
}
	  
