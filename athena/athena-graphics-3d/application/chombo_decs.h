  /* CHOMBO declarations */

  char compTmp[80];
  HDF5Handle handle;
  HDF5attributes attr;
  int r;
  herr_t ret;
  hid_t boxdataspace[100], boxdataset[100];
  hsize_t* FABoffsets;
  hsize_t flatdims[1];
  herr_t err;
  hsize_t count[1];
  hssize_t offset[1];
  hid_t attributeGroup, compattr, aid;
  hid_t dataspace[100], dataset[100], memdataspace, datatype;
  int comps;

#define Float Float
  ChPrecision precision = Double;

