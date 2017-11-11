/*
**   _______              __
**  / ___/ /  ___  __ _  / /  ___
** / /__/ _ \/ _ \/  ' \/ _ \/ _ \
** \___/_//_/\___/_/_/_/_.__/\___/ 
**
** chomboHDF5.h
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

/*  Serial programming interface */


#include <hdf5.h>
#include <string.h>

struct {
  int i;
  int j;
}typedef intvect2d ;

struct {
  int i;
  int j;
  int k;
}typedef intvect3d ;

struct {
  intvect2d lo;
  intvect2d hi;
}typedef box2d ;

struct{
  intvect3d lo;
  intvect3d hi;
} typedef box3d;

union 
{
  box2d b2; 
  box3d b3;
} typedef box;

int numPnts2(const box2d*);
int numPnts3(const box3d*);
void refine2(box2d*, int refinement);
void refine3(box3d*, int refinement);

static hid_t intvect2d_id, intvect3d_id, box2d_id, box3d_id ;

void initializeHDF5datatypes();  /* gives values to *_id objects */



/* information pointed at by 'data' in these classes
   is assumed to be packed like a fortran data array:

   int c, i, j, sizei, sizej, index;
   float d;
   FloatArray2D arrayBox;
   .
   .
   sizei = arrayBox.hi.i-arrayBox.lo.i + 1;
   sizej = arrayBox.hi.j-arrayBox.lo.j + 1;
   for(c=0; c<arrayBox.ncomps; ++c){
     for(j=0; j<sizej; ++j){
       for(i=0; i<sizei; ++i){
         index = c*(sizei*sizej) + j*sizei + i;
	 d = arrayBox.data[index];
       }
     }
   }

   */
struct
{
  float* data;
  int    ncomps;
  box2d  box;
}typedef  FloatArrayBox2D;

struct
{
  double* data;
  int     ncomps;
  box2d   box;
}typedef DoubleArrayBox2D ;

struct 
{
  float* data;
  int    ncomps;
  box3d  box;
}typedef FloatArrayBox3D ;

struct
{
  double* data;
  int     ncomps;
  box3d   box;
}typedef DoubleArrayBox3D ;

union
{
  FloatArrayBox2D f2;
  DoubleArrayBox2D d2;
  FloatArrayBox3D f3;
  DoubleArrayBox3D d3;
} typedef FAB;

void freeDataArray(FAB* ArrayBoxCArray,
		   int length); /* send in {Float|Double}ArrayBox{2|3}D* */

enum{Float, Double} typedef  ChPrecision;

struct
{
  hid_t file_ID;
  hid_t group_ID;
  int   dim;
  ChPrecision precision;
}typedef HDF5Handle;

/* three next functions return 0 on success */

/* returns negative number if file open failed
           1  if 'Chombo_global' couldn't be found and ChDim and ChPrecision not set
   accessMode is one of H5F_ACC_RDWR, H5F_ACC_RDONLY, H5F_ACC_CREATE
	   */
int Handleopen(HDF5Handle* handle, const char* filename, hid_t accessMode);
int HandlesetGroup(HDF5Handle* handle, const char* group);/* ie. "/level_0" */
int Handleclose(HDF5Handle handle);
int HandleCreateGroup(HDF5Handle* handle, const char* group);/* ie. "/level_0" */


#define ChTYPES 8
enum{INTEGER, FLOAT, DOUBLE, CHAR,  INTVECT2D, INTVECT3D, BOX2D, BOX3D}typedef ChattributeType;
static const char* TYPE_NAMES[ChTYPES];

struct
{
  void* data;
  char* name;
  int attributeType;
  void* next;
} typedef Attribute;

/*  struct that holds all the attributes of an HDF5 group.  

  each entry in the accessByType array is the head of a linked
  list of members terminated with a null.

  to access directly by type you can use code like this:

   HDF5attributes a;
   readHDF5attributes(&a, handle);
   int numChar = a.numByType(CHAR);
   Attribute* at = a.accessByType[CHAR];
   while(at != NULL)
   {
      printf("char attribute %s = %s",
      at->name,
      (char*)(at->data));
      at = at->next;
   }
   */
struct
{
  int numByType[ChTYPES];
  Attribute*  accessByType[ChTYPES]; 
  int   dim;
  ChPrecision precision;
} typedef HDF5attributes;

void freeHDF5attributes(HDF5attributes* attributes);
void printAttributes(HDF5attributes* attributes);


/* utility functions for appending attributes to an HDF5attributes struct
   attrib:  HDF5atributes struct */
Attribute* ADD_tmp;
Attribute* ADD_at;
#define ADD_ATTRIBUTE(attrib, Ttype, ChType, Sname, value)         \
ADD_tmp = attrib.accessByType[ChType] ;                            \
ADD_at = (Attribute*) malloc(sizeof(Attribute));                   \
attrib.accessByType[ChType] = ADD_at;                              \
ADD_at->name = malloc(strlen(Sname));                              \
sprintf(ADD_at->name, Sname);                                      \
ADD_at->data = malloc(sizeof(Ttype));                              \
ADD_at->attributeType = ChType;                                    \
*((Ttype *)ADD_at->data) = value;                                  \
  ADD_at->next = ADD_tmp;                                          \
attrib.numByType[ChType]++;                                        \

#define ADD_CHAR_ATTRIBUTE(attrib,  Sname, value)                  \
ADD_tmp = attrib.accessByType[CHAR] ;                              \
ADD_at = (Attribute*) malloc(sizeof(Attribute));                   \
attrib.accessByType[CHAR] = ADD_at;                                \
ADD_at->name = malloc(strlen(Sname));                              \
sprintf(ADD_at->name, Sname);                                      \
ADD_at->data = malloc(strlen(value));                              \
ADD_at->attributeType = CHAR;                                      \
sprintf((char*)ADD_at->data, value);                               \
  ADD_at->next = ADD_tmp;                                          \
attrib.numByType[CHAR]++;                                          \

/* next eight functions return 0 on success */

/* - users responsibility to make sure HDF5Handle is pointing at a valid group
   in an open file. 
   - User must use freeHDF5attributes when finished.
   - User must pass in a pointer to a valid HDF5attributes struct.*/
int readHDF5attributes(HDF5attributes* attr, HDF5Handle handle);

/* input: HDF5Handle
   users responsibility to cast *(boxes) to either box2d* 
   or box3d* and free it */
int readBoxes(box** boxes, int* length, HDF5Handle handle);

/* input: HDF5Handle, length, boxes
   users responsibility to cast *(arrayboxes) to one of
   {Float|Double}ArrayBox{2|3}D*, and to call freeDataArray and 
   free(*arrayboxes) when finished 
  
   returns negative on HDF5 error.
  
*/
int readData(FAB** arrayboxes, box* boxes, int length, HDF5Handle handle);

/* - users responsibility to make sure HDF5Handle is pointing at a valid group
   in an open file. 
   - User must use freeHDF5attributes when finished.
   - User must pass in a pointer to a valid HDF5attributes struct.

   returns 0 on success

   -1 attribute creation problems  */
int writeHDF5attributes(HDF5attributes* attr, HDF5Handle handle);

int writeBoxes(box* boxes,    int length, HDF5Handle handle);
int writeFABBoxes(FAB* FABArray, int length, HDF5Handle handle);

int writeData(FAB* FABArray, int length, HDF5Handle handle);


/* Top level C functions, simplest interface  */
 
/* everything is input.  

   returns 0 on success.

  refinementRatios :  refinement ratio at all levels
  (ith entry is refinement ratio between levels i and i + 1)

  returns

  -1: error opening file
  -2: error writing attributes to header
  -4: error writing boxes to data file
  -5: error writing FAB data to file
  -20: error closing file on exit

*/

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
		      int* refinementRatios);

/* input : filename

   all pointers must point to a valid allocated object whose value will be filled in

   User is responsible for calling 'freeDataArray' on all the 'levels' 
   data when finished with
   them and calling 'free' on all data members passed out of this routine.

   refinementRatios :  refinement ratio at all levels
   (ith entry is refinement ratio between levels i and i + 1)

   returns 0 on success

   -1: unable to open file
   -2: Unable to read attributes from group "/"
   -3: Number of levels not found in HDF5 file
   -4: Number of components not found in HDF5 filen
   -5: num_levels <= 0
   -6: num_components <= 0
   -7: missing component name
   -8: error reading boxes
   -20: error closing file on exit;
   -100 + readData error: error reading data.
   */
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
		     int** refinementRatios);
		      
