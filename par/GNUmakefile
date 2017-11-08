#--------------------------------------------------------------------------
#
# This Gmake file will compile the PARAMESH library and create a
# set of library files to which you can link. To use it, make sure
# it is in the PARAMESH root directory.  
# It works by running gmake on the Makefile.gnu files which appear
# in the headers, source and mpi_source sub-directories.
# To simply create these PARAMESH library files, type
#     gmake -f Makefile.gnu
# when in the PARAMESH root directory. The library files will
# then be found in a newly created sub-directory called libs.
#
# If you type 
#     gmake -f Makefile.gnu Tests
# it will also compile and link the test programs in the Tests
# sub-directory. There is a file called Makefile.gnu inside Tests
# which is used.
# 
# To compile and link application files in a sub-directory called
# User_applic you could type
#     gmake -f Makefile.gnu User_applic
# provided you copy Makefile.gnu from Tests to User_applic, and modify
# it appropriately.
#
#
# Written : Ernest Mamikonyan        April 2002.
#
#--------------------------------------------------------------------------


# Define the fortran compiler
###export FC = mpif90.gm

export FC = mpif90

#export FC = f90

#-----------------------------------------------
 
# Set the desired compilation flags
# SGIs
#export FFLAGS = -cpp -O3 -r8 -DREAL8 -I../headers 

# NAG95 linux - you must use -float-store. 
#export FFLAGS = -O2 -r8 -DREAL8 -dusty -w -I../headers
#export FFLAGS = -O2 -r8 -I../headers
#export FFLAGS = -g -C -r8 -I../headers
#export FFLAGS = -O2 -float-store -nan -r8  -I../headers
#export FFLAGS = -O4 -float-store -r8  -I../headers
#export FFLAGS = -O4 -float-store -r8 -C=all -gline -I../headers

#-----------------------------------------------

# Additional libraries to link to. You do not need
# to add the shmem library. This is automatically added
# if you define SHMEM=1 below.
# export ADD_LIB = /usr/lib32/libmpi.so

# FOR ifc
#export FC = ifc
#export FFLAGS = -O3 -ip -i4 -r8 -I../headers -I/usr/local/mpich-intel/include
export FFLAGS = -O3  -unroll -align -ip -i4 -r8 -I../headers  -w -w90 -w95 -c -m32
export FFLAGS = -O3 -i4 -r8 -I../headers 
export FFLAGS = -O3 -I../headers -m32
#export ADD_LIB = -L/usr/local/mpich-intel/lib -lmpichf90 -lmpich -lPEPCF90

#-----------------------------------------------

# some compilers can generate make rules to stdout from the source files
# if you have such a compiler, provide the flags, otherwise comment it out
#export MY_CPP := gcc -E -MM -MG  # for the GNU C Preprocessor

#-----------------------------------------------

# SHMEM or MPI ?
# uncomment to use SHMEM
#export SHMEM = 1

#--------------------------------------------------------------------------


.PHONY: all
ifdef SHMEM
all: libs headers source
else
all: libs headers mpi_source source
endif

.PHONY: headers
headers:
	$(MAKE) -C $@ -f Makefile.gnu
	cp -f headers/libmodules.a libs

.PHONY: mpi_source
mpi_source: headers
	$(MAKE) -C $@ -f Makefile.gnu
	cp -f mpi_source/libmpi_paramesh.a libs

.PHONY: source
source: headers
	$(MAKE) -C $@ -f Makefile.gnu
	cp -f source/libparamesh.a libs

.PHONY: clean
clean:
	$(RM) -r *~ libs
	for dir in headers {mpi_,}source Tests amr_kernel_bmark; do \
	  $(MAKE) -C $$dir -f Makefile.gnu clean; \
	done

.PHONY: Tests
Tests: all
	$(MAKE) -C $@ -f Makefile.gnu

libs:
	mkdir $@

# An example target to match an application directory name other than Tests
# in which the users application files are located.
amr_kernel_bmark:	all
	$(MAKE) -C $@ -f Makefile.gnu
