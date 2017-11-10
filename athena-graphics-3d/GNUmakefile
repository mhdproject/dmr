# NAG Fortran90 compiler and flags

export F90 = mpif90
##export F90FLAGS = -O0 -w -132  -r8 -DREAL8

#export F90 = mpif90.gm

#export F90FLAGS = -O4 -w -132 -dusty -float-store -r8 -DREAL8
#export F90FLAGS = -132 -w -w90 -w95 -r8 -DREAL8  -g
export F90FLAGS = -O3 -align -ip -unroll -132 -w -w90 -w95 -r8 -DREAL8 -g 
export F90FLAGS =-O3 -fast -Mextend -r8 -DREAL8 -g
export F90FLAGS = -extend_source -macro_expand  -r8 -DREAL8 -OPT:Olimit=7917 -O3 
export F90FLAGS =  -O3  -ffixed-line-length-132 -cpp -m32 -x f95-cpp-input

###export F90FLAGS = -O0  -tpp7 -axW -xW  -132 -w -w90 -w95 -r8 -DREAL8 -g 

##export F90FLAGS = -O0 -g90 -w -132 -C=array -C=do -C=present -C=pointer  -dusty -float-store -r8 -DREAL8 -gline


#SGI commands
#export F90 = f90
#export F90FLAGS = -extend_source -macro_expand -cpp -64 -r8 -DREAL8 -DNDEBUG -g -C -DEBUG:trap_uninitialized=ON
#export F90FLAGS = -extend_source -macro_expand -cpp -64 -r8 -DREAL8 -OPT:Olimit=7917 -O3 

# C compiler and flags
export CC = gcc
# linux
export CFLAGS = -O3  -I /hdf5-1.6.5/hdf5/include -DH5_USE_16_API -m32
#-I/usr/local/ChomboVis-4.1.3/usr/include
#-I${HOME}/hdf5/include 
#-I/usr/local/ChomboVis-3.94.3/usr/include  
#-I${HOME}/include \
#-I/usr/local/include 
# for gauss
#export CFLAGS = -O3 -64 -I/usr/local/src/Mesa-3.5/include

# linux
#export LDFLAGS = -lglut -lGLU -lGL
export LDFLAGS = -lm  \
-L/hdf5-1.6.5/hdf5/lib/ -lhdf5 -m32
#-L/usr/local/ChomboVis-4.1.3/usr/lib -lhdf5
#-L${HOME}/hdf5/lib -lhdf5

# NVIDIA driver for linux doesn't have off_screen routines so we need to use
# the Mesa lib that was moved to xxx.libGL.a.RPMSAVE
# that we then moved to libdefaultGL.a to make usable
# uncomment only if using the NVIDIA driver
#export LDFLAGS += -ldl -ldefaultGL

# for gauss
#export LDFLAGS = -64 -lmpi -L/usr/local/lib/Mesa-3.4.2 -lglut -lGLU -lGL \
#		 -lX11 -lXmu -lm

# top-level paramesh directory
#export PARAMESH = /home/sdo/UMA1/amr_jet/pulsar/paramesh3.1
##export PARAMESH = ${HOME}/cooling/paramesh3.1
export PARAMESH = ../../par
#export PARAMESH = ${HOME}/amr_jet/paramesh3.1
#export PARAMESH = /home/gauss_6/people/spicer/Fortran/paramesh
all:
	$(MAKE) -C application 

clean:
	$(RM) *~
	$(MAKE) -C application clean
