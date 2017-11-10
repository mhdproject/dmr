FROM 32bit/ubuntu:16.04

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gcc \
    gcc-multilib \
    gfortran \
    gzip \
    libopenmpi-dev \
    make \
    ssh \
    tar \
    vim \
    wget

#RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.8 10
#RUN ln -s /usr/bin/cpp-4.8 /usr/bin/cpp
#RUN wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.6/hdf5-1.6.5.tar.gz

COPY hdf5-1.6.5.tar.gz /

# RUN wget https://www.open-mpi.org/nightly/v3.0.x/openmpi-v3.0.x-201711040323-888fac7.tar.gz
#RUN tar xzvf openmpi-v3.0.x-201711040323-888fac7.tar.gz
#WORKDIR /openmpi-v3.0.x-201711040323-888fac7
#RUN ./configure  CFLAGS=-m32 CXXFLAGS=-m32 FFLAGS=-m32 FCFLAGS=-m32
#RUN make all install 
WORKDIR /
RUN tar xzvf hdf5-1.6.5.tar.gz

WORKDIR hdf5-1.6.5

RUN ./configure CFLAGS="-w -msse2 -m32" LDFLAGS="-m32"
RUN make
RUN make install

COPY par.tgz /
WORKDIR /
RUN tar xzvf par.tgz 
WORKDIR /par
RUN make


COPY ath.tgz /
WORKDIR /
RUN tar xzvf ath.tgz 
WORKDIR /athena-graphics-3d/application
RUN mpif90 -c -I../../par/headers -ffixed-line-length-132  -x f95-cpp-input physcons.F
RUN mpif90 -c -I../../par/headers -ffixed-line-length-132 -x f95-cpp-input module_grid.F 
RUN mpif90 -c -I../../par/headers -ffixed-line-length-132 -x f95-cpp-input module_hstrongfield.F
RUN mpif90 -c -I../../par/headers -ffixed-line-length-132 -x f95-cpp-input module_athena_specific_parameters.F 
RUN mpif90 -c -I../../par/headers -ffixed-line-length-132 -x f95-cpp-input module_interior_gbc.F
WORKDIR /athena-graphics-3d/
RUN make
WORKDIR /athena-graphics-3d/application
RUN apt-get install -y openmpi-bin



