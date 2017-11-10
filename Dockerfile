FROM i686/ubuntu

RUN apt-get update && apt-get install -y \
    make \
    tar \
    build-essential \
    libopenmpi-dev \
    gzip \
    gfortran \
    gcc \
    gcc-multilib \
    wget \
    ssh \
    vim

#RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.8 10
#RUN ln -s /usr/bin/cpp-4.8 /usr/bin/cpp
COPY par.tgz /
COPY ath.tgz /
#RUN wget https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.6/hdf5-1.6.5.tar.gz

#COPY hdf5-1.6.5.tar.gz /

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

WORKDIR /
RUN tar xzvf par.tgz 
WORKDIR /par
RUN make


WORKDIR /
RUN tar xzvf ath.tgz 
WORKDIR /athena-graphics-3d/application
RUN mpif90 -c physcons.F -I../../par/headers -ffixed-line-length-132
RUN mpif90 -c module_grid.F -I../../par/headers -ffixed-line-length-132
RUN mpif90 -c module_hstrongfield.F -I../../par/headers
RUN mpif90 -c module_athena_specific_parameters.F -I../../par/headers -ffixed-line-length-132
RUN mpif90 -c module_interior_gbc.F -I../../par/headers -ffixed-line-length-132
WORKDIR /athena-graphics-3d/
#RUN make

