sources := \
paramesh_dimensions.F \
paramesh_interfaces.F \
physicaldata.F \
prolong_arrays.F \
timings.F \
tree.F \
workspace.F \
io.F \
constants.F

ifndef SHMEM
sources += mpi_morton.F paramesh_mpi_interfaces.F
endif

objects := $(sources:.F=.o)

libmodules.a: $(objects)
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
	$(RM) libmodules.a *.o *.mod *.d *~ GNUmakefile.include
