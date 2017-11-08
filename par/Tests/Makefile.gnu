#LD_SPARSKIT = \
#  ../Unsupported/SPARSKIT2/UNSUPP/BLAS1/blas1.o \
#  ../Unsupported/distdot.o \
#  -L../Unsupported/SPARSKIT2 \
#  -lskit 

LD_SPARSKIT = 


LDFLAGS += -L../libs
ifndef SHMEM
LDFLAGS += -lmpi_paramesh
endif

# normal LDFLAGS
ifndef T3E
LDFLAGS += -lparamesh -lmodules
endif

LDFLAGS += $(LD_SPARSKIT)

ifndef SHMEM
LDFLAGS += -lmpi_paramesh
endif

ifdef T3E
# T3E LDFLAGS
LDFLAGS += -lparamesh
endif

ifdef SHMEM
ifdef SGI
# SGI LDFLAGS
LDFLAGS += -lsma
endif
ifdef COMPAQ
# Compaq LDFLAGS
LDFLAGS += -lshmem
endif
endif

LDFLAGS += $(ADD_LIB)

# List all application source files
sources := \
 amr_1blk_bcset.F \
 check_data.F \
 zero_guardcells.F 
# divb_test_supp.F

objects := $(sources:.F=.o)

# tests are assumed to be name test_*.F
# so this strips test_ and .F from all files that match test_*.F
# tests := $(patsubst %.F,%,$(patsubst test_%,%,$(wildcard test_*.F)))
tests := \
guardcell_1blk \
prolong_1blk \
c_to_f_1blk \
c_to_f_1blk_2 \
c_to_f_1blk_3 \
c_to_f_1blk_4 \
c_to_f_1blk_quad \
c_to_f_1blk_divb \
1blk_guardcell_icoord \
1blk_guardcell_nlayers \
1blk_guardcell_big \
prolong_multilevel_1blk \
checkpoint \
checkpoint1 \
restrict_1blk \
flux_conserve_1blk \
edges_1blk \
multi_level_1 \
multi_level_2 \
force_consist \
gcell_on


# compiles all tests as defined above
.PHONY: all
all: $(wildcard test_*.F) $(objects)
	for test in $(tests); do $(MAKE) -f Makefile.gnu $$test; done

# compiles one particular test
# this is the target that does the work for the one above
%: test_%.F $(objects)
	$(FC) $(FFLAGS) -o test_$@ $^ $(LDFLAGS)

$(objects): test_defs.fh

.PHONY: clean
clean:
	$(RM) *.o *.i *~
	for test in $(tests); do $(RM) $(addprefix test_,$$test); done
