F90FLAGS += -I$(PARAMESH)/headers 

modules := physcons.F $(wildcard module_*.F)
sources := $(filter-out $(modules), $(wildcard *.F))
objects := $(modules:.F=.o) $(sources:.F=.o)
headers := $(wildcard *.fh)

all: athena athena_graphics

athena: $(objects)  CHOMBO.o chomboHDF5.o
	$(F90) -o $@ $^ -L$(PARAMESH)/libs -lmpi_paramesh -lparamesh \
	  -lmodules $(LDFLAGS)

athena_graphics: athena_graphics.o 
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS) 

$(objects): $(headers)

.F.o:
	$(F90) $(F90FLAGS) -c $<

.PHONY: clean
clean:
	$(RM) core* *.o *.mod *~
