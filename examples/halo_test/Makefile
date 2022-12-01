#include ../../src/Makefile.inc

#INCLUDE = -I../../include
FFLAGS := $(subst $(MODFLAG),$(MODFLAG)../../,$(FFLAGS))
FFLAGS := $(patsubst -I%,-I../../%,$(FFLAGS))
LIBS = -L../../ -l$(LIBDECOMP) $(LFLAGS)

NP ?= 1
MPIRUN ?= mpirun

all: halo_test exchange_test

halo_test: halo_test.o
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -o $@ $^ $(LIBS)

exchange_test: exchange_test.o
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -o $@ $^ $(LIBS)

ifeq ($(PARAMOD),gpu)
check:
	$(MPIRUN) -n $(NP) ./bind.sh ./halo_test
else
check:
	$(MPIRUN) -n $(NP) ./halo_test
endif

clean:
	rm -f *.o halo_test exchange_test  *.log

%.o : %.f90
	$(FC) $(FFLAGS) $(OPT) $(DEFS) $(INC) -c $< -o $@
