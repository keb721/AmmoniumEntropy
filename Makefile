#-*- mode: makefile; mode: font-lock; vc-back-end: RCS -*-

SHELL = /bin/sh

# Compiler and settings

F90       = gfortran
LD        = gfortran
FFLAGS    = -O3 -g 
INCLUDE   =

.PRECIOUS: %.o
.PHONY:  clean

%: %.o
%.o: %.f90
	$(F90) $(FFLAGS) $(INCLUDE) -c -o $@ $<

%: %.o
	$(F90) $(FFLAGS) $(INCLUDE)  -o $@ $^

all :  config_objects refine_objects

CONFIG_OBJECTS = random.o ion_ghost.o run_ion.o
REFINE_OBJECTS = refine_bonds.o test_gro.o

config_objects : $(CONFIG_OBJECTS)
	$(LD) -o make_configs $(CONFIG_OBJECTS)

refine_objects : $(REFINE_OBJECTS)
	$(LD) -o make_gro $(REFINE_OBJECTS)


clean :
	rm -f *.mod *.d *.il *.o work.* make_configs make_gro













