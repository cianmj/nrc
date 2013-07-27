goal: ghwp1D.x

.SUFFIXES: .f90 .o

F90 = ifc -FR -w95 -static -O0 -g -p -C
# -O3 (remove all after -O0)
EXTRALIB = -L /opt/intel/mkl/lib/32 -lmkl_lapack -lmkl_p3 -lguide -lpthread

OBJECTS = ghwp_options.o potential_xxx.o ghwp_math.o ghwp_pe.o lapack.o \
         ghwp_properties.o main.o ghwp.o

ghwp1D.x : $(OBJECTS)
	$(F90) $(OBJECTS) $(EXTRALIB) -o ghwp1D.x

main.o: main.f90 ghwp.o ghwp_options.o ghwp_pe.o potential_xxx.o

ghwp.o: ghwp.f90 ghwp_options.o ghwp_pe.o lapack.o

ghwp_pe.o: ghwp_pe.f90 ghwp_options.o ghwp_math.o potential_xxx.o

ghwp_properties.o : ghwp_properties.f90 ghwp_options.o ghwp_math.o

ghwp_math.o: ghwp_math.f90 ghwp_options.o

potential_xxx.o: potential_xxx.f90 ghwp_options.o

ghwp_options.o: ghwp_options.f90

lapack.o : lapack.f90

.f90.o:
	$(F90) -c $<

clean:
	-rm *.x *.mod *.out *.o
