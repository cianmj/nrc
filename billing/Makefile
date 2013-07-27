goal: bilrun.x

F90 = ifc -g -FR -C -w95 -static -p -O0 -r16 -d16
EXTRALIB = -L /opt/intel/mkl/lib/32 -lmkl_lapack -lmkl_p3 -lguide -lpthread
BILLING = modules.o billing.f dapprx.o DPoten.o Potential.o Hamdiag.o Hzeros.o hermite.o AngPot.o

bilrun.x : $(BILLING)
	$(F90) $(BILLING) $(EXTRALIB) -o bilrun.x

dapprx.o : dapprx.f
	$(F90) -c dapprx.f -o dapprx.o

modules.o : modules.f
	$(F90) -c modules.f -o modules.o

Potential.o : Potential.f
	$(F90) -c Potential.f -o Potential.o

AngPot.o : AngPot.f
	$(F90) -c AngPot.f -o AngPot.o

DPoten.o : DPoten.f
	$(F90) -c DPoten.f -o DPoten.o

Hamdiag.o : Hamdiag.f
	$(F90) -c Hamdiag.f -o Hamdiag.o

hermite.o : hermite.f
	$(F90) -c hermite.f -o hermite.o

Hzeros.o : Hzeros.f
	$(F90) -c Hzeros.f -o Hzeros.o

clean:
	-rm *.x *.mod *.out
