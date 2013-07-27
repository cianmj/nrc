goal: run.x bilrun.x nrgmin.x testp.x

F90 = ifc -FR -w95 -static -r16 
# -O0 -g -p -C -d16
# -O3 (remove all after -O0)
EXTRALIB = -L /opt/intel/mkl/lib/32 -lmkl_lapack -lmkl_p3 -lguide -lpthread
OBJECTS = modules.o Potential.o DPoten.o Hamdiag.o Hzeros.o hermite.o dapprx.o linslv.o AngPot.o norm.o
PROPGTN = propgtn.f calcrate.o calcquad.o calcsing.o SVD.o condcheck.o evalwf.o
TESTP = testp.f oldPot.o
BILLING = billing.f evalwf.o calcbil.o
NRGMIN = nrgmin.f

nrgmin.x : $(OBJECTS) $(NRGMIN)
	$(F90) $(OBJECTS) $(NRGMIN) $(EXTRALIB) -o nrgmin.x

bilrun.x : $(OBJECTS) $(BILLING)
	$(F90) $(OBJECTS) $(BILLING) $(EXTRALIB) -o bilrun.x

run.x : $(OBJECTS) $(PROPGTN)
	$(F90) $(OBJECTS) $(PROPGTN) $(EXTRALIB) -o run.x

testp.x : $(OBJECTS) $(TESTP)
	$(F90) $(OBJECTS) $(TESTP) $(EXTRALIB) -o testp.x

calcbil.o : calcbil.f
	$(F90) -c calcbil.f -o calcbil.o

norm.o : norm.f
	$(F90) -c norm.f -o norm.o

evalwf.o : evalwf.f
	$(F90) -c evalwf.f -o evalwf.o

condcheck.o : condcheck.f
	$(F90) -c condcheck.f -o condcheck.o

oldPot.o : oldPot.f
	$(F90) -c oldPot.f -o oldPot.o

SVD.o : SVD.f
	$(F90) -c SVD.f -o SVD.o

AngPot.o : AngPot.f
	$(F90) -c AngPot.f -o AngPot.o

modules.o : modules.f
	$(F90) -c modules.f -o modules.o

Potential.o : Potential.f
	$(F90) -c Potential.f -o Potential.o

DPoten.o : DPoten.f
	$(F90) -c DPoten.f -o DPoten.o

Hamdiag.o : Hamdiag.f
	$(F90) -c Hamdiag.f -o Hamdiag.o

linslv.o : linslv.f
	$(F90) -c linslv.f -o linslv.o

hermite.o : hermite.f
	$(F90) -c hermite.f -o hermite.o

Hzeros.o : Hzeros.f
	$(F90) -c Hzeros.f -o Hzeros.o

dapprx.o : dapprx.f
	$(F90) -c dapprx.f -o dapprx.o

calcrate.o : calcrate.f
	$(F90) -c calcrate.f -o calcrate.o

calcquad.o : calcquad.f
	$(F90) -c calcquad.f -o calcquad.o

calcsing.o : calcsing.f
	$(F90) -c calcsing.f -o calcsing.o

clean:
	-rm *.x *.mod *.out
