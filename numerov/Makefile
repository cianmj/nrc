goal: he3-c60.x he4-c60.x h2-S.x h2-L.x d2-S.x hd-S.x xxx.x

# FLAGS = -C -g -O0
FLAGS = -O2 -tpp6 -ip -xK -mp

he3-c60.x:	he3-c60.f90 process.f90 numerov.f90
	ifc $(FLAGS) -o he3-c60.x numerov.f90 he3-c60.f90

he4-c60.x:	he4-c60.f90 process.f90 numerov.f90
	ifc $(FLAGS) -o he4-c60.x numerov.f90 he4-c60.f90

h2-S.x:	h2-S.f90 process.f90 numerov.f90
	ifc $(FLAGS) -o h2-S.x numerov.f90 h2-S.f90

d2-S.x:	d2-S.f90 process.f90 numerov.f90
	ifc $(FLAGS) -o d2-S.x numerov.f90 d2-S.f90

hd-S.x:	hd-S.f90 process.f90 numerov.f90
	ifc $(FLAGS) -o hd-S.x numerov.f90 hd-S.f90

h2-L.x:	h2-L.f90 process.f90 numerov.f90
	ifc $(FLAGS) -o h2-L.x numerov.f90 h2-L.f90

xxx.x:  xxx.f90 process.f90 numerov.f90
	ifc $(FLAGS) -o xxx.x numerov.f90 xxx.f90

to-office:
	rsync -urtv -e ssh . ps@demon.sims.nrc.ca:numerov/

from-office:
	rsync -urtv -e ssh ps@demon.sims.nrc.ca:numerov/ .

clean:
	-rm *.out *.x *.mod *.dat *.gpl
