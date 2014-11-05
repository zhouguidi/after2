prefix=/usr/local/bin
FC = ifort
#FC = gfortran# -std=f2003
FFLAGS = -standard-semantics -std03
CDIINC=/home/gzhou/Fakeroot/usr/local/include
CDILIB=/home/gzhou/Fakeroot/usr/local/lib
LIBS = -lcdi

after2: after2.o cdiio.o config.o util.o
	$(FC) -o $@ $^ -L$(CDILIB) $(LIBS)
	
after2.o: after2.f90 cdiio.o config.o util.o
	$(FC) -c $< $(FFLAGS) -I$(CDIINC)

cdiio.o: cdiio.f90 config.o util.o
	$(FC) -c $< $(FFLAGS) -I$(CDIINC)

config.o: config.f90 util.o
	$(FC) -c $< $(FFLAGS) -I$(CDIINC)

util.o: util.f90
	$(FC) -c $< $(FFLAGS) -I$(CDIINC)

install:
	cp after2 after2.vardef $(prefix)

clean:
	rm *.o *.mod
