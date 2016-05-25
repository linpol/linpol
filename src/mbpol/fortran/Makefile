FC = ifort
#FFLAGS = -O2 -openmp -lstdc++
FFLAGS = -O2 -lstdc++

all: test-mbpol

test-mbpol: test-mbpol.f ../mbpol/libmbpol.a
	$(FC) -o test-mbpol test-mbpol.f ../mbpol/libmbpol.a $(FFLAGS)

clean:
	rm -f test-mbpol

.PHONY: all clean
