CXX = icpc
CXXFLAGS = -O2
CPPFLAGS = -I. -I../mbpol
LDFLAGS = -O2

all: test-mbpol

test-mbpol: test-mbpol.o io-xyz.o xyz-water-utils.o
	$(CXX) -o test-mbpol test-mbpol.o io-xyz.o xyz-water-utils.o \
                  ../mbpol/libmbpol.a $(LDFLAGS)

clean:
	rm -f test-mbpol *.o

.PHONY: all clean
