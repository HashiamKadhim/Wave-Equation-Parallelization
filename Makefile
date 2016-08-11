#
#  
# Modules required
#
#module load gcc/4.8.1
#module load openmpi/gcc/1.6.4
#module load pgplot/5.2.2-gcc


CXX=mpic++
CPPFLAGS=-I${SCINET_PGPLOT_INC}
CXXFLAGS=-g -O0 -std=c++11
LDLIBS=-lcpgplot -lpgplot -lX11 -lgfortran -lpng 
LDFLAGS=-L${SCINET_PGPLOT_LIB}

all: wave1d

# Makefile for wave equation

wave1d: wave1d.o ticktock.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

wave1d.o: wave1d.cc ticktock.h
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) -o $@ $<

ticktock.o: ticktock.cc ticktock.h
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) -o $@ $<

clean:
	rm -f *.o
