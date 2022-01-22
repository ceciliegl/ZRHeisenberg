LIBS   = -lm -lstdc++
LIBDIR = -L$(LD_LIBRARY_PATH)
#INCDIR = -I$(LD_LIBRARY_PATH) -I$(EBROOTEIGEN)/include
INCDIR = -I$(EBROOTEIGEN)/include

# C++ compiler
CCC= $(CXX)
CCOPTS= $(CFLAGS) -g -O3 -std=c++11 -Wno-deprecated-register
LINKOPTS =  -static -Wall -O3
CCFLAGS = $(CCOPTS) $(INCDIR)


RUNSOURCES =  run.cpp Makefile Solver.hpp indexstate.hpp BasicFunctions.hpp Quantities.hpp PrintFunctions.hpp ReadInputFiles.hpp


program: run.o
	$(CCC) $(LINKOPTS) -o $@ $^ $(LIBDIR) $(LIBS)

run.o:   $(RUNSOURCES)
				 $(CCC) $(CCFLAGS) -c -o $@ $<



spotless:
	make clean

clean   :
	rm -f *.dvi
	rm -f *.aux
	rm -f *~
	rm -f core
	rm -f *.o
	rm -f *.exec
	rm -f program
