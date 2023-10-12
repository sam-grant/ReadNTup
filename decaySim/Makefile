C++ = g++

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)

Target = Plotter.out SanityPlotter.out

all:$(Target)

Plotter.out: Plotter.C
	${C++} -o $@ Plotter.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

SanityPlotter.out: SanityPlotter.C
	${C++} -o $@ SanityPlotter.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)
	
clean:
	rm *.out
