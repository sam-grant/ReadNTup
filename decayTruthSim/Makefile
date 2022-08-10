C++ = g++

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)

Target = Plotter.exe Acceptance.exe Count.exe

all:$(Target)

Plotter.exe: Plotter.C
	${C++} -o $@ Plotter.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

Acceptance.exe: Acceptance.C
	${C++} -o $@ Acceptance.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

Count.exe: Count.C
	${C++} -o $@ Count.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

clean:
	rm *exe
