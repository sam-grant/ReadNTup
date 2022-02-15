C++ = g++

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)

Target = Plotter.exe OrthogonalPlotter.exe Acceptance.exe

all:$(Target)

Plotter.exe: Plotter.C
	${C++} -o $@ Plotter.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

OrthogonalPlotter.exe: OrthogonalPlotter.C
	${C++} -o $@ OrthogonalPlotter.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

Acceptance.exe: Acceptance.C
	${C++} -o $@ Acceptance.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)
	
clean:
	rm *exe
