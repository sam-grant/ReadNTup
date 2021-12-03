C++ = g++

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)

Target = Plotter.exe FastRotation.exe VerticalOffset.exe OrthogonalPlotter.exe test.exe

all:$(Target)

Plotter.exe: Plotter.C
	${C++} -o $@ Plotter.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

FastRotation.exe: FastRotation.C
	${C++} -o $@ FastRotation.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

VerticalOffset.exe: VerticalOffset.C
	${C++} -o $@ VerticalOffset.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

OrthogonalPlotter.exe: OrthogonalPlotter.C
	${C++} -o $@ OrthogonalPlotter.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

test.exe: test.C
	${C++} -o $@ test.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)


clean:
	rm *exe
