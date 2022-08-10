C++ = g++

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)

Target = Plotter.exe Acceptance.exe Count.exe PlotVerticalPosition.exe VerticalOffset.exe MomentumDist.exe UniqueTracks.exe SanityPlots.exe

all:$(Target)

Plotter.exe: Plotter.C
	${C++} -o $@ Plotter.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

MomentumDist.exe: MomentumDist.C
	${C++} -o $@ MomentumDist.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

Acceptance.exe: Acceptance.C
	${C++} -o $@ Acceptance.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

Count.exe: Count.C
	${C++} -o $@ Count.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

PlotVerticalPosition.exe: PlotVerticalPosition.C
	${C++} -o $@ PlotVerticalPosition.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

VerticalOffset.exe: VerticalOffset.C
	${C++} -o $@ VerticalOffset.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

UniqueTracks.exe: UniqueTracks.C
	${C++} -o $@ UniqueTracks.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

SanityPlots.exe: SanityPlots.C
	${C++} -o $@ SanityPlots.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

clean:
	rm *exe
