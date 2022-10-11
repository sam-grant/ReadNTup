C++ = g++

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)

Target = Plotter.exe FastRotation.exe VerticalOffset.exe VerticalOffsetVsTime.exe OrthogonalPlotter.exe test.exe MomentumDist.exe Acceptance.exe PlotVerticalPosition.exe Count.exe ThetaYvsMomentum.exe CountEverything.exe

all:$(Target)

Plotter.exe: Plotter.C
	${C++} -o $@ Plotter.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

Count.exe: Count.C
	${C++} -o $@ Count.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

FastRotation.exe: FastRotation.C
	${C++} -o $@ FastRotation.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

VerticalOffset.exe: VerticalOffset.C
	${C++} -o $@ VerticalOffset.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

OrthogonalPlotter.exe: OrthogonalPlotter.C
	${C++} -o $@ OrthogonalPlotter.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

test.exe: test.C
	${C++} -o $@ test.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

MomentumDist.exe: MomentumDist.C
	${C++} -o $@ MomentumDist.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

Acceptance.exe: Acceptance.C
	${C++} -o $@ Acceptance.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

PlotVerticalPosition.exe: PlotVerticalPosition.C
	${C++} -o $@ PlotVerticalPosition.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

VerticalOffsetVsTime.exe: VerticalOffsetVsTime.C
	${C++} -o $@ VerticalOffsetVsTime.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

ThetaYvsMomentum.exe: ThetaYvsMomentum.C 
	${C++} -o $@ ThetaYvsMomentum.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

CountEverything.exe: CountEverything.C
	${C++} -o $@ CountEverything.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

clean:
	rm *exe
