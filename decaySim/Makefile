C++ = g++

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)

Target = Plotter.exe Acceptance.exe Count.exe SanityPlots.exe ThetaYvsMomentum.exe 

all:$(Target)

Plotter.exe: Plotter.C
	${C++} -o $@ Plotter.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

Acceptance.exe: Acceptance.C
	${C++} -o $@ Acceptance.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

Count.exe: Count.C
	${C++} -o $@ Count.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

SanityPlots.exe: SanityPlots.C
	${C++} -o $@ SanityPlots.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)
	
ThetaYvsMomentum.exe: ThetaYvsMomentum.C
	${C++} -o $@ ThetaYvsMomentum.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

clean:
	rm *exe
