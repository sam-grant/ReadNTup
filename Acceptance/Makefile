C++ = g++

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)

Target = Plotter.out # PlotterMisalignment.exe #PlotterTest.exe

all:$(Target)

Plotter.out: Plotter.C
	${C++} -o $@ Plotter.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

# PlotterMisalignment.exe: PlotterMisalignment.C
# 	${C++} -o $@ PlotterMisalignment.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

# PlotterTest.exe: PlotterTest.C
# 	${C++} -o $@ PlotterTest.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

clean:
	rm *.out