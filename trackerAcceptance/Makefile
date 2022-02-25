C++ = g++

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)

Target = Plotter.exe #PlotterTest.exe

all:$(Target)

Plotter.exe: Plotter.C
	${C++} -o $@ Plotter.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

PlotterTest.exe: PlotterTest.C
	${C++} -o $@ PlotterTest.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

clean:
	rm *exe