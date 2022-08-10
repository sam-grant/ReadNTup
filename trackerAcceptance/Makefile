C++ = g++

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)

<<<<<<< HEAD
Target = Plotter.exe PlotterMisalignment.exe #PlotterTest.exe
=======
Target = Plotter.exe #PlotterTest.exe
>>>>>>> 61e4c6dfc811565a0487218eb445720b32e0a372

all:$(Target)

Plotter.exe: Plotter.C
	${C++} -o $@ Plotter.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

<<<<<<< HEAD
PlotterMisalignment.exe: PlotterMisalignment.C
	${C++} -o $@ PlotterMisalignment.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)
=======
PlotterTest.exe: PlotterTest.C
	${C++} -o $@ PlotterTest.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)
>>>>>>> 61e4c6dfc811565a0487218eb445720b32e0a372

clean:
	rm *exe