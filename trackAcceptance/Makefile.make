C++ = g++

ROOTFLAGS = $(shell root-config --cflags)
ROOTLIBS  = $(shell root-config --libs)

Target = CheckFiles.exe

all:$(Target)

CheckFiles.exe: CheckFiles.C
	${C++} -o $@ CheckFiles.C $(Objects) $(ROOTFLAGS) $(ROOTLIBS)

clean:
	rm *exe
