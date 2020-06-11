

all: ptf_bfield.exe

ptf_bfield.exe: ptf_bfield.o
	g++ -std=c++11 -Wall -O2 ptf_bfield.o `root-config --libs` -o ptf_bfield.exe

ptf_bfield.o: ptf_bfield.cpp t2kstyle.h
	g++ -std=c++11 -Wall -O2 -I`root-config --incdir` -c ptf_bfield.cpp


clean:
	rm *.o *.exe

.PHONY: all clean
