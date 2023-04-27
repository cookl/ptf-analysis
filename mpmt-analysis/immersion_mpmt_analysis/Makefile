#**********************************************************************
# Author: Hassan Elshabasy
# Email:  helshabasy@triumf.ca
# Date:   2023/01/16
#********************************************************************** 
INCDIR= $(shell pwd)/../../include
SRCDIR= $(shell pwd)/../../src
OBJDIR= $(shell pwd)/../../obj
BINDIR= $(shell pwd)/../../bin

VPATH = $(SRCDIR)

CFLAGS=-c -g -Wall `root-config --cflags` -I${INCDIR}
LDFLAGS=`root-config --glibs` -lHistPainter -lMinuit -L${ROOTSYS}/lib 

all: mpmt_immersion_analysis.exe mpmt_time_calibration.exe

mpmt_immersion_analysis.exe:  mpmt_immersion_analysis.o WaveformFitResult.o
	CPATH=/usr/local/include $(CXX) $^ -o $@ $(LDFLAGS)

mpmt_immersion_analysis.o: mpmt_immersion_analysis.cpp
	$(CXX) $(CFLAGS) $< -o $@
	
mpmt_time_calibration.exe:  mpmt_time_calibration.o WaveformFitResult.o
	CPATH=/usr/local/include $(CXX) $^ -o $@ $(LDFLAGS)

mpmt_time_calibration.o: mpmt_time_calibration.cpp
	$(CXX) $(CFLAGS) $< -o $@


WaveformFitResult.o: ${SRCDIR}/WaveformFitResult.cpp
	$(CXX) $(CFLAGS) $< -o $@

clean:
	- $(RM) *.exe *.o
