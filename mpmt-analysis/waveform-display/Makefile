#**********************************************************************
# Author: John Walker
# Email:  jmgwalker@triumf.ca
# Date:   2019/09/25
#********************************************************************** 
INCDIR= $(shell pwd)/../../include
SRCDIR= $(shell pwd)/../../src
OBJDIR= $(shell pwd)/../../obj
BINDIR= $(shell pwd)/../../bin

VPATH = $(SRCDIR)

CFLAGS=-c -g -Wall `root-config --cflags` -I${INCDIR}
LDFLAGS=`root-config --glibs` -lHistPainter -lMinuit -L${ROOTSYS}/lib

all: waveform_plotting.exe

waveform_plotting.exe: waveform_plotting.o 
	CPATH=/usr/local/include $(CXX) $^ -o $@ $(LDFLAGS)

waveform_plotting.o: waveform_plotting.cpp
	$(CXX) $(CFLAGS) $< -o $@

wrapper.o: waveform_plotting
	$(CXX) $(CFLAGS) $< -o $@

clean:
	- $(RM) *.exe *.o