#**********************************************************************
# Author: John Walker
# Email:  jmgwalker@triumf.ca
# Date:   2019/09/25
#********************************************************************** 
INCDIR= $(shell pwd)/include
SRCDIR= $(shell pwd)/src
OBJDIR= $(shell pwd)/obj
BINDIR= $(shell pwd)/bin

VPATH = $(SRCDIR)

CFLAGS=-c -g -Wall `root-config --cflags` -I${INCDIR}
LDFLAGS=`root-config --glibs` -lHistPainter

TARGET1=field_to_csv.cpp
TARGET2=ptf_analysis.cpp
TARGET3=ptf_ttree_analysis.cpp
TARGET4=ptf_qe_analysis.cpp
TARGET5=ptf_field_analysis.cpp
TARGET6=ptf_charge_analysis.cpp
TARGET7=ptf_timing_analysis.cpp

EXECUTABLE1=$(TARGET1:%.cpp=$(BINDIR)/%.app)
EXECUTABLE2=$(TARGET2:%.cpp=$(BINDIR)/%.app)
EXECUTABLE3=$(TARGET3:%.cpp=$(BINDIR)/%.app)
EXECUTABLE4=$(TARGET4:%.cpp=$(BINDIR)/%.app)
EXECUTABLE5=$(TARGET5:%.cpp=$(BINDIR)/%.app)
EXECUTABLE6=$(TARGET6:%.cpp=$(BINDIR)/%.app)
EXECUTABLE7=$(TARGET7:%.cpp=$(BINDIR)/%.app)

FILES= $(wildcard $(SRCDIR)/*.cpp)
SOURCES=$(FILES)

OBJECTS = $(FILES:$(SRCDIR)/%.cpp=${OBJDIR}/%.o)

OBJ1=$(TARGET1:%.cpp=${OBJDIR}/%.o) $(OBJECTS)
OBJ2=$(TARGET2:%.cpp=${OBJDIR}/%.o) $(OBJECTS)
OBJ3=$(TARGET3:%.cpp=${OBJDIR}/%.o) $(OBJECTS)
OBJ4=$(TARGET4:%.cpp=${OBJDIR}/%.o) $(OBJECTS)
OBJ5=$(TARGET5:%.cpp=${OBJDIR}/%.o) $(OBJECTS)
OBJ6=$(TARGET6:%.cpp=${OBJDIR}/%.o) $(OBJECTS)
OBJ7=$(TARGET7:%.cpp=${OBJDIR}/%.o) $(OBJECTS)

all: MESSAGE $(EXECUTABLE1) $(EXECUTABLE2) $(EXECUTABLE3) $(EXECUTABLE4) $(EXECUTABLE5) $(EXECUTABLE6) $(EXECUTABLE7)

MESSAGE:
	@echo '**********************************************************************'
	@echo '* Compiling ptf-analysis-2 programs:                                 *'
	@echo '*   - field_to_csv                                                   *'
	@echo '*   - ptf_analysis                                                   *'
	@echo '*   - ptf_ttree_analysis                                             *'
	@echo '*   - ptf_qe_analysis                                                *'
	@echo '*   - ptf_field_analysis                                             *'
	@echo '*   - ptf_charge_analysis                                            *'
	@echo '*   - ptf_timing_analysis                                            *'
	@echo '**********************************************************************'

$(EXECUTABLE1): $(OBJECTS) $(OBJ1)
	$(CXX) $^ -o $@ $(LDFLAGS)

$(EXECUTABLE2): $(OBJECTS) $(OBJ2)
	$(CXX) $^ -o $@ $(LDFLAGS)

$(EXECUTABLE3): $(OBJECTS) $(OBJ3)
	$(CXX) $^ -o $@ $(LDFLAGS)

$(EXECUTABLE4): $(OBJECTS) $(OBJ4)
	$(CXX) $^ -o $@ $(LDFLAGS)

$(EXECUTABLE5): $(OBJECTS) $(OBJ5)
	$(CXX) $^ -o $@ $(LDFLAGS)

$(EXECUTABLE6): $(OBJECTS) $(OBJ6)
	$(CXX) $^ -o $@ $(LDFLAGS)

$(EXECUTABLE7): $(OBJECTS) $(OBJ7)
	$(CXX) $^ -o $@ $(LDFLAGS)

$(OBJDIR)/%.o: %.cpp
	$(CXX) $(CFLAGS) $< -o $@

clean:
	- $(RM) $(BINDIR)/* $(OBJDIR)/*
