.PHONY: all
all: bin/r2

CC=gcc
CXX=g++

CXXFLAGS= -std=c++11 -O2 $(VERSION)
CFLAGS = -O2 $(VERSION)

LFLAGS = -lz -lm -lpthread

# spdlog
SPDLOGDIR=external/spdlog
IFLAGS += -Iexternal/ -I$(SPDLOGDIR) 

#htslib stuff
HTSDIR=external/htslib-1.12
include $(HTSDIR)/htslib.mk
HTSLIB = $(HTSDIR)/libhts.a
IFLAGS += -I$(HTSDIR)

#special builds for debugging and profiling
debug: CXXFLAGS = -std=c++11 -g -O0 -Wall $(VERSION)
debug: CFLAGS = -g -O0 $(VERSION)
debug: all

profile: CXXFLAGS = -std=c++11 -pg $(VERSION)
profile: CFLAGS =  -pg $(VERSION)
profile: all

OBJS=CorrelationCalculator.o
%.o: %.cpp %.hh
	$(CXX) $(CXXFLAGS) $(IFLAGS) $(TESTFLAGS) -c -o $@ $<

bin/r2: r2.cpp $(OBJS) $(HTSLIB)
	mkdir -p bin/
	$(CXX) $(CXXFLAGS) -o $@ $< $(OBJS) $(HTSLIB) $(IFLAGS) $(LFLAGS)
clean:
	rm bin/r2 *.o
	test/observed*.txt
