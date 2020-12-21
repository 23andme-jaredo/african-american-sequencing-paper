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
HTSDIR=external/htslib-1.10.2
include $(HTSDIR)/htslib.mk
HTSLIB = $(HTSDIR)/libhts.a
IFLAGS += -I$(HTSDIR)

#special builds for debugging and profiling
debug: CXXFLAGS = -std=c++11 -g -O1 -Wall $(VERSION)
debug: CFLAGS = -g -O1 $(VERSION)
debug: all

profile: CXXFLAGS = -std=c++11 -pg $(VERSION)
profile: CFLAGS =  -pg $(VERSION)
profile: all

bin/r2: r2.cpp $(OBJS) $(HTSLIB)
	mkdir -p bin/
	$(CXX) $(CXXFLAGS) -o $@ $< $(HTSLIB) $(IFLAGS) $(LFLAGS)
clean:
	rm bin/r2
