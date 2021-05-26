LIBDIR=
LIBFLAGS=
WHICHSYS := $(shell uname -n)
ifeq ($(WHICHSYS), damic.uchicago.edu)
ROOTLIBDIR=-L/usr/local/lib64
ROOTINCDIR=-I/usr/local/include/root
else
ifeq ($(WHICHSYS), opus.uchicago.edu)
ROOTLIBDIR=-L/usr/local/lib64
ROOTINCDIR=-I/usr/local/include/root
LIBDIR=-L/usr/local/gcc620/lib64/
LIBFLAGS=-Wl,-rpath=/usr/local/gcc620/lib64
else
ifeq ($(WHICHSYS), zev.uchicago.edu)
ROOTLIBDIR=-L/usr/local/lib64
ROOTINCDIR=-I/usr/local/include/root
else
ROOTLIBDIR=-L/usr/lib/x86_64-linux-gnu 
ROOTINCDIR=-I/usr/include/root
endif
endif
endif
ROOTLIBS=-lCore -lHist -lRIO
#CC=g++ 

#LD=G++
CFLAGS=-t -g -O0 -std=c++0x -D_FILE_OFFSET_BITS=64 -fmax-errors=5 -Wall
LIBS=-lcfitsio
LIBS+=$(ROOTLIBDIR) $(ROOTLIBS)


SRCDIR=./
SRC=$(wildcard $(SRCDIR)/*.cpp)
INCDIR=./
INCLUDE=-Iinclude -I$(INCDIR) $(ROOTINCDIR)
OBJS=$(SRC:.cpp=.o)
MAIN=process_data

all: $(MAIN) 


debug: CFLAGS += -DDEBUG 
debug: clean all



$(MAIN): $(OBJS)
	$(CXX) $(OBJS) $(INCLUDE) $(CFLAGS) $(LIBDIR) $(LIBS) $(LIBFLAGS) -o $@

.cpp.o:
	$(CXX) $(INCLUDE) $(CFLAGS) -c $< -o $@

clean: 
	rm -f $(OBJS) $(MAIN)


# DO NOT DELETE
