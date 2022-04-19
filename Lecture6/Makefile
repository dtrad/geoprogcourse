SHELL= /bin/sh
OPTC= -O3 -Wall
CWPROOT=$(HOME)/cwp
IDIR=-I$(CWPROOT)/include -I. -I../Lecture4/.

include $(CWPROOT)/src/Makefile.config

D = $L/libcwp.a $L/libpar.a $L/libsu.a $L/libgeoprog.a

LFLAGS= $(PRELFLAGS) -L$L -lsu -lpar -lcwp -lm -lgeoprog $(POSTLFLAGS)

PROG=suradonfreq

FLAG= -Wno-deprecated -Wno-deprecated-declarations -Wno-write-strings -g -O0

OBJECTS = $(PROG).o save_vector.o cleansegy.o radonclass.o fft.o radonsolver.o radon_toeplitz.o rcdot.o

CPP= g++
CC= gcc

$(PROG): $(OBJECTS)
	$(CPP) $(OBJECTS) $(IDIR) $(D) $(LFLAGS) -o $@
%.o: %.cpp 
	$(CPP)  $(IDIR) $(FLAG) -c $<
%.o: %.c
	$(CC)  $(IDIR) $(FLAG) -c $<
clean: 
	rm -rf $(PROG) *.o *~ data model.su p1 p2 p4 p pp










