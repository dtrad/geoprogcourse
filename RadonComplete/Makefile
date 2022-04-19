SHELL= /bin/sh
OPTC= -O3 -Wall
CWPROOT=/home/dtrad/cwp
IDIR=-I$(CWPROOT)/include -I.

include $(CWPROOT)/src/Makefile.config

D = $L/libcwp.a $L/libpar.a $L/libsu.a $L/libsudan.a

LFLAGS= $(PRELFLAGS) -L$L -lsu -lpar -lcwp -lm -lsudan $(POSTLFLAGS)

PROG=suradonfreq

FLAG= -Wno-deprecated -Wno-deprecated-declarations -Wno-write-strings -g -O3

OBJECTS = $(PROG).o radonsolver0_subtract.o radonsolver.o modelweight.o \
		 radon_matrix.o nmo.o interpovv.o getvelocities.o radoninv.o \
	       	 wtcgls_solvers.o radon_cgfft2.o radon_toeplitz.o radon_cgtoep.o \
		 Mute.o



LIBC=./clibrary
CPP= g++

$(PROG): $(OBJECTS)
	g++ $(OBJECTS) $(IDIR) $(D) $(LFLAGS) -o $@
%.o: %.cpp 
	g++  $(IDIR) $(FLAG) -c $<
%.o: %.c
	g++  $(IDIR) $(FLAG) -c $<
clean: 
	rm -rf $(PROG) *.o *.~ *.d primaries.su dataout.su model.su landcdpin.su landcdpout.su data.demul.su

#need to add dependency list.
#include $(subst .o,.d,$(OBJECTS))
#%.d: %.cpp
#	$(CPP) -MM -std=c++0x $< > $@

# dependence 
#../bin/$(PROG): $(OBJ) $D
#	${CPP} -o ../bin/$(PROG) $(OBJ) $(CFLAGS) $(LFLAGS) 
#$(PROG).o: $(PROG).cpp
#	${CPP} $(FLAG) $(OPTC) $(ENDIANFLAG) $? $(IDIR)








