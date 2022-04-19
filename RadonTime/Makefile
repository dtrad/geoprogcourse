SHELL= /bin/sh
OPTC= -O3 -Wall
CWPROOT=/home/dtrad/cwp
IDIR=-I$(CWPROOT)/include -I.

include $(CWPROOT)/src/Makefile.config

D = $L/libcwp.a $L/libpar.a $L/libsu.a

LFLAGS= $(PRELFLAGS) -L$L -lsu -lpar -lcwp -lm $(POSTLFLAGS)

PROG=suradontime

FLAG= -Wno-deprecated -Wno-deprecated-declarations -Wno-write-strings -g -O0

OBJECTS = $(PROG).o  radontd_win.o wpcgnr2.o mytaper.o contran.o getvelocities.o \
		     interpovv.o smoothing.o read_ascii_file.o get_wavelet.o     \
	             subtract_multiples.o smute.o modelweight.o dot.o save_vector.o \
		     plotgather.o

#LIBC=./clibrary
CPP= g++

$(PROG): $(OBJECTS)
	g++ $(OBJECTS) $(IDIR) $(D) $(LFLAGS) -o $@
%.o: %.cpp 
	g++  $(IDIR) $(FLAG) -c $<
%.o: %.c
	g++  $(IDIR) $(FLAG) -c $<
clean: 
	rm -rf $(PROG) *.o *~ *.d

#need to add dependency list.
#include $(subst .o,.d,$(OBJECTS))
#%.d: %.cpp
#	$(CPP) -MM -std=c++0x $< > $@

# dependence 
#../bin/$(PROG): $(OBJ) $D
#	${CPP} -o ../bin/$(PROG) $(OBJ) $(CFLAGS) $(LFLAGS) 
#$(PROG).o: $(PROG).cpp
#	${CPP} $(FLAG) $(OPTC) $(ENDIANFLAG) $? $(IDIR)








