SHELL= /bin/sh
CWPROOT=/home/dtrad/cwp
IDIR=-I$(CWPROOT)/include  -I../Lecture4 -I.

include $(CWPROOT)/src/Makefile.config

CC=gcc
CPP=g++

OPTC+= -g -O0 -Wno-comment -Wno-write-strings
OPTCPP=-g -O0 -Wno-comment -Wno-write-strings

D = $L/libcwp.a $L/libpar.a $L/libsu.a $L/libgeoprog.a

LFLAGS= -L$L -lsu -lpar -lcwp -lm -lgeoprog 

PROG=suadapdiff

FLAG=-c -O3 -Wall

OBJ=  $(PROG).o save_vector.o cleansegy.o smoothing.o fftgo.o fftback.o cdot.o

LIBC=.

# dependence 
$(PROG): $(OBJ)
	$(CPP) $(OBJ) $(CFLAGS) $(LFLAGS) -o $@ 
%.o: %.c
	$(CC) $(OPTC) $(ENDIANFLAG) $? $(IDIR) -c $<
%.o: %.cpp
	$(CPP) $(OPTCPP) $(ENDIANFLAG) $? $(IDIR) -c $<
clean:
	rm -rf $(PROG) *.o *~ *.su



