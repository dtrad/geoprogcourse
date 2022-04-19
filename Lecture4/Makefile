#building static library
SHELL= /bin/sh
OPTC= -O3 -Wall

IDIR=-I$(CWPROOT)/include  -I=. -I$(CWPROOT)/src/Complex/include

include $(CWPROOT)/src/Makefile.config

OPTC=-std=c++0x -g -O0 -Wno-comment -Wno-write-strings

D = $L/libcwp.a $L/libpar.a $L/libsu.a 

LFLAGS= -L$L -lsu -lpar -lcwp -lm  

PROG=libgeoprog.a


FLAG=-c -O -Wall 

OBJ= dot.o conv.o Dcomplex.o Complex.o Atimesx.o 

CPP= g++

# dependence 
$L/$(PROG): $(OBJ)
	ar r  $L/$(PROG) $(OBJ) 
%.o:  %.cpp
	${CPP} $(FLAG) $(OPTC) $(ENDIANFLAG) $(IDIR) -c $<
clean:
	rm -rf $L/$(PROG) *.o *~ 











