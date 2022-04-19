# building Shared Library
SHELL= /bin/sh
OPTC= -O3 -Wall

IDIR=-I$(CWPROOT)/include  -I=. -I$(CWPROOT)/src/Complex/include

include $(CWPROOT)/src/Makefile.config

OPTC+=-fPIC
OPTCPP=-std=c++0x -g -O0 -Wno-comment -Wno-write-strings -fPIC

D = $L/libcwp.a $L/libpar.a $L/libsu.a 

LFLAGS= -L$L -lsu -lpar -lcwp -lm  

PROG=libgeoprog.so


FLAG=-c -O -Wall 

OBJ= dot.o conv.o Dcomplex.o Complex.o 

CPP= g++
CC= gcc

# dependence 
$(PROG): $(OBJ)
	$(CPP) -shared -fPIC -o $(PROG) $(OBJ) 
%.o:  %.cpp
	${CPP} $(FLAG) $(OPTCPP) $(ENDIANFLAG) $(IDIR) -c $<
%.o:  %.c
	${CC} $(FLAG) $(OPTC) $(ENDIANFLAG) $(IDIR) -c $<

clean:
	rm -rf $(PROG) *.o *~ 











