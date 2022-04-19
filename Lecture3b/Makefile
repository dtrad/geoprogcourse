#Using shared library
CWPROOT=$(HOME)/cwp
include $(CWPROOT)/src/Makefile.config
NAME=susynth

CPP=g++
D = $L/libcwp.a $L/libpar.a $L/libsu.a  
IDIR=-I$(CWPROOT)/include -I.  
LD_LIBRARY_PATH=$(HOME)/geoprog/Lecture4b

LFLAGS= -L$L  -lsu -lpar -lcwp -lm 
LFLAGS+= -L$(LD_LIBRARY_PATH) -lgeoprog  

OPTC=-std=c++0x -g -O0 -Wno-comment -Wno-write-strings -Wno-pedantic

OBJECTS = $(NAME).o
# dependence 
$(NAME): $(OBJECTS)
	$(CPP) $(OBJECTS) $(CFLAGS) $(LFLAGS) -o $@
%.o: %.cpp
	$(CPP) $(OPTC) $(ENDIANFLAG) $(IDIR) -c $<
clean:
	rm -rf $(NAME) *.o *~ core p.su p2.su




