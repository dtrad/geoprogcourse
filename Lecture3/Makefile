CWPROOT=$(HOME)/cwp
include $(CWPROOT)/src/Makefile.config
NAME=susynth

CPP=g++
D = $L/libcwp.a $L/libpar.a $L/libsu.a 
IDIR=-I$(CWPROOT)/include -I.  

LFLAGS= -L$L -lsu -lpar -lcwp -lm 
OPTC-=std
OPTC+=-std=c++0x 
OPTC+=-g -O0 -Wno-comment -Wno-write-strings -Wno-pedantic

OBJECTS = $(NAME).o Complex.o
# dependence 
$(NAME): $(OBJECTS)
	$(CPP) $(OBJECTS) $(CFLAGS) $(LFLAGS) -o $@
%.o: %.cpp
	$(CPP) $(OPTC) $(ENDIANFLAG) $(IDIR) -c $<
clean:
	rm -rf $(NAME) *.o *~ core p.su p2.su




