CWPROOT=$(HOME)/cwp
include $(CWPROOT)/src/Makefile.config
NAME=susynthdt

CPP=g++
D = $L/libcwp.a $L/libpar.a $L/libsu.a 
IDIR=-I$(CWPROOT)/include  


LFLAGS= -L$L -lsu -lpar -lcwp -lm 
OPTC-=std
OPTC+=-std=c++0x 
OPTC+=-g -O0 -Wno-comment -Wno-write-strings

OBJECTS = $(NAME).o Complex.o


# dependence 
$(NAME): $(OBJECTS)
	$(CPP) $(OBJECTS) $(CFLAGS) $(LFLAGS) -o $@
%.o: %.cpp
	$(CPP) $(OPTC) $(ENDIANFLAG) $(IDIR) -c $<
clean:
	rm -rf $(NAME) *.o *~ core




