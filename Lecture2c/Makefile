CWPROOT=$(HOME)/cwp
include $(CWPROOT)/src/Makefile.config
NAME=susynth

CPP=gcc
D = $L/libcwp.a $L/libpar.a $L/libsu.a 
IDIR=-I$(CWPROOT)/include  


LFLAGS= -L$L -lsu -lpar -lcwp -lm 
OPTC-=std
OPTC+=-std=c++0x 
OPTC+=-g -O0 -Wno-comment -Wno-write-strings

OBJECTS = $(NAME).o 


# dependence 
$(NAME): $(OBJECTS)
	$(CPP) $(OBJECTS) $(CFLAGS) $(LFLAGS) -o $@
%.o: %.c
	$(CPP) $(OPTC) $(ENDIANFLAG) $(IDIR) -c $<
clean:
	rm -rf $(NAME) *.o *~ core p.su p2.su




