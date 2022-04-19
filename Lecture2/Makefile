SHELL= /bin/sh
include $(CWPROOT)/src/Makefile.config

IDIR=-I$(CWPROOT)/include  
D = $L/libcwp.a $L/libpar.a $L/libsu.a 

LFLAGS= -L$L -lsu -lpar -lcwp -lm 
PROG=supick

OPTC+=-g -O0

OBJ=  $(PROG).o 


# dependence 
$(PROG): $(OBJ)
	gcc -o $(PROG) $(OBJ) $(LFLAGS) 
$(PROG).o: $(PROG).c
	gcc -c $(PROG).c $(OPTC) $(CFLAGS) $(IDIR)
clean:
	rm -rf $(PROG) *.o *~ core list pp




