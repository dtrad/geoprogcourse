SHELL= /bin/sh
include $(CWPROOT)/src/Makefile.config
IDIR=-I$(CWPROOT)/include -I.


D = $L/libcwp.a $L/libpar.a $L/libsu.a 
LFLAGS= $(PRELFLAGS) -L$L -lsu -lpar -lcwp -lm $(POSTLFLAGS)

PROG=sumpitest

OPTC+= -g -O0
FLAG= -c -g -Wall -O0
OBJ=  $(PROG).o 

# dependence 
$(PROG): $(OBJ)
	mpicc -o $(PROG) $(OBJ) $(CFLAGS) $(LFLAGS) 
$(PROG).o: $(PROG).c
	mpicc -c $(OPTC) $(ENDIANFLAG) $(PROG).c $(FLAG) $(IDIR)
clean:
	rm -rf $(PROG) *.o 
