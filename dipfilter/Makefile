SHELL= /bin/sh
OPTC= -O -Wall -fopenmp
IDIR=-I$(CWPROOT)/include  

include $(CWPROOT)/src/Makefile.config

D = $L/libcwp.a $L/libpar.a $L/libsu.a 

#$L/Complex.a

LFLAGS= $(PRELFLAGS) -L$L -lsu -lpar -lcwp -lm -lgomp $(POSTLFLAGS)

PROG=sudipfilt0

FLAG=-c -O -Wall

OBJ=  $(PROG).o 

LIBC=./clibrary

# dependence 
$(PROG): $(OBJ)
	gcc -o $(PROG) $(OBJ) $(CFLAGS) $(LFLAGS) -lgomp
$(PROG).o: $(PROG).c
	gcc -c $(OPTC) $(ENDIANFLAG) $(PROG).c $(IDIR) -fopenmp
clean:
	rm -f *.o $(PROG)




