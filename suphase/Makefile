SHELL= /bin/sh


include $(CWPROOT)/src/Makefile.config
OPTC = -g -Wall -pedantic -Wno-long-long  -m64 -Wno-write-strings
D = $L/libcwp.a $L/libpar.a $L/libsu.a 
#OPTC -= -Wmultichar
IDIR=-I$(CWPROOT)/include  -I/. -I/data/home/dtrad/work/suphase


#$L/Complex.a
LFLAGS= $(PRELFLAGS) -L$L -lsu -lpar -lcwp -lm $(POSTLFLAGS) 

PROG=suphase

DEBUG= -g -O0 

#DEBUG=

OBJ=  $(PROG).o Complex.o

#LIBC=./clibrary

# dependence 
$(PROG): $(OBJ)
	g++ -o $(PROG) $(OBJ) $(CFLAGS) $(LFLAGS) 
$(PROG).o: $(PROG).cpp
	g++ -c $(OPTC) $(ENDIANFLAG) $(PROG).cpp $(IDIR) $(DEBUG)
Complex.o: Complex.cpp
	g++ -c $(OPTC) $(ENDIANFLAG) $? $(IDIR) $(DEBUG)
clean : 
	rm $(OBJ) $(PROG)



