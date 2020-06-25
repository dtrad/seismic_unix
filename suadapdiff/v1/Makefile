SHELL= /bin/sh
OPTC=

IDIR=-I$(CWPROOT)/include  -I.


include $(CWPROOT)/src/Makefile.config

D = $L/libcwp.a $L/libpar.a $L/libsu.a 

LFLAGS= $(PRELFLAGS) -L$L -lsu -lpar -lcwp -lm $(POSTLFLAGS)

PROG=suadapdiff

FLAG=-c -O3 -Wall

OBJ=  $(PROG).o dot.o save_vector.o cleansegy.o smoothing.o fftgo.o fftback.o cdot.o

LIBC=.

# dependence 
$(HOME)/bin/$(PROG): $(OBJ)
	g++ -o $(HOME)/bin/$(PROG) $(OBJ) $(CFLAGS) $(LFLAGS) 
$(PROG).o: $(PROG).cpp
	g++ -c $(OPTC) $(ENDIANFLAG) $? $(IDIR)
dot.o: $(LIBC)/dot.cpp
	g++ -c $(OPTC) $(ENDIANFLAG) $? $(IDIR)
save_vector.o: $(LIBC)/save_vector.cpp
	g++ -c $(OPTC) $(ENDIANFLAG) $? $(IDIR)
cleansegy.o: $(LIBC)/cleansegy.c
	g++ -c $(OPTC) $(ENDIANFLAG) $? $(IDIR)
smoothing.o: $(LIBC)/smoothing.cpp
	g++ -c $(OPTC) $(ENDIANFLAG) $? $(IDIR)
fftgo.o: $(LIBC)/fftgo.cpp
	g++ -c $(OPTC) $(ENDIANFLAG) $? $(IDIR)
fftback.o: $(LIBC)/fftback.cpp
	g++ -c $(OPTC) $(ENDIANFLAG) $? $(IDIR)
cdot.o: $(LIBC)/cdot.cpp
	g++ -c $(OPTC) $(ENDIANFLAG) $? $(IDIR)
clean:
	rm $(PROG) *~ *.o
