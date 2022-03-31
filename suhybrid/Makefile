SHELL= /bin/sh

CWPROOT=$(HOME)/cwp
IDIR=-I$(CWPROOT)/include  -I./

include $(CWPROOT)/src/Makefile.config

D = $L/libcwp.a $L/libpar.a $L/libsu.a

LFLAGS= $(PRELFLAGS) -L$L -lsu -lpar -lcwp -lm   $(POSTLFLAGS)

PROG=suradonhybrid

FLAG=-c -Wall -Wno-write-strings -g -Wno-deprecated -Wno-deprecated-declarations 
OPTC= -O3 

OBJ= $(PROG).o radonhybrid.o radon_wtcgls_2op.o radoninv_2op.o radon_param.o \
		radon_param_2op.o interpovv.o modelweight.o nmo.o \
		radon_matrix.o interval.o \
		taper.o wtcgls_solvers.o getvelocities.o save_vector.o cleansegy.o \
		xtimesy.o xequaly.o Atimesx.o dot.o cdot.o fft.o fft2.o freqweight.o window.o \
		Complex.o rcdot.o plotgather.o filters.o
LIBC=./
CPP= g++

# dependence 
../bin/$(PROG): $(OBJ) $D
	${CPP} -o $(PROG) $(OBJ) $(CFLAGS) $(LFLAGS) 
$(PROG).o: $(PROG).cpp
	${CPP} $(FLAG) $(OPTC) $(ENDIANFLAG) $? $(IDIR)
radonhybrid.o: $(LIBC)/radonhybrid.cpp
	${CPP} $(FLAG) $(OPTC) $(ENDIANFLAG) $? $(IDIR)
radon_wtcgls_2op.o: $(LIBC)/radon_wtcgls_2op.cpp
	${CPP} $(FLAG) $(OPTC) $(ENDIANFLAG) $? $(IDIR)
radoninv_2op.o: $(LIBC)/radoninv_2op.cpp
	${CPP} $(FLAG) $(OPTC) $(ENDIANFLAG) $? $(IDIR)
radon_param.o: $(LIBC)/radon_param.cpp
	${CPP} $(FLAG) $(OPTC) $(ENDIANFLAG) $? $(IDIR)
radon_param_2op.o: $(LIBC)/radon_param_2op.cpp
	${CPP} $(FLAG) $(OPTC) $(ENDIANFLAG) $? $(IDIR)
interpovv.o: $(LIBC)/interpovv.cpp
	${CPP} $(FLAG) $(OPTC) $(ENDIANFLAG) $? $(IDIR)
modelweight.o: $(LIBC)/modelweight.cpp
	${CPP} $(FLAG) $(OPTC) $(ENDIANFLAG) $? $(IDIR)
nmo.o: $(LIBC)/nmo.cpp
	${CPP} $(FLAG) $(OPTC) $(ENDIANFLAG) $? $(IDIR)
radon_matrix.o: $(LIBC)/radon_matrix.cpp
	${CPP} $(FLAG) $(OPTC) $(ENDIANFLAG) $? $(IDIR)
interval.o: $(LIBC)/interval.cpp
	${CPP} $(FLAG) $(OPTC) $(ENDIANFLAG) $? $(IDIR)
taper.o: $(LIBC)/taper.c
	${CPP} $(FLAG) $(OPTC) $(ENDIANFLAG) $? $(IDIR)
wtcgls_solvers.o: $(LIBC)/wtcgls_solvers.cpp
	${CPP} $(FLAG) $(OPTC) $(ENDIANFLAG) $? $(IDIR)
getvelocities.o: $(LIBC)/getvelocities.cpp
	${CPP} $(FLAG) $(OPTC) $(ENDIANFLAG) $? $(IDIR)
save_vector.o: $(LIBC)/save_vector.cpp
	${CPP} $(FLAG) $(OPTC) $(ENDIANFLAG) $? $(IDIR)
cleansegy.o: $(LIBC)/cleansegy.c
	${CPP} $(FLAG) $(OPTC) $(ENDIANFLAG) $? $(IDIR)
xtimesy.o: $(LIBC)/xtimesy.cpp
	${CPP} $(FLAG) $(OPTC) $(ENDIANFLAG) $? $(IDIR)
xequaly.o: $(LIBC)/xequaly.cpp
	${CPP} $(FLAG) $(OPTC) $(ENDIANFLAG) $? $(IDIR)
cdot.o: $(LIBC)/cdot.cpp
	${CPP} $(FLAG) $(OPTC) $(ENDIANFLAG) $? $(IDIR)
Atimesx.o: $(LIBC)/Atimesx.cpp
	${CPP} $(FLAG) $(OPTC) $(ENDIANFLAG) $? $(IDIR)
dot.o: $(LIBC)/dot.cpp
	${CPP} $(FLAG) $(OPTC) $(ENDIANFLAG) $? $(IDIR)
rcdot.o: $(LIBC)/rcdot.cpp
	${CPP} $(FLAG) $(OPTC) $(ENDIANFLAG) $? $(IDIR)
fft.o: $(LIBC)/fft.cpp
	${CPP} $(FLAG) $(OPTC) $(ENDIANFLAG) $? $(IDIR)
fft2.o: $(LIBC)/fft2.cpp
	${CPP} $(FLAG) $(OPTC) $(ENDIANFLAG) $? $(IDIR)
window.o: $(LIBC)/window.cpp
	${CPP} $(FLAG) $(OPTC) $(ENDIANFLAG) $? $(IDIR)
freqweight.o: $(LIBC)/freqweight.cpp
	${CPP} $(FLAG) $(OPTC) $(ENDIANFLAG) $? $(IDIR)
Complex.o: $(LIBC)/Complex.cpp
	${CPP} $(FLAG) $(OPTC) $(ENDIANFLAG) $? $(IDIR)
plotgather.o: $(LIBC)/plotgather.cpp
	${CPP} $(FLAG) $(OPTC) $(ENDIANFLAG) $? $(IDIR)
filters.o: $(LIBC)/filters.cpp
	${CPP} $(FLAG) $(OPTC) $(ENDIANFLAG) $? $(IDIR)
clean::
	rm -f *.o *~ p1 p2 pp file* 	







