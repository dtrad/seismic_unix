# ----------------------- PROGRAM NAME:

program = simptst

#------------------------ OBJECT FILES:

srcs = $(program).o nfft_routs.o fft_gen.o precision.o

# ----------------------- DESTINATION:

dest = .

# ----------------------- COMPILER:

fc = f90

# ----------------------- COMPILER OPTIONS:

opts = -O 

# ----------------------- LIBRARIES:

libs = -lsunperf

# ----------------------- THE MAKING OF THE PROGRAM ITSELF:

$(program) : $(srcs) 
	$(fc) $(opts) $(srcs) $(libs) -o $(dest)/$(program)

precision.o : precision.f90
	$(fc) $(opts) -c precision.f90
nfft_routs.o : nfft_routs.f90 precision.o fft_gen.o
	$(fc) $(opts) -c nfft_routs.f90
fft_gen.o : fft_gen.f90 precision.o
	$(fc) $(opts) -c fft_gen.f90

simptst.o : simptst.f90 nfft_routs.o
	$(fc) $(opts) -c simptst.f90

