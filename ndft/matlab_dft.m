% matlab_dft.m  test program dft. See simptst

% define same data as in simptst:
dsig=[1+i -1+i .8+.3i 10-i 2+2i -2+2i .4+.6i 20-2i].';
ind_irr=[ 1.1, 1.9, 2.9, 4., 5.1, 6.3, 6.9, 8.1];

M=8;
Xmax=2*pi/.78;

%----- DFT matrix -----
N=length(dsig);

m=(-(M/2):(-1+M/2))';
delta_kx=2*pi/Xmax;

Airr=exp(i*m*ind_irr*delta_kx);

%----- Actual DFT -----
fsig=Airr*dsig


