suximage d2=5 < moddata2.rms.su
suximage d2=5 < moddata2.rms.su colorbar=1
suximage d2=5 < moddata2.rms.su bar=1
ximage
suximage d2=5 < moddata2.rms.su legend=1
bg
lt
less unisamb.900
sed < unisamb.800 's/cdp=800//' > pp; mv pp unisamb.800
less unisamb.800
emacs unisamb.900
emacs unisamb.800 
emacs unisamb.700 
emacs unisamb.600 
emacs unisamb.500 
emacs unisamb.400 
emacs unisamb.300 
emacs unisamb.200 
emacs unisamb.100 
Xplotvel 900
Xplotvel 800
Xplotvel 700
Xplotvel 600
Xplotvel 500
Xplotvel 400
Xplotvel 300
Xplotvel 200
Xplotvel 100
lt \*.ps
lt ../cseg/\*.ps
Xplotvelan
h | grep Velan0
h | grep Xsukmig1c
df
ls *.ps
h | head -10
h | tail -10
Velan0 moddata2.filt.su 0.3 100 900 100  
Velan0 moddata2.filt.su 0.6 100 900 100  
clean
clean
Xsukmig1c 100 900 2 41 10 0.02 20 1 100 
lt | head -2
cp moddata2.csp.su moddata2.csp.B.su
Xsukmig1c 100 900 4 41 10 0.02 20 1 100 
lt | head -2
cp moddata2.csp.su moddata2.csp.LSM.su
mv ../ps/fig3cseg.ps ../ps/fig3seg2000n.ps
Xplotvelan
ls m*.ps
Xplotvelan
h | grep Velan0
 Velan0 moddata2.csp.B.su 0.5 100 900 100 
clean
 Velan0 moddata2.csp.B.su 0.3 100 900 100 
clean
 Velan0 moddata2.csp.LSM.su 0.5 100 900 100 
clean
Xplotvelan
 Velan0 moddata2.csp.LSM.su 0.7 100 900 100 
clean
Xplotvelan
Xplotvelan
Xplotvelan
ls *.ps
mv ../cseg ../seg2000
ls *.ps | more
mv *.ps ../seg2000/.
Xplotvelan
 Velan0 moddata2.csp.B.su 0.3 100 900 100 
clean
mv *.ps ../seg2000/.
Xplotvelan
Xsukmig1c 100 900 2 41 10 0.02 20 1 150 
 Velan0 moddata2.csp.B.su 0.3 100 900 100 
clean
mv *.ps ../seg2000/.
Xplotvelan
 Velan0 moddata2.csp.B.su 0.3 100 900 100 
mv *.ps ../seg2000/.
Xplotvelan
clean; Velan0 moddata2.csp.LSM.su 0.7 100 900 100 
mv *.ps ../seg2000/.
Xplotvelan
clean; Velan0 moddata2.csp.LSM.su 0.8 100 900 100 
clean;Xplotvelan
Xsukmig1c 100 900 2 41 10 0.02 20 1 100 
h | grep cp
 cp moddata2.csp.su moddata2.csp.B.stkvel.su 
Xsukmig1c 100 900 4 41 10 0.02 20 1 150 
clean
 cp moddata2.csp.su moddata2.csp.LSM.stkvel.su 
clean; Velan0 moddata2.csp.LSM.stkvel.su 0.7 100 900 100 
clean; Velan0 moddata2.csp.LSM.stkvel.su 0.3 100 900 100 
clean
clean; Velan0 moddata2.csp.B.stkvel.su 0.1 100 900 100 
clean; Velan0 moddata2.csp.LSM.stkvel.su 0.1 100 900 100 
clean; Xplotvelan
mv *.ps ../seg2000/.
ls ../seg2000
