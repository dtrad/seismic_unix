function [scaleb]=normalize_energy(a,b)
% Given two data sets, a and b, it computes the scale factor to use in wigb.
% The first array, a, must be plotted with maximum 1.
% The second array b, must be plotted with maximum scaleb,
% Daniel Trad, UBC. 1999
a=normalize(a);
b=normalize(b);
ena=(sum(sum(a.^2))).^0.5;
enb=(sum(sum(b.^2))).^0.5;
[ma,na]=size(a);
[mb,nb]=size(b);
ena=ena/ma/na;
enb=enb/mb/nb;
scaleb=ena/enb;

