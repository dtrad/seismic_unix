% Program to plot all Tau_p programs hr_ss, hr_vs, hr_parab, (from Mauricio Sacchi).
% Daniel Trad - UBC - Canada. 6/11/98
nmo=1; %0, 1 nmo.  
% Original data set
nt=512;
nh=40;
load cdp3.txt %data
data=reshape(cdp3,nt,nh);
load cdp3.h   %offset
h0=cdp3;

[v,av,bv]=read_hb('vel_gather.out');
[d,ad,bd]=read_hb('rec_data.out');

if nmo==1
[dnmo,adnmo,bdnmo]=read_hb('datanmo.out');
dnmo=reshape(dnmo,adnmo(2),adnmo(1));
h=0:ad(1)-1;
h=h*bd(1)-bd(3);
t=0:ad(2)-1;
t=t*bd(2)-bd(4);
figure,wigb(dnmo,1,h,t);title('NMO corrected data'),
xlabel('offset (m)');ylabel('time (s)');
end

h=0:ad(1)-1;
h=h*bd(1)-bd(3);
t=0:ad(2)-1;
t=t*bd(2)-bd(4);

q=0:av(1)-1;
q=q*bv(1)+bv(3);
tau=0:av(2)-1;
tau=tau*bv(2)-bv(4);

figure,

subplot(221),wigb(data,1,h,t)
title('Data');xlabel('p');ylabel('tau')

v=reshape(v,av(2),av(1));
subplot(222),wigb(v,1,q,tau)
title('Velocity Gather');xlabel('p');ylabel('tau')

max_data=max(max(abs(data)));

d=reshape(d,ad(2),ad(1));

max_d=max(max(abs(d)));
d=d.*max_data/max_d;

subplot(223),wigb(d,1,h,t)
title('Data recovered');xlabel('offset');ylabel('t')

r=data-d;
subplot(224),wigb(r,max_data,h,t)
title('Residuals');xlabel('offset');ylabel('t')

