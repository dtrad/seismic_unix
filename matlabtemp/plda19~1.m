eq.-1  freq-offset
c         stx(nt,nx)      - if index.eq. 1  time-offset
c         nfft            - lenght of the transform
c
c       Notes:
c
c        The input/output changes according to index:
c
c          index = -1 TX ----> FX
c          index =  1 FX ----> TX

        complex * 16    sfx(nfft,nx), aux(nfft)
        real    * 8     stx(nfft,nx)

        if(index.eq.-1)  then

        do 140 ix=1,nx
		do 130 it=1,nt
	        aux(it)=dcmplx(stx(it,ix),0.d0)
130     continue
      );
end
[d,ad,bd]=read_hb('rec_data.out');

if nmo==1
[dnmo,adnmo,bdnmo]=read_hb('datanmo.out');
dnmo=reshape(dnmo,adnmo(2),adnmo(1));
t=0:ad(2)-1;
t=t*bd(2)+bd(4);
figure,wigb(dnmo,1,h0,t);title('NMO corrected data'),
xlabel('offset (m)');ylabel('time (s)');
end

h=0:ad(1)-1;
h=h*bd(1)+bd(3);
t=0:ad(2)-1;
t=t*bd(2)-bd(4);

q=0:av(1)-1;
q=q*bv(1)+bv(3);
tau=0:av(2)-1;
tau=tau*bv(2)-bv(4);


figure,
data=normalize(data);
scaled=normalize_energy(data,d);
subplot(221),wigb(data,1,h0,t)
title('Data');xlabel('offset');ylabel('t (sec)')

v=reshape(v,av(2),av(1));
subplot(222),wigb(v,1,q,tau)
title('Velocity Gather');xlabel('q');ylabel('\tau')

max_data=max(max(abs(data)));

d=reshape(d,ad(2),ad(1));
max_d=max(max(abs(d)));

d=normalize(d);
if (ss==1)
	subplot(223),wigb(d,scaled,h0,t)
elseif ss~=1
   subplot(223),wigb(d,scaled,h,t)
end
title('Data recovered');xlabel('offset');ylabel('t')
if (mute==1)
   subplot(224),wigb(v2,1,q,tau);
   title('Muted Radon domain');xlabel('p');ylabel('\tau')
end

if (residuals==1)
r=data-d;
subplot(224),wigb(r,max_data,h,t)
title('Residuals');xlabel('offset');ylabel('t')
end
if (ss==1)
	figure,wigb(d,scaled,h0,t)
elseif ss~=1
   figure,wigb(d,scaled,h,t)
end

title('Data recovered');xlabel('offset (m)');ylabel('t (s)')
