load c:\daniel\mauricio\data.dat
[v,av,bv]=read_hb2('c:\daniel\mauricio\hr_parab\debug\vel_gather.out');
[d,ad,bd]=read_hb2('c:\daniel\mauricio\hr_parab\debug\rec_data.out');




figure,

subplot(221),wigb(data)
title('Data');xlabel('p');ylabel('tau')

v=reshape(v,av(2),av(1));
subplot(222),wigb(v)
title('Velocity Gather');xlabel('p');ylabel('tau')

max_data=max(max(abs(data)));

d=reshape(d,ad(2),ad(1));

max_d=max(max(abs(d)));
d=d.*max_data/max_d;

subplot(223),wigb(d)
title('Data recovered');xlabel('offset');ylabel('t')


r=data-d;
subplot(224),wigb(r,max_data)
title('Residuals');xlabel('offset');ylabel('t')
