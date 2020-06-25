load c:\daniel\mauricio\data.in
load c:\daniel\mauricio\hr_ss\debug\vel_gather.out
load c:\daniel\mauricio\hr_ss\debug\rec_data.out
figure,

subplot(221),wigb(data)
title('Data');xlabel('p');ylabel('tau')

p=reshape(vel_gather,512,80);
subplot(222),wigb(p)
title('Velocity Gather');xlabel('p');ylabel('tau')

max_data=max(max(abs(data)));

d=reshape(rec_data,512,64);
max_d=max(max(abs(d)));
d=d.*max_data/max_d;
subplot(223),wigb(d)
title('Data recovered');xlabel('offset');ylabel('t')


r=data-d;
subplot(224),wigb(r,max_data)
title('Residuals');xlabel('offset');ylabel('t')
