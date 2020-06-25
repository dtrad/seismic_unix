load cdps.dat
fid1 = fopen('tnmopar','r');
fid2 = fopen('vnmopar','r');

fid3 = fopen('parfile','w');

%fid4 = fopen('cdps','w');
%for i=1:length(cdps)
%  fprintf(fid4,'%d,',cdps(i));
%end
%fprintf(fid4,'\n');
%fclose(fid4)

fid4 = fopen('cdps','r');

S0=fgets(fid4);
fprintf(fid3,'cdps=',S0);
for i=1:45
 S1 = fgets(fid1);
 S2 = fgets(fid2);
 fprintf(fid3,'%s%s',S1,S2);
end 
fclose(fid1);
fclose(fid2);
fclose(fid3);
fclose(fid4);