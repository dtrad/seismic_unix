function [p,a,b]=read_hb(name)
fid=fopen(name,'rt');
a=fread(fid,2,'int32');
b=fread(fid,4,'float64');
p=fread(fid,inf,'float64');
fclose(fid);
