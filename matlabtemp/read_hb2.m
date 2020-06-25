function [p,a,b]=read_hb(name)
fid=fopen(name,'r');
a=fread(fid,2,'int32');
b=fread(fid,4,'float32');
p=fread(fid,inf,'float32');
fclose(fid);
