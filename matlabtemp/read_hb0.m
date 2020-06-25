function [p]=read_hb0(name)
fid=fopen(name,'r');
p=fread(fid,inf,'float32');
fclose(fid);
