function [p]=read_bin(name)
fid=fopen(name,'r');
p=fread(fid,inf,'float32');
fclose(fid);
