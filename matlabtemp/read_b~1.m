function [p]=read_bin_double(name)
fid=fopen(name,'r');
p=fread(fid,inf,'float64');
fclose(fid);
