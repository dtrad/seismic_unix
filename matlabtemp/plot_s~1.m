function plot_ss_td(file)
p=read_bin(file);
nx=max(size(p))/2048;
 pp=reshape(p,2048,nx);
 figure,wigb(pp)