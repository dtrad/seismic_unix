void elliptical_taper(float **Wd, int np, int nt, int ntaper, int ntapertime, int taperopt, 
		      float *p, float *t, float dt)

{
  int centre;
  float min=fabs(p[0]);
  
  if (taper==1||taper==3){ 
    for (ip=0;ip<np;ip++){
      if(fabs(p[ip])<min){
	min=fabs(p[ip]);
	centre=ip;
      }
    }
    tapercentre(Wd,np,nt,ntaper,ntapertime,centre);
    if (taper==1){
      save_gather(Wd,np,nt,dt,"Wd");
      system("suxwigb < Wd clip=1.5  title=\"plotgather\" &");
    }
  }
