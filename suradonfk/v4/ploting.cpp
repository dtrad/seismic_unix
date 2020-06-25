#include "radonfk.h"

void xplotgather(float **d, int nh, int nt, float dt, char *s, int num, char *s2)
{
  char buf[180];
  save_gather(d,nh,nt,dt,s);
  //sprintf(buf,"suximage < %s title=%s%d curve=curve1 npair=5 %s\n",s,s,num,s2);
  
  sprintf(buf,"suximage < %s title=%s%d legend=1 wbox=500 %s \n",s,s,num,s2);
  fprintf(stderr,buf);  
  system(buf);
  return;
}

void xplotdiff(float **d1, float **d2, int nh, int nt, float dt, char *s,  char *s2)
{
  char buf[180];
  int ih,it;
  float **d3=0;
  d3= ealloc2float(nt,nh);
  for (ih=0;ih<nh;ih++) for(it=0;it<nt;it++) d3[ih][it]=d1[ih][it]-d2[ih][it];

  save_gather(d3,nh,nt,dt,s);
  sprintf(buf,"suximage < %s title=%s legend=1 wbox=500 %s \n",s,s,s2);
  system(buf);
  
  free2float(d3);

  return;
}

/* test for static variable used in fputtr */
void plot_after_stretch(float **d, int nh, int nt, int dt, char *s1, char *s2)
{
  xplotgather(d,nh,nt,dt,s1,0,s2);
  return;
}

void ximageplotgather(float **d, int nh, int nt, float dt, char *s, int num, char *s2)
{
  char buf[180];
  save_gather(d,nh,nt,dt,s);
  sprintf(buf,"suximage < %s title=%s%d curve=curve1 npair=5 hbox=900 wbox=700 %s legend=1 \n",s,s,num,s2);
  system(buf);
  return;
}

void plotimage(float **d, int n1, int n2, const char *s, float perc)
{
  FILE* fp;

  char buf[80];
  int i2;
  fp=efopen(s,"w");

  for (i2=0;i2<n2;i2++)  
    efwrite(d[i2],sizeof(float),n1,fp);
  
  fclose(fp);
  sprintf(buf,"ximage < %s n1=%d n2=%d  perc=%f title=%s legend=1 & \n",s,n1,n2,perc,s);
  system(buf);

  return;
  
}



void save2d(float **d, int nh, int nt, float dt, char* name)
{
  segy tr2;
  
  int  itr;
  FILE* fp;
  
  if ((fp=fopen(name,"w"))==NULL){ 
    warn("Cannot open fp\n");
    return;
  }
  fprintf(stderr,"nh=%d,nt=%d,dt=%f\n",nh,nt,dt);
  for (itr=0;itr<nh;itr++){
      memcpy((void *) tr2.data,
	     (const void *) d[itr],nt*sizeof(float));
      tr2.tracl=itr+1;
      tr2.dt=(int) (dt*1e6);
      tr2.ns=nt;
      tr2.ntr=nh;
      fprintf(stderr,"==>itr=%d\n",itr);
      fputtr(fp,&tr2);
      
      //for (int it=0;it<nt;it++) fprintf(stderr,"tr2.data[%d]=%f\n",it,tr2.data[it]); 
  }    

  fflush(fp);
  fclose(fp);
  
  return;
  
}



// FK spectrum using standard format (symmetric w, non symmetric k)
void plotAmpSpec(float **p, int nx, int nt, float dt, char *c)
{
  int nw, nk, ntfft, nxfft;
  int it,ix;
  complex czero;czero.r=czero.i=0;
  int sign=1;
  char buf[12];
  float scale;
  
  ntfft = npfar(nt);
  nxfft = npfa(nx);
  nw=ntfft/2+1;
  nk=nxfft;
  scale=nxfft*ntfft;
  //df=1./(ntfft*dt);
  //dk=1./(nxfft*dx);
  float**   pfft  = ealloc2float(ntfft,nxfft);
  complex** cpfft = ealloc2complex(nw,nk);

  /* copy data from input to FFT array and pad with zeros */
  for (ix=0;ix<nx;ix++){
    /* if ix odd, negate to center transform of dimension 2 */
    if (ISODD(ix)) sign=-1; else sign=1;
    for (it=0; it<nt; it++) pfft[ix][it]=sign*p[ix][it];
    for (it=nt; it< ntfft;it++) pfft[ix][it] = 0.0;
  }
  for (ix=nx;ix<nxfft;ix++) for(it=0;it<ntfft;it++) pfft[ix][it] = 0.0;

  /* Fourier transform t to w */
  pfa2rc(1,1,ntfft,nx,pfft[0],cpfft[0]);
  pfa2cc(-1,2,nw,nxfft,cpfft[0]);
  for (ix=0;ix<nk;ix++) for(it=0;it<nw;it++) pfft[ix][it] = abs(cpfft[ix][it]);

  sprintf(buf,"fkspec%s",c);
  if (1) plotimage(pfft,nw,nk,buf,99.0);
  free2complex(cpfft);
  free2float(pfft);
  
}

// FK spectrum using Stolt format (symmetric k, non symmetric w)
void plotAmpSpec(complex **p, int nk, int nt, float dt, char *c)
{
  int nw;
  float dw;
  float fw;
  complex *pp;
  int it,ik;
  complex czero;czero.r=czero.i=0;
  
  float **aspec; // Amplitude spectrum
  char buf[12];
  
  nw = (int) (nt*2/0.6);
  nw = npfao(nw,nw*2);
  dw = 2.0*PI/(nw*dt);
  fw = -PI/dt;

  fprintf(stderr,"Inside plotAmpSpec nw=%d nk=%d\n",nw,nk);

  aspec=ealloc2float(nw,nk);
  pp = ealloc1complex(nw);
  
  for (ik=0;ik<nk;ik++){
    for (it=0; it<nt; it+=2)
      pp[it] = p[ik][it];
    for (it=1; it<nt; it+=2) {
      pp[it].r = -p[ik][it].r;
      pp[it].i = -p[ik][it].i;
    }
    for (it=nt; it<nw; ++it) pp[it] = czero;
    pfacc(1,nw,pp);
    /* zero -Nyquist frequency for symmetry */
    pp[0] = czero;
    for (it=0; it<nw; it++) aspec[ik][it]=abs(pp[it]);
      
  }

  if (0){
    save2d(aspec,nk,nw,dw,"temp");
    system("suxwigb < temp perc=100 title=temp &"); 
  }
  sprintf(buf,"fkspec%s",c);
  if (1) plotimage(aspec,nw,nk,buf,99.0);
  free1complex(pp);
  free2float(aspec);
  
}




