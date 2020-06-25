#include "su.h"
#include "stdio.h"
#include "segy.h"
#define TRACE fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__)
void plotgather(float **d, int nh, int nt, float dt)
{
  segy tr;
  
  int  itr;
  FILE* fp;
  
  if ((fp=popen("suximage legend=1 title=\"plotgather\"","w"))==NULL){ 
    warn("Cannot open fp\n");
    return;
  }
  fprintf(stderr,"nh=%d,nt=%d,dt=%f\n",nh,nt,dt);
  for (itr=0;itr<nh;itr++){
      memcpy((void *) tr.data,
	     (const void *) d[itr],nt*sizeof(float));
      tr.offset=(int) itr;//h[itr];
      tr.tracl=itr;
      tr.dt=(int) (dt*1e6);
      tr.ns=(int) nt;
      tr.ntr=(int) nh;
      //fprintf(stderr,"==>itr=%d\n",itr);
      fflush(fp);
      fputtr(fp,&tr);

      //for (int it=0;it<nt;it++) fprintf(stderr,"tr.data[%d]=%f\n",it,tr.data[it]); 
  }    


  pclose(fp);

  return;
  
}

void plotgather(float **d, int nh, int nt, float dt, const char *s)
{
  segy tr;
  
  int  itr;
  FILE* fp;
  
  if ((fp=popen(s,"w"))==NULL){ 
    warn("Cannot open fp\n");
    return;
  }
  fprintf(stderr,"nh=%d,nt=%d,dt=%f\n",nh,nt,dt);
  for (itr=0;itr<nh;itr++){
      memcpy((void *) tr.data,
	     (const void *) d[itr],nt*sizeof(float));
      tr.offset=(int) itr;//h[itr];
      tr.tracl=itr;
      tr.dt=(int) (dt*1e6);
      tr.ns=(int) nt;
      tr.ntr=(int) nh;
      //fprintf(stderr,"==>itr=%d\n",itr);
      fflush(fp);
      fputtr(fp,&tr);

      //for (int it=0;it<nt;it++) fprintf(stderr,"tr.data[%d]=%f\n",it,tr.data[it]); 
  }    


  pclose(fp);

  return;
  
}

void gather_pipe(float **d, int nh, int nt, const char *s)
{
  FILE* fp;
  int itr;
  if ((fp=popen(s,"w"))==NULL){ 
    warn("Cannot open fp\n");
    return;
  }

  //fprintf(stderr,"In plotgather nh=%d,nt=%d\n",nh,nt);
  for (itr=0;itr<nh;itr++)  
    efwrite(d[itr],sizeof(float),nt,fp);

  pclose(fp);

  return;
  
}

void plotvector(float *d, int nt, const char *s)
{
  FILE* fp;
  int it;
  if ((fp=fopen("wavelet","w"))==NULL){ 
    warn("Cannot open fp\n");
    return;
  }

  for (it=0;it<nt;it++) fprintf(stderr,"d(%d)=%f\n",it,d[it]);

  efwrite(d,sizeof(float),nt,fp);
  fclose(fp);
  system(s);
  


  return;
  
}

void plotcurves(float *d, int n1, int n2, const char *s)
{
  FILE* fp;

  char buf[80];
  fp=efopen(s,"w");
  efwrite(d,sizeof(float),n1*n2,fp);
  fclose(fp);

  sprintf(buf,"xgraph < %s n=%d nplot=%d pairs=2 d1=1 style=normal title=%s & \n",s,n1,n2,s);
  system(buf);
  


  return;
  
}

void plotcurves(float **d, int n1, int n2, int f2, int l2, const char *s)
{
  FILE* fp;
  int i;
  
  char buf[80];
  fp=efopen(s,"w");
  for (i=f2;i<l2;i++) efwrite(d[i],sizeof(float),n1,fp);
  fclose(fp);

  sprintf(buf,"xgraph < %s n=%d nplot=%d pairs=2 d1=1 style=normal title=%s & \n",
	  s,n1,l2-f2,s);

  system(buf);

  return;
  
}

void plotgather(float *d, int n1, int n2, const char *s)
{
  FILE* fp;

  char buf[80];
  fp=efopen(s,"w");
  efwrite(d,sizeof(float),n1*n2,fp);
  fclose(fp);

  sprintf(buf,"xwigb < %s n1=%d n2=%d  perc=90 title=%s & \n",s,n1,n2,s);
  system(buf);
  
  fclose(fp);

  return;
  
}

void plotgather(float **d, int n1, int n2, const char *s, float perc)
{
  FILE* fp;

  char buf[80];
  int i2;
  fp=efopen(s,"w");

  for (i2=0;i2<n2;i2++)  
    efwrite(d[i2],sizeof(float),n1,fp);
  
  fclose(fp);
  sprintf(buf,"xwigb < %s n1=%d n2=%d  perc=%f title=%s & \n",s,n1,n2,perc,s);
  system(buf);

  return;
  
}



void plotgather_pipe(float *d, int n2, int n1, const char *s)
{
  FILE* fp;

  char buf[80];

  sprintf(buf,"xwigb n1=%d n2=%d perc=98  title=%s  \n",n1,n2,s);

  if ((fp=popen(buf,"w"))==NULL){ 
    warn("Cannot open fp\n");
    return;
  }

  efwrite(d,sizeof(float),n1*n2,fp);

  pclose(fp);

  return;
  
}

void plotgather_pipe(float **d, int n2, int n1, const char *s)
{
  FILE* fp;
  int i2;
  char buf[80];

  sprintf(buf,"xwigb n1=%d n2=%d perc=98  title=%s  \n",n1,n2,s);

  if ((fp=popen(buf,"w"))==NULL){ 
    warn("Cannot open fp\n");
    return;
  }

  //fprintf(stderr,"In plotgather nh=%d,nt=%d\n",nh,nt);
  for (i2=0;i2<n2;i2++)  
    efwrite(d[i2],sizeof(float),n1,fp);

  pclose(fp);

  return;
  
}





