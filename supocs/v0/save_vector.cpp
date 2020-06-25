#include "su.h"
//#include "stdio.h"
#include "segy.h"
#define TRACE fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__)
segy cleansegy(segy tr);
void save_vector(float *d, int n, char* name)
{
  int nw;
  FILE* fp;

  if ((fp=fopen(name,"w"))==NULL){ 
    warn("Cannot open fp\n");
    return;
  }
  nw=efwrite(d,sizeof(float),n,fp);
  fprintf(stderr,"nw=%d,n=%d\n",nw,n);
  if (nw!=n) warn("writing only %d of %d\n",nw,n);
  fclose(fp);
  return;
  
}

void save_gather(float **d, int nh, float *h, int nt, float dt, char* name)
{
  
  segy tr;
  tr=cleansegy(tr);
  int type=1; // Default = t - offset gather
              // option  type=2 : t - q Radon gather
  int itr;
  FILE* fp;
  if (fabs(h[1]-h[0]) < 0.1 ) type=2;

  if ((fp=fopen(name,"w"))==NULL){ 
    warn("Cannot open fp\n");
    return;
  }
  fprintf(stderr,"nh=%d,nt=%d,dt=%f\n",nh,nt,dt);
  for (itr=0;itr<nh;itr++){
      memcpy((void *) tr.data,
	     (const void *) d[itr],nt*sizeof(float));

      if (type==1) tr.offset=(int) h[itr];
      else if (type==2) tr.f2=h[itr];

      tr.tracl=itr;
      tr.dt=(int) (dt*1e6);
      tr.ns=(int) nt;
      tr.ntr=(int) nh;
      tr.delrt=0;
      tr.trid=1;
      //fprintf(stderr,"==>itr=%d\n",itr);
      fvputtr(fp,&tr);

      //for (int it=0;it<nt;it++) fprintf(stderr,"tr.data[%d]=%f\n",it,tr.data[it]); 
  }    

  fflush(fp);
  fclose(fp);

  return;
  
}

void save_gather(float ***d, int nx, float *vel, int nh, float *h, int nt, float dt, char* name)
{
  /* Save a 3 dimensional gather including the header words for all axes */

  segy tr;
  tr=cleansegy(tr);

  int ix, ih;
  FILE* fp;

  if ((fp=fopen(name,"w"))==NULL){ 
    warn("Cannot open fp\n");
    return;
  }

  for (ix=0;ix<nx;ix++){
    for (ih=0;ih<nh;ih++){
      memcpy((void *) tr.data,(const void *) d[ix][ih],nt*sizeof(float));
      tr.offset=(int) h[ih];
      tr.f2=vel[ix];
      tr.tracl=ix*nh+ih;
      tr.dt=(int) (dt*1e6);
      tr.ns=(int) nt;
      tr.ntr=(int) nh*nx;
      //tr.delrt=(int) (t[0]*1e6);
      tr.trid=1;
      fvputtr(fp,&tr);
    }
  }    

  fflush(fp);
  fclose(fp);

  return;
  
}


void save_gather(float ***d, int nx, float *vel, int nh, float *h, int nt, float *t, char* name)
{
  /* Save a 3 dimensional gather including the header words for all axes */

  segy tr;
  tr=cleansegy(tr);

  int ix, ih;
  float dt=t[1]-t[0];
  FILE* fp;

  if ((fp=fopen(name,"w"))==NULL){ 
    warn("Cannot open fp\n");
    return;
  }

  for (ix=0;ix<nx;ix++){
    for (ih=0;ih<nh;ih++){
      memcpy((void *) tr.data,(const void *) d[ix][ih],nt*sizeof(float));
      tr.offset=(int) h[ih];
      tr.f2=vel[ix];
      tr.tracl=ix*nh+ih;
      tr.dt=(int) (dt*1e6);
      tr.ns=(int) nt;
      tr.ntr=(int) nh*nx;
      tr.delrt=(int) (t[0]*1e6);
      tr.trid=1;
      fvputtr(fp,&tr);
    }
  }    

  fflush(fp);
  fclose(fp);

  return;
  
}



void save_gather(float *d, int nh, float *h, int nt, float dt, char* name)
{
  
  segy trr;
  trr=cleansegy(trr);
  int type=1; // Default = t - offset gather
              // option  type=2 : t - q Radon gather
  int itr;
  FILE* fp;
  if (fabs(h[1]-h[0]) < 0.1 ) type=2;

  if ((fp=fopen(name,"w"))==NULL){ 
    warn("Cannot open fp\n");
    return;
  }
  fprintf(stderr,"nh=%d,nt=%d,dt=%f\n",nh,nt,dt);
  for (itr=0;itr<nh;itr++){
      memcpy((void *) trr.data,
	     (const void *) &d[itr*nt],nt*sizeof(float));

      if (type==1) trr.offset=(int) h[itr];
      else if (type==2) trr.f2=h[itr];

      trr.tracl=itr;
      trr.dt=(int) (dt*1e6);
      trr.ns=(int) nt;
      trr.ntr=(int) nh;
      trr.delrt=0;
      trr.trid=1;
      fvputtr(fp,&trr);

  }    

  fflush(fp);
  fclose(fp);

  return;
  
}

void save_gather_tx(float **d, int nh, float *h, int nt, float dt, char* name)
{
  /* The same as before but the gather is d[t][h] instead of d[h][t] */
  segy tr;
  tr=cleansegy(tr);
  int type=1; // Default = t - offset gather
  int it;            // option  type=2 : t - q Radon gather
  int itr;
  FILE* fp;
  if (fabs(h[1]-h[0]) < 0.1 ) type=2;

  if ((fp=fopen(name,"w"))==NULL){ 
    warn("Cannot open fp\n");
    return;
  }
  fprintf(stderr,"nh=%d,nt=%d,dt=%f\n",nh,nt,dt);
  for (itr=0;itr<nh;itr++){
      for (it=0;it<nt;it++) tr.data[it]=d[it][itr];
      if (type==1) tr.offset=(int) h[itr];
      else if (type==2) tr.f2=h[itr];

      tr.tracl=itr;
      tr.dt=(int) (dt*1e6);
      tr.ns=(int) nt;
      tr.ntr=(int) nh;
      //fprintf(stderr,"==>itr=%d\n",itr);
      fvputtr(fp,&tr);

      //for (int it=0;it<nt;it++) fprintf(stderr,"tr.data[%d]=%f\n",it,tr.data[it]); 
  }    

  fflush(fp);
  fclose(fp);

  return;
  
}



void save_gather(float **d, int nh, int nt, float dt, char* name)
{

  int  itr;
  FILE* fp;
  segy trr;
  trr=cleansegy(trr);

 
  if ((fp=fopen(name,"w"))==NULL){ 
    warn("Cannot open fp\n");
    return;
  }
  //fprintf(stderr,"nh=%d,nt=%d,dt=%f\n",nh,nt,dt);
  for (itr=0;itr<nh;itr++){
      memcpy((void *) trr.data,
	     (const void *) d[itr],nt*sizeof(float));
      trr.tracl=itr+1;
      //fprintf(stderr,"trr.dt=%d, dt=%f\n",trr.dt,dt);
      trr.dt=(int) (dt*1e6);
      trr.ns=nt;
      trr.ntr=nh;
      //fprintf(stderr,"==>itr=%d\n",itr);
      fvputtr(fp,&trr);

      //for (int it=0;it<nt;it++) fprintf(stderr,"tr.data[%d]=%f\n",it,tr.data[it]); 
  }    

  //fflush(fp);
  fclose(fp);
  
  return;
  
}


void save_gather(float ***d, int nx, int nh, int nt, float dt, char* name)
{

  int  ih, ix;
  FILE* fp;
  segy trr;
  trr=cleansegy(trr);

  TRACE;
  if ((fp=fopen(name,"w"))==NULL){ 
    warn("Cannot open fp\n");
    return;
  }

  for (ix=0;ix<nx;ix++){
    for (ih=0;ih<nh;ih++){
      memcpy((void *) trr.data,(const void *) d[ix][ih],nt*sizeof(float));
      trr.tracl=ih+1;
      trr.dt=(int) (dt*1e6);
      trr.ns=nt;
      trr.ntr=nh*nx;
      fvputtr(fp,&trr);
    }
  }    
  TRACE;

  fclose(fp);
  
  return;
  
}


void save_gather(float *d, int nh, int nt, float dt, char* name)
{
  segy tr;
  tr=cleansegy(tr);
  int  itr, it;
  FILE* fp;
  
  if ((fp=fopen(name,"w"))==NULL){ 
    warn("Cannot open fp\n");
    return;
  }
  fprintf(stderr,"nh=%d,nt=%d,dt=%f\n",nh,nt,dt);
  for (itr=0;itr<nh;itr++){
    for (it=0;it<nt;it++) tr.data[it]=d[itr*nt+it];
    tr.offset=(int) itr;//h[itr];
    tr.tracl=itr;
    tr.dt=(int) (dt*1e6);
    tr.ns=(int) nt;
    tr.ntr=(int) nh;
    //fprintf(stderr,"==>itr=%d\n",itr);
    fvputtr(fp,&tr);
      //for (int it=0;it<nt;it++) fprintf(stderr,"tr.data[%d]=%f\n",it,tr.data[it]); 
  }    

  fflush(fp);
  fclose(fp);

  return;
  
}


void save_gather(complex **d, int nh, float *h, int nt, float dt, int nfft, char* name)
{
  segy trf;
  
  int  itr, it;
  FILE* fpf;
  
  if ((fpf=fopen(name,"w+"))==NULL){ 
    err("Cannot open fp\n");
    return;
  }
  fprintf(stderr,"nh=%d,nt=%d,dt=%f,nfft=%d\n",nh,nt,dt,nfft);
  for (itr=0;itr<nh;itr++){
      for (it=0;it<nfft/2;it++) trf.data[it]=abs(d[it][itr]);
      for (it=nfft/2+1;it<nt;it++) trf.data[it]=0;
      trf.offset=(int) h[itr];
      trf.d1=1./(nfft*dt);
      //tr.dt=(int) 1000;
      trf.f1=0;
      trf.ns=(int) nfft/2+1;
      //trf.ns=(int) nt;
      trf.ntr=(int) nh;
      trf.tracl=itr;
      
      fvputtr(fpf,&trf);

      fprintf(stderr,"itr=%d\n",itr);
      //for (int it=0;it<nt;it++) fprintf(stderr,"tr.data[%d]=%f\n",it,tr.data[it]); 
  }    
      fprintf(stderr,"itr=%d\n",itr);
      //   fvputtr(fpf,&trf);
  fclose(fpf);

  return;
  
}

void save_gather(complex **d, int nh, int nt, float dt, float dx, char* name)
{
  segy trf;
  
  int  itr, it;
  FILE* fpf;
  
  if ((fpf=fopen(name,"w+"))==NULL){ 
    err("Cannot open fp\n");
    return;
  }
  fprintf(stderr,"nh=%d,nt=%d,dt=%f\n",nh,nt,dt);
  for (itr=0;itr<nh;itr++){
      for (it=0;it<nt;it++) trf.data[it]=abs(d[itr][it]);
      trf.d1=dt;
      trf.d2=dx;
      trf.f1=0;
      trf.f2=-dx*nh/2.;
      trf.ns=(int) nt;
      trf.ntr=(int) nh;
      trf.tracl=itr;
      fvputtr(fpf,&trf);
  }    
  fprintf(stderr,"itr=%d\n",itr);

  fclose(fpf);

  return;
  
}


void save1dfile(float *d, int n1, const char *s)
{
  FILE* fp;

  if ((fp=fopen(s,"w"))==NULL){ 
    warn("Cannot open fp\n");
    return;
  }

  efwrite(d,sizeof(float),n1,fp);

  fclose(fp);

  return;
  
}

void save2dfile(float **d, int nh, int nt, float dt, const char *s)
{
  FILE* fp;
  int itr;
  if ((fp=fopen(s,"w"))==NULL){ 
    warn("Cannot open fp\n");
    return;
  }

  fprintf(stderr,"In plotgather nh=%d,nt=%d,dt=%f\n",nh,nt,dt);
  
  for (itr=0;itr<nh;itr++)  
    efwrite(d[itr],sizeof(float),nt,fp);

  fclose(fp);

  return;
  
}

void save2dfile(float **d, int n1, int n2, const char *s, int transpose)
{
  FILE* fp;
  int i1,i2;
  if ((fp=fopen(s,"w"))==NULL){ 
    warn("Cannot open fp\n");
    return;
  }

  fprintf(stderr,"In plotgather n1=%d,n2=%d\n",n1,n2);
  // Normal saving 
  if (!transpose) 
    for (i2=0;i2<n2;i2++)  
      efwrite(d[i2],sizeof(float),n1,fp);
  else // transposes the matrix 
    for (i1=0;i1<n1;i1++)
     for (i2=0;i2<n2;i2++)    
       efwrite(&d[i2][i1],sizeof(float),1,fp);

  fclose(fp);

  return;
  
}


void save3dfile(complex ***d, int n1, int n2, int n3, const char *s)
{
  FILE* fp;
  int i2, i3;
  if ((fp=fopen(s,"w"))==NULL){ 
    warn("Cannot open fp\n");
    return;
  }

  fprintf(stderr,"In plotgather n1=%d,n2=%d\n",n1,n2);
  // Normal saving 
   
  for (i3=0;i3<n3;i3++)
    for (i2=0;i2<n2;i2++)  
      efwrite(d[i3][i2],sizeof(complex),n1,fp);


  fclose(fp);

  return;
  
}


void save2dfile(complex **d, int n1, int n2, const char *s)
{
  FILE* fp;
  int i2;
  if ((fp=fopen(s,"w"))==NULL){ 
    warn("Cannot open fp\n");
    return;
  }

  fprintf(stderr,"In plotgather n1=%d,n2=%d\n",n1,n2);
  // Normal saving 
   
  
  for (i2=0;i2<n2;i2++)  
    efwrite(d[i2],sizeof(complex),n1,fp);


  fclose(fp);

  return;
  
}



void read1dfile(float *d, int n1, const char *s)
{
  FILE* fp;
  if ((fp=fopen(s,"r"))==NULL){ 
    warn("Cannot open fp\n");
    return;
  }

  efread(d,sizeof(float),n1,fp);

  fclose(fp);

  return;
  
}



void read2dfile(complex **d, int n1, int n2, const char *s)
{
  FILE* fp;
  int i2;
  if ((fp=fopen(s,"r"))==NULL){ 
    warn("Cannot open fp\n");
    return;
  }

  fprintf(stderr,"In plotgather n1=%d,n2=%d\n",n1,n2);

  for (i2=0;i2<n2;i2++)  
    efread(d[i2],sizeof(complex),n1,fp);

  fclose(fp);

  return;
  
}

void read3dfile(complex ***d, int n1, int n2, int n3, const char *s)
{
  FILE* fp;
  int i2, i3;
  if ((fp=fopen(s,"r"))==NULL){ 
    warn("Cannot open fp\n");
    return;
  }

  fprintf(stderr,"In plotgather n1=%d,n2=%d\n",n1,n2);
  // Normal saving 
   
  for (i3=0;i3<n3;i3++)
    for (i2=0;i2<n2;i2++)  
      efread(d[i3][i2],sizeof(complex),n1,fp);


  fclose(fp);

  return;
  
}








