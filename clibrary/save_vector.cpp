#include "su.h"
//#include "stdio.h"
#include "segy.h"
#define TRACE fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__)
segy cleansegy(segy tr){


  fprintf(stderr,"cleaning segy..\n");
  tr.tracl=0;
  tr.tracr=0;	/* trace sequence number within reel */
  tr.fldr=0;	/* field record number */
  tr.tracf=0;	/* trace number within field record */
  tr.ep=0;	/* energy source potr.number */
  tr.cdp=0;	/* CDP ensemble number */
  tr.cdpt=0;	/* trace number within CDP ensemble */
  tr.trid=0;	/* trace identification code:*/
  tr.nvs=0;	/* number of vertically summed traces (see vscode
			   in bhed structure) */
  tr.nhs=0;	/* number of horizontally summed traces (see vscode
			   in bhed structure) */
  tr.duse=0;	/* data use:
		   1 = production
				2 = test */
  tr.offset=0;	/* distance from source point to receiver
		   group (negative if opposite to direction
		   in which the line was shot) */
  tr.gelev=0;	/* receiver group elevation from sea level
			   (above sea level is positive) */
  tr.selev=0;	/* source elevation from sea level
		   (above sea level is positive) */
  tr.sdepth=0;	/* source depth (positive) */
  tr.gdel=0;	/* datum elevation at receiver group */
  tr.sdel=0;	/* datum elevation at source */
  tr.swdep=0;	/* water depth at source */
  tr.gwdep=0;	/* water depth at receiver group */
  tr.scalel=0;	/* scale factor for previous 7 entries
		   with value plus or minus 10 to the
		   power 0, 1, 2, 3, or 4 (if positive,
			   multiply, if negative divide) */
  tr.scalco=0;	/* scale factor for next 4 entries
		   with value plus or minus 10 to the
		   power 0, 1, 2, 3, or 4 (if positive,
		   multiply, if negative divide) */
  tr. sx=0;	/* X source coordinate */
  tr. sy=0;	/* Y source coordinate */
  tr. gx=0;	/* X group coordinate */
  tr. gy=0;	/* Y group coordinate */
  tr.counit=0;	/* coordinate units code:
		   for previous four entries
		   1 = length (meters or feet)
		   2 = seconds of arc (in this case, the
		   X values are longitude and the Y values
		   are latitude, a positive value designates
		   the number of seconds east of Greenwich
		   or north of the equator */
  tr.wevel=0;	/* weathering velocity */
  tr.swevel=0;	/* subweathering velocity */
  tr.sut=0;	/* uphole time at source */
  tr.gut=0;	/* uphole time at receiver group */
  tr.sstat=0;	/* source static correction */
  tr.gstat=0;	/* group static correction */
  tr.tstat=0;	/* total static applied */
  tr.laga=0;	/* lag time A, time in ms between end of 240-
		   byte trace identification header and time
		   break, positive if time break occurs after
		   end of header, time break is defined as
		   the initiation pulse which maybe recorded
		   on an auxiliary trace or as otherwise
		   specified by the recording system */
  
  tr.lagb=0;	/* lag time B, time in ms between the time break
		   and the initiation time of the energy source,
		   may be positive or negative */
  
  tr.delrt=0;	/* delay recording time, time in ms between
		   initiation time of energy source and time
		   when recording of data samples begins
		   (for deep water work if recording does not
		   start at zero time) */
  
  tr.muts=0;	/* mute time--start */
  tr.mute=0;	/* mute time--end */
  tr.ns=0;	/* number of samples in this trace */
  tr.dt=0;	/* sample interval=0; in micro-seconds */
  tr.gain=0;	/* gain type of field instruments code:
		   1 = fixed
		   2 = binary
		   3 = floating point
				4 ---- N = optional use */
  tr.igc=0;	/* instrument gain constant */
  tr.igi=0;	/* instrument early or initial gain */
  tr.corr=0;	/* correlated:
		   1 = no
		   2 = yes */
  tr.sfs=0;	/* sweep frequency at start */
  tr.sfe=0;	/* sweep frequency at end */
  tr.slen=0;	/* sweep length in ms */
  tr.styp=0;	/* sweep type code:
		   1 = linear
		   2 = cos-squared
		   3 = other */
  tr.stas=0;	/* sweep trace length at start in ms */
  tr.stae=0;	/* sweep trace length at end in ms */
  tr.tatyp=0;	/* taper type: 1=linear, 2=cos^2, 3=other */
  tr.afilf=0;	/* alias filter frequency if used */
  tr.afils=0;	/* alias filter slope */
  tr.nofilf=0;	/* notch filter frequency if used */
  tr.nofils=0;	/* notch filter slope */
  tr.lcf=0;	/* low cut frequency if used */
  tr.hcf=0;	/* high cut frequncy if used */
  tr.lcs=0;	/* low cut slope */
  tr.hcs=0;	/* high cut slope */
  tr.year=0;	/* year data recorded */
  tr.day=0;	/* day of year */
  tr.hour=0;	/* hour of day (24 hour clock) */
  tr.minute=0;	/* minute of hour */
  tr.sec=0;	/* second of minute */
  tr.timbas=0;	/* time basis code:
		   1 = local
		   2 = GMT
		   3 = other */
  tr.trwf=0;	/* trace weighting factor, defined as 1/2^N
		   volts for the least sigificant bit */
  tr.grnors=0;	/* geophone group number of roll switch
		   position one */
  tr.grnofr=0;	/* geophone group number of trace one within
		   original field record */
  tr.grnlof=0;	/* geophone group number of last trace within
		   original field record */
  tr.gaps=0;	/* gap size (total number of groups dropped) */
  tr.otrav=0;	/* overtravel taper code:
		   1 = down (or behind)
		   2 = up (or ahead) */
  /* local assignments */
  tr.d1=0;	/* sample spacing for non-seismic data */
  tr.f1=0;	/* first sample location for non-seismic data */
  tr.d2=0;	/* sample spacing between traces */
  tr.f2=0;	/* first trace location */
  tr.ungpow=0;	/* negative of power used for dynamic
		   range compression */
  tr.unscale=0;	/* reciprocal of scaling factor to normalize
		   range */
  tr.ntr=0; 	/* number of traces */
  tr.mark=0;	/* mark selected traces */
  
  tr.shortpad=0; /* alignment padding */
  return(tr);
}






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








