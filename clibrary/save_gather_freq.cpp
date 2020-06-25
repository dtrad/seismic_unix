#include "su.h"
#include "segy.h"
#include "clibrary.h"

void save_gather_freq(complex **d, int nh, float *h, int nt, float dt, int nfft, char* name)
{
  segy trf;
  
  int itr, it;
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
      
      fputtr(fpf,&trf);

      fprintf(stderr,"itr=%d\n",itr);
      //for (int it=0;it<nt;it++) fprintf(stderr,"tr.data[%d]=%f\n",it,tr.data[it]); 
  }    
      fprintf(stderr,"itr=%d\n",itr);
      //   fputtr(fpf,&trf);
  fclose(fpf);

  return;
  
}
