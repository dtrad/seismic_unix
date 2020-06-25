#include "su.h"
int taper(float **data, int nt, int nh, int ntaper,int flag)
     /*  Given a data gather apply a taper to the edge ntaper traces     */
     /*  The taper is an \"ntaper\" point sine-squared taper 		 */
     /*  symmetrically applied at each end of the data set.		 */
     /* 0 both ends */
     /* 1 beginning */
     /* 2 end       */
     /* see sutaper.c for original function */
     /* Daniel Trad - UBC - June 20, 2000 */
{
  float *taper;
  int k;
  float s;
  int it, ih;
  if (ntaper<=0) return EXIT_SUCCESS;
  taper = ealloc1float(ntaper);
  for (k = 0; k < ntaper; ++k) {
    s = sin(k*PI/(2*ntaper));
    //taper[k] = s*s;
    taper[k]= pow(s,100.0);
    if (0) fprintf(stderr,"s=%f,taper[%d]=%f\n",s,k,taper[k]);
  }

  if(flag==0){
    for(ih=0;ih<ntaper;ih++)
      for (it=0;it<nt;it++) data[ih][it]*=taper[ih];
    for(ih=nh-ntaper;ih<nh;ih++)
      for (it=0;it<nt;it++) data[ih][it]*=taper[nh-ih-1]; 
  }


  if(flag==1)
    for(ih=0;ih<ntaper;ih++)
      for (it=0;it<nt;it++) data[ih][it]*=taper[ih];

  if(flag==2)
    for(ih=nh-ntaper;ih<nh;ih++)
      for (it=0;it<nt;it++) data[ih][it]*=taper[nh-ih-1]; 

  if(flag==4){
    for(ih=0;ih<ntaper;ih++)
      for (it=0;it<nt;it++) data[ih][it]*=0.;//taper[ih];
    for(ih=nh-ntaper;ih<nh;ih++)
      for (it=0;it<nt;it++) data[ih][it]*=0.;//taper[nh-ih-1]; 
  }

  if(flag==5){
    for(ih=0;ih<ntaper;ih++)
      for (it=0;it<nt;it++) data[ih][it]*=0.;//taper[ih];
    for(ih=nh-ntaper;ih<nh;ih++)
      for (it=0;it<nt;it++) data[ih][it]*=0.;//taper[nh-ih-1]; 
  }

  free1float(taper);

  if (flag==3)
    for (ih=ntaper;ih<nh;ih++){
      float scalemute=pow(1.0-((float) (ih-ntaper)/(nh-1-ntaper)),2.0);
      fprintf(stderr,"scalemute=%f\n",scalemute);
      for (it=0;it<nt;it++){
	data[ih][it]*=scalemute;
      }
    } 
  
  
  return EXIT_SUCCESS;
}

int taper(float *data, int nt, int nh, int ntaper,int flag)
     /*  Given a data gather apply a taper to the edge ntaper traces     */
     /*  The taper is an \"ntaper\" point sine-squared taper 		 */
     /*  symmetrically applied at each end of the data set.		 */
     /* 0 both ends */
     /* 1 beginning */
     /* 2 end       */
     /* see sutaper.c for original function */
     /* Daniel Trad - UBC - June 20, 2000 */
{
  float *taper;
  int k;
  float s;
  int it, ih;

  taper = ealloc1float(ntaper);
  for (k = 0; k < ntaper; ++k) {
    s = sin(k*PI/(2*ntaper));
    taper[k] = pow(s,100.0);
    if (0) fprintf(stderr,"s=%f,taper[%d]=%f\n",s,k,taper[k]);
  }

  if(flag==0){
    for(ih=0;ih<ntaper;ih++)
      for (it=0;it<nt;it++) data[ih*nt+it]*=taper[ih];
    for(ih=nh-ntaper;ih<nh;ih++)
      for (it=0;it<nt;it++) data[ih*nt+it]*=taper[nh-ih-1]; 
  }


  if(flag==1)
    for(ih=0;ih<ntaper;ih++)
      for (it=0;it<nt;it++) data[ih*nt+it]*=taper[ih];

  if(flag==2)
    for(ih=nh-ntaper;ih<nh;ih++)
      for (it=0;it<nt;it++) data[ih*nt+it]*=taper[nh-ih-1]; 

  free1float(taper);

  if (flag==3)
    for (ih=ntaper;ih<nh;ih++){
      float scalemute=pow(1.0-((float) (ih-ntaper)/(nh-1-ntaper)),2.0);
      fprintf(stderr,"scalemute=%f\n",scalemute);
      for (it=0;it<nt;it++){
	data[ih*nt+it]*=scalemute;
      }
    } 

  if(flag==5)
    for(ih=nh-ntaper;ih<nh;ih++)
      for (it=0;it<nt;it++) data[ih*nt+it]*=0;   
  
  return EXIT_SUCCESS;
}


int taper(float *data, int nh, int ntaper)
{
  float *taper;
  int k, ih;
  float s;
  
  taper=ealloc1float(ntaper);

  for (k=0;k<ntaper;++k){ 
    s=sin(k*PI/(2*ntaper));
    taper[k]=s;
  }  
  taper[0]=1e-5; 
  /* Taper at the left end of the data set */
  for (ih = 0; ih < ntaper; ++ih) data[ih] /= taper[ih];

  /* Taper at the right end of the data set */
  for (ih=nh - ntaper ; ih < nh; ++ih)  data[ih] /= taper[nh - ih - 1];

  return EXIT_SUCCESS;

}










