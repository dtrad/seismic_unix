#include "dan.h"
#include "su.h"
#define SLOPE 2

int taper(float **data, int nt, int nh, int ntaper, int flag);

/* Function to filter non hyperbolic events in CSP gathers
   Because the forward operator can take different shapes I need to write different 
   routines, but the only difference is the forward operator used. 
   The name rest_multiples does not make sense but I keep it for compatibility
   with the previous program suradonlinetd that I use as the main body of this
   denoising algorithm. 

   Some experimental variables.
   plot: =1 plot gathers before and after a
   itm:  First multiple is below this time

   Daniel Trad UBC - June 2001

*/

void rest_multiples(void  (*oper) (float *, float *, unsigned int **, int, int, 
				   int, int, int, float *, int), 
		    float **model, float **data, float **Wd, unsigned int **index, 
		    int adj, int nt, int nh, int nq, int nsparse, float *wavelet, 
		    int nw, float parmute, float *q, float *h, int itm, int plot)
{
  int side=1;  // side 2 --> right mute ; side 1 left mute
  float **datarec=0; // multiples
  float dpd;
  float dpdp;
  float scale; 
  float **modeltemp=0;
  int iq,ih,it;
  int nmute=0;
  int slope=SLOPE; // Defines a slope in the muting

  modeltemp=ealloc2float(nt,nq);
  datarec=ealloc2float(nt,nh);
  

  side=1;
  iq=0; while(q[iq]<parmute) iq++; 
  if (side==2) nmute=nq-iq;
  else if (side==1) nmute=iq;
  fprintf(stderr,"MUTING 2 at nmute=%d************************\n",nmute);

  for (iq=0;iq<nq;iq++)
    memcpy(modeltemp[iq],model[iq],nt*FSIZE);

  if ((plot==2)||(plot==3)){
    //system("rm multiples.su");
    save_gather(data,nh,h,nt,0.004,"csp.su");
    system("suxwigb < csp.su clip=1 key=offset  title=csp xbox=0 &");
  }
  if ((plot==2)||(plot==3)){
    //system("rm mutedmodel.su");
    save_gather(modeltemp,nq,q,nt,0.004,"nomutedmodel.su");
    system("suxwigb < nomutedmodel.su perc=99 key=f2 title=nomutedmodel &");
  }
  
  if (0) taper(modeltemp,nt,nq,nmute,side);
  else{
    for (iq=0;iq<nq;iq++) 
      for (it=0;it<(itm-(iq-nmute)*slope);it++)  // 0 should be it0 
	modeltemp[iq][it]=0;
  }

  // Temporal change to leave primaries untouched at the shallow times
  for (iq=0;iq<nq;iq++) 
    for (it=0;it<(itm-(iq-nmute)*slope);it++) 
      modeltemp[iq][it]=0;   

  if ((plot==2)||(plot==3)){
    //system("rm mutedmodel.su");
    save_gather(modeltemp,nq,q,nt,0.004,"mutedmodel.su");
    system("suxwigb < mutedmodel.su perc=99 key=f2 title=mutedmodel xbox=600 ");
  }

  oper(modeltemp[0],datarec[0],index,0,nt,nh,nq,nsparse,wavelet,nw);
  /* save the predicted data before applying back the data weighting */
  if ((plot==2)||(plot==3)){
    //system("rm multiples.su");
    save_gather(datarec,nh,h,nt,0.004,"predicted.su");
    system("suxwigb < predicted.su clip=1 key=offset  title=predicted xbox=600 &");
  }
  /* Save the predicted data afted putting back the amplitudes through Wd */
  for (ih=0;ih<nh;ih++) for (it=0;it<nt;it++) datarec[ih][it]*=Wd[ih][it];
  if ((plot==2)||(plot==3)){
    //system("rm multiples.su");
    save_gather(datarec,nh,h,nt,0.004,"multiples.su");
    system("suxwigb < multiples.su clip=1 key=offset  title=multiples xbox=600 &");
  }

  fprintf(stderr,"substracting multiples \n");
  dpd=dot(nh*nt,data[0],datarec[0]);
  dpdp=dot(nh*nt,datarec[0],datarec[0]);
  scale=dpd/dpdp;
  fprintf(stderr,"scale 2===>%f\n",scale);

  for (ih=0;ih<nh;ih++) 
    for (it=0;it<nt;it++)
      data[ih][it]=data[ih][it]-scale*datarec[ih][it];

  if ((plot==2)||(plot==3)){
    //system("rm multiples.su");
    save_gather(data,nh,h,nt,0.004,"primaries.su");
    system("suxwigb < primaries.su clip=1 key=offset  title=primaries ");
  }  
    
  
  free2float(modeltemp);
  free2float(datarec);
  return;
}




void rest_multiples(void  (*oper) (float *, float *, unsigned int **, int, int, 
				   int, int, int), 
		    float **model, float **data, float **Wd, unsigned int **index, 
		    int adj, int nt, int nh, int nq, int nsparse, float parmute, 
		    float *q, float *h, int itm, int plot)
{
  int side=1;  // side 2 --> right mute ; side 1 left mute
  float **datarec=0; // multiples
  float dpd;
  float dpdp;
  float scale;
  float **modeltemp=0;
  int iq,ih,it;
  int nmute=0;
  int slope=SLOPE; // Defines a slope in the muting
  
  modeltemp=ealloc2float(nt,nq);
  datarec=ealloc2float(nt,nh);
  

  side=1;
  iq=0; while(q[iq]<parmute) iq++; 
  if (side==2) nmute=nq-iq;
  else if (side==1) nmute=iq;
  fprintf(stderr,"MUTING 2 at nmute=%d************************\n",nmute);

  for (iq=0;iq<nq;iq++)
    memcpy(modeltemp[iq],model[iq],nt*FSIZE);

  if ((plot==2)||(plot==3)){
    //system("rm multiples.su");
    save_gather(data,nh,h,nt,0.004,"csp.su");
    system("suxwigb < csp.su clip=1 key=offset  title=csp xbox=0 &");
  }
  if ((plot==2)||(plot==3)){
    //system("rm mutedmodel.su");
    save_gather(modeltemp,nq,q,nt,0.004,"nomutedmodel.su");
    system("suxwigb < nomutedmodel.su perc=99 key=f2 title=nomutedmodel &");
  }

  
  if (0) taper(modeltemp,nt,nq,nmute,side);
  else{
    for (iq=0;iq<nmute;iq++) 
      for (it=0;it<nt;it++)  // 0 should be it0 
      modeltemp[iq][it]=0; 
  } 

 
  // Temporal change to leave primaries untouched at the shallow times
  for (iq=0;iq<nq;iq++) 
    for (it=0;it<(itm-(iq-nmute)*slope);it++)  // 0 should be it0 
      modeltemp[iq][it]=0; 

  if ((plot==2)||(plot==3)){
    //system("rm mutedmodel.su");
    save_gather(modeltemp,nq,q,nt,0.004,"mutedmodel.su");
    system("suxwigb < mutedmodel.su perc=99 key=f2 title=mutedmodel xbox=600 ");
  }

  oper(modeltemp[0],datarec[0],index,0,nt,nh,nq,nsparse);
  /* save the predicted data before applying back the data weighting */
  if ((plot==2)||(plot==3)){
    //system("rm multiples.su");
    save_gather(datarec,nh,h,nt,0.004,"predicted.su");
    system("suxwigb < predicted.su clip=1 key=offset  title=predicted xbox=600 &");
  }
  /* Save the predicted data afted putting back the amplitudes through Wd */

  for (ih=0;ih<nh;ih++) for (it=0;it<nt;it++) datarec[ih][it]*=Wd[ih][it];
  if ((plot==2)||(plot==3)){
    //system("rm multiples.su");
    save_gather(datarec,nh,h,nt,0.004,"multiples.su");
    system("suxwigb < multiples.su clip=1 key=offset  title=multiples xbox=600 &");
  }

  fprintf(stderr,"substracting multiples \n");
  dpd=dot(nh*nt,data[0],datarec[0]);
  dpdp=dot(nh*nt,datarec[0],datarec[0]);
  scale=dpd/dpdp;
  fprintf(stderr,"scale 2===>%f\n",scale);

  for (ih=0;ih<nh;ih++) 
    for (it=0;it<nt;it++)
      data[ih][it]=data[ih][it]-scale*datarec[ih][it];

  if ((plot==2)||(plot==3)){
    //system("rm multiples.su");
    save_gather(data,nh,h,nt,0.004,"primaries.su");
    system("suxwigb < primaries.su clip=1 key=offset  title=primaries ");
  }  
  
  free2float(modeltemp);
  free2float(datarec);
  return;
}


/*****************************************************************************/
/* Functions to predict and plot data corresponding to the model at any stage */
/* With convolution */
void show_predicted_data(void  (*oper) (float *, float *, unsigned int **, int, int, 
				   int, int, int, float *, int), 
		    float **model, float **Wd, unsigned int **index, 
		    int nt, int nh, int nq, int nsparse, float *wavelet, 
		    int nw, float *q, float *h)
{

  float **data=0; 
  int ih,it;

  data=ealloc2float(nt,nh);
  /* Plot the model */
  save_gather(model,nq,q,nt,0.004,"model.su");
  system("suxwigb < model.su key=f2  title=model  xbox=0 &");

  /* Predict the data */
  oper(model[0],data[0],index,0,nt,nh,nq,nsparse,wavelet,nw);
  /* save the predicted data before applying back the data weighting */
  save_gather(data,nh,h,nt,0.004,"predicted0.su");
  system("suxwigb < predicted0.su clip=1 key=offset  title=predicted0 xbox=400 &");
  /* Save the predicted data afted putting back the amplitudes through Wd */
  for (ih=0;ih<nh;ih++) for (it=0;it<nt;it++) data[ih][it]*=Wd[ih][it];
  save_gather(data,nh,h,nt,0.004,"predicted1.su");
  system("suxwigb < predicted1.su clip=1 key=offset  title=predicted1 xbox=800 &");
  
  free2float(data);

  return;
}

/* Without convolution */
void show_predicted_data(void  (*oper) (float *, float *, unsigned int **, int,
					int, int, int, int), 
			 float **model, float **Wd, unsigned int **index, 
			 int nt, int nh, int nq, int nsparse, float *q, float *h)
{

  float **data=0; 
  int ih,it;

  data=ealloc2float(nt,nh);
  /* Plot the model */
  save_gather(model,nq,q,nt,0.004,"model.su");
  system("suxwigb < model.su key=f2  title=model  xbox=0 &");
  /* Predict the data */
  oper(model[0],data[0],index,0,nt,nh,nq,nsparse);
  /* save the predicted data before applying back the data weighting */
  save_gather(data,nh,h,nt,0.004,"predicted0.su");
  system("suxwigb < predicted0.su clip=1 key=offset  title=predicted0 xbox=400 &");
  /* Save the predicted data afted putting back the amplitudes through Wd */
  for (ih=0;ih<nh;ih++) for (it=0;it<nt;it++) data[ih][it]*=Wd[ih][it];
  save_gather(data,nh,h,nt,0.004,"predicted1.su");
  system("suxwigb < predicted1.su clip=1 key=offset  title=predicted1 xbox=800 &");
  
  free2float(data);

  return;
}


