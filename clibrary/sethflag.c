/*********************** self documentation **********************/
/*****************************************************************************
SETHFLAG   Given two data sets, one with offset h, another with offset h2, this function sets the offset flag hflag = 1 if there is already a trace with offset h2[ih]. This means that the original trace with offset h[ih] will be output instead of the new trace h2[ih].  

void sethflag(float *hflag, int nhflag, float *h, int nh, float *h2, int nh2)

 Daniel Trad - June 2000 - UBC 

******************************************************************************/

void sethflag(int *hflag, int nhflag, float *h, int nh, float *h2, int nh2)
{
  int ih,ih2;
  float under=0.999;
  float over=1.001;

  for (ih2=0;ih2<nh2;ih2++){
    hflag[ih2]=0;
    for (ih=0;ih<nh;ih++){
      if (h[ih]>=0) 
	if ((h2[ih2]>under*h[ih])&&(h2[ih2]<over*h[ih])) hflag[ih2]=1;
      if (h[ih]<0) 
	if ((h2[ih2]<under*h[ih])&&(h2[ih2]>over*h[ih])) hflag[ih2]=1;
    }
  }
}

void sethflag2(int *hflag, int nhflag, float *h, int nh, float *h2, int nh2, float under, float over)
{
  int ih,ih2;

  for (ih2=0;ih2<nh2;ih2++){
    hflag[ih2]=0;
    for (ih=0;ih<nh;ih++){
      if (h[ih]>=0) 
	if ((h2[ih2]>under*h[ih])&&(h2[ih2]<over*h[ih])) hflag[ih2]=1;
      if (h[ih]<0) 
	if ((h2[ih2]<under*h[ih])&&(h2[ih2]>over*h[ih])) hflag[ih2]=1;
    }
  }
}

