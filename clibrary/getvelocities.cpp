#include "su.h"

 
void getvelocities(float dt, int nt, int ncdp, float *cdp, float **ovv)
{
  /* This function reads the par file in the command line and computes the vector
  ovv. This vecotr can then be used with interpovv to get the interpolated velocities at the desired cdp
 input :
    float *cdp;	         array[ncdp] of cdps 
    int nt
    float dt
    int ncdp;	 number of cdps specified 
    
 output :
    float **ovv;	 array[ncdp][nt] of sloth (1/velocity^2) functions 

The calling function requires this:

  float **ovv;
  float *cdp;	
  int ncdp;	
  float *velint;

  ncdp = countparval("cdp");
  cdp = ealloc1float(ncdp);
  ovv = ealloc2float(nt,ncdp);
  velint=ealloc1float(nt);

  getvelocities(dt,nt,ncdp,cdp,ovv);

  interpovv(nt,ncdp,cdp,ovv,cdpgather,velint);

  free1float(cdp);
  free2float(ovv);
  free1float(velint);

  */
  int it;
  int icdp;	/* index into cdp array */
  int jcdp;	/* index into cdp array */
  int nvnmo;	/* number of vnmos specified */
  float *vnmo;	/* array[nvnmo] of vnmos */
  int ntnmo;	/* number of tnmos specified */
  float *tnmo;	/* array[ntnmo] of tnmos */
  float tn;     /* temporary time */
  float acdp;	/* temporary used to sort cdp array */
  float *aovv;	/* temporary used to sort ovv array */

  ///////////////////////////////////////////////////////////////
  
  //float **oa1;	/* array[ncdp][nt] of anis1 functions */
  //float **oa2;	/* array[ncdp][nt] of anis2 functions */
  //float **ovv;	/* array[ncdp][nt] of sloth (1/velocity^2) functions */
  //float *cdp;	        /* array[ncdp] of cdps */  
  ///////////////////////////////////////////////////////////////

  /* get velocity functions, linearly interpolated in time */
  if (!getparfloat("cdp",cdp)) cdp[0] = 1;
  //fprintf(stderr,"countparname=%d\n",countparname("vnmo"));
  if (ncdp>0) {
    if (countparname("vnmo")!=ncdp)
      err("a vnmo array must be specified for each cdp");
    if (countparname("tnmo")!=ncdp)
      err("a tnmo array must be specified for each cdp");
  } else {
    ncdp = 1;
    if (countparname("vnmo")>1)
      err("only one (or no) vnmo array must be specified");
    if (countparname("tnmo")>1)
      err("only one (or no) tnmo array must be specified");
  }

  for (icdp=0; icdp<ncdp; ++icdp) {
    nvnmo = countnparval(icdp+1,"vnmo");
    ntnmo = countnparval(icdp+1,"tnmo");
    if (nvnmo!=ntnmo && !(ncdp==1 && nvnmo==1 && ntnmo==0))
      err("number of vnmo and tnmo values must be equal");
    if (nvnmo==0) nvnmo = 1;
    if (ntnmo==0) ntnmo = nvnmo;
    /* equal numbers of parameters vnmo, tnmo, anis1, anis2 */
    vnmo = ealloc1float(nvnmo);
    tnmo = ealloc1float(nvnmo);
    if (!getnparfloat(icdp+1,"vnmo",vnmo)) vnmo[0] = 2000.0;
    if (!getnparfloat(icdp+1,"tnmo",tnmo)) tnmo[0] = 0.0;
    for (it=1; it<ntnmo; ++it)
      if (tnmo[it]<=tnmo[it-1]){
	fprintf(stderr,"Error for #cdp  %d\n",icdp);
	err("tnmo values must increase monotonically");
      }
    for (it=0,tn=0; it<nt; ++it,tn+=dt) 
      intlin(ntnmo,tnmo,vnmo,vnmo[0],vnmo[nvnmo-1],1,&tn,&ovv[icdp][it]);
    
    free1float(vnmo);
    free1float(tnmo);
  }
  
  /* sort (by insertion) sloth and anis functions by increasing cdp */
  for (jcdp=1; jcdp<ncdp; ++jcdp) {
    acdp = cdp[jcdp];
    aovv = ovv[jcdp];
    for (icdp=jcdp-1; icdp>=0 && cdp[icdp]>acdp; --icdp) {
      cdp[icdp+1] = cdp[icdp];
      ovv[icdp+1] = ovv[icdp];
    }
    cdp[icdp+1] = acdp;
    ovv[icdp+1] = aovv;
  } 
}





void getvelocities(float dt, int nt, int ncdp, float *cdp, float **ovv, float **oa1,
		   float **oa2)
{
  /* This function reads the par file in the command line and computes the vectors
  ovv, oa1, and oa2. These 3 vecotrs can then be sued with interpovv to get the interpolated velocities at the desired cdp
 input :
    float **oa1;	 array[ncdp][nt] of anis1 functions 
    float **oa2;	 array[ncdp][nt] of anis2 functions 
    float **ovv;	 array[ncdp][nt] of sloth (1/velocity^2) functions 
    float *cdp;	         array[ncdp] of cdps 
    int nt
    float dt
    int ncdp;	 number of cdps specified 
    
 output :
    float *velint; array[nt] of vel for a particular trace 
    float *oa1t;	 array[nt] of anis1 for a particular trace 
    float *oa2t;	 array[nt] of anis2 for a particular trace 

The calling function requires this:

  float **oa1;
  float **oa2;
  float **ovv;
  float *cdp;	
  int ncdp;	
  float *velint;
  float *oa1t;	
  float *oa2t;	

  ncdp = countparval("cdp");
  cdp = ealloc1float(ncdp);
  ovv = ealloc2float(nt,ncdp);
  oa1 = ealloc2float(nt,ncdp);
  oa2 = ealloc2float(nt,ncdp);
  velint=ealloc1float(nt);
  oa1t = ealloc1float(nt);
  oa2t = ealloc1float(nt);

  getvelocities(dt,nt,ncdp,cdp,ovv,oa1,oa2);

  interpovv(nt,ncdp,cdp,ovv,oa1,oa2,cdpgather,velint,oa1t,oa2t);

  free1float(cdp);
  free2float(ovv);
  free2float(oa1);
  free2float(oa2);
  free1float(oa1t);
  free1float(oa2t);
  free1float(velint);

  */

  int icdp;	/* index into cdp array */
  int jcdp;	/* index into cdp array */
  int nvnmo;	/* number of vnmos specified */
  float *vnmo;	/* array[nvnmo] of vnmos */
  int ntnmo;	/* number of tnmos specified */
  float *tnmo;	/* array[ntnmo] of tnmos */
  int nanis1;	/* number of anis1's specified */
  int nanis2;	/* number of anis2's specified */
  float *anis1;	/* array[nanis1] of anis1's */
  float *anis2;	/* array[nanis2] of anis2's */
  float tn;     /* temporary time */
  float acdp;	/* temporary used to sort cdp array */
  float *aovv;	/* temporary used to sort ovv array */
  float *aoa1;	/* temporary used to sort oa1 array */
  float *aoa2;	/* temporary used to sort oa2 array */
  int i, it;

  ///////////////////////////////////////////////////////////////
  
  //float **oa1;	/* array[ncdp][nt] of anis1 functions */
  //float **oa2;	/* array[ncdp][nt] of anis2 functions */
  //float **ovv;	/* array[ncdp][nt] of sloth (1/velocity^2) functions */
  //float *cdp;	        /* array[ncdp] of cdps */  
  ///////////////////////////////////////////////////////////////


  /* get velocity functions, linearly interpolated in time */
  if (!getparfloat("cdp",cdp)) cdp[0] = 1;
  if (ncdp>0) {
    if (countparname("vnmo")!=ncdp)
      err("a vnmo array must be specified for each cdp");
    if (countparname("tnmo")!=ncdp)
      err("a tnmo array must be specified for each cdp");
    if (countparname("anis1")!=ncdp &&
	countparname("anis1")!=0)
      err("an anis1 array must be specified for each cdp, "
	  "or omitted at all");
    if (countparname("anis2")!=ncdp &&
	countparname("anis2")!=0)
      err("an anis2 array must be specified for each cdp, "
	  "or omitted at all");
  } else {
    ncdp = 1;
    if (countparname("vnmo")>1)
      err("only one (or no) vnmo array must be specified");
    if (countparname("tnmo")>1)
      err("only one (or no) tnmo array must be specified");
    if (countparname("anis1")>1)
      err("only one (or no) anis1 array must be specified");
    if (countparname("anis2")>1)
      err("only one (or no) anis2 array must be specified");    
  }

  for (icdp=0; icdp<ncdp; ++icdp) {
    nvnmo = countnparval(icdp+1,"vnmo");
    ntnmo = countnparval(icdp+1,"tnmo");
    nanis1 = countnparval(icdp+1,"anis1");
    nanis2 = countnparval(icdp+1,"anis2");
    if (nvnmo!=ntnmo && !(ncdp==1 && nvnmo==1 && ntnmo==0))
      err("number of vnmo and tnmo values must be equal");
    if (nanis1!=nvnmo && nanis1 != 0)
      err("number of vnmo and anis1 values must be equal");
    if (nanis2!=nvnmo && nanis2 != 0)
      err("number of vnmo and anis2 values must be equal");
    if (nvnmo==0) nvnmo = 1;
    if (ntnmo==0) ntnmo = nvnmo;
    if (nanis1==0) nanis1 = nvnmo;
    if (nanis2==0) nanis2 = nvnmo;
    /* equal numbers of parameters vnmo, tnmo, anis1, anis2 */
    vnmo = ealloc1float(nvnmo);
    tnmo = ealloc1float(nvnmo);
    anis1 = ealloc1float(nvnmo);
    anis2 = ealloc1float(nvnmo);
    if (!getnparfloat(icdp+1,"vnmo",vnmo)) vnmo[0] = 2000.0;
    if (!getnparfloat(icdp+1,"tnmo",tnmo)) tnmo[0] = 0.0;
    if (!getnparfloat(icdp+1,"anis1",anis1)) 
      for (i=0; i<nvnmo; i++) anis1[i] = 0.0;
    if (!getnparfloat(icdp+1,"anis2",anis2))
      for (i=0; i<nvnmo; i++) anis2[i] = 0.0;
    for (it=1; it<ntnmo; ++it)
      if (tnmo[it]<=tnmo[it-1]){
	fprintf(stderr,"Error for #cdp  %d\n",icdp);
	err("tnmo values must increase monotonically");
      }
    for (it=0,tn=0; it<nt; ++it,tn+=dt) 
      intlin(ntnmo,tnmo,vnmo,vnmo[0],vnmo[nvnmo-1],1,&tn,&ovv[icdp][it]);
    
    for (it=0,tn=0; it<nt; ++it,tn+=dt) 
      intlin(ntnmo,tnmo,anis1,anis1[0],anis1[nanis1-1],1,&tn,&oa1[icdp][it]);
    
    for (it=0,tn=0; it<nt; ++it,tn+=dt) 
      intlin(ntnmo,tnmo,anis2,anis2[0],anis2[nanis2-1],1,&tn,&oa2[icdp][it]);
    
    free1float(vnmo);
    free1float(tnmo);
    free1float(anis1);
    free1float(anis2);
  }

  /* sort (by insertion) sloth and anis functions by increasing cdp */
  for (jcdp=1; jcdp<ncdp; ++jcdp) {
    acdp = cdp[jcdp];
    aovv = ovv[jcdp];
    aoa1 = oa1[jcdp];
    aoa2 = oa2[jcdp];
    for (icdp=jcdp-1; icdp>=0 && cdp[icdp]>acdp; --icdp) {
      cdp[icdp+1] = cdp[icdp];
      ovv[icdp+1] = ovv[icdp];
      oa1[icdp+1] = oa1[icdp];
      oa2[icdp+1] = oa2[icdp];
    }
    cdp[icdp+1] = acdp;
    ovv[icdp+1] = aovv;
    oa1[icdp+1] = aoa1;
    oa2[icdp+1] = aoa2;
  } 

}








