/* Copyright (c) University of British Columbia, 1999.*/
/* All rights reserved.                       */
/* SUEOMIG:  $Date: March 1999  */

#include "su.h"
#include "segy.h"
#include "clibrarytd.h"
#include "header.h"
#include <math.h>

/*********************** self documentation **********************/
char *sdoc[] = {
  " 	   								",
  " SUREADVELI - writes veloctiy sufile from parfile                    ",
  "                                                                     ", 
  "	                                                          	",
  " 	   								",
  " suredceli< stdin > stdout [optional parameters]          		",
  " 									",
  " cdpmin=0        Fisrt CDP in meters                                 ",
  " cdpmax=100      Last CDP 						",
  " par=            file with stacking velocities and time as obtained  ",
  "                 from Velan                                          ",
  "                                                                     ",
  " verbose=0       =1  Extra Information                               ",
  "                                      				",
  NULL};
/* Credits:
 *	Daniel Trad.
 * Trace header fields accessed: ns, cdp, dt, offset
/**************** end self doc ***********************************/
// Function Prototypes

void interpovv(int nt, int ncdp, float *cdp, float **ovv, float **oa1, 
	       float **oa2, float cdpt, float *ovvt, float *oa1t, float *oa2t);

//////////////////////////////////////////////////////////////////////


int main(int argc, char **argv)
{
  segy tr; 
  int  nt;
  int  nx;         // number of output cdps
  float dt,dx;
  int j, i, k, ix; // General counters 
  register int it;
  float *x;      /* axis for CDPs */
  float  *t;     // time axis for input and output

  /* Velocity */
  int ncdp;	/* number of cdps specified */
  float *cdp;	/* array[ncdp] of cdps */
  int icdp;	/* index into cdp array */
  int jcdp;	/* index into cdp array */
  int nvnmo;	/* number of vnmos specified */
  float *vnmo;	/* array[nvnmo] of vnmos */
  int ntnmo;	/* number of tnmos specified */
  float *tnmo;	/* array[ntnmo] of tnmos */
  float **ovv;	/* array[ncdp][nt] of sloth (1/velocity^2) functions */
  float *velint;/* array[nt] of vel for a particular trace */
  int nanis1;	/* number of anis1's specified */
  int nanis2;	/* number of anis2's specified */
  float *anis1;	/* array[nanis1] of anis1's */
  float *anis2;	/* array[nanis2] of anis2's */
  float **oa1;	/* array[ncdp][nt] of anis1 functions */
  float **oa2;	/* array[ncdp][nt] of anis2 functions */
  float *oa1t;	/* array[nt] of anis1 for a particular trace */
  float *oa2t;	/* array[nt] of anis2 for a particular trace */
  float tn;     /* temporary time */
  float acdp;	/* temporary used to sort cdp array */
  float *aovv;	/* temporary used to sort ovv array */
  float *aoa1;	/* temporary used to sort oa1 array */
  float *aoa2;	/* temporary used to sort oa2 array */
  float  cdpmin;          //limit for the output
  float  cdpmax;          //limit for the output
  int verbose;		/* flag for echoing info		*/

  // Initialize 
  initargs(argc, argv);
  requestdoc(1);

  // Get parameters 
  
  if (!getparfloat("cdpmin", &cdpmin))  cdpmin = 0;
  if (!getparfloat("cdpmax", &cdpmax))  cdpmax = 1000;
  if (!getparfloat("dxcdp", &dxcdp))  dxcdp = 10;
  if (!getparint("verbose", &verbose)) verbose = 0;  
  if (!getparfloat("dt",&dt)) dt=0.004;
  if (!getparint("nt", &nt)) nt = 0; 
 
  nx=(int) floor((cdpmax-cdpmin)/dxcdp + 1);

  // Because this implementation computes ncsp csps together
  // we need to modify nx such that mod(nx/ncsp) =0 

  if (verbose) fprintf(stderr,"***sureadveli => par file -> sufile\n");

  /* get velocity functions, linearly interpolated in time */
  ncdp = countparval("cdp");
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

  cdp = ealloc1float(ncdp);
  if (!getparfloat("cdp",cdp)) cdp[0] = tr.cdp;
  ovv = ealloc2float(nt,ncdp);
  oa1 = ealloc2float(nt,ncdp);
  oa2 = ealloc2float(nt,ncdp);
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
      if (tnmo[it]<=tnmo[it-1])
	err("tnmo values must increase monotonically");

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

  /* allocate workspace */
 
  oa1t = ealloc1float(nt);
  oa2t = ealloc1float(nt);

  // Allocate memory for data and model
    
  if ((x=ealloc1float(nx))==NULL)
    fprintf(stderr,"***Sorry, space for x could not be allocated\n");
  
  if ((velint=ealloc1float(nt))==NULL)
    fprintf(stderr,"***Space for velint could not be allocated\n");

  /* Create axis for output cpds, equivalent offset and time */
   
  for(ix=0;ix<nx;ix++) x[ix]=cdpmin+dx*ix;

  /* Loop to compute one CSP gather */
  for (ix=0;ix<nx;ix+=ncsp){
    fprintf(stderr,"cdp=%f\n",x[ix]);   
    /* compute new square slowness and anis function */
    interpovv(nt,ncdp,cdp,ovv,oa1,oa2,x[ix],velint,oa1t,oa2t); 
    /* Map this trace to the current csp */
    memcpy((void *) tr.data,(const void *) velint,nt*sizeof(float));

    tr.cdp=(int) x[ix]; // front of the CSP
    tr.dt=(int) (dt*1e6);       
    tr.ntr=nx;
    tr.ns=nt;
    tr.tracl=ix;
    tr.tracr=ix;
    tr.sx=(int) x[ix];
    tr.gx=(int) x[ix];
    puttr(&tr);    
  }


  free1float(cdp);
  free2float(ovv);
  free2float(oa1);
  free2float(oa2);
  free1float(oa1t);
  free1float(oa2t);
  free1float(velint);
  free1float(x);
  return EXIT_SUCCESS;
}


void interpovv (int nt, int ncdp, float *cdp, float **ovv, float **oa1, 
	float **oa2, float cdpt, float *ovvt, float *oa1t, float *oa2t)
{
	static int indx=0;
	int it;
	float a1,a2;

	/* if before first cdp, constant extrapolate */
	if (cdpt<=cdp[0]) {
		for (it=0; it<nt; ++it) {
			ovvt[it] = ovv[0][it];
			oa1t[it] = oa1[0][it];
			oa2t[it] = oa2[0][it];
		      };
	
	/* else if beyond last cdp, constant extrapolate */
	} else if (cdpt>=cdp[ncdp-1]) {
		for (it=0; it<nt; ++it) {
			ovvt[it] = ovv[ncdp-1][it];
			oa1t[it] = oa1[ncdp-1][it];
			oa2t[it] = oa2[ncdp-1][it];
		      };
	
	/* else, linearly interpolate */
	} else {
		xindex(ncdp,cdp,cdpt,&indx);
		a1 = (cdp[indx+1]-cdpt)/(cdp[indx+1]-cdp[indx]);
		a2 = (cdpt-cdp[indx])/(cdp[indx+1]-cdp[indx]);
		for (it=0; it<nt; ++it) {
			ovvt[it] = a1*ovv[indx][it]+a2*ovv[indx+1][it];
			oa1t[it] = a1*oa1[indx][it]+a2*oa1[indx+1][it];
			oa2t[it] = a1*oa2[indx][it]+a2*oa2[indx+1][it];
		      };
	}
}





















