/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUVCAT: $Revision: 1.18 $ ; $Date: 2011/11/16 23:09:52 $	*/

#include "su.h"
#include "segy.h"
#include <math.h>

/*********************** self documentation **********************/
char *sdoc[] = {
" 								",
" SUCREATEDATA - add data to a geometry                         ",
" sucreatedata < stdin > stdout					",
" if offsetxy=0                                                 ",
" writes in key stae the actual azimuth,   			",
"      		stas the binned azimuth (dazim)			",
" 		offset the true offset.				",
" if offsetxy=1 						",
" writes in stae   = azimuth               			",
"           offset = offsetx                                    ",
" 	    stas   = offsety					",
"                                                               ",
" if cdp=0 it writes a file with midpointx, midpointy and cdp   ",
"                                                               ",
" Required parameters:						",
"        none							",
"                                                               ",
" Example: sucreatedata < input.su > output.su dazim=45         ",
"   sucreatedata < geometry.su | suconv sufile=wavelet.su  ...  ",
"                   | suresamp nt=1501 > data.su    ",
NULL};

/* Credits:
 */
/**************** end self doc ***********************************/

segy tr;
//float myPI = acos(-1.);
void calcOffsetAz(int sx, int sy, int gx, int gy, float* offset, float* azim){
    float offx = (float) (sx - gx);
    float offy = (float) (sy - gy);
    *offset = (float) (sqrt(offx * offx + offy * offy));

    double small = 1.e-15;    
    if (fabs(offx) < small) offx = small;
    *azim = (float) (180./PI*atan2(offy, offx));
    if ((fabs(*azim-180.) < 2)&&(offx<0)) *azim=-*azim;// correction for close to 0.
    //fprintf(stderr,"azim=%f offx=%f offy=%f offset=%f\n",*azim,offx,offy,*offset);
    
    return;
}
int setCdp(int sx, int sy, int gx, int gy, int nx, float dx, float dy){
    // SU uses integers. 
    float midx = (sx + gx)/2.;
    float midy = (sy + gy)/2.;
    int inl = (int) ((midx-0.)/dy) + 1;
    int xl  = (int) ((midy-0.)/dx) + 1;
    int cdp = ((inl-1)*nx)+xl;
    
    return cdp;
}

int createtrace(segy* tr, int sample, int makezero, float amp, float azim){
  
  if (makezero) memset(tr->data,0,tr->ns*sizeof(float)); 
  tr->data[sample]=amp;

  return 1;
}

int deviation(float azim, float azimaxis, float scale){
  int dev = (int) (cos( 2.*(azim+azimaxis) * PI/180. ) * scale);
  return (dev);
}

int main(int argc, char **argv){

  	FILE *fp1 = NULL;	/* file pointer for first file		*/
	int itr = 0;	/* number of trace being processed	*/
	/* Initialize */
	initargs(argc, argv);
	requestdoc(1); /* two file args required */
	int doffset=50;
	float dazim=30;
	int offsetxy=0; // if 1 apply offsetxy binning instead;
	int verbose = 0;
	float azdir = 0; // anisotropy direction
	float scale= 1; // scale of perturbation
	if (!getparint("doffset",&doffset)) doffset=0;
	if (!getparfloat("dazim",&dazim)) dazim=45;
	if (!getparint("offsetxy",&offsetxy)) offsetxy=0;
	if (!getparint("verbose",&verbose)) verbose=0;
	if (!getparfloat("azdir",&azdir)) azdir=110.;
	if (!getparfloat("scale",&scale)) scale=5.;
        checkpars();

	/* Open two files given as arguments for reading */
	fp1 = efopen("list", "w"); 
	/* Loop over the traces */
	float azim, offset;
        int az;

	gettr(&(tr));
	// write a value on the data to plot azimuth
	int sample1 = (int) ((0.9 * 1000000) / tr.dt );
	int sample2 = (int) ((1.05 * 1000000) / tr.dt );
	erewind(stdin);
        while (gettr(&tr)){
            calcOffsetAz(tr.sx,tr.sy,tr.gx,tr.gy,&offset,&azim);
            az = (int) azim;
            tr.stae = az;           

	    if (!offsetxy){
	      if ( doffset){
		float offsetnew = (int) ((int) ((offset+doffset/2)/doffset)*doffset);
		if (verbose) fprintf(stderr,"%f %f ",offsetnew,offset);
		tr.offset = offsetnew;
		
	      }
	      else tr.offset=offset;
	      // to bin in azimuth with decimal degree multiply all by 10
	      int dazimb = (int) (dazim*10);
	      int azimb= (int) (azim*10);

	      tr.stas   = (int) (0.1* ((int) ((azimb-dazimb/2)/dazimb)*dazimb));
	      
	      if (verbose) fprintf(stderr,"azim=%d binned=%d  ",tr.stae,tr.stas);
	    }
	    else{
	      tr.offset=tr.sx-tr.gx; // use same keywords to store values. 
	      tr.stas=tr.sy-tr.gy;
	    }
	    // calculate a cmp number 
	    //if (!tr.cdp){
	    if (1){
	      float mx = (tr.sx + tr.gx)/2.;
	      float my = (tr.sy + tr.gy)/2.;
	      tr.cdp = (int) (my / 30) * 1000 + (int) (mx / 30) + 1;
	      if (verbose) fprintf(fp1,"%f %f %d\n",mx,my,tr.cdp);
	    }




	    int dev = deviation(azim,azdir,scale);
	    //	    fprintf (stderr,"azim=%f dev=%f \n",azim,fabs(azim-45));
	    if (verbose) fprintf(stderr,"moveout = %d %d dev = %d \n",sample1,sample2, dev);

	    //createtrace(&tr,sample,1,1.);
	    createtrace(&tr,sample1+dev,1,1.,azim);
	    createtrace(&tr,sample2+dev/2,0,1.,azim);

	    
	    //tr.data[sample+dev] = 1;//az/180.; 
	    //tr.data[sample] = 1;
            //fprintf(stderr,"azin=%d azout=%d \n",tr.stae,tr.stas);
            puttr(&tr); 
            ++itr;
	}
	fclose(fp1);

	for (int i=-180;i<=180;i+=10){
	  float azimuth=i;
	  int dev = deviation(azimuth,azdir,scale);
	  //	  int dev = (int) (((fabs(azimuth+110))/(180.))*20.);
	  fprintf(stderr,"azimuth - deviation = %f %d\n",azimuth,dev);
	}

	return(CWP_Exit());
}
