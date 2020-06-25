/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUVCAT: $Revision: 1.18 $ ; $Date: 2011/11/16 23:09:52 $	*/

#include "su.h"
#include "segy.h"
#include <math.h>

/*********************** self documentation **********************/
char *sdoc[] = {
    " 								",
    " SUCAT -  append one data set  				",
    " sugeomfix < data1 > data2 				",
    " 								",
    " Required parameters:					",
    "        none						",
    " shift=1 shift to origin 0,0                               ",
    " rotation=1 rotate parallel to axis x                      ",
    " first=1    first tracl number at beginning of rcvr line   ",
    " last=100   last  tracl number at end of rcvr line         ",
    " 								",
    NULL
};

/* Credits:
 */
/**************** end self doc ***********************************/

segy tr;
float myPI;
void select(int* px, int* py, int x, int y, int id, int value) {
    if (id == value) {
        *px = x;
        *py = y;
        fprintf(stderr, "found id=%d \n",id);
    }
}

float findAngle(int* x, int* y, int* id, int ntr, int first, int last){

    int gxlb,gylb,gxle,gyle;
    gxlb=gylb=gxle=gyle=0;
    for (int itr=0;itr<ntr;itr++){        
        select(&gxlb,&gylb,x[itr],y[itr],id[itr],first);
        select(&gxle,&gyle,x[itr],y[itr],id[itr],last);
    }

    float deltax = gxle - gxlb;
    float deltay = gyle - gylb;
    float angle = atan2(deltay, deltax);
    fprintf(stderr, "deltax=%f, deltay=%f\n", deltax, deltay);
    return angle;
}

void applyShift(int* sx, int* sy, int* gx, int* gy, int ntr){
    int xmin = sx[0];
    int ymin = sy[0];
    for (int itr=0;itr<ntr;itr++){
        if (xmin > sx[itr]) xmin = sx[itr];
        if (ymin > sy[itr]) ymin = sy[itr];
        if (xmin > gx[itr]) xmin = gx[itr];
        if (ymin > gy[itr]) ymin = gy[itr];
    }
    for (int itr=0;itr<ntr;itr++){
        sx[itr] -= xmin;
        sy[itr] -= ymin;
        gx[itr] -= xmin;
        gy[itr] -= ymin;
    }
    fprintf(stderr, "read=%i xmin=%i ymin=%i \n", ntr, xmin, ymin);
 
}

void rotate(int* x, int* y, int ntr, float sina, float cosa){
    int xb, yb;
    for (int i=0;i<ntr;i++){
        xb=x[i];
        yb=y[i];
        x[i]=(xb*cosa - yb*sina);
        y[i]=(xb*sina + yb*cosa);
    }
}


int main(int argc, char **argv) {
    myPI=acos(-1.);
    int itr = 0; /* number of trace being processed	*/
    /* Initialize */
    initargs(argc, argv);
    requestdoc(1); /* two file args required */
    int doffset = 50;
    int dazim = 30;
    int offsetxy = 0; // if 1 apply offsetxy binning instead;
    int rotation = 0;
    int shift = 0;
    
    int first; // first tracl number to estimate angle
    int last;  // last tracl number to estimate angle;
    int maxtracl; // max tracl to process;
    float angle;  // angle can be given or calculated
    int scalecoord; // scale shot/receiver
    if (!getparint("doffset", &doffset)) doffset = 50;
    if (!getparint("dazim", &dazim)) dazim = 30;
    if (!getparint("offsetxy", &offsetxy)) offsetxy = 0;
    if (!getparint("rotation", &rotation)) rotation = 0;
    if (!getparint("shift", &shift)) shift = 0;
    if (!getparint("first",&first)) first=1;
    if (!getparint("last",&last))   last=100;
    if (!getparint("maxtracl",&maxtracl)) maxtracl=100;
    if (!getparfloat("angle",&angle)) angle = 0;
    if (!getparint("scalecoord",&scalecoord)) scalecoord= 1;
    checkpars();
    if (!gettr(&tr)) fprintf(stderr, "can't open stdin file\n");
    int ntr=tr.ntr;
    if (!ntr) CWP_Exit();
    int* sx = ealloc1int(ntr);
    int* sy = ealloc1int(ntr);
    int* gx = ealloc1int(ntr);
    int* gy = ealloc1int(ntr);
    int* id = ealloc1int(ntr);
    if ((shift)||(rotation)){
        rewind(stdin);
        itr=0;
        while(gettr(&tr)){            
            sx[itr]=tr.sx;
            sy[itr]=tr.sy;
            gx[itr]=tr.gx;
            gy[itr]=tr.gy;
            id[itr]=tr.tracl; // hardcoded to tracl, can be changed here
            itr++;
        }
        if (itr<ntr) tr.ntr = itr;
        else if (itr>ntr) err("tr.ntr set too small\n");        
    }
    
    fprintf(stderr,"ntr=%d\n",ntr);
    
    /* Open two files given as arguments for reading */
    /* Loop over the traces */
    // first scan for sx,sy,gx,gy

    if (shift) applyShift(sx,sy,gx,gy,ntr);
    
    if (rotation) {
      if (!angle) findAngle(gx,gy,id,ntr,first,last);
      else fprintf(stderr,"angle given = %f\n",angle);
      angle=angle*myPI/180.;
      float sina = sin(-angle);
      float cosa = cos(-angle);
        
      fprintf(stderr, "angle =%f sina=%f cosa=%f \n", 180. / myPI*angle, sina, cosa);
      rotate(sx,sy,ntr,sina,cosa);
      rotate(gx,gy,ntr,sina,cosa);
      // check angle after rotation
      angle= findAngle(gx,gy,id,ntr,first,last);
      fprintf(stderr, "angle =%f sina=%f cosa=%f \n", 180. / myPI*angle, sina, cosa);
    }
    
    // if shift again repeat first step;
    if (shift) applyShift(sx,sy,gx,gy,ntr);

    // output traces with new headers;
    rewind(stdin);
    itr = 0;
    while (gettr(&tr)) {
        tr.sx=sx[itr]/scalecoord;
        tr.sy=sy[itr]/scalecoord;
        tr.gx=gx[itr]/scalecoord;
        tr.gy=gy[itr]/scalecoord;

        puttr(&tr);
        itr++;
    }
    free1int(sx);
    free1int(sy);
    free1int(gx);
    free1int(gy);
    free1int(id);
    return (CWP_Exit());
}
