/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUVCAT: $Revision: 1.18 $ ; $Date: 2011/11/16 23:09:52 $	*/

#include "su.h"
#include "segy.h"
#include <math.h>

/*********************** self documentation **********************/
char *sdoc[] = {
    " 								",
    " SUBIN4D -  bin shot receiver coordinates for plotting     ",
    " subin4d < data1 > data2 		         		",
    " 								",
    " Required parameters:					",
    "        none						",
    " binsx=                                                    ",
    " binsy=                                                    ",
    " bingx=                                                    ",
    " bingy=                                                    ",
    " 								",
    NULL
};

/* Credits:
 */
/**************** end self doc ***********************************/

segy tr;
int main(int argc, char **argv) {


    /* Initialize */
    initargs(argc, argv);
    requestdoc(1); /* two file args required */
    int binsx, binsy, bingx, bingy;
    int verb=0;
    if (!getparint("binsx", &binsx)) binsx = 1;
    if (!getparint("binsy", &binsy)) binsy = 1;
    if (!getparint("bingx", &bingx)) bingx = 1;
    if (!getparint("bingy", &bingy)) bingy = 1;
    if (!getparint("verb", &verb)) verb = 0;

    checkpars();
    int a,b,c,d;
    while(gettr(&tr)){
      a=tr.sx;b=tr.sy;c=tr.gx;d=tr.gy;
      if (verb) fprintf(stderr,"before %d %d %d %d",tr.sx,tr.sy,tr.gx,tr.gy);
      tr.sx=((tr.sx+binsx/2)/binsx)*binsx;
      tr.sy=((tr.sy+binsy/2)/binsy)*binsy;
      tr.gx=((tr.gx+bingx/2)/bingx)*bingx;
      tr.gy=((tr.gy+bingy/2)/bingy)*bingy;
      if (verb) fprintf(stderr,"--> %d %d %d %d \n",tr.sx,tr.sy,tr.gx,tr.gy);
      if (verb) fprintf(stderr,"%d %d %d %d \n",a-tr.sx,b-tr.sy,c-tr.gx,d-tr.gy);
      puttr(&tr);
    }
    

    return (CWP_Exit());
}
