void radon_param(double fmax, double xmin,double xmax,double dx,
double qmin, double qmaxt, double qmax, double dq, int nq, int  rtmethod)
/* Given the field geometry (dx, xmin, xmax and fmax) and the chosen qmin,    
   it computes dq, maximum allowable qmax, nq.
    Hence, the NMO must be adjust such that q < qmax.	
    The dx can be the average or maximum, it is under research.
    rtmethod=1 PRT
    rtmethod=2 LRT 
    Daniel Trad- UBC- 16-2-99
*/
{   
   if (rtmethod==1) {
	  dq= 1/(fmax*(xmax-xmin)*(xmax-xmin)); 
	  qmaxt = 1/(2*fmax*(xmax-xmin)*dx);
	  qmax = 1.5*qmaxt;
	  nq= (int) ((qmax-qmin)/dq) + 1;
    }
    else if(rtmethod==2) {
	  dq= 1/(fmax*(xmax-xmin)); 
	  qmaxt = 1/(fmax*dx);
	  qmax = 1.5*qmaxt;
	  nq= (int) ((qmax-qmin)/dq) + 1;
    }  
    return;
}

	





