#include <iostream.h>
#include "su.h"

#define TRACE fprintf(stderr,"%s @ %u\n",__FILE__,__LINE__)

class TGather {
private:
  int nx;
  int nt;
  float **d;
  float *t;
  float *x;
  float dt;
  float dx;
  
public:
  TGather(int n1, int n2, float dt, float dx, float t0, float x0);
  ~TGather();
  void Display(){fprintf(stderr,"nt=%d,nx=%d\n",nt,nx);}
  void getsize(int &n1, int &n2){n1=nt;n2=nx;}
  void setx(float x_ix, int ix){x[ix]=x_ix;}
  float getd(int ix, int it) { return(d[ix][it]);}
  float getx(int ix) { return(x[ix]);}
  float gett(int it) { return(t[it]);}
  void plot(char *name, char *key);
  void gettaxis(float *to);
  void getxaxis(float *xo);
  void getdata(float **data);
  void setdata(float **data);
  void gettrace(float *data, int ix);
  void settrace(float *data, int ix);
  void save_gather0(float **d, int nh, float *h, int nt, float dt, char* name);
};



  




