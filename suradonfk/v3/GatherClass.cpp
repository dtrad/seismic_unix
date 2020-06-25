#include "GatherClass.hpp"

#include <vector>
#include "Alloc.hpp"
//#include "HeapSortUtl.hpp"

#ifndef MARK
#define MARK printf("%s @ %u\n",__FILE__,__LINE__)
#endif


#define debug 1

// Note Gather class has indexes from 0 to Nx-1.
// Since vega attributes are from 1 to N it may be important to
// change gatherClass as well. Not decided.

using namespace std;


GatherClass::GatherClass(){
  cout << "Null Gather Constructor\n";
  myNx=0;
  myNt=0;
  myXaxis=0;
  myData=0;
  myOwnData=0;
	  
}
     
/** Constructors and destructors */
GatherClass::GatherClass(int Nt, int Nx){
  if (debug){
    cout << "Gather Constructor with ";
    cout << "Nx=" << Nx << ", " ;
    cout << "Nt=" << Nt << endl;
  }
  myNx=Nx;
  myNt=Nt;
  myOwnData=1;
  new1d(myXaxis,myNx);
  new2d(myData,myNx,myNt);

}

/** Constructor the duplicates memmory */
GatherClass::GatherClass(int Nt, int Nx, float Dt, float T0, float* Xaxis){

  if (debug) cout << " Gather constructor own data. (X:float)  \n ";
  myNt=Nt;
  myNx=Nx;
  myT0=T0;
  myDt=Dt;

  new1d(myXaxis,myNx);
  setXaxis(Xaxis);

  new2d(myData,myNx,myNt);
  setDataToZero();

  myOwnData=1;
}

/** The same but receives a pointer to double for offset */
GatherClass::GatherClass(int Nt, int Nx, float Dt, float T0, double* Xaxis){
  if (debug) cout << " Gather constructor own data. (X:double) \n ";
  myNt=Nt;
  myNx=Nx;
  myT0=T0;
  myDt=Dt;


  if (0) cout << "myNt=" << myNt << ", myNx=" << myNx << ", myT0=" << myT0 << endl; 
  new1d(myXaxis,myNx);
  setXaxis(Xaxis);
  new2d(myData,myNx,myNt);
  setDataToZero();
  myOwnData=1;
}
     
GatherClass::~GatherClass(){
  if (myOwnData){
    if (debug) cout << "Gather Destructor" << endl;
    del2d(myData);
    del1d(myXaxis);
	       
  }
  else {
    if (debug) cout << "Gather Destructor without data " << endl; }
	  
}
     
/* Copy Constructor */
GatherClass& GatherClass::operator=(const GatherClass &rhs){
	  
  cout << " Copy constructor \n" ;
  if (this==&rhs) return *this;

  myNx=rhs.getNx();
  myNt=rhs.getNt();
  myDt=rhs.getDt();
  myT0=rhs.getT0();
  myOwnData=1;

  setData(rhs.getData());
  setXaxis(rhs.getXaxis());
	  
  return *this;
}
     
/** Data accessors */
void GatherClass::setDataToZero(){
	  
  for (int ix=0;ix<myNx;ix++)
    memset( (void *) myData[ix], (int) '\0', myNt *sizeof(float) );
}

void GatherClass::setData(float** Data){
  for (int ix=0; ix<myNx; ix++) 
    memcpy(myData[ix],Data[ix],myNt*sizeof(float));
}

//      void GatherClass::setData(float** Data, int buffer){
// 	  // This member function adds a buffer with zeroes at
// 	  // the beginning and the end to get rid of wrap
// 	  // Note that in this case myNt = nDataSamples + 2 * buffer
// 	  for (int ix=0; ix<myNx; ix++){ 
// 	       memset(myData[ix],(int) '\0',buffer*sizeof(float));
// 	       memcpy(&myData[ix][buffer],Data[ix],(myNt-buffer)*sizeof(float));

// 	  }
//      }


     
void GatherClass::setTrace(float* Trace, int ix, int ftime){
  memcpy(myData[ix],&Trace[ftime],myNt*sizeof(float));
}

void GatherClass::setXaxis(float* Xaxis){
  memcpy(myXaxis,Xaxis,myNx*sizeof(float));
}
     
void GatherClass::setXaxis(double* Xaxis){
  for (int ix=0;ix<myNx;ix++){
    //cout << ix << " = " << Xaxis[ix] << endl
    myXaxis[ix]=(float) Xaxis[ix];
  }
}
void GatherClass::setXaxis(float Xaxis, int ix){
  myXaxis[ix]=Xaxis;
}
     
/** copy one trace to a 1d array */
void GatherClass::copyTraceTo(float* Trace, int ix, int ftime ){
  if (Trace==myData[ix]) return;
  else memcpy(&Trace[ftime],myData[ix],myNt*sizeof(float));
}


/** Miscelaneous: saving, plotting, debugging etc. */
void GatherClass::saveData(const char *s){
  FILE* fp;
  if ((fp=fopen(s,"w"))==NULL){ 
    fprintf(stderr,"Cannot open fp");
    return;
  }
  for (int ix=0;ix< myNx ;ix++ )  
    fwrite(myData[ix],sizeof( float ), myNt, fp);
	  
  fclose(fp);
	  
  return;
}
/** Miscelaneous: saving, plotting, debugging etc. */
void GatherClass::saveSortedData(const char *s, int *index){
  // If the data are not sorted in the x axis
  // save them by using an index vector corresponding to
  // the sorted elements.
	  
  FILE* fp;
	  
  if ((fp=fopen(s,"w"))==NULL){ 
    fprintf(stderr,"Cannot open fp");
    return;
  }
  for (int ix=0;ix< myNx ;ix++ )  
    fwrite(myData[index[ix]],sizeof( float ), myNt, fp);
	  
  fclose(fp);
	  
  return;
}
void GatherClass::plotData(char *filename){
  char buf[80];
	  
  saveData(filename);
  sprintf(buf,"ximage < %s n1=%d n2=%d perc=99 title=%s\n",filename,myNt,myNx,filename);
  system(buf);
  return;
}

void GatherClass::plotSortedData(char *filename){
  // If the data are not sorted in the x axis
  // save them by using an index vector corresponding to
  // the sorted elements.
	   
  char * buf  = new char[myNx*10];
  char * buf2 = new char[myNx*10];
  char * buf3 = new char[10];
  buf[0]  = '\0';
  buf2[0] = '\0';
  buf3[0] = '\0';
  int itemp;
  int i;
	  
  vector<int> index(myNx);
  vector<float> xaxis(myNx);
	  
  for (i=0; i < myNx; i++){
    index[i]=i;
    xaxis[i]=myXaxis[i];
  }
	  
  //vUtl::heapSort(myNx,&xaxis.front(),&index.front());

  // write a string with the offsets
  for (i=0; i < myNx-1; i++){
    itemp = (int) xaxis[i];
    sprintf(buf3,"%d,",itemp);
    strcat(buf2,buf3);
  }
  itemp = (int) xaxis[myNx-1];
  sprintf(buf3,"%d",itemp);
  strcat(buf2,buf3);
	  
  //saveSortedData(filename,&index.front());
  sprintf(buf,"xwigb < %s n1=%d n2=%d perc=99 title=%s x2=%s\n",filename,
	  myNt,myNx,filename,buf2);
	  
  system(buf);
  cout << buf ;
	  
  delete buf;
  delete buf2;
  delete buf3;
  return;
}
      
void GatherClass::Display(){
  cout << " nt = " << myNt << ", nx = " << myNx << endl ;
}

      

