#ifndef GatherClass_HPP
#define GatherClass_HPP

#include <iostream>

class GatherClass {
	  
	  
public:
  /** Constructors and destructors */
  /** Null gather constructor */
  GatherClass();

  /** Simple gather constructor */
  GatherClass(int Nt, int Nx);

  /** Constructor that allocates new memory and set data to zero */
  GatherClass(int Nt, int Nx, float dt, float t0, float* x);
	  
  /** The same constructor but the offset is given as double (vega uses double) */
  GatherClass(int Nt, int Nx, float dt, float t0, double* x);
	  
  /** Destructor */
  ~GatherClass();

  /** Data accessors */

  /** gets */
  float** getData()const {return myData;}
  float* getXaxis()const {return myXaxis;}
  float getXaxis(int ix)const{return myXaxis[ix];}
  int getNx()const {return myNx;} 
  int getNt()const {return myNt;}
  bool getOwnData()const {return myOwnData;}
  float getDt()const {return myDt;}
  float getT0()const {return myT0;}
  float* getTrace(int ix){return myData[ix];}
  void copyTraceTo(float* Trace, int ix, int ftime );
	  
  /** sets */
  void setDataToZero();
  void setData(float** Data);
  void setTrace(float* Trace, int ix, int ftime  );
  void setXaxis(float* Xaxis);
  void setXaxis(double* Xaxis);
  void setXaxis(float  Xaxis, int i);

  /** Miscelaneous: saving, plotting, debugging etc. */
  void saveData(const char *s);
  void saveSortedData(const char *s, int * index);
  void plotData(char *filename);
  void plotSortedData(char *filename);
	  

  /** Possible additions */
  void Display();
	  
  void getSize(int &n1, int &n2){n1=myNt;n2=myNx;}
  float getDataAt(int ix, int it) { return(myData[ix][it]);}

  GatherClass& operator=(const GatherClass &rhs);
  //GatherClass& operator=(const GatherClass *rhs);
	  
  /** Add operator=
   *  Add operator+
   *  Add operator-
   *  Add operator*
   *  Add operator<<
   */

private:
  int myNx;
  int myNt;
  float **myData;
  float *myXaxis;
  float myDt;
  float myT0;
  bool  myOwnData;
};


#endif
