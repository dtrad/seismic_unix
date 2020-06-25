#ifndef dataClass_HPP
#define dataClass_HPP

class dataClass{
public:
  dataClass();
  dataClass(int nt, int nh, int dt);
private:
  bool  myOwnData;
  int   myNt;
  int   myNx;
  float myDt;
  float **myData;
  float *myXaxis;
  
  
};

#endif
