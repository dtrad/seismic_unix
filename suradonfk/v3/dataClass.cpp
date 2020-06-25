#include  <iostream>
#include "dataClass.hpp"
#include "alloc.hpp"

dataClass::dataClass(){
  
  std::cout << "Null Gather Constructor\n";
  myNx=0;
  myNt=0;
  myDt=0;
  myXaxis=0;
  myData=0;
  myOwnData=0;
}

dataClass::dataClass(int nt, int nh, int dt){
  
}
