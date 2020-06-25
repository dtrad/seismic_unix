#include <iostream>
#include <stdlib.h>
#include <stdio.h>

using namespace std;
int main(int argc, char **argv){
  FILE * tfp1  = fopen("marmsmooth.time.bin","r");
  FILE * tfp2  = fopen("ppout","w");
  FILE * tfp3  = fopen("ray","w");
  int n2 = 488;
  float* a = new float[n2];
  for (int i=0;i<argc;i++) cout << argv[i] << endl;
  int n1 = atoi(argv[1]);
  for (int i = 0; i< n1;i++){
    fread(a,sizeof(float)*n2,1,tfp1);
    for (int j =0; j< n2;j++) if (a[j] == 999) a[j] = 0;
    fwrite(a,sizeof(float)*n2,1,tfp2);
    if (i%10) cout << a[0] << endl;
    

  }
  fclose(tfp1);
  fclose(tfp2);
  fclose(tfp3);
  delete []a;
  


  return 0;
}
