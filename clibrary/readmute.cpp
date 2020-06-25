#include <stdio.h>
#include <stdlib.h>
int main()
{
  char* file="filempicks";
  char* file1="tfileb";
  char* file2="xfileb";
  float x,t;
  int i,k;

  FILE *filep,*filep1,*filep2;
  
  if((filep=fopen(file,"r"))==NULL)
      printf("cannot open file=%s\n",file);
  
  if((filep1=fopen(file1,"w"))==NULL)
    printf("cannot open file=%s\n",file1);
  
  if((filep2=fopen(file2,"w"))==NULL)
    printf("cannot open file=%s\n",file2);
  k=0;
  while((fscanf(filep,"%f%e",&t,&x))==2){
    //for (i=0;i<3;i++){
    //aa=fscanf(filep,"%f%e",&t,&x);
    //printf("aa=%d\n",aa);
    //printf("t=%f\n",t);
    //printf("x=%e\n",x);
    fwrite(&t,sizeof(float),1,filep1);
    fwrite(&x,sizeof(float),1,filep2);
    k+=1;
  }
  fclose(filep);
  fclose(filep1);
  fclose(filep2);
  printf("%d\n",k);
  return(k);

}




