#include <iostream.h>
#include <fstream.h>
void main(char file[100])
{
  char* file1="tfileb";
  char* file2="xfileb";
  float x,t;
  int i;
  ifstream fin;
  ofstream fout1;
  ofstream fout2;

  fin.open(filename);
  fout1.open("tfile",ios:bin);
  fout2.open("xfile",ios:bin);
  do{
    fin >> t;
    fin >> x;
    fout << 
  
  
FILE *filep,*filep1,*filep2;

 if((filep=fopen(file,"r"))==NULL)
      printf("cannot open file=%s\n",file);

if((filep1=fopen(file1,"w"))==NULL)
      printf("cannot open file=%s\n",file1);

if((filep1=fopen(file2,"w"))==NULL)
      printf("cannot open file=%s\n",file2);

 while((fscanf(filep,"&t,&x"))==2){
   fwrite(t,sizeof(float),1,filep1);
   fwrite(x,sizeof(float),1,filep2);
 }
fclose(filep);
fclose(file1);
fclose(file2);
}

