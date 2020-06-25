/*Este programa lee un archivo de datos binario obtenido con el equipo
de medicion de MAGNETOTELURICA EMI y lo escribe en ASCII para su posterior
tratamiento matematico.
Estos archivos tienen un header ASCII donde estan especificadas las caracteris-
ticas de registracion (ganancias, frecuencias, filtros, etc.)
Luego siguen los datos agrupados por ventanas de 512 valores con 5
canales cada una.Cada valor esta formado por dos bytes.
Al final de cada canal hay un alimentador de linea y un retorno de carro.
Despues de cada ventana hay un ok para con alimentador de linea y retorno
de carro.*/
/*#include <alloc.h>*/
#include <stdio.h>
#include <stdlib.h>
#include "nrutil.h"

#define win 512

void read_header(int);
void read_data(int *ival);
void write_data(int count3,int lwin,int *ival);
char res;
int count4=0;
FILE *fpt,*fout; /*punteros de archivos de entrada y salida*/
int main(int argc, char **argv)
{
  char inputfile[30],outputfile[30];
  int f,nwin,lwin,nlin,count3;
  int *ival=NULL;
  unsigned char c1,c2;
  unsigned int i1,i2,i,j;

  /*printf("What long of windows do you like?\n"); time series continuos*/
  /*scanf("%d",&lwin);                They can grouped in different long*/

  printf("maximum number of windows\n");
  scanf("%d",&nwin);
  printf("input file? \n");
  scanf("%s",inputfile);
  printf("output file? \n");
  scanf("%s",outputfile);
  printf("number of lines in the header\n");
  scanf("%d",&nlin);
  printf("Do you like with header :y or n\n"); /*answer n if you like plot the series*/
  scanf("%s",&res);                            /*answer y if you like process them*/
  if ((ival = (int*)malloc(sizeof(int)*win*5*1)) == NULL) /*But it is not needed dinamic allocation*/
    {printf("Not enough memory to allocate array\n");exit(1);}
  fpt=fopen(inputfile,"rb");  /*bynary file*/
  fout=fopen(outputfile,"w"); /*text file*/
  if(fpt==NULL) exit(1);
  read_header(nlin);
  for(count3=1;count3<=nwin;++count3)
    {
      read_data(ival);
      waveletfilter(ival,lwin);
      write_data(count3,lwin,ival);
      printf("window %d",count3);
    }
  fclose(fpt);
  fclose(fout);
  return 0;
}

void read_header(int nlin)
/*function for reading the header and put it in text file*/
{
  int count=0;
  char c1;
  while(count<nlin)
    {
      c1=fgetc(fpt);
      printf("%c",c1);
      if ((c1!='\r')&&((res=='y')||(res=='Y'))) fprintf(fout,"%c",c1);
      if (c1=='\n') ++count;
    }
  return;
}


void write_data(int count3,int lwin,int *ival)
     /*function for writing data to text file. Too it puts the number of window
       at top of each one*/
{
  int i,j;
  /*if ((count4%lwin==0)&&((res=='y')||(res=='Y'))) fprintf(fout,"%d\n",count3);*/
  for(i=0;i<win;++i){
      count4++;
      for(j=0;j<5;++j){    /*write five columns */
	fprintf(fout,"%d\t",ival[i+j*win]);
	if(j==4) fprintf(fout,"\n");
      }
  }
  return ;
}


void read_data(int *ival)
     /*function for reading binary data. */
{
  int f,count=0,count2=0;
  unsigned char c1,c2;
  unsigned i1,i2,i,j;
  do {
    ++count;
    if((((count-1)%win)==0)&&(count!=1))   /*end of channel*/
      do{
	c1=fgetc(fpt);
	printf("%c channel",c1);     /*it jumps \n and \r */
      } while (c1!='\n');
    if((((count-1)%(5*win))==0)&&(count!=1)){  /*end of window*/
      do {
	c1=fgetc(fpt);
	printf("%c",c1);
      } while (c1!='\n');
      printf("end of window\n");
      break;
    }
    c1=fgetc(fpt);
    c2=fgetc(fpt);
    i1=(int)c1;
    i2=(int)c2;                        /*cast*/
    ival[count2]=i1*256+i2-32768;      /*computing the integer value*/
    ++count2;
    f=feof(fpt);         /*end of file*/
    if(f){
      fclose(fpt);
      fclose(fout);
      exit(1);
    }
  }while ((!feof(fpt))&&(count2<=(512*5)));
  return ;
}

void waveletfilter(int *ival,lwin)
{
  int i;
  int nw=5;
  float *data;
  data=(float *) malloc(lwin*sizeof(float));
 
  for (i=0;i<nw;i++){
    for (j=0; j< lwin; j++) data[j]=ival[j+nw*i];
    
    
    
  }
  free((float *) data);
  return ;
}






