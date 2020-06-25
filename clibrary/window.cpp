#include "su.h"
#include "dan.h"
/************************************************************************** 

point a submatrix to a[a2f..a2f][a1f..a1l] from Press et al. 1988
However, this submatrix CANNOT be used in functions that expect a vector
(as the original matrix can)
 
Hence use submatrix only for array operations
Daniel Trad January 2001
***************************************************************************/

float **window(float **a, int a2f, int a2l, int a1f, int a1l)
/* point a submatrix to a[a2f..a2l][a1f..a1l] */
{
	int i,j,n2=a2l-a2f+1,n1=a1f;

	float **m=NULL;

	/* allocate array of pointers to rows */
	m=(float **) malloc((size_t) ((n2)*sizeof(float*)));
	if (!m) err("allocation failure in window()");

	/* set pointers to rows */
	for(i=a2f,j=0;i<=a2l;i++,j++) m[j]=a[i]+n1;

	/* return pointer to array of pointers to rows */
	return m;
}

complex **window(complex **a, int a2f, int a2l, int a1f, int a1l)

{
	int i,j,n2=a2l-a2f+1,n1=a1f;

	complex **m=NULL;

	/* allocate array of pointers to rows */
	m=(complex **) malloc((size_t) ((n2)*sizeof(complex*)));
	if (!m) err("allocation failure in window()");

	/* set pointers to rows */
	for(i=a2f,j=0;i<=a2l;i++,j++) m[j]=a[i]+n1;

	/* return pointer to array of pointers to rows */
	return m;
}


  // Test for window
  //plotgather(data[0],nt,nh,"data");
  //float **dataw=window(data,41,60,101,400);
  //for (ih=2;ih<5;ih++) for (it=0;it<300;it++) dataw[ih][it]*=0;
  //zero_array(dataw,3,300);
  //plotgather(dataw[0],nt,20,"dataw");
  //plotgather(data[0],nt,nh,"data_after");


void window_cpy(float **a, int a2f, int a2l, int a1f, int a1l, float **m, int go)
{
  
  int i1,i2; // indexes for data
  int j2,j1; // indexes for window

  
  if (go){
    for (i2=a2f,j2=0;i2<a2l;i2++,j2++)
      for (i1=a1f,j1=0;i1<a1l;i1++,j1++)
	m[j2][j1]=a[i2][i1];
  }
  else{
    for (i2=a2f,j2=0;i2<a2l;i2++,j2++)
      for (i1=a1f,j1=0;i1<a1l;i1++,j1++)
	a[i2][i1]=m[j2][j1];
  } 
  
  return;

}


float **eallocwindow(int a2f, int a2l, int a1f, int a1l)
{
  float **a;
  a=ealloc2float(a1l-a1f+1,a2l-a2f+1);
  return(a);

}

