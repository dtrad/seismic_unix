#ifndef ALLOC_HPP
#define ALLOC_HPP


#include <iostream>



// Template functions to allocate/delete 1d/2d/3d buffer
//

template<class A>
void new1d( A* &p1d, unsigned int size )
{
     p1d = new A[size];
}

template<class A>
void del1d( A* &p1d )
{
    delete []p1d;
    p1d = 0;
}


template<class A>
void new2d( A** &p2d, unsigned int size1,  unsigned int size2 )
{
     p2d = 0; p2d = new A*[size1];
     p2d[0] = new A[size1*size2];
     for( unsigned int i = 1; i < size1; ++i )
	  p2d[i] = p2d[i-1] + size2;

}

template<class A>
void del2d( A** &p2d )
{
    if( p2d != 0 ) delete []p2d[0];
    delete []p2d;
    p2d = 0;
}

template<class A> 
void new3d( A*** &p3d, unsigned int size1, unsigned int size2, unsigned int size3 )
{
     p3d = 0; p3d = new A**[size1];
     p3d[0] = 0; p3d[0] = new A*[size1*size2];
     p3d[0][0] = new A[size1*size2*size3];
     for( unsigned int i = 0, offset = 0; i < size1; ++i ) {
	  p3d[i] = p3d[0] + size2*i;
	  for( unsigned int j = 0; j < size2; ++j, offset += size3 ) {
	       p3d[i][j] = p3d[0][0] + offset;
	  }
     }
     
}

template<class A> 
void del3d( A*** &p3d )
{
    if( p3d !=0 ) {
	 delete []p3d[0][0];
	 delete []p3d[0];        
	 delete []p3d;
    }
    p3d = 0;
}

// delete a singular object, assign NULL to pointer
//
template< class A >
inline void MDEL( A& aa )
{
     if( aa != NULL ) {
	  delete aa;
	  aa = NULL;
     }
}

// delete an array pointer, assign NULL to pointer
//
template< class A >
inline void MDELARR( A& aa )
{
     if( aa != NULL ) {
	  delete [] aa;
	  aa = NULL;
     }
}




#endif
