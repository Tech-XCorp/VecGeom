/*
 * Vector3D.h
 *
 *  Created on: Nov 10, 2013
 *      Author: swenzel
 */

#ifndef VECTOR3D_H_
#define VECTOR3D_H_

#include "GlobalDefs.h"
#include "mm_malloc.h"
#include <vector>

struct Vector3D
{
	double x,y,z;
	Vector3D(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {};
	Vector3D() : x(0), y(0), z(0) {};

	/*
	double operator[](int i) const
	{
		return &x+i;
	}
*/

	Vector3D & operator=(Vector3D const & rhs)
	{
		this->x=rhs.x;
		this->y=rhs.y;
		this->z=rhs.z;
		return *this;
	}
};

struct Vectors3DSOA
{
private:
  double * xvec;
  double * yvec;
  double * zvec;

public:
  int size; // the size
  double *x; // these are just "views" to the real data in xvec ( with the idea that we can set x to point to different start locations in xvec )
  double *y;
  double *z;

  Vector3D getAsVector(int index) const
  {
	  return Vector3D(x[index],y[index],z[index]);
  }

  void getAsVector(int index, Vector3D & v) const
  {
	  v.x=x[index];
  	  v.y=y[index];
  	  v.z=z[index];
  }

  void toStructureOfVector3D( std::vector<Vector3D> & v )
  {
	  for(auto i=0;i<size;i++)
	  {
		  getAsVector(i, v[i]);
	  }
  }

  void setFromVector(int index, Vector3D const &v) const
  {
  	  x[index]=v.x;
  	  y[index]=v.y;
  	  z[index]=v.z;
  }

  void setstartindex(int index)
  {
    x=&xvec[index];
    y=&yvec[index];
    z=&zvec[index];
  }

  // checks if index will be aligned and if not go back to an aligned address
  void setstartindex_aligned(int /*index*/)
  {
    // TOBEDONE
  }

  void fill(double const * onedvec)
  {
    for(int i=0;i<size;++i)
      {
    	xvec[i]=onedvec[3*i];
    	yvec[i]=onedvec[3*i+1];
    	zvec[i]=onedvec[3*i+2];
      }
  }

  void alloc(int n)
  {
    this->size=n;
    xvec=(double*)_mm_malloc(sizeof(double)*size,ALIGNMENT_BOUNDARY); // aligned malloc (32 for AVX )
    yvec=(double*)_mm_malloc(sizeof(double)*size,ALIGNMENT_BOUNDARY);
    zvec=(double*)_mm_malloc(sizeof(double)*size,ALIGNMENT_BOUNDARY);
    x=xvec;y=yvec;z=zvec;
  }

  void dealloc()
  {
    _mm_free(xvec);
    _mm_free(yvec);
    _mm_free(zvec);
    x=0;y=0;z=0;
  }
};

 #endif /* VECTOR3D_H_ */
