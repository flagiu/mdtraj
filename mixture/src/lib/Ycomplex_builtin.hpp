/*
	Spherical harmonics
https://en.wikipedia.org/wiki/Table_of_spherical_harmonics

https://en.cppreference.com/w/cpp/numeric/special_functions/sph_legendre
*/

#ifndef _YCOMPLEX_BUILTIN_H_
#define _YCOMPLEX_BUILTIN_H_
#include<cstdlib>
#include<complex>
#include<cmath>
#include "utility.hpp"
using namespace std;

template<class ntype>
complex<ntype> Y(int l, int m, myvec<ntype,3> R)
{
  return Y(l,m,R,R.norm());
}

template<class ntype>
complex<ntype> Y(int l, int m, myvec<ntype,3> R, ntype r)
{
  ntype pi=M_PI, x=R[0], y=R[1], z=R[2];
  unsigned int ul=l; // degree of the spherical harmonic
  unsigned int um=(m>0?m:-m); // order of the spherical harmonic
  if(!(l>=0)) { cerr<<"[ERROR: negative angular momentum selected in Y(): l = " << l << "]\n"; exit(1); }
  if(!(l<=127)) { cerr<<"[ERROR: l>127 gives unpredictable Y(): l = " << l << "]\n"; exit(1); }
  if(!(ul>=um)) { cerr<<"[ERROR: condition l>=|lz| not satisfied in Y(): l = " << l << ", lz = " << m << "]\n"; exit(1); }
  if(r==0.0) { cerr<<"[ERROR: cannot compute Y() of a null direction]\n"; exit(1); }
  ntype theta=acos(z/r); // polar angle
  ntype phi=atan2(y,x); // azimuthal angle
  return sph_legendre(ul,um,theta) * exp(1i*(m*phi));
}
# endif
