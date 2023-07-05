#ifndef _MYVEC_H_
#define _MYVEC_H_

#include<cstdlib>
#include<iostream>
#include<cmath>
#include "randnumgen.hpp"
using namespace std;
template <class ntype, int N>
class myvec
{
  ntype v[N];
  unsigned long int idx=0;
public:
  inline myvec<ntype,N>& operator,(const ntype& val)
    {
      idx++;
      if (idx >= N)
        {
          cout << "Too many elements in comma initialization!\n";
          exit(-1);
        }
      v[idx]=val;
      return (*this);
    }
  inline myvec<ntype,N>& operator<<(const ntype& val)
    {
      v[0] = val;
      idx=0;
      return (*this);
    } 
  inline ntype& operator[](const int& i)
  {
    return v[i];
  }
  inline const ntype& operator[](const int& i) const
  {
    return v[i];
  }
  int length(void) const
    {
      return N;
    }
  void set(int i, ntype val)
    {
      v[i] = val;
    }
  void show(void) const
    {
      cout << "(";
      for (int i = 0; i < length(); i++) {
        cout << v[i];
        if (i < length()-1)
          cout << ",";
        else
          cout << ")\n";
      }
    }
  void show(const char *str) const
  {
    cout << str << ": ";
    show();
  }

  template <class ntype2>
  void to_array(ntype2 arr[N]) const
  {
    for (int i=0;i<length();i++)
	arr[i] = (ntype2)(v[i]);
  }

  // overloaded operations returning THIS myvec
  //// operations with vectors - returning *this
  
  inline myvec<ntype,N>& operator=(const myvec<ntype,N>& other) //assignment
    {
      for (int i=0; i < length(); i++)
        v[i] = other[i]; //cast to non-const
      return (*this);
    }
  inline myvec<ntype,N>& operator += (const myvec<ntype,N>& other)
  {
    int i;
    for (i=0; i < length(); i++)
      v[i] += other[i];
    return *this;
  }
  inline myvec<ntype,N>& operator -= (const myvec<ntype,N>& other)
  {
    int i;
    for (i=0; i < length(); i++)
      v[i] -= other[i];
    return *this;
  }

  void normalize()
  {
    ntype n=norm();
    if (n!=0)
      (*this)/=n;
  }
  
  ////operations with scalars - returning *this
  inline myvec<ntype,N>& operator *= (const ntype& sc)
  {
    for (int i=0; i < length(); i++)
      v[i] *= sc;
    return *this;
  }
  inline myvec<ntype,N>& operator *= (const int& sc)
  {
    for (int i=0; i < length(); i++)
      v[i] *= sc;
    return *this;
  }
  inline myvec<ntype,N>& operator += (const ntype& sc)
  {
    for (int i=0; i < length(); i++)
      v[i] += sc;
    return *this;
  }
  inline myvec<ntype,N>& operator /= (const ntype& sc)
  {
    for (int i=0; i < length(); i++)
      v[i] /= sc;
    return *this;
  }
  inline myvec<ntype,N>& operator -= (const ntype& sc)
  {
    for (int i=0; i < length(); i++)
      v[i] -= sc;
    return *this;
  }
  inline myvec<ntype,N>& operator ^= (const ntype& sc)
  {
    for (int i=0; i < length(); i++)
      v[i] = pow(v[i],sc);
    return *this;
  }
  inline myvec<ntype,N>& exp_from_base(const ntype& sc)
  {
    for (int i=0; i < length(); i++)
      v[i] = pow(sc,v[i]);
    return *this;
  }

  inline ntype operator*(myvec<ntype,N>& other) //dot product
  {
    ntype sum=0;
    for (int i=0;i<length();i++)
	sum += v[i]*other[i];
    return sum;
  }

  inline ntype sq(void) //squared norm
  {
    return (*this)*(*this);
  }

  inline ntype norm(void)
  {
    return sqrt(sq());
  }

  inline myvec<ntype,N> project_on(myvec<ntype,N>& other)
  {
    myvec<ntype,N> versor;
    ntype proj;
    versor=other/other.norm();
    proj=(*this)*versor;
    return versor*proj;
  }
  inline myvec<ntype,N> cross(const myvec<ntype,N>& other)
  {
    myvec<ntype,N> vv;
    if (length()==3 && other.length()==3)
      {
	vv[0] = v[1]*other[2]-v[2]*other[1];
	vv[1] = v[2]*other[0]-v[0]*other[2];
	vv[2] = v[0]*other[1]-v[1]*other[0];
      }
    else
      cout << "[Error: vector product is valid only in 3D! Returning null vector.]\n";
    return vv;
  }
  
  void set_zero()
  {
    for (int i=0;i<length();i++)
      set(i, 0);
  }

  void ranf()
  {
    for (int i=0;i<length();i++)
      set(i,randnum.ranf());
  }

  void rand_box()
  {
    for (int i=0;i<length();i++)
      set(i,randnum.ranf()-0.5);
  }
  myvec<ntype,N> mic(myvec<ntype,N> L) const
  {
    myvec<ntype,N> vnew;
    for (auto i=0;i<length();i++)
      vnew[i] = v[i] - L[i]*rint(v[i]/L[i]);
    return vnew;
  }
  myvec<ntype,N> mic(ntype L) const
  {
    myvec<ntype,N> vnew;
    for (auto i=0;i<length();i++)
      vnew[i] = v[i] - L*rint(v[i]/L);
    return vnew;
  }

  ntype mean() const
  {
    ntype s=0.0;
    for (int i=0;i<length();i++)
      s += v[i];
    return s/(ntype)(length());
  }

  ntype var() const
  {
    if (length()==1) return 0.0;
    ntype s=0.0;
    ntype m=mean();
    for (int i=0;i<length();i++)
      s += (v[i]-m)*(v[i]-m);
    return s/(ntype)(length()-1);
  }

  ntype std() const
  {
    return sqrt(var());
  }

  ntype prod() const
  {
    ntype p = 1.0;
    for (int i=0;i<length();i++)
      p *= v[i];
    return p;
  }

  myvec<ntype,N-1> pop(const int& a) const
  {
    myvec<ntype,N-1> vec;
    int i=0,j=0;
    while (i<length())
      {
	if (i==a) i++;
	vec[j]= v[i];
	i++; j++;
      }
    return vec;
  }
  
  myvec<ntype,N> abs(void) const
  {
    myvec<ntype,N> vec;
    for(int i=0;i<length();i++)
	vec[i]= fabs(v[i]);
    return vec;
  }
  
  ntype min(void) const
  {
    ntype m=v[0];
    for(int i=1;i<length();i++)
	if( v[i]<m )
	  m = v[i];
    return m;
  }
  
  ntype max(void) const
  {
    ntype m=v[0];
    for(int i=1;i<length();i++)
	if( v[i]>m )
	  m = v[i];
    return m;
  }
  
  myvec()
    {}
  ~myvec()
    {}
};

template<class ntype, int N>
inline myvec<ntype,N>& operator+(myvec<ntype,N> v1, myvec<ntype,N> v2)
{ // binary vector sum
  return v1+=v2;
}
template<class ntype, int N>
inline myvec<ntype,N>& operator-(myvec<ntype,N> v1, myvec<ntype,N> v2)
{ // binary vector sum
  return v1-=v2;
}

template<class ntype, int N>
inline myvec<ntype,N>& operator*(myvec<ntype,N> v1, const ntype& sc)
{ // binary scalar multiplication from right
  return v1*=sc;
}
template<class ntype, int N>
inline myvec<ntype,N>& operator*(const ntype& sc,  myvec<ntype,N> v1)
{ // binary scalar multiplication from left
  return v1*=sc;
}
template<class ntype, int N>
inline myvec<ntype,N>& operator+(myvec<ntype,N> v1, const ntype& sc)
{ // binary scalar sum from right
  return v1+=sc;
}
template<class ntype, int N>
inline myvec<ntype,N>& operator+(const ntype& sc,  myvec<ntype,N> v1)
{ // binary scalar sum from left
  return v1+=sc;
}
template<class ntype, int N>
inline myvec<ntype,N>& operator/(myvec<ntype,N> v1, const ntype& sc)
{ // binary scalar division from right
  return v1/=sc;
}
template<class ntype, int N>
inline myvec<ntype,N>& operator/(const ntype& sc,  myvec<ntype,N> v1)
{ // binary scalar division from left
  return v1/=sc;
}
template<class ntype, int N>
inline myvec<ntype,N>& operator-(myvec<ntype,N> v1, const ntype& sc)
{ // binary scalar subtraction from right
  return v1-=sc;
}
template<class ntype, int N>
inline myvec<ntype,N>& operator-(const ntype& sc,  myvec<ntype,N> v1)
{ // binary scalar subtraction from left
  return v1-=sc;
}

template<class ntype, int N>
inline myvec<ntype,N>& operator^(myvec<ntype,N> v1, const ntype& sc)
{ // binary scalar power (only from right)
  return v1^=sc;
}
template<class ntype, int N>
inline myvec<ntype,N>& operator^(const ntype& sc,  myvec<ntype,N> v1)
{ // power the base on the left to the vector
  return v1.exp_from_base(sc);
}

template<class ntype, int N>
inline myvec<ntype,N>& operator*(myvec<ntype,N> v1, const int& sc)
{ // scalar multiplication from right
  return v1*=sc;
}
template<class ntype, int N>
inline myvec<ntype,N>& operator*(const int& sc,  myvec<ntype,N> v1)
{ // scalar multiplication from left
  return v1*=sc;
}
template<class ntype, int N>
myvec<ntype,N> round(myvec<ntype,N> v1)
{ // round the whole vector
  myvec<ntype,N> u;
  for(auto i=0;i<N;i++)
    u[i] = (ntype)round(v1[i]);
  return u;
}
//using vec2d=myvec<double,2>;
//using vec3d=myvec<double,3>;
#endif
