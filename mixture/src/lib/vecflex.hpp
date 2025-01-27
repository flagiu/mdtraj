#ifndef _VECFLEX_H_
#define _VECFLEX_H_

#include "myvec.hpp"
using namespace std;
template <class ntype>
class vecflex
{
  ntype* v;
  int N;
  unsigned long int idx=0;
public:
  //---- Constructors ----//
  vecflex(){ N=0; }
  ~vecflex() {} // QUESTO DAREBBE ERRORE munmap_chunk(): invalid pointer --> { delete[] v; }
  //------- Operations ----------//
  void resize(int NN)
    {
      N = NN;
      v = new ntype[N];
    }
  inline vecflex<ntype>& operator,(const ntype& val)
    {
      idx++;
      if (idx >= N)
        {
          cerr << "Too many elements in comma initialization!\n";
          exit(-1);
        }
      v[idx]=val;
      return (*this);
    }
  inline vecflex<ntype>& operator<<(const ntype& val)
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

  void show(ostream& o)
  {
    o << "(";
    for (int i = 0; i < length(); i++) {
      o << v[i];
      if (i < length()-1)
        o << ", ";
      else
        o << ")\n";
    }
  }
  void show(ostream& o, const char *myname) const
  {
    o << myname << ": ";
    show(o);
  }

  void show(void) const
    {
      show(cerr);
    }
  void show(const char *myname) const
  {
    show(cerr, myname);
  }

  template<class ntype2>
  void to_array(ntype2* arr, int N2) const
  {
    if(length()!=N2) { cerr << "[ERROR: trying to convert a vecflex to a different-sized array]\n"; exit(1); }
    for (int i=0;i<length();i++)
	arr[i] = (ntype2)(v[i]);
  }

  // overloaded operations returning THIS vecflex
  //// operations with vectors - returning *this

  inline vecflex<ntype>& operator=(const vecflex<ntype>& other) //assignment
    {
      if(length()!=other.length()) { cerr << "[ERROR: different-sized vecflex]\n"; exit(1); }
      for (int i=0; i < length(); i++)
        v[i] = other[i]; //cast to non-const
      return (*this);
    }
  inline vecflex<ntype>& operator += (const vecflex<ntype>& other)
  {
    if(length()!=other.length()) { cerr << "[ERROR: different-sized vecflex]\n"; exit(1); }
    int i;
    for (i=0; i < length(); i++)
      v[i] += other[i];
    return *this;
  }
  inline vecflex<ntype>& operator -= (const vecflex<ntype>& other)
  {
    if(length()!=other.length()) { cerr << "[ERROR: different-sized vecflex]\n"; exit(1); }
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
  inline vecflex<ntype>& operator *= (const ntype& sc)
  {
    for (int i=0; i < length(); i++)
      v[i] *= sc;
    return *this;
  }
  inline vecflex<ntype>& operator *= (const int& sc)
  {
    for (int i=0; i < length(); i++)
      v[i] *= sc;
    return *this;
  }
  inline vecflex<ntype>& operator += (const ntype& sc)
  {
    for (int i=0; i < length(); i++)
      v[i] += sc;
    return *this;
  }
  inline vecflex<ntype>& operator /= (const ntype& sc)
  {
    for (int i=0; i < length(); i++)
      v[i] /= sc;
    return *this;
  }
  inline vecflex<ntype>& operator -= (const ntype& sc)
  {
    for (int i=0; i < length(); i++)
      v[i] -= sc;
    return *this;
  }
  inline vecflex<ntype>& operator ^= (const ntype& sc)
  {
    for (int i=0; i < length(); i++)
      v[i] = pow(v[i],sc);
    return *this;
  }
  inline vecflex<ntype>& exp_from_base(const ntype& sc)
  {
    for (int i=0; i < length(); i++)
      v[i] = pow(sc,v[i]);
    return *this;
  }

  inline ntype operator*(vecflex<ntype>& other) //dot product
  {
    if(length()!=other.length()) { cerr << "[ERROR: different-sized vecflex]\n"; exit(1); }
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

  inline vecflex<ntype> project_on(vecflex<ntype>& other)
  {
    if(length()!=other.length()) { cerr << "[ERROR: different-sized vecflex]\n"; exit(1); }
    vecflex<ntype> versor;
    versor.resize( length() );
    ntype proj;
    versor=other/other.norm();
    proj=(*this)*versor;
    return versor*proj;
  }
  inline vecflex<ntype> cross(const vecflex<ntype>& other)
  {
     if(length()!=other.length()) { cerr << "[ERROR: different-sized vecflex]\n"; exit(1); }
    vecflex<ntype> vv;
    vv.resize( length() );
    if (length()==3 && other.length()==3)
      {
	vv[0] = v[1]*other[2]-v[2]*other[1];
	vv[1] = v[2]*other[0]-v[0]*other[2];
	vv[2] = v[0]*other[1]-v[1]*other[0];
      }
    else
      cerr << "Error: vector product is valid only in 3D!\nReturning null vector.\n";
    return vv;
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
  vecflex<ntype> mic(vecflex<ntype> L)
  {
    vecflex<ntype> vnew;
    vnew.resize(length());
    for (auto i=0;i<length();i++)
      vnew[i] = v[i] - L[i]*rint(v[i]/L[i]);
    return vnew;
  }
  vecflex<ntype> mic(ntype L)
  {
    vecflex<ntype> vnew;
    vnew.resize(length());
    for (auto i=0;i<length();i++)
      vnew[i] = v[i] - L*rint(v[i]/L);
    return vnew;
  }

  ntype mean()
  {
    ntype s=0.0;
    for (int i=0;i<length();i++)
      s += v[i];
    return s/(ntype)(length());
  }

  ntype var()
  {
    if (length()==1) return 0.0;
    ntype s=0.0;
    ntype m=mean();
    for (int i=0;i<length();i++)
      s += (v[i]-m)*(v[i]-m);
    return s/(ntype)(length()-1);
  }

  ntype std()
  {
    return sqrt(var());
  }

  ntype prod()
  {
    ntype p = 1.0;
    for (int i=0;i<length();i++)
      p *= v[i];
    return p;
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

  vecflex<ntype> pop(const int& a) const
  {
    vecflex<ntype> vec;
    vec.resize( length()-1 );
    int i=0,j=0;
    while (i<length())
      {
	if (i==a) i++;
	vec[j]= v[i];
	i++; j++;
      }
    return vec;
  }
};

template<class ntype>
inline vecflex<ntype>& operator+(vecflex<ntype> v1, vecflex<ntype> v2)
{ // binary vector sum
  if(v1.length()!=v2.length()) { cerr << "[ERROR: different-sized vecflex]\n"; exit(1); }
  return v1+=v2;
}
template<class ntype>
inline vecflex<ntype>& operator-(vecflex<ntype> v1, vecflex<ntype> v2)
{ // binary vector sum
  if(v1.length()!=v2.length()) { cerr << "[ERROR: different-sized vecflex]\n"; exit(1); }
  return v1-=v2;
}

template<class ntype>
inline vecflex<ntype>& operator*(vecflex<ntype> v1, const ntype& sc)
{ // binary scalar multiplication from right
  return v1*=sc;
}
template<class ntype>
inline vecflex<ntype>& operator*(const ntype& sc,  vecflex<ntype> v1)
{ // binary scalar multiplication from left
  return v1*=sc;
}
template<class ntype>
inline vecflex<ntype>& operator+(vecflex<ntype> v1, const ntype& sc)
{ // binary scalar sum from right
  return v1+=sc;
}
template<class ntype>
inline vecflex<ntype>& operator+(const ntype& sc,  vecflex<ntype> v1)
{ // binary scalar sum from left
  return v1+=sc;
}
template<class ntype>
inline vecflex<ntype>& operator/(vecflex<ntype> v1, const ntype& sc)
{ // binary scalar division from right
  return v1/=sc;
}
template<class ntype>
inline vecflex<ntype>& operator/(const ntype& sc,  vecflex<ntype> v1)
{ // binary scalar division from left
  return v1/=sc;
}
template<class ntype>
inline vecflex<ntype>& operator-(vecflex<ntype> v1, const ntype& sc)
{ // binary scalar subtraction from right
  return v1-=sc;
}
template<class ntype>
inline vecflex<ntype>& operator-(const ntype& sc,  vecflex<ntype> v1)
{ // binary scalar subtraction from left
  return v1-=sc;
}

template<class ntype>
inline vecflex<ntype>& operator^(vecflex<ntype> v1, const ntype& sc)
{ // binary scalar power (only from right)
  return v1^=sc;
}
template<class ntype>
inline vecflex<ntype>& operator^(const ntype& sc,  vecflex<ntype> v1)
{ // power the base on the left to the vector
  return v1.exp_from_base(sc);
}

template<class ntype>
inline vecflex<ntype>& operator*(vecflex<ntype> v1, const int& sc)
{ // power the base on the left to the vector
  return v1*=sc;
}
template<class ntype>
inline vecflex<ntype>& operator*(const int& sc,  vecflex<ntype> v1)
{ // power the base on the left to the vector
  return v1*=sc;
}
template<class ntype>
vecflex<ntype> rint(vecflex<ntype> v1)
{ // round to integer the whole vector
  vecflex<ntype> u;
  u.resize( v1.length() );
  for(auto i=0;i<v1.length();i++)
    u[i] = (ntype)rint(v1[i]);
  return u;
}
#endif
