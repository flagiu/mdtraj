#ifndef _MYMATRIX_H_
#define _MYMATRIX_H_
#include "myvec.hpp"

template <class ntype, int Nr, int Nc>
class mymatrix
{
  using mv=myvec<ntype,Nc>;
  mv M[Nr]; //row vectors
  unsigned long int idx=0;

public:

  inline mymatrix<ntype,Nr,Nc>& operator,(const ntype& val)
    {
      idx++;
      if (idx >= Nr*Nc)
        {
          cout << "Too many elements in comma initialization!\n";
          exit(-1);
        }
      M[(int)idx/Nc][idx%Nc]=val;
      return (*this);
    }
  inline mymatrix<ntype,Nr,Nc>& operator<<(const ntype& val)
    {
      M[0][0] = val;
      idx=0;
      return (*this);
    } 
  inline mv& operator[](const int& i)
  {
    return M[i];
  }
  inline const mv& operator[](const int& i) const
  {
    return M[i];
  }
  void set(const int i, const int j, const ntype& val)
  {
    M[i].set(j,val);
  }
  void set(int i, mv& vv)
  {
    for (int j=0;j<cols();j++)
      set(i,j,vv[j]);
  }
  int rows() const
  {
    return Nr;
  }
  int cols() const
  {
    return Nc;
  }
  void show() const
  {
    cout<<"[\n";
    for (int i=0;i<rows();i++)
	M[i].show();
    cout<<"]\n";
  }
  void show(const char *str) const
  {
    cout<<str<<": ";
    show();
  }
  template <class ntype2>
  void to_array(ntype2 arr[Nr][Nc]) const
  {
    for (int i=0;i<rows();i++)
	M[i].to_array(arr[i]);
  }
  // Methods returning *this:
  inline mymatrix<ntype,Nr,Nc>& operator=(const mymatrix<ntype,Nr,Nc>& other)
  {
    for (int i=0;i<rows();i++)
      M[i]=other[i]; //cast to non-const
    return (*this);
  }
  inline mymatrix<ntype,Nr,Nc>& operator += (const mymatrix<ntype,Nr,Nc>& other)
  {//matrix addition
    int i;
    for (i=0; i < rows(); i++)
      M[i] += other[i];
    return *this;
  }
  inline mymatrix<ntype,Nr,Nc>& operator += (const ntype& sc)
  {// scalar addition
    int i;
    for (i=0; i < rows(); i++)
      M[i] += sc;
    return *this;
  }
  inline mymatrix<ntype,Nr,Nc>& operator *= (const ntype& sc)
  {// scalar multiplication
    int i;
    for (i=0; i < rows(); i++)
      M[i] *= sc;
    return *this;
  }
  inline mymatrix<ntype,Nr,Nc>& operator -= (const ntype& sc)
  {// scalar subtraction
    int i;
    for (i=0; i < rows(); i++)
      M[i] -= sc;
    return *this;
  }
  inline mymatrix<ntype,Nr,Nc>& operator /= (const ntype& sc)
  {// scalar division
    int i;
    for (i=0; i < rows(); i++)
      M[i] /= sc;
    return *this;
  }

  mymatrix<ntype,Nc,Nr> T()
  {
    mymatrix<ntype,Nc,Nr> mat;
    for (auto i=0;i<rows();i++)
      for (auto j=0;j<cols();j++)
	mat[j][i] = M[i][j];
    return mat;
  }

  const mymatrix<ntype,Nr-1,Nc-1> minore(const int& a, const int& b) const
  {
    mymatrix<ntype,Nr-1,Nc-1> mat;
    int i=0,j=0;
    while(i<rows()&&j<rows()-1)
      {
	if (i==a) i++;
	mat[j] = M[i].pop(b);
	i++; j++;
      }
    return mat;
  }

  ntype det() const
  {
    ntype d=0.0;
    if (rows()!=cols())
      {
	cout << "Error: determinant is defined only for square matrices. Returning zero."<<endl;
      }
    else if (Nr==1) M[0][0];
    else
      {
	int sign=1;
	mymatrix<ntype,Nr-1,Nc-1> m;
	for (auto i=0;i<rows();i++)
	  {
	    m=minore(i,0);
	    d += sign*M[i][0]*m.det();
	    sign=-sign;
	  }
      }
    return d;
  }
    
  //
  void ranf()
  {
    for (auto i=0;i<rows();i++)
      M[i].ranf();
  }
};
template<class ntype,int Nr,int Nc>
myvec<ntype,Nr> operator*(mymatrix<ntype,Nr,Nc>& mat, myvec<ntype,Nc>& vec)
{//vector multiplication from right
  myvec<ntype,Nr> vv;
  for (auto i=0;i<Nr;i++)
    vv[i] = mat[i]*vec;
  return vv;
}
template<class ntype,int Nr,int Nc>
myvec<ntype,Nc> operator*( myvec<ntype,Nr>& vec, mymatrix<ntype,Nr,Nc>& mat)
{//vector multiplication from left
  mymatrix<ntype,Nc,Nr> matT=mat.T();
  return matT*vec;
}
template<class ntype,int Nr,int Nc>
mymatrix<ntype,Nr,Nr> operator*(mymatrix<ntype,Nr,Nc>& m1, mymatrix<ntype,Nc,Nr>& m2)
{//matrix multiplication
  mymatrix<ntype,Nr,Nr> mat;
  for (auto i=0;i<Nr;i++)
    mat[i] = m1[i]*m2;
  return mat;
}
using matrix3d=mymatrix<double,3,3>;

template <class ntype> class mymatrix<ntype,1,1>
{
  using mv=myvec<ntype,1>;
  mv M[1];
public:
  inline mv& operator[](const int& i)
  {
    return M[i];
  }
  const inline mv& operator[](const int& i) const
  {
    return M[i];
  }
  void show()
  {
    cout<<"["<<M[0][0]<<"]"<<endl;
  }
  void show(const char *str)
  {
    cout<<str<<": ";
    show();
  }
  ntype det()
  {
    return M[0][0];
  }
  mymatrix()
  {
  }
  ~mymatrix()
  {}
};

#endif
