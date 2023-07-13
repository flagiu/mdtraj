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
  void write(fstream& o) const
  {
    for (auto i=0;i<rows();i++)
    	for (auto j=0;j<cols();j++)
    		{
    		  o << setprecision(20) << M[i][j];
    		  if(i==rows()-1 && j==cols()-1) o << '\n';
    		  else				 o << ' ';
    		}
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
  inline mymatrix<ntype,Nr,Nc>& operator -= (const mymatrix<ntype,Nr,Nc>& other)
  {//matrix diff
    int i;
    for (i=0; i < rows(); i++)
      M[i] -= other[i];
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
  
  myvec<ntype,Nr> operator*(myvec<ntype,Nc> vec)
  {//vector multiplication from right
    myvec<ntype,Nr> vv;
    for (auto i=0;i<Nr;i++) {
      vv[i] = M[i]*vec;
    }
    return vv;
  }
  mymatrix<ntype,Nr,Nr> operator*(const mymatrix<ntype,Nc,Nr>& other)
  {//matrix multiplication
    mymatrix<ntype,Nr,Nr> mat;
    for (auto i=0;i<Nr;i++)
      mat[i] = M[i]*other;
    return mat;
  }

  //----------------------------------
  mymatrix<ntype,Nc,Nr> T() const
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
  
  mymatrix<ntype,Nr,Nr> square() const
  {
    
    mymatrix<ntype,Nc,Nr> mat = (*this);
    return mat*(this->T());
  }
  
  ntype tr() const
  {
    ntype t=0.0;
    if (rows()!=cols()) { cout << "Error: trace is defined only for square matrices."<<endl; exit(1); }
    for(auto i=0;i<rows();i++) t+=M[i][i];
    return t;
  }

  ntype det() const
  {
    ntype d=0.0;
    if (rows()!=cols()) { cout << "Error: determinant is defined only for square matrices."<<endl; exit(1);}
    if (Nr==1) d=M[0][0];
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
  void set_zero()
  {
    for (auto i=0;i<rows();i++)
      for (auto j=0;j<rows();j++)
        set(i,j, 0.0);
    return;
  }
  void set_identity()
  {
    if (rows()!=cols()) { cout << "Error: identity is defined only for square matrices."<<endl; exit(1);}
    for (auto i=0;i<rows();i++) {
      for (auto j=0;j<rows();j++) {
        if(i==j) set(i,j, 1.0);
        else     set(i,j, 0.0);
      }
    }
    return;
  }
  
  void set_diag(myvec<ntype,Nr> vec)
  {
    if (rows()!=cols()) { cout << "Error: diagonal is defined only for square matrices."<<endl; exit(1);}
    for (auto i=0;i<rows();i++) {
      for (auto j=0;j<rows();j++) {
        if(i==j) set(i,j, vec[i]);
        else     set(i,j,    0.0);
      }
    }
    return;
  }

  myvec<ntype,Nr> diag() const
  {
    myvec<ntype,Nr> vec;
    if (rows()!=cols()) { cout << "Error: diagonal is defined only for square matrices."<<endl; exit(1);}
    for (auto i=0;i<rows();i++) vec[i] = M[i][i];
    return vec;
  }
  
  mymatrix<ntype,Nr,Nc> inverse() const
  {
    mymatrix<ntype,Nr,Nc> m;
    ntype d;
    if (rows()==2 && cols()==2) {
       d = M[0][0]*M[1][1] - M[0][1]*M[1][0];
       if( d==0.0) { cout << "[ERROR: this matrix is not invertible.]"<<endl; show(); exit(1); }
       m[0][0] =  M[1][1]/d;
       m[0][1] = -M[0][1]/d;
       m[1][0] = -M[1][0]/d;
       m[1][1] =  M[0][0]/d;
    }
    else if (rows()==3 && cols()==3) {
       d = det();
       if( d==0.0) { cout << "[ERROR: this matrix is not invertible.]"<<endl; show(); exit(1); }
       ntype trace, traceSq;
       mymatrix<ntype,Nr,Nc> self, sq, eye;
       self = (*this);
       sq = this->square();
       eye.set_identity();
       trace = tr();
       traceSq = sq.tr();
       m = ( ( 0.5*(trace*trace - traceSq)*eye ) - ( trace*self ) + sq ) / d;
    }
    else { cout << "[ERROR: inverse() not yet implemented for " << rows() << "x" << cols() << " matrices.]"<<endl; exit(1); }
    return m;
  }
};

template<class ntype,int Nr,int Nc>
inline mymatrix<ntype,Nr,Nc>& operator+(mymatrix<ntype,Nr,Nc>& m1, mymatrix<ntype,Nr,Nc>& m2)
{//pointwise sum
  return m1+=m2;
}
template<class ntype,int Nr,int Nc>
inline mymatrix<ntype,Nr,Nc>& operator-(mymatrix<ntype,Nr,Nc>& m1, mymatrix<ntype,Nr,Nc>& m2)
{//pointwise diff
  return m1-=m2;
}

template<class ntype,int Nr,int Nc>
inline mymatrix<ntype,Nr,Nc>& operator*(const ntype sc, mymatrix<ntype,Nr,Nc>& m)
{//scalar mult from left
  return m*=sc;
}

template<class ntype,int Nr,int Nc>
inline mymatrix<ntype,Nr,Nc>& operator*(mymatrix<ntype,Nr,Nc>& m, const ntype sc)
{//scalar mult from right
  return m*=sc;
}

template<class ntype,int Nr,int Nc>
inline mymatrix<ntype,Nr,Nc>& operator/(mymatrix<ntype,Nr,Nc>& m, const ntype sc)
{//scalar div from right
  return m/=sc;
}
/*
template<class ntype,int Nr,int Nc>
myvec<ntype,Nr> operator*(const mymatrix<ntype,Nr,Nc>& mat, const myvec<ntype,Nc>& vec)
{//vector multiplication from right
  myvec<ntype,Nr> vv;
  for (auto i=0;i<Nr;i++)
    vv[i] = mat[i]*vec;
  return vv;
}

template<class ntype,int Nr,int Nc>
mymatrix<ntype,Nr,Nr>& operator*(const mymatrix<ntype,Nr,Nc>& m1, const mymatrix<ntype,Nc,Nr>& m2)
{//matrix multiplication
  mymatrix<ntype,Nr,Nr> mat;
  for (auto i=0;i<Nr;i++)
    mat[i] = m1[i]*m2;
  return mat;
}
*/
template<class ntype,int Nr,int Nc>
myvec<ntype,Nc> operator*(myvec<ntype,Nr>& vec, const mymatrix<ntype,Nr,Nc>& mat)
{//vector multiplication from left
  mymatrix<ntype,Nc,Nr> matT=mat.T();
  return matT*vec;
}
//vector multiplication from right is implemented inside the class!

template<class ntype, int N>
myvec<ntype,N> mic(mymatrix<ntype,N,N> box, mymatrix<ntype,N,N>& boxInv, myvec<ntype,N>& vec)
{ // apply MIC for non-orthorombic boxes (if you already know the inverse)
  return vec - box * round(boxInv*vec);
}
template<class ntype, int N>
myvec<ntype,N> mic(mymatrix<ntype,N,N>& box, myvec<ntype,N>& vec)
{ // apply MIC for non-orthorombic boxes
  return mic(box, box.inverse(), vec);
}

template<class ntype, int N>
mymatrix<ntype,N,N> outer_product( myvec<ntype,N> vec1, myvec<ntype,N> vec2 )
{
  mymatrix<ntype,N,N> mat;
  for(auto i=0;i<N;i++)
    for(auto j=0;j<N;j++)
      mat[i][j] = vec1[i]*vec2[j];
  return mat;
}
template<class ntype, int N>
mymatrix<ntype,N,N> cross_product_matrix(myvec<ntype,N> vec)
{
  mymatrix<ntype,N,N> mat;
  myvec<ntype,N> e;
  mat.set_zero();
  e.set_zero();
  for(auto i=0;i<N;i++)
  {
      e[i] = 1.0;
      mat += outer_product(vec.cross(e), e);
      e[i] = 0.0;
  }
  return mat;
}

template<class ntype, int N>
mymatrix<ntype,N,N> rotation_matrix_axis_cossin( myvec<ntype,N> u, ntype c, ntype s )
{
  u.normalize();
  mymatrix<ntype,N,N> eye;
  eye.set_identity();
  return (eye*=c) + (cross_product_matrix(u)*=s) + (outer_product(u,u)*=(1.0-c));
}
using matrix3d=mymatrix<double,3,3>;

//-----------------------------------------------------------------//
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
  const mymatrix<ntype,1,1> inverse() const
  {
    mymatrix<ntype,1,1> m;
    if( M[0][0]==0.0) { cout << "[ERROR: this matrix is not invertible.]"<<endl; show(); exit(1); }
    m[0][0] = 1.0/M[0][0];
    return m;
  }
  mymatrix()
  {
  }
  ~mymatrix()
  {}
};

#endif
