/*
	Spherical harmonics
https://en.wikipedia.org/wiki/Table_of_spherical_harmonics
*/

#ifndef _YCOMPLEX_H_
#define _YCOMPLEX_H_
#include<cstdlib>
#include<complex>
#include<cmath>
#include "utility.hpp"
using namespace std;

template<class ntype>
complex<ntype> Y(int l, int m, myvec<ntype,3> R)
{
  ntype r=R.norm();
  return Y(l,m,R,r);
}

template<class ntype>
complex<ntype> Y(int l, int m, myvec<ntype,3> R, ntype r)
{
  ntype pi=M_PI, x=R[0], y=R[1], z=R[2];
  if(!(l>=0)) { cout << "[ERROR: negative angular momentum selected in Y(): l = " << l << "]\n"; exit(1); }
  if(!(l<=6)) { cout << "[ERROR: l>6 not yet implemented in Y(): l = " << l << "]\n"; exit(1); }
  if(!(l>=abs(m))) { cout << "[ERROR: condition l>=|lz| not satisfied in Y(): l = " << l << ", lz = " << m << "]\n"; exit(1); }
  if(r==0.0) { cout << "[ERROR: cannot compute Y() of a null direction]\n"; exit(1); }
  ntype c = z/r, c2,c3,c4,c5,c6; // cos(theta)
  complex<ntype> sp = (x + y*1i)/r; // exp(i*phi)*sin(theta)
  complex<ntype> sm = conj(sp); // exp(-i*phi)*sin(theta)
  switch(l)
  {
    case 0: return 0.5/sqrt(pi);
    case 1:
      switch(m)
      {
        case -1: return 0.5*sqrt(0.5*3/pi) * sm;
        case  0: return 0.5*sqrt(3/pi)     * c;
        case  1: return 0.5*sqrt(0.5*3/pi) * sp;
      }
    case 2:
      switch(m)
      {
        case -2: return  0.25*sqrt(0.5*15/pi) * sm*sm;
        case -1: return  0.5*sqrt(0.5*15/pi)  * sm*c;
        case  0: return  0.25*sqrt(5/pi)      * (3*c*c - 1);
        case  1: return -0.5*sqrt(0.5*15/pi)  * sp*c;
        case  2: return  0.25*sqrt(0.5*15/pi) * sp*sp;
      }
    case 3:
      switch(m)
      {
        case -3: return  0.125*sqrt(35/pi)     * sm*sm*sm;
        case -2: return  0.25*sqrt(0.5*105/pi) * sm*sm*c;
        case -1: return  0.125*sqrt(21/pi)     * sm*(5*c*c - 1);
        case  0: return  0.25*sqrt(7/pi)       * (5*c*c*c - 3*c);
        case  1: return -0.125*sqrt(21/pi)     * sp*(5*c*c - 1.0);
        case  2: return  0.25*sqrt(0.5*105/pi) * sp*sp*c;
        case  3: return -0.125*sqrt(35/pi)     * sp*sp*sp;
      }
    case 4:
      c2 = c*c;
      c3 = c2*c;
      c4 = c3*c;
      switch(m)
      {
        case -4: return  0.0625*3*sqrt(0.5*35/pi) * sm*sm*sm*sm;
        case -3: return  0.125*3*sqrt(35/pi)      * sm*sm*sm*c;
        case -2: return  0.125*3*sqrt(0.5*5/pi)   * sm*sm*(7*c2 - 1);
        case -1: return  0.125*3*sqrt(5/pi)       * sm*(7*c3 - 3*c);
        case  0: return  0.0625*3/sqrt(pi)        * (35*c4 - 30*c2 + 3);
        case  1: return -0.125*3*sqrt(5/pi)       * sp*(7*c3 - 3*c);
        case  2: return  0.125*3*sqrt(0.5*5/pi)   * sp*sp*(7*c2 - 1);
        case  3: return -0.125*3*sqrt(35/pi)      * sp*sp*sp*c;
        case  4: return  0.0625*3*sqrt(0.5*35/pi) * sp*sp*sp*sp;
      }
    case 5:
      c2 = c*c;
      c3 = c2*c;
      c4 = c3*c;
      c5 = c4*c;
      switch(m)
      {
        case -5: return  0.03125*3*sqrt(77/pi)     * sm*sm*sm*sm*sm;
        case -4: return  0.0625*3*sqrt(0.5*385/pi) * sm*sm*sm*sm*c;
        case -3: return  0.0625*sqrt(385/pi)       * sm*sm*sm*(9*c2 - 1);
        case -2: return  0.125*sqrt(0.5*1155/pi)   * sm*sm*(3*c3 - c);
        case -1: return  0.0625*sqrt(0.5*165/pi)   * sm*(21*c4 - 14*c2 + 1);
        case  0: return  0.0625*sqrt(11/pi)        * (63*c5 - 70*c3 + 15*c);
        case  1: return -0.0625*sqrt(0.5*165/pi)   * sp*(21*c4 - 14*c2 + 1);
        case  2: return  0.125*sqrt(0.5*1155/pi)   * sp*sp*(3*c3 - c);
        case  3: return -0.0625*sqrt(385/pi)       * sp*sp*sp*(9*c2 - 1);
        case  4: return  0.0625*3*sqrt(0.5*385/pi) * sp*sp*sp*sp*c;
        case  5: return -0.03125*3*sqrt(77/pi)     * sp*sp*sp*sp*sp;
      }
    case 6:
      c2 = c*c;
      c3 = c2*c;
      c4 = c3*c;
      c5 = c4*c;
      c6 = c5*c;
      switch(m)
      {
        case -6: return  0.015625*sqrt(3003/pi)    * sm*sm*sm*sm*sm*sm;
        case -5: return  0.03125*3*sqrt(1001/pi)   * sm*sm*sm*sm*sm*c;
        case -4: return  0.03125*3*sqrt(0.5*91/pi) * sm*sm*sm*sm*(11*c2 - 1);
        case -3: return  0.03125*sqrt(1365/pi)     * sm*sm*sm*(11*c3 - 3*c);
        case -2: return  0.015625*sqrt(1365/pi)    * sm*sm*(33*c4 - 18*c2 + 1);
        case -1: return  0.0625*sqrt(0.5*273/pi)   * sm*(33*c5 - 30*c3 + 5*c);
        case  0: return  0.03125*sqrt(13/pi)       * (231*c6 - 315*c4 + 105*c2 - 5);
        case  1: return -0.0625*sqrt(0.5*273/pi)   * sp*(33*c5 - 30*c3 + 5*c);
        case  2: return  0.015625*sqrt(1365/pi)    * sp*sp*(33*c4 - 18*c2 + 1);
        case  3: return -0.03125*sqrt(1365/pi)     * sp*sp*sp*(11*c3 - 3*c);
        case  4: return  0.03125*3*sqrt(0.5*91/pi) * sp*sp*sp*sp*(11*c2 - 1);
        case  5: return -0.03125*3*sqrt(1001/pi)   * sp*sp*sp*sp*sp*c;
        case  6: return  0.015625*sqrt(3003/pi)    * sp*sp*sp*sp*sp*sp;
       }
    default: break;
  }
  cout << "[WARNING: you should not reach this line. There's a bug in the code]\n";
  return 10000.0; // useless return
}
# endif
