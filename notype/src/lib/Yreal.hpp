/*
	Spherical harmonics
*/

#ifndef _Y_REAL_H_
#define _Y_REAL_H_
#include<cstdlib>
#include "particle.hpp"
#include "utility.hpp"

template<class ntype>
ntype Y(int l, int m, myvec<ntype,3> r)
{
  if(!(l>=0)) { cout << "[ERROR: negative angular momentum selected in Y(): l = " << l << "]\n"; exit(1); }
  if(!(l<=4)) { cout << "[ERROR: l>4 not yet implemented in Y(): l = " << l << "]\n"; exit(1); }
  if(!(l>=abs(m))) { cout << "[ERROR: condition l>=|lz| not satisfied in Y(): l = " << l << ", lz = " << m << "]\n"; exit(1); }
  if(l==0) { cout << "[Warning: bond orientation order parameter with l=0 is always equal to 1]"; }
  ntype pi=M_PI, x=r[0], y=r[1], z=r[2], r1,r2,r3,r4, r5;
  switch(l)
  {
    case 0: return 0.5/sqrt(pi);
    case 1:
      switch(m)
      {
        case -1: return 0.5*sqrt(3/pi) * y/r.norm();
        case  0: return 0.5*sqrt(3/pi) * z/r.norm();
        case  1: return 0.5*sqrt(3/pi) * x/r.norm();
      }
    case 2:
      switch(m)
      {
        case -2: return  0.5*sqrt(15/pi) * x*y/r.sq();
        case -1: return  0.5*sqrt(15/pi) * y*z/r.sq();
        case  0: return  0.25*sqrt(5/pi) * (3*z*z/r.sq() - 1.0);
        case  1: return  0.5*sqrt(15/pi) * x*z/r.sq();
        case  2: return 0.25*sqrt(15/pi) * (x*x - y*y)/r.sq();
      }
    case 3:
      r2 = r.sq();
      r1 = sqrt(r2);
      r3 = r1*r2;
      switch(m)
      {
        case -3: return 0.25*sqrt(0.5*35/pi) * y*(3*x*x - y*y)/r3;
        case -2: return 0.5*sqrt(105/pi)     * x*y*z/r3;
        case -1: return 0.25*sqrt(0.5*21/pi) * y*(5*z*z/r3 - 1.0);
        case  0: return 0.25*sqrt(7/pi)      * z*(5*z*z/r3 - 3/r1);
        case  1: return 0.25*sqrt(0.5*21/pi) * x*(5*z*z/r3 - 1.0);
        case  2: return 0.5*sqrt(105/pi)     * (x*x - y*y)*z/r3;
        case  3: return 0.25*sqrt(0.5*35/pi) * x*(x*x - 3*y*y)/r3;
      }
    case 4:
      r2 = r.sq();
      r4 = r2*r2;
      switch(m)
      {
        case -4: return 0.25*3*sqrt(35/pi)     * x*y*(x*x - y*y)/r4;
        case -3: return 0.25*3*sqrt(0.5*35/pi) * y*(3*x*x - y*y)*z/r4;
        case -2: return 0.25*3*sqrt(5/pi)      * x*y*(7*z*z - r2)/r4;
        case -1: return 0.25*3*sqrt(0.5*5/pi)  * y*(7*z*z - 3*r2)*z/r4;
        case  0: return 0.0625*3/sqrt(pi)      * (z*z*(35*z*z - 30*r2)/r4 + 3.0);
        case  1: return 0.25*3*sqrt(0.5*5/pi)  * (x*x - y*y)*z*(7*z*z - r2)/r4;
        case  2: return 0.125*3*sqrt(5/pi)     * (x*x - y*y)*(7*z*z - r2)/r4;
        case  3: return 0.25*3*sqrt(0.5*35/pi) * x*(x*x - 3*y*y)*z/r4;
        case  4: return 0.0625*3*sqrt(35/pi)   * (x*x*(x*x - 3*y*y) - y*y*(3*x*x - y*y))/r4;
      }
    case 5:
      r2 = r.sq();
      r1 = sqrt(r2);
      r3 = r1*r2;
      r5 = r3*r2;
      switch(m)
      {
        case -5: return 0.03125*3*sqrt(77/pi)     * x*y*(x*x - y*y)/r4;
      }
    default: break;
  }
  return 10000.0; // useless return
}
#endif
