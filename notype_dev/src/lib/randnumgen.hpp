#ifndef _RNG_
#define _RNG_

#include<stdlib.h>
#include<random>
#include<chrono>
using namespace std;
// default engine: marsenne-twister
template <class ntype=double, typename engine=mt19937>
class randnumgen
{
  using RNGtype=engine;
  using unidist = uniform_real_distribution<ntype>;
  RNGtype rng;
  unidist *distro;  
public:
  void seed(unsigned int s)
    {
      seed_seq seq{s, s+1, s+2, s+3, s+4, s+5, s+6, s+7};  
      rng.seed(seq);
    }
  void rseed(void)
    {
      // crea una sequenza di seeding random utilizzamndo la classe random_device
      random_device r;
      typename engine::result_type s0 = r() ^ (
            (typename engine::result_type)
            chrono::duration_cast<chrono::seconds>(
                chrono::system_clock::now().time_since_epoch()
                ).count() +
            (typename engine::result_type)
            chrono::duration_cast<chrono::microseconds>(
                chrono::high_resolution_clock::now().time_since_epoch()
                ).count() );
      typename engine::result_type s[8] = {s0, r(), r(), r(), r(), r(), r(), r()};
      seed_seq seq(s,s+8);
      rng.seed(seq);
    }
  void discard( unsigned long n)
    {
      rng.discard(n);
    }
  ntype ranf(void)
    {
      return (*distro)(rng);
    }

  randnumgen()
  {
    rng.seed(0);
    distro = new unidist(0.0,1.0);
  }
  randnumgen(unsigned int s)
    {
      rng.seed(s);
      distro = new unidist(0.0,1.0);
    }
  ~randnumgen()
  {
    delete distro;
  }
};
static randnumgen<double> randnum;
#endif
