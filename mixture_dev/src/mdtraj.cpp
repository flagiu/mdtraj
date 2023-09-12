#include "all_headers.hpp"
using namespace std;
using ntype=double;
using ptype=particle<ntype>;

//-----------------Main----------------------------------------//

int main(int argc, char* argv[])
{
  randnum.rseed();
  Trajectory<ntype,ptype> traj;
  traj.init(argc, argv);
  traj.run();
}

