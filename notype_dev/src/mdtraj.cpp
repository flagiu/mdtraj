#include "lib/mdtraj.hpp"
#include "args.cpp"
#include "input.cpp"
#include "statics.cpp"
#include "dynamics.cpp"
using namespace std;
using ntype=double;
using ptype=particle<ntype>;

//-----------------Main----------------------------------------//

int main(int argc, char* argv[])
{
  randnum.rseed();
  Trajectory<ntype,ptype> traj;
  traj.init(argc, argv);
//  traj.set_L(15.09, 15.68, 15.80);
  traj.run();
}

