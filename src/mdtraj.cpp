#include "lib/mdtraj.hpp"

//-----------------Main----------------------------------------//

int main(int argc, char* argv[])
{
  randnum.rseed();
  Trajectory<ntype,ptype> traj;
  traj.init(argc, argv);
//  traj.set_L(15.09, 15.68, 15.80);
  traj.run();
}

