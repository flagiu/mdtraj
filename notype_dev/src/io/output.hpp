using namespace std;
//------------ Priting output trajectories of various formats --------------------//

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
init_box() {
  ss.str(std::string()); ss << s_box << tag << ".dat"; fout.open(ss.str(), ios::out);
  if(!out_alphanes) fout << "#Ax Bx Cx Ay By Cy Az Bz Cz\n";
  fout.close();
}

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
print_box() {
  if(out_alphanes)
  {
    fout.open("box.dat", ios::app);
    // Ax, Bx, Cx, By, Cy, Cz (the only non-zero entries)
    fout<< box[0][0] <<" "<< box[0][1] <<" "<< box[0][2] <<" "<< box[1][1] <<" "<< box[1][2] <<" "<< box[2][2] <<endl;
  }
  else
  {
    ss.str(std::string()); ss << s_box << tag << ".dat"; fout.open(ss.str(), ios::app);
    box.write(fout);
  }
  fout.close();
}

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
init_out_xyz() {
  ss.str(std::string()); ss << s_out << tag << ".xyz"; fout.open(ss.str(), ios::out);
  fout.close();
}

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
print_out_xyz() {
  ss.str(std::string()); ss << s_out << tag << ".xyz"; fout.open(ss.str(), ios::app);
  fout << N << endl;
  fout << "Atoms. Timestep: " << timestep << endl;
  for(auto &p : ps)
  {
    if(pbc_out) p.r = mic(box, boxInv, p.r); // apply PBC to the position
    p.write_xyz(fout);
  }
  fout.close();
}


template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
init_out_alphanes() {
  fout.open("box.dat", ios::out);
  fout.close();
  fout.open("pos.dat", ios::out);
  fout.close();
  /*

  fout.open("force.dat", ios::out);
  fout.close();

  fout.open("energy.dat", ios::out);
  fout.close();
  */
}

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
print_out_alphanes() {
  print_box();
  fout.open("pos.dat", ios::app);
  for(auto i=0;i<ps.size();i++) {
    fout << ps[i].r[0] <<" "<< ps[i].r[1] <<" "<< ps[i].r[2];
    if(i<ps.size()-1) fout << " ";
    else fout << endl;
  }
  fout.close();
  /* Forces and energies are not yet implemented!
  fout.open("force.dat", ios::app);
  for(auto &p : ps) {
    fout << p.f[0] <<" "<< p.f[1] <<" "<< p.f[2] <<endl;
  }
  fout.close();

  fout.open("energy.dat", ios::app);
  fout << energy_per_atom << endl;
  fout.close();
  */
}
