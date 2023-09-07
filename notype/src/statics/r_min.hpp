using namespace std;
//------------ Minimum atomic distance ----------------//
template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
init_rmin()
{
  ss.str(std::string()); ss << s_rmin << tag << ".dat"; fout.open(ss.str(), ios::out);
  fout << "# Minimum atomic distance. # cutoff1 = " << cutoff[0]<< endl;
  fout.close();
}

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
compute_rmin()
{
  ntype rSq, rminSq = cutoffSq[0];
  int i,j, k;
  if(debug) cout << "*** RMIN computation for timestep " << timestep << " STARTED ***\n";
  for(i=0;i<N;i++){
    for(k=0;k<ps[i].neigh_list[0].size();k++){
      j = ps[i].neigh_list[0][k];
      if(j>i) continue; // avoid double counting!
      rSq = ps[i].rijSq_list[0][k];
      if(rSq<rminSq) rminSq=rSq;
    }
  }
  ss.str(std::string()); ss << s_rmin << tag << ".dat"; fout.open(ss.str(), ios::app);
  fout << sqrt(rminSq) << endl;
  fout.close();
  if(debug) cout << "*** RMIN computation for timestep " << timestep << " ENDED ***\n";
  return;
}
