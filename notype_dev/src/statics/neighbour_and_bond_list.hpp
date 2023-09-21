using namespace std;

//---- "init" functions: prepare arrays, averages and headers for output files
//---- "compute" functions: compute quantities for each frame and total average

//------------------------------ Neighbour and Bond list --------------------------------------------//

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
init_neigh(){
  int u;
  for(u=0;u<maxshell;u++){
    if(neigh[u].length()!=N) neigh[u].resize(N);
  }
}

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
build_neigh() {
  int u,i,j;
  vec rij, rij_mic;
  ntype rijSq, rijSq_mic;
  if(debug) cout << "\n*** Neighbours computation STARTED ***\n";
  //---- Reset counters and lists ----//
  for(u=0;u<maxshell;u++) {
    bond_list[u].clear();
  }
  for(i=0;i<N;i++){
      for(u=0;u<maxshell;u++){
        ps[i].neigh_list[u].clear();
        ps[i].rij_list[u].clear();
        ps[i].rijSq_list[u].clear();
        neigh[u].set(i, 0.0);
      }
  }
  if(debug) cout << " * Reset counters and lists DONE\n";
    //---- Build neighbour list (and save rij vectors) ----//
    for(i=0;i<N;i++){
      for(j=i+1;j<N;j++){
        rij = ps[j].r - ps[i].r;
        rijSq = rij.sq();
        rij_mic = mic(box, boxInv, rij); // first periodic image
        rijSq_mic = rij_mic.sq();
        if(rijSq_mic < rijSq){ // if closer, choose first periodic image
          rijSq = rijSq_mic;
          rij = rij_mic;
        }
        for(u=0;u<maxshell;u++){
          if(rijSq <= cutoffSq[u]){
            bond_list[u].push_back( ij2int(i,j,N) );

            ps[i].neigh_list[u].push_back(j);
            ps[j].neigh_list[u].push_back(i);

            ps[i].rij_list[u].push_back(   rij);
            ps[j].rij_list[u].push_back(-1*rij);

            ps[i].rijSq_list[u].push_back(rijSq);
            ps[j].rijSq_list[u].push_back(rijSq);
          }
        }
      }
    }
  if(debug) cout << " * Build neighbour and bond lists DONE\n";
  if(debug) cout << "*** Neighbour computation for timestep " << timestep << " ENDED ***\n\n";
  return;
}

//------------------------------------- Coordination Number ---------------------------------------//
template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
init_coordnum()
{
      if(neigh[0].length()!=N) neigh[0].resize(N);
      ss.str(std::string()); ss << s_coordnum << tag << ".dat"; fout.open(ss.str(), ios::out);
      fout << "#Timestep, Particle index, Coordination number. # cutoffs = ";
      for(auto  i=0;i<Nshells;i++) fout << cutoff[i]<<", ";
      fout << endl;
      fout.close();

      ss.str(std::string()); ss << s_coordnum << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "#Timestep, average coordination number, fluctuations. # cutoffs = ";
      for(auto  i=0;i<Nshells;i++) fout << cutoff[i]<<", ";
      fout << endl;
      fout.close();
}

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
compute_coordnum() {
  int i, j, k;
  ntype fval, rijSq;
  vec rij;
  for(i=0;i<N;i++) neigh[0].set(i,0.0); // start counter
  for(i=0;i<N;i++){
    for(k=0;k<ps[i].neigh_list[0].size();k++){ // search in 1st shell neighbour list
      j = ps[i].neigh_list[0][k];
      if(j>i) continue; // avoid double counting!
      rij = ps[i].rij_list[0][k];
      rijSq = ps[i].rijSq_list[0][k];
      // fval = fcut( rijSq/cutoffSq[0], p1half, p2half );    // smooth
      fval = ( rijSq <= cutoffSq[0] ? 1.0 : 0.0 );         // sharp
      neigh[0][i] += fval;
      neigh[0][j] += fval;
    }
  }
  ss.str(std::string()); ss << s_coordnum << tag << ".dat"; fout.open(ss.str(), ios::app);
  for(i=0;i<N;i++) fout << timestep << " " << i << " " << neigh[0][i] << endl;
  fout.close();
  ss.str(std::string()); ss << s_coordnum << tag << ".ave"; fout.open(ss.str(), ios::app);
  fout << timestep << " " << neigh[0].mean() << " " << neigh[0].std()/sqrt(N) << endl;
  fout.close();
}
