using namespace std;
//------------------------------- Radial Distribution Function ----------------------------------//
template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
init_rdf()
{
  int t1,t2,a;
  ntype r, shell1, shell2, normalization;
  rdf_nbins = int(floor( 0.5*L.abs().min() / rdf_binw ));
  //rdf_binw = 0.5*(L.abs().min()) / (rdf_nbins-1) ;
  rdf_bins.resize(rdf_nbins);
  shell1 = 0.0;
  ss.str(std::string()); ss << s_rdf << tag << ".traj"; fout.open(ss.str(), ios::out);
  fout << "# First block: r; Other blocks: g_00(r), g_01(r), ... for each frame.\n";
  for( auto i=0; i<rdf_nbins; i++){
    rdf_bins[i] = (i+0.5)*rdf_binw; // take the center of the bin for histogram
    fout << rdf_bins[i] << endl;
  }

  rdf_norm.resize(nTypePairs);
  rdf.resize(nTypePairs);
  rdf_ave.resize(nTypePairs);
  rdf2_ave.resize(nTypePairs);

  for(t1=0;t1<nTypes;t1++)
  {
    for(t2=t1;t2<nTypes;t2++)
    {
      a = types2int(t1,t2,nTypes);
      if(debug) cout << " Inizializzo la coppia (" << t1 <<","<<t2<<")-->"<<  a << "\n";
      rdf_norm[a].resize(rdf_nbins);
      rdf[a].resize(rdf_nbins);
      rdf_ave[a].resize(rdf_nbins);
      rdf2_ave[a].resize(rdf_nbins);
      if(t1!=t2) normalization =     Nt[t2]/V * 4.0*M_PI/3.0 * Nt[t1]/2.0 * 2.0; // multiply by 2 because it's off-diagonal
      else       normalization = (Nt[t2]-1)/V * 4.0*M_PI/3.0 * Nt[t1]/2.0;
      shell1 = 0.0;
      for( auto i=0; i<rdf_nbins; i++)
      {
        r = (i+1)*rdf_binw; // take the upper side of the bin for shells
        shell2 = r*r*r;
        rdf_norm[a][i] = normalization*(shell2 - shell1); // normalize to ideal gas radial distribution
        shell1 = shell2;
        rdf_ave[a][i] = 0.0;
        rdf2_ave[a][i] = 0.0;
      }
    }
  }
  fout.close();
  ss.str(std::string()); ss << s_rdf << tag << ".ave"; fout.open(ss.str(), ios::out);
  fout << "# r | g_00(r), g_01(r), ... | error for each g(r).\n";
  fout.close();
}

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
compute_rdf(int frameidx)
{
  int i,j, k, tp, ti,tj;
  vec rij;
  ntype r,r_mic;
  for(tp=0;tp<nTypePairs;tp++){
    for(k=0; k<rdf_nbins; k++){
      rdf[tp][k] = 0.0;
    }
  }
  if(debug) cout << "*** RDF computation for timestep " << timestep << " STARTED ***\n";
  for(i=0;i<N-1;i++){
    for(j=i+1;j<N;j++){// i<j
      rij = ps[j].r - ps[i].r; // real distance
      r = rij.norm();
      rij = mic(box, boxInv, rij); // first periodic image
      r_mic = rij.norm();
      if( r_mic < r ) r = r_mic; // if closer, choose first periodic image
      // if(debug) cout << "  PBC for (i,j)=(" <<i<<","<<j<<") : " << r_true << " --> " << r << endl;
      k = int(floor( r/rdf_binw));
      //if(debug) cout << "  r="<<r<<", k="<<k<<endl;
      if(k<rdf_nbins){
        ti = ps[i].label;
        tj = ps[j].label;
        tp = types2int(ti,tj,nTypes);
        //if(debug) cout << "  updatingr g("<<ti<<","<<tj<<",r)\n";
        rdf[tp][k] += 1.0;
      }
    }
  }
  ss.str(std::string()); ss << s_rdf << tag << ".traj"; fout.open(ss.str(), ios::app);
  fout << endl; // start new block
  for(k=0; k<rdf_nbins; k++)
  {
    for(tp=0;tp<nTypePairs;tp++)
    {
      rdf[tp][k] /= rdf_norm[tp][k];
      fout << rdf[tp][k];
      rdf_ave[tp][k] += rdf[tp][k];
      rdf2_ave[tp][k] += rdf[tp][k]*rdf[tp][k];
      if(tp<nTypePairs-1) fout << " ";
      else fout << endl;
    }
  }
  fout.close();
  if(frameidx == (nframes-1)){
      ss.str(std::string()); ss << s_rdf << tag << ".ave"; fout.open(ss.str(), ios::app);
      for(k=0; k<rdf_nbins; k++)
      {
        fout << rdf_bins[k] << " "; // print bin
        for(tp=0;tp<nTypePairs;tp++)
        {
          rdf_ave[tp][k] /= (frameidx+1);
          fout << rdf_ave[tp][k] << " "; // print average g(r)
        }
        for(tp=0;tp<nTypePairs;tp++)
        {
          rdf2_ave[tp][k] /= (frameidx+1);
          fout << sqrt( (rdf2_ave[tp][k]-rdf_ave[tp][k]*rdf_ave[tp][k])/(frameidx+1 -1) ); // print st.dev. of average g(r)
          if(tp<nTypePairs-1) fout << " ";
          else fout << endl;
        }
      }
      fout.close();
      if(debug) cout << "Average RDF printed to file\n";
  }
  if(debug) cout << "*** RDF computation for timestep " << timestep << " ENDED ***\n\n";
  return;
}
