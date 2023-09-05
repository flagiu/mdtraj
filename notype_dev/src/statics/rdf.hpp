using namespace std;
//------------------------------- Radial Distribution Function ----------------------------------//
template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
init_rdf()
{
      rdf_nbins = int(floor( 0.5*L.abs().min() / rdf_binw ));
      //rdf_binw = 0.5*(L.abs().min()) / (rdf_nbins-1) ;
      rdf_norm.resize(rdf_nbins);
      rdf_bins.resize(rdf_nbins);
      rdf.resize(rdf_nbins);
      rdf_ave.resize(rdf_nbins);
      rdf2_ave.resize(rdf_nbins);
      ntype r, shell1, shell2, normalization = (N-1)/V * 4.0*M_PI/3.0 * N/2.0;
      shell1 = 0.0;
      ss.str(std::string()); ss << s_rdf << tag << ".traj"; fout.open(ss.str(), ios::out);
      fout << "# First block: radial distance; other blocks: g(r)\n";
      for( auto i=0; i<rdf_nbins; i++){
        rdf_bins[i] = (i+0.5)*rdf_binw; // take the center of the bin for histogram
        fout << rdf_bins[i] << endl;
        r = (i+1)*rdf_binw; // take the upper side of the bin for shells
        shell2 = r*r*r;
        rdf_norm[i] = normalization*(shell2 - shell1); // normalize to ideal gas radial distribution
        shell1 = shell2;
        rdf_ave[i] = 0.0;
        rdf2_ave[i] = 0.0;
      }
      fout.close();
      ss.str(std::string()); ss << s_rdf << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "# Radial distance, g(r), g(r) error.\n";
      fout.close();
}

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
compute_rdf(int frameidx)
{
    int i,j, k;
    vec rij;
    ntype r,r_mic;
    for(k=0; k<rdf_nbins; k++) rdf[k] = 0.0;
    if(debug) cout << "*** RDF computation for timestep " << timestep << " STARTED ***\n";
    for(i=0;i<N;i++){
      for(j=i+1;j<N;j++){ // i<j
         rij = ps[j].r - ps[i].r; // real distance
         r = rij.norm();
         rij = mic(box, boxInv, rij); // first periodic image
         r_mic = rij.norm();
         if( r_mic < r) r = r_mic; // if closer, choose first periodic image
         // if(debug) cout << "  PBC for (i,j)=(" <<i<<","<<j<<") : " << r_true << " --> " << r << endl;
         k = int(floor( r/rdf_binw ));
         if(k<rdf_nbins){
           rdf[k] += 1.0;
         }
      }
    }
    ss.str(std::string()); ss << s_rdf << tag << ".traj"; fout.open(ss.str(), ios::app);
    fout << endl; // start new block
    for(k=0; k<rdf_nbins; k++){
      rdf[k] /= rdf_norm[k];
      fout << rdf[k] << endl;
      rdf_ave[k] += rdf[k];
      rdf2_ave[k] += rdf[k]*rdf[k];
    }
    fout.close();
    if(frameidx == (nframes-1)){
        ss.str(std::string()); ss << s_rdf << tag << ".ave"; fout.open(ss.str(), ios::app);
        for(k=0; k<rdf_nbins; k++){
          rdf_ave[k] /= (frameidx+1);
          rdf2_ave[k] /= (frameidx+1);
          fout << rdf_bins[k] << " " << rdf_ave[k] << " " << sqrt( (rdf2_ave[k]-rdf_ave[k]*rdf_ave[k])/(frameidx+1 -1) ) << endl;
        }
        fout.close();
        if(debug) cout << "Average RDF printed to file\n";
    }
    if(debug) cout << "*** RDF computation for timestep " << timestep << " ENDED ***\n\n";
    return;
}
