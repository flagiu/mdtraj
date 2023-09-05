using namespace std;
//--------------------- Angular Distribution Function (distribution of cos(theta_ijk) ----------------------------------//
// must be preceded by computation of nearest neighbours
template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
init_adf()
{
      adf_nbins = int(floor( 2.0 / adf_binw )); // -1 < cos(angle) < 1
      adf_bins.resize(adf_nbins);
      adf.resize(adf_nbins);
      adf_ave.resize(adf_nbins);
      adf2_ave.resize(adf_nbins);
      ss.str(std::string()); ss << s_adf << tag << ".traj"; fout.open(ss.str(), ios::out);
      fout << "# First block: cos(angle); other blocks: ADF for each frame )\n";
      for( auto i=0; i<adf_nbins; i++){
        adf_bins[i] = -1.0 + (i+0.5)*adf_binw; // take the center of the bin for histogram
        fout << adf_bins[i] << endl;
        adf_ave[i] = 0.0;
        adf2_ave[i] = 0.0;
      }
      fout.close();
      ss.str(std::string()); ss << s_adf << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "# cos(angle), ADF, ADF error.\n";
      fout.close();
}

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
compute_adf(int frameidx)
{
    int i,j,k, a,b, bin, counts=0;
    vec rij, rik;
    ntype rijSq, rikSq, costheta;
    for(a=0; a<adf_nbins; a++) adf[a] = 0.0;
    if(debug) cout << "*** ADF computation for timestep " << timestep << " STARTED ***\n";
    for(i=0;i<N;i++){
      for(a=1;a<ps[i].neigh_list[0].size();a++){
          j = ps[i].neigh_list[0][a];
          if(j>i) continue; // avoid double counting!
          rij = ps[i].rij_list[0][a];
          rijSq = ps[i].rijSq_list[0][a];
          for(b=0;b<a;b++){
            k = ps[i].neigh_list[0][b]; //useless?
            rik = ps[i].rij_list[0][b];
            rikSq = ps[i].rijSq_list[0][b];
            costheta = (rij*rik) / sqrt(rijSq*rikSq);
            bin = int(floor( (costheta+1.0)/adf_binw)); // min value is -1.0 !
            if(bin<adf_nbins){
              adf[bin] += 1.0;
              counts++;
            }
          }
      }
    }
    ss.str(std::string()); ss << s_adf << tag << ".traj"; fout.open(ss.str(), ios::app);
    fout << endl; // start new block
    for(k=0; k<adf_nbins; k++){
      adf[k] /= (counts*adf_binw); // ci sta la larghezza del bin?
      fout << adf[k] << endl;
      adf_ave[k] += adf[k];
      adf2_ave[k] += adf[k]*adf[k];
    }
    fout.close();
    if(frameidx == (nframes-1)){
        ss.str(std::string()); ss << s_adf << tag << ".ave"; fout.open(ss.str(), ios::app);
        for(k=0; k<adf_nbins; k++){
          adf_ave[k] /= (frameidx+1);
          adf2_ave[k] /= (frameidx+1);
          fout << adf_bins[k] << " " << adf_ave[k] << " " << sqrt( (adf2_ave[k]-adf_ave[k]*adf_ave[k])/(frameidx+1 -1) ) << endl;
        }
        fout.close();
        if(debug) cout << "Average ADF printed to file\n";
    }
    if(debug) cout << "*** ADF computation for timestep " << timestep << " ENDED ***\n\n";
    return;
}
