using namespace std;
//--------------------- ALTBC ----------------------------------//
// must be preceded by computation of nearest neighbours (?)
template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
init_altbc()
{
      altbc_nbins = int(floor( (cutoff[0]-altbc_rmin)/altbc_binw ));
      altbc_bins.resize(altbc_nbins);
      altbc.resize(altbc_nbins*altbc_nbins);
      altbc_ave.resize(altbc_nbins*altbc_nbins);
      altbc2_ave.resize(altbc_nbins*altbc_nbins);
      altbc_cos = cos( altbc_angle/180. * M_PI ); // should be close to 1.0 for 'collinear bonds'
      if(debug) cout << "ALTBC cos threshold = " << altbc_cos<<endl;
      ss.str(std::string()); ss << s_altbc << tag << ".traj"; fout.open(ss.str(), ios::out);
      fout << "# 1st block: radial distance; other blocks: ALTBC 2D matrix for each frame; # deviation >= " << altbc_angle << " degrees\n";
      for( auto i=0; i<altbc_nbins; i++){
        altbc_bins[i] = altbc_rmin + (i+0.5)*altbc_binw; // take the center of the bin for histogram
        fout << altbc_bins[i] << endl;
      }
      fout << endl; // end of block
      for( auto i=0; i<altbc.length(); i++){
        altbc_ave[i] = 0.0;
        altbc2_ave[i] = 0.0;
      }
      fout.close();
      ss.str(std::string()); ss << s_altbc << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "# 1st block: radial distance; 2nd block: ALTBC 2D matrix; # deviation >= " << altbc_angle << " degrees\n";
      fout.close();
}

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
compute_altbc(int frameidx)
{
    int i,j,k,k0,k1, a,b, bin0, bin1, counts=0;
    vec rij, rik;
    ntype rijSq, rikSq, costheta, rijNorm, rikNorm, norm_factor= V*V/(N*(N-1)*(N-2));
    for(a=0; a<altbc.length(); a++) altbc[a] = 0.0;
    if(debug) cout << "*** altbc computation for timestep " << timestep << " STARTED ***\n";
    for(i=0;i<N;i++){
      for(a=1;a<ps[i].neigh_list[0].size();a++){
          j = ps[i].neigh_list[0][a];
          //if(j>i) continue; // avoid double counting! WRONG
          rij = ps[i].rij_list[0][a];
          rijSq = ps[i].rijSq_list[0][a];
          rijNorm = sqrt(rijSq);
          if(rijNorm < altbc_rmin) continue;
          for(b=0;b<a;b++){
            k = ps[i].neigh_list[0][b]; //useless?
            rik = ps[i].rij_list[0][b];
            rikSq = ps[i].rijSq_list[0][b];
            rikNorm = sqrt(rikSq);
            if(rikNorm < altbc_rmin) continue;
            costheta = (rij*rik) / (rijNorm*rikNorm);
            bin0 = int(floor( (rijNorm-altbc_rmin)/altbc_binw ));
            bin1 = int(floor( (rikNorm-altbc_rmin)/altbc_binw ));
            if(bin0<altbc_nbins && bin1<altbc_nbins && fabs(costheta) >= altbc_cos){
              altbc[bin0 + altbc_nbins*bin1] += 1.0;
              counts++;
            }
          }
      }
    }
    ss.str(std::string()); ss << s_altbc << tag << ".traj"; fout.open(ss.str(), ios::app);
    fout << endl; // start new block
    for(k0=0; k0<altbc_nbins; k0++){
      for(k1=0; k1<altbc_nbins; k1++){
        k = k0 + altbc_nbins*k1;
        if(counts>0) altbc[k] /= (counts*altbc_binw*altbc_binw); // ci sta l'area del bin??
        altbc[k] *= norm_factor;
        fout << altbc[k] << " "; // 2D matrix
        altbc_ave[k] += altbc[k];
        altbc2_ave[k] += altbc[k]*altbc[k];
      }
      fout << endl;
    }
    fout << endl; // end of block
    fout.close();
    if(frameidx == (nframes-1)){
        ss.str(std::string()); ss << s_altbc << tag << ".ave"; fout.open(ss.str(), ios::app);
        for(k0=0; k0<altbc_nbins; k0++){
          fout << altbc_bins[k0] << endl; // bins
        }
        fout << endl; // end of block
        for(k0=0; k0<altbc_nbins; k0++){
          for(k1=0; k1<altbc_nbins; k1++){
            k = k0 + altbc_nbins*k1;
            altbc_ave[k] /= (frameidx+1);
            altbc2_ave[k] /= (frameidx+1);
            fout << altbc_ave[k] << " "; // average
          }
          fout << endl;
        }
        /*
        fout << endl; // end of block
        for(k0=0; k0<altbc_nbins; k0++){
          for(k1=0; k1<altbc_nbins; k1++){
            k = k0 + altbc_nbins*k1;
            fout << sqrt( (altbc2_ave[k]-altbc_ave[k]*altbc_ave[k])/(frameidx+1 -1) ) << " "; // error
          }
          fout << endl;
        }
        */
        fout.close();
        if(debug) cout << "Average ALTBC printed to file\n";
    }
    if(debug) cout << "*** ALTBC computation for timestep " << timestep << " ENDED ***\n\n";
    return;
}
