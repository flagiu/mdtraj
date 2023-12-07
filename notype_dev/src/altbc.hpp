#ifndef _ALTBC_H_
#define _ALTBC_H_
using namespace std;
//--------------------- ALTBC ----------------------------------//
// must be preceded by computation of nearest neighbours (?)

template <class ntype, class ptype>
class ALTBC_Calculator
{
  using vec=myvec<ntype,3>;
  using mat=mymatrix<ntype,3,3>;
  private:
    ntype r_min, binw, r_max, angle_th, cos_th;
    int nbins, nbins2;
    vecflex<ntype> bins, norm, value, ave, ave2;
    string string_out, myName, tag;
    fstream fout;
    stringstream ss;
  public:
    ALTBC_Calculator(){
      myName = "ALTBC";
    }
    virtual ~ALTBC_Calculator(){}

    void init(ntype r_min_, ntype binw_, ntype r_max_, ntype angle_th_, int N, ntype V, string string_out_, string tag_, bool debug)
    {
      int i,j, idx;
      ntype shell1, shell2;
      string_out = string_out_;
      tag = tag_;
      r_min = r_min_;
      binw = binw_;
      r_max = r_max_;
      nbins = int(floor( (r_max-r_min)/binw ));
      nbins2 = nbins*nbins;
      angle_th = angle_th_;
      cos_th = cos( angle_th/180. * M_PI ); // should be close to 1.0 for 'collinear bonds'
      if(debug) cout << myName <<": cosine threshold = " << cos_th<<endl;

      bins.resize(nbins);
      norm.resize(nbins2);
      value.resize(nbins2);
      ave.resize(nbins2);
      ave2.resize(nbins2);

      ss.str(std::string()); ss << string_out << tag << ".traj"; fout.open(ss.str(), ios::out);
      fout << "# 1st block: radial distance; other blocks: ALTBC 2D matrix for each frame; # deviation >= " << angle_th << " degrees\n";
      for(i=0; i<nbins; i++)
      {
        bins[i] = r_min + (i+0.5)*binw; // take the center of the bin for histogram
        fout << bins[i] << endl;
      }
      fout << endl; // end of block
      fout.close();

      ntype NVfactor = ((ntype)N) *(((ntype)(N-1))/V) * (((ntype)(N-2))/V);
      for(i=0; i<nbins; i++)
      {
        shell1 = 4.0 * M_PI * bins[i]*bins[i] * binw;
        for(j=0;j<nbins;j++)
        {
          shell2 = 4.0 * M_PI * bins[j]*bins[j] * binw;
          idx = nbins*i + j;
          norm[idx] = shell1*shell2 * NVfactor ; // prepare normalization
          ave[idx]=ave2[idx]=0.0; // reset averages
        }
      }

      ss.str(std::string()); ss << string_out << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "# 1st block: radial distance; 2nd block: ALTBC 2D matrix; # deviation >= " << angle_th << " degrees\n";
      fout.close();
    }

    void compute(int frameidx, int nframes, int timestep, vector<ptype> ps, bool debug)
    {
      const int N=ps.size();
      int i,j,k,k0,k1, a,b, bin0, bin1, counts;
      vec rij, rik;
      ntype rijSq, rikSq, costheta, rijNorm, rikNorm;
      for(a=0; a<nbins2; a++) value[a] = 0.0;
      if(debug) cout << "*** "<<myName<<" computation for timestep " << timestep << " STARTED ***\n";

      counts=0;
      for(i=0;i<N;i++)
      {
        for(a=1;a<ps[i].neigh_list[0].size();a++)
        {
          j = ps[i].neigh_list[0][a];
          //if(j>i) continue; // avoid double counting! WRONG
          rij = ps[i].rij_list[0][a];
          rijSq = ps[i].rijSq_list[0][a];
          rijNorm = sqrt(rijSq);
          if(rijNorm < r_min) continue;
          for(b=0;b<a;b++)
          {
            k = ps[i].neigh_list[0][b]; //useless?
            rik = ps[i].rij_list[0][b];
            rikSq = ps[i].rijSq_list[0][b];
            rikNorm = sqrt(rikSq);
            if(rikNorm < r_min) continue;
            costheta = (rij*rik) / (rijNorm*rikNorm);
            bin0 = int(floor( (rijNorm-r_min)/binw ));
            bin1 = int(floor( (rikNorm-r_min)/binw ));
            if(bin0<nbins && bin1<nbins && fabs(costheta) >= cos_th)
            {
              value[bin0 + nbins*bin1] += 1.0;
              counts++;
              value[bin1 + nbins*bin0] += 1.0;
              counts++;
            }
          }
        }
      }
      // Print current histogram and add it to the average
      ss.str(std::string()); ss << string_out << tag << ".traj"; fout.open(ss.str(), ios::app);
      fout << endl; // start new block
      for(k0=0; k0<nbins; k0++)
      {
        for(k1=0; k1<nbins; k1++)
        {
          k = k0 + nbins*k1;
          if(counts>0) value[k] /= (counts*binw*binw); // ci sta l'area del bin??
          value[k] /= norm[k];
          fout << value[k] << " "; // 2D matrix
          ave[k] += value[k];
          ave2[k] += value[k]*value[k];
        }
        fout << endl;
      }
      fout << endl; // end of block
      fout.close();

      if(frameidx == (nframes-1)) // if last frame: print final results
      {
        ss.str(std::string()); ss << string_out << tag << ".ave"; fout.open(ss.str(), ios::app);
        for(k0=0; k0<nbins; k0++){
          fout << bins[k0] << endl; // bins
        }
        fout << endl; // end of block
        for(k0=0; k0<nbins; k0++){
          for(k1=0; k1<nbins; k1++){
            k = k0 + nbins*k1;
            ave[k] /= nframes;
            ave2[k] /= nframes;
            fout << ave[k] << " "; // average
          }
          fout << endl;
        }
        /* // do NOT print error bars... who does even look at them?
        fout << endl; // end of block
        for(k0=0; k0<nbins; k0++){
          for(k1=0; k1<nbins; k1++){
            k = k0 + nbins*k1;
            fout << sqrt( (ave2[k]-ave[k]*ave[k])/(nframes -1) ) << " "; // error
          }
          fout << endl;
        }
        */
        fout.close();
        if(debug) cout << "Average "<<myName<<" printed to file\n";
      }
      if(debug) cout << "*** "<<myName<<" computation for timestep " << timestep << " ENDED ***\n\n";
      return;
    }
};
#endif
