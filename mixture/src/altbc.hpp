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
    int nbins, nbins2, max_nTypes, nTypes, nTypePairs;
    vecflex<ntype> bins, norm;
    vecflex< vecflex<ntype> > value, ave, ave2;
    string string_out, myName, tag;
    fstream fout;
    stringstream ss;
    bool debug, verbose;
    struct CollinearBonds {
      /*
       j--i----k
       - r_short = |rij|
       - ratio_shortlong = rij/rik
       - cos = cos(angle(jik))
      */
      ntype r_short, ratio_shortlong, cos;
    };
    struct PeierlsDistortion {
      static const int maxCollBonds=3;
      int numCollBonds;
      CollinearBonds collBonds[maxCollBonds];
    };
    vector<PeierlsDistortion> peierlsDists;

  public:
    ALTBC_Calculator(){
      myName = "ALTBC";
    }
    virtual ~ALTBC_Calculator(){}

    void init(ntype r_min_, ntype binw_, ntype r_max_, ntype angle_th_,
              int nTypes_, int N, ntype V, string string_out_,
              string tag_, bool debug_, bool verbose_)
    {
      int i,j, idx, ti;
      ntype shell1, shell2;
      string_out = string_out_;
      tag = tag_;
      debug = debug_;
      verbose = verbose_;
      nTypes = nTypes_; // in case of dynamic_types, this should be intended as max_nTypes
      nTypePairs = nTypes*(nTypes+1)/2;
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
      value.resize(nTypes);
      ave.resize(nTypes);
      ave2.resize(nTypes);
      for(j=0;j<nTypes;j++){
        value[j].resize(nbins2);
        ave[j].resize(nbins2);
        ave2[j].resize(nbins2);
      }

      for(i=0; i<nbins; i++){
        bins[i] = r_min + (i+0.5)*binw; // take the center of the bin for histogram
      }

      if(verbose)
      {
        for(ti=0;ti<nTypes;ti++){
          ss.str(std::string()); ss << string_out << tag << "_type"<<ti<< ".traj"; fout.open(ss.str(), ios::out);
          fout << "# 1st block: radial distance; other blocks: ALTBC 2D matrix for each frame; # deviation >= " << angle_th << " degrees\n";
          for(i=0; i<nbins; i++){
            fout << bins[i] << endl;
          }
          fout << endl; // end of block
          fout.close();
        }
      }

      // same normalization for all types ...
      ntype NVfactor = ((ntype)N) *(((ntype)(N-1))/V) * (((ntype)(N-2))/V);
      for(i=0; i<nbins; i++)
      {
        shell1 = 4.0 * M_PI * bins[i]*bins[i] * binw;
        for(j=0;j<nbins;j++)
        {
          shell2 = 4.0 * M_PI * bins[j]*bins[j] * binw;
          idx = nbins*i + j;
          norm[idx] = shell1*shell2 * NVfactor; // prepare normalization
          for(ti=0;ti<nTypes;ti++){
            ave[ti][idx]=ave2[ti][idx]=0.0; // reset averages
          }
        }
      }

      for(ti=0;ti<nTypes;ti++){
        ss.str(std::string()); ss << string_out << tag << "_type"<<ti<< ".ave"; fout.open(ss.str(), ios::out);
        fout << "# 1st block: radial distance; 2nd block: ALTBC 2D matrix; # deviation >= " << angle_th << " degrees\n";
        fout.close();

        ss.str(std::string()); ss << string_out << "_r1r2" << tag << "_type"<<ti<< ".ave"; fout.open(ss.str(), ios::out);
        fout << "# Timestep | <r_1/r_2> (with r_1<=r_2) | fluctuations\n";
        fout.close();
      }

      if(verbose)
      {
        for(ti=0;ti<nTypes;ti++){
          /*
          ss.str(std::string()); ss << string_out << "_r1r2" << tag << "_type"<<ti<< ".traj"; fout.open(ss.str(), ios::out);
          fout << "# r_1/r_2 (with r_1<=r_2) for each count; one line per frame\n";
          fout.close();
          */
          ss.str(std::string()); ss << string_out << "_peierls" << tag << "_type"<<ti<< ".traj"; fout.open(ss.str(), ios::out);
          fout << "# one line per frame: for each atom: { numCollBonds(integer) for each collBond: {rshort, ratio_shortlong, cos } }\n";
          fout.close();
        }
      }
    }

    void compute(int frameidx, int nframes, int timestep, vector<ptype> ps)
    {
      const int N=ps.size();
      int i,j,k,k0,k1, a,b, bin0, bin1, ti;
      vec rij, rik;
      ntype rijSq, rikSq, costheta, rijNorm, rikNorm, r1r2;
      vector<int> counts;
      vecflex<ntype> r1r2_ave, r1r2_ave2;
      CollinearBonds collBond;
      if(debug) cout << "*** "<<myName<<" computation for timestep " << timestep << " STARTED ***\n";

      peierlsDists.resize(N);
      counts.resize(nTypes);
      r1r2_ave.resize(nTypes);
      r1r2_ave2.resize(nTypes);
      // reset counters and values
      for(ti=0;ti<nTypes;ti++){
        for(k=0; k<nbins2; k++){ value[ti][k] = 0.0; }
        counts[ti]=0;
        r1r2_ave[ti]=r1r2_ave2[ti]=0.0;
      }

      // check angles of the type j---i---k (with i at center)
      for(i=0;i<N;i++)
      {
        peierlsDists[i].numCollBonds=0;
        ti=ps[i].label;
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
              value[ti][bin0 + nbins*bin1] += 1.0;
              value[ti][bin1 + nbins*bin0] += 1.0;
              counts[ti]+=2;
              collBond.r_short = rijNorm>rikNorm ? rikNorm : rijNorm;
              collBond.ratio_shortlong = rijNorm>rikNorm ? rikNorm/rijNorm: rijNorm/rikNorm;
              collBond.cos = costheta;
              if(peierlsDists[i].numCollBonds==peierlsDists[i].maxCollBonds) {
                if(debug){
                  cerr << "Warning: numCollBonds for atom i="<<i<<" exceeded maxCollBonds="<<peierlsDists[i].maxCollBonds<<endl;
                  cerr << "         I will ignore this bond in peierlsDists and in ALTBC.\n";
                  cerr << "                 r_short = "<<collBond.r_short<<endl;
                  cerr << "         ratio_shortlong = "<<collBond.ratio_shortlong<<endl;
                  cerr << "                     cos = "<<collBond.cos<<endl<<endl;
                }
                continue;
              }
              peierlsDists[i].collBonds[peierlsDists[i].numCollBonds] = collBond;
              peierlsDists[i].numCollBonds++;
              r1r2 = collBond.ratio_shortlong;
              r1r2_ave[ti] += 2 * r1r2; // 2x because we are increasing count twice!
              r1r2_ave2[ti] += 2 * r1r2*r1r2;
              //
            }
          }
        }
      }
      // print <r1/r2>, per type
      for(ti=0;ti<nTypes;ti++){
        ss.str(std::string()); ss << string_out << "_r1r2" << tag << "_type"<<ti<< ".ave"; fout.open(ss.str(), ios::app);
        if(counts[ti]>0)
        {
          r1r2_ave[ti] /= counts[ti];
          r1r2_ave2[ti] /= counts[ti];
          if(counts[ti]>1) r1r2_ave2[ti] = sqrt( (r1r2_ave2[ti] - r1r2_ave[ti]*r1r2_ave[ti])/(counts[ti]-1) );
        }
        fout << timestep << " " << r1r2_ave[ti] << " " << r1r2_ave2[ti] << endl;
        fout.close();
      }
      // print perierlsDists, per type
      if(verbose)
      {
        for(ti=0;ti<nTypes;ti++){
          ss.str(std::string()); ss << string_out << "_peierls" << tag << "_type"<<ti<< ".traj"; fout.open(ss.str(), ios::app);
          for(i=0;i<N;i++){
            if(ps[i].label!=ti) continue;
            fout << peierlsDists[i].numCollBonds << " ";
            for(j=0;j<peierlsDists[i].numCollBonds; j++) {
                fout << peierlsDists[i].collBonds[j].r_short << " ";
                fout << peierlsDists[i].collBonds[j].ratio_shortlong << " ";
                fout << peierlsDists[i].collBonds[j].cos << " ";
            }
          }
          fout << endl;
          fout.close();
        }
      }
      // Add current histogram to the average...
      for(ti=0;ti<nTypes;ti++)
      {
        for(k=0; k<nbins2; k++)
        {
          if(counts[ti]>0) value[ti][k] /= (counts[ti]*binw*binw); // ci sta l'area del bin??
          value[ti][k] /= norm[k];
          ave[ti][k] += value[ti][k];
          ave2[ti][k] += value[ti][k]*value[ti][k];
        }
        // ... and print it, if verbose
        if(verbose)
        {
          ss.str(std::string()); ss << string_out << tag << "_type"<<ti<< ".traj"; fout.open(ss.str(), ios::app);
          fout << endl; // start new block
          for(k0=0; k0<nbins; k0++){ // print as 2D matrix
            for(k1=0; k1<nbins; k1++){
              k = k0 + nbins*k1;
              fout << value[ti][k] << " ";
            }
            fout << endl;
          }
          fout << endl; // end of block
          fout.close();
        }
      }

      if(frameidx == (nframes-1)) // if last frame: print final results
      {
        for(ti=0;ti<nTypes;ti++)
        {
          ss.str(std::string()); ss << string_out << tag << "_type"<<ti<< ".ave"; fout.open(ss.str(), ios::app);
          for(k0=0; k0<nbins; k0++){
            fout << bins[k0] << endl; // bins
          }
          fout << endl; // end of block
          for(k0=0; k0<nbins; k0++){
            for(k1=0; k1<nbins; k1++){
              k = k0 + nbins*k1;
              ave[ti][k] /= nframes;
              ave2[ti][k] /= nframes;
              fout << ave[ti][k] << " "; // average
            }
            fout << endl;
          }
          /* // do NOT print errors... who does even look at them?
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
        }
        if(debug) cout << "Average "<<myName<<" printed to file\n";
      }
      if(debug) cout << "*** "<<myName<<" computation for timestep " << timestep << " ENDED ***\n\n";
      return;
    }
};
#endif
