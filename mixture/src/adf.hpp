#ifndef _ADF_H_
#define _ADF_H_
using namespace std;
//--------------------- Angular Distribution Function (distribution of cos(theta_ijk) ----------------------------------//
// must be preceeded by computation of nearest neighbours

template <class ntype, class ptype>
class ADF_Calculator
{
  using vec=myvec<ntype,3>;
  using mat=mymatrix<ntype,3,3>;
  private:
    ntype binw;
    int nbins, nTypes, nTypePairs;
    vecflex<ntype> bins;
    vecflex< vecflex< vecflex<ntype> >> value, ave, ave2;
    string string_out, myName, tag;
    fstream fout;
    stringstream ss;
    bool debug, verbose;

    int types2int(int ti, int tj){ // map type pairs (ti,tj) in 0,1,...,nTypes-1 to integer index 0,1,...,nTypePairs
      if (ti>tj) return types2int(tj,ti); // map to ti<=tj
      if(ti==tj) return ti*nTypes - int(ti*(ti-1)/2);
      else return 1 + types2int(ti,tj-1);
    }

  public:
    ADF_Calculator(){
      myName = "ADF";
    }
    virtual ~ADF_Calculator(){}

    void init(ntype binw_, int nTypes_, string string_out_, string tag_, bool debug_, bool verbose_)
    {
      nTypes = nTypes_;
      nTypePairs = nTypes*(nTypes+1)/2;
      string_out = string_out_;
      tag = tag_;
      debug = debug_;
      verbose = verbose_;
      binw = binw_;
      nbins = int(floor( 2.0 / binw )); // -1 < cos(angle) < 1
      bins.resize(nbins);
      int i,t1,t2,tp;
      value.resize(nTypes);
      ave.resize(nTypes);
      ave2.resize(nTypes);
      if(verbose)
      {
        ss.str(std::string()); ss << string_out << tag << ".traj"; fout.open(ss.str(), ios::out);
        fout << "# First block: cos(angle); other blocks: ADF for each frame )\n";
      }
      for(i=0; i<nbins; i++)
      {
        bins[i] = -1.0 + (i+0.5)*binw; // take the center of the bin for histogram
        if(verbose) fout << bins[i] << endl;
      }
      if(verbose) fout.close();

      for(t1=0;t1<nTypes;t1++)
      {
        value[t1].resize(nTypePairs);
        ave[t1].resize(nTypePairs);
        ave2[t1].resize(nTypePairs);
        for(tp=0;tp<nTypePairs;tp++)
        {
          //if(debug) cout << " Inizializzo la coppia (" << t1 <<","<<t2<<")-->"<<  a << "\n";
          value[t1][tp].resize(nbins);
          ave[t1][tp].resize(nbins);
          ave2[t1][tp].resize(nbins);
          for(i=0; i<nbins; i++)
          {
            ave[t1][tp][i] = ave2[t1][tp][i] = 0.0;
          }
        }
      }

      ss.str(std::string()); ss << string_out << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "# cos(angle) | ADF's: 0-0-0 0-0-1 0-0-2 ... 1-0-1 1-0-2 ... 0-1-0 0-1-1 ... | ADF errors.\n";
      fout.close();
    }

    void compute(int frameidx, int nframes, int timestep, vector<ptype> ps)
    {
      const int N=ps.size();
      int i,j,k, a,b, bin, counts=0, ti,tj,tk,tp;
      vec rij, rik;
      ntype rijSq, rikSq, costheta;
      for(ti=0;ti<nTypes;ti++){
        for(tp=0;tp<nTypePairs;tp++){
          for(bin=0; bin<nbins; bin++){
            //if(debug) cout << ti << " " << tp << " " << bin << endl;
            value[ti][tp][bin] = 0.0;
          }
        }
      }
      if(debug) cout << "*** "<<myName<<" computation for timestep " << timestep << " STARTED ***\n";
      for(i=0;i<N;i++){
        ti = ps[i].label;
        for(a=1;a<ps[i].neigh_list[0].size();a++){
            j = ps[i].neigh_list[0][a];
            //if(j>i) continue; // avoid double counting!
            tj = ps[j].label;
            rij = ps[i].rij_list[0][a];
            rijSq = ps[i].rijSq_list[0][a];
            for(b=0;b<a;b++){
              k = ps[i].neigh_list[0][b];
              tk = ps[k].label;
              rik = ps[i].rij_list[0][b];
              rikSq = ps[i].rijSq_list[0][b];
              costheta = (rij*rik) / sqrt(rijSq*rikSq);
              bin = int(floor( (costheta+1.0)/binw)); // min value is -1.0 !
              if(bin<nbins){
                tp = types2int(tj,tk);
                value[ti][tp][bin] += 1.0; //  j--i--k angle
                counts++;
              }
            }
        }
      }
      if(verbose)
      {
        ss.str(std::string()); ss << string_out << tag << ".traj"; fout.open(ss.str(), ios::app);
        fout << endl; // start new block
      }
      for(k=0; k<nbins; k++){
        for(ti=0;ti<nTypes;ti++){
          for(tp=0;tp<nTypePairs;tp++){
            if(counts>0) value[ti][tp][k] /= (counts*binw); // ci sta la larghezza del bin?
            ave[ti][tp][k] += value[ti][tp][k];
            ave2[ti][tp][k] += value[ti][tp][k]*value[ti][tp][k];
            if(verbose)
            {
              fout << value[ti][tp][k];
              if(ti==nTypes-1 && tp==nTypePairs-1) fout << endl;
              else fout << " ";
            }
          }
        }
      }
      if(verbose) fout.close();

      if(frameidx == (nframes-1))
      {
        ss.str(std::string()); ss << string_out << tag << ".ave"; fout.open(ss.str(), ios::app);
        for(k=0; k<nbins; k++)
        {
          fout << bins[k] << " ";
          for(ti=0;ti<nTypes;ti++){
            for(tp=0;tp<nTypePairs;tp++){
              ave[ti][tp][k] /= nframes;
              fout << ave[ti][tp][k] << " ";
            }
          }
          for(ti=0;ti<nTypes;ti++){
            for(tp=0;tp<nTypePairs;tp++){
              ave2[ti][tp][k] /= nframes;
              if(nframes>1) fout << sqrt( (ave2[ti][tp][k]-ave[ti][tp][k]*ave[ti][tp][k])/(nframes -1) );
              else fout << 0.0;
              if(ti==nTypes-1 && tp==nTypePairs-1) fout << endl;
              else fout << " ";
            }
          }
        }
        fout.close();
        if(debug) cout << "Average "<<myName<<" printed to file\n";
      }
      if(debug) cout << "*** "<<myName<<" computation for timestep " << timestep << " ENDED ***\n\n";
      return;
    }
};

#endif
