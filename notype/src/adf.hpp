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
    int nbins;
    vecflex<ntype> bins, value, ave, ave2;
    string string_out, myName, tag;
    fstream fout;
    stringstream ss;
  public:
    ADF_Calculator(){
      myName = "ADF";
    }
    virtual ~ADF_Calculator(){}

    void init(ntype binw_, ntype rcut1, string string_out_, string tag_)
    {
      string_out = string_out_;
      tag = tag_;
      binw = binw_;
      nbins = int(floor( 2.0 / binw )); // -1 < cos(angle) < 1
      bins.resize(nbins);
      value.resize(nbins);
      ave.resize(nbins);
      ave2.resize(nbins);
      ss.str(std::string()); ss << string_out << tag << ".traj"; fout.open(ss.str(), ios::out);
      fout << "# First block: cos(angle); other blocks: ADF for each frame. # cutoff = "<<rcut1<<endl;
      for( auto i=0; i<nbins; i++){
        bins[i] = -1.0 + (i+0.5)*binw; // take the center of the bin for histogram
        fout << bins[i] << endl;
        ave[i] = ave2[i] = 0.0;
      }
      fout.close();
      ss.str(std::string()); ss << string_out << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "# cos(angle), ADF, ADF error. # cutoff = "<<rcut1<<endl;
      fout.close();
    }

    void compute(int frameidx, int nframes, int timestep, vector<ptype> ps, bool debug)
    {
      const int N=ps.size();
      int i,j,k, a,b, bin, counts=0;
      vec rij, rik;
      ntype rijSq, rikSq, costheta;
      for(a=0; a<nbins; a++) value[a] = 0.0;
      if(debug) cout << "*** "<<myName<<" computation for timestep " << timestep << " STARTED ***\n";
      for(i=0;i<N;i++){
        if(debug) cout << "ADF: i="<<i<<endl;
        for(a=0;a<ps[i].neigh_list[0].size();a++){
            j = ps[i].neigh_list[0][a];
            if(debug) cout << "ADF:  a="<<a<<" --> j="<<j<<endl;
            //if(j>i) continue; // avoid double counting! NONSENSE!!!!!!
            rij = ps[i].rij_list[0][a];
            rijSq = ps[i].rijSq_list[0][a];
            for(b=a+1;b<ps[i].neigh_list[0].size();b++){
              k = ps[i].neigh_list[0][b]; //useless?
              if(debug) cout << "ADF:   b="<<b<<" --> k="<<k<<endl;
              rik = ps[i].rij_list[0][b];
              rikSq = ps[i].rijSq_list[0][b];
              costheta = (rij*rik) / sqrt(rijSq*rikSq);
              bin = int(floor( (costheta+1.0)/binw)); // min value is -1.0 !
              if(debug) cout << "ADF:     (j--i--k)=("<<j<<"--"<<i<<"--"<<k<<") --> cos(jik)="<<costheta<<" --> bin="<<bin<<" out of nbins="<<nbins<<endl;
              if(bin<nbins){
                value[bin] += 1.0;
                counts++;
              }
            }
        }
      }
      ss.str(std::string()); ss << string_out << tag << ".traj"; fout.open(ss.str(), ios::app);
      fout << endl; // start new block
      for(k=0; k<nbins; k++){
        value[k] /= (counts*binw); // ci sta la larghezza del bin?
        fout << value[k] << endl;
        ave[k] += value[k];
        ave2[k] += value[k]*value[k];
      }
      fout.close();

      if(frameidx == (nframes-1)){
          ss.str(std::string()); ss << string_out << tag << ".ave"; fout.open(ss.str(), ios::app);
          for(k=0; k<nbins; k++){
            ave[k] /= nframes;
            ave2[k] /= nframes;
            fout << bins[k] << " " << ave[k] << " " << sqrt( (ave2[k]-ave[k]*ave[k])/(nframes -1) ) << endl;
          }
          fout.close();
          if(debug) cout << "Average "<<myName<<" printed to file\n";
      }
      if(debug) cout << "*** "<<myName<<" computation for timestep " << timestep << " ENDED ***\n\n";
      return;
    }
};

#endif
