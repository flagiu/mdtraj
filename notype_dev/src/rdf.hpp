#ifndef _RDF_H_
#define _RDF_H_
using namespace std;
//------------------------------- Radial Distribution Function ----------------------------------//

template <class ntype, class ptype>
class RDF_Calculator
{
  using vec=myvec<ntype,3>;
  using mat=mymatrix<ntype,3,3>;
  private:
    ntype binw;
    int nbins;
    vecflex<ntype> bins, norm, value, ave, ave2;
    string string_out, myName, tag;
    fstream fout;
    stringstream ss;
  public:
    RDF_Calculator(){
      myName = "RDF";
    }
    virtual ~RDF_Calculator(){}

    void init(ntype binw_, vec box_diagonal, int N, ntype V, string string_out_, string tag_)
    {
      string_out = string_out_;
      tag = tag_;
      binw = binw_;
      nbins = int(floor( 0.5*box_diagonal.abs().min() / binw ));
      bins.resize(nbins);
      norm.resize(nbins);
      value.resize(nbins);
      ave.resize(nbins);
      ave2.resize(nbins);
      ntype r, shell1, shell2, normalization = (N-1)/V * 4.0*M_PI/3.0 * N/2.0;
      // normalization:
      //    N*(N-1)/2 / V = average density of pairs i<j
      //    4*pi/3*((r_2)^3 - (r_1)^3) = volume of the shell
      shell1 = 0.0;
      ss.str(std::string()); ss << string_out << tag << ".traj"; fout.open(ss.str(), ios::out);
      fout << "# First block: radial distance; other blocks: g(r)\n";
      for( auto i=0; i<nbins; i++){
        bins[i] = (i+0.5)*binw; // take the center of the bin for histogram
        fout << bins[i] << endl;
        r = (i+1)*binw; // take the upper side of the bin for shells
        shell2 = r*r*r;
        norm[i] = normalization*(shell2 - shell1); // normalize to ideal gas radial distribution
        shell1 = shell2;
        ave[i] = 0.0;
        ave2[i] = 0.0;
      }
      fout.close();
      ss.str(std::string()); ss << string_out << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "# Radial distance, g(r), g(r) error.\n";
      fout.close();
    }

    void compute(int frameidx, int nframes, int timestep, vector<ptype> ps, mat box, mat boxInv, bool debug)
    {
      int i,j, k;
      const int N=ps.size();
      vec rij;
      ntype r,r_mic;
      for(k=0; k<nbins; k++) value[k] = 0.0;
      if(debug) cout << "*** "<<myName<<" computation for timestep " << timestep << " STARTED ***\n";
      for(i=0;i<N;i++){
        for(j=i+1;j<N;j++){ // i<j
          rij = ps[j].r - ps[i].r; // real distance
          r = rij.norm();
          r_mic = (rij - box*round(boxInv*rij)).norm(); // first periodic image
          if( r_mic < r) r = r_mic; // if closer, choose first periodic image
          // if(debug) cout << "  PBC for (i,j)=(" <<i<<","<<j<<") : " << r_true << " --> " << r << endl;
          k = int(floor( r/binw ));
          if(k<nbins){
            value[k] += 1.0;
          }
        }
      }
      ss.str(std::string()); ss << string_out << tag << ".traj"; fout.open(ss.str(), ios::app);
      fout << endl; // start new block
      for(k=0; k<nbins; k++){
        value[k] /= norm[k];
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
