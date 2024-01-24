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
    int image_convention=1; // 0: no images (cluster) ; 1: 1 replica (minimum image) ; -1: all replicas (crystal)
    ntype binw;
    int nbins, nTypes, nTypePairs;
    vecflex<ntype> bins;
    vector< vecflex<ntype> > norm, value, ave, ave2;
    string string_out, myName, tag;
    fstream fout;\
    stringstream ss;
    bool debug, verbose;

    int types2int(int ti, int tj){ // map type pairs (ti,tj) in 0,1,...,nTypes-1 to integer index 0,1,...,nTypePairs
      if (ti>tj) return types2int(tj,ti); // map to ti<=tj
      if(ti==tj) return ti*nTypes - int(ti*(ti-1)/2);
      else return 1 + types2int(ti,tj-1);
    }
    void int2types(int t, int *t1, int *t2){
      int x,low,high;
      for(x=0; x<nTypes; x++){
        low = x * nTypes - int(x*(x-1)/2);
        high = (x+1) * nTypes - int((x+1)*x/2) - 1;
        if(t>=low && t<=high){
          *t1 = x;
          *t2 = x + (t-low);
          break;
        }
      }
      return;
    }

  public:
    ntype rmax;

    RDF_Calculator(){
      myName = "RDF";
    }
    virtual ~RDF_Calculator(){}

    void init(ntype binw_, ntype rmax_, mat box, int N, ntype V, int nTypes_, int* Nt_, string string_out_, string tag_, bool debug_, bool verbose_)
    {
      string_out = string_out_;
      tag = tag_;
      debug =debug_;
      verbose = verbose_;
      binw = binw_;
      nTypes = nTypes_;
      nTypePairs = nTypes*(nTypes+1)/2;
      int k,t1,t2,tp;
      ntype r, r_limit, shell1, shell2, normalization;
      // r_limit = 1/2 * min( |ax+bx+cx|, |ay+by+cy|, |az+bz+cz| ) for minimum image convention
      r_limit = 0.5* min( min( fabs(box[0][0]+box[0][1]+box[0][2]), fabs(box[1][0]+box[1][1]+box[1][2]) ), fabs(box[2][0]+box[2][1]+box[2][2]) );
      if(rmax_<=0) rmax = r_limit;
      else
      {
        rmax = rmax_;
        if(rmax<binw) { cout << "[ERROR: bin_width > r_max for g(r) ! ]\n"; exit(1);}
        if(rmax>r_limit) {
          cout << "WARNING: r_max="<<rmax<<" with r_limit="<<r_limit<<" may result in artefacts.\n";
        }
      }
      nbins = int(floor( rmax / binw ));
      bins.resize(nbins);

      shell1 = 0.0;
      if(verbose) {
          ss.str(std::string()); ss << string_out << tag << ".traj"; fout.open(ss.str(), ios::out);
          fout << "# First block: r; Other blocks: g_00(r), g_01(r), ... for each frame.\n";
      }
      for( k=0; k<nbins; k++){
        bins[k] = (k+0.5)*binw; // take the center of the bin for histogram
        if(verbose) fout << bins[k] << endl;
      }
      if(verbose) fout.close();

      norm.resize(nTypePairs);
      value.resize(nTypePairs);
      ave.resize(nTypePairs);
      ave2.resize(nTypePairs);
      for(t1=0;t1<nTypes;t1++)
      {
        for(t2=t1;t2<nTypes;t2++)
        {
          tp = types2int(t1,t2);
          //if(debug) cout << " Inizializzo la coppia (" << t1 <<","<<t2<<")-->"<<  a << "\n";
          norm[tp].resize(nbins);
          value[tp].resize(nbins);
          ave[tp].resize(nbins);
          ave2[tp].resize(nbins);
          //if(t1!=t2) normalization =     Nt_[t1]/V * 4.0*M_PI/3.0 * Nt_[t2]/2.0; // * 2.0; // multiply by 2 because it's off-diagonal
          //else       normalization = (Nt_[t1]-1)/V * 4.0*M_PI/3.0 * Nt_[t2]/2.0;
          // normalization:
          //    N[t1]*N[t2]/2 / V = average density of ordered pairs i(t1) < j(t2)
          //    N[t1]*(N[t1]-1)/2 /V = average density of ordered pairs i(t1) < j(t1)
          //    4*pi/3*((r_2)^3 - (r_1)^3) = volume of the shell
          //    the off-diagonal part must be further divided by 2 because of the symmetry between t1 and t2
          if(t1!=t2) normalization =       Nt_[t1]*Nt_[t2]/2 /V * 4.0*M_PI/3.0 *2.0;
          else  {
            if(Nt_[t1]>1) normalization = Nt_[t1]*(Nt_[t1]-1)/2 /V * 4.0*M_PI/3.0;
            else          normalization = 2*Nt_[t1]*(2*Nt_[t1]-1)/2 /(2*V) * 4.0*M_PI/3.0;
          }
          shell1 = 0.0;
          for(k=0; k<nbins; k++){
            bins[k] = (k+0.5)*binw; // take the center of the bin for histogram
            r = (k+1)*binw; // take the upper side of the bin for shells
            shell2 = r*r*r;
            norm[tp][k] = normalization*(shell2 - shell1); // normalize to ideal gas radial distribution
            shell1 = shell2;
            ave[tp][k] = 0.0;
            ave2[tp][k] = 0.0;
          }
        }
      }

      ss.str(std::string()); ss << string_out << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "# r | g_00(r), g_01(r), ... | error for each g(r).\n";
      fout.close();
    }

    void compute(int frameidx, int nframes, int timestep, vector<ptype> ps, PBC<ntype> *pbc)
    {
      int i,j, k, tp, ti,tj, num_images,img;
      const int N=ps.size();
      vec rij, r_images[MAX_N_IMAGES];
      ntype r,r_mic;
      for(tp=0;tp<nTypePairs;tp++){
        for(k=0; k<nbins; k++){
          value[tp][k] = 0.0;
        }
      }
      if(debug) cout << "*** "<<myName<<" computation for timestep " << timestep << " STARTED ***\n";
      for(i=0;i<N;i++){
        for(j=i;j<N;j++){
          //if(i==j) continue;
          rij = ps[j].r - ps[i].r; // real distance
          num_images = pbc->apply(rij, r_images);
          /*
          r = rij.norm();
          r_mic = (rij - box*round(boxInv*rij)).norm(); // first periodic image
          if( r_mic < r) r = r_mic; // if closer, choose first periodic image
          */
          // if(debug) cout << "  PBC for (i,j)=(" <<i<<","<<j<<") : " << r_true << " --> " << r << endl;
          for(img=0;img<num_images;img++){
            r = r_images[img].norm();
            if(r==0.0) continue;
            k = int(floor( r/binw ));
            if(k<nbins){
              ti = ps[i].label;
              tj = ps[j].label;
              tp = types2int(ti,tj);
              value[tp][k] += 1.0;
            }
          }
        }
      }

      if(verbose)
      {
        ss.str(std::string()); ss << string_out << tag << ".traj"; fout.open(ss.str(), ios::app);
        fout << endl; // start new block
      }
      for(k=0; k<nbins; k++)
      {
        for(tp=0;tp<nTypePairs;tp++)
        {
          value[tp][k] /= norm[tp][k];
          ave[tp][k] += value[tp][k];
          ave2[tp][k] += value[tp][k]*value[tp][k];
          if(verbose)
          {
            fout << value[tp][k];
            if(tp<nTypePairs-1) fout << " ";
            else fout << endl;
          }
        }
      }
      if(verbose) fout.close();

      if(frameidx == (nframes-1)){
          ss.str(std::string()); ss << string_out << tag << ".ave"; fout.open(ss.str(), ios::app);
          for(k=0; k<nbins; k++)
          {
            fout << bins[k] << " "; // print bin
            for(tp=0;tp<nTypePairs;tp++)
            {
              ave[tp][k] /= nframes;
              fout << ave[tp][k] << " ";
            }
            for(tp=0;tp<nTypePairs;tp++)
            {
              ave2[tp][k] /= nframes;
              if(nframes>1) fout << sqrt( (ave2[tp][k]-ave[tp][k]*ave[tp][k])/(nframes -1) );
              else  fout << 0.0;
              if(tp<nTypePairs-1) fout << " ";
              else fout << endl;
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
