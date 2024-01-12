#ifndef _NEIGH_AND_BOND_LIST_H_
#define _NEIGH_AND_BOND_LIST_H_
#include<algorithm>
using namespace std;

//---- "init" functions: prepare arrays, averages and headers for output files
//---- "compute" functions: compute quantities for each frame and total average

//------------------------------ Neighbour and Bond list --------------------------------------------//

#define MAX_NSPHERE 3

template <class ntype, class ptype>
class Neigh_and_Bond_list
{
  using vec=myvec<ntype,3>;
  using mat=mymatrix<ntype,3,3>;
  private:
    int p1half, p2half, p1,p2;
    string string_cn_out, string_rmin_out, string_rmax_out, log_file, myName, tag;
    fstream fout;
    stringstream ss;
    bool debug, verbose;

    ntype fcut_generic(ntype x, int pow1, int pow2)
    {
      if(x==1.0) return ((ntype)pow1) / ((ntype)pow2); // avoid divergence
      ntype x1=1.0, x2=1.0;
      for(auto i=0; i<pow2; i++) {
        x2 *= x;
        if(i==pow1-1) x1=x2;
      }
      return (1.0-x1) / (1.0-x2); // x1 = x^p1, x2 = x^p2
    }

  public:
    ntype rmaxSq;
    int Nsphere, N, nTypes, nTypePairs;
    vecflex<ntype> rcut[MAX_NSPHERE], rcutSq[MAX_NSPHERE];
    vector< vecflex<ntype> > neigh[MAX_NSPHERE]; // Nspheres X nTypePairs X N
    vecflex<ntype> neigh_anytype[MAX_NSPHERE]; // Nspheres X N (agnostic of types)
    vector<int> bond_list[MAX_NSPHERE];  // stores all bonds (i,j), encoded into an integer through ij2int()

    Neigh_and_Bond_list(){
      myName = "NEIGH & BOND List";
    }
    virtual ~Neigh_and_Bond_list(){}

    int ij2int(int i, int j, int N){
      return (i<j ? N*i+j : N*j+i); // i<j = 0,...,N-1
    }
    int int2i(int x, int N){ return x/N; }
    int int2j(int x, int N){ return x%N; }
    int get_bond_index(int i, int j, int N, int shell) {
      int bond_encoding = ij2int(i,j,N);
      int bond_idx = indexOf<int>(&(bond_list[shell]), bond_encoding); // indexOf is defined in "utility.hpp"
      if(bond_idx==-1) {
        cout << "ERROR: bond not found with (i,j)=("<<i<<","<<j<<")\n";
        exit(1);
      }
      return bond_idx;
    }
    ntype fcut(ntype xSq) { return fcut_generic(xSq,p1half,p2half); }
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

    void read_rcut_from_file(string s_rcut, vecflex<ntype> cutoffDefault[MAX_NSPHERE])
    {
      int u,t;
      string line, x;
      if(s_rcut=="__NOT_DEFINED__") {
        cout << "*** Using default RCUT for all type pairs ***\n";
        for(u=0;u<Nsphere;u++) {
          for(t=0;t<nTypePairs;t++){
            rcut[u][t]=cutoffDefault[u][0];
            rcutSq[u][t]=SQUARE(rcut[u][t]);
          }
        }
      }
      else {
        if(debug) cout << "*** Reading RCUT from file "<<s_rcut<<" ***\n";
        fout.open(s_rcut, ios::in);
        if(fout.fail()) // checks to see if file opended
        {
          cout << "ERROR: could not open RCUT file '"<<s_rcut<<"'\n\n";
          exit(1);
        }
        for(u=0;u<Nsphere;u++)
        {
          getline(fout,line);
          ss << line;
          for(t=0;t<nTypePairs;t++)
          {
            ss >> x;
            rcut[u][t] = stof(x);
            rcutSq[u][t] = SQUARE(rcut[u][t]);
          }
          ss.str(std::string()); ss.clear(); // clear the string stream!
        }
        fout.close();
        if(debug) cout << "*** Reading RCUT file DONE ***\n";
      }
      cout << "#-------- RCUT summary --------#\n";
      for(u=0;u<Nsphere;u++)
      {
        cout << "SPHERE u="<<u<<" : ";
        rcut[u].show();
      }
      cout << "#------------------------------#\n";
    }

    void print_bond_summary(vector<ptype> ps)
    {
      int i,u,ii;
      ss.str(std::string()); ss << log_file << tag; fout.open(ss.str(), ios::app);
      fout << "#---------------------- BOND SUMMARY --------------------------#\n";
      for(u=0;u<Nsphere;u++)
      {
        for(i=0;i<N;i++)
        {
          fout << "  Sphere u="<<u<<" of particle i="<<i<<" of type="<<ps[i].label<<" contains Nc="<<ps[i].neigh_list[u].size()<<" neighbours:\n   ";
          for(ii=0;ii<ps[i].neigh_list[u].size();ii++) fout<<ps[i].neigh_list[u][ii]<<" ";
          fout << "\n   ";
          for(ii=0;ii<ps[i].neigh_list[u].size();ii++) fout<<sqrt(ps[i].rijSq_list[u][ii])<<" ";
          fout << endl;
        }
        fout << "> System has "<<bond_list[u].size()<<" bonds within rcut=";
        for(int t=0;t<nTypePairs;t++) fout <<rcut[u][t]<< " ";
        fout << endl;
      }
      fout << "#--------------------------------------------------------------#\n";
      fout.close();
      //cout << " Summary of neighbour and bond lists saved into log file: "<<log_file <<tag <<endl;
    }

    void sort_by_distance_3vec(vector<ntype>* dist, vector<int>* a, vector<vec>* b)
    {
      // assuming same length for all inputs!
      vector<size_t> p(dist->size()); // vector of 0,1,...,size-1
      for(auto i=0;i<dist->size();i++) p[i]=i;
      // sort p according to rSq_list
      std::sort(p.begin(), p.end(), [&](std::size_t i, std::size_t j){ return (*dist)[i] < (*dist)[j]; } );
      // now p is a permutation vector!
      vector<ntype> dist_copy = (*dist);
      vector<int> a_copy = (*a);
      vector<vec> b_copy = (*b);
      for(auto i=0;i<dist->size();i++)
      {
        (*dist)[i] = dist_copy[p[i]];
        (*a)[i] = a_copy[p[i]];
        (*b)[i] = b_copy[p[i]];
      }
      return;
    }
    void sort_by_distance_2vec(vector<ntype>* dist, vector<int>* a)
    {
      // assuming same length for all inputs!
      vector<size_t> p(dist->size()); // vector of 0,1,...,size-1
      for(auto i=0;i<dist->size();i++) p[i]=i;
      // sort p according to rSq_list
      std::sort(p.begin(), p.end(), [&](std::size_t i, std::size_t j){ return (*dist)[i] < (*dist)[j]; } );
      // now p is a permutation vector!
      vector<ntype> dist_copy = (*dist);
      vector<int> a_copy = (*a);
      for(auto i=0;i<dist->size();i++)
      {
        (*dist)[i] = dist_copy[p[i]];
        (*a)[i] = a_copy[p[i]];
      }
      return;
    }

    void init(string s_rcut, vecflex<ntype> cutoffDefault[MAX_NSPHERE], int Nsphere_,
      int p1half_, int N_, int nTypes_, string log_file_, string tag_, bool debug_, bool verbose_)
    {
      log_file = log_file_;
      tag=tag_;
      debug=debug_;
      verbose=verbose_;
      Nsphere=Nsphere_;
      N=N_;
      nTypes=nTypes_;
      nTypePairs=nTypes*(nTypes+1)/2;
      for(int u=0;u<Nsphere;u++)
      {
        rcut[u].resize(nTypePairs);
        rcutSq[u].resize(nTypePairs);
      }
      read_rcut_from_file(s_rcut, cutoffDefault);
      p1half=p1half_;
      p2half = 2*p1half;
      p1 = 2*p1half;
      p2 = 2*p2half;
      for(int u=0;u<Nsphere;u++){
        if(neigh[u].size()!=nTypePairs) neigh[u].resize(nTypePairs);
        for(int j=0;j<nTypePairs;j++){
          if(neigh[u][j].length()!=N) neigh[u][j].resize(N);
        }
        if(neigh_anytype[u].length()!=N) neigh_anytype[u].resize(N);
      }
    }

    void build(int timestep, vector<ptype>& ps, mat box, mat boxInv)
    { // NOTA BENE: ps deve essere passato con &, perche' dobbiamo modificare le sue variabili
      int u,i,j,t;
      vec rij, rij_mic;
      ntype rijSq, rijSq_mic;
      vector<ntype> bond_list_rijSq[MAX_NSPHERE];
      if(debug) cout << "\n*** "<<myName<<" computation STARTED ***\n";

      //---- Reset counters and lists ----//
      for(u=0;u<Nsphere;u++)
      {
        bond_list[u].clear();
        bond_list_rijSq[u].clear();
      }

      for(i=0;i<N;i++)
      {
        for(u=0;u<Nsphere;u++)
        {
          ps[i].neigh_list[u].clear();
          ps[i].rij_list[u].clear();
          ps[i].rijSq_list[u].clear();
          for(t=0;t<nTypePairs;t++) neigh[u][t][i]=0.;
        }
      }
      if(debug) cout << " * Reset counters and lists DONE\n";

      //---- Build neighbour list (and save rij vectors) and compute rmax ----//
      rmaxSq=0.0;
      for(i=0;i<N;i++)
      {
        for(j=i+1;j<N;j++)
        {
          t = types2int( ps[i].label, ps[j].label);
          rij = ps[j].r - ps[i].r;
          rijSq = rij.sq();
          rij_mic = rij - box*round(boxInv*rij); // first periodic image
          rijSq_mic = rij_mic.sq();
          if(rijSq_mic < rijSq){ // if closer, choose first periodic image
            rijSq = rijSq_mic;
            rij = rij_mic;
            //if(rijSq>rmaxSq) rmaxSq=rijSq;
          }
          //else if(rijSq_mic>rmaxSq) rmaxSq=rijSq_mic; // update rmaxSq with max(rijSq,rijSq_mic)
          if(rijSq>rmaxSq) rmaxSq=rijSq;

          for(u=0;u<Nsphere;u++)
          {
            if(rijSq <= rcutSq[u][t])
            {
              bond_list[u].push_back( ij2int(i,j,N) );
              bond_list_rijSq[u].push_back(rijSq); // this is just for ordering bond_list later

              ps[i].neigh_list[u].push_back(j);
              ps[j].neigh_list[u].push_back(i);

              ps[i].rij_list[u].push_back(   rij);
              ps[j].rij_list[u].push_back(-1*rij);

              ps[i].rijSq_list[u].push_back(rijSq);
              ps[j].rijSq_list[u].push_back(rijSq);
            }
          }
        }
      }
      // sort neigh_list and rij_list according to rijSq_list, for each particle
      for(i=0;i<N;i++)
      {
        for(u=0;u<Nsphere;u++)
        {
          sort_by_distance_3vec( &ps[i].rijSq_list[u], &ps[i].neigh_list[u], &ps[i].rij_list[u] );
        }
      }
      // sort bond_list according to bond_list_rijSq
      for(u=0;u<Nsphere;u++) sort_by_distance_2vec( &bond_list_rijSq[u], &bond_list[u]);

      if(debug || verbose) print_bond_summary(ps); // to the log file
      if(debug)
      {
        cout << " * Build neighbour and bond lists DONE\n";
        cout << "*** "<<myName<<" computation for timestep " << timestep << " ENDED ***\n\n";
      }
      return;
    }

    //----------------- Coordination Number (done on 1st sphere) ------------//
    void init_coordnum(string string_cn_out_)
    {
      if(debug) cout<<"*** Initializing COORDNUM within "<<myName<<"***\n";
      string_cn_out = string_cn_out_;
      if(verbose)
      {
        ss.str(std::string()); ss << string_cn_out << tag << ".dat"; fout.open(ss.str(), ios::out);
        fout << "#Timestep | Particle idx | Coordination number for each type pair: 00 | 01 | 02 ... . # cutoffs = ";
        for(int t=0;t<nTypePairs;t++) fout<<rcut[0][t]<<" ";
        fout<<endl;
        fout.close();
      }

      ss.str(std::string()); ss << string_cn_out << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "#Timestep | <coordination number> for each type pair | Fluctuations for each. # cutoffs = ";
      for(int t=0;t<nTypePairs;t++) fout<<rcut[0][t]<<" ";
      fout<<endl;
      fout.close();

      const int u=0;
      if(neigh[u].size()!=nTypePairs) neigh[u].resize(nTypePairs);
      for(int j=0;j<nTypePairs;j++){
        if(neigh[u][j].length()!=N) neigh[u][j].resize(N);
      }
      //for(int i=0;i<N;i++) neigh[0][i]=0.0 ; // start counter
      if(debug) cout<<"*** Initialization completed ***\n";
    }

    void compute_coordnum(int timestep, vector<ptype> ps)
    {
      if(debug) cout << "*** COORDNUM computation for timestep " << timestep << " STARTED ***\n";
      int i, j, k, t, u;
      ntype fval, rijSq;
      vec rij;

      for(i=0;i<N;i++){
        for(t=0;t<nTypePairs;t++) neigh[0][t][i] = 0.; // start counters
        neigh_anytype[0][i] = 0.; // start counters
      }

      u=0; // 1st sphere neighbours
      for(i=0;i<N;i++){
        for(k=0;k<ps[i].neigh_list[u].size();k++){ // search in neigh list
          j = ps[i].neigh_list[u][k];
          if(j>i) continue; // avoid double counting!
          rij = ps[i].rij_list[u][k];
          rijSq = ps[i].rijSq_list[u][k];
          // fval = fcut( rijSq/rcutSq[0][t] );         // smooth
          fval = ( rijSq <= rcutSq[u][t] ? 1.0 : 0.0 ); // sharp
          t = types2int( ps[i].label, ps[j].label);
          neigh[u][t][i] += fval;
          neigh_anytype[u][i] += fval;
          t = types2int(ps[j].label, ps[i].label); // take into account i<j
          neigh[u][t][j] += fval;
          neigh_anytype[u][j] += fval;
        }
      }

      if(verbose)
      {
        ss.str(std::string()); ss << string_cn_out << tag << ".dat"; fout.open(ss.str(), ios::app);
        for(i=0;i<N;i++)
        {
          fout << timestep << " " << i;
          for(t=0;t<nTypePairs;t++) fout << " " << neigh[0][t][i];
          fout << endl;
        }
        fout.close();
      }

      ss.str(std::string()); ss << string_cn_out << tag << ".ave"; fout.open(ss.str(), ios::app);
      fout << timestep;
      for(t=0;t<nTypePairs;t++) fout << " " << neigh[u][t].mean();
      for(t=0;t<nTypePairs;t++) fout << " " << neigh[u][t].std()/sqrt(N);
      fout << endl;
      fout.close();
      if(debug) cout << "*** COORDNUM computation for timestep " << timestep << " ENDED ***\n";
    }

    //---------------------- Minimum atomic distance ---------------------------------//
    void init_rmin(string string_rmin_out_)
    {
      if(debug) cout<<"*** Initializing RMIN within "<<myName<<"***\n";
      string_rmin_out = string_rmin_out_;
      ss.str(std::string()); ss << string_rmin_out << tag << ".dat"; fout.open(ss.str(), ios::out);
      fout << "# Minimum atomic distance. # cutoffs = ";
      for(int t=0;t<nTypePairs;t++) fout<<rcut[0][t]<<" ";
      fout<<endl;
      fout.close();
    }

    void compute_rmin(int timestep, vector<ptype> ps)
    {
      ntype rSq, rminSq = rcut[0].max();
      int i,j, k;
      if(debug) cout << "*** RMIN computation for timestep " << timestep << " STARTED ***\n";
      for(i=0;i<N;i++){
        for(k=0;k<ps[i].neigh_list[0].size();k++){
          j = ps[i].neigh_list[0][k];
          if(j>i) continue; // avoid double counting!
          rSq = ps[i].rijSq_list[0][k];
          if(rSq<rminSq) rminSq=rSq;
        }
      }
      ss.str(std::string()); ss << string_rmin_out << tag << ".dat"; fout.open(ss.str(), ios::app);
      fout << sqrt(rminSq) << endl;
      fout.close();
      if(debug) cout << "*** RMIN computation for timestep " << timestep << " ENDED ***\n";
      return;
    }

    //-------------------- Maximum atomic distance --------------------------------------//
    void init_rmax(string string_rmax_out_)
    {
      if(debug) cout<<"*** Initializing RMAX within "<<myName<<"***\n";
      string_rmax_out = string_rmax_out_;
      ss.str(std::string()); ss << string_rmax_out << tag << ".dat"; fout.open(ss.str(), ios::out);
      fout << "# Maximum atomic distance within PBC. # cutoffs = ";
      for(int t=0;t<nTypePairs;t++) fout<<rcut[0][t]<<" ";
      fout<<endl;
      fout.close();
    }
    void print_rmax(int timestep)
    {
      if(debug) cout << "*** PRINT RMAX for timestep " << timestep << " STARTED ***\n";
      ss.str(std::string()); ss << string_rmax_out << tag << ".dat"; fout.open(ss.str(), ios::app);
      fout << sqrt(rmaxSq) << endl;
      fout.close();
      if(debug) cout << "*** PRINT RMAX for timestep " << timestep << " ENDED ***\n";
      return;
    }
};
#endif
