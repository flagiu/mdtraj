#ifndef _NEIGH_AND_BOND_LIST_H_
#define _NEIGH_AND_BOND_LIST_H_
#include<algorithm>
using namespace std;

//---- "init" functions: prepare arrays, averages and headers for output files
//---- "compute" functions: compute quantities for each frame and total average

//------------------------------ Neighbour and Bond list --------------------------------------------//

#define MAX_NSHELL 3

template <class ntype, class ptype>
class Neigh_and_Bond_list
{
  using vec=myvec<ntype,3>;
  using mat=mymatrix<ntype,3,3>;
  private:
    int p1half, p2half, p1,p2;
    string string_cn_out, string_rmin_out, string_rmax_out, myName, tag;
    fstream fout;
    stringstream ss;

    ntype fcut_generic(ntype x, int pow1, int pow2)
    {
      ntype x1=1.0, x2=1.0;
      for(auto i=0; i<pow2; i++) {
        x2 *= x;
        if(i==pow1-1) x1=x2;
      }
      return (1.0-x1) / (1.0-x2); // x1 = x^p1, x2 = x^p2
    }

  public:
    int Nshell, N;
    ntype rcut[MAX_NSHELL], rcutSq[MAX_NSHELL], rmaxSq;
    vecflex<ntype> neigh[MAX_NSHELL];
    vector<int> bond_list[MAX_NSHELL];  // stores all bonds encoded into an integer through ij2int()

    Neigh_and_Bond_list(){
      myName = "NEIGH & BOND List";
    }
    virtual ~Neigh_and_Bond_list(){}

    int ij2int(int i, int j, int N){
      return (i<j ? N*i+j : N*j+i); // i<j = 0,...,N-1
    }
    int int2i(int x, int N){ return x/N; }
    int int2j(int x, int N){ return x%N; }
    ntype fcut(ntype xSq) { return fcut_generic(xSq,p1half,p2half); }

    void print_bond_summary(vector<ptype> ps)
    {
      int i,u,ii;
      for(u=0;u<Nshell;u++)
      {
        for(i=0;i<N;i++)
        {
          cout << "  Shell u="<<u<<" of particle i="<<i<<" contains "<<ps[i].neigh_list[u].size()<<" neighbours:\n   ";
          for(ii=0;ii<ps[i].neigh_list[u].size();ii++) cout<<ps[i].neigh_list[u][ii]<<" ";
          cout << endl;
        }
        cout << "> System has "<<bond_list[u].size()<<" bonds within rcut="<<rcut[u]<<endl;
      }
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

    void init(int Nshell_, ntype *rcut_, int p1half_, int N_)
    {
      Nshell=Nshell_;
      N=N_;
      for(int u=0;u<Nshell;u++)
      {
        rcut[u]=rcut_[u];
        rcutSq[u] = rcut[u]*rcut[u];
      }
      p1half=p1half_;
      p2half = 2*p1half;
      p1 = 2*p1half;
      p2 = 2*p2half;
    }

    void build(int timestep, vector<ptype>& ps, mat box, mat boxInv, bool debug)
    { // NOTA BENE: ps deve essere passato con &, perche' dobbiamo modificare le sue variabili
      int u,i,j, a;
      vec rij, rij_mic;
      ntype rijSq, rijSq_mic;
      vector<ntype> bond_list_rijSq[MAX_NSHELL];
      if(debug) cout << "\n*** "<<myName<<" computation STARTED ***\n";

      //---- Reset counters and lists ----//
      for(u=0;u<Nshell;u++)
      {
        if(neigh[u].length()!=N) neigh[u].resize(N);
        bond_list[u].clear();
        bond_list_rijSq[u].clear();
      }

      for(i=0;i<N;i++)
      {
        for(u=0;u<Nshell;u++)
        {
          ps[i].neigh_list[u].clear();
          ps[i].rij_list[u].clear();
          ps[i].rijSq_list[u].clear();
          neigh[u][i]=0.0;
        }
      }
      if(debug) cout << " * Reset counters and lists DONE\n";

      //---- Build neighbour list (and save rij vectors) and compute rmax ----//
      rmaxSq=0.0;
      for(i=0;i<N;i++)
      {
        for(j=i+1;j<N;j++)
        {
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
          
          for(u=0;u<Nshell;u++)
          {
            if(rijSq <= rcutSq[u])
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
        for(u=0;u<Nshell;u++)
        {
          /*
          if(debug) {
            cout << "\nUnsorted arrrays:\n  ";
            for (auto x: ps[i].rijSq_list[u]) std::cout << x << ' ';
            cout << "\n  ";
            for (auto x: ps[i].neigh_list[u]) std::cout << x << ' ';
            cout << "\n  ";
          } */
          sort_by_distance_3vec( &ps[i].rijSq_list[u], &ps[i].neigh_list[u], &ps[i].rij_list[u] );
          /*
          if(debug) {
            cout << "\nSorted arrrays:\n  ";
            for (auto x: ps[i].rijSq_list[u]) std::cout << x << ' ';
            cout << "\n  ";
            for (auto x: ps[i].neigh_list[u]) std::cout << x << ' ';
            cout << "\n  ";
          }
          */
        }
      }
      // sort bond_list according to bond_list_rijSq
      for(u=0;u<Nshell;u++) sort_by_distance_2vec( &bond_list_rijSq[u], &bond_list[u]);
      if(debug)
      {
        cout << " * Build & Sort the Neighbour list & the Bond list DONE\n";
        //print_bond_summary(ps);
        cout << "*** "<<myName<<" computation for timestep " << timestep << " ENDED ***\n\n";
      }
      return;
    }

    //----------------- Coordination Number (done on 1st shell) ------------//
    void init_coordnum(string string_cn_out_, string tag_, bool debug)
    {
      if(debug) cout<<"*** Initializing COORDNUM within "<<myName<<"***\n";
      string_cn_out = string_cn_out_;
      tag = tag_;
      ss.str(std::string()); ss << string_cn_out << tag << ".dat"; fout.open(ss.str(), ios::out);
      fout << "#Timestep, Particle index, Coordination number. # cutoff = "<<rcut[0]<<endl;
      fout.close();

      ss.str(std::string()); ss << string_cn_out << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "#Timestep, average coordination number, fluctuations. # cutoff = "<<rcut[0]<<endl;
      fout.close();

      if(neigh[0].length()!=N) neigh[0].resize(N);
      for(int i=0;i<N;i++) neigh[0][i]=0.0 ; // start counter
      if(debug) cout<<"*** Initialization completed ***\n";
    }

    void compute_coordnum(int timestep, vector<ptype> ps, bool debug)
    {
      if(debug) cout << "*** COORDNUM computation for timestep " << timestep << " STARTED ***\n";
      int i, j, k;
      ntype fval, rijSq;
      vec rij;
      for(i=0;i<N;i++){
        for(k=0;k<ps[i].neigh_list[0].size();k++){ // search in 1st shell neighbour list
          j = ps[i].neigh_list[0][k];
          if(j>i) continue; // avoid double counting!
          rij = ps[i].rij_list[0][k];
          rijSq = ps[i].rijSq_list[0][k];
          // fval = fcut( rijSq/rcutSq[0] );    // smooth
          fval = ( rijSq <= rcutSq[0] ? 1.0 : 0.0 );         // sharp
          neigh[0][i] += fval;
          neigh[0][j] += fval;
        }
      }
      ss.str(std::string()); ss << string_cn_out << tag << ".dat"; fout.open(ss.str(), ios::app);
      for(i=0;i<N;i++) fout << timestep << " " << i << " " << neigh[0][i] << endl;
      fout.close();
      ss.str(std::string()); ss << string_cn_out << tag << ".ave"; fout.open(ss.str(), ios::app);
      fout << timestep << " " << neigh[0].mean() << " " << neigh[0].std()/sqrt(N) << endl;
      fout.close();
      if(debug) cout << "*** COORDNUM computation for timestep " << timestep << " ENDED ***\n";
    }

    //---------------------- Minimum atomic distance ---------------------------------//
    void init_rmin(string string_rmin_out_, string tag_, bool debug)
    {
      if(debug) cout<<"*** Initializing RMIN within "<<myName<<"***\n";
      string_rmin_out = string_rmin_out_;
      tag = tag_;
      ss.str(std::string()); ss << string_rmin_out << tag << ".dat"; fout.open(ss.str(), ios::out);
      fout << "# Minimum atomic distance. # cutoff = " << rcut[0]<< endl;
      fout.close();
    }

    void compute_rmin(int timestep, vector<ptype> ps, bool debug)
    {
      ntype rSq, rminSq = rcut[0];
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
    void init_rmax(string string_rmax_out_, string tag_, bool debug)
    {
      if(debug) cout<<"*** Initializing RMAX within "<<myName<<"***\n";
      string_rmax_out = string_rmax_out_;
      tag = tag_;
      ss.str(std::string()); ss << string_rmax_out << tag << ".dat"; fout.open(ss.str(), ios::out);
      fout << "# Maximum atomic distance within PBC. # cutoff = " << rcut[0]<< endl;
      fout.close();
    }
    void print_rmax(int timestep, bool debug)
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
