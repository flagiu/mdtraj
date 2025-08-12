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
    string string_cn_out, string_nna_out,string_nnd_out, string_rmin_out, string_rmax_out;
    string log_file, myName, tag;
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
    int Nsphere, N, nTypes, nTypePairs, max_num_nn, max_num_nna,max_num_nnd;
    vecflex<ntype> rcut[MAX_NSPHERE], rcutSq[MAX_NSPHERE];
    vector< vecflex<ntype> > neigh[MAX_NSPHERE]; // Nspheres X nTypePairs X N
    vecflex<ntype> neigh_anytype[MAX_NSPHERE]; // Nspheres X N (agnostic of types)
    vector<int> bond_list[MAX_NSPHERE];  // stores all bonds (i,j), encoded into an integer through ij2int()

    vector<int> cluster_of_particle, head_of_cluster, next_of_particle, size_of_cluster;
    int sphere_for_clustering, num_clusters, maxClusterSize;
    vector<std::size_t> cluster_permutation_by_size;
    vecflex<ntype> rcut_clustering, rcutSq_clustering;

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
        cerr << "ERROR: bond not found with (i,j)=("<<i<<","<<j<<")\n";
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
        cerr << "*** Using default RCUT for all type pairs ***\n";
        for(u=0;u<Nsphere;u++) {
          for(t=0;t<nTypePairs;t++){
            rcut[u][t]=cutoffDefault[u][0];
            rcutSq[u][t]=SQUARE(rcut[u][t]);
          }
        }
      }
      else {
        if(debug) cerr << "*** Reading RCUT from file "<<s_rcut<<" ***\n";
        fout.open(s_rcut, ios::in);
        if(fout.fail()) // checks to see if file opended
        {
          cerr << "ERROR: could not open RCUT file '"<<s_rcut<<"'\n\n";
          exit(1);
        }
        for(u=0;u<Nsphere;u++)
        {
          getline(fout,line);
          ss.str(std::string()); ss.clear(); // clear the string stream!
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
        if(debug) cerr << "*** Reading RCUT file DONE ***\n";
      }
      // output r_cut to the log file
      ss.str(std::string()); ss << log_file << tag; fout.open(ss.str(), ios::app);
      fout << "#-------- RCUT summary --------#\n";
      for(u=0;u<Nsphere;u++)
      {
        fout << "SPHERE u="<<u<<" : ";
        rcut[u].show(fout);
      }
      fout << "#------------------------------#\n";
      fout.close();
    }

    void read_rcut_clustering_from_file(string s_rcut, vecflex<ntype> cutoffDefault)
    {
      int u,t;
      string line, x;
      if(s_rcut=="__NOT_DEFINED__") {
        cerr << "*** Using default RCUT_CLUSTERING for all type pairs ***\n";
        for(t=0;t<nTypePairs;t++){
          rcut_clustering[t]=cutoffDefault[0];
          rcutSq_clustering[t]=SQUARE(rcut_clustering[t]);
        }
      }
      else {
        if(debug) cerr << "*** Reading RCUT_CLUSTERING from file "<<s_rcut<<" ***\n";
        fout.open(s_rcut, ios::in);
        if(fout.fail()) // checks to see if file opended
        {
          cerr << "ERROR: could not open RCUT_CLUSTERING file '"<<s_rcut<<"'\n\n";
          exit(1);
        }
        getline(fout,line);
        ss.str(std::string()); ss.clear(); // clear the string stream!
        ss << line;
        for(t=0;t<nTypePairs;t++)
        {
          ss >> x;
          rcut_clustering[t] = stof(x);
          rcutSq_clustering[t] = SQUARE(rcut_clustering[t]);
        }
        ss.str(std::string()); ss.clear(); // clear the string stream!
        fout.close();
        if(debug) cerr << "*** Reading RCUT_CLUSTERING file DONE ***\n";
      }
      // output to the log file
      ss.str(std::string()); ss << log_file << tag; fout.open(ss.str(), ios::app);
      fout << "#--- RCUT_CLUSTERING summary ---#\n";
      rcut_clustering.show(fout);
      fout << "#-------------------------------#\n";
      fout.close();
    }

    void print_bond_summary(vector<ptype> ps)
    {
      int i,u,ii;
      ss.str(std::string()); ss << log_file << tag; fout.open(ss.str(), ios::app);
      fout << "#---------------------- BOND SUMMARY --------------------------#\n";
      //for(u=0;u<Nsphere;u++)
      for(u=0;u<1;u++)
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
      //cerr << " Summary of neighbour and bond lists saved into log file: "<<log_file <<tag <<endl;
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
      if(debug) cerr << "\n*** "<<myName<<" computation STARTED ***\n";

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
      if(debug) cerr << " * Reset counters and lists DONE\n";

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
        cerr << " * Build neighbour and bond lists DONE\n";
        cerr << "*** "<<myName<<" computation for timestep " << timestep << " ENDED ***\n\n";
      }
      return;
    }

    //----------------- Coordination Number (done on 1st sphere) ------------//
    void init_coordnum(string string_cn_out_)
    {
      if(debug) cerr<<"*** Initializing COORDNUM within "<<myName<<"***\n";
      string_cn_out = string_cn_out_;
      if(verbose)
      {
        ss.str(std::string()); ss << string_cn_out << tag << ".dat"; fout.open(ss.str(), ios::out);
        fout << "#Timestep | Particle idx | Particle type | Coordination number for each neighbor type ... . # cutoffs = ";
        for(int t=0;t<nTypePairs;t++) fout<<rcut[0][t]<<" ";
        fout<<endl;
        fout.close();
      }

      ss.str(std::string()); ss << string_cn_out << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "#Timestep | <coordination number> and its fluctuations for each type pair. # cutoffs = ";
      for(int t=0;t<nTypePairs;t++) fout<<rcut[0][t]<<" ";
      fout<<endl;
      fout.close();

      const int u=0;
      if(neigh[u].size()!=nTypePairs) neigh[u].resize(nTypePairs);
      for(int j=0;j<nTypePairs;j++){
        if(neigh[u][j].length()!=N) neigh[u][j].resize(N);
      }
      //for(int i=0;i<N;i++) neigh[0][i]=0.0 ; // start counter
      if(debug) cerr<<"*** Initialization completed ***\n";
    }

    void compute_coordnum(int timestep, vector<ptype> ps)
    {
      if(debug) cerr << "*** COORDNUM computation for timestep " << timestep << " STARTED ***\n";
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
          fout << timestep << " " << i << " " << ps[i].label;
          for(int tj=0;tj<nTypes;tj++) {
		  t=types2int(ps[i].label, tj);
		  fout << " " << neigh[0][t][i];
	  }
          fout << endl;
        }
        fout.close();
      }

      ss.str(std::string()); ss << string_cn_out << tag << ".ave"; fout.open(ss.str(), ios::app);
      fout << timestep;
      ntype x1,x2,count;
      for(int ti=0;ti<nTypes;ti++){
	      for(int tj=ti;tj<nTypes;tj++){
		      t=types2int(ti,tj);
		      x1=x2=count=0.0;
		      for(i=0;i<N;i++){
			if(ps[i].label==ti){
				x1 += neigh[u][t][i];
				x2 += SQUARE(neigh[u][t][i]);
				count += 1;
			}
		      }
		      x1/=count;
		      x2/=count;
		      fout << " " << x1 << " " << sqrt( (x2-x1*x1)/(count-1));
	      }
      }
//      for(t=0;t<nTypePairs;t++) fout << " " << neigh[u][t].mean();
//      for(t=0;t<nTypePairs;t++) fout << " " << neigh[u][t].std()/sqrt(N-1);
      fout << endl;
      fout.close();
      if(debug) cerr << "*** COORDNUM computation for timestep " << timestep << " ENDED ***\n";
    }

    //----------------- Nearest Neighbours distances (done on 1st sphere) ------------//
    void init_nearest_neigh_dists(string string_nnd_out_, int max_num_nnd_)
    {
      max_num_nnd=max_num_nnd_;
      if(debug) cerr<<"*** Initializing NearestNeighDistances within "<<myName<<"***\n";
      cerr<<"WARNING: NearestNeighDistances produces a large output file\n";
      string_nnd_out = string_nnd_out_;
      ss.str(std::string()); ss << string_nnd_out << tag << ".dat"; fout.open(ss.str(), ios::out);
      fout << "#TimeStep | ParticleIdx | ParticleType | (NeighbourType, DistanceSquared) for the first "<<max_num_nnd<<" neighs, sorted by distance | # cutoffs = ";
      for(int t=0;t<nTypePairs;t++) fout<<rcut[0][t]<<" ";
      fout<<endl;
      fout.close();
      if(debug) cerr<<"*** Initialization completed ***\n";
    }

    void compute_nearest_neigh_dists(int timestep, vector<ptype> ps)
    {
      if(debug) cerr << "*** NearestNeighDistances computation for timestep " << timestep << " STARTED ***\n";
      int i, j, k, u;
      ntype rijSq;

      ss.str(std::string()); ss << string_nnd_out << tag << ".dat"; fout.open(ss.str(), ios::app);

      u=0; // 1st sphere neighbours
      for(i=0;i<N;i++){
        fout << timestep << " " << i << " " << ps[i].label;
        for(k=0;k<ps[i].neigh_list[u].size()&&k<max_num_nnd;k++){ // search in neigh list
          j = ps[i].neigh_list[u][k];
          fout << " " << ps[j].label << " " << ps[i].rijSq_list[u][k];
        }
        fout << endl;
      }

      fout.close();
      if(debug) cerr << "*** NearestNeighDistances computation for timestep " << timestep << " ENDED ***\n";
    }

    //----------------- Nearest Neighbours angles (done on 1st sphere) ------------//
    void init_nearest_neigh_angles(string string_nna_out_, int max_num_nna_)
    {
      max_num_nna=max_num_nna_;
      if(debug) cerr<<"*** Initializing NearestNeighAngles within "<<myName<<"***\n";
      cerr<<"WARNING: NearestNeighAngles produces a large output file\n";
      string_nna_out = string_nna_out_;
      ss.str(std::string()); ss << string_nna_out << tag << ".dat"; fout.open(ss.str(), ios::out);
      fout << "#TimeStep | Particle0Idx | Particle0Type | (Neighbour1Type, Neighbour2Type, cos(102)) for the first "<<max_num_nna<<" neighs, sorted by r1<=r2 | # cutoffs = ";
      for(int t=0;t<nTypePairs;t++) fout<<rcut[0][t]<<" ";
      fout<<endl;
      fout.close();
      if(debug) cerr<<"*** Initialization completed ***\n";
    }

    void compute_nearest_neigh_angles(int timestep, vector<ptype> ps)
    {
      if(debug) cerr << "*** NearestNeighAngles computation for timestep " << timestep << " STARTED ***\n";
      int i, j,jj, k,kk, u;
      vec rij, rik;
      ntype rijSq, rikSq, cos;

      ss.str(std::string()); ss << string_nna_out << tag << ".dat"; fout.open(ss.str(), ios::app);

      u=0; // 1st sphere neighbours
      for(i=0;i<N;i++){
        fout << timestep << " " << i << " " << ps[i].label;
        for(jj=0;jj<ps[i].neigh_list[u].size()&&jj<max_num_nna;jj++){ // search in neigh list
          j = ps[i].neigh_list[u][jj];
          rij = ps[i].rij_list[u][jj];
          rijSq = ps[i].rijSq_list[u][jj];
          for(kk=jj+1;kk<ps[i].neigh_list[u].size()&&kk<max_num_nna;kk++){ // search in neigh list
            k = ps[i].neigh_list[u][kk];
            rik = ps[i].rij_list[u][kk];
            rikSq = ps[i].rijSq_list[u][kk];
            cos = (rij*rik) / sqrt(rijSq*rikSq);
            fout << " " << ps[j].label << " " << ps[k].label << " " << cos;
          }
        }
        fout << endl;
      }

      fout.close();
      if(debug) cerr << "*** NearestNeighAngles computation for timestep " << timestep << " ENDED ***\n";
    }

    //---------------------- Minimum atomic distance ---------------------------------//
    void init_rmin(string string_rmin_out_)
    {
      if(debug) cerr<<"*** Initializing RMIN within "<<myName<<"***\n";
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
      if(debug) cerr << "*** RMIN computation for timestep " << timestep << " STARTED ***\n";
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
      if(debug) cerr << "*** RMIN computation for timestep " << timestep << " ENDED ***\n";
      return;
    }

    //-------------------- Maximum atomic distance --------------------------------------//
    void init_rmax(string string_rmax_out_)
    {
      if(debug) cerr<<"*** Initializing RMAX within "<<myName<<"***\n";
      string_rmax_out = string_rmax_out_;
      ss.str(std::string()); ss << string_rmax_out << tag << ".dat"; fout.open(ss.str(), ios::out);
      fout << "# Maximum atomic distance within PBC. # cutoffs = ";
      for(int t=0;t<nTypePairs;t++) fout<<rcut[0][t]<<" ";
      fout<<endl;
      fout.close();
    }
    void print_rmax(int timestep)
    {
      if(debug) cerr << "*** PRINT RMAX for timestep " << timestep << " STARTED ***\n";
      ss.str(std::string()); ss << string_rmax_out << tag << ".dat"; fout.open(ss.str(), ios::app);
      fout << sqrt(rmaxSq) << endl;
      fout.close();
      if(debug) cerr << "*** PRINT RMAX for timestep " << timestep << " ENDED ***\n";
      return;
    }

    //-------------------- Clusterize --------------------------------------//
    void init_clusterize(string string_cluster_out, string s_rcut_clusters, vecflex<ntype> cutoffDefault){
      rcut_clustering.resize(nTypePairs);
      rcutSq_clustering.resize(nTypePairs);
      read_rcut_clustering_from_file(s_rcut_clusters, cutoffDefault);
      // decide which neighbour list sphere fits the specified cutoff:
      sphere_for_clustering=-1;
      int bad_tp=-1;
      bool ok=false;
      for(int sphere=0;(sphere<Nsphere)&&(!ok);sphere++){
        ok=true;
        for(int tp=0;tp<(nTypePairs)&&(ok);tp++){
          if(rcut_clustering[tp]>rcut[sphere][tp]) {
            ok=false;
            bad_tp=tp;
          }
        }
        if(ok) sphere_for_clustering=sphere;
      }

      if(!ok) {
        cerr<<"ERROR: rcut_clustering[tp] > rcut[sphere][tp] for all spheres";
        cerr<<"(at least for type pair tp="<<bad_tp<<")!\n";
        cerr<<"       Please specify a smaller rcut_clustering or a larger rcut\n";
        exit(1);
      }

      ss.str(std::string()); ss<<string_cluster_out<<".maxsize"<<tag<<".dat"; fout.open(ss.str(), ios::out);
      fout<<"# Timestep, max_cluster_size\n";
      fout.close();

      ss.str(std::string()); ss<<string_cluster_out<<".num"<<tag<<".dat"; fout.open(ss.str(), ios::out);
      fout<<"# Timestep, num_clusters\n";
      fout.close();

      ss.str(std::string()); ss<<string_cluster_out<<".size"<<tag<<".dat"; fout.open(ss.str(), ios::out);
      fout<<"# Timestep, cluster_size_for_each_cluster\n";
      fout.close();
    }

    void clusterize(int timestep, vector<ptype> particles, string string_cluster_out,
                    vecflex<ntype> particle_observable, ntype p_o_threshold){
      // clustering criterion: particle_observable[i]>p_o_threshold AND |r[i]-r[j]|<rcutSq_clustering
      if(sphere_for_clustering<0) {
        cerr<<"ERROR: bad initialization of sphere_for_clustering in init_clusterize() !\n";
        exit(1);
      }
      const int num_bonds=bond_list[sphere_for_clustering].size(); //ordered list i<j
      if(particles.size()!=N){
        cerr<<"ERROR: size of particles != N\n";
        cerr<<"       size of particles = "<<particles.size()<<endl;
        cerr<<"                       N = "<<N<<endl;
        exit(1);
      }
      if(particle_observable.length()!=N){
        cerr<<"ERROR: size of particle_observable < size of particles\n";
        cerr<<"       size of particle_observable = "<<particle_observable.length()<<endl;
        cerr<<"       size of particles           = "<<N<<endl;
        exit(1);
      }
      int i,j,k, typePair, ci,cj,ck, bond_idx,bond_encoded;
      bool found;

      if(debug) cerr<<"*** STARTED computation of CLUSTERS ***\n";
      if(debug) cerr<<"  string_cluster_out = "<<string_cluster_out<<endl;

      // cluster to which each particle belongs to
      // (N elements, whose value is in [0,num_clusters-1] when assigned, -1 otherwise)
      if(cluster_of_particle.size()<N) cluster_of_particle.resize(N);
      // first particle in this cluster
      // (num_clusters<=N elements, whose value is in [0,N-1] when assigned, -1 otherwise)
      if(head_of_cluster.size()<N) head_of_cluster.resize(N);
      // next particle in the cluster of particle i
      if(next_of_particle.size()<N) next_of_particle.resize(N);
      // size of cluster (num_clusters<=N elements, whose value is non-negative when assigned, 0 otherwise)
      if(size_of_cluster.size()<N) size_of_cluster.resize(N);
      if(debug) cerr<<"  Resize check Done\n";

      // start with 1 single-particle-cluster for each crystalline particle
      num_clusters=0;
      for(i=0;i<N;i++){
        next_of_particle[i]=-1;   // particle i don't have any next
        if(particle_observable[i]>p_o_threshold){
          //if(debug) cerr<<"  New cluster from 1 particle: i,Oi = "<<i<<" "<<particle_observable[i]<<endl;
          cluster_of_particle[i]=num_clusters;
          head_of_cluster[num_clusters]=i;
          size_of_cluster[num_clusters]=1;
          num_clusters++;
        } else {
          cluster_of_particle[i]=-1; // the particle doesn't belong to any cluster
          head_of_cluster[i]=-1;    // the clusters don't have any head (they don't exist)
          size_of_cluster[i]=0;
        }
      }

      if(debug) cerr<<"  Single-particle clusters initialization Done\n";

      for(bond_idx=0;bond_idx<num_bonds;bond_idx++){
        bond_encoded=bond_list[sphere_for_clustering][bond_idx];
        i = int2i(bond_encoded,N);
        j = int2j(bond_encoded,N);
        ci=cluster_of_particle[i];
        cj=cluster_of_particle[j];

        // Requirement 1: Both particles' observable > threshold,
        if(particle_observable[i]<p_o_threshold) continue;
        if(particle_observable[j]<p_o_threshold) continue;

        /*
        if(debug) {
          cerr<<"  bond_idx,i,j,ci,cj,Oi,Oj = "<<bond_idx<<" "<<i<<" "<<j<<" "<<" "<<ci<<" "<<cj<<
              " "<<particle_observable[i]<<" "<<particle_observable[j]<<endl;
        }
        */

        // Requirement 2: Closer than a type-dependent cutoff:
        typePair = types2int(particles[i].label, particles[j].label);
        // // 1) find j in i's neighbour list
        found=false;
        for(k=0;(!found) && k<particles[i].neigh_list[sphere_for_clustering].size();k++){
          if(j==particles[i].neigh_list[sphere_for_clustering][k]) found=true;
        }
        if(!found) {
          cerr<<"ERROR: i's list does not contain neighbour j. i="<<i<<" j="<<j<<endl;
          exit(1);
        }
        // // 2) check if they are close enough for clustering
        if(particles[i].rijSq_list[sphere_for_clustering][k]>rcutSq_clustering[typePair]) continue;

        // Cases (update: actually only the case ci>=0 && cj>=0 will occurr):
        if(ci<0 && cj<0) {
          // new cluster of 2 elements (not possible!)
          //if(debug) cerr<<"  New cluster from 2 independent particles\n";
          cluster_of_particle[i]=cluster_of_particle[j]=num_clusters;
          head_of_cluster[num_clusters]=i; //take arbitrarily i as head...
          next_of_particle[i]=j; //... and j as next of i
          size_of_cluster[num_clusters]=2;
          num_clusters++; // num of clusters increases
        } else if(ci>=0 && cj<0) {
          // prepend particle j to cluster of particle i
          //if(debug) cerr<<"  Prepending to cluster: ci<--j\n";
          cluster_of_particle[j]=ci;
          next_of_particle[j]=head_of_cluster[ci]; // add j as head...
          head_of_cluster[ci]=j; //... and set the former head as j's next
          size_of_cluster[ci]++;
          //num_clusters+=0;
        } else if(ci<0 && cj>=0) {
          // same with i<-->j
          //if(debug) cerr<<"  Prepending to cluster: i-->cj\n";
          cluster_of_particle[i]=cj;
          next_of_particle[i]=head_of_cluster[cj];
          head_of_cluster[cj]=i;
          size_of_cluster[cj]++;
          //num_clusters+=0;
        } else if(ci>=0 && cj>=0 && ci!=cj){
          // merge 2 DIFFERENT clusters into one (move j's into i's)
          //if(debug) cerr<<"  Merging clusters: ci<--cj\n";
          // prepend particles in j's cluster to i's cluster
          int head_of_cluster_i=head_of_cluster[ci]; // store the head of i's cluster
          k=head_of_cluster[cj]; // take j's head...
          head_of_cluster[ci]=k; // ... move it to i's head...
          cluster_of_particle[k]=ci; // ...and assign it to i's cluster...
          size_of_cluster[ci]++; // ...and increase the counter
          while(next_of_particle[k]>=0) {
            k=next_of_particle[k]; // loop over j's cluster...
            cluster_of_particle[k]=ci; //... and assign it to i's cluster...
            size_of_cluster[ci]++; // ...and increase the counter
          }
          next_of_particle[k] = head_of_cluster_i; // move i's head as next of last particle in j's cluster

          // shift by -1 the clusters ck with ck>cj
          for(ck=cj+1;ck<num_clusters;ck++){
            head_of_cluster[ck-1]=head_of_cluster[ck]; // copy head to previous cluster
            size_of_cluster[ck-1]=size_of_cluster[ck]; // copy size to previous cluster
            k=head_of_cluster[ck];
            cluster_of_particle[k]=ck-1;  // shift by -1 the cluster for each particle
            while(next_of_particle[k]>=0) {
              k=next_of_particle[k];
              cluster_of_particle[k]=ck-1;
            }
          }
          head_of_cluster[num_clusters-1]=-1; // reset the last cluster head
          size_of_cluster[num_clusters-1]=0; // reset the last cluster size
          num_clusters--; // decrease num of clusters by 1
        } else {
          // ci=cj>=0: do nothing if particles already belong to same cluster
        }

      }

      if(debug) cerr<<"  Clusterization Done\n";

      sort_clusters_by_size_and_get_maxsize();

      if(debug) cerr<<"  Sorting by size Done\n";

      ss.str(std::string()); ss<<string_cluster_out<<".size"<<tag<<".dat"; fout.open(ss.str(), ios::app);
      fout<<timestep<<" ";
      for(ck=0;ck<num_clusters;ck++){
        fout<<size_of_cluster[cluster_permutation_by_size[ck]]<<" ";
      }
      fout<<endl;
      fout.close();

      ss.str(std::string()); ss<<string_cluster_out<<".num"<<tag<<".dat"; fout.open(ss.str(), ios::app);
      fout<<timestep<<" "<<num_clusters<<endl;
      fout.close();

      ss.str(std::string()); ss<<string_cluster_out<<".maxsize"<<tag<<".dat"; fout.open(ss.str(), ios::app);
      fout<<timestep<<" "<<maxClusterSize<<endl;
      fout.close();

      if(debug) cerr<<"  Output Done\n";

      if(debug) cerr<<"*** COMPLETED computation of CLUSTERS ***\n";

    }

    int get_cluster_size(int cluster_index){
      if(cluster_index>=num_clusters) {
        cerr<<"ERROR: cluster_index >= num_clusters\n";
        cerr<<"       cluster_index = "<<cluster_index<<endl;
        cerr<<"       num_clusters  = "<<num_clusters<<endl;
        exit(1);
      }
      int k, count=0;
      //if(debug) { cerr<<"  Cluster "<<cluster_index<<" : ";}
      k=head_of_cluster[cluster_index];
      if(k>=0) {
        count++;
        //if(debug) { cerr<<k<<" ";}
        while(next_of_particle[k]>=0) {
          k=next_of_particle[k];
          count++;
          //if(debug) { cerr<<k<<" ";}
        }
      }
      //if(debug) { cerr<<endl;}
      if(count!=size_of_cluster[cluster_index]) {
        cerr<<"ERROR: unexpected size of cluster "<<cluster_index<<":\n";
        cerr<<"       from get_cluster_size() : "<<count<<endl;
        cerr<<"       from size_of_cluster    : "<<size_of_cluster[cluster_index]<<endl;
        exit(1);
      }
      return count;
    }

    void sort_clusters_by_size_and_get_maxsize()
    {
      int foo;
      cluster_permutation_by_size.resize(num_clusters); // vector of 0,1,...,size-1
      for(std::size_t i=0;i<num_clusters;i++) {
        cluster_permutation_by_size[i]=i;
        foo=get_cluster_size((int)i);
      }
      // sort this vector according to cluster size
      std::sort(
        cluster_permutation_by_size.begin(),
        cluster_permutation_by_size.end(),
        [&](std::size_t i, std::size_t j){
          return size_of_cluster[i] > size_of_cluster[j];
        }
      );

      maxClusterSize=( num_clusters==0?0:size_of_cluster[cluster_permutation_by_size[0]] );

      return;
    }
};
#endif
