#ifndef _MDTRAJ_H_
#define _MDTRAJ_H_

/*
 TO DO: 
 Optimization ideas:

*/

#include<cstdlib>
#include<cstring>
#include<complex>
#include<cmath>
#include "utility.hpp"
#include "vecflex.hpp"
#include "particle.hpp"
#include "Ycomplex.hpp"
using ntype=double;
using ptype=particle<ntype>;

template<class ntype, class ptype>
class Trajectory {
public:
  using vec=myvec<ntype,3>;
  int N; // number of particles
  vec L; // dimensions of an orthorombic simulation box ( ??? centered in 0: -Lx/2 < x < Lx/2 )
  ntype V, mdens, ndens; // volume, mass density, nuerical density
  vector<ptype> ps, ps_new; // vector of particles
  int nframes, timestep, period, rdf_nbins;
  bool c_coordnum, c_bondorient, c_msd, c_rdf; // compute or not
  string s_xyz, tag, s_coordnum, s_bondorient, s_bondcorr, s_nxtal, s_msd, s_ngp, s_rdf; // for file naming
  int l; //angular momentum  for bond orientational parameter
  static const int Nshells=3;
  int p1, p2; // parameters for fcut
  ntype cutoff[Nshells]; // cutoff radii
  vecflex<ntype> neigh[Nshells], ql, Cl_ij, ql_dot, Ql_dot;
  vector< vecflex< complex<ntype> > > qlm; // a collection of l_deg vectors of local average qlm=<Ylm> with a cutoff function
  vector<int> bond_list[Nshells]; // records all bonds encoded into an integer through ij2int()
  vecflex<ntype> r2,r4, r2CM; // for MSD
  vecflex<ntype> rdf_bins, rdf_norm, rdf, rdf_ave, rdf2_ave; // for RDF
  vector<vec> rs;
  bool msdAverageOverTime0;
  bool debug, verbose;

private:
  bool l_is_odd;
  int nlines, t0frame, dtframe, l_deg, p1half, p2half, periodIdx, Nperiods;
  fstream fin, fout;
  stringstream ss;
  ntype cutoffSq[Nshells], invN, qldot_th, rdf_binw;

  int ij2int(int i, int j, int N){
    return (i<j ? N*i+j : N*j+i); // i<j = 0,...,N-1
  }
  int int2i(int x, int N){
    return x/N;
  }
  int int2j(int x, int N){
    return x%N;
  }

public:
  Trajectory(){}
  ~Trajectory(){}

  void print_state(){
    cout << "Summary of parameters:\n";
    cout << " debug = \t " << debug << endl;
    cout << " verbose = \t " << debug << endl;
    cout << " c_coordnum = \t " << c_coordnum << endl;
    cout << " c_bondorient = \t " << c_bondorient << endl;
    cout << " c_msd = \t " << c_msd << endl;
    cout << " c_rdf = \t " << c_rdf << endl;
    cout << " l = \t " << l << endl;
    cout << " L = \t "; L.show();
    cout << " p1 = \t " << p1 << endl;
    cout << " p2 = \t " << p2 << endl;
    cout << " period = \t " << period << endl;
    cout << " rdf_nbins = \t " << rdf_nbins << endl;
    cout << " qldot_th = \t " << qldot_th << endl;
    cout << " tag = \t " << tag << endl;
    cout << " V = \t " << V << endl;
    for(auto  i=0;i<Nshells;i++) cout << " rcut" << i << " = \t " << cutoff[i] << endl;
    cout << " s_xyz = \t " << s_xyz << endl;
    cout << endl;
  }

  void init(int argc, char** argv){
    debug = false;
    verbose = false;
    c_coordnum = false;
    c_bondorient = false;
    c_msd = false;
    c_rdf = false;
    s_xyz="dump.run.xyz";
    s_coordnum="coordnum";
    s_bondorient="boo";
    s_bondcorr="boc";
    s_nxtal="nc";
    s_msd="msd";
    s_ngp="ngp";
    s_rdf="RDF";
    tag="";
    period = -1; // default: don't average over t0 for MSD
    cutoff[0] = 3.75; // 3.6 in glass, 3.75-3.89 in xtal
    cutoff[1] = 5.15;
    cutoff[2] = 8.8;
    l = 4;
    L << 1.0,1.0,1.0;
    p1=2*6;
    p2=2*p1; //2*12;
    qldot_th = 0.65;
    rdf_nbins=0;
    args(argc, argv);
    compute_volume();
    for(auto  i=0;i<Nshells;i++) cutoffSq[i] = cutoff[i]*cutoff[i];
    p1half = p1/2;
    p2half = p2/2;
    l_is_odd = (l%2!=0);
    l_deg = 2*l+1;
    if(debug) { cout << "State after reading args():\n"; print_state(); }
    if(l==0) cout << "[WARNING: bond orientation order parameter with l=0 is always equal to 1.0]\n";
  }
  
  void set_L(ntype Lx, ntype Ly, ntype Lz) {
    L << Lx, Ly, Lz;
    compute_volume();
  }

  void compute_volume() {
    V = L.prod(); // valid only for orthorombic box
  }

  void run()
  {
    PrintProgress printProgress;
    if(verbose || debug) {cout << "Begin of run():\n"; print_state();}
    nlines = getLineCount(s_xyz);
    if(debug) cout << "Read " << nlines << " lines in file " << s_xyz << ". Opening again for reading trajectory.\n";
    fin.open(s_xyz, ios::in);
    read_xyz_frame(fin, true);
    invN = 1.0/N;
    nframes = nlines / (N+2);
    t0frame = timestep;
    if(debug) cout << "Read first frame. Set N = " << N << ", t0frame = " << t0frame << ".\n";
    if(debug) cout << "Deduced nframes = " << nframes << " (assuming N is constant).\n";
    read_xyz_frame(fin, false);
    dtframe = timestep - t0frame;
    if(debug) cout << "Read second frame. Set dtframe = " << dtframe << " (assumed to be constant).\n";
    init_computations();
    if(debug) cout << "Initialized arrays for computations.\n";
    fin.close();
    // Restart reading
    fin.open(s_xyz, ios::in);
    printProgress.init( 0.0 );
    for(int i=0; i<nframes; i++)
    {
      read_xyz_frame(fin, false);
      if(N != ps.size()) { cout << "[ERROR: N has changed]\n"; exit(1); }
      if( (timestep - t0frame)%dtframe != 0) {
        cout << "[ERROR: timestep interval has changed]\n";
        cout << "[t0frame = "<<t0frame<<", dtframe = "<<dtframe<<", timestep = "<<timestep<<"]\n";
        exit(1);
      }
      if(c_coordnum) compute_coord_num();
      if(c_bondorient) compute_boo_boc();
      if(c_msd) compute_msd(i);
      if(c_rdf) compute_rdf(i);
      if(verbose && !debug) printProgress.update( (i+1)/(float)nframes *100.0 );
    }
    fin.close();
    if(debug) cout << "Closed input file.\n";
    if(c_msd) print_msd();
    if(debug || verbose) cout << "\n\nExecution completed.\n\n";
  }

  void read_xyz_frame(fstream &i, bool resetN)
  {
    string line, a,b,c,d;
    getline(i,line); // first line
    istringstream(line) >> a;
    N = stoi(a);
    if(debug) cout << "  Line 1: " << N << " atoms\n";
    getline(i, line); // second line
    istringstream(line) >> b >> c >> d;
    timestep = stoi(d);
    if(debug) cout << "  Line 2: Timestep " << timestep << endl;
    if(resetN) ps.resize(N);
    for(auto &p: ps) p.read_xyz(i); // N particle lines
    if(debug) cout << "  Read particle positions DONE.\n";
  }

  //-------------Trajectory analysis---------------------------//

  void init_computations() {
    ndens = N/V;
    if(c_coordnum)
    {
      neigh[0].resize(N);
      ss.str(std::string()); ss << s_coordnum << tag << ".dat"; fout.open(ss.str(), ios::out);
      fout << "#Timestep, Particle index, Coordination number. # cutoffs = ";
      for(auto  i=0;i<Nshells;i++) fout << cutoff[i]<<", ";
      fout << endl;
      fout.close();
      ss.str(std::string()); ss << s_coordnum << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "#Timestep, average coordination number, fluctuations. # cutoffs = ";
      for(auto  i=0;i<Nshells;i++) fout << cutoff[i]<<", ";
      fout << endl;
      fout.close();
    }
    if(c_bondorient)
    {
      for(auto i=0;i<Nshells;i++) neigh[i].resize(N);
      qlm.resize(l_deg);
      for(auto a=0;a<l_deg;a++) qlm[a].resize(N);
      ql.resize(N);
      ql_dot.resize(N);
      Ql_dot.resize(N);
      ss.str(std::string()); ss << s_bondorient << ".l" << l << tag << ".dat"; fout.open(ss.str(), ios::out);
      fout << "#Timestep, Particle index, Bond order orientation parameter q_l. # l = " << l << ", cutoffs = ";
      for(auto  i=0;i<Nshells;i++) fout << cutoff[i]<<", ";
      fout << endl;
      fout.close();
      ss.str(std::string()); ss << s_bondorient << ".l" << l << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "#Timestep, <q_l>, fluctuations. # l = " << l << ", cutoffs = ";
      for(auto  i=0;i<Nshells;i++) fout << cutoff[i]<<", ";
      fout << endl;
      fout.close();
      ss.str(std::string()); ss << s_bondcorr << ".l" << l << tag << ".dat"; fout.open(ss.str(), ios::out);
      fout << "#Timestep, Particle index, Bond order correlation parameter q_l_dot. # l = " << l << ", cutoffs = ";
      for(auto  i=0;i<Nshells;i++) fout << cutoff[i]<<", ";
      fout << endl;
      fout.close();
      ss.str(std::string()); ss << s_nxtal << ".l" << l << tag << ".dat"; fout.open(ss.str(), ios::out);
      fout << "#Timestep, Nc = number of crystalline particles # l = " << l << ", cutoffs = ";
      for(auto  i=0;i<Nshells;i++) fout << cutoff[i]<<", ";
      fout << endl;
      fout.close();
      ss.str(std::string()); ss << s_bondcorr << ".l" << l << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "#Timestep, <q_l_dot>, fluctuations. # l = " << l << ", cutoffs = ";
      for(auto  i=0;i<Nshells;i++) fout << cutoff[i]<<", ";
      fout << endl;
      fout.close();
      ss.str(std::string()); ss << s_bondcorr << ".l" << l << tag << ".xyz"; fout.open(ss.str(), ios::out);
      fout.close();
      ss.str(std::string()); ss << s_bondcorr << ".l" << l << tag << ".local.ave"; fout.open(ss.str(), ios::out);
      fout << "#Timestep, Particle index, local <q_l_dot(i)>. # l = " << l << ", cutoffs = ";
      for(auto  i=0;i<Nshells;i++) fout << cutoff[i]<<", ";
      fout << endl;
      fout.close();
    }
    if(c_msd){
      if(period>0 && period<nframes*dtframe){
        periodIdx = ceil(period/dtframe);
      } else{
        periodIdx = nframes;
      } // periodIdx = length of the trajectory, in index units
      Nperiods = nframes/periodIdx; // Nperiods = n. of trajectories with different t0
      if(debug) cout << "Set periodIdx = " << periodIdx << ", Nperiods = " << Nperiods << endl;
      rs.resize((N+1)*periodIdx); // store all coordinates at all times + center of mass at all time
      r2.resize(periodIdx-1); // < |r(t) - r(t0)|^2 >
      r4.resize(periodIdx-1); // < |r(t) - r(t0)|^4 >
      r2CM.resize(periodIdx-1); // < |rCM(t) - rCM(t0)|^2 >
      ss.str(std::string()); ss << s_msd << tag << ".traj"; fout.open(ss.str(), ios::out);
      fout << "#Delta timestep, then |r(t)-r(t0)|^2 # block file for different t0\n";
      fout.close();
      ss.str(std::string()); ss << s_msd << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "#Delta timestep, MSD, MSD error, MSD of c.o.m.\n";
      fout.close();
      ss.str(std::string()); ss << s_ngp << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "#Delta timestep, NGP (-0.4: homog. displ., 0.0: Brownian displ., >0: heterog. displ.)\n";
      fout.close();
    }
    if(c_rdf){
      rdf_norm.resize(rdf_nbins);
      rdf_bins.resize(rdf_nbins);
      rdf.resize(rdf_nbins);
      rdf_ave.resize(rdf_nbins);
      rdf2_ave.resize(rdf_nbins);
      rdf_binw = 0.5*(L.abs().min()) / (rdf_nbins-1) ;
      ntype r, shell1, shell2, normalization = (N-1)/V * 4.0*M_PI/3.0 * N/2.0;
      shell1 = 0.0;
      ss.str(std::string()); ss << s_rdf << tag << ".traj"; fout.open(ss.str(), ios::out);
      fout << "# First block: radial distance; other blocks: g(r)\n";
      for( auto i=0; i<rdf_nbins; i++){
        rdf_bins[i] = (i+0.5)*rdf_binw; // take the center of the bin for histogram
        fout << rdf_bins[i] << endl;
        r = (i+1)*rdf_binw; // take the upper side of the bin for shells
        shell2 = r*r*r;
        rdf_norm[i] = normalization*(shell2 - shell1); // normalize to ideal gas radial distribution
        shell1 = shell2;
        rdf_ave[i] = 0.0;
        rdf2_ave[i] = 0.0;
      }
      fout.close();
      ss.str(std::string()); ss << s_rdf << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "# Radial distance, g(r), g(r) error.\n";
      fout.close();
    }
  }

  void compute_coord_num() { //----- Coordination Number -----//
    ntype fval;
    vec rij;
    for(auto i=0;i<N;i++) neigh[0].set(i,0.0); // start counter
    for(auto i=0;i<N-1;i++){
      for(auto j=i+1;j<N;j++){
        rij = ps[j].r - ps[i].r;
        rij = ps[i].mic(rij, L); // minimum image convention
        fval = fcut( rij.sq()/cutoffSq[0], p1half, p2half );
        neigh[0][i] += fval;
        neigh[0][j] += fval;
//        if( ps[j].is_inside(ps[i], cutoff[0], L) ){
//          neigh[0][i]++;
//        }
      }
    }
    ss.str(std::string()); ss << s_coordnum << tag << ".dat"; fout.open(ss.str(), ios::app);
    for(auto i=0;i<N;i++) fout << timestep << " " << i << " " << neigh[0][i] << endl;
    fout.close();
    ss.str(std::string()); ss << s_coordnum << tag << ".ave"; fout.open(ss.str(), ios::app);
    fout << timestep << " " << neigh[0].mean() << " " << neigh[0].std()/sqrt(N) << endl;
    fout.close();
  }

  void compute_boo_boc() { //---- Bond Orientational Parameter -------//
    ntype fval, rijSq, Q,Q2, ql_factor = sqrt(4*M_PI/l_deg), re,im;
    complex<ntype> Yval;
    vec rij;
    int i,j,k,a,m,u;
    if(debug) cout << "\n*** BOO,BOC computation STARTED ***\n";
    //---- Reset counters and lists ----//
    for(u=0;u<Nshells;u++) bond_list[u].clear();
    for(i=0;i<N;i++){
        for(u=0;u<Nshells;u++){
          ps[i].neigh_list[u].clear();
          ps[i].rij_list[u].clear();
          ps[i].rijSq_list[u].clear();
          neigh[u].set(i, 0.0);
        }
        ql.set(i, 0.0);
        ql_dot.set(i, 0.0);
        Ql_dot.set(i, 0.0);
        for(a=0;a<l_deg; a++) qlm[a].set(i, 0.0);
    }
    if(debug) cout << " * Reset counters and lists DONE\n";
    //---- Build neighbour list (and save rij vectors) ----//
    for(i=0;i<N;i++){
      for(j=i+1;j<N;j++){
        rij = ps[j].r - ps[i].r;
        rij = ps[i].mic(rij, L); // minimum image convention
        rijSq = rij.sq();
        for(u=0;u<Nshells;u++){
          if(rijSq <= cutoffSq[u]){
            bond_list[u].push_back( ij2int(i,j,N) );
            
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
    if(debug) cout << " * Build neighbour lists DONE\n";
    //---- Compute Y_lm's and neighs ----//
    for(i=0;i<N;i++){
      for(u=0;u<Nshells;u++){
        for(k=0;k<ps[i].neigh_list[u].size();k++){
          j = ps[i].neigh_list[u][k];
          if(j>i) continue; // avoid double counting!
          rij = ps[i].rij_list[u][k];
          rijSq = ps[i].rijSq_list[u][k];
//          fval = fcut( rijSq/cutoffSq[u], p1half, p2half );
          if(u==0){ //Ider uses a step function for u==0?
            fval = ( rijSq <= cutoffSq[u] ? 1.0 : 0.0 );
          } else {
            fval = fcut( rijSq/cutoffSq[u], p1half, p2half );
          }
          neigh[u][i] += fval;
          neigh[u][j] += fval;
          if(u!=0) continue; // Y(l,m) are computed on first shell neighbours only
          for(a=0;a<l_deg; a++){
              m=a-l; // -l <= m <= l
              Yval = Y(l,m, rij, sqrt(rijSq) );
              qlm[a][i] += fval*Yval;
              if(l_is_odd) Yval=-Yval; // Y(-r) = (-1)^l * Y(r)
              qlm[a][j] += fval*Yval;
          }
        }
      }
    }
    if(debug) cout << " * Compute ql(m,i) and neigh(i) DONE (but ql(m,i) has yet to be normalized)\n";
    //---- Compute ql~qlm(i)*qlm(i) ----//
    Q = 0.0; //  <ql>
    Q2 = 0.0; // <ql^2>
    ss.str(std::string()); ss << s_bondorient << ".l" << l << tag << ".dat"; fout.open(ss.str(), ios::app);
    for(i=0;i<N;i++){
        for(a=0;a<l_deg; a++) {
           if(neigh[0][i]>0) qlm[a][i] /= neigh[0][i]; // qlm(i) = <Ylm>(i) completed
           re=real(qlm[a][i]);
           im=imag(qlm[a][i]);
           ql[i] += re*re + im*im;
        }
        ql[i] = ql_factor * sqrt(ql[i]);
        fout << timestep << " " << i << " " << ql[i] << endl;
        Q += ql[i] / N;
        Q2 += ql[i]*ql[i] / N;
    }
    fout.close();
    ss.str(std::string()); ss << s_bondorient << ".l" << l << tag << ".ave"; fout.open(ss.str(), ios::app);
    fout << timestep << " " << Q << " " << sqrt((Q2-Q*Q)/N) << endl;
    fout.close();
    if(debug) cout << " * Compute ql(i) (BOO) DONE\n";
    //---- Compute Cl_ij~qlm(i)*qlm(j) (I take the real part) and ql_dot(i) (BOC) ----//
    Cl_ij.resize( bond_list[1].size() );
    for(k=0;k<bond_list[1].size();k++){ // Cij and BOC are computed on 2nd shell neighbours
      Cl_ij.set(k, 0.0);
      i = int2i( bond_list[1][k], N );
      j = int2j( bond_list[1][k], N );
      for(a=0;a<l_deg; a++) {
         Cl_ij[k] += ( real(qlm[a][i])*real(qlm[a][j]) + imag(qlm[a][i])*imag(qlm[a][j]) );
      }
      Cl_ij[k] /= (ql[i]*ql[j]); // Cl_ij[k] done
      a = indexOf<int>( ps[i].neigh_list[1], j ); // find index of j in i's neighbour list
      rijSq = ps[i].rijSq_list[1][a]; // and use it to recover the radius
      fval = fcut( rijSq/cutoffSq[1], p1half, p2half );
      ql_dot[i] += fval*Cl_ij[k] / neigh[1][i];
      ql_dot[j] += fval*Cl_ij[k] / neigh[1][j];
    }
    if(debug) cout << " * Compute C_l(i,j) and ql_dot(i) (BOC) DONE\n";
    //---- Compute global average of ql_dot(i) ----//
    ss.str(std::string()); ss << s_bondcorr << ".l" << l << tag << ".dat"; fout.open(ss.str(), ios::app);
    Q = 0.0; //  <ql_dot>
    Q2 = 0.0; // <ql_dot^2>
    for(i=0;i<N;i++){
/*      for(k=0;k<ps[i].neigh_list[1].size();k++){
        j = ps[i].neigh_list[1][k];
        rijSq = ps[i].rijSq_list[1][k];
        fval = fcut( rijSq/cutoffSq[1], p1half, p2half );
        ql_dot[i] += fval*Cl_ij[ a++ ] / neigh[1][i];
      }*/
      fout << timestep << " " << i << " " << ql_dot[i] << endl;
      Q += ql_dot[i] / N;
      Q2 += ql_dot[i]*ql_dot[i] / N;
    }
    fout.close();
    ss.str(std::string()); ss << s_bondcorr << ".l" << l << tag << ".ave"; fout.open(ss.str(), ios::app);
    fout << timestep << " " << Q << " " << sqrt((Q2-Q*Q)/N) << endl;
    fout.close();
    if(debug) cout << " * Compute global average of BOC DONE\n";
    //---- Compute a Locally Confined version of the global BOC ~ average <ql_dot(i)> around i ----//
    ss.str(std::string()); ss << s_bondcorr << ".l" << l << tag << ".local.ave"; fout.open(ss.str(), ios::app);
    for(i=0;i<N;i++){
      for(k=0;k<ps[i].neigh_list[2].size();k++){
        j = ps[i].neigh_list[2][k];
        rijSq = ps[i].rijSq_list[2][k];
        fval = fcut( rijSq/cutoffSq[2], p1half, p2half );
        Ql_dot[i] += fval*ql_dot[j] / neigh[2][i];
      }
      fout << timestep << " " << i << " " << Ql_dot[i] << endl;
    }
    fout.close();
    if(debug) cout << " * Compute local average of BOC DONE\n";
    ss.str(std::string()); ss << s_bondcorr << ".l" << l << tag << ".xyz"; fout.open(ss.str(), ios::app);
    fout << N << endl << "-1: crystal, -2: noncrystal. Timestep = " << timestep << endl;
    m=0; // number of crystalline particles
    for(i=0;i<N;i++){
      if(ql_dot[i] > qldot_th){
        ps[i].label =  -1;
        m++;
      } else{
        ps[i].label =  -2;
      }
      ps[i].write_xyz(fout);
    }
    fout.close();
    ss.str(std::string()); ss << s_nxtal << ".l" << l << tag << ".dat"; fout.open(ss.str(), ios::app);
    fout << timestep << " " << m << endl;
    fout.close();
    if(debug) cout << "*** BOO,BOC computation for timestep " << timestep << " ENDED ***\n\n";
  }

  void compute_msd(int frameidx) { //---- Mean Squared Displacement -------//
    int i, idx, idx_old, dframe, nperiod;
    vec dr;
    ntype dr2, r2t0;
    if(debug) cout << "*** MSD computation for timestep " << timestep << " STARTED ***\n";
    nperiod = frameidx / periodIdx;
    dframe = frameidx % periodIdx;
    ss.str(std::string()); ss << s_msd << tag << ".traj"; fout.open(ss.str(), ios::app);
    if(frameidx==0){
      for(i=0;i<periodIdx-1;i++){ // reset averages
        r2[i]=0.0;
        r4[i]=0.0;
        rs[ (N+1)*i + N ] << 0.0,0.0,0.0;
        r2CM[i]=0.0;
        fout << (i+1)*dtframe << endl;
      }
      rs[ (N+1)*(periodIdx-1) + N ] << 0.0,0.0,0.0;
      fout << endl; // end of the first block (delta timesteps)
    }
    if(dframe==0){
      r2t0 = 0.0; // average r^2 over atoms at fixed t,t0
      for(i=0;i<N;i++){
        rs[ (N+1)*0 + i] = ps[i].r; //store initial positions for this sample trajectory
        rs[ (N+1)*0 + N ] += ( ps[i].r * invN ); // store initial center of mass for '''
      }
    }else{
      for(i=0;i<N;i++)
      {
        idx = (N+1)*dframe + i;
        idx_old = idx - (N+1);       // previous frame
        dr = ps[i].r - rs[idx_old];  // displacement since previous frame
        dr = dr.mic(L);              // correct for periodic boundary conditions
        rs[idx] = rs[idx_old] + dr;  // add correction
//        if(i==0) cout << timestep << " " << rs[idx][0]-rs[(N+1)*0+i][0] << endl;
        rs[ (N+1)*dframe + N ] += ( rs[idx] * invN ); // compute center of mass
        
        idx_old = (N+1)*0 + i;       // t0 for this sample trajectory
        dr = rs[idx] - rs[idx_old];  // displacement r(t)-r(t0)
        dr2 = dr.sq();
        r2t0 += dr2 * invN;          // sample trajectory
        r2[dframe-1] += dr2;         // average over t0 and atoms
        r4[dframe-1] += dr2*dr2;
      }
      fout << r2t0 << endl; // output sample trajectory
      if(dframe == periodIdx-1) fout << endl; // end of the block (sample trajectory)
      // after evaluating all particles, do the same for the center of mass
      dr = rs[ (N+1)*dframe + N ] - rs[ (N+1)*0 + N ];
      dr2 = dr.sq();
      r2CM[dframe-1] += dr2;
    }
    fout.close();
    if(debug) cout << "*** MSD computation for timestep " << timestep << " ENDED ***\n\n";
  }

  void print_msd() {
    ss.str(std::string()); ss << s_msd << tag << ".ave"; fout.open(ss.str(), ios::app);
    for(auto i=0; i<periodIdx-1; i++) {
      r2[i] /= (Nperiods*N);
      r4[i] /= (Nperiods*N);
      r2CM[i] /= Nperiods;
      fout << (i+1)*dtframe << " " << r2[i] << " " << sqrt( (r4[i] - r2[i]*r2[i])/(Nperiods*N) ) << " " << r2CM[i] << endl;
    }
    fout.close();
    if(debug) cout << "MSD printed to file\n";
    ss.str(std::string()); ss << s_ngp << tag << ".ave"; fout.open(ss.str(), ios::app);
    for(auto i=0; i<periodIdx-1; i++) {
      fout << (i+1)*dtframe << " " << 0.6*(r4[i]/(r2[i]*r2[i]))-1.0 << endl;
    }
    fout.close();
    if(debug) cout << "NGP printed to file\n";
  }
  
  void compute_rdf(int frameidx) { //---- Radial Distribution Function -------//
    int i,j, k;
    vec rij;
    ntype r;
    for(k=0; k<rdf_nbins; k++)
      rdf[k] = 0.0;
    if(debug) cout << "*** RDF computation for timestep " << timestep << " STARTED ***\n";
    for(i=0;i<N;i++){
      for(j=i+1;j<N;j++){
         rij = ( ps[i].r - ps[j].r ).mic(L);
         r = rij.norm();
         k = int(floor( r/rdf_binw));
         if(k<rdf_nbins){
           rdf[k] += 1.0;
         }
      }
    }
    ss.str(std::string()); ss << s_rdf << tag << ".traj"; fout.open(ss.str(), ios::app);
    fout << endl; // start new block
    for(k=0; k<rdf_nbins; k++){
      rdf[k] /= rdf_norm[k];
      fout << rdf[k] << endl;
      rdf_ave[k] += rdf[k];
      rdf2_ave[k] += rdf[k]*rdf[k];
    }
    fout.close();
    if(frameidx == (nframes-1)){
        ss.str(std::string()); ss << s_rdf << tag << ".ave"; fout.open(ss.str(), ios::app);
        for(k=0; k<rdf_nbins; k++){
          rdf_ave[k] /= (frameidx+1);
          rdf2_ave[k] /= (frameidx+1);
          fout << rdf_bins[k] << " " << rdf_ave[k] << " " << sqrt( (rdf2_ave[k]-rdf_ave[k]*rdf_ave[k])/(frameidx+1 -1) ) << endl;
        }
        fout.close();
        if(debug) cout << "Average RDF printed to file\n";
    }
    if(debug) cout << "*** RDF computation for timestep " << timestep << " ENDED ***\n\n";
    return;
  }

  ntype fcut(ntype x, int pow1, int pow2) {
    ntype x1=1.0, x2=1.0;
    for(auto i=0; i<pow2; i++) {
      x2 *= x;
      if(i==pow1-1) x1=x2;
    }
    return (1.0-x1) / (1.0-x2); // x1 = x^p1, x2 = x^p2
  }

//------------ Input reading and interaction --------------------//

void print_usage(char argv0[])
{
  fprintf(stderr, "Usage: %s [-d -h -v] [-bo] [-cn] [-in] [-l] [-L] [-msd] [-p/--period] [-rcut1] [-rcut2] [-rcut3] [-rdf] [-tag]\n", argv0);
}
void print_summary()
{
  const char *a_xyz = s_xyz.c_str();
  fprintf(stderr, "\nThis program computes some statistical quantities over a MD trajectory.\n");
  fprintf(stderr, "\n -h/--help \t Print this summary.\n -bo \t Compute the bond order orientation (BOO) and correlation (BOC) parameters. Angular momentum is defined by the option -l.\n -cn \t Compute the coordination number.\n -d/--debug \t Open in debug mode.\n -in/--input \t Input .xyz trajectory file [default: %s].\n -l \t Angular momentum for the computed bond order parameters [default %d].\n -L \t Lx,Ly,Lz sizes of the orthorombic supercell [default %.2f %.2f %.2f].\n -msd \t Compute Mean Squared Displacement.\n -p/--period \t Average over initial time t0 every 'period' (in timesteps units) when computing MSD. If negative, don't. [default %d].\n -rcut1 \t Cutoff radius for cutoff functions in 1st shell [default %.2f].\n -rcut2 \t Cutoff radius for cutoff functions in 2nd shell [default %.2f].\n -rcut3 \t Cutoff radius for cutoff functions in 3rd shell [default %.2f].\n -rdf \t Compute Radial Distribution Function using the given number of bins.\n -tag \t Add this text tag inside output files' name [default none].\n -v/--verbose \t Print a lot of outputs during execution.\n\n", a_xyz, l, L[0],L[1],L[2], period, cutoff[0],cutoff[1],cutoff[2]);
}

void args(int argc, char** argv)
{
  int i=1;
  // strcmp returns 0 if the strings are identical
  if (argc == 2 && ( !strcmp(argv[i], "-h") || !strcmp(argv[i],"--help") ))
    {
      print_usage(argv[0]);
      print_summary();
      exit(-1);
    }
  else if (argc >1 )
    {
      while (i<argc)
	{
	  if ( !strcmp(argv[i], "-d") || !strcmp(argv[i],"--debug") )
	    debug = true;
	  else if ( !strcmp(argv[i], "-v") || !strcmp(argv[i],"--verbose") )
	    verbose = true;
	  else if ( !strcmp(argv[i], "-bo") )
	      c_bondorient = true;
	  else if ( !strcmp(argv[i], "-cn") )
	      c_coordnum = true;
	  else if ( !strcmp(argv[i], "-in") || !strcmp(argv[i],"--input") )
	    {
	      i++;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: '-in' must be followed by file name!\n");
		  exit(-1);
		}
	      s_xyz = string(argv[i]);
	    }
	  else if ( !strcmp(argv[i], "-l") )
	    {
	      i++;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: '-l' must be followed by angular momentum value!\n");
		  exit(-1);
		}
	      l = atoi(argv[i]);
	    }
	  else if ( !strcmp(argv[i], "-L") )
	    {
	        for(int j=0;j<3;j++)
	            {
		      i++;
		      if (i == argc)
			{
			  fprintf(stderr, "ERROR: '-L' must be followed by the 3 box sizes!\n");
			  exit(-1);
			}
		      L[j] = atof(argv[i]);
		    }
	    }
	  else if ( !strcmp(argv[i], "-msd") )
	      c_msd = true;
	  else if ( !strcmp(argv[i], "-p") || !strcmp(argv[i], "--period") )
	    {
	      i++;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: '-p/--period' must be followed by a positive integer value!\n");
		  exit(-1);
		}
	      period = atoi(argv[i]);
	    }
	  else if ( !strcmp(argv[i], "-rcut1") )
	    {
	      i++;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: '-rcut1' must be followed by cutoff value!\n");
		  exit(-1);
		}
	      cutoff[0] = atof(argv[i]);
	    }
	  else if ( !strcmp(argv[i], "-rcut2") )
	    {
	      i++;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: '-rcut2' must be followed by cutoff value!\n");
		  exit(-1);
		}
	      cutoff[1] = atof(argv[i]);
	    }
	  else if ( !strcmp(argv[i], "-rcut3") )
	    {
	      i++;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: '-rcut3' must be followed by cutoff value!\n");
		  exit(-1);
		}
	      cutoff[2] = atof(argv[i]);
	    }
	  else if ( !strcmp(argv[i], "-rdf") )
	    {
	      c_rdf = true;
	      i++;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: '-rdf' must be followed by number of bins!\n");
		  exit(-1);
		}
	      rdf_nbins = atof(argv[i]);
	      if(rdf_nbins < 2){
	        fprintf(stderr, "ERROR: too few bins for RDF!\n");
	        exit(1);
	      }
	    }
	  else if ( !strcmp(argv[i], "-tag") )
	    {
	      i++;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: '-tag' must be followed by some text!\n");
		  exit(-1);
		}
	      tag = argv[i];
	      tag.insert(0, 1, '.'); // add a dot to the beginning of the tag
	    }
	  else
	    {
	      fprintf(stderr, "ERROR: Invalid argumet!\n");
	      exit(-1);
	    }
      	  i++;
	}
    }
}
};
#endif
