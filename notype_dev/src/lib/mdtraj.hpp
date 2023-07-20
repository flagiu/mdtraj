#ifndef _MDTRAJ_H_
#define _MDTRAJ_H_

#include<cstdlib>
#include<cstring>
#include<complex>
#include<cmath>
#include "utility.hpp"
#include "vecflex.hpp"
#include "particle.hpp"
#include "mymatrix.hpp"
#include "Ycomplex.hpp"
using namespace std;

enum class FileType {
  XYZ, XYZ_CP2K, CONTCAR, XDATCAR, XDATCARV, ALPHANES, ALPHANES9, JMD
};

template<class ntype, class ptype>
class Trajectory {
public:
  using vec=myvec<ntype,3>;
  using mat=mymatrix<ntype,3,3>;
  int N; // number of particles
  vec L; // dimensions of an orthorombic simulation box ( ??? centered in 0: -Lx/2 < x < Lx/2 )
  mat box, boxInv; // most general simulation box
  bool remove_rot_dof; // remove the 3 rotational degrees of freedom from the box?
  ntype V, mdens, ndens; // volume, mass density, nuerical density
  vector<ptype> ps, ps_new; // vector of particles
  int nframes, timestep, period, l, rdf_nbins, adf_nbins, altbc_nbins;
  bool c_coordnum, c_bondorient, c_msd, c_rdf, c_adf, c_rmin, c_altbc; // compute or not
  string s_in, s_out, tag, s_box, s_ndens, s_coordnum, s_bondorient, s_bondcorr, s_nxtal, s_msd, s_ngp, s_rdf, s_adf, s_rmin, s_tbc, s_altbc; // for file naming
  static const int Nshells=3;
  ntype cutoff[Nshells], altbc_rmin, altbc_angle;
  vecflex<ntype> neigh[Nshells], ql, Cl_ij, ql_dot, Ql_dot;
  vector< vecflex< complex<ntype> > > qlm; // a collection of l_deg vectors of local average qlm=<Ylm> with a cutoff function
  vector<int> bond_list[Nshells]; // records all bonds encoded into an integer through ij2int()
  vecflex<ntype> r2,r4, r2CM; // for MSD
  vecflex<ntype> rdf_bins, rdf_norm, rdf, rdf_ave, rdf2_ave; // for RDF
  vecflex<ntype> adf_bins, adf, adf_ave, adf2_ave; // for ADF
  vecflex<ntype> altbc_bins, altbc, altbc_ave, altbc2_ave; // for ALTBC
  vector<vec> rs;
  bool msdAverageOverTime0, out_box, out_xyz, out_alphanes;
  bool debug, verbose;

private:
  bool l_is_odd;
  int nlines, t0frame, dtframe, l_deg, periodIdx, Nperiods, maxshell;
  int p1half, p2half, p1, p2; // parameters for fcut (only p1half is free)
  fstream fin, fout;
  stringstream ss;
  ntype cutoffSq[Nshells], invN, qldot_th, rdf_binw, adf_binw, altbc_binw, altbc_cos;
  FileType filetype;

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

//------------ Input reading and interaction --------------------//
// implemented in other ../args.cpp

  void print_usage(char argv0[]);
  void print_summary();
  void args(int argc, char** argv);

  void print_state(){
    cout << "Summary of parameters:\n";
    cout << " debug = \t " << debug << endl;
    cout << " verbose = \t " << debug << endl;
    cout << " c_coordnum = \t " << c_coordnum << endl;
    cout << " c_bondorient = \t " << c_bondorient << endl;
    cout << " c_msd = \t " << c_msd << endl;
    cout << " c_rdf = \t " << c_rdf << endl;
    cout << " c_adf = \t " << c_adf << endl;
    cout << " c_rmin = \t " << c_rmin << endl;
    cout << " c_altbc = \t " << c_altbc << endl;
    cout << " angular momentum for ql: l = \t " << l << endl;
    cout << " qldot threshold = \t " << qldot_th << endl;
    cout << " box (a|b|c) = \t "; box.show();
    cout << " total volume V = \t " << V << endl;
    cout << " box inverse = \t "; boxInv.show();
    cout << " L (length of each box vector) = \t "; L.show();
    cout << " p1 = \t " << p1 << endl;
    cout << " p2 = \t " << p2 << endl;
    cout << " MSD period = \t " << period << endl;
    cout << " rdf_nbins = \t " << rdf_nbins << endl;
    cout << " adf_nbins = \t " << adf_nbins << endl;
    cout << " altbc_nbins = \t " << altbc_nbins << endl;
    cout << " altbc_rmin = \t " << altbc_rmin << endl;
    cout << " altbc_angle = \t " << altbc_angle << endl;
    cout << " remove rotational degrees of freedom = \t " << remove_rot_dof << endl;
    cout << " out_box = \t " << out_box << endl;
    cout << " out_xyz = \t " << out_xyz << endl;
    cout << " out_alphanes = \t " << out_alphanes << endl;
    cout << " tag = \t " << tag << endl;
    for(auto  i=0;i<Nshells;i++) cout << " rcut" << i << " = \t " << cutoff[i] << endl;
    cout << " s_in = \t " << s_in << endl;
    cout << endl;
  }

  void init(int argc, char** argv){
    // Default values for primitive parameters:
    debug = false;
    verbose = false;
    c_coordnum = false;
    c_bondorient = false;
    c_msd = false;
    c_rdf = false;
    c_adf = false;
    c_rmin = false;
    c_altbc = false;
    filetype=FileType::XYZ;
    s_ndens="ndens";
    s_box="box";
    s_coordnum="coordnum";
    s_bondorient="boo";
    s_bondcorr="boc";
    s_nxtal="nc";
    s_msd="msd";
    s_ngp="ngp";
    s_rdf="rdf";
    s_adf="adf";
    s_rmin="rmin";
    s_altbc="altbc";
    tag="";
    s_out="traj";
    period = -1; // default: don't average over t0 for MSD
    remove_rot_dof = false;
    out_box = false;
    out_xyz = false;
    out_alphanes = false;
    // default: nonsense box
    L << 0.0, 0.0, 0.0;
    for(auto i=0;i<3;i++) {
      box[i] << 0.0, 0.0, 0.0;
      boxInv[i] << 0.0, 0.0, 0.0;
    }
    V=0.0;
    cutoff[0] = 3.75; // 3.6 in glass, 3.75-3.89 in xtal
    cutoff[1] = 5.15;
    cutoff[2] = 8.8;
    l = 4;
    p1half=6;
    qldot_th = 0.65;
    adf_nbins=0;
    rdf_nbins=0;
    altbc_nbins=0;
    altbc_rmin=0.0;
    altbc_angle=-1.0;
    //-------- Update parameters with input arguments: -----------//
    args(argc, argv);
    // Compute non-primitive parameters:
    for(auto  i=0;i<Nshells;i++) cutoffSq[i] = cutoff[i]*cutoff[i];
    l_is_odd = (l%2!=0);
    l_deg = 2*l+1;
    p2half = 2*p1half;
    p1 = 2*p1half;
    p2 = 2*p2half;
    if(c_bondorient)             maxshell=Nshells; // init all neigh shells
    else if(c_coordnum || c_adf || c_rmin || c_altbc) maxshell=1;       // init only first neigh shell
    else                         maxshell=0;       // do not init any
    // Print a recap:
    if(debug) { cout << "State after reading args():\n"; print_state(); }
    if(l==0) cout << "[WARNING: bond orientation order parameter with l=0 is always equal to 1.0]\n";
  }

  void set_box_from_L() {
    box.set_diag(L); // diagonal box
    boxInv = box.inverse();
    compute_volume();
  }

  void set_L_from_box() {
    boxInv=box.inverse();
    L[0] = box.T()[0].norm(); // in general, L[i] is the length of the i-th box vector. Is it useful? idk
    L[1] = box.T()[1].norm();
    L[2] = box.T()[2].norm();
    compute_volume();
  }

  void compute_volume() {
    V = box.det(); //determinant
    if(V<0.) {
      cout << "[ Warning: det(box)="<<V<<" follows left-hand rule. ]\n";
      V=-V;
    }
    else if (V==0.) {
      box.show();
      cout << "[ Error: det(box)=0.0 not supported. ]\n";
      exit(1);
    }
  }

  void read_frame(fstream &i, bool resetN);
  void removeRotDof();
  void read_contcar_frame(fstream &i, bool resetN);
  void read_xdatcar_frame(fstream &i, bool resetN, bool constantBox);
  void read_xyz_frame(fstream &i, bool resetN);
  void read_xyz_cp2k_frame(fstream &i, bool resetN);
  void read_alphanes_frame(fstream &i, bool resetN);
  void read_alphanes9_frame(fstream &i, bool resetN);
  void read_jmd_frame(fstream &i, bool resetN);

  //------- COMPUTE things ---------------//

  void run()
  {
    PrintProgress printProgress;
    if(verbose || debug) {cout << "Begin of run():\n"; print_state();}
    nlines = getLineCount(s_in);
    if(debug) cout << "Read " << nlines << " lines in file " << s_in << ". Opening again for reading trajectory.\n";
    fin.open(s_in, ios::in);

    if(filetype==FileType::CONTCAR || filetype==FileType::ALPHANES || filetype==FileType::ALPHANES9) {
      timestep=-1; dtframe = 1;
    } // set manual time for CONTCAR, ALPHANES, ALPHANES9 file format
    read_frame(fin, true);
    t0frame = timestep;
    if(debug) cout << "Read first frame. Set N = " << N << " (assumed to beconstant), t0frame = " << t0frame << ".\n";
    if(debug) cout << "Deduced nframes = " << nframes << ".\n";

    read_frame(fin, false);
    if(filetype!=FileType::CONTCAR && filetype!=FileType::ALPHANES && filetype!=FileType::ALPHANES) dtframe = timestep - t0frame;
    if(debug) cout << "Read second frame. Set dtframe = " << dtframe << " (assumed to be constant).\n";
    init_computations();
    if(debug) cout << "Initialized arrays for computations.\n";
    fin.close();

    // Restart reading
    fin.open(s_in, ios::in);
    printProgress.init( nframes );
    if(filetype==FileType::CONTCAR || filetype==FileType::ALPHANES) {timestep=-1; } // set manual time
    string junk_line;
    if(filetype==FileType::XDATCARV) { for(int i=0;i<7;i++) getline(fin, junk_line); } // skip first 7 lines (so that you can use resetN=false)
    for(int i=0; i<nframes; i++)
    {
      read_frame(fin, false);
      if(N != ps.size()) { cout << "[Warning: N has changed]\n"; exit(1);}
      if( (timestep - t0frame)%dtframe != 0) {
        cout << "[Warning: timestep interval has changed]\n";
        cout << "[t0frame = "<<t0frame<<", dtframe = "<<dtframe<<", timestep = "<<timestep<<"]\n";
        exit(1);
      }
      do_computations_and_output(i);
      printProgress.update( i+1 );
    }
    printProgress.end();
    fin.close();
    if(debug) cout << "Closed input file.\n";
    if(c_msd) print_msd();
    if(debug || verbose) cout << "\nExecution completed.\n\n";
  }

  void init_computations() {
    init_density();
    if(out_box) init_box();
    if(out_xyz) init_out_xyz();
    if(out_alphanes) init_out_alphanes();
    if(maxshell>0) init_neigh();
    if(c_coordnum) init_coordnum();
    if(c_bondorient) init_bondorient();
    if(c_msd) init_msd();
    if(c_rdf) init_rdf();
    if(c_adf) init_adf();
    if(c_rmin) init_rmin();
    if(c_altbc) init_altbc();
  }

  void do_computations_and_output(int i) {
    compute_density();
    if(out_box) print_box();
    if(out_xyz) print_out_xyz();
    if(out_alphanes) print_out_alphanes();

    if(maxshell>0) build_neigh();
    if(c_coordnum) compute_coordnum();
    if(c_bondorient) compute_bondorient();
    if(c_msd) compute_msd(i);
    if(c_rdf) compute_rdf(i);
    if(c_adf) compute_adf(i);
    if(c_rmin) compute_rmin();
    if(c_altbc) compute_altbc(i);
  }

  void init_density() {
    compute_volume();
    ndens = N/V;
    ss.str(std::string()); ss << s_ndens << tag << ".dat"; fout.open(ss.str(), ios::out);
    fout << "#Timestep, Number density\n";
    fout.close();
  }
  void compute_density() {
    compute_volume();
    ndens = N/V;
    ss.str(std::string()); ss << s_ndens << tag << ".dat"; fout.open(ss.str(), ios::app);
    fout << timestep << " " << ndens << endl;
    fout.close();
  }
  
//-------------Trajectory output, implemented in ../output.cpp -----------------//
  void init_box();
  void print_box();
  void init_out_xyz();
  void print_out_xyz();
  void init_out_alphanes();
  void print_out_alphanes();

  //-------------Trajectory analysis, implemented in ../statics.cpp and ../dynamics.cpp -----------------//
  void init_neigh();
  void build_neigh();
  void init_coordnum();
  void compute_coordnum();
  void init_bondorient();
  void compute_bondorient();
  void init_rdf();
  void compute_rdf(int frameidx);
  void init_adf();
  void compute_adf(int frameidx);
  void init_rmin();
  void compute_rmin();
  void init_altbc();
  void compute_altbc(int frameidx);
  void init_msd();
  void compute_msd(int frameidx);
  void print_msd();

  ntype fcut(ntype x, int pow1, int pow2) {
    ntype x1=1.0, x2=1.0;
    for(auto i=0; i<pow2; i++) {
      x2 *= x;
      if(i==pow1-1) x1=x2;
    }
    return (1.0-x1) / (1.0-x2); // x1 = x^p1, x2 = x^p2
  }

};

#endif
