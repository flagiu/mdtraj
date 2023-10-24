#ifndef _MDTRAJ_H_
#define _MDTRAJ_H_

#include<cstdlib>
#include<cstring>
#include<complex>
#include<cmath>
#include "lib/utility.hpp"
#include "lib/vecflex.hpp"
#include "lib/particle.hpp"
#include "lib/mymatrix.hpp"
#include "lib/Ycomplex.hpp"
using namespace std;

const string root_path="/home/flavio/programmi/mdtraj/mixture_dev";
#define MAX_N_TYPES 5
enum class FileType {
  XYZ, XYZ_CP2K, CONTCAR, XDATCAR, XDATCARV, ALPHANES, ALPHANES9, JMD, LAMMPSTRJ, YUHAN
};

template<class ntype, class ptype>
class Trajectory {
public:
  using vec=myvec<ntype,3>;
  using mat=mymatrix<ntype,3,3>;
  int N, nTypes, nTypePairs, Nt[MAX_N_TYPES]; // number of particles, number of types, num of particles for each type
  string type_names[MAX_N_TYPES]; // name for each type
  vec L; // dimensions of an orthorombic simulation box ( ??? centered in 0: -Lx/2 < x < Lx/2 )
  mat box, boxInv; // most general simulation box
  bool remove_rot_dof; // remove the 3 rotational degrees of freedom from the box?
  bool pbc_out; // print output with PBC?
  ntype V, mdens, ndens; // volume, mass density, nuerical density
  vector<ptype> ps, ps_new; // vector of particles
  int nframes, timestep;
  bool c_coordnum, c_bondorient, c_msd, c_rdf, c_adf, c_rmin, c_rmax, c_altbc, c_sq, c_sqt, c_edq; // compute or not
  string s_in, s_out, s_rcut, tag, s_atom_label, s_box, s_ndens, s_coordnum, s_bondorient, s_bondcorr, s_nxtal, s_msd, s_ngp, s_rdf, s_adf, s_rmin, s_tbc, s_altbc, s_sq, s_sqt, s_log, s_rmax, s_edq; // for file naming
  bool out_box, out_xyz, out_alphanes;
  bool debug, verbose;
  //
  int maxshell; // <= MAX_NSHELL
  vecflex<ntype> cutoff[MAX_NSHELL], cutoffSq[MAX_NSHELL];
  Neigh_and_Bond_list<ntype,ptype> *n_b_list;
  //
  int l;
  ntype qldot_th;
  Bond_Parameters<ntype,ptype> *bond_parameters;
  ED_Bond_Parameter<ntype,ptype> *ed_q_calculator;
  //
  ntype rdf_binw;
  RDF_Calculator<ntype,ptype> *rdf_calculator;
  //
  int qmodmin,qmodmax,qmodstep;
  SQ_Calculator<ntype,ptype> *sq_calculator;
  //
  ntype adf_binw;
  ADF_Calculator<ntype,ptype> *adf_calculator;
  //
  int period;
  bool msdAverageOverTime0;
  MSD_Calculator<ntype,ptype> *msd_calculator;
  //
  ntype altbc_rmin, altbc_binw, altbc_angle_th;
  ALTBC_Calculator<ntype,ptype> *altbc_calculator;
  //
  SQT_Calculator<ntype,ptype> *sqt_calculator;

private:
  bool timings;
  int nlines, t0frame, dtframe, nskip0, nskip1, nframes_original;
  int p1half, p2half, p1, p2; // parameters for fcut (only p1half is free)
  float fskip0, fskip1;
  fstream fin, fout;
  stringstream ss;
  FileType filetype;
  Timer timer, sq_timer;

  int types2int(int ti, int tj, int nTypes){ // map type pairs (ti,tj) in 0,1,...,nTypes-1 to integer index 0,1,...,nTypePairs
    if (ti>tj) return types2int(tj,ti,nTypes); // map to ti<=tj
    if(ti==tj) return ti*nTypes - int(ti*(ti-1)/2);
    else return 1 + types2int(ti,tj-1,nTypes);
  }
  void int2types(int t, int nTypes, int *t1, int *t2){
    int x,low,high;
    for(x=0; x<nTypes; x++){
      low = x * nTypes - x*(x-1)/2;
      high = (x+1) * nTypes - x*(x-1)/2 - 1;
      if(t>=low && t<=high){
        *t1 = x;
        *t2 = t-low;
      }
    }
    return;
  }

public:
  Trajectory(){}
  ~Trajectory(){}

//------------ Input reading and interaction --------------------//

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
    cout << " c_rmax = \t " << c_rmax << endl;
    cout << " c_altbc = \t " << c_altbc << endl;
    cout << " c_sq = \t " << c_sq << endl;
    cout << " c_sqt = \t " << c_sqt << endl;
    cout << " c_edq = \t " << c_edq << endl;
    cout << " angular momentum for ql: l = \t " << l << endl;
    cout << " qldot threshold = \t " << qldot_th << endl;
    cout << " box (a|b|c) = \t "; box.show();
    cout << " total volume V = \t " << V << endl;
    cout << " box inverse = \t "; boxInv.show();
    cout << " L (length of each box vector) = \t "; L.show();
    cout << " p1 = \t " << p1 << endl;
    cout << " p2 = \t " << p2 << endl;
    cout << " MSD period = \t " << period << endl;
    cout << " rdf_binw = \t " << rdf_binw << endl;
    cout << " adf_binw = \t " << adf_binw << endl;
    cout << " q_mod_min,q_modmax,q_mod_step = \t " << qmodmin << ", " << qmodmax << ", "<< qmodstep << endl;
    cout << " remove rotational degrees of freedom = \t " << remove_rot_dof << endl;
    cout << " out_box = \t " << out_box << endl;
    cout << " out_xyz = \t " << out_xyz << endl;
    cout << " out_alphanes = \t " << out_alphanes << endl;
    cout << " pbc_out = \t " << pbc_out << endl;
    cout << " fskip_from_beginning = \t " << fskip0 << endl;
    cout << " fskip_from_end = \t " << fskip1 << endl;
    cout << " timings = \t " << timings << endl;
    cout << " tag = \t " << tag << endl;
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
    c_rmax = false;
    c_altbc = false;
    c_sq = false;
    c_sqt = false;
    c_edq = false;
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
    s_rmax="rmax";
    s_altbc="altbc";
    s_sq="sq";
    s_sqt="sqt";
    s_edq="ed_q";
    tag="";
    s_out="traj";
    s_log="log";
    s_atom_label="labels";
    ss.str(std::string()); ss << s_log << tag; fout.open(ss.str(), ios::out);
    fout << "LOG SUMMARY"<<endl;
    fout.close();
    period = -1; // default: don't average over t0 for MSD
    remove_rot_dof = false;
    out_box = false;
    out_xyz = false;
    out_alphanes = false;
    pbc_out = false;
    timings = false;
    // default: nonsense box
    L << 0.0, 0.0, 0.0;
    for(auto i=0;i<3;i++) {
      box[i] << 0.0, 0.0, 0.0;
      boxInv[i] << 0.0, 0.0, 0.0;
    }
    V=0.0;
    nTypes=nTypePairs=1;
    cutoff[0].resize(1); // 1st shell
    cutoff[0][0] = 3.75; // Antimony: 3.6 in glass, 3.75-3.89 in xtal
    cutoff[1].resize(1); // 2nd shell
    cutoff[1][0] = 5.15;
    cutoff[2].resize(1); // 3rd shell
    cutoff[2][0] = 8.8;
    l = 4;
    p1half=6;
    qldot_th = 0.65;
    rdf_binw=0.0;
    qmodmin=2;
    qmodmax=100;
    qmodstep=1;
    adf_binw=0.0;
    altbc_binw=0.0;
    altbc_rmin=0.0;
    altbc_angle_th=-1.0;
    fskip0=fskip1=0.0;
    //-------- Update parameters with input arguments: -----------//
    args(argc, argv);
    // Compute non-primitive parameters:
    p2half = 2*p1half;
    p1 = 2*p1half;
    p2 = 2*p2half;
    if(c_bondorient)             maxshell=MAX_NSHELL; // init all neigh shells
    else if(c_coordnum || c_adf || c_rmin || c_rmax || c_altbc || c_edq) maxshell=1;       // init only first neigh shell
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
    //L[0] = box.T()[0].norm(); // in general, L[i] is the length of the i-th box vector. Is it useful? idk
    //L[1] = box.T()[1].norm();
    //L[2] = box.T()[2].norm();
    for(int i=0;i<3;i++) L[i] = box[i][i]; // in general, L[i] is the length of the i-th box vector. Is it useful? idk
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

  void read_frame(fstream &i, bool resetN, int frameIdx);
  void removeRotDof();
  void read_contcar_frame(fstream &i, bool resetN);
  void read_xdatcar_frame(fstream &i, bool resetN, bool constantBox);
  void read_xyz_frame(fstream &i, bool resetN);
  void read_xyz_cp2k_frame(fstream &i, bool resetN);
  void read_alphanes_frame(fstream &i, bool resetN);
  void read_alphanes9_frame(fstream &i, bool resetN);
  void read_jmd_frame(fstream &i, bool resetN);
  void read_lammpstrj_frame(fstream &i, bool resetN);
  void read_yuhan_frame(fstream &i, bool resetN, bool isFirstFrame);

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
    read_frame(fin, true, 0);
    t0frame = timestep;
    if(debug) cout << "Read first frame. Set N = " << N << " (assumed to beconstant), t0frame = " << t0frame << ".\n";
    if(debug) cout << "Deduced nframes = " << nframes << ".\n";
    print_types();

    nskip0=int(fskip0*nframes);
    nskip1=int(fskip1*nframes);
    nframes_original = nframes;
    nframes = nframes - nskip0 - nskip1;
    if(nframes<2) { cout << "[ Error: skipped too many frames.\n  Total: "<<nframes_original<<"; Skipped: "<<nskip0<<"+"<<nskip1<<"; Remaining: "<<nframes<<" ]\n\n"; exit(1);}

    read_frame(fin, false, 1);
    if(filetype!=FileType::CONTCAR && filetype!=FileType::ALPHANES && filetype!=FileType::ALPHANES9) dtframe = timestep - t0frame;
    if(debug) cout << "Read second frame. Set dtframe = " << dtframe << " (assumed to be constant).\n";
    init_computations();
    if(debug) cout << "Initialized arrays for computations.\n";
    fin.close();

    // Restart reading
    if(debug || verbose) cout << "#------- MAIN LOOP ------#\n";
    fin.open(s_in, ios::in);
    printProgress.init( nframes, 2000 ); // update % every 2000 ms
    if(filetype==FileType::CONTCAR || filetype==FileType::ALPHANES || filetype==FileType::ALPHANES9) {timestep=-1; } // set manual time
    string junk_line;
    if(filetype==FileType::XDATCARV) { for(int i=0;i<7;i++) getline(fin, junk_line); } // skip first 7 lines (so that you can use resetN=false)
    timer.go();
    for(int i=0; i<nframes_original; i++)
    {
      read_frame(fin, false, i);
      if(N != ps.size()) { cout << "[Warning: N has changed]\n"; exit(1);}
      if( (timestep - t0frame)%dtframe != 0) {
        cout << "[Warning: timestep interval has changed]\n";
        cout << "[t0frame = "<<t0frame<<", dtframe = "<<dtframe<<", timestep = "<<timestep<<"]\n";
        exit(1);
      }
      if(i+1>nskip0 && i<nframes_original-nskip1) do_computations_and_output(i-nskip0);
      printProgress.update( i+1-nskip0, timer.lap() );
    }
    printProgress.end();
    fin.close();
    if(debug) cout << "Closed input file.\n";
    print_final_computations();
    if(debug || verbose) cout << "\nExecution completed.\n\n";
  }

  void print_types() {
    nTypePairs = nTypes*(nTypes+1)/2;
    if(debug) cout << " Found nTypes =" << nTypes << ", nTypePairs =" << nTypePairs << endl;
    if(debug)
     for(int j=0;j<nTypes;j++)
        cout << "   n. atoms of type " << j << " = " << Nt[j] << " ("<< setw(3) << Nt[j]/(float)N*100 <<"%)" << endl;

    ss.str(std::string());  ss << s_atom_label << tag << ".dat"; fout.open(ss.str(), ios::out);
    for(auto a=0;a<nTypes;a++) {
      if(filetype==FileType::XYZ_CP2K) fout << type_names[a] <<" "<<Nt[a]<<endl;
      else                             fout << a <<" "<<Nt[a]<<endl;
    }
    fout.close();
  }

  void init_computations() {
    init_density();
    if(out_box) init_box();
    if(out_xyz) init_out_xyz();
    if(out_alphanes) init_out_alphanes();
    if(maxshell>0)
    {
      n_b_list = new Neigh_and_Bond_list<ntype,ptype>();
      n_b_list->init(s_rcut, maxshell, p1half, N, nTypes, s_log, tag, debug);
    }
    if(c_coordnum) n_b_list->init_coordnum(s_coordnum, tag, debug);
    if(c_rmin) n_b_list->init_rmin(s_rmin, tag, debug);
    if(c_rmax) n_b_list->init_rmax(s_rmax, tag, debug);
    if(c_bondorient)
    {
      bond_parameters = new Bond_Parameters<ntype,ptype>();
      bond_parameters->init(n_b_list, l, qldot_th, s_bondorient, s_bondcorr, s_nxtal, tag);
    }
    if(c_edq)
    {
      ed_q_calculator = new ED_Bond_Parameter<ntype,ptype>();
      ed_q_calculator->init(n_b_list, s_edq, tag);
    }
    if(c_msd) {
      msd_calculator = new MSD_Calculator<ntype,ptype>();
      msd_calculator->init(dtframe,nframes,period, N, s_msd,s_ngp,tag,debug);
    }
    if(c_rdf) {
      rdf_calculator = new RDF_Calculator<ntype,ptype>();
      rdf_calculator->init(rdf_binw, L, N, V, nTypes, Nt, s_rdf, tag);
    }
    if(c_adf) {
      adf_calculator = new ADF_Calculator<ntype,ptype>();
      adf_calculator->init(adf_binw, s_adf, tag);
    }
    if(c_altbc) {
      altbc_calculator = new ALTBC_Calculator<ntype,ptype>();
      altbc_calculator->init(altbc_rmin, altbc_binw, n_b_list->rcut[0][0], altbc_angle_th, N, V, s_altbc, tag, debug);
    }
    if(c_sq) {
      sq_calculator = new SQ_Calculator<ntype,ptype>();
      sq_calculator->init(qmodmin, qmodmax, qmodstep, L, s_sq, tag);
    }
    if(c_sqt) {
      sqt_calculator = new SQT_Calculator<ntype,ptype>();
      sqt_calculator->init(qmodmin, qmodmax, qmodstep, dtframe,nframes,period, L, s_sqt, tag, debug); // init_rdf();
    }
  }

  void do_computations_and_output(int i) {
    compute_density();
    if(out_box) print_box();
    if(out_xyz) print_out_xyz();
    if(out_alphanes) print_out_alphanes();

    if(maxshell>0) n_b_list->build(timestep, ps, box,boxInv, debug);
    if(c_coordnum) n_b_list->compute_coordnum(timestep, ps, debug);
    if(c_rmin) n_b_list->compute_rmin(timestep, ps, debug);
    if(c_rmax) n_b_list->print_rmax(timestep, debug);
    if(c_bondorient) bond_parameters->compute(timestep, ps, debug);
    if(c_edq) ed_q_calculator->compute(timestep, ps, debug);
    if(c_msd) msd_calculator->compute(i,timestep,ps,box,boxInv,debug);
    if(c_rdf) rdf_calculator->compute(i,nframes,timestep,ps,box,boxInv,debug);
    if(c_adf) adf_calculator->compute(i,nframes,timestep,ps,debug);
    if(c_altbc) altbc_calculator->compute(i,nframes,timestep,ps,debug);
    if(c_sq)
    {
      if(timings) sq_timer.go();
      sq_calculator->compute(i,nframes,timestep,ps,debug);
      if(timings) timing_log( "sq_timing(ms): ", sq_timer.lap() );
    }
    if(c_sqt) sqt_calculator->compute(i,nframes,timestep,ps,debug);
  }

  void print_final_computations() {
    if(c_msd) msd_calculator->print(dtframe,N,debug);
  }

  void timing_log(string comment, float time)
  {
    ss.str(std::string()); ss << s_log << tag; fout.open(ss.str(), ios::app);
    fout << comment << time <<endl;
    fout.close();
  }ntype fcut(ntype x, int pow1, int pow2) {
    ntype x1=1.0, x2=1.0;
    for(auto i=0; i<pow2; i++) {
      x2 *= x;
      if(i==pow1-1) x1=x2;
    }
    return (1.0-x1) / (1.0-x2); // x1 = x^p1, x2 = x^p2
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

//-------------Trajectory output, implemented in io/output.cpp -----------------//
  void init_box();
  void print_box();
  void init_out_xyz();
  void print_out_xyz();
  void init_out_alphanes();
  void print_out_alphanes();

// apply MIC for general periodic boxes
  vec mic(vec& r)
  {
    return r - box * round(boxInv*r);
  }

};

#endif
