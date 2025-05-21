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
#include "lib/Ycomplex_builtin.hpp"
#include "lib/pbc.hpp"
#include "lib/logtimesteps.hpp"
using namespace std;

const string root_path="/home/flavio/programmi/mdtraj/mixture";
#define MAX_N_TYPES 10
#define MAX_N_ANGMOM 5
enum class FileType {
  NONE, XYZ, XYZ_CP2K, CONTCAR, POSCAR, XDATCAR, XDATCARV, ALPHANES, ALPHANES9, JMD, LAMMPSDATA, LAMMPSTRJ, YUHAN, RUNNER
};

template<class ntype, class ptype>
class Trajectory {
public:
  using vec=myvec<ntype,3>;
  using mat=mymatrix<ntype,3,3>;
  int N, nTypes, nTypePairs, Nt[MAX_N_TYPES]; // number of particles, number of types, num of particles for each type
  int max_nTypes; // in case of dynamic types (e.g. clusters)
  string type_names[MAX_N_TYPES]; // name for each type
  vec L; // length of box vectors ( ??? centered in 0: -Lx/2 < x < Lx/2 )
  mat box, boxInv; // most general simulation box
  int image_convention; // 1 (Minimum Image), 0 (Cluster), -1 (All within the cutoff radius)
  bool remove_rot_dof; // remove the 3 rotational degrees of freedom from the box?
  bool pbc_out; // print output with PBC?
  ntype V, mdens, ndens; // volume, mass density, nuerical density
  vector<ptype> ps, ps_new; // vector of particles
  int nframes, timestep;
  bool c_coordnum, c_nna,c_nnd, c_bondorient, c_msd, c_rdf, c_adf, c_rmin, c_rmax;
  bool c_altbc, c_sq, c_sqt, c_edq, c_clusters, c_pmp, c_oct; // compute or not
  bool dynamic_types;
  string s_in, s_out, s_rcut, s_rcut_clusters, tag, s_logtime, s_atom_label, s_box, s_ndens, s_coordnum, s_clusters;
  string s_nna,s_nnd, s_bondorient, s_bondcorr, s_nxtal, s_msd, s_ngp, s_overlap, s_rdf, s_adf;
  string s_rmin, s_tbc, s_altbc, s_sq, s_sqt, s_log, s_rmax, s_edq, s_pmp, s_oct; // for file naming
  bool ignore_double_frames, logtime, nodynamics, out_box, out_xyz, out_alphanes, out_lammpsdump;
  bool debug, verbose;
  //
  LogTimesteps logt;
  //
  int maxsphere,max_num_nna, max_num_nnd; // <= MAX_NSPHERE
  vecflex<ntype> defaultCutoff[MAX_NSPHERE];
  Neigh_and_Bond_list<ntype,ptype> *n_b_list;
  //
  int l, num_l, l_list[MAX_N_ANGMOM];
  string s_l_list[MAX_N_ANGMOM];
  ntype qldot_th;
  Bond_Parameters<ntype,ptype>* bond_parameters[MAX_N_ANGMOM];
  ED_Bond_Parameter<ntype,ptype>* ed_q_calculator;
  PatternMatchingParameters<ntype,ptype>* pmp_calculator;
  OctahedralParameter<ntype,ptype>* oct_calculator;
  AtomicTimeCorrelator<ntype> *oct_correlator;
  //
  ntype rdf_binw, rdf_rmax;
  PBC<ntype> *pbc;
  RDF_Calculator<ntype,ptype> *rdf_calculator;
  //
  int qmodmin,qmodmax,qmodstep;
  SQ_Calculator<ntype,ptype> *sq_calculator;
  //
  ntype adf_binw;
  ADF_Calculator<ntype,ptype> *adf_calculator;
  //
  int period; // this is in timestep units, not in number of frames
  ntype Qoverlap_cutoff;
  bool msdAverageOverTime0;
  MSDU_Calculator<ntype,ptype> *msd_calculator;
  //
  ntype altbc_rmin, altbc_binw, altbc_angle_th;
  ALTBC_Calculator<ntype,ptype> *altbc_calculator;
  //
  SQT_Calculator<ntype,ptype> *sqt_calculator;

private:
  string myName="mdtraj";
  bool timings;
  int nlines, t0frame, dtframe,last_dtframe,last_timestep, nskip0, nskip1, nframes_original;
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
  Trajectory(){}
  ~Trajectory(){}

//------------ Input reading and interaction --------------------//

  void print_usage(char argv0[]);
  void print_summary();
  void args(int argc, char** argv);

  void print_state(){
    cerr << "Summary of parameters:\n";
    cerr << " debug = \t " << debug << endl;
    cerr << " verbose = \t " << debug << endl;
    cerr << " c_coordnum = \t " << c_coordnum << endl;
    cerr << " c_bondorient = \t " << c_bondorient << endl;
    cerr << " c_msd = \t " << c_msd << endl;
    cerr << " c_rdf = \t " << c_rdf << endl;
    cerr << " c_adf = \t " << c_adf << endl;
    cerr << " c_rmin = \t " << c_rmin << endl;
    cerr << " c_rmax = \t " << c_rmax << endl;
    cerr << " c_altbc = \t " << c_altbc << endl;
    cerr << " c_sq = \t " << c_sq << endl;
    cerr << " c_sqt = \t " << c_sqt << endl;
    cerr << " c_edq = \t " << c_edq << endl;
    cerr << " c_nna = \t " << c_nna << endl;
    cerr << " c_nnd = \t " << c_nnd << endl;
    cerr << " c_pmp = \t " << c_pmp << endl;
    cerr << " c_oct = \t " << c_oct << endl;
    cerr << " c_clusters = \t " << c_clusters << endl;
    cerr << " dynamic_types = \t " << dynamic_types << endl;
    cerr << " angular momentum for ql: l = \t " << l << endl;
    cerr << " qldot threshold = \t " << qldot_th << endl;
    cerr << " box (a|b|c) = \t "; box.show();
    cerr << " total volume V = \t " << V << endl;
    cerr << " box inverse = \t "; boxInv.show();
    cerr << " L (length of each box vector) = \t "; L.show();
    cerr << " p1 = \t " << p1 << endl;
    cerr << " p2 = \t " << p2 << endl;
    cerr << " period for MSD & NGP & S(q,t) = \t " << period << endl;
    cerr << " ignore_double_frames = \t " << ignore_double_frames << endl;
    cerr << " logtime = \t " << logtime << endl;
    cerr << " s_logtime = \t " << s_logtime << endl;
    cerr << " rdf_binw = \t " << rdf_binw << endl;
    cerr << " rdf_rmax = \t " << rdf_rmax << endl;
    cerr << " adf_binw = \t " << adf_binw << endl;
    cerr << " q_mod_min,q_modmax,q_mod_step = \t " << qmodmin << ", " << qmodmax << ", "<< qmodstep << endl;
    cerr << " remove rotational degrees of freedom = \t " << remove_rot_dof << endl;
    cerr << " nodynamics = \t " << nodynamics << endl;
    cerr << " out_box = \t " << out_box << endl;
    cerr << " out_xyz = \t " << out_xyz << endl;
    cerr << " out_lammpsdump = \t " << out_lammpsdump << endl;
    cerr << " out_alphanes = \t " << out_alphanes << endl;
    cerr << " pbc_out = \t " << pbc_out << endl;
    cerr << " fskip_from_beginning = \t " << fskip0 << endl;
    cerr << " fskip_from_end = \t " << fskip1 << endl;
    cerr << " timings = \t " << timings << endl;
    cerr << " tag = \t " << tag << endl;
    cerr << " s_in = \t " << s_in << endl;
    cerr << endl;
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
    c_nna = false;
    c_nnd = false;
    c_pmp = false;
    c_oct = false;
    c_clusters = false;
    dynamic_types = false;

    filetype=FileType::NONE;
    s_in="__NOT_DEFINED__";

    s_ndens="ndens";
    s_box="box";
    s_coordnum="coordnum";
    s_bondorient="boo";
    s_bondcorr="boc";
    s_nxtal="nc";
    s_msd="msd";
    s_ngp="ngp";
    s_overlap="Qoverlap";
    s_nna="nna";
    s_nnd="nnd";
    s_rdf="rdf";
    s_adf="adf";
    s_rmin="rmin";
    s_rmax="rmax";
    s_altbc="altbc";
    s_sq="sq";
    s_sqt="sqt";
    s_edq="ed_q";
    s_pmp="q_pmp";
    s_oct="q_oct";
    s_clusters="clusters";
    tag="";
    s_out="traj";
    s_log="log";
    s_atom_label="labels";
    s_rcut="__NOT_DEFINED__";
    s_rcut_clusters="__NOT_DEFINED__";
    ss.str(std::string()); ss << s_log << tag; fout.open(ss.str(), ios::out);
    fout << "LOG SUMMARY"<<endl;
    fout.close();

    period = -1; // default: don't average over t0 for MSD
    image_convention = 1; // Minimum Image Convention
    remove_rot_dof = false;
    nodynamics = false; // ignore timestep values?
    out_box = false;
    out_xyz = false;
    out_lammpsdump = false;
    ignore_double_frames = false;
    logtime = false;
    s_logtime="";
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
    defaultCutoff[0].resize(1); // 1st sphere
    defaultCutoff[0][0] = 3.75; // Antimony: 3.6 in glass, 3.75-3.89 in xtal
    defaultCutoff[1].resize(1); // 2nd sphere
    defaultCutoff[1][0] = 5.15;
    defaultCutoff[2].resize(1); // 3rd sphere
    defaultCutoff[2][0] = 8.8;
    //l=4;
    num_l=1;
    l_list[0]=4;
    s_l_list[0]="l.4";
    p1half=6;
    qldot_th = 0.65;
    rdf_binw=0.0;
    rdf_rmax=0.0;
    qmodmin=2;
    qmodmax=100;
    qmodstep=1;
    adf_binw=0.0;
    altbc_binw=0.0;
    altbc_rmin=0.0;
    altbc_angle_th=-1.0;
    fskip0=fskip1=0.0;
    Qoverlap_cutoff = 2.0; // default 2 angstrom is ok for antimony

    //-------- Update parameters with input arguments: -----------//
    args(argc, argv);

    // Compute non-primitive parameters:
    get_angular_momentum_list();
    p2half = 2*p1half;
    p1 = 2*p1half;
    p2 = 2*p2half;
    if(c_bondorient) maxsphere=MAX_NSPHERE; // init all neigh spheres
    else if(c_coordnum || c_nna || c_nnd || c_adf || c_rmin || c_rmax || c_altbc || c_edq || c_pmp || c_oct) maxsphere=1;       // init only first neigh sphere
    else                         maxsphere=0;       // do not init any
    if(c_clusters && !c_bondorient) {
      cerr << "ERROR: cannot compute clusters without computing BOC parameters!\n";
      exit(1);
    }

    // Print a recap:
    if(debug) { cerr << "State after reading args():\n"; print_state(); }

    if(filetype==FileType::NONE) {
      cerr << "ERROR: input file was not defined. Use -h for help\n";
      exit(1);
    }

    // check conflicts
    if(nodynamics){
      cerr<<"WARNING: '-nodynamics' may produce junk ensemble averages.\n";
      if(logtime && nodynamics){
        cerr<<"WARNING: combination of '-nodynamics' and '-logtime' will LIKELY produce junk.\n";
      }
    }
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
      cerr << "[ "<<myName<<" Warning ] det(box)="<<V<<" follows left-hand rule. ]\n";
      V=-V;
    }
    else if (V==0.) {
      box.show();
      cerr << "[ "<<myName<<" Error ] det(box)=0.0 not supported. ]\n";
      exit(1);
    }
  }

  void read_frame(fstream &i, bool resetN, bool reset_nTypes, int frameIdx);
  void removeRotDof();
  void read_contcar_frame(fstream &i, bool resetN);
  void read_poscar_frame(fstream &i, bool resetN);
  void read_xdatcar_frame(fstream &i, bool resetN, bool constantBox);
  void read_xyz_frame(fstream &i, bool resetN);
  void read_xyz_cp2k_frame(fstream &i, bool resetN);
  void read_alphanes_frame(fstream &i, bool resetN);
  void read_alphanes9_frame(fstream &i, bool resetN);
  void read_jmd_frame(fstream &i, bool resetN);
  void read_lammpsdata_frame(fstream &i, bool resetN, bool reset_nTypes);
  void read_lammpstrj_frame(fstream &i, bool resetN, bool reset_nTypes);
  void read_yuhan_frame(fstream &i, bool resetN, bool isFirstFrame);
  void read_runner_frame(fstream &i, bool resetN);

  //------- COMPUTE things ---------------//

  void run()
  {
    PrintProgress printProgress;
    if(verbose || debug) {cerr << "#Begin of run():\n"; print_state();}
    nlines = getLineCount(s_in);
    if(debug) cerr << "#Read " << nlines << " lines in file " << s_in << ". Opening again for reading trajectory.\n";
    fin.open(s_in, ios::in);

    // set manual time for the following file formats:
    bool set_manual_time=(filetype==FileType::CONTCAR ||
                          filetype==FileType::POSCAR ||
                          filetype==FileType::ALPHANES ||
                          filetype==FileType::ALPHANES9);
    if(set_manual_time) { timestep=-1; dtframe=1; }

    //---------- Read 1st frame -------------//
    read_frame(fin, true, true, 0);
    t0frame = timestep;
    nframes_original = nframes;
    if(debug) cerr << "# I read first frame: set N = " << N << " (I assume it's constant) and t0frame = " << t0frame << "\n";
    if(debug) cerr << "#  From N, I deduce nframes = " << nframes << "\n";
    print_types();

    // this is valid for trajectories with linear timesteps
    nskip0=int(fskip0*nframes_original);
    nskip1=int(fskip1*nframes_original);
    nframes = nframes_original - nskip0 - nskip1;

    // read log timesteps full schedule and deduce nskip's
    if(logtime) {
      logt.deduce_fromfile(s_logtime, debug);
      nframes = logt.apply_fskip(fskip0,fskip1);
      nskip0 = logt.ncyc_skip0*logt.npc;
      if(debug) logt.print_summary();
      if(t0frame!=logt.get_first_timestep_noskip()){
        cerr<<"ERROR: first timestep="<<t0frame<<" =/= first scheduled timestep="
            <<logt.get_first_timestep_noskip()<<endl;
            exit(1);
      }
      // nskip0 is correct, but nskip1 may not be if trajectory was early-stopped,
      // so compute it by subtracting nframes and nskip0 to the actual (original) n. of frames:
      nskip1 = nframes_original - nskip0 - nframes;
    }
    if(debug||verbose){
      if(logtime) cerr << "# Log-skipping applied:    nskip0="<<nskip0<<" nframes="<<nframes<<" => nskip1="<<nskip1<<endl;
      else        cerr << "# Linear-skipping applied: nskip0="<<nskip0<<" nskip1="<<nskip1<<" => nframes="<<nframes<<endl;
    }

    if(nframes<1) { cerr << "[ "<<myName<<" Error ] skipped too many frames.\n  Total: "<<nframes_original<<"; Skipped: "<<nskip0<<"+"<<nskip1<<"; Remaining: "<<nframes<<" ]\n\n"; exit(1);}

    //---------- Read 2nd frame (if it exists) -------------//
    try
    {
      read_frame(fin, false, dynamic_types, 1);
      if(!set_manual_time){ dtframe=timestep-t0frame; }
      if(debug) cerr << "# I read the second frame: dtframe = "<<dtframe<<endl;
    } catch (...) {
      cerr << "WARNING: only 1 frame in trajectory.\n";
      dtframe=last_dtframe=1; // this is meaningless, but it avoids nonsense later
    }
    fin.close();

    // Restart reading!!

    // check max nTypes
    max_nTypes=nTypes;
    if(dynamic_types) {
      if(debug||verbose) { cerr << "# dynamic_types: checking for max number of types in trajectory\n"; }
      if(debug) { cerr << "# dynamic_types: # frame nTypes\n"; }
      fin.open(s_in, ios::in);
      for(int i=0; i<nframes_original; i++)
      {
        if(i+1>nskip0 && i<nframes_original-nskip1){
          read_frame(fin, false, true, i); // reset ntypes but don't reset N and nframes
          max_nTypes = max(max_nTypes, nTypes);
          if(debug) { cerr << "# dynamic_types: "<< i << " " << nTypes << endl; }
        } else{
          read_frame(fin, false, false, i);
        }
        if(N != ps.size()) { cerr << "[Error: N has changed]\n"; exit(1);}
      }
      fin.close();
      if(debug||verbose) { cerr << "# dynamic_types: found max "<< max_nTypes << " types\n"; }
    }

    if(debug) cerr << "Initialization of Computations STARTED\n";
    init_computations();
    if(debug) cerr << "Initialization of Computations COMPLETED\n";

    if(debug || verbose) cerr << "#------- MAIN LOOP ------#\n";
    fin.open(s_in, ios::in);
    printProgress.init( nframes, 2000 ); // update % every 2000 ms
    if(filetype==FileType::CONTCAR || filetype==FileType::CONTCAR ||
      filetype==FileType::ALPHANES || filetype==FileType::ALPHANES9) {timestep=-1; } // set manual time
    string junk_line;
    if(filetype==FileType::XDATCARV) { for(int i=0;i<7;i++) getline(fin, junk_line); } // skip first 7 lines (so that you can use resetN=false)
    timer.go();
    for(int i=0; i<nframes_original; i++)
    {
      if(i+1>nskip0 && i<nframes_original-nskip1)
        read_frame(fin, false, dynamic_types, i);
      else
        read_frame(fin, false, false, i);

      if(N != ps.size()) { cerr << "[Error: N has changed]\n"; exit(1);}

      if(!set_manual_time){
        if(i==0) dtframe=0;
        else     dtframe=timestep-last_timestep;
      }

      if(i>1 && dtframe==0) {
        if(ignore_double_frames){
          cerr<<"["<<myName<<" Warning] Skipping timestep= "<<timestep<<" because dt==0 and ignore_double_frames=="<<ignore_double_frames<<"\n";
          continue; //!!!
        }else{
          cerr<<"["<<myName<<" Error] dt==0 at timestep= "<<timestep<<" and ignore_double_frames=="<<ignore_double_frames<<"\n";
          exit(1);
        }
      }

      if(i>1 && !nodynamics && !logtime && dtframe!=last_dtframe) {
        cerr<<"["<<myName<<" Error] nodynamics==false and logtime==false, but timestep interval is NOT linearly spaced:\n";
        cerr<<"               last_dt = "<<last_dtframe<<", dt = "<<dtframe<<", at timestep = "<<timestep<<"\n";
        exit(1);
      }

      if(i+1>nskip0 && i<nframes_original-nskip1)
        do_computations_and_output(i-nskip0); // use index 0,1,...,nframes-1

      printProgress.update( i+1-nskip0, timer.lap() );
      last_dtframe=dtframe;
      last_timestep=timestep;
    }
    printProgress.end();
    fin.close();
    if(debug) cerr << "Closed input file.\n";
    print_final_computations();
    if(debug || verbose) cerr << "\nExecution completed.\n\n";
  }

  void print_types() {
    if(nTypes<=0) {
      cerr << "ERROR: nTypes = "<<nTypes<<" is not admitted.\n";
      exit(1);
    }
    nTypePairs = nTypes*(nTypes+1)/2;
    if(debug) cerr << " Found nTypes =" << nTypes << ", nTypePairs =" << nTypePairs << endl;
    if(debug)
     for(int j=0;j<nTypes;j++)
        cerr << "   n. atoms of type " << j << " = " << Nt[j] << " ("<< setw(3) << Nt[j]/(float)N*100 <<"%)" << endl;

    ss.str(std::string());  ss << s_atom_label << tag << ".dat"; fout.open(ss.str(), ios::out);
    for(auto a=0;a<nTypes;a++) {
      if(filetype==FileType::XYZ_CP2K || filetype==FileType::XYZ ||
        filetype==FileType::CONTCAR || filetype==FileType::POSCAR ||
        filetype==FileType::XDATCAR || filetype==FileType::XDATCARV ||
        filetype==FileType::YUHAN) fout << type_names[a] <<" "<<Nt[a]<<endl;
      else                             fout << a <<" "<<Nt[a]<<endl;
    }
    fout.close();
  }

  void init_computations() {
    init_density();
    if(out_box) init_box();
    if(out_xyz) init_out_xyz();
    if(out_lammpsdump&&!c_clusters) init_out_lammpsdump(s_out);
    if(out_alphanes) init_out_alphanes();

    if(maxsphere>0 || c_rdf || c_msd)
    {
      pbc = new PBC<ntype>();
      // rmax>0 is currently implemented only for g(r)
      if(c_rdf) pbc->init(image_convention,rdf_rmax,box,debug);
      else      pbc->init(image_convention,      -1,box,debug);
    }
    if(maxsphere>0)
    {
      n_b_list = new Neigh_and_Bond_list<ntype,ptype>();
      n_b_list->init(s_rcut, defaultCutoff, maxsphere, p1half, N, max_nTypes, s_log,
        tag, debug, verbose);
    }
    if(c_coordnum) n_b_list->init_coordnum(s_coordnum);
    if(c_nna) n_b_list->init_nearest_neigh_angles(s_nna,max_num_nna);
    if(c_nnd) n_b_list->init_nearest_neigh_dists(s_nnd,max_num_nnd);
    if(c_rmin) n_b_list->init_rmin(s_rmin);
    if(c_rmax) n_b_list->init_rmax(s_rmax);
    if(c_bondorient)
    {
      for(auto l_=0;l_<num_l;l_++){
        bond_parameters[l_] = new Bond_Parameters<ntype,ptype>();
        bond_parameters[l_]->init(n_b_list, l_list[l_], qldot_th, s_bondorient, s_bondcorr,
          s_nxtal, tag, debug, verbose);
        if(c_clusters) {
          ss.str(std::string()); ss<<s_clusters<<s_l_list[l_];
          n_b_list->init_clusterize( ss.str(),
                                    s_rcut_clusters, defaultCutoff[0]);
          if(out_lammpsdump) {
            ss.str(std::string()); ss<<s_clusters<<s_l_list[l_];
            init_out_lammpsdump(ss.str());
          }
        }
      }
    }
    if(c_edq)
    {
      ed_q_calculator = new ED_Bond_Parameter<ntype,ptype>();
      ed_q_calculator->init(n_b_list, s_edq, tag, debug, verbose);
    }
    if(c_pmp) {
      pmp_calculator = new PatternMatchingParameters<ntype,ptype>();
      pmp_calculator->init(n_b_list, s_pmp, tag, debug, verbose);
    }
    if(c_oct) {
      oct_calculator = new OctahedralParameter<ntype,ptype>();
      oct_calculator->init(n_b_list, s_oct, tag, debug, verbose);
      oct_correlator = new AtomicTimeCorrelator<ntype>();
      oct_correlator->init(dtframe, nframes, period, N, logtime, &logt,
        s_oct+"_time", tag,debug,verbose);
    }
    if(c_msd) {
      msd_calculator = new MSDU_Calculator<ntype,ptype>();
      msd_calculator->init(dtframe,nframes,period, Qoverlap_cutoff, N,nTypes,Nt,
        logtime,logt, s_msd,s_ngp,s_overlap,tag, debug,verbose);
    }
    if(c_rdf) {
      rdf_calculator = new RDF_Calculator<ntype,ptype>();
      rdf_calculator->init(rdf_binw, rdf_rmax, box, N, V, nTypes, Nt, logtime,logt,
        s_rdf, tag, debug, verbose);
    }
    if(c_adf) {
      adf_calculator = new ADF_Calculator<ntype,ptype>();
      adf_calculator->init(adf_binw, nTypes, s_adf, tag, debug, verbose);
    }
    if(c_altbc) {
      altbc_calculator = new ALTBC_Calculator<ntype,ptype>();
      altbc_calculator->init(altbc_rmin, altbc_binw, n_b_list->rcut[0][0],
        altbc_angle_th, max_nTypes, N, V, s_altbc, tag, debug, verbose);
    }
    if(c_sq) {
      sq_calculator = new SQ_Calculator<ntype,ptype>();
      sq_calculator->init(qmodmin, qmodmax, qmodstep, L, logtime,logt,
        s_sq, tag, debug, verbose);
    }
    if(c_sqt) {
      sqt_calculator = new SQT_Calculator<ntype,ptype>();
      sqt_calculator->init(qmodmin, qmodmax, qmodstep, dtframe, nframes,
        period, L, nTypes, Nt, logtime,logt, s_sqt, tag, debug, verbose);
    }
  }

  void do_computations_and_output(int i) {
    compute_density();
    if(maxsphere>0) n_b_list->build(timestep, ps, box,boxInv);
    if(c_coordnum) n_b_list->compute_coordnum(timestep, ps);
    if(c_nna) n_b_list->compute_nearest_neigh_angles(timestep, ps);
    if(c_nnd) n_b_list->compute_nearest_neigh_dists(timestep, ps);
    if(c_rmin) n_b_list->compute_rmin(timestep, ps);
    if(c_rmax) n_b_list->print_rmax(timestep);
    if(c_bondorient) {
      for(int l_=0;l_<num_l;l_++){
        bond_parameters[l_]->compute(timestep, ps);
        if(c_clusters) {
          ss.str(std::string()); ss<<s_clusters<<s_l_list[l_];
          n_b_list->clusterize(timestep, ps, ss.str(),
                              bond_parameters[l_]->ql_dot, bond_parameters[l_]->qldot_th);
          if(out_lammpsdump){ // visualization of clusters
            int c;
            vector<int> original_labels;
            original_labels.resize(ps.size());
            for(auto i=0;i<ps.size();i++){
              c=n_b_list->cluster_of_particle[i];
              original_labels[i]=ps[i].label;
              // label 1 for particles belonging to no cluster, label 2,3,... for clusters of decreasing size
              if(c>=0) ps[i].label = 2+(int)(n_b_list->cluster_permutation_by_size[c]);
              else     ps[i].label = 1;
            }
            ss.str(std::string()); ss<<s_clusters<<s_l_list[l_];
            print_out_lammpsdump(ss.str());
            for(auto i=0;i<ps.size();i++) ps[i].label=original_labels[i]; //!!!
          }
        }
      }
    }
    if(c_edq) ed_q_calculator->compute(timestep, ps);
    if(c_pmp) pmp_calculator->compute(timestep, ps);
    if(c_oct) {
      oct_calculator->compute(timestep, ps);
      oct_correlator->correlate(i,timestep,oct_calculator->my_q_oct);
    }
    if(c_msd) msd_calculator->compute(i,timestep,ps,pbc);
    if(c_rdf) rdf_calculator->compute(i,nframes,timestep,ps,pbc);
    if(c_adf) adf_calculator->compute(i,nframes,timestep,ps);
    if(c_altbc) altbc_calculator->compute(i,nframes,timestep,ps);
    if(c_sq)
    {
      if(timings) sq_timer.go();
      sq_calculator->compute(i,nframes,timestep,ps);
      if(timings) timing_log( "sq_timing(ms): ", sq_timer.lap() );
    }
    if(c_sqt) sqt_calculator->compute(i,nframes,timestep,ps);

    if(out_box) print_box();
    if(out_xyz && !c_clusters) print_out_xyz();
    if(out_lammpsdump && !c_clusters) print_out_lammpsdump(s_out);
    if(out_alphanes) print_out_alphanes();
  }

  void print_final_computations() {
    if(c_oct) {
      oct_correlator->print(dtframe);
    }
    if(c_msd) msd_calculator->print(dtframe);
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

  void get_angular_momentum_list() {
    // get tag name from l_list for each l_
    for(int i=0;i<num_l;i++){
      ss.str(std::string());
      ss<<".l"<<std::to_string(l_list[i]);
      s_l_list[i]=ss.str();
      ss.clear();
    }
    /*
    // split the input INTEGER into digits, e.g. l=46 into a vector containing 4,6
    string segment;
    stringstream sss;
    ss.str(std::string()); ss<<std::to_string(l);
    int i=0, ic;
    char c;
    while (ss.get(c) && i<MAX_N_ANGMOM){
      // vector of integers
      ic = c - '0';
      l_list[i] = ic;//stoi(segment);
      // vector of file tags
      sss.str(std::string()); sss<<".l"<<ic;//stoi(segment);
      s_l_list[i] = sss.str();
      i++;
    }
    ss.clear(); //!!!
    num_l=i;
    */

    /* // old method
    if(l>0&&l<=6){
      num_l=1;
      l_list[0]=l;
    } else if(l==46 || l==64) {
      num_l=2;
      l_list[0]=6;
      l_list[1]=4;
    } else {
      cerr<<"[ "<<myName<<" Error ] l="<<l<<" not supported. Multi-angular-momentum is still experimental.\n";
      exit(1);
    }
    */

    if(debug) {
      cerr<<"[ "<<myName<<" ] Set angular momentum list:";
      for(int l_=0;l_<num_l;l_++){
        cerr<<" "<<l_list[l_];
      }
      cerr<<endl;
      cerr<<"[ "<<myName<<" ] with the following file tags:";
      for(int l_=0;l_<num_l;l_++){
        cerr<<" "<<s_l_list[l_];
      }
      cerr<<endl;
    }
    return;
  }

//-------------Trajectory output, implemented in io/output.cpp -----------------//
  void init_box();
  void print_box();

  void init_out_xyz();
  void print_out_xyz();

  void init_out_alphanes();
  void print_out_alphanes();

  void init_out_lammpsdump(string);
  void print_out_lammpsdump(string);

};

#endif
