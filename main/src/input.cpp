#include "lib/mdtraj.hpp"
using namespace std;
//------------ Reading input files of various formats --------------------//

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
read_frame(fstream &i, bool resetN)
  {
    switch(filetype)
    {
      case FileType::XYZ:
        read_xyz_frame(i, resetN);
        break;
      case FileType::XDATCAR:
        read_xdatcar_frame(i, resetN);
        break;
      case FileType::CONTCAR:
        read_contcar_frame(i, resetN);
        timestep+=dtframe;
        break;
      default:
        cout << "[ERROR: filetype = " << static_cast<int>(filetype) << " not recognized]\n";
        exit(1);
    }
    if(debug) cout << "  Read frame at timestep " << timestep << " DONE.\n";
  }

//-----------------------------------------------------------------------

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
read_xyz_frame(fstream &i, bool resetN)
  {
    string line, a,b,c,d;
    getline(i,line); // first line
    istringstream(line) >> a;
    N = stoi(a);
    if(debug) cout << "\n  Line 1: " << N << " atoms\n";
    getline(i,line); // second line
    istringstream(line) >> b >> c >> d;
    timestep = stoi(d);
    if(debug) cout << "  Line 2: Timestep " << timestep << endl;
    if(resetN) {
      ps.resize(N);
      invN = 1.0/N;
      nframes = nlines / (N+2);
    }
    for(auto &p: ps) p.read_xyz(i); // N particle lines
  }
  
template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
read_contcar_frame(fstream &i, bool resetN)
{
    string line, a,b,c,d;
    ntype s;
    int format;
    
    getline(i,line); // line 1
    if(debug) cout << "\n  Line 1: " << line << endl;
    
    getline(i,line); // line 2
    istringstream(line) >> a;
    s=stof(a);
    if(debug) cout << "  Line 2: scaling factor = " << s << endl;
    
    getline(i,line); // line 3
    istringstream(line) >> a >> b >> c;
    box[0] << s*stof(a), s*stof(b), s*stof(c);
    getline(i,line); // line 4
    istringstream(line) >> a >> b >> c;
    box[1] << s*stof(a), s*stof(b), s*stof(c);
    getline(i,line); // line 5
    istringstream(line) >> a >> b >> c;
    box[2] << s*stof(a), s*stof(b), s*stof(c);
    box = box.T(); // !!!!!!!!!!!
    boxInv = box.inverse();
    if(debug) cout << "  Line 3-5: lattice vectors\n";
    if(debug) box.show();
    if(debug) boxInv.show();
    
    getline(i,line); // line 6
    if(debug) cout << "  Line 6: species = " << line << endl;
    
    getline(i,line); // line 7
    istringstream(line) >> a >> b; // MULTI-SPECIES TO BE IMPLEMENTED
    N = stoi(a);
    if(debug) cout << "  Line 7: number of species 1 = " << N << endl;

    getline(i,line); // line 8
    if(debug) cout << "  Line 8: coordinates format = " << line << endl;
    if(line=="Direct" || line=="direct") format=0; // lattice units
    else if(line=="Cartesian" || line=="cartesian") format=1; // cartesian units
    else { cout << "[ERROR: coordinates format not recognized: " << line << "]\n"; exit(1); }
    
    if(resetN) {
      ps.resize(N);
      invN = 1.0/N;
      nframes = nlines / (8+N+ 9+N+ 4+3*N); // N positions, N velocities, 3*N predictor-corrector
    }
    for(auto &p: ps) {
      p.read_3cols(i); // N particle lines
      if(format==0) p.r = box*p.r; // x' = (x*ax + y*bx + z*cx) --> r' = r*Box
      else if(format==1) p.r*=s;
    }

    for(auto j=0;j<9+N;j++) getline(i,line); // velocities
    for(auto j=0;j<4+3*N;j++) getline(i,line); // predictor-corrector
}

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
read_xdatcar_frame(fstream &i, bool resetN)
{
    string line, a,b,c,d;
    ntype s;
    int format;
    
    getline(i,line); // line 1
    if(debug) cout << "\n  Line 1: " << line << endl;
    
    getline(i,line); // line 2
    istringstream(line) >> a;
    s=stof(a);
    if(debug) cout << "  Line 2: scaling factor = " << s << endl;
    
    getline(i,line); // line 3
    istringstream(line) >> a >> b >> c;
    box[0] << s*stof(a), s*stof(b), s*stof(c);
    getline(i,line); // line 4
    istringstream(line) >> a >> b >> c;
    box[1] << s*stof(a), s*stof(b), s*stof(c);
    getline(i,line); // line 5
    istringstream(line) >> a >> b >> c;
    box[2] << s*stof(a), s*stof(b), s*stof(c);
    box = box.T(); // box = (a|b|c) where a,b,c = lattice vectors as columns !!!!!!!!!!!
    boxInv = box.inverse();
    if(debug) cout << "  Line 3-5: lattice vectors\n";
    
    getline(i,line); // line 6
    if(debug) cout << "  Line 6: species = " << line << endl;
    
    getline(i,line); // line 7
    istringstream(line) >> a >> b; // MULTI-SPECIES TO BE IMPLEMENTED
    N = stoi(a);
    if(debug) cout << "  Line 7: number of species 1 = " << N << endl;

    getline(i,line); // line 8
    istringstream(line) >> a >> b >> c;
    timestep=stoi(c);
    if(debug) cout << "  Line 8: coordinates format = " << a << ", timestep = " << timestep << endl;
    if(a=="Direct" || a=="direct") format=0; // lattice units
    else if(a=="Cartesian" || a=="cartesian") format=1; // cartesian units
    else { cout << "[ERROR: coordinates format not recognized: " << a << "]\n"; exit(1); }
    
    if(resetN) {
      ps.resize(N);
      invN = 1.0/N;
      nframes = nlines / (8+N); // 8 system info + N position lines
    }
    for(auto &p: ps) {
      p.read_3cols(i); // N particle lines
      // from direct to cartesian
      if(format==0) p.r = box*p.r; // x' = (x*ax + y*bx + z*cx) etc.
      // already cartesian: must only be scaled by s
      else if(format==1) p.r*=s;
    }
}
