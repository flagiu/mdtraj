using namespace std;
//------------ Reading input files of various formats --------------------//

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
read_frame(fstream &i, bool resetN, int frameIdx)
  {
    switch(filetype)
    {
      case FileType::XYZ:
        read_xyz_frame(i, resetN);
        break;
      case FileType::XYZ_CP2K:
        read_xyz_cp2k_frame(i, resetN);
        break;
      case FileType::XDATCAR:
        read_xdatcar_frame(i, resetN, false);
        break;
      case FileType::XDATCARV:
        read_xdatcar_frame(i, resetN, true);
        break;
      case FileType::CONTCAR:
        read_contcar_frame(i, resetN);
        timestep+=dtframe; // artificial time
        break;
      case FileType::ALPHANES:
        read_alphanes_frame(i, resetN);
        timestep+=dtframe; // artificial time
        break;
      case FileType::ALPHANES9:
        read_alphanes9_frame(i, resetN);
        timestep+=dtframe; // artificial time
        break;
      case FileType::JMD:
        read_jmd_frame(i, resetN);
        break;
      case FileType::LAMMPSTRJ:
        read_lammpstrj_frame(i, resetN);
        break;
      case FileType::YUHAN:
        read_yuhan_frame(i, resetN, frameIdx==0);
        break;
      default:
        cout << "[ERROR: filetype = " << static_cast<int>(filetype) << " not recognized]\n";
        exit(1);
    }
    if(remove_rot_dof) removeRotDof();
    boxInv = box.inverse();
    set_L_from_box();
    if(c_sq || c_sqt)
    {
      for(auto &p: ps) p.s = boxInv * p.r; // pre-compute fractional coordinates
    }
    if(debug) cout << "  Read frame at timestep " << timestep << " DONE.\n";
  }

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
removeRotDof()
{
  mat R1,R2,Rtot;
  vec u, xdir;
  ntype a_xy_radius, c, s, t, angle;
  if(debug) cout <<"\n*------ Removing the 3 rotational degrees of freedom ------*\n";
  if(debug) cout <<"\n[ Remember that box=(a|b|c) ]\n";
  // let a, b, c be the columns of the box matrix.
  // 1) Align a to x-axis
  // u=(0,az,-ay) or (0,-az,ay) is the axis of rotation (to be normalized)
  if(debug) { cout << "Raw box = "; box.show(); }
  if(debug) { cout << "Raw volume (with sign)= "<< box.det()<<endl; }
  u << 0.0, box[2][0], -1*box[1][0];
  if(debug) { cout << "  axis1 = "; u.show(); }
  // c = cos(angle) = ax/|a| ; s = sin(angle)
  c = box[0][0] / sqrt( box[0][0]*box[0][0] + box[1][0]*box[1][0] + box[2][0]*box[2][0] );
  s = sqrt( 1.0 - c*c ); // + or - ? with the choice of u=(0,az,-ay) it should be +
  if(debug) { cout << "  cos1 = "<< c << endl; }
  if(debug) { cout << "  sin1 = "<< s << endl; }
  // build the rotation matrix
  R1 = rotation_matrix_axis_cossin( u, c, s ); // declared in lib/matrix.hpp
  if(box.det()<0) R1 *= -1.; // invert to go back to right-hand rule
  if(debug) { cout << "Rotation matrix R1 = "; R1.show(); }
  box = R1*box;
  if(debug) { cout << "box (after R1) = "; box.show(); }
  // 2) Rotate around the x-axis to bring the new b within the x-y plane
  xdir << 1.0, 0.0, 0.0;
  if(debug) { cout << "  axis2 (x-axis) = "; xdir.show(); }
  // tan(angle) = -bz / by
  angle = -atan2( box[2][1], box[1][1] );
  if(debug) { cout << "  angle2 (deg) = "<<angle*180/M_PI<<endl; }
  R2 = rotation_matrix_axis_cossin( xdir, cos(angle), sin(angle) );
  if(debug) { cout << "Rotation matrix R2 = "; R2.show(); }
  box = R2*box;
  if( abs(box[1][1])<1e-15 ) box[1][1]=0.0; // remove rounding errors
  if( abs(box[1][1])<1e-15 ) box[2][1]=0.0;
  if( abs(box[2][2])<1e-15 ) box[2][2]=0.0;
  if(debug) { cout << "box (after R2) = "; box.show(); }
  // 3) Apply both rotations to each particle
  Rtot=R2*R1;
  for(auto &p: ps) p.r = Rtot*p.r;
  // Final result: A=(Ax,0,0)^T, B=(Bx,By,0)^T, C generic ==> box is upper diagonal
  if(debug) { cout <<"\n*------------------------------------------------------------*\n"; }
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
  read_xyz_cp2k_frame(fstream &i, bool resetN)
    {
      string line, a,b,c,d,e;
      getline(i,line); // first line
      istringstream(line) >> a;
      N = stoi(a);
      if(debug) cout << "\n  Line 1: " << N << " atoms\n";
      getline(i,line); // second line
      istringstream(line) >> b >> c >> d >> e;
      if(debug) cout << "d="<<d<<endl;
      timestep = stoi(d);
      if(debug) cout << "  Line 2: Timestep " << timestep << endl;
      if(resetN) {
        ps.resize(N);
        invN = 1.0/N;
        nframes = nlines / (N+2);
      }
      for(auto &p: ps)
      {
        getline(i, line);
        istringstream(line) >> a >> b >> c >> d;
        p.r[0] = stof(b);
        p.r[1] = stof(c);
        p.r[2] = stof(d);
      }
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

    for(auto j=0;j<9+N;j++) getline(i,line); // skip velocities
    for(auto j=0;j<4+3*N;j++) getline(i,line); // skip predictor-corrector
}

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
read_xdatcar_frame(fstream &i, bool resetN, bool constantBox)
{
    string line, a,b,c,d;
    ntype s;
    int format;
    if(resetN || !constantBox) {
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
	    if(debug) cout << "  Line 3-5: lattice vectors\n";

	    getline(i,line); // line 6
	    if(debug) cout << "  Line 6: species = " << line << endl;

	    getline(i,line); // line 7
	    istringstream(line) >> a >> b; // MULTI-SPECIES TO BE IMPLEMENTED
	    N = stoi(a);
	    if(debug) cout << "  Line 7: number of species 1 = " << N << endl;
    }
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
      if(constantBox) nframes = (nlines - 7) / (1+N);
      else            nframes = nlines / (8+N); // 8 system info + N position lines
    }
    for(auto &p: ps) {
      p.read_3cols(i); // N particle lines
      // from direct to cartesian
      if(format==0) p.r = box*p.r; // x' = (x*ax + y*bx + z*cx) etc.
      // already cartesian: must only be scaled by s
      else if(format==1) p.r*=s;
    }
}

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
read_alphanes_frame(fstream &i, bool resetN)
{
    string line;
    stringstream ss;
    ntype x;
    int ncols=0, natom, nxyz;

    getline(i,line); // all frame (box+positions) is contained in one line
    if(debug) cout << "\n  Line: " << line << endl;

    ss << line;
    while( ss >> x )
    {
      if(ncols==0)      box[0][0] = x; // ax
      else if(ncols==1) box[0][1] = x; // bx
      else if(ncols==2) box[0][2] = x; // cx
      else if(ncols==3) box[1][1] = x; // by
      else if(ncols==4) box[1][2] = x; // cy
      else if(ncols==5) box[2][2] = x; // cz
      else if(!resetN) { // assign positions only after they are allocated
        natom=int((ncols-6)/3);
        nxyz=(ncols-6)%3;
        ps[natom].r[nxyz] = x;
      }
      ncols++;
    }
    box[1][0]=box[2][0]=box[2][1] = 0.0; //  ay=az=bz=0.0
    natom=int((ncols-6)/3);
    N = natom;

    if(resetN) {
      ps.resize(N);
      invN = 1.0/N;
      nframes = nlines;
    }
}

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
read_alphanes9_frame(fstream &i, bool resetN)
{
    string line;
    stringstream ss;
    ntype x;
    int ncols=0, natom, nxyz;

    getline(i,line); // all frame (box+positions) is contained in one line
    if(debug) cout << "\n  Line: " << line << endl;

    ss << line;
    while( ss >> x )
    {
      if(ncols==0)      box[0][0] = x; // ax
      else if(ncols==1) box[0][1] = x; // bx
      else if(ncols==2) box[0][2] = x; // cx
      else if(ncols==3) box[1][0] = x; // ay
      else if(ncols==4) box[1][1] = x; // by
      else if(ncols==5) box[1][2] = x; // cy
      else if(ncols==6) box[2][0] = x; // az
      else if(ncols==7) box[2][1] = x; // bz
      else if(ncols==8) box[2][2] = x; // cz
      else if(!resetN) { // assign positions only after they are allocated
        natom=int((ncols-9)/3);
        nxyz=(ncols-9)%3;
        ps[natom].r[nxyz] = x;
      }
      ncols++;
    }
    natom=int((ncols-9)/3);
    N = natom;

    if(resetN) {
      ps.resize(N);
      invN = 1.0/N;
      nframes = nlines;
    }
}

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
read_jmd_frame(fstream &i, bool resetN)
{
    string line, x;
    stringstream ss;
    int ncols=0;

    getline(i,line);
    if(debug) cout << "\n  Line 1: " << line << endl;

    ss << line;
    while( ss >> x )
    {
      if     (ncols==0) timestep = stoi(x);
      else if(ncols==1) N = stoi(x);
      else if(ncols==2) box[0][0] = stof(x); // Lx
      else if(ncols==3) box[1][1] = stof(x); // Ly
      else if(ncols==4) box[2][2] = stof(x); // Lz
      ncols++;
    }
    box[0][1]=box[0][2]=box[1][0]=box[1][2]=box[2][0]=box[2][1] = 0.0; // orthorombic
    set_L_from_box();

    if(resetN) {
      ps.resize(N);
      invN = 1.0/N;
      nframes = nlines / (N+1);
    }
    for(auto &p: ps) p.read_3cols(i); // N particle lines
}

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
read_lammpstrj_frame(fstream &i, bool resetN)
{
  enum class LAMMPS_ATOM_ENTRIES {
    ID, TYPE, X,Y,Z, XS,YS,ZS, IX,IY,IZ, VX,VY,VZ
  };
  string line, x, a,b,c;
  LAMMPS_ATOM_ENTRIES entries[15];
  stringstream ss;
  int ncols, num_atomic_entries, x_count, xs_count, pi_count, v_count;
  ntype xlo,ylo,zlo, xlob,ylob;
  ntype xhi,yhi,zhi, xhib,yhib;
  ntype xy,xz,yz;

  getline(i,line); // ITEM: TIMESTEP
  if(debug) cout << "\n  Line 1: " << line << endl;

  getline(i,line);
  timestep = stoi(line);
  if(debug) cout << timestep << endl;

  getline(i,line); // ITEM: NUMBER OF ATOMS
  if(debug) cout << "\n  Line 3: " << line << endl;
  getline(i,line);
  N = stoi(line);
  if(debug) cout << N << endl;

  // ITEM: BOX BOUNDS ...
  getline(i,line);
  ncols=0;
  if(debug) cout << "\n  Line 5: " << line << endl;
  ss << line;
  while( ss >> x )
  {
    if(ncols>=3) if(debug) cout << x << endl;
    ncols++;
  }
  ss.str(std::string()); ss.clear(); // clear the string stream!

  xy=xz=yz=0.0;
  switch (ncols) {
    case 6:
      getline(i,line); istringstream(line) >> a >> b;
      xlo = stof(a);
      xhi = stof(b);

      getline(i,line); istringstream(line) >> a >> b;
      ylo = stof(a);
      yhi = stof(b);

      getline(i,line); istringstream(line) >> a >> b;
      zlo = stof(a);
      zhi = stof(b);
      break;
    case 9:
      getline(i,line); istringstream(line) >> a >> b >> c;
      xlob = stof(a);
      xhib = stof(b);
      xy = stof(c);

      getline(i,line); istringstream(line) >> a >> b >> c;
      ylob = stof(a);
      yhib = stof(b);
      xz = stof(c);

      getline(i,line); istringstream(line) >> a >> b >> c;
      zlo = stof(a);
      zhi = stof(b);
      yz = stof(c);

      xlo = xlob - min(min(0.0,xy), min(xz,xy+xz) );
      xhi = xhib - max( max(0.0,xy), max(xz,xy+xz) );
      ylo = ylob - min(0.0,yz);
      yhi = yhib - max(0.0,yz);
      break;
    default:
      cout << "[ Error: box format not recognized with "<<ncols<<" columns in ITEM]\n\n";
      exit(1);
  }
  box[0][0] = xhi - xlo;
  box[1][1] = yhi - ylo;
  box[2][2] = zhi - zlo;
  box[0][1] = xy;
  box[0][2] = xz;
  box[1][2] = yz;
  box[1][0] = box[2][0] = box[2][1] = 0.0;
  set_L_from_box();

  if(resetN) {
    ps.resize(N);
    invN = 1.0/N;
    nframes = nlines / (N+9);
  }

  getline(i,line); // ITEM: ATOMS ...
  if(debug) cout << "\n  Line 9: " << line << endl;
  ncols=0;
  x_count=xs_count=pi_count=v_count=0;
  ss << line;
  while( ss >> x )
  {
    if(ncols==1 && x!="ATOMS")
    {
      cout << " [ERROR: I expected to find atoms but I found the following line:]\n"<<line<<endl;
      exit(1);
    }
    if(ncols>=2)
    {
      if(debug) cout << x << endl;
      if      (x=="id")   { entries[ncols-2]=LAMMPS_ATOM_ENTRIES::ID; }
      else if (x=="type") { entries[ncols-2]=LAMMPS_ATOM_ENTRIES::TYPE; }

      else if (x=="x") { entries[ncols-2]=LAMMPS_ATOM_ENTRIES::X; x_count++; }
      else if (x=="y") { entries[ncols-2]=LAMMPS_ATOM_ENTRIES::Y; x_count++; }
      else if (x=="z") { entries[ncols-2]=LAMMPS_ATOM_ENTRIES::Z; x_count++; }

      else if (x=="xs") { entries[ncols-2]=LAMMPS_ATOM_ENTRIES::XS; xs_count++; }
      else if (x=="ys") { entries[ncols-2]=LAMMPS_ATOM_ENTRIES::YS; xs_count++; }
      else if (x=="zs") { entries[ncols-2]=LAMMPS_ATOM_ENTRIES::ZS; xs_count++; }

      else if (x=="ix") { entries[ncols-2]=LAMMPS_ATOM_ENTRIES::IX; pi_count++; }
      else if (x=="iy") { entries[ncols-2]=LAMMPS_ATOM_ENTRIES::IY; pi_count++; }
      else if (x=="iz") { entries[ncols-2]=LAMMPS_ATOM_ENTRIES::IZ; pi_count++; }

      else if (x=="vx") { entries[ncols-2]=LAMMPS_ATOM_ENTRIES::VX; v_count++; }
      else if (x=="vy") { entries[ncols-2]=LAMMPS_ATOM_ENTRIES::VY; v_count++; }
      else if (x=="vz") { entries[ncols-2]=LAMMPS_ATOM_ENTRIES::VZ; v_count++; }

      else
      {
        cout << " [ERROR: lammps entry '"<<x<<"' for atoms not recognized.]\n";
        exit(1);
      }
    }
    ncols++;
  }
  ss.str(std::string()); ss.clear(); // clear the string stream!
  num_atomic_entries=ncols-2;

  for(auto &p: ps)
  {
    getline(i,line);
    ss << line;
    //if(debug) cout << line << endl;
    for(int j=0;j<num_atomic_entries;j++)
    {
      ss >> x;
      //if(debug) cout << "j="<<j<<" ; x="<<x << endl;
      switch( entries[j] )
      {
        case LAMMPS_ATOM_ENTRIES::ID : break; // atom number
        case LAMMPS_ATOM_ENTRIES::TYPE : break; // atom type
        case LAMMPS_ATOM_ENTRIES::X : p.r[0] = stof(x); break; // x cartesian
        case LAMMPS_ATOM_ENTRIES::Y : p.r[1] = stof(x); break; // y cartesian
        case LAMMPS_ATOM_ENTRIES::Z : p.r[2] = stof(x); break; // z cartesian
        case LAMMPS_ATOM_ENTRIES::XS : p.s[0] = stof(x); break; // x scaled
        case LAMMPS_ATOM_ENTRIES::YS : p.s[1] = stof(x); break; // y scaled
        case LAMMPS_ATOM_ENTRIES::ZS : p.s[2] = stof(x); break; // z scaled
        case LAMMPS_ATOM_ENTRIES::IX : p.pi[0] = stof(x); break; // x periodic image
        case LAMMPS_ATOM_ENTRIES::IY : p.pi[1] = stof(x); break; // y periodic image
        case LAMMPS_ATOM_ENTRIES::IZ : p.pi[2] = stof(x); break; // z periodic image
        case LAMMPS_ATOM_ENTRIES::VX : break;
        case LAMMPS_ATOM_ENTRIES::VY : break;
        case LAMMPS_ATOM_ENTRIES::VZ : break;
        default:
          cout << "[ERROR: lammps entry at column '"<<j<<"' not recognized while reading atom.]\n";
          exit(1);
      }
    }
    if(x_count)
    {
      p.r[0] -= xlo; // CORRETTO traslare per xl0,ylo,zlo tutto ???
      p.r[1] -= ylo;
      p.r[2] -= zlo;
    }
    if(x_count && !xs_count) p.s = boxInv * p.r; // compute scaled coordinate from Cartesian
    if(xs_count && !x_count) p.r = box * p.r; // compute Cartesian coordinate from scaled
    //if(debug) p.show();
    ss.str(std::string()); ss.clear(); // clear the string stream!
  }
}

//------------------------------------------------------------//
template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
read_yuhan_frame(fstream &i, bool resetN, bool isFirstFrame)
{
  string line, x, a,b,c;
  stringstream ss;
  int ncols;

  if(isFirstFrame) // if first frame:
  {
    getline(i,line); // CONFIG
    if(debug) cout << "\n  Line: " << line << endl;

    getline(i,line); // CONFIG details
    if(debug) cout << "\n  Line: " << line << endl;
  }

  getline(i,line); // timestep timestep_value N  
  if(debug) cout << "\n  Line: " << line << endl;
  ss << line;
  ncols=0;
  while( ss >> x )
  {
    //if(debug) cout << "x="<<x<< endl;
    switch(ncols)
    {
      case 0:
        if(x!="timestep") {
          cout<<"[ERROR: expected 'timestep' at beginning of line; instead: "<<x<<"]\n";
          exit(1);
        }
        break;
      case 1:
        timestep = stoi(x);
        if(debug) cout << "timestep="<<timestep<<endl;
        break;
      case 2:
        N = stoi(x);
        if(debug) cout << "N="<<N<<endl;
      case 3: break; // I don't know the meaning
      case 4: break; // I don't know the meaning
      case 5: break; // I don't know the meaning
      case 6: break; // I don't know the meaning
      default:
        cout<<"[ERROR: I did not expect more than 7 columns in Line 3]\n";
        exit(1);
    }
    ncols++;
  }
  ss.str(std::string()); ss.clear(); // clear the string stream!

  getline(i,line); // lattice vector 'a'
  istringstream(line) >> a >> b >> c;
  box[0][0] = stof(a);
  box[1][0] = stof(b);
  box[2][0] = stof(c);

  getline(i,line); // lattice vector 'ab
  istringstream(line) >> a >> b >> c;
  box[0][1] = stof(a);
  box[1][1] = stof(b);
  box[2][1] = stof(c);

  getline(i,line); // lattice vector 'c'
  istringstream(line) >> a >> b >> c;
  box[0][2] = stof(a);
  box[1][2] = stof(b);
  box[2][2] = stof(c);

  set_L_from_box();

  if(resetN) {
    ps.resize(N);
    invN = 1.0/N;
    nframes = (nlines-2) / (4*N+4);
  }

  for(auto &p: ps)
  {
    getline(i,line); // 1st line: type index ... ... ==> useless for notype
    getline(i,line); // 2nd line: x y z
    istringstream(line) >> a >> b >> c;
    p.r[0] = stof(a);
    p.r[1] = stof(b);
    p.r[2] = stof(c);

    getline(i,line); // 3rd line: vx vy vz
    getline(i,line); // 4th line: fx fy fz
  }
}