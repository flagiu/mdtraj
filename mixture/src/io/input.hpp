using namespace std;
//------------ Reading input files of various formats --------------------//

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
read_frame(fstream &i, bool resetN, bool reset_nTypes, int frameIdx)
{
  if(resetN || reset_nTypes){
    nTypes=0;
    for(auto a=0;a<MAX_N_TYPES;a++) Nt[a]=0;
  }
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
    case FileType::POSCAR:
      read_poscar_frame(i, resetN);
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
      read_lammpstrj_frame(i, resetN, reset_nTypes);
      break;
    case FileType::YUHAN:
      read_yuhan_frame(i, resetN, frameIdx==0);
      break;
    case FileType::RUNNER:
      read_runner_frame(i, resetN);
      break;
    default:
      cerr << "[ERROR: filetype = " << static_cast<int>(filetype) << " not recognized]\n";
      exit(1);
  }
  if(remove_rot_dof) removeRotDof();
  boxInv = box.inverse();
  set_L_from_box();
  if(c_sq || c_sqt)
  {
    for(auto &p: ps) p.s = boxInv * p.r; // pre-compute fractional coordinates
  }
  if(debug) cerr << "  Read frame at timestep " << timestep << " DONE.\n";
}

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
removeRotDof()
{
  mat R1,R2,Rtot;
  vec u, xdir;
  ntype a_xy_radius, c, s, t, angle;
  if(debug) cerr <<"\n*------ Removing the 3 rotational degrees of freedom ------*\n";
  if(debug) cerr <<"  Pss: Remember that box=(a|b|c)\n";
  // Remove rounding errors
  if( abs(box[1][0])<1e-15 ) box[1][0]=0.0;
  if( abs(box[2][0])<1e-15 ) box[2][0]=0.0;
  if( abs(box[2][1])<1e-15 ) box[2][1]=0.0;
  if(debug) { cerr << "  Raw box = "; box.show(); }
  if( box[1][0]==box[2][0]==box[2][1]==0.0 ) {
    if(debug) {
      cerr <<"  Box has already the correct <=6 degrees of freedom!\n";
      cerr <<"*----------------------------------------------------------*\n";
    }
    return;
  }
  // let a, b, c be the columns of the box matrix.
  // 1) Align a to x-axis
  // u=(0,az,-ay) or (0,-az,ay) is the axis of rotation (to be normalized)

  if(debug) { cerr << "  Raw volume (with sign)= "<< box.det()<<endl; }
  u << 0.0, box[2][0], -1*box[1][0];
  if(debug) { cerr << "    axis1 = "; u.show(); }
  // c = cos(angle) = ax/|a| ; s = sin(angle)
  c = box[0][0] / sqrt( box[0][0]*box[0][0] + box[1][0]*box[1][0] + box[2][0]*box[2][0] );
  s = sqrt( 1.0 - c*c ); // + or - ? with the choice of u=(0,az,-ay) it should be +
  if(debug) { cerr << "    cos1 = "<< c << endl; }
  if(debug) { cerr << "    sin1 = "<< s << endl; }
  // build the rotation matrix
  R1 = rotation_matrix_axis_cossin( u, c, s ); // declared in lib/matrix.hpp
  if(box.det()<0) R1 *= -1.; // invert to go back to right-hand rule
  if(debug) { cerr << "  Rotation matrix R1 = "; R1.show(); }
  box = R1*box;
  if(debug) { cerr << "  box (after R1) = "; box.show(); }
  // 2) Rotate around the x-axis to bring the new b within the x-y plane
  xdir << 1.0, 0.0, 0.0;
  if(debug) { cerr << "    axis2 (x-axis) = "; xdir.show(); }
  // tan(angle) = -bz / by
  angle = -atan2( box[2][1], box[1][1] );
  if(debug) { cerr << "    angle2 (deg) = "<<angle*180/M_PI<<endl; }
  R2 = rotation_matrix_axis_cossin( xdir, cos(angle), sin(angle) );
  if(debug) { cerr << "  Rotation matrix R2 = "; R2.show(); }
  box = R2*box;
  if( abs(box[1][1])<1e-15 ) box[1][1]=0.0; // remove rounding errors
  if( abs(box[1][1])<1e-15 ) box[2][1]=0.0;
  if( abs(box[2][2])<1e-15 ) box[2][2]=0.0;
  if(debug) { cerr << "  box (after R2) = "; box.show(); }
  // 3) Apply both rotations to each particle
  Rtot=R2*R1;
  for(auto &p: ps) p.r = Rtot*p.r;
  // Final result: A=(Ax,0,0)^T, B=(Bx,By,0)^T, C generic ==> box is upper diagonal
  if(debug) { cerr <<"\n*------------------------------------------------------------*\n"; }
}

//-----------------------------------------------------------------------

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
read_xyz_frame(fstream &i, bool resetN)
{
  // any string labels, but must be ordered by type
  string line, a,b,c,d, cur_type_str;
  getline(i,line); // first line
  istringstream(line) >> a;
  N = stoi(a);
  if(debug) cerr << "\n  Line 1: " << N << " atoms\n";
  getline(i,line); // second line
  istringstream(line) >> b >> c >> d;
  timestep = stoi(d);
  if(debug) cerr << "  Line 2: Timestep " << timestep << endl;
  if(resetN) {
    ps.resize(N);
    nframes = nlines / (N+2);
  }
  int cur_type=0;
  for(int j=0;j<N;j++)
  {
    getline(i, line);
    istringstream(line) >> a >> b >> c >> d;
    //try { p.label = stoi(a); } except { cerr << "[Warning: could not convert type to integer.]\n"; }
    ps[j].r[0] = stof(b);
    ps[j].r[1] = stof(c);
    ps[j].r[2] = stof(d);

    if(j==0)
    {
      cur_type_str=a;
      type_names[cur_type]=cur_type_str;
    }
    else if(a!=cur_type_str)
    {
      if(debug) cerr << " Updating " << cur_type_str << "-->" << a << endl;
      cur_type_str=a;
      cur_type++;
      if(cur_type>=MAX_N_TYPES) {
        cerr << "[ Error: exceeded max number of allowed types ("<<MAX_N_TYPES<<"). Change MAX_N_TYPES and recompile if you need more. ]\n\n";
        exit(1);
      }
      type_names[cur_type]=cur_type_str;
    }
    ps[j].label = cur_type;
    if(resetN) Nt[ps[j].label]++;
  }

  if(resetN) nTypes = cur_type+1;
}

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
read_xyz_cp2k_frame(fstream &i, bool resetN)
{
  // assumes particles are ordered by type !!
  string line, a,b,c,d,e, cur_type_str;
  stringstream ss;
  int cur_type=0;
  getline(i,line); // first line
  istringstream(line) >> a;
  N = stoi(a);
  if(debug) cerr << "\n  Line 1: " << N << " atoms\n";
  getline(i,line); // second line
  ss << line;
  do {
    try {
      ss >> a;
      timestep = stoi(a);
      break;
    }
    catch (...) {
      if(debug) cerr << "a="<<a <<endl;
      continue;
    }
  } while(true);
  if(debug) cerr << "  Line 2: Timestep " << timestep << endl;
  if(resetN) {
    ps.resize(N);
    nframes = nlines / (N+2);
  }
  for(int j=0;j<N;j++)
  {
    getline(i, line);
    istringstream(line) >> a >> b >> c >> d;
    //if(debug) cerr << "Read particle: " << a << " @ " << b << " @ " << c << " @ " << d << endl;
    if(j==0)
    {
      cur_type_str=a;
      type_names[cur_type]=cur_type_str;
    }
    else if(a!=cur_type_str)
    {
      if(debug) cerr << " Updating " << cur_type_str << "-->" << a << endl;
      cur_type_str=a;
      cur_type++;
      if(cur_type>=MAX_N_TYPES) {
        cerr << "[ Error: exceeded max number of allowed types ("<<MAX_N_TYPES<<"). Change MAX_N_TYPES and recompile if you need more. ]\n\n";
        exit(1);
      }
      type_names[cur_type]=cur_type_str;
    }
    ps[j].label = cur_type;
    ps[j].r[0] = stof(b);
    ps[j].r[1] = stof(c);
    ps[j].r[2] = stof(d);
    if(resetN) Nt[cur_type]++;
  }
  if(resetN) nTypes = cur_type+1;
}

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
read_contcar_frame(fstream &i, bool resetN)
{
  string line, a,b,c,d;
  stringstream ss;
  ntype s;
  int format;

  getline(i,line); // line 1
  if(debug) cerr << "\n  Line 1: " << line << endl;

  getline(i,line); // line 2
  istringstream(line) >> a;
  s=stof(a);
  if(debug) cerr << "  Line 2: scaling factor = " << s << endl;

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
  if(debug) cerr << "  Line 3-5: lattice vectors\n";
  if(debug) box.show();
  if(debug) boxInv.show();

  getline(i,line); // line 6
  if(debug) cerr << "  Line 6: species = " << line << endl;
  if(resetN)
  {
    ss<<line;
    int cur_type=0;
    while(ss>>a)
    {
      type_names[cur_type]=a;
      cur_type++;
    }
    nTypes = cur_type;
    ss.str(std::string()); ss.clear(); // clear the string stream!
  }

  getline(i,line); // line 7
  if(debug) cerr << "  Line 7: number of atoms per species = " << line << endl;
  if(resetN)
  {
    ss<<line;
    N=0;
    for(int cur_type=0;cur_type<nTypes;cur_type++)
    {
      ss>>a;
      Nt[cur_type]=stoi(a);
      N += Nt[cur_type];
    }
  }

  getline(i,line); // line 8
  if(debug) cerr << "  Line 8: coordinates format = " << line << endl;
  if(line=="Direct" || line=="direct") format=0; // lattice units
  else if(line=="Cartesian" || line=="cartesian") format=1; // cartesian units
  else { cerr << "[ERROR: coordinates format not recognized: " << line << "]\n"; exit(1); }

  if(resetN) {
    ps.resize(N);
    nframes = nlines / (8+N+ 9+N+ 4+3*N); // N positions, N velocities, 3*N predictor-corrector
  }
  int cur_type=0;
  int cumulative_Nt=Nt[cur_type];
  for(int j=0;j<N;j++)
  {
    if(j==cumulative_Nt)
    {
      cur_type++;
      cumulative_Nt+=Nt[cur_type];
    }
    ps[j].read_3cols(i); // N particle lines
    ps[j].label = cur_type;
    // from direct to cartesian
    if(format==0) ps[j].r = box*ps[j].r; // x' = (x*ax + y*bx + z*cx) etc.
    // already cartesian: must only be scaled by s
    else if(format==1) ps[j].r*=s;
  }

  for(auto j=0;j<9+N;j++) getline(i,line); // skip velocities
  for(auto j=0;j<4+3*N;j++) getline(i,line); // skip predictor-corrector
}

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
read_poscar_frame(fstream &i, bool resetN)
{
  string line, a,b,c,d;
  stringstream ss;
  ntype s;
  int format;

  getline(i,line); // line 1
  if(debug) cerr << "\n  Line 1: " << line << endl;

  getline(i,line); // line 2
  istringstream(line) >> a;
  s=stof(a);
  if(debug) cerr << "  Line 2: scaling factor = " << s << endl;

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
  if(debug) cerr << "  Line 3-5: lattice vectors\n";
  if(debug) box.show();
  if(debug) boxInv.show();

  getline(i,line); // line 6
  if(debug) cerr << "  Line 6: species = " << line << endl;
  if(resetN)
  {
    ss<<line;
    int cur_type=0;
    while(ss>>a)
    {
      type_names[cur_type]=a;
      cur_type++;
    }
    nTypes = cur_type;
    ss.str(std::string()); ss.clear(); // clear the string stream!
  }

  getline(i,line); // line 7
  if(debug) cerr << "  Line 7: number of atoms per species = " << line << endl;
  if(resetN)
  {
    ss<<line;
    N=0;
    for(int cur_type=0;cur_type<nTypes;cur_type++)
    {
      ss>>a;
      Nt[cur_type]=stoi(a);
      N += Nt[cur_type];
    }
  }

  getline(i,line); // line 8
  if(debug) cerr << "  Line 8: coordinates format = " << line << endl;
  if(line=="Direct" || line=="direct") format=0; // lattice units
  else if(line=="Cartesian" || line=="cartesian") format=1; // cartesian units
  else { cerr << "[ERROR: coordinates format not recognized: " << line << "]\n"; exit(1); }

  if(resetN) {
    ps.resize(N);
    nframes = nlines / (8+N);
  }
  int cur_type=0;
  int cumulative_Nt=Nt[cur_type];
  for(int j=0;j<N;j++)
  {
    if(j==cumulative_Nt)
    {
      cur_type++;
      cumulative_Nt+=Nt[cur_type];
    }
    ps[j].read_3cols(i); // N particle lines
    ps[j].label = cur_type;
    // from direct to cartesian
    if(format==0) ps[j].r = box*ps[j].r; // x' = (x*ax + y*bx + z*cx) etc.
    // already cartesian: must only be scaled by s
    else if(format==1) ps[j].r*=s;
  }
}

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
read_xdatcar_frame(fstream &i, bool resetN, bool constantBox)
{
  string line, a,b,c,d;
  stringstream ss;
  ntype s;
  int format;
  if(resetN || !constantBox) {
    getline(i,line); // line 1
    if(debug) cerr << "\n  Line 1: " << line << endl;

    getline(i,line); // line 2
    istringstream(line) >> a;
    s=stof(a);
    if(debug) cerr << "  Line 2: scaling factor = " << s << endl;

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
    if(debug) cerr << "  Line 3-5: lattice vectors\n";

    getline(i,line); // line 6
    if(debug) cerr << "  Line 6: species = " << line << endl;
    if(resetN)
    {
      ss<<line;
      int cur_type=0;
      while(ss>>a)
      {
        type_names[cur_type]=a;
        cur_type++;
      }
      nTypes = cur_type;
      ss.str(std::string()); ss.clear(); // clear the string stream!
    }

    getline(i,line); // line 7
    if(debug) cerr << "  Line 7: number of atoms per species = " << line << endl;
    if(resetN)
    {
      ss<<line;
      N=0;
      for(int cur_type=0;cur_type<nTypes;cur_type++)
      {
        ss>>a;
        Nt[cur_type]=stoi(a);
        N += Nt[cur_type];
      }
    }
  }

  getline(i,line); // line 8
  ss.str(std::string()); ss.clear(); // clear the string stream!
  ss<<line;
  // get first segment separed by a whitespace
  ss>>a;
  // then get 2nd segment separed by a '='
  for(int foo=0;foo<2;foo++){
    if(!getline(ss, b, '=')) {
      cerr << "ERROR during reading of timestep in a xdatcar/xdatcarV frame.\n";
      exit(1);
    }
  }
  timestep=stoi(b);
  if(debug) cerr << "  Line 8: coordinates format = " << a << " , timestep = " << timestep << endl;
  if(a=="Direct" || a=="direct") format=0; // lattice units
  else if(a=="Cartesian" || a=="cartesian") format=1; // cartesian units
  else { cerr << "[ERROR: coordinates format not recognized: " << a << "]\n"; exit(1); }

  if(resetN) {
    ps.resize(N);
    if(constantBox) nframes = (nlines - 7) / (1+N);
    else            nframes = nlines / (8+N); // 8 system info + N position lines
  }
  int cur_type=0;
  int cumulative_Nt=Nt[cur_type];
  for(int j=0;j<N;j++) {
    if(j==cumulative_Nt)
    {
      cur_type++;
      cumulative_Nt+=Nt[cur_type];
    }
    ps[j].read_3cols(i); // N particle lines read into ps[j].r
    ps[j].label = cur_type;
    // from direct to cartesian
    if(format==0) {
      ps[j].s = ps[j].r; // we did read the scaled coordinate
      ps[j].r = box*ps[j].s; // x' = (x*ax + y*bx + z*cx) etc.
    }
    // already cartesian: must only be scaled by s
    else if(format==1) {
      ps[j].r*=s;
      ps[j].s = boxInv*ps[j].r;
    }
    // coordinates are wrapped: you can tell nothing about the unwrapped ones
    ps[j].su=ps[j].s; ps[j].ru=ps[j].r;
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
  if(debug) cerr << "\n  Line: " << line << endl;

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
    nframes = nlines;
    nTypes=1;
    Nt[0]=N;
    cerr << "WARNING: only monospecies input is implemented for alphanes format.\n";
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
  if(debug) cerr << "\n  Line: " << line << endl;

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
    nframes = nlines;
    nTypes=1;
    Nt[0]=N;
  }
}

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
read_jmd_frame(fstream &i, bool resetN)
{
  string line, x;
  stringstream ss;
  int ncols, ncol;

  nTypes = 1; // ONLY MONO-SPECIES IMPLEMENTED!!!

  getline(i,line);
  if(debug) cerr << "\n  Line 1: " << line << endl;

  ss.str(std::string()); ss.clear(); // clear the string stream!
  ss << line;
  // count the number of columns
  ncols=0;
  while( ss >> x ) {ncols++;};
  ss.str(std::string()); ss.clear(); // clear the string stream!
  ss << line; // reset the string string stream

  if(ncols==5) // box is written as Lx Ly Lz
  {
    ncol=0;
    while( ss >> x )
    {
      if     (ncol==0) timestep = stoi(x);
      else if(ncol==1) N = stoi(x);
      else if(ncol==2) box[0][0] = stof(x); // Lx
      else if(ncol==3) box[1][1] = stof(x); // Ly
      else if(ncol==4) box[2][2] = stof(x); // Lz
      ncol++;
    }
    box[0][1]=box[0][2]=box[1][0]=box[1][2]=box[2][0]=box[2][1] = 0.0;
  }
  else if(ncols==8) // box is written as Ax Bx, Cx, By, Cy, Cz
  {
    ncol=0;
    while( ss >> x )
    {
      if     (ncol==0) timestep = stoi(x);
      else if(ncol==1) N = stoi(x);
      else if(ncol==2) box[0][0] = stof(x); // Ax
      else if(ncol==3) box[0][1] = stof(x); // Bx
      else if(ncol==4) box[0][2] = stof(x); // Cx
      else if(ncol==5) box[1][1] = stof(x); // By
      else if(ncol==6) box[1][2] = stof(x); // Cy
      else if(ncol==7) box[2][2] = stof(x); // Cz
      ncol++;
    }
    box[1][0]=box[2][0]=box[2][1] = 0.0; // Ay,Az,Bz=0
  }
  else if(ncols>0)
  {
    cerr << "ERROR: ncols="<<ncols<<" not accepted for Line 1 of JMD-formatted frame.\n";
    exit(1);
  }

  set_L_from_box();

  if(resetN) {
    ps.resize(N);
    nframes = nlines / (N+1);
    Nt[0]=N; // ONLY MONOSPECIES IMPLEMENTED
  }

  for(auto &p: ps) p.read_3cols(i); // N particle lines
  return;
}

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
read_lammpstrj_frame(fstream &i, bool resetN, bool reset_nTypes)
{
  enum class LAMMPS_ATOM_ENTRIES {
    ID, TYPE, X,Y,Z, XS,YS,ZS, IX,IY,IZ, XU,YU,ZU, VX,VY,VZ
  };
  // Assumes that particles are labelled as 1,2,...,ntypes
  // This must be mapped to our convention: 0,1,...,ntypes-1
  string line, x, a,b,c;
  LAMMPS_ATOM_ENTRIES entries[15];
  stringstream ss;
  int ncols, num_atomic_entries, x_count, xs_count, xu_count, pi_count, v_count;
  ntype xlo,ylo,zlo, xlob,ylob;
  ntype xhi,yhi,zhi, xhib,yhib;
  ntype xy,xz,yz;

  getline(i,line); // ITEM: TIMESTEP
  if(debug) cerr << "\n  Line 1: " << line << endl;

  getline(i,line);
  timestep = stoi(line);
  if(debug) cerr << timestep << endl;

  getline(i,line); // ITEM: NUMBER OF ATOMS
  if(debug) cerr << "\n  Line 3: " << line << endl;
  getline(i,line);
  N = stoi(line);
  if(debug) cerr << N << endl;

  // ITEM: BOX BOUNDS ...
  getline(i,line);
  ncols=0;
  if(debug) cerr << "\n  Line 5: " << line << endl;
  ss << line;
  while( ss >> x )
  {
    if(ncols>=3) if(debug) cerr << x << endl;
    ncols++;
  }
  ss.str(std::string()); ss.clear(); // clear the string stream!

  xy=xz=yz=0.0;
  switch (ncols) {
    case 3:
      cerr << "WARNING: dump file has non-periodic cubic box\n";
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
      cerr << "[ Error: box format not recognized with "<<ncols<<" columns in ITEM]\n\n";
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
    nframes = nlines / (N+9);
  }

  getline(i,line); // ITEM: ATOMS ...
  if(debug) cerr << "\n  Line 9: " << line << endl;
  ncols=0;
  x_count=xs_count=pi_count=xu_count=v_count=0;
  ss << line;
  while( ss >> x )
  {
    if(ncols==1 && x!="ATOMS")
    {
      cerr << " [ERROR: I expected to find atoms but I found the following line:]\n"<<line<<endl;
      exit(1);
    }
    if(ncols>=2)
    {
      if(debug) cerr << x << endl;
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

      else if (x=="xu") { entries[ncols-2]=LAMMPS_ATOM_ENTRIES::XU; xu_count++; }
      else if (x=="yu") { entries[ncols-2]=LAMMPS_ATOM_ENTRIES::YU; xu_count++; }
      else if (x=="zu") { entries[ncols-2]=LAMMPS_ATOM_ENTRIES::ZU; xu_count++; }

      else if (x=="vx") { entries[ncols-2]=LAMMPS_ATOM_ENTRIES::VX; v_count++; }
      else if (x=="vy") { entries[ncols-2]=LAMMPS_ATOM_ENTRIES::VY; v_count++; }
      else if (x=="vz") { entries[ncols-2]=LAMMPS_ATOM_ENTRIES::VZ; v_count++; }

      else
      {
        cerr << " [ERROR: lammps entry '"<<x<<"' for atoms not recognized.]\n";
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
    //if(debug) cerr << line << endl;
    for(int j=0;j<num_atomic_entries;j++)
    {
      ss >> x;
      //if(debug) cerr << "j="<<j<<" ; x="<<x << endl;
      switch( entries[j] )
      {
        case LAMMPS_ATOM_ENTRIES::ID : break; // atom number
        case LAMMPS_ATOM_ENTRIES::TYPE : p.label = stoi(x)-1; break; // atom type
        case LAMMPS_ATOM_ENTRIES::X : p.r[0] = stof(x); break; // x cartesian
        case LAMMPS_ATOM_ENTRIES::Y : p.r[1] = stof(x); break; // y cartesian
        case LAMMPS_ATOM_ENTRIES::Z : p.r[2] = stof(x); break; // z cartesian
        case LAMMPS_ATOM_ENTRIES::XS : p.s[0] = stof(x); break; // x scaled
        case LAMMPS_ATOM_ENTRIES::YS : p.s[1] = stof(x); break; // y scaled
        case LAMMPS_ATOM_ENTRIES::ZS : p.s[2] = stof(x); break; // z scaled
        case LAMMPS_ATOM_ENTRIES::IX : p.pi[0] = stof(x); break; // x periodic image
        case LAMMPS_ATOM_ENTRIES::IY : p.pi[1] = stof(x); break; // y periodic image
        case LAMMPS_ATOM_ENTRIES::IZ : p.pi[2] = stof(x); break; // z periodic image
        case LAMMPS_ATOM_ENTRIES::XU : p.ru[0] = stof(x); break; // x cartesian unwrapped
        case LAMMPS_ATOM_ENTRIES::YU : p.ru[1] = stof(x); break; // y cartesian unwrapped
        case LAMMPS_ATOM_ENTRIES::ZU : p.ru[2] = stof(x); break; // z cartesian unwrapped
        case LAMMPS_ATOM_ENTRIES::VX : break;
        case LAMMPS_ATOM_ENTRIES::VY : break;
        case LAMMPS_ATOM_ENTRIES::VZ : break;
        default:
          cerr << "[ERROR: lammps entry at column '"<<j<<"' not recognized while reading atom.]\n";
          exit(1);
      }
    }
    if(x_count>0)
    {
      p.r[0] -= xlo; // CORRETTO traslare per xl0,ylo,zlo tutto ???
      p.r[1] -= ylo;
      p.r[2] -= zlo;
    }
    if(xu_count>0)
    {
      p.ru[0] -= xlo; // CORRETTO traslare per xl0,ylo,zlo tutto ???
      p.ru[1] -= ylo;
      p.ru[2] -= zlo;
    }
    // given only the wrapped cartesian:
    if(x_count>0 && xs_count==0 && pi_count==0 && xu_count==0) {
      p.s = boxInv * p.r;
      // you can tell nothing about the unwrapped ones
      p.su=p.s; p.ru=p.r; p.pi<<0,0,0;
    }
    // given only the wrapped fractional:
    if(xs_count>0 && x_count==0 && pi_count==0 && xu_count==0) {
      p.r = box * p.s;
      // you can tell nothing about the unwrapped ones
      p.su=p.s; p.ru=p.r; p.pi<<0,0,0;
    }
    // given only the wrapped fractionals and the periodic images:
    if(xs_count>0 && pi_count>0 && x_count==0 && xu_count==0) {
      p.r = box*p.s;    // get wrapped cartesian
      p.su = p.s + p.pi; // get unwrapped fractional
      p.ru = box*p.su;   // get unwrapped cartesian
    }
    // given only unwrapped cartesian coordinates:
    if(xu_count>0 && x_count==0 && xs_count==0 && pi_count==0) {
      p.su = boxInv * p.ru; // get unwrapped fractional
      p.pi = round(p.su);   // get periodic image
      p.s  = p.su - p.pi;   // get wrapped fractional in (-0.5,0.5]
      p.r = box * p.s;      // get wrapped cartesian
    }
    //if(debug) p.show();
    ss.str(std::string()); ss.clear(); // clear the string stream!
    if(resetN || reset_nTypes)
    {
      nTypes = max(nTypes,p.label+1);
      if(nTypes>MAX_N_TYPES) {
        cerr << "[ Error: exceeded max number of allowed types ("<<MAX_N_TYPES<<"). Change MAX_N_TYPES and recompile if you need more. ]\n\n";
        exit(1);
      }
      Nt[p.label]++;
    }
  }
}

//------------------------------------------------------------//
template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
read_yuhan_frame(fstream &i, bool resetN, bool isFirstFrame)
{
  string line, x, a,b,c, cur_type_str;
  int cur_type;
  stringstream ss;
  int ncols;

  if(isFirstFrame) // if first frame:
  {
    getline(i,line); // CONFIG
    if(debug) cerr << "\n  Line: " << line << endl;

    getline(i,line); // CONFIG details
    if(debug) cerr << "\n  Line: " << line << endl;
  }

  getline(i,line); // timestep timestep_value N
  if(debug) cerr << "\n  Line: " << line << endl;
  ss << line;
  ncols=0;
  while( ss >> x )
  {
    //if(debug) cerr << "x="<<x<< endl;
    switch(ncols)
    {
      case 0:
        if(x!="timestep") {
          cerr<<"[ERROR: expected 'timestep' at beginning of line; instead: "<<x<<"]\n";
          exit(1);
        }
        break;
      case 1:
        timestep = stoi(x);
        if(debug) cerr << "timestep="<<timestep<<endl;
        break;
      case 2:
        N = stoi(x);
        if(debug) cerr << "N="<<N<<endl;
      case 3: break; // I don't know the meaning
      case 4: break; // I don't know the meaning
      case 5: break; // I don't know the meaning
      case 6: break; // I don't know the meaning
      default:
        cerr<<"[ERROR: I did not expect more than 7 columns in Line 3]\n";
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

  getline(i,line); // lattice vector 'b'
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
    nframes = (nlines-2) / (4*N+4);
  }

  cur_type=0;
  for(int j=0;j<N;j++)
  {
    getline(i,line); // 1st line: type index ... ... and other stuff
    ss << line;
    ss >> a;
    ps[j].label = cur_type;
    if(j==0)
    {
      cur_type_str=a;
      type_names[cur_type]=cur_type_str;
    }
    else if(a!=cur_type_str)
    {
      if(debug) cerr << " Updating " << cur_type_str << "-->" << a << endl;
      cur_type_str=a;
      cur_type++;
      if(cur_type>=MAX_N_TYPES) {
        cerr << "[ Error: exceeded max number of allowed types ("<<MAX_N_TYPES<<"). Change MAX_N_TYPES and recompile if you need more. ]\n\n";
        exit(1);
      }
      type_names[cur_type]=cur_type_str;
    }
    Nt[cur_type]++;

    ss.str(std::string()); ss.clear(); // clear the string stream!
    getline(i,line); // 2nd line: x y z
    istringstream(line) >> a >> b >> c;
    ps[j].r[0] = stof(a);
    ps[j].r[1] = stof(b);
    ps[j].r[2] = stof(c);

    getline(i,line); // 3rd line: vx vy vz
    getline(i,line); // 4th line: fx fy fz
  }

  if(resetN) nTypes = cur_type+1;
}

//------------------------------------------------------------//
template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
read_runner_frame(fstream &i, bool resetN)
{
  string line, x, a,b,c, cur_type_str;
  int cur_type;
  stringstream ss;
  int ncols;

  getline(i,line); // begin
  if(debug) cerr << "\n  Line: " << line << endl;

  getline(i,line); // comment
  if(debug) cerr << "\n  Line: " << line << endl;

  getline(i,line); // lattice vector 'a'
  istringstream(line) >> x >> a >> b >> c;
  box[0][0] = stof(a);
  box[1][0] = stof(b);
  box[2][0] = stof(c);

  getline(i,line); // lattice vector 'b'
  istringstream(line) >> x >> a >> b >> c;
  box[0][1] = stof(a);
  box[1][1] = stof(b);
  box[2][1] = stof(c);

  getline(i,line); // lattice vector 'c'
  istringstream(line) >> x >> a >> b >> c;
  box[0][2] = stof(a);
  box[1][2] = stof(b);
  box[2][2] = stof(c);

  set_L_from_box();

  if(resetN) {
    ps.resize(N);
    nframes = (nlines-2) / (4*N+4);
  }
  /*  to be completed */
  return;
}
