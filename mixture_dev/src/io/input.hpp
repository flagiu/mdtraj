using namespace std;
//------------ Reading input files of various formats --------------------//
#define MAX_N_TYPES 5

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
read_frame(fstream &i, bool resetN)
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
      default:
        cout << "[ERROR: filetype = " << static_cast<int>(filetype) << " not recognized]\n";
        exit(1);
    }
    if(remove_rot_dof) removeRotDof();
    boxInv = box.inverse();
    set_L_from_box();
    if(c_sq)
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
    string line, a,b,c,d,e, cur_type_str;
    int cur_type=0, Ntype[5];
    for(auto a=0;a<5;a++) Ntype[a]=0;
    getline(i,line); // first line
    istringstream(line) >> a;
    N = stoi(a);
    if(debug) cout << "\n  Line 1: " << N << " atoms\n";
    getline(i,line); // second line
    istringstream(line) >> b >> c >> d >> e;
    timestep = stoi(d);
    if(debug) cout << "  Line 2: Timestep " << timestep << endl;
    if(resetN) {
      ps.resize(N);
      invN = 1.0/N;
      nframes = nlines / (N+2);
    }
    for(int j=0;j<N;j++)
    {
      getline(i, line);
      istringstream(line) >> a >> b >> c >> d;
      //if(debug) cout << "Read particle: " << a << " @ " << b << " @ " << c << " @ " << d << endl;
      if(j==0) cur_type_str=a;
      else if(a!=cur_type_str)
      {
        if(debug) cout << " Updating " << cur_type_str << "-->" << a << endl;
        cur_type_str=a;
        cur_type++;
      }
      ps[j].label = cur_type;
      ps[j].r[0] = stof(b);
      ps[j].r[1] = stof(c);
      ps[j].r[2] = stof(d);
      Ntype[cur_type]++;
    }
    if(debug) cout << " Found " << cur_type+1 << " types"<< endl;
    if(resetN)
    {
      nTypes = cur_type+1;
      Nt.resize(nTypes);
      for(auto a=0;a<nTypes;a++) {
        Nt[a]=Ntype[a];
        if(debug) cout << "   n. atoms of type " << a << " = " << Nt[a] << endl;
      }
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

    for(auto j=0;j<9+N;j++) getline(i,line); // velocities
    for(auto j=0;j<4+3*N;j++) getline(i,line); // predictor-corrector
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
      if(ncols==0) box[0][0] = x; // ax
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
  // Assumes that particles are labelled as 1,2,...,ntypes
  // This must be mapped to our convention: 0,1,...,ntypes-1
    string line, x, a,b,c;
    stringstream ss;
    int ncols=0, Nt_array[MAX_N_TYPES];;
    ntype xlo,ylo,zlo, xlob,ylob;
    ntype xhi,yhi,zhi, xhib,yhib;
    ntype xy,xz,yz;
    if(resetN)
    {
      for(int j=0;j<MAX_N_TYPES;j++) Nt_array[j]=0;
      nTypes=0;
    }

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
    ss << line;
    if(debug) cout << "\n  Line 5: " << line << endl;
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

    getline(i,line); // ITEM: ATOMS id type xs ys zs
    for(auto &p: ps)
    {
      getline(i,line);
      ss << line;
      ss >> x; // atom number
      ss >> x; p.label = stoi(x)-1; // atom type
      ss >> x; p.s[0] = stof(x); // x scaled
      ss >> x; p.s[1] = stof(x); // y scaled
      ss >> x; p.s[2] = stof(x); // z scaled
      p.r = box * p.s; // real coordinates
      p.r[0] -= xlo; // CORRETTO???
      p.r[1] -= ylo;
      p.r[2] -= zlo;
      //if(debug) p.show();
      ss.str(std::string()); ss.clear(); // clear the string stream!
      if(resetN)
      {
        nTypes = max(nTypes,p.label+1);
        Nt_array[p.label]++;
      }
    }
    if(debug) cout << " Found " << nTypes << " types"<< endl;
    if(resetN)
    {
      Nt.resize(nTypes);
      for(int j=0;j<nTypes;j++) {
        Nt[j]=Nt_array[j];
        if(debug) cout << "   n. atoms of type " << j+1 << " = " << Nt[j] << endl;
      }
    }
}
