using namespace std;
//------------ Input reading and interaction --------------------//

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::print_usage(char argv0[])
{
  fprintf(stderr, "\nUsage: %s [-d -h -v] [-alphanes -alphanes9 -contcar -jmd -lammpstrj -poscar -xdatcar -xdatcarV -xyz -xyz_cp2k -yuhan]"
  				  " [-box1 -box3 .box6 -box9 -image_convention -remove_rot_dof] [-outxyz] [-adf -altbc -bo -cn -edq -l -msd -rdf -rmin -rmax -sq -sqt]"
				  " [-rcut -p1half -period] [-out_xyz -out_alphanes -pbc_out -fskip -tag -timings]\n", argv0);
}

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::print_summary()
  {
  fprintf(stderr, "\n#---------------------------------------------------------------------------------------#.\n");
  fprintf(stderr, "  Computes some statistical quantities over the MD trajectory of a multi-species system.\n");
  fprintf(stderr, "#---------------------------------------------------------------------------------------#.\n");
  fprintf(stderr, " !!! WARNING: maybe I forgot to implement some functions !!!\n");
  fprintf(stderr, "#---------------------------------------------------------------------------------------#.\n");
  fprintf(stderr, "\n -h/--help \t Print this summary.");
  fprintf(stderr, "\n -d/--debug \t Open in debug mode.");
  fprintf(stderr, "\n -v/--verbose \t Print a lot of outputs during execution.");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n INPUT FILES (only one of the following must be selected, followed by the file name):");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n -alphanes \t alpha_nes format. It expects a [paste -d ' ' box.dat pos.dat] file (one row for each timestep; 6 columns for box + 3N columns for coordinates).");
  fprintf(stderr, "\n -alphanes9 \t like -alphanes but with 9 columns for box.");
  fprintf(stderr, "\n -contcar \t Concatenation of CONTCAR files containing lattice, positions, velocities, lattice velocities and gear-predictor-corrector data.");
  fprintf(stderr, "\n -poscar \t Concatenation of POSCAR files.");
  fprintf(stderr, "\n -jmd \t John Russo's Molecular Dynamics format. It expects a [rm tmp; ls pos_* | sort -V | while read el; do cat $el >> tmp; done] file (first row: time N Lx Ly Lz; then N rows for coordinates; repeat).");
  fprintf(stderr, "\n -lammpstrj \t LAMMPS format. It expects the output of a 'dump atom' command.");
  fprintf(stderr, "\n -xdatcar \t XDATCAR format.");
  fprintf(stderr, "\n -xdatcarV \t XDATCAR format, with constant box.");
  fprintf(stderr, "\n -xyz \t .xyz format. Box size is supplied via -box. Particle labels must be 0,1,2,... .");
  fprintf(stderr, "\n -xyz_cp2k \t .xyz format from CP2K. Box size is supplied via -box. Particles must be ordered by label.");
  fprintf(stderr, "\n -yuhan \t Hi Yuhan, this is the format you gave me.");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n BOX (Note: it will be overwritten if present in the input file):");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n -box1 \t INPUT: L. Size of the cubic box [this is the default, with L=%.2f].",L[0]);
  fprintf(stderr, "\n -box3 \t INPUT: Lx,Ly,Lz. Sizes of the orthorombic box");
  fprintf(stderr, "\n -box6 \t INPUT: Ax,Bx,Cx,By,Cy,Cz. Components of an upper-diagonalized box");
  fprintf(stderr, "\n -box9 \t INPUT: Ax,Bx,Cx,Ay,By,Cy,Az,Bz,Cz. Components of a general box");
  fprintf(stderr, "\n -image_convention \t How many box images? [ONLY FOR g(r)]. 1 (Minimum Image), 0 (Cluster), -1 (All within the cutoff). [default 1].");
  fprintf(stderr, "\n -remove_rot_dof \t Rotate the positions and the box to upper-diagonalize it. [default don't].");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n STATISTICAL ANALYSIS:");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n -adf \t Compute Angular Distribution Function within the 1st shell. INPUT: bin_width (in terms of the cosine). OUTPUT: %s.{traj,ave}.", s_adf.c_str() );
  fprintf(stderr, "\n -altbc \t Compute Angular-Limited Three-Body Correlation. INPUT: bin_width rmin maxangle. Uses the given bin width for each dimension, with rmin <= bond length <= rcut1 and |180°- bond angle|<=maxangle. OUTPUT: %s.{traj,ave}.", s_altbc.c_str() );
  fprintf(stderr, "\n -bo \t Compute the bond order orientation (BOO) and correlation (BOC) parameters. Angular momentum is defined by the option -l.  OUTPUT: %s.l*.{dat,ave}, %s.l*.{dat,ave,local.ave,.xyz}, %s.l*.dat.", s_bondorient.c_str(), s_bondcorr.c_str(), s_nxtal.c_str());
  fprintf(stderr, "\n -cn \t Compute the coordination number, i.e., the number of neighbours in the 1st shell, weighted by a cutoff function. OUTPUT: %s.{dat,ave}.", s_coordnum.c_str());
  fprintf(stderr, "\n -edq \t Compute the Eddington-Debenedetti 'q' bond order parameter.  OUTPUT: %s.{dat,ave,_classes.dat}.", s_edq.c_str() );
  fprintf(stderr, "\n -msd \t Compute the Mean Squared Displacement and the Non-Gaussianity Parameter. OUTPUT: %s.{traj,ave,ngp}.", s_msd.c_str() );
  fprintf(stderr, "\n -rdf \t Compute the Radial Distribution Function g(r). INPUT: bin_width, max_distance. OUTPUT: %s.{traj,ave}.", s_rdf.c_str() );
  fprintf(stderr, "\n -rmin \t Compute the minimum distance between atoms. OUTPUT: %s.dat.", s_rmin.c_str() );
  fprintf(stderr, "\n -rmax \t Compute the maximum distance between atoms. OUTPUT: %s.dat.", s_rmax.c_str() );
  fprintf(stderr, "\n -sq \t Compute the Static Structure Factor S(q). ONLY FOR CUBIC BOXES. INPUT: q_mod_min q_mod_max q_mod_step. OUTPUT: %s.{traj,ave}. [default %d %d %d]", s_sq.c_str(), qmodmin,qmodmax,qmodstep );
  fprintf(stderr, "\n -sqt \t Compute the Dynamic Structure Factor S(q,t). ONLY FOR CUBIC BOXES. INPUT: q_mod_min q_mod_max q_mod_step. OUTPUT: %s.{traj,ave}. [default %d %d %d]", s_sqt.c_str(), qmodmin,qmodmax,qmodstep);
  fprintf(stderr, "\n");
  fprintf(stderr, "\n -l \t Angular momentum for the computed bond order parameters [default %d].", l);
  fprintf(stderr, "\n -rcut \t File containing %d blocks of cutoff radii for each pair of atom types, in a symmetric matrix form. They will be used for cutoff functions in 1st,2nd,3rd shells [default %.2f,%.2f,%.2f for every pair].", MAX_NSHELL, cutoff[0][0], cutoff[1][0], cutoff[2][0]);
  fprintf(stderr, "\n -p1half \t Half the power for the radial cutoff function f(x) = (1-x^p1)/(1-x^p2) with p2=2*p1, p1=2*p1half. Must be integer [default %d].", p1half);
  fprintf(stderr, "\n -period \t Average over initial time t0 every 'period' (in timesteps units) when computing MSD and S(q,t). If negative, don't. [default %d].", period);
  fprintf(stderr, "\n");
  fprintf(stderr, "\n OTHER PARAMETERS:");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n -out_xyz \t Produces an output traj.xyz file.");
  fprintf(stderr, "\n -out_alphanes \t [TO BE COMPLETED] Produce the following self-explaining files: type.dat, box.dat, pos.dat. Box is rotated if -remove_rot_dof is activated. No tag is addded.");
  fprintf(stderr, "\n -pbc_out \t Apply PBC to the output trajectory files. [default don't].");
  fprintf(stderr, "\n -fskip \t Skip the given fraction of frames from beginning and from end. INPUT: fskip_from_beginning fskip_from_end. [default: 0.0 0.0].");
  fprintf(stderr, "\n -tag \t Add this text tag inside output files' name [default none].");
  fprintf(stderr, "\n -timings \t Measure time for computation of S(q) into the log file: %s[default don't].", s_log.c_str());
  fprintf(stderr, "\n");
  fprintf(stderr, "\n TIPS:");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n - Convert a .traj file to a column file with mdtraj/shell/traj2nxy.sh; then plot it with mdtraj/*/python/{rdf,msd,...}.traj.py");
  fprintf(stderr, "\n - Plot the ALTBC with mdtraj/*/python/plot.altbc.py");
  fprintf(stderr, "\n\n");
  }

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::args(int argc, char** argv)
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
	  else if ( !strcmp(argv[i], "-adf") )
	    {
	      c_adf = true;
	      i++;
	      if (i == argc) { fprintf(stderr, "ERROR: '-adf' must be followed by bin width!\n"); exit(-1); }
	      adf_binw = atof(argv[i]);
        if(adf_binw<=0 || adf_binw>=1) { fprintf(stderr, "ERROR: adf bin width value must be in (0,1)!\n"); exit(1); }
	    }
	  else if ( !strcmp(argv[i], "-bo") )
	      c_bondorient = true;
	  else if ( !strcmp(argv[i], "-cn") )
	      c_coordnum = true;
	  else if ( !strcmp(argv[i], "-edq") )
	      c_edq = true;
	  else if ( !strcmp(argv[i], "-box1") )
	    {
		    i++;
  		  if (i == argc){
  		    fprintf(stderr, "ERROR: '-box1' must be followed by the box size!\n");
  		    exit(-1);
  		  }
	      for(int j=0;j<3;j++) L[j] = atof(argv[i]);
		    set_box_from_L();
	    }
	  else if ( !strcmp(argv[i], "-box3") )
	    {
	      for(int j=0;j<3;j++){
		      i++;
		      if (i == argc){
			      fprintf(stderr, "ERROR: '-box3' must be followed by the 3 box sizes!\n");
			      exit(-1);
			    }
		      L[j] = atof(argv[i]);
		    }
		    set_box_from_L();
	    }
	  else if ( !strcmp(argv[i], "-box6") )
	    {
	      for(int j=0;j<6;j++){
		      i++;
		      if (i == argc){
			      fprintf(stderr, "ERROR: '-box6' must be followed by the 6 box sizes!\n");
			      exit(-1);
  			  }
  			  switch(j)
  			  {
  			    case 0: box[0][0] = atof(argv[i]); break;
  			    case 1: box[0][1] = atof(argv[i]); break;
  			    case 2: box[0][2] = atof(argv[i]); break;
  			    case 3: box[1][1] = atof(argv[i]); break;
  			    case 4: box[1][2] = atof(argv[i]); break;
  			    case 5: box[2][2] = atof(argv[i]); break;
  			  }
		    }
  			box[1][0]=box[2][0]=box[2][1]=0.0;
  			set_L_from_box();
	    }
	  else if ( !strcmp(argv[i], "-box9") )
	    {
	        for(int j=0;j<9;j++){
		      i++;
		      if (i == argc){
			    fprintf(stderr, "ERROR: '-box9' must be followed by the 6 box sizes!\n");
			    exit(-1);
			  }
			  switch(j)
			  {
			    case 0: box[0][0] = atof(argv[i]); break;
			    case 1: box[0][1] = atof(argv[i]); break;
			    case 2: box[0][2] = atof(argv[i]); break;
			    case 3: box[1][0] = atof(argv[i]); break;
			    case 4: box[1][1] = atof(argv[i]); break;
			    case 5: box[1][2] = atof(argv[i]); break;
			    case 6: box[2][0] = atof(argv[i]); break;
			    case 7: box[2][1] = atof(argv[i]); break;
			    case 8: box[2][2] = atof(argv[i]); break;
			  }
		    }
			set_L_from_box();
	    }
	  else if ( !strcmp(argv[i], "-image_convention") )
    {
      i++;
      if (i == argc) { fprintf(stderr, "ERROR: '-image_convention' must be followed by an integer value!\n"); exit(-1); }
	    image_convention = atoi(argv[i]);
      if (image_convention!=1 && image_convention!=0 && image_convention!=-1) {
        fprintf(stderr, "ERROR: '-image_convention' must be followed by an integer value among {1,0,-1}!\n");
        exit(-1);
      }
    }
    else if ( !strcmp(argv[i], "-remove_rot_dof") )
        remove_rot_dof = true;
	  else if ( !strcmp(argv[i], "-l") )
	    {
	      i++;
	      if (i == argc) { fprintf(stderr, "ERROR: '-l' must be followed by angular momentum value!\n"); exit(-1); }
	      l = atoi(argv[i]);
	    }
	  else if ( !strcmp(argv[i], "-msd") )
	      c_msd = true;
	  else if ( !strcmp(argv[i], "-period") )
	    {
	      i++;
	      if (i == argc) { fprintf(stderr, "ERROR: '-period' must be followed by a positive integer value!\n"); exit(-1); }
	      period = atoi(argv[i]);
	      if(period<=0) { fprintf(stderr, "ERROR: '-period must be a positive integer value!\n"); exit(-1); }
	    }
	  else if ( !strcmp(argv[i], "-p1half") )
	    {
	      i++;
	      if (i == argc) { fprintf(stderr, "ERROR: '-p1half must be followed by a positive integer value!\n"); exit(-1); }
	      p1half = atoi(argv[i]);
	      if(p1half<=0) { fprintf(stderr, "ERROR: '-p1half must be a positive integer value!\n"); exit(-1); }
	    }
	  else if ( !strcmp(argv[i], "-rcut") )
	  	{
	      i++;
	      if (i == argc)
		  {
		  	fprintf(stderr, "ERROR: '-rcut' must be followed by file name!\n");
		  	exit(-1);
		  }
	      s_rcut = string(argv[i]);
		}
	  else if ( !strcmp(argv[i], "-rdf") )
	    {
	      c_rdf = true;
	      i++;
	      if (i == argc) { fprintf(stderr, "ERROR: '-rdf' must be followed by (bin width) and (max distance)!\n"); exit(-1); }
	      rdf_binw = atof(argv[i]);
	      i++;
	      if (i == argc) { fprintf(stderr, "ERROR: '-rdf' must be followed by (bin width) and (max distance)!\n"); exit(-1); }
	      rdf_rmax = atof(argv[i]);
	    }
    else if ( !strcmp(argv[i], "-rmin") )
       c_rmin = true;
     else if ( !strcmp(argv[i], "-rmax") )
        c_rmax = true;
 	  else if ( !strcmp(argv[i], "-sq") )
 	    {
        c_sq = true;
 	      i++;
 	      if (i == argc) { fprintf(stderr, "ERROR: '-sq' must be followed by 3 integer values!\n"); exit(-1); }
 	      qmodmin = atoi(argv[i]);
 	      if (qmodmin<2) { fprintf(stderr, "ERROR: q_mod_min must be >=2 !\n"); exit(-1); }
 	      i++;
 	      if (i == argc) { fprintf(stderr, "ERROR: '-sq' must be followed by 3 integer values!\n"); exit(-1); }
 	      qmodmax = atoi(argv[i]);
 	      if (qmodmax<qmodmin || qmodmax>500) { fprintf(stderr, "ERROR: q_mod_max must be > q_mod_min and <=500 !\n"); exit(-1); }
 	      i++;
 	      if (i == argc) { fprintf(stderr, "ERROR: '-sq' must be followed by 3 integer values!\n"); exit(-1); }
 	      qmodstep = atoi(argv[i]);
 	      if (qmodstep<1) { fprintf(stderr, "ERROR: q_mod_step must be >=1 !\n"); exit(-1); }
 	    }

 	  else if ( !strcmp(argv[i], "-sqt") )
 	    {
        c_sqt = true;
 	      i++;
 	      if (i == argc) { fprintf(stderr, "ERROR: '-sqt' must be followed by 3 integer values!\n"); exit(-1); }
 	      qmodmin = atoi(argv[i]);
 	      if (qmodmin<2) { fprintf(stderr, "ERROR: q_mod_min must be >=2 !\n"); exit(-1); }
 	      i++;
 	      if (i == argc) { fprintf(stderr, "ERROR: '-sqt' must be followed by 3 integer values!\n"); exit(-1); }
 	      qmodmax = atoi(argv[i]);
 	      if (qmodmax<qmodmin || qmodmax>500) { fprintf(stderr, "ERROR: q_mod_max must be > q_mod_min and <=500 !\n"); exit(-1); }
 	      i++;
 	      if (i == argc) { fprintf(stderr, "ERROR: '-sqt' must be followed by 3 integer values!\n"); exit(-1); }
 	      qmodstep = atoi(argv[i]);
 	      if (qmodstep<1) { fprintf(stderr, "ERROR: q_mod_step must be >=1 !\n"); exit(-1); }
 	    }
	  else if ( !strcmp(argv[i], "-altbc") )
	    {
	      c_altbc = true;
	      i++;
	      if (i == argc) { fprintf(stderr, "ERROR: '-altbc' must be followed by [bin_width rmin maxangle]!\n"); exit(-1); }
	      altbc_binw = atof(argv[i]);
	      i++;
	      if (i == argc) { fprintf(stderr, "ERROR: '-altbc' must be followed by [bin_width rmin maxangle]!\n"); exit(-1); }
	      altbc_rmin = atof(argv[i]);
	      i++;
	      if (i == argc) { fprintf(stderr, "ERROR: '-altbc' must be followed by [bin_width rmin maxangle]!\n"); exit(-1); }
	      altbc_angle_th = atof(argv[i]);
	      if(altbc_angle_th<0.0 || altbc_angle_th>180.0){ fprintf(stderr, "ERROR: ALTBC angular limit must be within 0° and 180°!\n"); exit(1); }
	    }
	  else if ( !strcmp(argv[i], "-fskip") )
	    {
	      i++;
	      if (i == argc) { fprintf(stderr, "ERROR: '-fskip' must be followed by 2 values in [0.0,1.0) !\n"); exit(-1); }
	      fskip0 = atof(argv[i]);
	      if (fskip0<0.0 || fskip0>=1.0) { fprintf(stderr, "ERROR: '-fskip' must be followed by 2 values in [0.0,1.0) !\n"); exit(-1); }
	      i++;
	      if (i == argc) { fprintf(stderr, "ERROR: '-fskip' must be followed by 2 values in [0.0,1.0) !\n"); exit(-1); }
	      fskip1 = atof(argv[i]);
	      if (fskip1<0.0 || fskip1>=1.0) { fprintf(stderr, "ERROR: '-fskip' must be followed by 2 values in [0.0,1.0) !\n"); exit(-1); }
	    }
	  else if ( !strcmp(argv[i], "-tag") )
	    {
	      i++;
	      if (i == argc) { fprintf(stderr, "ERROR: '-tag' must be followed by some text!\n"); exit(-1); }
	      tag = argv[i];
	      tag.insert(0, 1, '.'); // add a dot to the beginning of the tag
	    }
	  else if ( !strcmp(argv[i], "-alphanes") )
	    {
	      i++;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: '-alphanes' must be followed by file name!\n");
		  exit(-1);
		}
          filetype = FileType::ALPHANES;
	      s_in = string(argv[i]);
	    }
	  else if ( !strcmp(argv[i], "-alphanes9") )
	    {
	      i++;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: '-alphanes9' must be followed by file name!\n");
		  exit(-1);
		}
          filetype = FileType::ALPHANES9;
	      s_in = string(argv[i]);
	    }
	  else if ( !strcmp(argv[i], "-contcar") )
    {
      i++;
      if (i == argc)
  		{
  		  fprintf(stderr, "ERROR: '-contcar' must be followed by file name!\n");
  		  exit(-1);
  		}
      filetype = FileType::CONTCAR;
	    s_in = string(argv[i]);
	  }

    else if ( !strcmp(argv[i], "-poscar") )
    {
      i++;
      if (i == argc)
  		{
  		  fprintf(stderr, "ERROR: '-poscar' must be followed by file name!\n");
  		  exit(-1);
  		}
      filetype = FileType::POSCAR;
	    s_in = string(argv[i]);
	  }

    else if ( !strcmp(argv[i], "-jmd") )
	  {
	    i++;
	    if (i == argc)
      {
    	  fprintf(stderr, "ERROR: '-jmd' must be followed by file name!\n");
    	  exit(-1);
    	}
      filetype = FileType::JMD;
	    s_in = string(argv[i]);
	  }

	  else if ( !strcmp(argv[i], "-lammpstrj") )
    {
	    i++;
	    if (i == argc)
    	{
    	  fprintf(stderr, "ERROR: '-lammpstrj' must be followed by file name!\n");
    	  exit(-1);
    	}
      filetype = FileType::LAMMPSTRJ;
	    s_in = string(argv[i]);
	  }

    else if ( !strcmp(argv[i], "-xdatcar") )
	  {
	    i++;
	    if (i == argc)
		  {
        fprintf(stderr, "ERROR: '-xdatcar' must be followed by file name!\n");
        exit(-1);
      }
      filetype = FileType::XDATCAR;
	    s_in = string(argv[i]);
	  }

    else if ( !strcmp(argv[i], "-xdatcarV") )
	    {
	      i++;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: '-xdatcarV' must be followed by file name!\n");
		  exit(-1);
		}
              filetype = FileType::XDATCARV;
	      s_in = string(argv[i]);
	    }
	  else if ( !strcmp(argv[i], "-xyz") )
	    {
	      i++;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: '-xyz' must be followed by file name!\n");
		  exit(-1);
		}
              filetype = FileType::XYZ;
	      s_in = string(argv[i]);
	    }
	  else if ( !strcmp(argv[i], "-xyz_cp2k") )
	    {
	      i++;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: '-xyz_cp2k' must be followed by file name!\n");
		  exit(-1);
		}
        filetype = FileType::XYZ_CP2K;
	      s_in = string(argv[i]);
	    }
	  else if ( !strcmp(argv[i], "-yuhan") )
	    {
	      i++;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: '-yuhan' must be followed by file name!\n");
		  exit(-1);
		}
        filetype = FileType::YUHAN;
	      s_in = string(argv[i]);
	    }
	  else if ( !strcmp(argv[i], "-out_xyz") )
	    out_xyz = true;
	  else if ( !strcmp(argv[i], "-out_alphanes") )
	  {
	    out_alphanes = true;
		remove_rot_dof = true; // important!
	  }
  	else if ( !strcmp(argv[i], "-pbc_out") )
    	pbc_out = true;
  	else if ( !strcmp(argv[i], "-timings") )
    	timings = true;
	  else
	    {
	      fprintf(stderr, "ERROR: Invalid argumet!\n");
	      exit(-1);
	    }
      	  i++;
	}
    }
}
