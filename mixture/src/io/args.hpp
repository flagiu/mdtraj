using namespace std;
//------------ Input reading and interaction --------------------//

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::print_usage(char argv0[])
{
  fprintf(stderr, "\nUsage: %s [-d -h -v] [-alphanes -alphanes9 -contcar -jmd -lammpsdata -lammpstrj -poscar -xdatcar -xdatcarV -xyz -xyz_cp2k -yuhan]"
  				  " [-box1 -box3 -box6 -box9 -image_convention -remove_rot_dof] [-adf -altbc -bo -cn -edq -l -msd -nna -nnd -oct -pmp -Qself -rdf -rmin -rmax -sq -sqt]"
				  " [-rcut -p1half -period] [ -dynamic_types -nodynamics -out_lammpsdump -out_xyz -pbc_out -fskip -tag -timings]\n", argv0);
}

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::print_summary()
  {
  fprintf(stderr, "\n#---------------------------------------------------------------------------------------#.\n");
  fprintf(stderr, "  Computes some statistical quantities over the MD trajectory of a multi-species system.\n");
  fprintf(stderr, "#---------------------------------------------------------------------------------------#.\n");
  fprintf(stderr, " !!! WARNING: maybe I forgot to convert some functions from monospecies to multispecies !!!\n");
  fprintf(stderr, "#---------------------------------------------------------------------------------------#.\n");
  fprintf(stderr, "\n -h/--help \t Print this summary.");
  fprintf(stderr, "\n -d/--debug \t Run in debug mode.");
  fprintf(stderr, "\n -v/--verbose \t Print a lot of outputs during execution (both on-screen-text, log file and data files).");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n INPUT FILES (only one of the following must be selected, followed by the file name):");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n -alphanes \t alpha_nes format. It expects a [paste -d ' ' box.dat pos.dat] file");
  fprintf(stderr, "\n           \t (one row for each timestep; 6 columns for box + 3N columns for coordinates). ONLY MONOSPECIES.");
  fprintf(stderr, "\n -alphanes9 \t like -alphanes but with 9 columns for box. ONLY MONOSPECIES.");
  fprintf(stderr, "\n -contcar \t Concatenation of CONTCAR files containing lattice, positions, velocities,");
  fprintf(stderr, "\n          \t lattice velocities and gear-predictor-corrector data.");
  fprintf(stderr, "\n -poscar \t Concatenation of POSCAR files.");
  fprintf(stderr, "\n -jmd \t John Russo's Molecular Dynamics format. It expects a ");
  fprintf(stderr, "\n      \t [rm tmp; ls pos_* | sort -V | while read el; do cat $el >> tmp; done]");
  fprintf(stderr, "\n      \t file (first row: time N Lx Ly Lz; then N rows for coordinates; repeat). ONLY MONOSPECIES.");
  fprintf(stderr, "\n -lammpsdata \t LAMMPS 'data' format in 'atomic' mode. It expects the output of a 'write_data' command.");
  fprintf(stderr, "\n -lammpstrj \t LAMMPS 'dump' format. It expects the output of a 'dump atom' command.");
  fprintf(stderr, "\n -xdatcar \t XDATCAR format.");
  fprintf(stderr, "\n -xdatcarV \t XDATCAR format, with constant box.");
  fprintf(stderr, "\n -xyz \t .xyz format. Box size is supplied via -box. Particle labels must be 0,1,2,... .");
  fprintf(stderr, "\n -xyz_cp2k \t .xyz format from CP2K. Box size is supplied via -box. Particles must be ordered by label.");
  fprintf(stderr, "\n -yuhan \t Hi Yuhan, this is the format you gave me.");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n BOX (Note: it will be overwritten if present in the input file):");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n -box1 \t INPUT: L. Size of a cubic box [this is the default, with L=%.2f].",L[0]);
  fprintf(stderr, "\n -box3 \t INPUT: Lx,Ly,Lz. Sizes of an orthorombic box");
  fprintf(stderr, "\n -box6 \t INPUT: Ax,Bx,Cx,By,Cy,Cz. Components of an upper-diagonalized box");
  fprintf(stderr, "\n -box9 \t INPUT: Ax,Bx,Cx,Ay,By,Cy,Az,Bz,Cz. Components of a general box");
  fprintf(stderr, "\n -image_convention \t How many box images? [ONLY FOR g(r)]. 1 (Minimum Image),");
  fprintf(stderr, "\n                   \t 0 (Cluster), -1 (All within the cutoff). [default 1].");
  fprintf(stderr, "\n -remove_rot_dof \t Rotate the positions and the box to upper-diagonalize it. [default don't].");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n STATISTICAL ANALYSIS:");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n -adf \t Compute Angular Distribution Function within the 1st sphere.");
  fprintf(stderr, "\n      \t INPUT: bin_width (in terms of the cosine). OUTPUT: %s.{traj,ave}.", s_adf.c_str() );
  fprintf(stderr, "\n -altbc \t Compute Angular-Limited Three-Body Correlation.");
  fprintf(stderr, "\n        \t INPUT: bin_width rmin maxangle. Uses the given bin width for each dimension,");
  fprintf(stderr, "\n        \t with rmin <= bond length <= rcut1 and |180°- bond angle|<=maxangle. OUTPUT: %s.{traj,ave}.", s_altbc.c_str() );
  fprintf(stderr, "\n -bo \t Compute the bond order orientation (BOO) and correlation (BOC) parameters.");
  fprintf(stderr, "\n     \t Defines crystalline particles by a threshold of %.2f on the BOC.",qldot_th);
  fprintf(stderr, "\n     \t Orbital angular momentum is provided by the option -l.");
  fprintf(stderr, "\n     \t OUTPUT: %s.l*.{dat,ave}, %s.l*.{dat,ave,local.ave}, %s.l*.{indexes,dat}.", s_bondorient.c_str(), s_bondcorr.c_str(), s_nxtal.c_str());
  fprintf(stderr, "\n -clusters \t Clusterize the particles if qdot>qdot_th AND distance<cutoff. Must be followed by a file where cutoff is specified: one row; one column for each type pair. OUTPUT: %s{.dat,_size.dat}.", s_clusters.c_str());
  fprintf(stderr, "\n -cn \t Compute the coordination number, i.e., the number of neighbours in the 1st sphere, weighted by a cutoff function. OUTPUT: %s.{dat,ave}.", s_coordnum.c_str());
  fprintf(stderr, "\n -edq \t Compute the Errington-Debenedetti 'q' bond order parameter.  OUTPUT: %s.{dat,ave,_classes.dat}.", s_edq.c_str() );
  fprintf(stderr, "\n -msd \t Compute the Mean Squared Displacement and the Non-Gaussianity Parameter. OUTPUT: %s.{traj,ave,ngp}.", s_msd.c_str() );
  fprintf(stderr, "\n -nna \t Compute Nearest Neighbour Angles for neighbours in the 1st sphere. INPUT: max_number_of_neighbours OUTPUT: %s.{dat,ave}.", s_nna.c_str());
  fprintf(stderr, "\n -nnd \t Compute Nearest Neighbour Distances for neighbours in the 1st sphere. INPUT: max_number_of_neighbours OUTPUT: %s.{dat,ave}.", s_nnd.c_str());
  fprintf(stderr, "\n -oct \t Compute a custom Octahedral order parameter on exactly 6 nearest neighbours");
  fprintf(stderr, "\n      \t q_oct(i)=SUM_{j!=k}[ 1/3 H(angle-135°)(cos(angle))^2 + 1/12 H(135°-angle)(cos(angle)+1)^2 ]");
  fprintf(stderr, "\n      \t with H=heavyside function.  OUTPUT: %s.dat.", s_oct.c_str() );
  fprintf(stderr, "\n -pmp \t Compute Pattern Matching Parameter 'q_tetrahedral' and 'q_octahedral' https://doi.org/10.3389/fmats.2017.00034 .  OUTPUT: %s.{dat,ave}.", s_pmp.c_str() );
  fprintf(stderr, "\n -Qself \t Compute the self-overlap parameter Q_s(t) and its susceptibility (uses the same routine for MSDU). INPUT: cutoff. OUTPUT: %s.ave.", s_overlap.c_str() );
  fprintf(stderr, "\n -rdf \t Compute the Radial Distribution Function g(r). INPUT: bin_width, max_distance. OUTPUT: %s.{traj,ave}.", s_rdf.c_str() );
  fprintf(stderr, "\n -rmin \t Compute the minimum distance between atoms. OUTPUT: %s.dat.", s_rmin.c_str() );
  fprintf(stderr, "\n -rmax \t Compute the maximum distance between atoms. OUTPUT: %s.dat.", s_rmax.c_str() );
  fprintf(stderr, "\n -sq \t Compute the Static Structure Factor S(q). ONLY FOR CUBIC BOXES. INPUT: q_mod_min q_mod_max q_mod_step. OUTPUT: %s.{traj,ave}. [default %d %d %d]", s_sq.c_str(), qmodmin,qmodmax,qmodstep );
  fprintf(stderr, "\n -sqt \t Compute the Dynamic Structure Factor S(q,t). ONLY FOR CUBIC BOXES. CONSIDERS ALL ATOMS AS SAME TYPE. INPUT: q_mod_min q_mod_max q_mod_step. OUTPUT: %s.{traj,ave}. [default %d %d %d]", s_sqt.c_str(), qmodmin,qmodmax,qmodstep);
  fprintf(stderr, "\n");
  fprintf(stderr, "\n -l \t One or more (separed by spaces) Angular momentum values for the computation");
  fprintf(stderr, "\n    \t of bond order parameters (0<=l<=127). [default l=%d].", l);
  fprintf(stderr, "\n -rcut \t File containing <=%d lines of cutoff radii for each pair of atom types, each line ordered by type pair (e.g. for 3 types: r00 r01 r02 r11 r12 r22).", MAX_NSPHERE);
  fprintf(stderr, "\n       \t They will be used for cutoff functions in neighbour-spheres. No need to specify higher order spheres if not required for calculation. [default %.2f;%.2f;%.2f for every pair].", defaultCutoff[0][0], defaultCutoff[1][0], defaultCutoff[2][0]);
  fprintf(stderr, "\n -p1half \t Half the power for the radial cutoff function f(x) = (1-x^p1)/(1-x^p2) with p2=2*p1, p1=2*p1half. Must be integer [default %d].", p1half);
  fprintf(stderr, "\n -period \t Average over initial time t0 every 'period' (in timesteps units, not number of frames!),");
  fprintf(stderr, "\n         \t when computing MSD and S(q,t). If negative, don't. [default %d].", period);
  fprintf(stderr, "\n -qdot_th \t Threshold for the qdot bond-order parameter. A particle is crystalline if qdot>qdot_th. Must be in [0,1] [default %f].", qldot_th);
  fprintf(stderr, "\n");
  fprintf(stderr, "\n OTHER PARAMETERS:");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n -dynamic_types \t (Dangerous!) Allow for number of types and number of particles per type to change during simulation.");
  fprintf(stderr, "\n -fskip \t Skip the given fraction of frames from beginning and from end. ");
  fprintf(stderr, "\n        \t Note: if logtime is set, the skip refers to the time schedule");
  fprintf(stderr, "\n        \t and is internally converted into the trajectory's fraction of frames.");
  fprintf(stderr, "\n        \t INPUT: fskip_from_beginning fskip_from_end. [default: 0.0 0.0].");
  fprintf(stderr, "\n -ignore_double_frames \t Skip a frame if timestep is the same as before. [default: raise an error]");
  fprintf(stderr, "\n -logtime \t Input has log-linear timesteps specified by the following file,");
  fprintf(stderr, "\n          \t containing the full schedule. Note: trajectory timesteps must");
  fprintf(stderr, "\n          \t match at the beginning, but can be truncated at the end [default: %s]",s_logtime.c_str());
  fprintf(stderr, "\n -nodynamics \t Do not raise error if timesteps are not equally spaced.");
  fprintf(stderr, "\n             \t WARNING: dynamic averages will be junk.");
  fprintf(stderr, "\n -out_lammpsdump \t Produces an output *.dump file.");
  fprintf(stderr, "\n                 \t If combined with -clusters, the particle type will reflect crystallinity.");
  fprintf(stderr, "\n                 \t If combined with -oct, the per-atom variable q_oct will be printed next to x,y,z.");
  fprintf(stderr, "\n -out_xyz \t Produces an output traj.xyz file.");
  fprintf(stderr, "\n -out_alphanes \t [TO BE COMPLETED] Produce the following self-explaining");
  fprintf(stderr, "\n               \t files: type.dat, box.dat, pos.dat. Box is rotated if");
  fprintf(stderr, "\n               \t -remove_rot_dof is activated. No tag is addded.");
  fprintf(stderr, "\n -pbc_out \t Apply PBC to the output trajectory files. [default don't].");
  fprintf(stderr, "\n -tag \t Add this text tag inside output files' name [default none].");
  fprintf(stderr, "\n -timings \t Measure time for computation of S(q) into the ");
  fprintf(stderr, "\n          \t log file: %s. [default don't].", s_log.c_str());
  fprintf(stderr, "\n");
  fprintf(stderr, "\n TIPS:");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n - Convert a .traj file to a column file with mdtraj/shell/traj2nxy.sh;");
  fprintf(stderr, "\n     then plot it with mdtraj/*/python/plot_{rdf,adf,msd,...}_trajectory.py");
  fprintf(stderr, "\n - Plot a statistical average for all type pairs with");
  fprintf(stderr, "\n     mdtraj/.../python/plot_{rdf,adf,msd,...}_average.py");
  fprintf(stderr, "\n - Plot the ALTBC with mdtraj/.../python/plot_altbc.py");
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
    else if ( !strcmp(argv[i], "-clusters") )
	  	{
	      i++;
	      if (i == argc)
		  {
		  	fprintf(stderr, "ERROR: '-clusters' must be followed by file name!\n");
		  	exit(-1);
		  }
      c_clusters = true;
	    s_rcut_clusters = string(argv[i]);
		}
	  else if ( !strcmp(argv[i], "-cn") )
	      c_coordnum = true;
    else if ( !strcmp(argv[i], "-nna") )
      {
        c_nna = true;
        i++;
        if (i == argc){
          fprintf(stderr, "ERROR: '-nna' must be followed by the max_number_of_neighbours!\n");
          exit(-1);
        }
        max_num_nna = atoi(argv[i]);
      }
      else if ( !strcmp(argv[i], "-nnd") )
        {
          c_nnd = true;
          i++;
          if (i == argc){
            fprintf(stderr, "ERROR: '-nnd' must be followed by the max_number_of_neighbours!\n");
            exit(-1);
          }
          max_num_nnd = atoi(argv[i]);
        }
	  else if ( !strcmp(argv[i], "-edq") )
	      c_edq = true;
    else if ( !strcmp(argv[i], "-oct") )
	      c_oct = true;
    else if ( !strcmp(argv[i], "-pmp") )
	      c_pmp = true;
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
	      if (i == argc) { fprintf(stderr, "ERROR: '-l' must be followed by (at least one) angular momentum value!\n"); exit(-1); }
	      num_l=0;
        l_list[num_l++] = atoi(argv[i++]);
        while(num_l<MAX_N_ANGMOM && i<argc && is_positive_integer(argv[i])){
          l_list[num_l++] = atoi(argv[i++]);
        }
        i--;
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
    else if ( !strcmp(argv[i], "-qdot_th") )
	    {
	      i++;
	      if (i == argc) { fprintf(stderr, "ERROR: '-qdot_th must be followed by a number in [0,1]!\n"); exit(-1); }
	      qldot_th = atof(argv[i]);
	      if(qldot_th<0||qldot_th>1) { fprintf(stderr, "ERROR: '--qdot_th must be followed by a number in [0,1]!\n"); exit(-1); }
	    }
    else if ( !strcmp(argv[i], "-Qself") )
	    {
	      c_msd = true;
	      i++;
	      if (i == argc) { fprintf(stderr, "ERROR: '-Qself' must be followed by (cutoff)!\n"); exit(-1); }
	      Qoverlap_cutoff = atof(argv[i]);
        if(Qoverlap_cutoff<=0) { fprintf(stderr, "ERROR: '--Qself must be followed by a positive number!\n"); exit(-1); }
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
	      if (i == argc) { fprintf(stderr, "ERROR: '-tag' must be followed by some text! (can be an empty string '')\n"); exit(-1); }
	      tag = argv[i];
	      if(tag.length()>0) tag.insert(0, 1, '.'); // add a dot to the beginning of the tag, if not empty
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

    else if ( !strcmp(argv[i], "-lammpsdata") )
    {
	    i++;
	    if (i == argc)
    	{
    	  fprintf(stderr, "ERROR: '-lammpsdata' must be followed by file name!\n");
    	  exit(-1);
    	}
      filetype = FileType::LAMMPSDATA;
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
	      if (i == argc){
    		  fprintf(stderr, "ERROR: '-yuhan' must be followed by file name!\n");
    		  exit(-1);
    		}
        filetype = FileType::YUHAN;
	      s_in = string(argv[i]);
	    }
    else if ( !strcmp(argv[i], "-dynamic_types") )
        dynamic_types = true;
    else if ( !strcmp(argv[i], "-ignore_double_frames") )
    	ignore_double_frames = true;
    else if ( !strcmp(argv[i], "-logtime") )
    {
       i++;
       if(i==argc){
         fprintf(stderr, "ERROR: '-logtime' must be followed by file name!\n");
         exit(-1);
       }
       logtime = true;
       s_logtime = string(argv[i]);
     }
   else if ( !strcmp(argv[i], "-nodynamics") )
	    nodynamics = true;
	  else if ( !strcmp(argv[i], "-out_xyz") )
	    out_xyz = true;
	  else if ( !strcmp(argv[i], "-out_alphanes") )
	  {
	    out_alphanes = true;
		remove_rot_dof = true; // important!
	  }
    else if ( !strcmp(argv[i], "-out_lammpsdump") )
	  {
	    out_lammpsdump = true;
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
