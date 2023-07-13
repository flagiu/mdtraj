#include "lib/mdtraj.hpp"
using namespace std;
//------------ Input reading and interaction --------------------//

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::print_usage(char argv0[])
{
  fprintf(stderr, "\nUsage: %s [-d -h -v] [-alphanes -contcar -jmd -xdatcar -xdatcarV -xyz -xyz_cp2k] [-box1 -box3 .box6 -box9 -general_box] [-outxyz] [-adf -altbc -bo -cn -l -msd -rdf -rmin] [-period -p1half -rcut1 -rcut2 -rcut3 -tag]\n", argv0);
}

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::print_summary()
  {
  fprintf(stderr, "\n#---------------------------------------------------------------------------------------#.\n");
  fprintf(stderr, "  Computes some statistical quantities over the MD trajectory of a mono-species system.\n");
  fprintf(stderr, "#---------------------------------------------------------------------------------------#.\n");
  fprintf(stderr, "\n -h/--help \t Print this summary.");
  fprintf(stderr, "\n -d/--debug \t Open in debug mode.");
  fprintf(stderr, "\n -v/--verbose \t Print a lot of outputs during execution.");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n INPUT FILES (only one of the following must be selected, followed by the file name):");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n -alphanes \t alpha_nes format. It expects a [paste -d ' ' box.dat pos.dat] file (one row for each timestep; 6 columns for box + 3N columns for coordinates).");
  fprintf(stderr, "\n -contcar \t Concatenation of CONTCAR files containing lattice, positions, velocities, lattice velocities and gear-predictor-corrector data.");
  fprintf(stderr, "\n -jmd \t John Russo's Molecular Dynamics format. It expects a [rm tmp; ls pos_* | sort -V | while read el; do cat $el >> tmp; done] file (first row: time N Lx Ly Lz; then N rows for coordinates; repeat).");
  fprintf(stderr, "\n -xdatcar \t XDATCAR format.");
  fprintf(stderr, "\n -xdatcarV \t XDATCAR format, with constant box.");
  fprintf(stderr, "\n -xyz \t .xyz format. Box size is supplied via -box.");
  fprintf(stderr, "\n -xyz_cp2k \t .xyz format from CP2K. Box size is supplied via -box.");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n BOX (Note: it will be overwritten if present in the input file):");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n -box1 \t L size of the cubic box [this is the default, with L=%.2f].",L[0]);
  fprintf(stderr, "\n -box3 \t Lx,Ly,Lz sizes of the orthorombic box");
  fprintf(stderr, "\n -box6 \t Ax,Bx,Cx,By,Cy,Cz components of an upper-diagonalized box");
  fprintf(stderr, "\n -box9 \t Ax,Bx,Cx,Ay,By,Cy,Az,Bz,Cz components of a general box");
  fprintf(stderr, "\n -general_box \t Keep all 9 box elements, instead of rotating the box to upper-diagonalize it.");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n STATISTICAL ANALYSIS:");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n -adf \t Compute Angular Distribution Function within the 1st shell. INPUT: nbins. OUTPUT: %s.{traj,ave}.", s_adf.c_str() );
  fprintf(stderr, "\n -bo \t Compute the bond order orientation (BOO) and correlation (BOC) parameters. Angular momentum is defined by the option -l.  OUTPUT: %s.l*.{dat,ave}, %s.l*.{dat,ave,local.ave,.xyz}, %s.l*.dat.", s_bondorient.c_str(), s_bondcorr.c_str(), s_nxtal.c_str());
  fprintf(stderr, "\n -cn \t Compute the coordination number, i.e., the number of neighbours in the 1st shell. OUTPUT: %s.{dat,ave}.", s_coordnum.c_str());
  fprintf(stderr, "\n -msd \t Compute Mean Squared Displacement and Non-Gaussianity Parameter. OUTPUT: %s.{traj,ave,ngp}.", s_msd.c_str() );
  fprintf(stderr, "\n -rdf \t Compute Radial Distribution Function. INPUT: nbins. OUTPUT: %s.{traj,ave}.", s_rdf.c_str() );
  fprintf(stderr, "\n -altbc \t Compute Angular-Limited Three-Body Correlation. INPUT: nbins rmin maxangle. Uses the given number of bins for each dimension, with tmin <= bond length <= rcut1 and |180°- bond angle|<=maxangle. OUTPUT: %s.{traj,ave}.", s_altbc.c_str() );
  fprintf(stderr, "\n -rmin \t Compute the minimum distance between atoms. OUTPUT: %s.dat.", s_rmin.c_str() );
  fprintf(stderr, "\n");
  fprintf(stderr, "\n OTHER PARAMETERS:");
  fprintf(stderr, "\n");
  fprintf(stderr, "\n -l \t Angular momentum for the computed bond order parameters [default %d].", l);
  fprintf(stderr, "\n -outxyz \t Print an output traj.xyz file. [default don't]");
  fprintf(stderr, "\n -period \t Average over initial time t0 every 'period' (in timesteps units) when computing MSD. If negative, don't. [default %d].", period);
  fprintf(stderr, "\n -p1half \t Half the power for the radial cutoff function f(x) = (1-x^p1)/(1-x^p2) with p2=2*p1, p1=2*p1half. Must be integer [default %d].", p1half);
  fprintf(stderr, "\n -rcut1 \t Cutoff radius for cutoff functions in 1st shell [default %.2f].", cutoff[0]);
  fprintf(stderr, "\n -rcut2 \t Cutoff radius for cutoff functions in 2nd shell [default %.2f].", cutoff[1]);
  fprintf(stderr, "\n -rcut3 \t Cutoff radius for cutoff functions in 3rd shell [default %.2f].", cutoff[2]);
  fprintf(stderr, "\n -tag \t Add this text tag inside output files' name [default none].");
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
	      if (i == argc) { fprintf(stderr, "ERROR: '-adf' must be followed by number of bins!\n"); exit(-1); }
	      adf_nbins = atoi(argv[i]);
	      if(adf_nbins < 2){ fprintf(stderr, "ERROR: too few bins for ADF!\n"); exit(1); }
	    }
	  else if ( !strcmp(argv[i], "-bo") )
	      c_bondorient = true;
	  else if ( !strcmp(argv[i], "-cn") )
	      c_coordnum = true;
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
	  else if ( !strcmp(argv[i], "-general_box") )
	      remove_rot_dof = false;
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
	  else if ( !strcmp(argv[i], "-rcut1") )
	    {
	      i++;
	      if (i == argc) { fprintf(stderr, "ERROR: '-rcut1' must be followed by cutoff value!\n"); exit(-1); }
	      cutoff[0] = atof(argv[i]);
	      if(cutoff[0]<=0) { fprintf(stderr, "ERROR: '-rcut1 must be a positive real value!\n"); exit(-1); }
	    }
	  else if ( !strcmp(argv[i], "-rcut2") )
	    {
	      i++;
	      if (i == argc) { fprintf(stderr, "ERROR: '-rcut2' must be followed by cutoff value!\n"); exit(-1); }
	      cutoff[1] = atof(argv[i]);
	      if(cutoff[1]<=0) { fprintf(stderr, "ERROR: '-rcut2 must be a positive real value!\n"); exit(-1); }
	    }
	  else if ( !strcmp(argv[i], "-rcut3") )
	    {
	      i++;
	      if (i == argc) { fprintf(stderr, "ERROR: '-rcut3' must be followed by cutoff value!\n"); exit(-1); }
	      cutoff[2] = atof(argv[i]);
	      if(cutoff[2]<=0) { fprintf(stderr, "ERROR: '-rcut3 must be a positive real value!\n"); exit(-1); }
	    }
	  else if ( !strcmp(argv[i], "-rdf") )
	    {
	      c_rdf = true;
	      i++;
	      if (i == argc) { fprintf(stderr, "ERROR: '-rdf' must be followed by number of bins!\n"); exit(-1); }
	      rdf_nbins = atoi(argv[i]);
	      if(rdf_nbins < 2){ fprintf(stderr, "ERROR: too few bins for RDF!\n"); exit(1); }
	    }
    else if ( !strcmp(argv[i], "-rmin") )
       c_rmin = true;
	  else if ( !strcmp(argv[i], "-altbc") )
	    {
	      c_altbc = true;
	      i++;
	      if (i == argc) { fprintf(stderr, "ERROR: '-altbc' must be followed by [nbins rmin maxangle]!\n"); exit(-1); }
	      altbc_nbins = atoi(argv[i]);
	      if(altbc_nbins < 2){ fprintf(stderr, "ERROR: too few bins for ALTBC!\n"); exit(1); }
	      i++;
	      if (i == argc) { fprintf(stderr, "ERROR: '-altbc' must be followed by [nbins rmin maxangle]!\n"); exit(-1); }
	      altbc_rmin = atof(argv[i]);
	      i++;
	      if (i == argc) { fprintf(stderr, "ERROR: '-altbc' must be followed by [nbins rmin maxangle]!\n"); exit(-1); }
	      altbc_angle = atof(argv[i]);
	      if(altbc_angle<0.0 || altbc_angle>180.0){ fprintf(stderr, "ERROR: ALTBC angular limit must be within 0° and 180°!\n"); exit(1); }
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
	  else if ( !strcmp(argv[i], "-outxyz") )
	    print_out_xyz = true;
	  else
	    {
	      fprintf(stderr, "ERROR: Invalid argumet!\n");
	      exit(-1);
	    }
      	  i++;
	}
    }
}
