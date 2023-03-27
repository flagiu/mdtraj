#include "lib/mdtraj.hpp"
using namespace std;
//------------ Input reading and interaction --------------------//

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::print_usage(char argv0[])
{
  fprintf(stderr, "Usage: %s [-d -h -v] [-adf] [-bo] [-cn] [-in] [-l] [-L] [-msd] [-p/--period] [-p1half] [-rcut1] [-rcut2] [-rcut3] [-rdf] [-tag]\n", argv0);
}

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::print_summary()
  {
  const char *a_xyz = s_xyz.c_str();
  fprintf(stderr, "\nThis program computes some statistical quantities over a MD trajectory.\n");
  fprintf(stderr, "\n -h/--help \t Print this summary.\n -adf \t Compute Angular Distribution Function within the radial cutoff rcut1, using the given number of bins.\n -bo \t Compute the bond order orientation (BOO) and correlation (BOC) parameters. Angular momentum is defined by the option -l.\n -cn \t Compute the coordination number.\n -d/--debug \t Open in debug mode.\n -in/--input \t Input .xyz trajectory file [default: %s].\n -l \t Angular momentum for the computed bond order parameters [default %d].\n -L \t Lx,Ly,Lz sizes of the orthorombic supercell [default %.2f %.2f %.2f].\n -msd \t Compute Mean Squared Displacement.\n -p/--period \t Average over initial time t0 every 'period' (in timesteps units) when computing MSD. If negative, don't. [default %d].\n -p1half \t Half the power for the radial cutoff function f(x) = (1-x^p1)/(1-x^p2) with p2=2*p1, p1=2*p1half. Must be integer [default %d].\n -rcut1 \t Cutoff radius for cutoff functions in 1st shell [default %.2f].\n -rcut2 \t Cutoff radius for cutoff functions in 2nd shell [default %.2f].\n -rcut3 \t Cutoff radius for cutoff functions in 3rd shell [default %.2f].\n -rdf \t Compute Radial Distribution Function using the given number of bins.\n -tag \t Add this text tag inside output files' name [default none].\n -v/--verbose \t Print a lot of outputs during execution.\n\n", a_xyz, l, L[0],L[1],L[2], period, p1half, cutoff[0],cutoff[1],cutoff[2]);
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
	  else if ( !strcmp(argv[i], "-in") || !strcmp(argv[i],"--input") )
	    {
	      i++;
	      if (i == argc)
		{
		  fprintf(stderr, "ERROR: '-in' must be followed by file name!\n");
		  exit(-1);
		}
	      s_xyz = string(argv[i]);
	    }
	  else if ( !strcmp(argv[i], "-l") )
	    {
	      i++;
	      if (i == argc) { fprintf(stderr, "ERROR: '-l' must be followed by angular momentum value!\n"); exit(-1); }
	      l = atoi(argv[i]);
	    }
	  else if ( !strcmp(argv[i], "-L") )
	    {
	        for(int j=0;j<3;j++)
	            {
		      i++;
		      if (i == argc)
			{
			  fprintf(stderr, "ERROR: '-L' must be followed by the 3 box sizes!\n");
			  exit(-1);
			}
		      L[j] = atof(argv[i]);
		    }
	    }
	  else if ( !strcmp(argv[i], "-msd") )
	      c_msd = true;
	  else if ( !strcmp(argv[i], "-p") || !strcmp(argv[i], "--period") )
	    {
	      i++;
	      if (i == argc) { fprintf(stderr, "ERROR: '-p/--period' must be followed by a positive integer value!\n"); exit(-1); }
	      period = atoi(argv[i]);
	      if(period<=0) { fprintf(stderr, "ERROR: '-p/--period must be a positive integer value!\n"); exit(-1); }
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
	  else if ( !strcmp(argv[i], "-tag") )
	    {
	      i++;
	      if (i == argc) { fprintf(stderr, "ERROR: '-tag' must be followed by some text!\n"); exit(-1); }
	      tag = argv[i];
	      tag.insert(0, 1, '.'); // add a dot to the beginning of the tag
	    }
	  else
	    {
	      fprintf(stderr, "ERROR: Invalid argumet!\n");
	      exit(-1);
	    }
      	  i++;
	}
    }
}
