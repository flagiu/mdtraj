using namespace std;

//---- "init" functions: prepare arrays, averages and headers for output files
//---- "compute" functions: compute dynamic averages
//---- "print" functions: complete the calculations at the end of trajectory and print out

//------------------------------------- Mean Squared Displacement and Non-Gaussian Parameter ---------------------------------------//
template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
init_msd()
{
      if(period>0 && period<nframes*dtframe){
        periodIdx = ceil(period/dtframe);
      } else{
        periodIdx = nframes;
      } // periodIdx = length of the trajectory, in index units
      Nperiods = nframes/periodIdx; // Nperiods = n. of trajectories with different t0
      if(debug) cout << "Set periodIdx = " << periodIdx << ", Nperiods = " << Nperiods << endl;
      rs.resize((N+1)*periodIdx); // store all coordinates at all times + center of mass at all time
      r2.resize(periodIdx-1); // < |r(t) - r(t0)|^2 >
      r4.resize(periodIdx-1); // < |r(t) - r(t0)|^4 >
      r2CM.resize(periodIdx-1); // < |rCM(t) - rCM(t0)|^2 >
      ss.str(std::string()); ss << s_msd << tag << ".traj"; fout.open(ss.str(), ios::out);
      fout << "#Delta timestep, then |r(t)-r(t0)|^2 # block file for different t0\n";
      fout.close();
      ss.str(std::string()); ss << s_msd << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "#Delta timestep, MSD, MSD error, MSD of c.o.m.\n";
      fout.close();
      ss.str(std::string()); ss << s_ngp << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "#Delta timestep, NGP (-0.4: homog. displ., 0.0: Brownian displ., >0: heterog. displ.)\n";
      fout.close();
}

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
compute_msd(int frameidx)
{
    int i, idx, idx_old, dframe, nperiod;
    vec dr;
    ntype dr2, r2t0;
    if(debug) cout << "*** MSD computation for timestep " << timestep << " STARTED ***\n";
    nperiod = frameidx / periodIdx;
    dframe = frameidx % periodIdx;
    ss.str(std::string()); ss << s_msd << tag << ".traj"; fout.open(ss.str(), ios::app);
    if(frameidx==0){
      for(i=0;i<periodIdx-1;i++){ // reset averages
        r2[i]=0.0;
        r4[i]=0.0;
        rs[ (N+1)*i + N ] << 0.0,0.0,0.0;
        r2CM[i]=0.0;
        fout << (i+1)*dtframe << endl;
      }
      rs[ (N+1)*(periodIdx-1) + N ] << 0.0,0.0,0.0;
      fout << endl; // end of the first block (delta timesteps)
    }
    if(dframe==0){
      r2t0 = 0.0; // average r^2 over atoms at fixed t,t0
      for(i=0;i<N;i++){
        rs[ (N+1)*0 + i] = ps[i].r; //store initial positions for this sample trajectory
        rs[ (N+1)*0 + N ] += ( ps[i].r * invN ); // store initial center of mass for '''
      }
    }else{
      for(i=0;i<N;i++)
      {
        idx = (N+1)*dframe + i;
        idx_old = idx - (N+1);       // previous frame
        dr = ps[i].r - rs[idx_old];  // displacement since previous frame
        //dr = dr.mic(L);              // correct for periodic boundary conditions
        dr = dr - box*round(boxInv*dr);
        rs[idx] = rs[idx_old] + dr;  // add correction
//        if(i==0) cout << timestep << " " << rs[idx][0]-rs[(N+1)*0+i][0] << endl;
        rs[ (N+1)*dframe + N ] += ( rs[idx] * invN ); // compute center of mass

        idx_old = (N+1)*0 + i;       // t0 for this sample trajectory
        dr = rs[idx] - rs[idx_old];  // displacement r(t)-r(t0)
        dr2 = dr.sq();
        r2t0 += dr2 * invN;          // sample trajectory
        r2[dframe-1] += dr2;         // average over t0 and atoms
        r4[dframe-1] += dr2*dr2;
      }
      fout << r2t0 << endl; // output sample trajectory
      if(dframe == periodIdx-1) fout << endl; // end of the block (sample trajectory)
      // after evaluating all particles, do the same for the center of mass
      dr = rs[ (N+1)*dframe + N ] - rs[ (N+1)*0 + N ];
      dr2 = dr.sq();
      r2CM[dframe-1] += dr2;
    }
    fout.close();
    if(debug) cout << "*** MSD computation for timestep " << timestep << " ENDED ***\n\n";
}

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
print_msd()
{
    ss.str(std::string()); ss << s_msd << tag << ".ave"; fout.open(ss.str(), ios::app);
    for(auto i=0; i<periodIdx-1; i++) {
      r2[i] /= (Nperiods*N);
      r4[i] /= (Nperiods*N);
      r2CM[i] /= Nperiods;
      fout << (i+1)*dtframe << " " << r2[i] << " " << sqrt( (r4[i] - r2[i]*r2[i])/(Nperiods*N) ) << " " << r2CM[i] << endl;
    }
    fout.close();
    if(debug) cout << "MSD printed to file\n";
    ss.str(std::string()); ss << s_ngp << tag << ".ave"; fout.open(ss.str(), ios::app);
    for(auto i=0; i<periodIdx-1; i++) {
      fout << (i+1)*dtframe << " " << 0.6*(r4[i]/(r2[i]*r2[i]))-1.0 << endl;
    }
    fout.close();
    if(debug) cout << "NGP printed to file\n";
}
