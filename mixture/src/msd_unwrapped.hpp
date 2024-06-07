#ifndef _MSD_H_
#define _MSD_H_
using namespace std;

//---- "init" functions: prepare arrays, averages and headers for output files
//---- "compute" functions: compute dynamic averages
//---- "print" functions: complete the calculations at the end of trajectory and print out

//------------------------------------- Mean Squared Displacement and Non-Gaussian Parameter ---------------------------------------//

template <class ntype, class ptype>
class MSD_Calculator
{
  using vec=myvec<ntype,3>;
  using mat=mymatrix<ntype,3,3>;
  private:
    int period_in_dt_units, num_periods, N, nTypes;
    vector<int> Nt;
    vector<vec> particle_rvec; // collection of all 3D coordinates at all times
    vector<vec> CM_rvec; // collection of Center-of-Mass 3D coordinates at all times
    vecflex< vecflex<ntype> > r2; // < |r(t) - r(t0)|^2 > averaged over t0 and particles, for each type
    vecflex< vecflex<ntype> > r4; // < |r(t) - r(t0)|^4 > averaged over t0 and particles, for each type
    vecflex<ntype> r2CM; // < |rCM(t) - rCM(t0)|^2 > averaged over t0, for all types
    string string_out_msd, string_out_ngp, myName, tag;
    fstream fout;
    stringstream ss;
    bool debug,verbose;
  public:
    MSD_Calculator(){
      myName = "MSD and NGP";
    }
    virtual ~MSD_Calculator(){}

    void init(int dtframe, int nframes, int period_in_real_units, int N_, int nTypes_, int *Nt_, string string_out_msd_, string string_out_ngp_, string tag_, bool debug_, bool verbose_)
    {
      if(debug) cout << myName << " Initialization STARTED\n";
      N=N_;
      nTypes=nTypes_;
      Nt.resize(nTypes);
      for(int t=0;t<nTypes;t++) Nt[t]=Nt_[t]; // number of atoms for each type
      string_out_msd = string_out_msd_;
      string_out_ngp = string_out_ngp_;
      tag = tag_;
      debug=debug_;
      verbose=verbose_;
      if(period_in_real_units>0 && period_in_real_units<nframes*dtframe){
        period_in_dt_units = ceil(period_in_real_units/dtframe);
      } else if(period_in_real_units>0){
        cout << "WARNING: the input period for MSD calculation exceeds maximum time interval.\n";
        cout << "         I will set it equal to the maximum time interval.\n";
        period_in_dt_units = nframes;
      }
      else{
        period_in_dt_units = nframes;
      }
      num_periods = nframes/period_in_dt_units;
      if(debug) cout << myName<<": Set period_in_dt_units = " << period_in_dt_units << ", num_periods = " << num_periods << endl;

      particle_rvec.resize(N*period_in_dt_units);
      // index = N*(current_time_index) + particle_index,
      //   with 0 <= current_time_index < period_in_dt_units,
      //   with 0 <= particle_index < N.
      CM_rvec.resize(period_in_dt_units);
      r2.resize(nTypes);
      r4.resize(nTypes);
      for(int t=0;t<nTypes;t++){
        r2[t].resize(period_in_dt_units-1); // < |r(t0+t) - r(t0)|^2 > for type t
        r4[t].resize(period_in_dt_units-1); // < |r(t0+t) - r(t0)|^4 > for type t
      }
      r2CM.resize(period_in_dt_units-1); // < |rCM(t0+t) - rCM(t0)|^2 >

      if(verbose)
      {
        ss.str(std::string()); ss << string_out_msd << tag << ".traj"; fout.open(ss.str(), ios::out);
        fout << "#Delta timestep, then |r(t)-r(t0)|^2 # block file for different t0\n";
      }
      for(int i=0;i<period_in_dt_units-1;i++){ // reset averages
        for(int t=0;t<nTypes;t++){
          r2[t][i]=0.0;
          r4[t][i]=0.0;
        }
        CM_rvec[i] << 0.0,0.0,0.0;
        r2CM[i]=0.0;
        if(verbose){ fout << (i+1)*dtframe << endl; }
      }
      CM_rvec[period_in_dt_units-1] << 0.0,0.0,0.0; // last time was not counted
      if(verbose){
        fout << endl; // end of the first block (delta timesteps)
        fout.close();
      }
      ss.str(std::string()); ss << string_out_msd << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "#Delta timestep | MSD for each type | MSD error for each type | MSD of c.o.m.\n";
      fout.close();
      ss.str(std::string()); ss << string_out_ngp << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "#Delta timestep | NGP for each type (-0.4: homog. displ., 0.0: Brownian displ., >0: heterog. displ.)\n";
      fout.close();
      if(debug) cout << myName << " Initialization COMPLETED\n";
    }

    void compute(int frameidx, int timestep, vector<ptype> ps, PBC<ntype> *pbc)
    {
      int i, idx, idx_old, dframe, nperiod, num_images, t;
      const ntype invN=1.0/(ntype)N;
      vec dr, dr_image;
      ntype dr2, r2t0;
      if(debug) cout << "*** MSD computation for timestep " << timestep << " STARTED ***\n";
      nperiod = frameidx / period_in_dt_units;
      dframe = frameidx % period_in_dt_units;

      if(dframe==0)
      {
        for(i=0;i<N;i++)
        {
          particle_rvec[ N*0 + i] = ps[i].r; //store initial positions for this sample trajectory
          CM_rvec[0] += ( ps[i].r * invN ); // store initial center of mass for '''
        }
      }
      else
      {
        if(verbose) r2t0 = 0.0; // average r^2 over atoms at fixed t,t0
        for(i=0;i<N;i++)
        {
          idx = N*dframe + i;
          idx_old = idx - N;       // previous frame

          dr = ps[i].r - particle_rvec[idx_old];  // displacement since previous frame
          num_images = pbc->apply(dr, &dr_image); // apply periodic boundary conditions

          particle_rvec[idx] = particle_rvec[idx_old] + dr_image;  // add periodic correction
  //        if(i==0) cout << timestep << " " << rs[idx][0]-rs[(N+1)*0+i][0] << endl; //debug
          CM_rvec[dframe] += ( particle_rvec[idx] * invN ); // compute center of mass

          idx_old = N*0 + i;       // t0 for this sample trajectory
          dr = particle_rvec[idx] - particle_rvec[idx_old];  // displacement r(t)-r(t0)
          dr2 = dr.sq();
          if(verbose) r2t0 += dr2 * invN;          // sample trajectory
          t=ps[i].label;                  // type: integer 0,1,...,nTypes-1
          r2[t][dframe-1] += dr2;         // average over t0 and atoms of same type
          r4[t][dframe-1] += dr2*dr2;
        }
        if(verbose){
          ss.str(std::string()); ss << string_out_msd << tag << ".traj"; fout.open(ss.str(), ios::app);
          fout << r2t0 << endl; // output sample trajectory
          if(dframe == period_in_dt_units-1) fout << endl; // end of the block (sample trajectory)
          fout.close();
        }
        // after evaluating all particles, do the same for the center of mass
        dr = CM_rvec[dframe] - CM_rvec[0];
        dr2 = dr.sq();
        r2CM[dframe-1] += dr2;
      }
      if(debug) cout << "*** "<<myName<<" computation for timestep " << timestep << " ENDED ***\n\n";
    }

    void print(int dtframe)
    {
      ss.str(std::string()); ss << string_out_msd << tag << ".ave"; fout.open(ss.str(), ios::app);
      ntype r2_error[nTypes];

      for(int i=0; i<period_in_dt_units-1; i++) {
        // normalize the averages
        r2CM[i] /= num_periods;
        for(int t=0;t<nTypes;t++){
          r2[t][i] /= (num_periods*Nt[t]);
          r4[t][i] /= (num_periods*Nt[t]);
          r2_error[t] = sqrt( (r4[t][i] - r2[t][i]*r2[t][i])/(num_periods*Nt[t]) );
        }
        // print MSD to file
        fout << (i+1)*dtframe << " ";
        for(int t=0;t<nTypes;t++) fout << r2[t][i] << " ";
        for(int t=0;t<nTypes;t++) fout << r2_error[t] << " ";
        fout << r2CM[i] << endl;
      }
      fout.close();

      // print NGP to file
      ss.str(std::string()); ss << string_out_ngp << tag << ".ave"; fout.open(ss.str(), ios::app);
      for(int i=0; i<period_in_dt_units-1; i++) {
        fout << (i+1)*dtframe << " ";
        for(int t=0;t<nTypes;t++) fout << 0.6*(r4[t][i]/(r2[t][i]*r2[t][i]))-1.0 << endl;
      }
      fout.close();
      if(debug) cout << myName<<" printed to file\n";
    }
};

#endif
