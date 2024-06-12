#ifndef _MSDU_H_
#define _MSDU_H_
using namespace std;

//---- "init" functions: prepare arrays, averages and headers for output files
//---- "compute" functions: compute dynamic averages
//---- "print" functions: complete the calculations at the end of trajectory and print out

//------------------------------------- Mean Squared Displacement and Non-Gaussian Parameter ---------------------------------------//

template <class ntype, class ptype>
class MSDU_Calculator
{
  using vec=myvec<ntype,3>;
  using mat=mymatrix<ntype,3,3>;
  private:
    int period_in_dt_units, num_periods, ntimes, N, nTypes;
    vector<int> Nt;
    vector<vec> particle_rvec; // collection of all 3D coordinates at all times
    vector<vec> CM_rvec; // collection of Center-of-Mass 3D coordinates at all times
    vecflex< vecflex<ntype> > r2; // < |r(t) - r(t0)|^2 > averaged over t0 and particles, for each type
    vecflex< vecflex<ntype> > r4; // < |r(t) - r(t0)|^4 > averaged over t0 and particles, for each type
    vecflex<ntype> r2CM; // < |rCM(t) - rCM(t0)|^2 > averaged over t0, for all types
    vecflex<ntype> num_avg, num_avg_predicted;
    string string_out_msd, string_out_ngp, myName, tag;
    fstream fout;
    stringstream ss;
    bool debug,verbose, logtime;
    LogTimesteps logt;
  public:
    MSDU_Calculator(){
      myName = "MSDU and NGP";
    }
    virtual ~MSDU_Calculator(){}

    void initLogTime()
    {
      ntimes = logt.get_time_window_size();

      num_avg.resize(ntimes-1);
      num_avg_predicted.resize(ntimes-1);
      for(int j=0;j<ntimes-1;j++) {
        num_avg[j]=0.0;
        num_avg_predicted[j]=logt.get_num_avg(j+1);
      }

      if(verbose)
      {
        ss.str(std::string()); ss << string_out_msd << tag << ".traj"; fout.open(ss.str(), ios::out);
        fout << "#Delta timestep, then |r(t)-r(t0)|^2 # block file for different t0\n";
        for(int i=0;i<ntimes-1;i++){ // reset averages
          for(int t=0;t<nTypes;t++){
            fout << logt.get_dt(i+1) << endl;
          }
        }
        fout << endl; // end of the first block (delta timesteps)
        fout.close();
      }

    }

    void initLinearTime(int period_in_real_units, int dtframe, int nframes)
    {
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
      if(debug) cout << myName<<": Set period_in_dt_units = " << period_in_dt_units
        << ", num_periods = " << num_periods << endl;

      ntimes=period_in_dt_units;

      num_avg.resize(ntimes-1);
      num_avg_predicted.resize(ntimes-1);
      for(int j=0;j<ntimes-1;j++) {
        num_avg[j]=0.0;
        num_avg_predicted[j]=num_periods;
      }

      if(verbose)
      {
        ss.str(std::string()); ss << string_out_msd << tag << ".traj"; fout.open(ss.str(), ios::out);
        fout << "#Delta timestep, then |r(t)-r(t0)|^2 # block file for different t0\n";
        for(int i=0;i<ntimes-1;i++){ // reset averages
          for(int t=0;t<nTypes;t++){
            fout << (i+1)*dtframe << endl;
          }
        }
        fout << endl; // end of the first block (delta timesteps)
        fout.close();
      }

    }

    void init(int dtframe, int nframes, int period_in_real_units, int N_,
      int nTypes_, int *Nt_, bool logtime_, LogTimesteps logt_,
      string string_out_msd_, string string_out_ngp_,
      string tag_, bool debug_, bool verbose_)
    {
      if(debug) cout << myName << " Initialization STARTED\n";
      N=N_;
      nTypes=nTypes_;
      Nt.resize(nTypes);
      for(int t=0;t<nTypes;t++) Nt[t]=Nt_[t]; // number of atoms for each type
      logtime = logtime_;
      logt = logt_;
      string_out_msd = string_out_msd_;
      string_out_ngp = string_out_ngp_;
      tag = tag_;
      debug=debug_;
      verbose=verbose_;

      if(logtime){
        initLogTime();
      } else {
        initLinearTime(period_in_real_units, dtframe, nframes);
      }

      particle_rvec.resize(N*ntimes);
      // index = N*(current_time_index) + particle_index,
      //   with 0 <= current_time_index < ntimes,
      //   with 0 <= particle_index < N.
      CM_rvec.resize(ntimes);
      r2.resize(nTypes);
      r4.resize(nTypes);
      for(int t=0;t<nTypes;t++){
        r2[t].resize(ntimes-1); // < |r(t0+dt) - r(t0)|^2 > , dt>0, for type t
        r4[t].resize(ntimes-1); // < |r(t0+dt) - r(t0)|^4 > , dt>0, for type t
      }
      r2CM.resize(ntimes-1); // < |rCM(t0+dt) - rCM(t0)|^2 > , dt>0

      for(int i=0;i<ntimes-1;i++){ // reset averages
        for(int t=0;t<nTypes;t++){
          r2[t][i]=0.0;
          r4[t][i]=0.0;
        }
        r2CM[i]=0.0;
      }
      for(int i=0;i<ntimes-1;i++) CM_rvec[i]<<0.,0.,0.;

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
      int i,j, idx, idx_old, dframe, nperiod, num_images, t;
      const ntype invN=1.0/(ntype)N;
      vec dr, dr_image;
      ntype dr2, r2t0;
      if(debug) cout << "*** "<<myName<<" computation for timestep " << timestep << " STARTED ***\n";

      if(logtime){
        nperiod = frameidx / logt.npc;
        dframe = frameidx % logt.npc; // log-spaced at short time
        if(dframe==0 && nperiod>0) dframe = logt.npc-1+nperiod; // linearly-spaced at long time
        /*
        frameidx 0 1 ... npc-1 npc npc+1 ... 2*npc-1 2*npc 2*npc+1 ...
        nperiod  0 0 ...     0  1    1   ...    1       2     2    ...
        dframe   0 1 ... npc-1 npc!! 1   ...  npc-1  npc+1!!  1    ...

        dframe   0 1 ... npc-1 npc npc+1 npc+2 ... npc-1+ncycles
        frameidx 0 1 ... npc-1 npc 2*npc 3*npc ... (ncycles-1+1)*npc
        nperiod  0 0 ...   0    1    2     3   ...     ncycles
        */
      } else {
        nperiod = frameidx / ntimes;
        dframe = frameidx % ntimes;
      }

      if(dframe==0) // either linear and every nperiod, or logtime and first absolute frame
      {
        CM_rvec[dframe]<<0.,0.,0.;
        for(i=0;i<N;i++)
        {
          particle_rvec[N*dframe + i] = ps[i].ru; //store initial positions for this sample trajectory
          CM_rvec[dframe] += ( particle_rvec[N*dframe + i] * invN ); // store initial center of mass for '''
        }
        if(debug) cout << "  Saved positions of initial frame\n";
      }
      else
      {
        num_avg[dframe-1]+=1;
        CM_rvec[dframe]<<0.,0.,0.;
        if(debug) cout << "  Computing dr^2 for dframe="<<dframe<<endl;
        if(verbose) r2t0 = 0.0; // average r^2 over atoms at fixed t,t0
        for(i=0;i<N;i++)
        {
          idx = N*dframe + i;
          //idx_old = idx - N;       // previous frame
          particle_rvec[idx] = ps[i].ru;
          CM_rvec[dframe] += ( particle_rvec[idx] * invN ); // compute center of mass

          // Calculate r(t0+dt) - r(t0):
          // - linear sampling: t0 is the beginning of each cycle
          // - linlog sampling: t0 depends:
          if(logtime){
            if(dframe<logt.npc){
              if(nperiod==0) idx_old=N*0 + i; // else if first period, t0 is the begining of the first npc log cycle
              else           idx_old=N*(logt.npc-1+nperiod)+i; // else t0 is the begining of the current npc log cycle
            }
            else             idx_old=N*0 + i; // for linear subsampling, first compare with first absolute snapshot
          } else idx_old=N*0 + i;

          dr = particle_rvec[idx] - particle_rvec[idx_old];  // displacement r(t)-r(t0)
          num_images = pbc->apply(dr, &dr_image); // apply periodic boundary conditions
          dr2 = dr.sq();
          if(verbose) r2t0 += dr2 * invN;          // sample trajectory
          t=ps[i].label;                  // type: integer 0,1,...,nTypes-1
          r2[t][dframe-1] += dr2;         // average over t0 and atoms of same type
          r4[t][dframe-1] += dr2*dr2;

          if(logtime && dframe>=logt.npc) {
            if(debug) cout << "  Out-of-log-cycle linear subsampling:\n";
            for(j=1;j<=nperiod-1;j++) { // subtract a decreasing time interval
              idx_old=N*(logt.npc-1+j)+i; // ... t0 increases with j
              dr = particle_rvec[idx] - particle_rvec[idx_old];  // displacement r(t)-r(t0)
              num_images = pbc->apply(dr, &dr_image); // apply periodic boundary conditions
              dr2 = dr.sq();
              t=ps[i].label;                  // type: integer 0,1,...,nTypes-1
              r2[t][dframe-1-j] += dr2;         // average over t0 and atoms of same type
              r4[t][dframe-1-j] += dr2*dr2; // -time interval decreases with j
              //num_avg[dframe-1-j]+=1*invN; // time interval decreases with j (DONE in CM for saving operations)
            }
          }

        }

        if(verbose){
          ss.str(std::string()); ss << string_out_msd << tag << ".traj"; fout.open(ss.str(), ios::app);
          fout << r2t0 << endl; // output sample trajectory
          if(logtime){
            if(dframe == logt.npc-1) fout << endl; // end of the log block (sample trajectory)
          }else{
            if(dframe == ntimes-1) fout << endl; // end of the block (sample trajectory)
          }
          fout.close();
        }

        // after evaluating all particles, do the same for the center of mass
        if(logtime){
          if(dframe<logt.npc){
            if(nperiod==0) idx_old=0; // else if first period, t0 is the begining of the first npc log cycle
            else           idx_old=logt.npc-1+nperiod; // else t0 is the begining of the current npc log cycle
          }
          else             idx_old=0; // for linear subsampling, first compare with first absolute snapshot
        } else idx_old=0;

        dr = CM_rvec[dframe] - CM_rvec[idx_old];
        num_images = pbc->apply(dr, &dr_image); // apply periodic boundary conditions
        dr2 = dr.sq();
        r2CM[dframe-1] += dr2;

        if(logtime && dframe>=logt.npc) {
          if(debug) cout << "  Out-of-log-cycle linear subsampling for CM:\n";
          for(j=1;j<=nperiod-1;j++) { // subtract a decreasing time interval
            idx_old=logt.npc-1+j; // ... t0 increases with j
            dr = CM_rvec[dframe] - CM_rvec[idx_old];  // displacement r(t)-r(t0)
            num_images = pbc->apply(dr, &dr_image); // apply periodic boundary conditions
            dr2 = dr.sq();
            r2CM[dframe-1-j] += dr2;
            num_avg[dframe-1-j]+=1;
          }
        }

      }

      if(debug) cout << "*** "<<myName<<" computation for timestep " << timestep << " ENDED ***\n\n";
    }

    void print(int dtframe)
    {
      if(debug){
        cout<<"Please check that num_avg is the same as the predicted one:\n";
        cout<<" index num_avg predicted\n";
        for(int j=0;j<ntimes-1;j++)
          cout << " "<<j<<" "<<num_avg[j]<<" "<<num_avg_predicted[j]<<endl;
      }

      ss.str(std::string()); ss << string_out_msd << tag << ".ave"; fout.open(ss.str(), ios::app);
      ntype r2_error[nTypes];

      for(int i=0; i<ntimes-1; i++) {
        // normalize the averages
        r2CM[i] /= num_avg[i];
        for(int t=0;t<nTypes;t++){
          r2[t][i] /= (num_avg[i]*Nt[t]);
          r4[t][i] /= (num_avg[i]*Nt[t]);
          r2_error[t] = sqrt( (r4[t][i] - r2[t][i]*r2[t][i])/(num_avg[i]*Nt[t]-1) );
        }
        // print MSD to file
        if(logtime) fout << logt.get_dt(i+1) << " ";
        else        fout << (i+1)*dtframe << " ";
        for(int t=0;t<nTypes;t++) fout << r2[t][i] << " ";
        for(int t=0;t<nTypes;t++) fout << r2_error[t] << " ";
        fout << r2CM[i] << endl;
      }
      fout.close();

      // print NGP to file
      ss.str(std::string()); ss << string_out_ngp << tag << ".ave"; fout.open(ss.str(), ios::app);
      for(int i=0; i<ntimes-1; i++) {
        if(logtime) fout << logt.get_dt(i+1) << " ";
        else        fout << (i+1)*dtframe << " ";
        for(int t=0;t<nTypes;t++) fout << 0.6*(r4[t][i]/(r2[t][i]*r2[t][i]))-1.0 << endl;
      }
      fout.close();

      if(debug) cout << myName<<" printed to file\n";
    }
};

#endif
