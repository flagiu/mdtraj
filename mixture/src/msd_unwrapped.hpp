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
    int period_in_dt_units, num_periods, ntimes, N, nTypes, nTypePairs;
    vector<int> Nt;
    vector<vec> particle_rvec; // collection of all 3D coordinates at all times
    vector<vec> CM_rvec; // collection of Center-of-Mass 3D coordinates at all times
    vecflex< vecflex<ntype> > r2; // < |r(t) - r(t0)|^2 > averaged over t0 and particles, for each type
    vecflex< vecflex<ntype> > r4; // < |r(t) - r(t0)|^4 > averaged over t0 and particles, for each type
    vecflex<ntype> r2CM; // < |rCM(t) - rCM(t0)|^2 > averaged over t0, for all types
    vecflex<ntype> num_avg, num_avg_predicted;
    vecflex< vecflex<ntype> > Qs,Qss, Qd,Qdd, Qsd; // self- and distinct- overlap parameter
    vecflex< ntype > Qs_sample, Qd_sample;
    ntype Q_cutoff, Q_cutoff2; // cutoff for the overlap parameter
    string string_out_msd, string_out_ngp, string_out_overlap, myName, tag;
    fstream fout;
    stringstream ss;
    bool debug,verbose, logtime;
    LogTimesteps logt;

    ntype w(ntype dr2) {
      return dr2>Q_cutoff2 ? 0. : 1.;
    }

    int types2int(int ti, int tj){ // map type pairs (ti,tj) in 0,1,...,nTypes-1 to symmetric integer index 0,1,...,nTypePairs
      if (ti>tj) return types2int(tj,ti); // map to ti<=tj
      if(ti==tj) return ti*nTypes - int(ti*(ti-1)/2);
      else return 1 + types2int(ti,tj-1);
    }

  public:
    MSDU_Calculator(){
      myName = "MSDU and NGP and Self-Overlap";
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
        cerr << "WARNING: the input period for MSD calculation exceeds maximum time interval.\n";
        cerr << "         I will set it equal to the maximum time interval.\n";
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

    void init(int dtframe, int nframes, int period_in_real_units, ntype Q_cutoff_, int N_,
      int nTypes_, int *Nt_, bool logtime_, LogTimesteps logt_,
      string string_out_msd_, string string_out_ngp_, string string_out_overlap_,
      string tag_, bool debug_, bool verbose_)
    {
      if(debug) cout << myName << " Initialization STARTED\n";
      Q_cutoff = Q_cutoff_;
      Q_cutoff2 = Q_cutoff*Q_cutoff;
      N=N_;
      nTypes=nTypes_;
      nTypePairs = nTypes*(nTypes+1)/2;
      Nt.resize(nTypes);
      for(int t=0;t<nTypes;t++) Nt[t]=Nt_[t]; // number of atoms for each type
      logtime = logtime_;
      logt = logt_;
      string_out_msd = string_out_msd_;
      string_out_ngp = string_out_ngp_;
      string_out_overlap = string_out_overlap_;
      tag = tag_;
      debug=debug_;
      verbose=verbose_;
      int i,t,tp;

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
      Qs.resize(nTypes);
      Qss.resize(nTypes);
      Qd.resize(nTypePairs);
      Qdd.resize(nTypePairs);
      Qsd.resize(nTypePairs);
      Qs_sample.resize(nTypes);
      Qd_sample.resize(nTypePairs);
      for(t=0;t<nTypes;t++){
        r2[t].resize(ntimes-1); // < |r(t0+dt) - r(t0)|^2 > , dt>0, for type t
        r4[t].resize(ntimes-1); // < |r(t0+dt) - r(t0)|^4 > , dt>0, for type t
        Qs[t].resize(ntimes-1);
        Qss[t].resize(ntimes-1);
      }
      for(tp=0;tp<nTypePairs;tp++){
        Qd[tp].resize(ntimes-1);
        Qdd[tp].resize(ntimes-1);
        Qsd[tp].resize(ntimes-1);
      }
      r2CM.resize(ntimes-1); // < |rCM(t0+dt) - rCM(t0)|^2 > , dt>0

      for(i=0;i<ntimes-1;i++){ // reset averages
        for(t=0;t<nTypes;t++){
          r2[t][i]=0.0;
          r4[t][i]=0.0;
          Qs[t][i]=0.0;
          Qss[t][i]=0.0;
        }
        for(tp=0;tp<nTypePairs;tp++){
          Qd[tp][i]=0.0;
          Qdd[tp][i]=0.0;
          Qsd[tp][i]=0.0;
        }
        r2CM[i]=0.0;
      }
      for(i=0;i<ntimes-1;i++) CM_rvec[i]<<0.,0.,0.;

      ss.str(std::string()); ss << string_out_msd << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "#Delta timestep | MSD for each type | MSD error for each type | MSD of c.o.m.\n";
      fout.close();
      ss.str(std::string()); ss << string_out_ngp << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "#Delta timestep | NGP for each type (-0.4: homog. displ., 0.0: Brownian displ., >0: heterog. displ.)\n";
      fout.close();
      ss.str(std::string()); ss << string_out_overlap << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "#Delta timestep | averaged qs(1) qss(1) qd(1,1) qdd(1,1) qsd(1,1) qd(1,2) qdd(1,2) qsd(1,2) qs(2) qss(2) qd(2,2) ...\n";
      fout.close();

      if(debug) cout << myName << " Initialization COMPLETED\n";
    }

    void compute(int frameidx, int timestep, vector<ptype> ps, PBC<ntype> *pbc)
    {
      int i,j,k, idx, idx_old, dframe, nperiod, num_images, ti,tj,tp, idx_dt;
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
        idx_dt = dframe-1; // = 0,...,num_dt-1 however dt>0 always
        num_avg[idx_dt]+=1;
        CM_rvec[dframe]<<0.,0.,0.;
        if(debug) cout << "  Computing dr^2 and Qoverlap for dframe="<<dframe<<endl;
        if(verbose) r2t0 = 0.0; // average r^2 over atoms at fixed t,t0
        for(ti=0;ti<nTypes;ti++){ Qs_sample[ti]=0.0; }
        for(tp=0;tp<nTypePairs;tp++){ Qd_sample[tp]=0.0; }
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
          // apply periodic boundary conditions
          num_images = pbc->apply(dr, &dr_image);
          dr = dr_image;
          particle_rvec[idx] = particle_rvec[idx_old] + dr;
          // MSD and NGP
          dr2 = dr.sq();
          if(verbose) r2t0 += dr2 * invN;          // sample trajectory
          ti=ps[i].label;                  // type: integer 0,1,...,nTypes-1
          r2[ti][idx_dt] += dr2;         // average over t0 and atoms of same type
          r4[ti][idx_dt] += dr2*dr2;
          // overlap parameter, normalized to N
          // -- self
          Qs_sample[ti] += w(dr2);
          // -- distinct
          for(j=0;j<i;j++){
            dr = particle_rvec[idx] - particle_rvec[idx_old-i+j];  // displacement ri(t)-rj(t0)
            num_images = pbc->apply(dr, &dr_image); // apply periodic boundary conditions
            tj = ps[j].label;
            Qd_sample[types2int(ti,tj)] += 2*w(dr.sq()); // i<j + j>i
          }
        }
        // normalization to N of 2-body terms, and running average of 2-body and 4-body terms
        for(ti=0;ti<nTypes;ti++){
          Qs_sample[ti] *= invN;
          Qs[ti][idx_dt] += Qs_sample[ti];
          Qss[ti][idx_dt] += Qs_sample[ti]*Qs_sample[ti];
          /* if(idx_dt==0) {
            cerr << "#---- Adding Qs*Qs to Qss ----#\n";
            cerr << "# idx_dt="<<idx_dt<<" i_type="<<ti<<" navg="<<num_avg[idx_dt]<<" Qs_sample="<<Qs_sample[ti]<<endl;
            cerr << "# Qs="<<Qs[ti][idx_dt]<<" Qss="<<Qss[ti][idx_dt]<<" (to be normalized to navg)\n";
          } */
        }
        // (first normalize every self term, then go to distinct terms)
        for(ti=0;ti<nTypes;ti++){
          for(tj=ti;tj<nTypes;tj++){
            tp=types2int(ti,tj);
            Qd_sample[tp] *= invN;
            Qd[tp][idx_dt] += Qd_sample[tp];
            Qdd[tp][idx_dt] += Qd_sample[tp]*Qd_sample[tp];
            Qsd[tp][idx_dt] += Qs_sample[ti]*Qd_sample[tp];
            if(tj!=ti){ Qsd[tp][idx_dt] += Qs_sample[tj]*Qd_sample[tp]; } // type exchange
          }
        }

        if(logtime && dframe>=logt.npc) {
          if(debug) cout << "  Out-of-log-cycle linear subsampling:\n";
          for(k=1;k<=nperiod-1;k++) { // subtract a decreasing time interval
            idx_dt = dframe-1-k;
            num_avg[idx_dt]+=1;
            for(ti=0;ti<nTypes;ti++){ Qs_sample[ti]=0.0; }
            for(tp=0;tp<nTypePairs;tp++){ Qd_sample[tp]=0.0; }
            for(i=0;i<N;i++){
              idx = N*dframe + i;
              idx_old=N*(logt.npc-1+k)+i; // ... t0 increases with k
              dr = particle_rvec[idx] - particle_rvec[idx_old];  // displacement r(t)-r(t0)
              num_images = pbc->apply(dr, &dr_image); // apply periodic boundary conditions
              dr2 = dr.sq();
              ti=ps[i].label;                  // type: integer 0,1,...,nTypes-1
              r2[ti][idx_dt] += dr2;         // average over t0 and atoms of same type
              r4[ti][idx_dt] += dr2*dr2; // -time interval decreases with k
              // overlap parameter, normalized to N
              // -- self
              Qs_sample[ti] += w(dr2);
              // -- distinct
              for(j=0;j<i;j++){
                dr = particle_rvec[idx] - particle_rvec[idx_old-i+j];  // displacement ri(t)-rj(t0)
                num_images = pbc->apply(dr, &dr_image); // apply periodic boundary conditions
                tj = ps[j].label;
                Qd_sample[types2int(ti,tj)] += 2*w(dr.sq()); // i<j + j>i
              }
            }
            // normalization to N of 2-body terms, and running average of 2-body and 4-body terms
            for(ti=0;ti<nTypes;ti++){
              Qs_sample[ti] *= invN;
              Qs[ti][idx_dt] += Qs_sample[ti];
              Qss[ti][idx_dt] += Qs_sample[ti]*Qs_sample[ti];
              /* if(idx_dt==0) {
                cerr << "#---- Adding Qs*Qs to Qss ----#\n";
                cerr << "# idx_dt="<<idx_dt<<" i_type="<<ti<<" navg="<<num_avg[idx_dt]<<" Qs_sample="<<Qs_sample[ti]<<endl;
                cerr << "# Qs="<<Qs[ti][idx_dt]<<" Qss="<<Qss[ti][idx_dt]<<" (to be normalized to navg)\n";
              } */
            }
            // (first normalize every self term, then go to distinct terms)
            for(ti=0;ti<nTypes;ti++){
              for(tj=ti;tj<nTypes;tj++){
                tp=types2int(ti,tj);
                Qd_sample[tp] *= invN;
                Qd[tp][idx_dt] += Qd_sample[tp];
                Qdd[tp][idx_dt] += Qd_sample[tp]*Qd_sample[tp];
                Qsd[tp][idx_dt] += Qs_sample[ti]*Qd_sample[tp];
                if(tj!=ti){ Qsd[tp][idx_dt] += Qs_sample[tj]*Qd_sample[tp]; } // type inversion
              }
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
          for(k=1;k<=nperiod-1;k++) { // subtract a decreasing time interval
            idx_dt = dframe-1-k;
            idx_old=logt.npc-1+k; // ... t0 increases with j
            dr = CM_rvec[dframe] - CM_rvec[idx_old];  // displacement r(t)-r(t0)
            num_images = pbc->apply(dr, &dr_image); // apply periodic boundary conditions
            dr2 = dr.sq();
            r2CM[idx_dt] += dr2;
          }
        }

      }

      if(debug) cout << "*** "<<myName<<" computation for timestep " << timestep << " ENDED ***\n\n";
    }

    void print(int dtframe)
    {
      int i,t, ti,tj,tp;
      if(debug){
        cout<<"Please check that num_avg is the same as the predicted one:\n";
        cout<<" index num_avg predicted\n";
        for(i=0;i<ntimes-1;i++)
          cout << " "<<i<<" "<<num_avg[i]<<" "<<num_avg_predicted[i]<<endl;
      }

      ss.str(std::string()); ss << string_out_msd << tag << ".ave"; fout.open(ss.str(), ios::app);
      ntype tmpvar[nTypes];

      for(i=0; i<ntimes-1; i++) {
        // normalize the averages
        r2CM[i] /= num_avg[i];
        for(int t=0;t<nTypes;t++){
          r2[t][i] /= (num_avg[i]*Nt[t]);
          r4[t][i] /= (num_avg[i]*Nt[t]);
          tmpvar[t] = sqrt( (r4[t][i] - r2[t][i]*r2[t][i])/(num_avg[i]*Nt[t]-1) );
        }
        // print MSD to file
        if(logtime) fout << logt.get_dt(i+1) << " ";
        else        fout << (i+1)*dtframe << " ";
        for(int t=0;t<nTypes;t++) fout << r2[t][i] << " ";
        for(int t=0;t<nTypes;t++) fout << tmpvar[t] << " ";
        fout << r2CM[i] << endl;
      }
      fout.close();

      // print NGP to file
      ss.str(std::string()); ss << string_out_ngp << tag << ".ave"; fout.open(ss.str(), ios::app);
      for(i=0; i<ntimes-1; i++) {
        if(logtime) fout << logt.get_dt(i+1) << " ";
        else        fout << (i+1)*dtframe << " ";
        for(int t=0;t<nTypes;t++) fout << 0.6*(r4[t][i]/(r2[t][i]*r2[t][i]))-1.0 << endl;
      }
      fout.close();

      ss.str(std::string()); ss << string_out_overlap << tag << ".ave"; fout.open(ss.str(), ios::app);
      for(i=0; i<ntimes-1; i++) {
        // normalize the averages and print to file
        if(logtime) fout << logt.get_dt(i+1) << " ";
        else        fout << (i+1)*dtframe << " ";
        for(ti=0;ti<nTypes;ti++){
          Qs[ti][i] /= num_avg[i];
          Qss[ti][i] /= num_avg[i];
          fout << Qs[ti][i] << " " << Qss[ti][i] << " ";
          for(tj=ti;tj<nTypes;tj++){
            tp=types2int(ti,tj);
            Qd[tp][i] /= num_avg[i];
            Qdd[tp][i] /= num_avg[i];
            Qsd[tp][i] *= 2./num_avg[i]; // 2x because of double product Qs*Qd + Qd*Qs
            fout << Qd[tp][i] << " " << Qdd[tp][i] << " " << Qsd[tp][i] << " ";
          }
        }
        fout << endl;
      }
      fout.close();

      if(debug) cout << myName<<" printed to file\n";
    }
};

#endif
