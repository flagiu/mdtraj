#ifndef _CORRELATE_ATOMIC_H_
#define _CORRELATE_ATOMIC_H_
using namespace std;

// Compute the time correlation of a vector of per-atom scalars

//------------------------------------- XXX ---------------------------------------//

template <class ntype>
class AtomicTimeCorrelator
{
  using vec=myvec<ntype,3>;
  using mat=mymatrix<ntype,3,3>;
  private:
    int period_in_dt_units, num_periods, ntimes, N, counts1;
    vector<int> num_avg, num_avg_predicted;
    ntype x1static, x2static, x4static;
    vector<ntype> x_history, x2, x4;
    string string_out, myName, tag;
    fstream fout;
    stringstream ss;
    bool debug,verbose, logtime;
    LogTimesteps* logt;

  public:
    AtomicTimeCorrelator(){
      myName = "AtomicTimeCorrelator";
    }
    virtual ~AtomicTimeCorrelator(){}

    void initLogTime()
    {
      ntimes = logt->get_time_window_size(); // (log+lin) = npc-1 + ncycles_not_skipped
      num_avg_predicted.resize(ntimes-1);
      num_avg.resize(ntimes-1);
      for(int j=0;j<ntimes-1;j++) {
        num_avg[j]=0;
        num_avg_predicted[j]=logt->get_num_avg(j+1);
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

      num_avg_predicted.resize(ntimes-1);
      num_avg.resize(ntimes-1);
      for(int j=0;j<ntimes-1;j++) {
        num_avg[j]=0;
        num_avg_predicted[j]=logt->get_num_avg(j+1);
      }

    }

    void init(int dtframe, int nframes, int period_in_real_units, int N_,
              bool logtime_, LogTimesteps* logt_, string string_out_,
              string tag_, bool debug_, bool verbose_)
    {
      N=N_;
      logtime = logtime_;
      logt = logt_;
      string_out = string_out_;
      tag = tag_;
      debug=debug_;
      verbose=verbose_;
      int i,t,tp;
      if(debug) cerr << myName << " Initialization STARTED\n";

      if(logtime){
        initLogTime();
      } else {
        initLinearTime(period_in_real_units, dtframe, nframes);
      }

      if(debug) cerr<<"Allocating "<<N*ntimes<<" elements for x_history...\n";
      x_history.resize(N*ntimes);
      // index = N*(current_time_index) + particle_index,
      //   with 0 <= current_time_index < ntimes,
      //   with 0 <= particle_index < N.
      if(debug) cerr<<"Mallocating "<<ntimes-1<<" elements for x2 and x4...\n";
      x1static=x2static=x4static=0;
      counts1=0;
      x2.resize(ntimes-1);
      x4.resize(ntimes-1);
      for(int i=0;i<ntimes-1;i++){
        x2[i]=x4[i]=0;
      }

      ss.str(std::string()); ss << string_out << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "#Delta timestep | <dx(t)*dx(0)>/<dx(0)*dx(0)>\n";
      fout.close();

      if(debug) cerr << myName << " Initialization COMPLETED\n";
    }

    void correlate(int frameidx, int timestep, vector<ntype> x)
    {
      int i,j,k, idx, idx_old, dframe, nperiod, num_images, idx_dt;
      ntype xx;
      if(debug) cerr << "*** "<<myName<<" computation for timestep " << timestep << " STARTED ***\n";

      if(logtime){
        nperiod = frameidx / logt->npc;
        dframe = frameidx % logt->npc; // log-spaced at short time
        if(dframe==0 && nperiod>0) dframe = logt->npc-1+nperiod; // linearly-spaced at long time
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
        for(i=0;i<N;i++)
        {
          x_history[N*dframe + i] = x[i]; //store initial value for this sample trajectory
        }
        if(debug) cerr << "  Saved values of initial frame\n";
      }
      else
      {
        idx_dt = dframe-1; // = 0,...,num_dt-1 however dt>0 always
        num_avg[idx_dt]+=1;
        for(i=0;i<N;i++)
        {
          idx = N*dframe + i;
          x_history[idx] = x[i];
          x1static += x[i];
          x2static += x[i]*x[i];
          x4static += x[i]*x[i]*x[i]*x[i];
          counts1 += 1;

          // Calculate x(t0+dt) * x(t0) w.r.t. REFERENCE frame t0 indexed by idx_old:
          // - linear sampling: t0 is the beginning of each cycle
          // - linlog sampling: t0 depends:
          if(logtime){
            if(dframe<logt->npc){
              if(nperiod==0) idx_old=N*0 + i; // if first period, t0 is the begining of the first npc log cycle
              else           idx_old=N*(logt->npc-1+nperiod)+i; // else t0 is the begining of the current npc log cycle
            }
            else             idx_old=N*0 + i; // for linear subsampling, first compare with first absolute snapshot
          } else idx_old=N*0 + i;

          xx = x_history[idx] * x_history[idx_old];
          x2[idx_dt] += xx;
          x4[idx_dt] += xx*xx;
        }

        if(logtime && dframe>=logt->npc) {
          if(debug) cerr << "  Out-of-log-cycle linear subsampling:\n";
          for(k=1;k<=nperiod-1;k++) { // subtract a decreasing time interval
            idx_dt = dframe-1-k;
            num_avg[idx_dt]+=1;
            for(i=0;i<N;i++){
              idx = N*dframe + i;
              idx_old=N*(logt->npc-1+k)+i; // ... t0 increases with k
              xx = x_history[idx] * x_history[idx_old];
              x2[idx_dt] += xx;
              x4[idx_dt] += xx*xx;
            }
          }
        }

      }

      if(debug) cerr << "*** "<<myName<<" computation for timestep " << timestep << " ENDED ***\n\n";
    }

    void print(int dtframe)
    {
      int i;
      ntype c;
      if(debug){
        cerr<<"Please check that num_avg is the same as the predicted one:\n";
        cerr<<" index num_avg predicted\n";
        for(i=0;i<ntimes-1;i++)
          cout << " "<<i<<" "<<num_avg[i]<<" "<<num_avg_predicted[i]<<endl;
      }

      ss.str(std::string()); ss << string_out << tag << ".ave"; fout.open(ss.str(), ios::app);

      x1static/=counts1;
      x2static/=counts1;
      x4static/=counts1;
      fout<<"#<x> = "<<x1static<<" <x^2> = "<<x1static<<" <x^4> = "<<x4static<<endl;
      for(i=0; i<ntimes-1; i++) {
        // normalize the averages
        x2[i] /= (num_avg[i]*N);
        x4[i] /= (num_avg[i]*N);
        c=(x2[i]-x1static*x1static)/(x2static-x1static*x1static);
        // print to file
        if(logtime) fout << logt->get_dt(i+1) << " ";
        else        fout << (i+1)*dtframe << " ";
        fout<<c<<endl;
      }
      fout.close();

      if(debug) cerr << myName<<" printed to file\n";
    }
};

#endif
