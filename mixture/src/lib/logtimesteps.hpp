#ifndef _LOGTIMESTEPS_H_
#define _LOGTIMESTEPS_H_
using namespace std;

class LogTimesteps
{
public:
  int npc=0, ncycles=0, delta=0, t0=0;
  int ncyc_skip0=0, ncyc_skip1=0;
  float alpha=0.0;
  int log_period, time_window_size, num_frames;
  bool setup=false, debug=false;

  int get_time_window_size(){
    time_window_size=npc-1+(ncycles-ncyc_skip0-ncyc_skip1);
    return time_window_size;
  }

  int get_num_frames(){
    num_frames = (ncycles-ncyc_skip0-ncyc_skip1)*npc;
    return num_frames;
  }

  int get_log_period(){
    log_period = t0 * round(pow(alpha,npc));
    return log_period;
  }

  int apply_fskip(float fskip0, float fskip1){
    if(!setup) {
      cerr << "ERROR: apply_skip() called before setup\n";
      exit(1);
    }
    // idea: you can only skip a whole cycle
    int num_frames_full = get_num_frames();
    int nskip0=int(fskip0*num_frames_full);
    int nskip1=int(fskip1*num_frames_full);
    ncyc_skip0=int(nskip0/(float)npc+0.5);
    ncyc_skip1=int(nskip1/(float)npc+0.5);
    if(debug){
      cerr << "# LogTimesteps\n";
      cerr << "# num_frames_full="<<num_frames_full<<endl;
      cerr << "# fskip0="<<fskip0<<" ==> nskip0="<<nskip0<<" ==> ncyc_skip0="<<ncyc_skip0<<"\n";
      cerr << "# fskip1="<<fskip1<<" ==> nskip1="<<nskip1<<" ==> ncyc_skip1="<<ncyc_skip1<<"\n";
    }
    num_frames = get_num_frames();
    if(num_frames<=0) {
      cerr << "ERROR: skipped too many log cycles:";
      cerr <<"\n      ncycles="<<ncycles<<" ncyc_skip0="<<ncyc_skip0<<" ncyc_skip1="<<ncyc_skip1<<endl;
      exit(1);
    }
    return num_frames;
  }

  bool is_linear_sampling(int timestep){
    return ((timestep % log_period)==delta);
  }

  int get_dt(int i){
    if(i==0) return 0; // dt=0 static properties
    else if(i<npc) return (int)(t0*pow(alpha,i+1)+0.5)-delta; // log sampling
    else if(i<time_window_size) { // linear sampling
      int icycle=i-npc+1;
      return (int)(t0*pow(alpha,npc)*icycle+0.5);
    }
    else {
      cout << "ERROR: index "<<i<<" exceeds time_window_size "<<time_window_size<<endl;
      exit(1);
    }
  }

  int get_num_avg(int i){
    if(i<npc) return ncycles-ncyc_skip0-ncyc_skip1; // it is constant for dt=0 and for log sampling
    else if(i<time_window_size) {
      int icycle=i-npc+1;
      return ncycles-ncyc_skip0-ncyc_skip1-icycle; // it linearly decreases for linear sampling
    }
    else {
      cout << "ERROR: index "<<i<<" exceeds time_window_size "<<time_window_size<<endl;
      exit(1);
    }
  }

  int get_timestep(int icycle, int i_incycle) {
    // round, in case alpha is not perfectly deduced
    return (int) (0.5+ t0 * (pow(alpha,npc)*icycle + pow(alpha,(i_incycle+1)) ) );
  }

  int get_first_timestep_noskip(){
    return get_timestep(0,0);
  }
  int get_last_timestep_noskip(){
    return get_timestep(ncycles-1,npc-1);
  }
  int get_first_timestep_withskip(){
    return get_timestep(ncyc_skip0,0);
  }
  int get_last_timestep_withskip(){
    return get_timestep(ncycles-ncyc_skip1-1,npc-1);
  }

  void print_all_timesteps(){
    int i,j, t;
    t0=(int)(delta/alpha);
    for(i=ncyc_skip0;i<ncycles-ncyc_skip1;i++){
      for(j=0;j<npc;j++){
        t = get_timestep(i,j);
        cout<<t<<endl;
      }
    }
  }

  void print_summary(){
    cout<<"#-----------------------------#\n";
    cout<<"LogTimesteps summary:\n";
    cout<<"           alpha="<<alpha<<endl;
    cout<<"           delta="<<delta<<endl;
    cout<<"             npc="<<npc<<endl;
    cout<<"         ncycles="<<ncycles<<endl;
    cout<<"              t0="<<t0<<endl;
    cout<<endl;
    cout<<"      ncyc_skip0="<<ncyc_skip0<<endl;
    cout<<"      ncyc_skip1="<<ncyc_skip1<<endl;
    cout<<"      log_period="<<log_period<<endl;
    cout<<"time_window_size="<<time_window_size<<endl;
    cout<<"#-----------------------------#\n\n";
  }

  void deduce_fromfile(string file, bool debug_){
    debug=debug_;
    ifstream i;
    string line;
    int t1,t2,dt,dt_old,count=0;
    float ratio,ratio_old;
    bool firstCycleEnded=false;
    const int maxcount=99999999;
    i.open(file, ios::in);
    if (!i.is_open()){
      cout << "ERROR: Unable to open file "<<file<<endl;
      exit(1);
    }
    // deduce npc, delta, alpha, t0 from a FULL schedule file
    while ( count<maxcount && getline(i,line)){
      t2=stoi(line);
      if(count==0) delta=t2;
      if(count>0){
        dt=t2-t1;
        if(count>1){
          ratio=(float)t2/(float)t1;
          // deduce alpha from the ratio btw the two highest numbers available
          // in the first cycle, i.e. the last two. Find them by checking then
          // the step difference is shorter than the previous one
          if(!firstCycleEnded && dt<dt_old){
            npc=count;
            count=maxcount+1;
            alpha=ratio_old;
            t0=(int)round(((float)delta)/alpha);
            firstCycleEnded=true;
          }
          ratio_old=ratio;
        }
        dt_old=dt;
      }
      t1=t2;
      count++;
    }
    if(count==maxcount){
      cout << "ERROR: too many lines in "<<file<<": "<<maxcount<<endl;
      exit(1);
    }
    // deduce n. of cycles
    count=npc+1;
    while ( count<maxcount && getline(i,line)){
      count++;
    }
    if(count<maxcount){
      ncycles=(int)(count/(float)npc); // there is some tolerance
    } else {
      cout << "ERROR: too many lines in "<<file<<": "<<maxcount<<endl;
      exit(1);
    }
    i.close();
    get_log_period();
    get_time_window_size();
    setup=true;
  }

  LogTimesteps& operator=(const LogTimesteps& other)
    {
      if(!other.setup) return (*this);
      npc=other.npc;
      ncycles=other.ncycles;
      delta=other.delta;
      t0=other.t0;
      alpha=other.alpha;
      ncyc_skip0=other.ncyc_skip0;
      ncyc_skip1=other.ncyc_skip1;
      log_period=other.log_period;
      time_window_size=other.time_window_size;
      setup=other.setup;
      return (*this);
    }

};

#endif
