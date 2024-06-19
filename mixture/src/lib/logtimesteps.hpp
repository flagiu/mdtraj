#ifndef _LOGTIMESTEPS_H_
#define _LOGTIMESTEPS_H_
using namespace std;

class LogTimesteps
{
public:
  int npc, ncycles, delta, t0;
  float alpha;
  int log_period, time_window_size;
  bool setup=false;

  int get_time_window_size(){
    time_window_size=npc-1+ncycles;
    return time_window_size;
  }

  int get_log_period(){
    log_period = t0 * round(pow(alpha,npc));
    return log_period;
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
    if(i<npc) return ncycles; // it is constant for dt=0 and for log sampling
    else if(i<time_window_size) {
      int icycle=i-npc+1;
      return ncycles-icycle; // it linearly decreases for linear sampling
    }
    else {
      cout << "ERROR: index "<<i<<" exceeds time_window_size "<<time_window_size<<endl;
      exit(1);
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
    cout<<"      log_period="<<log_period<<endl;
    cout<<"time_window_size="<<time_window_size<<endl;
    cout<<"#-----------------------------#\n\n";
  }

  void print_all_timesteps(){
    int i,j, t;
    t0=(int)(delta/alpha);
    for(i=0;i<ncycles;i++){
      for(j=0;j<npc;j++){
        t = t0 * (pow(alpha,npc)*i + pow(alpha,(j+1)) );
        cout<<t<<endl;
      }
    }
  }

  void deduce_fromfile(string file){
    ifstream i;
    string line;
    int t1,t2,dt,dt_old,count=0;
    float ratio,ratio_old;
    const int maxcount=99999999;
    i.open(file, ios::in);
    if (!i.is_open()){
      cout << "ERROR: Unable to open file "<<file<<endl;
      exit(1);
    }
    // deduce npc, delta, alpha, t0
    while ( count<maxcount && getline(i,line)){
      t2=stoi(line);
      if(count==0) delta=t2;
      if(count>0){
        dt=t2-t1;
        if(count>1){
          ratio=(float)t2/(float)t1;
          if(dt<dt_old){
            npc=count;
            count=maxcount+1;
            alpha=ratio_old;
            t0=(int)(((float)delta)/alpha);
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
      log_period=other.log_period;
      time_window_size=other.time_window_size;
      setup=other.setup;
      return (*this);
    }

};

#endif
