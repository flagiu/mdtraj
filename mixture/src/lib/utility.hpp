#ifndef _UTILITY_H_
#define _UTILITY_H_

#include<cstdlib>
#include<fstream>
#include<iostream>
#include<iomanip>
#include<vector>
using namespace std;

#define SQUARE(x) ((x)*(x))

// Number of lines in a file
static int getLineCount(string file)
{
  int i=0;
  string str;
  fstream f;
  f.open(file, ios::in);
  while( getline(f,str) )
    i++;
  f.close();
  return i;
}

// Find index of a given element (first one, if it is repeated)
template<class T>
int indexOf(vector<T>* vec, T element)
{
  for(int i=0, size=vec->size(); i<size; i++)
    if((*vec)[i]==element) return i;
  return -1; // not found
}

class PrintProgress
{
private:
  int perc; // 'perc' is in percentage units
  int total_steps, print_every_ms, current_block;
public:
  void init(int tot_steps, int print_every_ms_)
  {
    total_steps = tot_steps;
    print_every_ms = print_every_ms_;
    perc=0;
    current_block=0;
  }
  void update(int step, float run_time_ms)
  {
    int block = ((int)run_time_ms) / print_every_ms;
    if(block>current_block)
    {
      current_block = block;
      perc = floor( step/(float)total_steps * 100.0);
      cout << "\r[" << setfill(' ') << setw(3) << perc << "%]" << flush;
    }
  }
  void end()
  {
    cout << "\r[100%]\n";
  }
};

int intsign(int x)
{
  if (x > 0) return 1;
  if (x < 0) return -1;
  return 0;
}
int floatsign(float x)
{
  if (x > 0) return 1;
  if (x < 0) return -1;
  return 0;
}

class Timer
{
private:
  clock_t start, end;
public:
  void go()
  {
    start = clock();
  }
  float lap()
  {
    end = clock();
    return (float) (end-start) / CLOCKS_PER_SEC * 1000.0; // elapsed time (ms)
  }
};

class LogTimesteps
{
public:
  int npc, ncycles, delta, t0;
  float alpha;
  int time_window_size;

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
    cout<<"#--------------------#\n";
    cout<<"LogTimesteps summary:\n";
    cout<<"   alpha="<<alpha<<endl;
    cout<<"   delta="<<delta<<endl;
    cout<<"     npc="<<npc<<endl;
    cout<<" ncycles="<<ncycles<<endl;
    cout<<"      t0="<<t0<<endl;
    cout<<"#--------------------#\n\n";
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
    time_window_size=npc-1+ncycles;
  }

};

#endif
