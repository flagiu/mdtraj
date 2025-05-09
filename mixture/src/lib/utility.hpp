#ifndef _UTILITY_H_
#define _UTILITY_H_

#include<cstdlib>
#include<fstream>
#include<iostream>
#include<iomanip>
#include<vector>
using namespace std;

#define SQUARE(x) ((x)*(x))

/*
template<class T>
void alloc(T *ptr, int n) {
  if(n<1){
    cerr<<"Error: cannot allocate n="<<n<<" elements\n";
    exit(1);
  }
  ptr = new (nothrow) T [n];
  if (ptr == nullptr) {
    cerr<<"Error during allocation of n="<<n<<" elements of type T="<<typeid(T).name()<<"\n";
    exit(1);
  }
}

template<class T>
void malloc(T *ptr, int n) {
  alloc(ptr, n);
  for(int i=0;i<n;i++) ptr[i]=0;
}
*/

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

bool is_positive_integer(const std::string& s)
{
    std::string::const_iterator it = s.begin();
    while (it != s.end() && std::isdigit(*it)) ++it;
    return !s.empty() && it == s.end();
}

#endif
