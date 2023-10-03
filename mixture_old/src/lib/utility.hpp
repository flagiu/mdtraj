#ifndef _UTILITY_H_
#define _UTILITY_H_

#include<cstdlib>
#include<fstream>
#include<iostream>
#include<iomanip>
#include<vector>
using namespace std;

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
int indexOf(vector<T>& vec, T element)
{
  for(int i=0, size=vec.size(); i<size; i++)
    if(vec[i]==element) return i;
  return -1; // not found
}

class PrintProgress
{
private:
  int perc; // in percentage units
  int total_steps, nblocks, block_size, current_block;
public:
  void init(int tot)
  {
    total_steps = tot;
    nblocks = 10;
    if(total_steps<nblocks) nblocks=total_steps;
    block_size = total_steps/nblocks;
    current_block=0;
    perc=0;
  }
  void update(int step)
  {
    if( step-current_block*block_size >= block_size)
    {
      current_block++;
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
  float stop()
  {
    end = clock();
    return (float) (end-start) / CLOCKS_PER_SEC * 1000.0;
  }
};

#endif