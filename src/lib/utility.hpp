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
public:
  void init(float p)
  {
    perc = floor(p);
  }
  void update(float p)
  {
    if(p-perc >= 1.0)
    {
      perc = floor(p);
      cout << '\r' << setfill('0') << "[" << perc << "%]";
    }
  }
};

#endif
