#ifndef _PARAMS_H_
#define _PARAMS_H_
#define NAME(x) (#x)
#define INT 0
#define FLO 1
#define ULI 2

#include<cstdlib>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<list>
using namespace std;

class Parameter
{
public:
  string name;
  int int_val=111;
  float flo_val=222;
  unsigned long int uli_val=333;
  int datatype;
  
  Parameter(int type, const char *nam)
  {
    name=nam;
    datatype=type;
  }
  template<class dtype>
  void setVal(dtype val)
  {
    switch(datatype){
    case INT: int_val=val;
      break;
    case FLO: flo_val=val;
      break;
    case ULI: uli_val=val;
      break;
    }
  }
  template<class dtype>
  void getVal(dtype* ptr)
  {
    switch(datatype){
    case INT: (*ptr) = int_val;
      break;
    case FLO: (*ptr) = flo_val;
      break;
    case ULI: (*ptr) = uli_val;
      break;
    }
  }

  void valFromStr(string str)
  {
    switch(datatype){
    case INT: int_val = stoi(str);
      break;
    case FLO: flo_val = stof(str);
      break;
    case ULI: uli_val = stoul(str);
      break;
    }
  }
  
  void show()
  {
    switch(datatype){
    case INT: cout << name << " " << int_val << endl;
      break;
    case FLO: cout << name << " " << flo_val << endl;
      break;
    case ULI: cout << name << " " << uli_val << endl;
      break;
    }
  }
  void write(fstream& o)
  {
    switch(datatype){
    case INT: o << name << " " << int_val << endl;
      break;
    case FLO: o << name << " " << flo_val << endl;
      break;
    case ULI: o << name << " " << uli_val << endl;
      break;
    }
  }
};

template <class ntype>
class Params
{
public:
  int debug, verbose;
  string trash;
  list<Parameter*> params;
  Params()
  {
    debug=0;
    verbose=0;
  }
  virtual void show()
  {
    cout << NAME(debug) << " " << debug << endl;
    cout << NAME(verbose) << " " << verbose << endl;
    for(auto param : params)
      param->show();
  }
  virtual void read(fstream& in)
  {
  }
  Parameter *searchParam(const char *nam)
  {
    for(auto param : params)
      if((param->name)==nam)
	  return param;
    {cout<<"Error: parameter "<<nam<<" not found.\n"; exit(1);}
  }
  void addParam(Parameter *par)
  {
    params.push_back(par);
  }
  template<class dtype>
  void getParam(dtype* ptr, const char *nam)
  {
    Parameter* p=searchParam(nam);
    p->getVal(ptr);
  }
  template<class dtype>
  void setParam(dtype val, const char *nam)
  {
    Parameter* p=searchParam(nam);
    p->setVal(val);
  }
  virtual void readParam(string fname, const char* nam)
  {// reads a single parameter from file fname, in any order!!
    Parameter* p=searchParam(nam);
    ifstream file(fname);
    if(file.is_open()) {
      string line, name, val;
      while(getline(file,line))
	{//proper format: "name value" (otherwise ignore line)
	  if ((istringstream(line) >> name >> val))
	    if(name==p->name)
	      p->valFromStr(val);
	}
    }
  }
  virtual void read(string fname)
  {// reads all parameters from file fname, in any order!!
    ifstream file(fname);
    if(file.is_open()) {
      string line, name, val;
      while(getline(file,line))
	{//proper format: "name value" (otherwise ignore line)
	  if ((istringstream(line) >> name >> val))
	    for(auto param: params)
	      if(name==param->name)
		param->valFromStr(val);
	}
    }
  }
  
  virtual void write(fstream& o)
  {
    for(auto param : params)
      param->write(o);
  }
};
#endif
