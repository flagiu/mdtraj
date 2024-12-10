#ifndef _DIFFUSION_H_
#define _DIFFUSION_H_

#include<cstdlib>
#include<fstream>
#include<iostream>
#include<iomanip>
#include<vector>
using namespace std;

#define SQUARE(x) ((x)*(x))

void chec_null_ptr(void *ptr){
  if (ptr==NULL) exit (1);
  return;
}

template <class ntype>
void DiffusionFlux1d(ntype* vs, int len, int nt, int origin_interval, ntype delta, bool debug)
{
  /*
  INPUT:
               vs: pointer to a 1d array of flux "velocity"
              len: length of vs
               nt: window size for correlation (<=len)
  origin_interval: shift of origins for correlation (<=len)
            delta: time interval btw data points
  OUTPUT:
     t, vacf, rvcf, msd
  */
  if(nt>len || nt<=0){
    cerr<<"ERROR: invalid nt="<<nt<<" ; len="<<len<<endl;
    exit(1);
  }
  if(origin_interval>len || origin_interval<=0){
    cerr<<"ERROR: invalid origin_interval="<<origin_interval<<" ; len="<<len<<endl;
    exit(1);
  }
  ntype *rs_buffer, *msd, *rvcf, *vacf, *norm,*v0,*r0;
  ntype r,v;
  int *t0, n0, mk, t, nk, dt, k;
  int rs_buf_size,i1,i0;
  bool full;

  rs_buf_size=2; // for Euler integration you only need the previous time

  rs_buffer = (ntype*)calloc(rs_buf_size, sizeof(ntype));
  chec_null_ptr(rs_buffer);

  msd = (ntype*)calloc(nt+1, sizeof(ntype));
  chec_null_ptr(msd);
  rvcf = (ntype*)calloc(nt+1, sizeof(ntype));
  chec_null_ptr(rvcf);
  vacf = (ntype*)calloc(nt+1, sizeof(ntype));
  chec_null_ptr(vacf);
  norm = (ntype*)calloc(nt+1, sizeof(ntype));
  chec_null_ptr(norm);

  // number of time origins
  n0=nt/origin_interval+1;
  mk=-1;
  full=false;

  t0 = (int*)malloc(n0*sizeof(int));
  chec_null_ptr(t0);
  v0 = (ntype*)malloc(n0*sizeof(ntype));
  chec_null_ptr(v0);
  r0 = (ntype*)malloc(n0*sizeof(ntype));
  chec_null_ptr(r0);

  for(t=0;t<len;t++)
  {
    v=vs[t];

    // simple Euler integration
    i1=t%rs_buf_size;
    if(t==0){
      r=0.0;
    }else{
      i0=(i1-1+rs_buf_size)%rs_buf_size;
      r = rs_buffer[i0] + vs[t-1]*delta;
    }
    rs_buffer[i1]=r;

    // Test to store as time origin
    if(t%origin_interval == 0)
    {
      if(debug){ cerr<<"\rProcessing t="<<t<<"      "; }
      mk = mk + 1;
      if(mk >= n0){
        full = true;
        mk = mk - n0; //Overwrite older values
      }
      t0[mk] = t; // Store time origin
      r0[mk] = r; // Store position at time origin
      v0[mk] = v; // Store velocity at time origin
    }

    // Correlate with all time origins stored so far
    nk = full ? n0 : mk+1;
    for(k=0;k<nk;k++){ // Loop over time origins
        dt=t-t0[k];
        if(dt<0){
          cerr<<"ERROR: dt="<<dt<<endl;
          exit(1);
        }
        if(dt<=nt){ // Check that dt is in range
          msd[dt]  += SQUARE(r-r0[k]);        // Increment msd
          rvcf[dt] += (r-r0[k])*v0[k]; // Increment cross correlation function
          vacf[dt] += v*v0[k];             // Increment autocorrelation function
          norm[dt] += 1.; // Increment normalizing factor
        }
    }

  }

  // Normalize
  for(t=0;t<nt+1;t++){
    if(norm[t]<1){
      cerr<<"ERROR: wrong normalization\n";
      exit(1);
    }
    msd[t]  /= norm[t];
    rvcf[t] /= norm[t];
    vacf[t] /= norm[t];
  }

  // Output
  for(t=0;t<nt+1;t++){
    printf("%.6f %15.8f %15.8f %15.8f\n",
      t*delta, vacf[t], rvcf[t], msd[t]
    );
  }
}


template <class ntype>
void DiffusionFlux1d_FromFile(string file, int nt, int origin_interval, ntype delta, bool debug)
{
  ifstream i;
  string line;
  int count=0, len;
  ntype x, *arr;
  const int maxcount=99999999;

  i.open(file, ios::in);
  if (!i.is_open()){
    cout << "ERROR: Unable to open file "<<file<<endl;
    exit(1);
  }
  // count lines
  while ( count<maxcount && getline(i,line)){
    x=stof(line);
    count++;
  }
  i.close();

  len=count;
  arr=(ntype*)malloc(len*sizeof(ntype));
  chec_null_ptr(arr);

  count=0;
  i.open(file, ios::in);
  // get_data
  while ( count<len && getline(i,line)){
    arr[count++]=stof(line);
  }
  i.close();

  return DiffusionFlux1d(arr, len, nt, origin_interval, delta, debug);
}

using ntype=double;
int main(int argc, char* argv[])
{
  if(argc<5||argc>6){
    cerr<<"Input: <file> <nt> <oi> <delta> [debug=0]\n";
    exit(1);
  }
  int nt, oi;
  ntype delta;
  bool debug=false;
  nt=stoi(argv[2]);
  oi=stoi(argv[3]);
  delta=stof(argv[4]);
  if(argc>=6) debug=stoi(argv[5]);
  DiffusionFlux1d_FromFile(argv[1], nt, oi, delta, debug);
  return 0;
}
#endif
