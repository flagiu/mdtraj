#ifndef _LL_H_
#define _LL_H_
#include "molecule.hpp"
using namespace std;
// Linked cells List for a cubic box

template <class ntype, class ptype>
class linked_list
{
public:
  int debug=0;
  int N,M,C;
  ntype l;
  vector<int> llist, head;//,  ll_old, h_old;

  template<class ptype2>
  linked_list<ntype,ptype>&  operator=(const linked_list<ntype,ptype2>& other)
  {// copy all
    N = other.N;
    M = other.M;
    C = other.C;
    l = other.l;
    llist.resize(N);
    head.resize(C);
    if (N>C)
      for (auto i=0;i<N;i++)
	{
	  llist[i]= other.llist[i];
	  if (i<C) head[i]= other.head[i];
	}
    else
      for (auto i=0;i<C;i++)
	{
	  head[i]= other.head[i];
	  if(i<N) llist[i] = other.llist[i];
	}
  }

  void init_pars(int N_, ntype rc, ntype L, int dbg)
  {
    debug=dbg;
    if(debug) cout<<"Initializing LL ";
    N=N_;            // n. of particles
    M = floor(L/rc); // n. of cells per side
    C=M*M*M;         // n. of total cells
    l = L/(ntype)M;  // cell side: l>=rc
    llist.resize(N);
    head.resize(C);
    if(debug) cout<<"with M = "<<M<<", l = "<<l<<" ...";
    //ll_old.resize(N);
    //h_old.resize(C);
    if (N>C)
      for (auto i=0;i<N;i++)
	{
	  llist[i]=-1;
	  if (i<C) head[i]=-1;
	}
    else
      for (auto i=0;i<C;i++)
	{
	  head[i]=-1;
	  if(i<N) llist[i]=-1;
	}
    if(debug) cout<<" done.\n";
  }
  void show(fstream& o)
  {
    o<<"head:\n";
    for (auto i=0;i<C;i++)
      o<<i<<" "<<head[i]<<endl;
    o<<"####################\n";
    o<<"list:\n";
    for (auto i=0;i<N;i++)
      o<<i<<" "<<llist[i]<<endl;
  }
  void show(int k0)
  {
    myvec<unsigned int,3> c0=get_coord(k0);
    int k, i;
    for(int x=c0[0]-1;x<=c0[0]+1;x++)
      {
	for(int y=c0[1]-1;y<=c0[1]+1;y++)
	  {
	    for(int z=c0[2]-1;z<=c0[2]+1;z++)
	      {
		k=get_idx(x,y,z);
		cout<<"Cell "<<k<<" ("<<x<<","<<y<<","<<z<<")"<<endl;
		i=head[k];
		while(i>-1)
		  {
		    cout<<" i+1 = "<<i+1<<endl;
		    i=llist[i];
		  }
	      }
	  }
      }
  }
  int get_idx(myvec<ntype,3> r)
  {
    int i,j,k;
    // (int)x cuts the decimals: (int)x = floor(x) if x>0, else ceil(x)
    i = (int)(r[0]/l+0.5*M); // x-coordinate in the 3D grid
    j = (int)(r[1]/l+0.5*M); // y
    k = (int)(r[2]/l+0.5*M); // z
    return get_idx(i,j,k);
  }
  int get_idx(int i, int j, int k)
  {
    return M*M*((i+M)%M) + M*((j+M)%M) + ((k+M)%M); // account for i,j,k <0 or >M-1
  }
  myvec<unsigned int,3> get_coord(int k)
  {// return the cell coordinates (i.e. write k in basis M)
    myvec<unsigned int,3> x;
    int a,b,c;
    a = k/(M*M); //'centinaia'
    c = k-a*(M*M);
    b = c/M; //'decine'
    c -= b*M; //'unità'
    x << a,b,c;
    return x;
  }
  void remove(int i, int k)
  {// remove particle i from cell k
    int j = head[k], j_old;
    if (j==i) // if i is the first one:
	head[k] = llist[i]; // replace it with its follower
    else
      {// else (i is not the first of the list)
	do
	  {
	    //cout<<"j="<<j<<endl;
	    j_old=j;
	    j = llist[j_old];
	  }
	while (j!=i && j>-1); // search through the list of particles in k, until i
	if (j==i) llist[j_old] = llist[i]; // replace the follower of j with the follower of i
	else
	  {
	    cout<<"!!! Problem in LL: searching i="<<i<<" in k="<<k<<", reached end of the list !!!\n";
	    //show();
	    exit(EXIT_FAILURE);
	  }
      }
    llist[i] = -1; // set its follower to none (useless if remove is followed by insert)
  }
  void insert(int i, int k)
  {// insert particle i into cell k
    llist[i] = head[k]; // set the first one as the follower of i
    head[k] = i; // set the particle as the first one
  }
  void build(vector<ptype>& ps)
  {
    int k;
    if (debug) cout<<"Building LL...";
    for (auto n=0;n<N;n++)
      {
        k = get_idx(ps[n].r);
	ps[n].k = k;
	insert(n,k);
      }
    if (debug) cout<<" done.\n";
  }
  ~linked_list()
  {
  }
};
#endif
