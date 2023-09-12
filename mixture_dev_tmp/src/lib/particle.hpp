#ifndef _PARTICLE_H_
#define _PARTICLE_H_
#include "myvec.hpp"
#include "params.hpp"

using namespace std;

template <class ntype>
class params_pt: public Params<ntype>
{
public:
  using par=Params<ntype>;
  using par::addParam, par::getParam;
  int label; // integer >=0
  ntype m=1,sigma=1,rcut=1,v; //v=particle volume

  params_pt()
  {
    addParam( new Parameter(FLO,NAME(m)) );
    addParam( new Parameter(FLO,NAME(sigma)) );
    addParam( new Parameter(FLO,NAME(rcut)) );
  }
  virtual void read(string fname) override
  {
    par::read(fname);
    getParam(&sigma, NAME(sigma));
    getParam(&rcut, NAME(rcut));
    v=M_PI*sigma*sigma*sigma/6.0; //sphere
  }
};
template <class ntype>
class particle
{
public:
  using vec3d=myvec<ntype,3>;
  vec3d r;
  vec3d rold, last_dr;
  int last_move; // -1: none; -2 box-move; 0: translation; 1: rotation;
  ntype sigma, sigmaSq, rcut; // for LJ particles
  int idx; // index for array of particles
  int kold,k; // index for linked cells list
  ntype Uold,U; // internal energy
  ntype infty = std::numeric_limits<ntype>::infinity();
  using pt=particle<ntype>;
  int label; // integer >=0
  static const int Nshells=3;
  vector<int> neigh_list[Nshells]; // list of neighbours' index
  vector<vec3d> rij_list[Nshells]; // list of neighbours' distance vector
  vector<ntype> rijSq_list[Nshells]; // list of neighbours' |rij|^2

  particle()
  {
    r << 0,0,0;
    last_dr << 0,0,0;
    last_move=-1;
    sigma=0;
    U=0;
    for(auto u=0;u<Nshells;u++){
      neigh_list[u].clear();
      rij_list[u].clear();
      rijSq_list[u].clear();
    }
  }
  void init(params_pt<ntype>* p)
    {
      //simpars = p;
      sigma= p->sigma;
      sigmaSq = sigma*sigma;
      rcut= p->rcut;
   }
  particle<ntype>& operator=(const particle<ntype>& p1)
    {
      r=p1.r;
      sigma=p1.sigma;
      rcut=p1.rcut;
      sigmaSq=p1.sigmaSq;
      // N.B. index of particle does not change, it is always the same
      // hence it is wrong to copy it
      return (*this);
    }
  ntype ranf()
  {
    return randnum.ranf();
  }

  void show() const
  {
    r.show("r");
  }

  virtual vec3d& get_r()
  {
    return r;
  }
  virtual void set_r(const ntype x, const ntype y, const ntype z)
  {
    r<<x,y,z;
  }
  virtual void rand_box(const vec3d& L) //random pos. in a (-L/2,L/2) box
  {
    r << L[0]*(ranf()-0.5), L[1]*(ranf()-0.5), L[2]*(ranf()-0.5);
  }
  virtual void rand_box(const ntype& L)
  {
    r << L*(ranf()-0.5), L*(ranf()-0.5), L*(ranf()-0.5);
  }
  //// Monte Carlo moves //////////////////////////////
  virtual void traj_move(const ntype del) //random displacement within +-delta
  {
    ntype dx, dy, dz;
    dx = del*(ranf()-0.5);
    dy = del*(ranf()-0.5);
    dz = del*(ranf()-0.5);
    last_dr << dx, dy, dz;
    r += last_dr;
    last_move=0;
  }
  virtual void box_move(ntype Lfact) //rescale all directions
  {
    r*=Lfact;
    last_move=-2;
  }
  ///////////////////////////////////////////////////////////////////
  ///// Periodic Boundary Conditions & Minimum Image Convenction ////
  void pbc(vec3d L) //è corretto?
  {
    for (auto i=0;i<3;i++)
      {
	if(r[i]<-0.5*L[i]) r[i]+=L[i];
	else if(r[i]>0.5*L[i]) r[i]-=L[i];
      }
  }
  void pbc(ntype L)
  {
    for (auto i=0;i<3;i++)
      {
	if(r[i]<-0.5*L) r[i]+=L;
	else if(r[i]>0.5*L) r[i]-=L;
      }
  }
  vec3d mic(vec3d rij, vec3d L)
  {
    vec3d newrij=rij;
    for (auto i=0;i<3;i++)
      newrij[i] -= L[i]*rint(rij[i]/L[i]);
    return newrij;
  }
  vec3d mic(vec3d rij, ntype L)
  {
    vec3d newrij=rij;
    for (auto i=0;i<3;i++)
      newrij[i] -= L*rint(rij[i]/L);
    return newrij;
  }
  //////////////////////////////////////////////////////
  //// Interaction with other particles ////////////////
  virtual ntype pairwise(pt& pB, ntype L)
  {//HARD SPHERES pairwise interaction. Override for different interactions
    return (overlap(pB,L) ? infty : 0.0);
  }
  virtual bool overlap(pt& pB, ntype L)
  {
    ntype rij=(mic(r-pB.r, L)).norm(); // MIC
    return ( rij < (0.5*(sigma+pB.sigma)) );
  }
  bool is_inside(pt& pB, ntype radius, ntype L)
  {// am I inside the given radius from particle B?
    vec3d rij = mic(r-pB.r, L); //MIC
    return (rij.sq() <= radius*radius);
  }
  bool is_inside(pt& pB, ntype radius, vec3d L)
  {// am I inside the given radius from particle B?
    vec3d rij = mic(r-pB.r, L); //MIC
    return (rij.sq() <= radius*radius);
  }
  ///////////////////////////////////////////////////////
  ///// Data storage /////////////////////
  virtual void store_coord(void)
  {
    rold = r;
    kold = k;
  }
  virtual void restore_coord(void)
  {
    r = rold;
    k = kold;
  }
  virtual void store_U(void)
  {
    Uold = U;
  }
  virtual void restore_U(void)
  {
    U = Uold;
  }
  /////////////////////////////////////////
  /////// I/O from stream //////////////////
  virtual void read(fstream& i) // just the position
  {
    for (auto j=0;j<3;j++)
      i >> r[j];
  }
  virtual void write(fstream& o)
  {
    o << setprecision(20) << r[0] << " " << r[1] << " " << r[2] << "\n";
  }
  virtual void read_molgl(fstream& i) // .molgl format
  {
    char at;
    for (auto j=0;j<3;j++)
      i >> r[j];
    i >> at;
    i >> sigma;
  }
  virtual void write_molgl(fstream& o)
  {
    o << setprecision(20) << r[0] << " " << r[1] << " " << r[2] << " @ " << sigma << endl;
  }

  virtual void read_xyz(fstream& i) // .xyz format
  {
    string line, lab, a[3];
    getline(i, line);
    istringstream(line) >> lab >> a[0] >> a[1] >> a[2];
    label = stoi(lab);
    for(auto j=0;j<3;j++) r[j] = stof(a[j]);
  }
  virtual void write_xyz(fstream& o)
  {
    o << setprecision(20) << label << " " << r[0] << " " << r[1] << " " << r[2] << endl;
  }

  virtual void read_3cols(fstream& i) // .xyz format
  {
    string line, a[3];
    getline(i, line);
    istringstream(line) >> a[0] >> a[1] >> a[2];
    for(auto j=0;j<3;j++) r[j] = stof(a[j]);
  }
  virtual void write_3cols(fstream& o)
  {
    o << setprecision(20) << r[0] << " " << r[1] << " " << r[2] << endl;
  }


  virtual void write_pdb(fstream& o, int myIndex, ntype myVal)
  {
    o << setprecision(10) << "ATOM \t " << myIndex+1 << " X\tXXX X " << label << "\t" << r[0] << " " << r[1] << " " << r[2] << " 1.00 " << myVal << "\t X\n";
  }
  /////////////////////////////////////////
  virtual ~particle()
  {}
};
#endif
