#ifndef _MOLECULE_H_
#define _MOLECULE_H_
#include "particle.hpp"
#include "mymatrix.hpp"

template <class ntype>
class params_mt: public params_pt<ntype>
{
public:
  using ppt = params_pt<ntype>;
  using ppt::addParam, ppt::getParam, ppt::read;
  //using ppt::r, ppt::w, ppt::read, ppt::write, ppt::show;
  ntype a,b,c;
  int rotmove=3;
  params_mt()
  {
    addParam( new Parameter(FLO,NAME(a)) );
    addParam( new Parameter(FLO,NAME(b)) );
    addParam( new Parameter(FLO,NAME(c)) );
    addParam( new Parameter(INT,NAME(rotmove)) );
  }
  virtual void read(string fname) override
  {
    ppt::read(fname);
    getParam(&a, NAME(a));
    getParam(&b, NAME(b));
    getParam(&c, NAME(c));
    getParam(&rotmove, NAME(rotmove));
  }
};
template <class ntype>
class molecule: public particle<ntype>
{
public:
  using vec3d=myvec<ntype,3>;
  using vec2d=myvec<ntype,2>;
  using pt = particle<ntype>;
  using pt::r, pt::last_dr, pt::ranf, pt::show, pt::sigma, pt::rcut, pt::traj_move, pt::rand_box, pt::mic, pt::read, pt::write;
  matrix3d R,Rold; // rotation matrix: R[0]=x', R[1]=y', R[2]=z' trasformate di x,y,z dal s.d.r. della molecola al s.d.r. del laboratorio
  ntype a,b,c; // semiaxes (as hard ellipsoids)
  int rotmove; // 3: 3d rotation; 4: 4d rotation
  int dir1ax; // index of the direction: 0,1,2 if uniaxial (1ax); -1 if spherical; -2 otherwise.
  ntype infty = std::numeric_limits<ntype>::infinity();
  using myptype=molecule<ntype>;

  molecule()
    {
      R << 1,0,0,0,1,0,0,0,1;
      rotmove=3;
      dir1ax=-1;
      a=b=c=1;
    }
  void init(params_mt<ntype>* p) 
    {
      pt::init(p);
      a=p->a;
      b=p->b;
      c=p->c;
      if (a==b && b==c) dir1ax=-1; //0ax - spherically simmetric
      else if (a==b) dir1ax=2; //1ax - uniaxial
      else if (b==c) dir1ax=0; //1ax
      else if (c==a) dir1ax=1; //1ax
      else dir1ax=-2; // 2ax - others
      rotmove=p->rotmove;
    }
  molecule<ntype>& operator=(const molecule<ntype>& p1)
    {
      pt::operator=(p1);
      R=p1.R;
      a=p1.a;
      b=p1.b;
      c=p1.c;
      dir1ax=p1.dir1ax;
      return (*this);
    }
  void store_coord(void) override
    {
      pt::store_coord();
      Rold = R;
    }

  void restore_coord(void) override
    {
      pt::restore_coord();
      R = Rold;
    }
  ///////////////////////////////////////////
  ///// Operations on rotation matrix //////
  void normalize_R()
    {
      vec3d u;
      for (auto i=0;i<3;i++)
	{
	  u=R[i];
	  u.normalize();
	  R[i]=u;
	}
    }
  void orthonormalize_R(int start_idx)
    {
      int i=start_idx, j=(i+1)%3, k=(i+2)%3;
      vec3d u1=R[i], u2=R[j], u3;
      u1.normalize();
      u2 = u2 - (u2*u1)*u1;
      u2.normalize();
      u3 = u1.cross(u2);
      u3.normalize();
      R.set(i,u1);
      R.set(j,u2);
      R.set(k,u3);
    }
  //////////////////////////////////////////////
  //// Orientation management //////////////////
  void set_1ax(const ntype x, const ntype y, const ntype z)
  {
    vec3d u;
    u<<x,y,z;
    u.normalize();
    R.set(dir1ax,u);
    orthonormalize_R(dir1ax);
  }
  vec3d random_orient(void); //declare external function?? (Marsaglia algorithm)
  void rand_R() //generate a new random orientation
  {
    int axis;
    if (dir1ax>-1) axis=dir1ax; //1ax
    else if (dir1ax==-2) axis=(int)(3*ranf()); //2ax
    vec3d u=random_orient();
    R[axis]=u;
    orthonormalize_R(axis);
  }
  ////////////////////////////////////////////////////
  //// MC rotational move //////////////////////////
  virtual void rot_move(ntype delta)
  {
    if (dir1ax==-1) return;
    else if (rotmove==1 && dir1ax > -1)
      rot_move_1ax(delta);
    else if (rotmove==3)
      rot_move_3d(delta);
    else if (rotmove==4)
      rot_move_4d(delta);
    else
      cout << "Error: invalid rotational move value: rotmove="<<rotmove<<endl;
    pt::last_dr << 0,0,0;
    pt::last_move=1;
  }
  void rot_move_1ax(ntype delta)
  {
    vec3d du = random_orient(); //random versor with Marsaglia algorithm
    vec3d u=R[dir1ax];
    u += (delta*du); //add delta*du to the uniaxial direction
    R.set(dir1ax,u);
    orthonormalize_R(dir1ax); // orthonormalize R w.r.t the uniax. direction
  }
  void rot_move_3d(ntype delta)
  {
    ntype th=delta*(randnum.ranf()-0.5); //random angle
    vec3d u=random_orient(); //random direction
    matrix3d Om, Om2, M;
    Om << 0,-u[2],u[1], u[2],0,-u[0], -u[1],u[0],0;
    Om*=(-sin(th));
    Om2=(Om*Om);
    Om2*=(1-cos(th));
    M = Om;
    M+= Om2;
    R += (R*M); // è normalizzata? Numericamente sembra di no
  }
  void rot_move_4d(ntype delrot)
  {
    //use a unitary quaternion rotation
    cout<<"4D rotation NOT YET IMPLEMENTED\n";
  }
  //////////////////////////////////////////////
  ////// Placement in a box //////////////////
  void rand_box(const vec3d& L) override
  {
    pt::rand_box(L);
    if(dir1ax!=-1) rand_R();
  }
  void rand_box(const ntype& L) override
  {
    pt::rand_box(L);
    if(dir1ax!=-1) rand_R();
  }
  ///////////////////////////////////////////////
  
  vec3d get_sax(void)
  {
    vec3d sax;
    sax<<a,b,c;
    return sax;
  }
  virtual void show() const
  {
    pt::show();
    R.show("R");
  }
  virtual void read(fstream &in) override
    {
      pt::read(in);
      for (auto i=0; i < 3; i++)
        for (auto j=0; j < 3; j++)
          in >> R[i][j]; 
    }
  virtual void write(fstream &o) override
    {
      pt::write(o);
      for (auto i=0; i < 3; i++)
        for (auto j=0; j < 3; j++)
          {
            o << setprecision(20) << R[i][j];
            if (i == 2 && j == 2)
              o << "\n";
            else
              o << " ";
          }
    }
  virtual ~molecule()
  {}
};

template <class ntype>
myvec<ntype,3> molecule<ntype>::random_orient(void)
{// Marsaglia algorithm: sample a vector from unitary sphere
  myvec<ntype,2> x;
  ntype n,xi;
  myvec<ntype,3> u;
  do{
    x.rand_box();
    x*=2.0;
    n=x.norm();
  } while(n>=1.0);
  xi=sqrt(fabs(1.0-n));
  u << 2*x[0]*xi, 2*x[1]*xi, 1-2*n;
  return u;
}
#endif
