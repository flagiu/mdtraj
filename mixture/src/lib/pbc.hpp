#ifndef _PBC_H_
#define _PBC_H_
using namespace std;

#define MAX_N_IMAGES 300

template <class ntype>
class PBC
{
  using vec=myvec<ntype,3>;
  using mat=mymatrix<ntype,3,3>;
  private:
    int image_convention=1; // 0: no images (cluster) ; 1: 1 replica (minimum image) ; -1: all replicas (crystal)
    ntype rmax;
    mat box, boxInv;
    int xmax,ymax,zmax, num_vectors;
    vec ds_within_rmax[MAX_N_IMAGES];
    const string myName="PBC";

    // Minimum Image Convention for general periodic boxes
    vec mic(vec& r)
    {
      vec r_mic = r - box * round(boxInv*r);
      return r_mic.norm()<r.norm() ? r_mic : r;
    }

    int all_images(vec r, vec* r_images)
    {
      int x,y,z, count=0;
      vec s, ds, trial;
      s = boxInv*r;
      for(x=-xmax;x<=xmax && count<MAX_N_IMAGES;x++)
      {
        for(y=-ymax;y<=ymax && count<MAX_N_IMAGES;y++)
        {
          for(z=-zmax;z<=zmax && count<MAX_N_IMAGES;z++)
          {
            ds << x,y,z;
            trial = box*(s+ds);
            if(trial.norm()<rmax)
            {
              *(r_images+count) = trial;
              count++;
            }
          }
        }
      }
      return count;
    }

  public:
    PBC(){}
    virtual ~PBC(){}

    void init(int img_conv, ntype rmax_, mat box_, bool debug)
    {
    if(debug) cerr << myName << " Initialization STARTED\n";
      image_convention = img_conv;
      rmax = rmax_;
      box = box_;
      boxInv = box.inverse();
      if(image_convention==-1)
      {
        int x,y,z, count=0;
        vec trial;
        xmax = round(rmax/fabs(box[0][0]+box[0][1]+box[0][2]));
        ymax = round(rmax/fabs(box[1][0]+box[1][1]+box[1][2]));
        zmax = round(rmax/fabs(box[2][0]+box[2][1]+box[2][2]));
        if(debug) cerr << "ALL_IMAGES\n xmax="<<xmax<<"\n ymax="<<ymax<<"\n zmax="<<zmax<<endl;
      }
      if(debug) cerr << myName << " Initialization COMPLETED\n";
      return;
    }

    int apply(vec r, vec* r_images)
    { // r_images must be allocated before this call!!!
      int ret;
      switch (image_convention)
      {
        case 1:
          (*r_images) = mic(r);
          ret=1;
          break;
        case 0:
          (*r_images) = r;
          ret=1;
          break;
        case -1:
          ret=all_images(r,r_images);
          break;
      }
      return ret;
    }

    // maps r_unwrapped --> [-L/2,L/2]
    // Example: x=L*0.9 --> L*0.9 - L*round(L*0.9/L) = L*0.9 - L*1 = -L*0.1
    void wrap_inplace(vec* r){
      vec r_wrapped = *r - box * round(boxInv*(*r));
      *r = r_wrapped;
    }
    vec wrap_outplace(vec* r){
      vec r_wrapped = *r - box * round(boxInv*(*r));
      return r_wrapped;
    }

    // unwraps assuming that |x_true-x_old| < L/2
    // Example: x_old=-.4*L ; x=.4*L
    // --> dx=.8*L
    //  --> dx'=.8*L-L*round(.8)=-0.2*L
    //   --> x'=x_old+dx'=-0.6*L
    // Example: x_old=12*L ; x=15.2*L
    // --> dx=3*L
    //  --> dx'=L*(3.2-round(3.2))=.2*L
    //   --> x'=x_old+dx'=12.2*L
    void unwrap_inplace(vec* r, vec* r_old){
      vec dr = *r - *r_old;
      dr -= box * round(boxInv*dr);
      *r = *r_old+ dr;
    }
    vec unwrap_outplace(vec* r, vec* r_old){
      vec dr = *r - *r_old;
      dr -= box * round(boxInv*dr);
      return *r_old + dr;
    }
};

#endif
