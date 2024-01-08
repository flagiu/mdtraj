#ifndef _PBC_H_
#define _PBC_H_
using namespace std;
//------------------------------- Radial Distribution Function ----------------------------------//
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

    // apply Minimum Image Convention for general periodic boxes
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
        if(debug) cout << "ALL_IMAGES\n xmax="<<xmax<<"\n ymax="<<ymax<<"\n zmax="<<zmax<<endl;
      }
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
};

#endif
