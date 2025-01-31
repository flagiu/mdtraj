#ifndef _PATTERN_MATCHING_PARAMETERS_H_
#define _PATTERN_MATCHING_PARAMETERS_H_

#include "neighbour_and_bond_list.hpp"
using namespace std;
//---------------------- https://doi.org/10.3389/fmats.2017.00034 -------------------------------//
template <class ntype, class ptype>
class PatternMatchingParameters
{
  using vec=myvec<ntype,3>;
  using mat=mymatrix<ntype,3,3>;
  private:
    Neigh_and_Bond_list<ntype,ptype> *nb_list;
    vecflex<ntype> q_tetr, q_oct;
    ntype tetr_angle, dtheta1,dtheta2,thetaSouth, dt1Sq,dt2Sq, rad2deg;
    string string_out, myName, tag;
    fstream fout;
    stringstream ss;
    bool debug, verbose;

  public:
    PatternMatchingParameters(){
      myName = "PatternMatchingParameters";
      tetr_angle = 109.47; //degrees
      dtheta1 = 12.0; // degrees, size of gaussian smoothening
      dtheta2 = 10.0;
      thetaSouth = 160.0; //degrees, threshold for being considered in South Pole
      dt1Sq = dtheta1*dtheta1;
      dt2Sq = dtheta2*dtheta2;
      rad2deg = 180/M_PI;
    }
    virtual ~PatternMatchingParameters(){}

    void init(Neigh_and_Bond_list<ntype,ptype> *nb_list_, string string_out_,
      string tag_, bool debug_, bool verbose_)
    {
      string_out = string_out_;
      tag = tag_;
      debug = debug_;
      verbose = verbose_;
      nb_list = nb_list_;
      const int N=nb_list->N;
      q_tetr.resize(N);
      q_oct.resize(N);
      ss.str(std::string()); ss << string_out << tag << ".dat"; fout.open(ss.str(), ios::out);
      fout << "#Timestep | Particle index  | Particle type | Coordination number | q_tetr | q_oct . # cutoffs = ";
      for(int t=0;t<nb_list->nTypePairs;t++)
        fout <<nb_list->rcut[0][t]<<" ";
      fout<<endl;
      fout.close();
      /*
      ss.str(std::string()); ss << string_out << "_classes"<< tag << ".dat"; fout.open(ss.str(), ios::out);
      fout << "#Timestep, Particle index, E-D class. # cutoff = "<<nb_list->rcut[0]<<endl;
      fout.close();
      */
    }

    void compute(int timestep, vector<ptype> ps)
    {
      const int N=nb_list->N;
      const int u=0; // first shell bonds
      ntype rijSq,rikSq,rimSq, theta_k,theta_m,phi, sumt1,sumt2,sumo1,sumo2;
      vec rij,rik,rim;
      int i,jj,kk,mm, j,k,m, num_neigh;
      if(debug) cout << "\n*** "<<myName<<" computation STARTED ***\n";
      for(i=0;i<N;i++)
      {
        sumt1 = sumo1 = 0.0;
        num_neigh = ps[i].neigh_list[u].size();

        for(jj=0;jj<num_neigh;jj++)
        {
          //j = ps[i].neigh_list[u][jj];
          rij = ps[i].rij_list[u][jj];
          rijSq = ps[i].rijSq_list[u][jj];

          for(kk=0;kk<num_neigh;kk++)
          {
            sumt2=sumo2=0.0;
            if(kk==jj){ continue; }
            //k = ps[i].neigh_list[u][kk];
            rik = ps[i].rij_list[u][kk];
            rikSq = ps[i].rijSq_list[u][kk];
            // angle j-i-k (polar)
            theta_k = rad2deg*acos( (rij*rik) / sqrt(rijSq*rikSq) );

            for(mm=0;mm<num_neigh;mm++)
            {
              if(mm==jj||mm==kk){ continue; }
              //m = ps[i].neigh_list[u][mm];
              rim = ps[i].rij_list[u][mm];
              rimSq = ps[i].rijSq_list[u][mm];
              // angle j-i-m (azimuth)
              theta_m = rad2deg*acos( (rij*rim) / sqrt(rijSq*rimSq) );
              // angle m-i-k (azimuth)
              phi = rad2deg*acos( (rim*rik) / sqrt(rimSq*rikSq) );

              sumt2 += SQUARE(cos(1.5*phi))*exp(-SQUARE(theta_m-tetr_angle)/(2*dt1Sq));

              if(theta_k<thetaSouth && theta_m<thetaSouth){
                sumo2 += SQUARE(cos(2.*phi))*exp(-SQUARE(theta_m-90.)/(2*dt2Sq));
              }

            }

            sumt1 += sumt2 * exp(-SQUARE(theta_k-tetr_angle)/(2*dt1Sq));

            if(theta_k<thetaSouth) {
              sumo1 += sumo2;
            } else {
              sumo1 += 3*exp(-SQUARE(theta_k-180.)/(2*dt1Sq));
            }

          }
        }

        q_tetr[i] = sumt1 / (num_neigh*(num_neigh-1)*(num_neigh-2));
        q_oct[i] = sumo1 / (num_neigh*(3+(num_neigh-2)*(num_neigh-3)));

      }

      ss.str(std::string()); ss << string_out << tag << ".dat"; fout.open(ss.str(), ios::app);
      for(i=0;i<N;i++)
      {
        fout << timestep <<" "<< i <<" "<< ps[i].label <<" "<< ps[i].neigh_list[u].size() << " "<< q_tetr[i] << " "<< q_oct[i] << endl;
      }
      fout.close();
      /*
      // Assign classes to each tetrahedral environment
      ss.str(std::string()); ss << string_out << "_classes" << tag << ".dat"; fout.open(ss.str(), ios::app);
      for(i=0;i<N;i++)
      {
        k=0;
        for(j=0;j<6 && k==0;j++)
        {
          if(q[i]<class_boundaries[j])
          {
            fout << timestep << " " << i << " " << class_labels[j] << endl;
            k=1;
          }
        }
        if(k==0) fout << timestep << " " << i << " " << class_labels[6] << endl;
      }
      fout.close();
      */
      if(debug) cout << "*** "<<myName<<" computation DONE ***\n";
    }

};
#endif
