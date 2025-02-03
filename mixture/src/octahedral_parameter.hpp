#ifndef _QOCT_H_
#define _QOCT_H_

#include "neighbour_and_bond_list.hpp"
using namespace std;
//---------------------- Custom octahedral order parameter -------------------------------//
// For the ideal simple cubic: CoordinationNumber=6, q_oct=1
// For a different phase with CN=6: q_oct<1
template <class ntype, class ptype>
class OctahedralParameter
{
  using vec=myvec<ntype,3>;
  using mat=mymatrix<ntype,3,3>;
  private:
    Neigh_and_Bond_list<ntype,ptype> *nb_list;
    vecflex<ntype> my_q_oct;
    ntype my_cos_th;
    string string_out, myName, tag;
    fstream fout;
    stringstream ss;
    bool debug, verbose;

  public:
    OctahedralParameter(){
      myName = "OctahedralParameter";
      my_cos_th = cos(3./4.*M_PI); // negative! Threshold to decide if closer to 90° or 180°
    }
    virtual ~OctahedralParameter(){}

    void init(Neigh_and_Bond_list<ntype,ptype> *nb_list_, string string_out_,
      string tag_, bool debug_, bool verbose_)
    {
      string_out = string_out_;
      tag = tag_;
      debug = debug_;
      verbose = verbose_;
      nb_list = nb_list_;
      const int N=nb_list->N;
      my_q_oct.resize(N);
      ss.str(std::string()); ss << string_out << tag << ".dat"; fout.open(ss.str(), ios::out);
      fout << "#Timestep | Particle index  | Particle type | q_oct . # cutoffs = ";
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
      ntype rijSq,rikSq, sum_my_oct;
      ntype cos_k;
      vec rij,rik;
      int i,jj,kk, j,k, num_neigh;
      if(debug) cout << "\n*** "<<myName<<" computation STARTED ***\n";
      for(i=0;i<N;i++)
      {
        sum_my_oct = 0.0;
        num_neigh = ps[i].neigh_list[u].size();
        if(num_neigh<6) {
          cerr<<"ERROR: less than 6 nearest neighbors. Increase the cutoff radius.\n";
          exit(1);
        }

        for(jj=0;jj<6;jj++)
        {
          //j = ps[i].neigh_list[u][jj];
          rij = ps[i].rij_list[u][jj];
          rijSq = ps[i].rijSq_list[u][jj];

          for(kk=jj+1;kk<6;kk++)
          {
            //k = ps[i].neigh_list[u][kk];
            rik = ps[i].rij_list[u][kk];
            rikSq = ps[i].rijSq_list[u][kk];
            // angle j-i-k (polar)
            cos_k = (rij*rik) / sqrt(rijSq*rikSq);

            if(cos_k<my_cos_th) { // if closer to 180°
              sum_my_oct+=SQUARE(cos_k+1) / 3.0;
            } else { // else closer to 90°
              sum_my_oct+=SQUARE(cos_k) / 12.0;
            }

          }
        }

        my_q_oct[i] = 1.0-sum_my_oct;
      }

      ss.str(std::string()); ss << string_out << tag << ".dat"; fout.open(ss.str(), ios::app);
      for(i=0;i<N;i++)
      {
        fout << timestep <<" "<< i <<" "<< ps[i].label << " " << my_q_oct[i]<< endl;
      }
      fout.close();
      if(debug) cout << "*** "<<myName<<" computation DONE ***\n";
    }

};
#endif
