#ifndef _ED_BOND_PARAMETER_H_
#define _ED_BOND_PARAMETER_H_

#include "neighbour_and_bond_list.hpp"
using namespace std;
//---------------------- Errington-Debenedetti Bond order parameter q -------------------------------//
template <class ntype, class ptype>
class ED_Bond_Parameter
{
  using vec=myvec<ntype,3>;
  using mat=mymatrix<ntype,3,3>;
  private:
    Neigh_and_Bond_list<ntype,ptype> *nb_list;
    vecflex<ntype> q;
    string string_out, myName, tag;
    fstream fout;
    stringstream ss;
    ntype class_boundaries[6] = {1./6., 5./12., 9./16., 11./16., 13./16., 15./16.};
    int class_labels[7] = { -7, -6, -4, -3, -1, -2, -5 };

  public:
    ED_Bond_Parameter(){
      myName = "Errington-Debenedetti BOND PARAMETER";
    }
    virtual ~ED_Bond_Parameter(){}

    void init(Neigh_and_Bond_list<ntype,ptype> *nb_list_, string string_out_, string tag_)
    {
      string_out = string_out_;
      tag = tag_;
      nb_list = nb_list_;
      const int N=nb_list->N;
      q.resize(N);
      ss.str(std::string()); ss << string_out << tag << ".dat"; fout.open(ss.str(), ios::out);
      fout << "#Timestep, Particle index, Coordination number, E-D Bond parameter q. # cutoff = "<<nb_list->rcut[0]<<endl;
      fout.close();
      /*
      ss.str(std::string()); ss << string_out << "_classes"<< tag << ".dat"; fout.open(ss.str(), ios::out);
      fout << "#Timestep, Particle index, E-D class. # cutoff = "<<nb_list->rcut[0]<<endl;
      fout.close();
      */
    }

    void compute(int timestep, vector<ptype> ps, bool debug)
    {
      const int N=nb_list->N;
      const int u=0; // first shell bonds
      ntype rijSq, rikSq, costheta, q_unnormalized;
      vec rij, rik;
      int i,j,k,a,b;
      if(debug) cout << "\n*** "<<myName<<" computation STARTED ***\n";
      for(i=0;i<N;i++)
      {
        q_unnormalized = 0.0;
        for(a=0; a<ps[i].neigh_list[u].size() ;a++)
        {
          //j = ps[i].neigh_list[u][a];
          rij = ps[i].rij_list[u][a];
          rijSq = ps[i].rijSq_list[u][a];
          for(b=a+1;b<ps[i].neigh_list[u].size();b++)
          {
            //k = ps[i].neigh_list[u][b];
            rik = ps[i].rij_list[u][b];
            rikSq = ps[i].rijSq_list[u][b];

            costheta = (rij*rik) / sqrt(rijSq*rikSq);
            q_unnormalized += SQUARE(1./3. + costheta);
          }
        }
        q[i] = 1. - (3./8.)*q_unnormalized;
      }
      ss.str(std::string()); ss << string_out << tag << ".dat"; fout.open(ss.str(), ios::app);
      for(i=0;i<N;i++)
      {
        fout << timestep << " " << i << " " << ps[i].neigh_list[u].size() << " "<< q[i] << endl;
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
