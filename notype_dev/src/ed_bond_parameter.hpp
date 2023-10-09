#ifndef _ED_BOND_PARAMETER_H_
#define _ED_BOND_PARAMETER_H_

#include "neighbour_and_bond_list.hpp"
using namespace std;
//---------------------- Eddington-Debenedetti Bond order parameter q -------------------------------//
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

  public:
    ED_Bond_Parameter(){
      myName = "Eddington-Debenetetti BOND PARAMETER";
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
      fout << "#Timestep, Particle index, E-D Bond parameter q. # cutoff = "<<nb_list->rcut[0]<<endl;
      fout.close();
      ss.str(std::string()); ss << string_out << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "#Timestep, <q>, fluctuations. # cutoff = "<<nb_list->rcut[0]<<endl;
      fout.close();
    }

    void compute(int timestep, vector<ptype> ps, bool debug)
    {
      const int N=nb_list->N;
      const int u=0; // first shell bonds
      ntype rijSq, rikSq, costheta, q_unnormalized, Q,Q2;
      vec rij, rik;
      int i,j,k,a,b;
      if(debug) cout << "\n*** "<<myName<<" computation STARTED ***\n";
      for(i=0;i<N;i++)
      {
        q_unnormalized = 0.0;
        for(a=0;a<min(3,int(ps[i].neigh_list[u].size()-1) );a++)
        {
          //j = ps[i].neigh_list[u][a];
          rij = ps[i].rij_list[u][a];
          rijSq = ps[i].rijSq_list[u][a];
          for(b=a+1;b<min(4,int(ps[i].neigh_list[u].size()) );b++)
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
      Q = 0.0; //  <q>
      Q2 = 0.0; // <q^2>
      for(i=0;i<N;i++)
      {
        fout << timestep << " " << i << " " << q[i] << endl;
        Q += q[i];
        Q2 += q[i]*q[i];
      }
      fout.close();
      Q /= N;
      Q2 /= N;
      ss.str(std::string()); ss << string_out << tag << ".ave"; fout.open(ss.str(), ios::app);
      fout << timestep << " " << Q << " " << sqrt((Q2-Q*Q)/N) << endl;
      fout.close();
      if(debug) cout << "*** "<<myName<<" computation DONE ***\n";
    }

};
#endif
