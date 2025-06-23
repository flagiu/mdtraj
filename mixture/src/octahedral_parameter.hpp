#ifndef _QOCT_H_
#define _QOCT_H_

#include "neighbour_and_bond_list.hpp"
using namespace std;
//---------------------- Custom octahedral order parameter -------------------------------//
// For the ideal simple cubic: CoordinationNumber=6, q_oct=0
// For a different phase with CN=6: q_oct>0
template <class ntype, class ptype>
class OctahedralParameter
{
  using vec=myvec<ntype,3>;
  using mat=mymatrix<ntype,3,3>;
  private:
    Neigh_and_Bond_list<ntype,ptype> *nb_list;
    ntype rad2deg;
    string string_out, myName, tag;
    fstream fout;
    stringstream ss;
    bool debug, verbose;

  public:
    vector<ntype> my_q_oct, counts_collinear, ratio_sl;
    ntype my_cos_th, ALTBC_cos_th, normalization;
    bool continuous;

    OctahedralParameter(){
      myName = "OctahedralParameter";
      my_cos_th = -2./3.; cos(3./4.*M_PI); // negative! Threshold to decide if closer to 90°or 180°
      normalization = 324./(13.*15.); // = 108/65;  324/13 comes from the addend; 1/15 comes from te number of angles
      rad2deg = 180./M_PI;
      ALTBC_cos_th = cos((180.-25.)/rad2deg);
      continuous=false;
    }
    virtual ~OctahedralParameter(){}

    void init(Neigh_and_Bond_list<ntype,ptype> *nb_list_, string string_out_, bool continuous_,
      string tag_, bool debug_, bool verbose_)
    {
      string_out = string_out_;
      tag = tag_;
      debug = debug_;
      verbose = verbose_;
      continuous = continuous_;
      nb_list = nb_list_;
      const int N=nb_list->N;
      my_q_oct.resize(N);
      counts_collinear.resize(N);
      ratio_sl.resize(N);
      ss.str(std::string()); ss << string_out << tag << ".dat"; fout.open(ss.str(), ios::out);
      fout << "#Timestep | Particle index  | Particle type | q_oct | <r_short/r_long>_collinear. # cutoffs = ";
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
      ntype rijSq,rikSq, sum_my_oct, sum_ratio_sl, counts_ratio_sl;
      ntype cos_k; //, r1r2[3]={0,0,0};
      vec rij,rik;
      int i,jj,kk, j,k, num_neigh; //, m;
      if(debug) cout << "\n*** "<<myName<<" computation STARTED ***\n";
      for(i=0;i<N;i++)
      {
        sum_my_oct = sum_ratio_sl = 0.0;
        counts_collinear[i]=0.0;
//	m=0;
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
              if(!continuous){
                sum_my_oct+=SQUARE(cos_k+1) / 3.0;
              }
              sum_ratio_sl += (rijSq<rikSq ? sqrt(rijSq/rikSq) : sqrt(rikSq/rijSq) );
              counts_collinear[i] += 1;
	  //    r1r2[m++] = (rijSq<rikSq ? sqrt(rijSq/rikSq) : sqrt(rikSq/rijSq) );
            } else { // else closer to 90°
              if(!continuous){
                sum_my_oct+=SQUARE(cos_k) / 12.0;
              }
            }
            // NEW CONTINUOUS DEFINITION
            if(continuous){
              sum_my_oct += (1-SQUARE( 2*SQUARE(cos_k) -1 ))/15.0;
            }

          }
        }

        my_q_oct[i] = normalization * sum_my_oct;
        if(counts_collinear[i]==0){
          ratio_sl[i] = 0.0;
          /*
          cerr<<"ERROR: particle "<<i<<" at timestep "<<timestep<<" has no bonds > "<<rad2deg*acos(my_cos_th)<<" degrees\n";
          cerr<<"List of neighbors:";
          for(jj=0;jj<6;jj++)
          {
            j = ps[i].neigh_list[u][jj];
            cerr<<" "<<j;
          }
          cerr<<endl;
          cerr<<"List of distances:";
          for(jj=0;jj<6;jj++)
          {
            rij = ps[i].rij_list[u][jj];
            rijSq = ps[i].rijSq_list[u][jj];
            cerr<<" "<<sqrt(rijSq);
          }
          cerr<<endl;
          cerr<<"List of angles:";
          for(jj=0;jj<6;jj++)
          {
            rij = ps[i].rij_list[u][jj];
            rijSq = ps[i].rijSq_list[u][jj];
            for(kk=jj+1;kk<6;kk++)
            {
              rik = ps[i].rij_list[u][kk];
              rikSq = ps[i].rijSq_list[u][kk];
              // angle j-i-k (polar)
              cos_k = (rij*rik) / sqrt(rijSq*rikSq);
              cerr<<" "<<rad2deg*acos(cos_k);
            }
          }
          cerr<<endl;
          exit(1);
          */
        } else {
          ratio_sl[i] = sum_ratio_sl/counts_collinear[i];
        }

      }

      ss.str(std::string()); ss << string_out << tag << ".dat"; fout.open(ss.str(), ios::app);
      for(i=0;i<N;i++)
      {
        fout << timestep <<" "<< i <<" "<< ps[i].label << " " << my_q_oct[i]<< " "<<ratio_sl[i]<< endl;
      }
      fout.close();
      if(debug) cout << "*** "<<myName<<" computation DONE ***\n";
    }

};
#endif
