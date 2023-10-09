#ifndef _BOND_PARAMETERS_H_
#define _BOND_PARAMETERS_H_

#include "neighbour_and_bond_list.hpp"
using namespace std;
//---------------------- Bond Orientational parameters: BOO (q_lm) and BOC (q_lm^dot) -------------------------------//
template <class ntype, class ptype>
class Bond_Parameters
{
  using vec=myvec<ntype,3>;
  using mat=mymatrix<ntype,3,3>;
  private:
    Neigh_and_Bond_list<ntype,ptype> *nb_list;
    int l, l_deg; // angular momentum and its degeneracy
    bool l_is_odd;
    vecflex<ntype> ql, Cl_ij, ql_dot, Ql_dot;
    vector< vecflex< complex<ntype> > > qlm; // a collection of l_deg vectors of local average qlm=<Ylm> with a cutoff function
    ntype qldot_th;
    string string_bondcorr_out, string_bondorient_out, string_nxtal_out, myName, tag;
    fstream fout;
    stringstream ss;
  public:
    Bond_Parameters(){
      myName = "BOND PARAMETERS";
    }
    virtual ~Bond_Parameters(){}

    void init(Neigh_and_Bond_list<ntype,ptype> *nb_list_, int l_, ntype qldot_th_, string string_boo_out_, string string_boc_out_, string string_nxtal_out_, string tag_)
    {
      string_bondcorr_out = string_boc_out_;
      string_bondorient_out = string_boo_out_;
      string_nxtal_out = string_nxtal_out_;
      tag = tag_;
      nb_list = nb_list_;
      l = l_;
      qldot_th = qldot_th_;
      l_is_odd = (l%2!=0);
      l_deg = 2*l+1;
      qlm.resize(l_deg);
      const int N=nb_list->N;
      const int Nshell=nb_list->Nshell;
      for(int a=0;a<l_deg;a++) qlm[a].resize(N);
      ql.resize(N);
      ql_dot.resize(N);
      Ql_dot.resize(N);
      ss.str(std::string()); ss << string_bondorient_out << ".l" << l << tag << ".dat"; fout.open(ss.str(), ios::out);
      fout << "#Timestep, Particle index, Bond order orientation parameter q_l. # l = " << l << ", cutoffs = ";
      for(auto  i=0;i<Nshell;i++) fout << nb_list->rcut[i]<<", ";
      fout << endl;
      fout.close();
      ss.str(std::string()); ss << string_bondorient_out << ".l" << l << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "#Timestep, <q_l>, fluctuations. # l = " << l << ", cutoffs = ";
      for(auto  i=0;i<Nshell;i++) fout << nb_list->rcut[i]<<", ";
      fout << endl;
      fout.close();
      ss.str(std::string()); ss << string_bondcorr_out << ".l" << l << tag << ".dat"; fout.open(ss.str(), ios::out);
      fout << "#Timestep, Particle index, Bond order correlation parameter q_l_dot. # l = " << l << ", cutoffs = ";
      for(auto  i=0;i<Nshell;i++) fout << nb_list->rcut[i]<<", ";
      fout << endl;
      fout.close();
      ss.str(std::string()); ss << string_nxtal_out << ".l" << l << tag << ".dat"; fout.open(ss.str(), ios::out);
      fout << "#Timestep, Nc = number of crystalline particles # l = " << l << ", cutoffs = ";
      for(auto  i=0;i<Nshell;i++) fout << nb_list->rcut[i]<<", ";
      fout << endl;
      fout.close();
      ss.str(std::string()); ss << string_bondcorr_out << ".l" << l << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "#Timestep, <q_l_dot>, fluctuations. # l = " << l << ", cutoffs = ";
      for(auto  i=0;i<Nshell;i++) fout << nb_list->rcut[i]<<", ";
      fout << endl;
      fout.close();
      ss.str(std::string()); ss << string_bondcorr_out << ".l" << l << tag << ".xyz"; fout.open(ss.str(), ios::out);
      fout.close();
      ss.str(std::string()); ss << string_bondcorr_out << ".l" << l << tag << ".local.ave"; fout.open(ss.str(), ios::out);
      fout << "#Timestep, Particle index, local <q_l_dot(i)>. # l = " << l << ", cutoffs = ";
      for(auto  i=0;i<Nshell;i++) fout << nb_list->rcut[i]<<", ";
      fout << endl;
      fout.close();
    }

    void compute(int timestep, vector<ptype> ps, bool debug)
    {
      if(debug) nb_list->print_bond_summary(ps);
      const int N=nb_list->N;
      const int Nshell=nb_list->Nshell;
      ntype fval, rijSq, Q,Q2, ql_factor = sqrt(4.0*M_PI/(ntype)l_deg), re,im;
      complex<ntype> Yval;
      vec rij;
      int i,j,k,a,m,u;
      if(debug) cout << "\n*** "<<myName<<" computation STARTED ***\n";
      //---- Reset counters ----//
      for(i=0;i<N;i++){
          ql.set(i, 0.0);
          ql_dot.set(i, 0.0);
          Ql_dot.set(i, 0.0);
          for(a=0;a<l_deg; a++) qlm[a].set(i, 0.0);
      }
      if(debug) cout << " * Reset counters DONE\n";
      //---- Compute Y_lm's and neighs ----//
      for(i=0;i<N;i++){
        for(u=0;u<Nshell;u++){
          for(k=0;k<ps[i].neigh_list[u].size();k++){
            j = ps[i].neigh_list[u][k];
            if(j>i) continue; // avoid double counting!
            rij = ps[i].rij_list[u][k];
            rijSq = ps[i].rijSq_list[u][k];
  //          fval = fcut( rijSq/cutoffSq[u], p1half, p2half );
            if(u==0){ // Ider uses a step function for u==0?
              fval = ( rijSq <= nb_list->rcut[u] ? 1.0 : 0.0 );
            } else {
              fval = nb_list->fcut( rijSq/nb_list->rcutSq[u] );
            }
            nb_list->neigh[u][i] += fval;
            nb_list->neigh[u][j] += fval;
            if(u!=0) continue; // Y(l,m) are computed on first shell neighbours only
            for(a=0;a<l_deg; a++){
                m=a-l; // -l <= m <= l
                Yval = Y(l,m, rij, sqrt(rijSq) );
                qlm[a][i] += fval*Yval;
                if(l_is_odd) Yval=-Yval; // Y(-r) = (-1)^l * Y(r)
                qlm[a][j] += fval*Yval;
            }
          }
        }
      }
      if(debug) cout << " * Compute ql(m,i) and neigh(i) DONE (but ql(m,i) has yet to be normalized)\n";
      //---- Compute ql~qlm(i)*qlm(i) ----//
      Q = 0.0; //  <ql>
      Q2 = 0.0; // <ql^2>
      ss.str(std::string()); ss << string_bondorient_out << ".l" << l << tag << ".dat"; fout.open(ss.str(), ios::app);
      for(i=0;i<N;i++){
          for(a=0;a<l_deg; a++) {
             if(nb_list->neigh[0][i]>0) qlm[a][i] /= nb_list->neigh[0][i]; // qlm(i) = <Ylm>(i) completed
             re=real(qlm[a][i]);
             im=imag(qlm[a][i]);
             ql[i] += re*re + im*im;
          }
          ql[i] = ql_factor * sqrt(ql[i]);
          fout << timestep << " " << i << " " << ql[i] << endl;
          Q += ql[i] / N;
          Q2 += ql[i]*ql[i] / N;
      }
      fout.close();
      ss.str(std::string()); ss << string_bondorient_out << ".l" << l << tag << ".ave"; fout.open(ss.str(), ios::app);
      fout << timestep << " " << Q << " " << sqrt((Q2-Q*Q)/N) << endl;
      fout.close();
      if(debug) cout << " * Compute ql(i) (BOO) DONE\n";
      //---- Compute Cl_ij~qlm(i)*qlm(j) (I take the real part) and ql_dot(i) (BOC) ----//
      Cl_ij.resize( nb_list->bond_list[1].size() );
      for(k=0;k<nb_list->bond_list[1].size();k++){ // Cij and BOC are computed on 2nd shell neighbours
        Cl_ij.set(k, 0.0);
        i = nb_list->int2i( nb_list->bond_list[1][k], N );
        j = nb_list->int2j( nb_list->bond_list[1][k], N );
        for(a=0;a<l_deg; a++) {
          Cl_ij[k] += ( real(qlm[a][i])*real(qlm[a][j]) + imag(qlm[a][i])*imag(qlm[a][j]) );
        }
        if( ql[i]==0.0 || ql[j]==0.0) Cl_ij[k]=0.0; // safety condition (NOT JUSTIFIED!)
        else Cl_ij[k] /= (ql[i]*ql[j]); // Cl_ij[k] done
        a = indexOf<int>( ps[i].neigh_list[1], j ); // find index of j in i's neighbour list
        if(a<0) {
          cout << "[ERROR: could not find j="<<j<<" in the 2nd-shell-neighbour-list of i="<<i<<"]\n";
          cout << "  2nd-shell list for particle i="<<i<<" contains "<<ps[i].neigh_list[1].size()<<" neighbours:\n   ";
          for(int ii=0;ii<ps[i].neigh_list[1].size();ii++) cout<<ps[i].neigh_list[1][ii]<<" ";
          cout << endl;
          cout << "  2nd-shell list for particle j="<<j<<" contains "<<ps[i].neigh_list[1].size()<<" neighbours:\n   ";
          for(int ii=0;ii<ps[j].neigh_list[1].size();ii++) cout<<ps[j].neigh_list[1][ii]<<" ";
          cout << endl;
          exit(1);
        }
        rijSq = ps[i].rijSq_list[1][a]; // and use it to recover the radius
        fval = nb_list->fcut( rijSq/nb_list->rcutSq[1] );
        ql_dot[i] += fval*Cl_ij[k] / nb_list->neigh[1][i];
        ql_dot[j] += fval*Cl_ij[k] / nb_list->neigh[1][j];
      }
      if(debug) cout << " * Compute C_l(i,j) and ql_dot(i) (BOC) DONE\n";
      //---- Compute global average of ql_dot(i) ----//
      ss.str(std::string()); ss << string_bondcorr_out << ".l" << l << tag << ".dat"; fout.open(ss.str(), ios::app);
      Q = 0.0; //  <ql_dot>
      Q2 = 0.0; // <ql_dot^2>
      for(i=0;i<N;i++){
  /*      for(k=0;k<ps[i].neigh_list[1].size();k++){
          j = ps[i].neigh_list[1][k];
          rijSq = ps[i].rijSq_list[1][k];
          fval = fcut( rijSq/cutoffSq[1], p1half, p2half );
          ql_dot[i] += fval*Cl_ij[ a++ ] / neigh[1][i];
        }*/
        fout << timestep << " " << i << " " << ql_dot[i] << endl;
        Q += ql_dot[i] / N;
        Q2 += ql_dot[i]*ql_dot[i] / N;
      }
      fout.close();
      ss.str(std::string()); ss << string_bondcorr_out << ".l" << l << tag << ".ave"; fout.open(ss.str(), ios::app);
      fout << timestep << " " << Q << " " << sqrt((Q2-Q*Q)/N) << endl;
      fout.close();
      if(debug) cout << " * Compute global average of BOC DONE\n";
      //---- Compute a Locally Confined version of the global BOC ~ average <ql_dot(i)> around i ----//
      ss.str(std::string()); ss << string_bondcorr_out << ".l" << l << tag << ".local.ave"; fout.open(ss.str(), ios::app);
      for(i=0;i<N;i++){
        for(k=0;k<ps[i].neigh_list[2].size();k++){
          j = ps[i].neigh_list[2][k];
          rijSq = ps[i].rijSq_list[2][k];
          fval = nb_list->fcut( rijSq/nb_list->rcutSq[2] );
          Ql_dot[i] += fval*ql_dot[j] / nb_list->neigh[2][i];
        }
        fout << timestep << " " << i << " " << Ql_dot[i] << endl;
      }
      fout.close();
      if(debug) cout << " * Compute local average of BOC DONE\n";
      ss.str(std::string()); ss << string_bondcorr_out << ".l" << l << tag << ".xyz"; fout.open(ss.str(), ios::app);
      fout << N << endl << "-1: crystal, -2: noncrystal. Timestep = " << timestep << endl;
      m=0; // number of crystalline particles
      for(i=0;i<N;i++){
        if(ql_dot[i] > qldot_th){
          ps[i].label =  -1;
          m++;
        } else{
          ps[i].label =  -2;
        }
        ps[i].write_xyz(fout);
      }
      fout.close();
      ss.str(std::string()); ss << string_nxtal_out << ".l" << l << tag << ".dat"; fout.open(ss.str(), ios::app);
      fout << timestep << " " << m << endl;
      fout.close();
      if(debug) cout << "*** "<<myName<<" computation for timestep " << timestep << " ENDED ***\n\n";
  }

};
#endif
