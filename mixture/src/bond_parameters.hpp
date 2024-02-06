#ifndef _BOND_PARAMETERS_H_
#define _BOND_PARAMETERS_H_

#include "neighbour_and_bond_list.hpp"
using namespace std;
//------ Bond Orientational parameters: BOO (q_lm) and BOC (q_lm^dot) --------//
// Bond orientational order parameter:
// q_lm(i) is the local average of spherical harmonics:
//   q_lm(i) = <Y_lm(r_ij)>_{neighbours j} = Sum_j { Y_lm(r_ij) } / Sum_j
//   Note: it's a complex vector with index m=-l,...,l.
//   Note: here we define the neighbours by a cutoff radius rcut1 (1st sphere).
//   Note: here we use a smoothed version where the Sum_j is weighted by a
//     smooth cutoff function f(r_ij/rcut1).
//   This is saved into ${string_bondorient_out}.dat
// q_l(i) is a local rotationally invariant descriptor of the environment (Steinhardt parameter)
//   q_l(i) = { 4*pi/(2*l+1) * Sum_m |q_lm(i)|^2 }^(1/2)
//          = {4*pi/(2*l+1)}^(1/2) * Norm(q_lm(i))
//   Note: it is proportional to the norm of the vector q_lm(i) with index m.
//   This is saved into ${string_bondorient_out}.ave
//
// Bond order correlation parameter:
// C_l(i,j) is the cosine similarity of a bond <i,j> between particles i,j,
// in terms of the BOO. A bond is defined by being in the 1st sphere.
//   C_l(i,j) = DotProduct( q_lm(i), q_lm(j) ) / [ Norm( q_lm(i) ) * Norm( q_lm(j) ) ]
//   Note: if i,j belong to the same crystalline cluster: C~1 ; otherwise C~0.
//   Problem: C_l is complex... I take the real part. Should I take the norm? No, real part is better
//   Problem: what happens when Norm(q_l(i))==0.0 ? I define C_l=0, but this is arbitrary.
// q_l^dot(i) is the local average of C_l(i,j) over neighbours j.
//   Note: here we use a cutoff radius rcut2 (2nd sphere) which can be > rcut1.
//   Note: here we weight Sum_j by a smooth cutoff function f(r_ij/rcut2).
//   This is saved into ${string_bondcorr_out}.ave
//
// One then defines that particle i is crystalline if q_l^dot(i) > threshold (usually 0.55 or 0.65).
// This is saved into the string_nxtal_out files.
//
// Otherwise, you can observe the average of q_l^dot(i) around any particle i,
// within a large cutoff radius rcut3 (3rd sphere).
//----------------------------------------------------------------------------//
#define PRINT_CUTOFFS_HEADER(Nsph,nbl,fout) {for(auto  u=0;u<Nsph;u++){ for(int t=0;t<nbl->nTypePairs;t++) {fout << nbl->rcut[u][t]<<" ";} fout << "; "; }; fout << endl;}

template <class ntype, class ptype>
class Bond_Parameters
{
  using vec=myvec<ntype,3>;
  using mat=mymatrix<ntype,3,3>;
  private:
    Neigh_and_Bond_list<ntype,ptype> *nb_list;
    int l, l_deg; // angular momentum and its degeneracy
    bool l_is_odd;
    vecflex<ntype> ql, ql_ave, Cl_ij, ql_dot, Ql_dot;
    vector< vecflex< complex<ntype> > > qlm, qlm_ave; // a collection of l_deg vectors of local average qlm=<Ylm> with a cutoff function
    ntype qldot_th;
    string string_bondcorr_out, string_bondorient_out, string_nxtal_out, myName, tag;
    fstream fout;
    stringstream ss;bool debug, verbose;

  public:
    Bond_Parameters(){
      myName = "BOND PARAMETERS";
    }
    virtual ~Bond_Parameters(){}

    void init(Neigh_and_Bond_list<ntype,ptype> *nb_list_, int l_, ntype qldot_th_,
      string string_boo_out_, string string_boc_out_, string string_nxtal_out_,
      string tag_, bool debug_, bool verbose_)
    {
      string_bondcorr_out = string_boc_out_;
      string_bondorient_out = string_boo_out_;
      string_nxtal_out = string_nxtal_out_;
      tag = tag_;
      debug = debug_;
      verbose = verbose_;
      nb_list = nb_list_;
      l = l_;
      qldot_th = qldot_th_;
      l_is_odd = (l%2!=0);
      l_deg = 2*l+1;
      const int N=nb_list->N;
      const int Nsphere=nb_list->Nsphere;
      const int nTypes=nb_list->nTypes;
      const int nTypePairs=nb_list->nTypePairs;
      qlm.resize(l_deg);
      for(int a=0;a<l_deg;a++) qlm[a].resize(N);
      ql.resize(N);
      qlm_ave.resize(l_deg);
      for(int a=0;a<l_deg;a++) qlm_ave[a].resize(N);
      ql_ave.resize(N);
      ql_dot.resize(N);
      Ql_dot.resize(N);

      if(verbose)
      {
        ss.str(std::string()); ss << string_bondorient_out << ".l" << l << tag << ".dat"; fout.open(ss.str(), ios::out);
        fout << "#Timestep, Particle index, Bond order orientation parameter q_l. # l = " << l << ", cutoffs = ";
        PRINT_CUTOFFS_HEADER(Nsphere,nb_list,fout);
        fout.close();
        ss.str(std::string()); ss << string_bondorient_out << ".l" << l << tag << ".ave"; fout.open(ss.str(), ios::out);
        fout << "#Timestep, <q_l>, fluctuations. # l = " << l << ", cutoffs = ";
        PRINT_CUTOFFS_HEADER(Nsphere,nb_list,fout);
        fout.close();

        ss.str(std::string()); ss << string_bondorient_out << "_ave.l" << l << tag << ".dat"; fout.open(ss.str(), ios::out);
        fout << "#Timestep, Particle index, bar{q_l}: Local average of q_l within 1st-neigh-sphere. # l = " << l << ", cutoffs = ";
        PRINT_CUTOFFS_HEADER(Nsphere,nb_list,fout);
        fout.close();
        ss.str(std::string()); ss << string_bondorient_out << "_ave.l" << l << tag << ".ave"; fout.open(ss.str(), ios::out);
        fout << "#Timestep, <bar{q_l}>, fluctuations. # l = " << l << ", cutoffs = ";
        PRINT_CUTOFFS_HEADER(Nsphere,nb_list,fout);
        fout.close();
      }

      ss.str(std::string()); ss << string_bondcorr_out << ".l" << l << tag << ".dat"; fout.open(ss.str(), ios::out);
      fout << "#Timestep, Particle index, Bond order correlation parameter q_l_dot. # l = " << l << ", cutoffs = ";
      PRINT_CUTOFFS_HEADER(Nsphere,nb_list,fout);
      fout.close();
      ss.str(std::string()); ss << string_bondcorr_out << ".l" << l << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "#Timestep, <q_l_dot>, fluctuations. # l = " << l << ", cutoffs = ";
      PRINT_CUTOFFS_HEADER(Nsphere,nb_list,fout);
      fout.close();

      ss.str(std::string()); ss << string_nxtal_out << ".l" << l << tag << ".dat"; fout.open(ss.str(), ios::out);
      fout << "#Timestep, Nc = number of crystalline particles # l = " << l << ", cutoffs = ";
      PRINT_CUTOFFS_HEADER(Nsphere,nb_list,fout);
      fout.close();
      ss.str(std::string()); ss << string_nxtal_out << ".l" << l << tag << ".indexes"; fout.open(ss.str(), ios::out);
      fout << "#Timestep | Indexes of crystalline particles (if any) # l = " << l << ", cutoffs = ";
      PRINT_CUTOFFS_HEADER(Nsphere,nb_list,fout);
      fout.close();

      if(verbose)
      {
        ss.str(std::string()); ss << string_bondcorr_out << ".l" << l << tag << ".local_ave"; fout.open(ss.str(), ios::out);
        fout << "#Timestep, Particle index, atom-centered <q_l_dot(i)>. # l = " << l << ", cutoffs = ";
        PRINT_CUTOFFS_HEADER(Nsphere,nb_list,fout);;
        fout.close();
      }
    }


    void compute(int timestep, vector<ptype> ps)
    {
      if(debug) nb_list->print_bond_summary(ps);
      const int N=nb_list->N;
      const int Nsphere=nb_list->Nsphere;
      const int nTypes=nb_list->nTypes;
      const int nTypePairs=nb_list->nTypePairs;
      ntype fval, rijSq, Q,Q2, ql_factor = sqrt(4.0*M_PI/(ntype)l_deg), re,im;
      complex<ntype> Yval;
      vec rij;
      int i,j,k,a,b,m,u, t;
      if(debug) cout << "\n*** "<<myName<<" computation STARTED ***\n";
      //---- Reset counters ----//
      for(i=0;i<N;i++){
          ql.set(i, 0.0);
          ql_ave.set(i, 0.0);
          ql_dot.set(i, 0.0);
          Ql_dot.set(i, 0.0);
          for(a=0;a<l_deg; a++){
            qlm[a].set(i, 0.0);
            qlm_ave[a].set(i, 0.0);
          }
          for(u=0;u<Nsphere;u++){
            nb_list->neigh_anytype[u][i] = 0.;
            for(t=0;t<nTypePairs;t++) nb_list->neigh[u][t][i] = 0.;
          }
      }
      if(debug) cout << " * Reset counters DONE\n";

      //---- Compute <Y_lm>'s (on 1st sphere) and neighs (on all spheres) ----//
      for(i=0;i<N;i++){
        for(u=0;u<Nsphere;u++){
          for(k=0;k<ps[i].neigh_list[u].size();k++){
            j = ps[i].neigh_list[u][k];
            if(i>j) continue; // avoid double counting!
            t = nb_list->types2int(ps[i].label, ps[j].label);
            rij = ps[i].rij_list[u][k];
            rijSq = ps[i].rijSq_list[u][k];
            // Ider uses a step function for u==0?
            // if(u==0) fval = ( rijSq <= nb_list->rcutSq[u][t] ? 1.0 : 0.0 );
            fval = nb_list->fcut( rijSq/nb_list->rcutSq[u][t] );
            nb_list->neigh[u][t][i] += fval;
            nb_list->neigh_anytype[u][i] += fval;
            t = nb_list->types2int(ps[j].label, ps[i].label);
            nb_list->neigh[u][t][j] += fval;
            nb_list->neigh_anytype[u][j] += fval;
            if(u!=0) continue; // Y(l,m) are computed on first sphere neighbours only
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
      if(debug) cout << " * Compute qlm(i) (to be normalized) and neigh(i) DONE\n";

      // Normalize qlm
      u=0;
      for(i=0;i<N;i++){
        if(nb_list->neigh_anytype[u][i]>0){ // avoid divergence
          for(a=0;a<l_deg; a++) qlm[a][i] /= nb_list->neigh_anytype[u][i];
        }
      }
      if(debug) cout << " * Normalization of qlm(i) DONE\n";

      // compute a locally averaged version of qlm(i) within the 1st sphere
      u=0;
      for(i=0;i<N;i++){
        for(k=0;k<ps[i].neigh_list[u].size();k++){
          j = ps[i].neigh_list[u][k];
          t = nb_list->types2int(ps[i].label, ps[j].label);
          rij = ps[i].rij_list[u][k];
          rijSq = ps[i].rijSq_list[u][k];
          // Ider uses a step function for u==0?
          // if(u==0) fval = ( rijSq <= nb_list->rcutSq[u][t] ? 1.0 : 0.0 );
          fval = nb_list->fcut( rijSq/nb_list->rcutSq[u][t] );
          for(a=0;a<l_deg; a++){
            qlm_ave[a][i] += fval*qlm[a][j]/(1.0 + nb_list->neigh_anytype[u][i]);
          }
        }
        for(a=0;a<l_deg; a++){
          qlm_ave[a][i] += qlm[a][i]/(1.0 + nb_list->neigh_anytype[u][i]); // sum particle i as well
        }
        // qlm_ave(i) = (local average of qlm(i) ) COMPLETED
      }
      if(debug) cout << " * Compute ql(m,i), ql(m,i)_ave and neigh(i) DONE\n";

      //---- Compute ql~qlm(i)*qlm(i) ----//
      Q = 0.0; //  <ql>
      Q2 = 0.0; // <ql^2>
      if(verbose) {
        ss.str(std::string()); ss << string_bondorient_out << ".l" << l << tag << ".dat"; fout.open(ss.str(), ios::app);
      }
      for(i=0;i<N;i++)
      {
        for(a=0;a<l_deg; a++) {
          re=real(qlm[a][i]);
          im=imag(qlm[a][i]);
          ql[i] += re*re + im*im;
        }
        ql[i] = ql_factor * sqrt(ql[i]);
        if(verbose) fout << timestep << " " << i << " " << ql[i] << endl;
        Q += ql[i] / N;
        Q2 += ql[i]*ql[i] / N;
      }
      if(verbose) fout.close();
      if(verbose)
      {
        ss.str(std::string()); ss << string_bondorient_out << ".l" << l << tag << ".ave"; fout.open(ss.str(), ios::app);
        fout << timestep << " " << Q << " " << sqrt((Q2-Q*Q)/N) << endl;
        fout.close();
      }
      if(debug) cout << " * Compute ql(i) (BOO) DONE\n";

      //---- Compute ql_ave~qlm_ave(i)*qlm_ave(i) ----//
      Q = 0.0; //  <ql>
      Q2 = 0.0; // <ql^2>
      if(verbose) {
        ss.str(std::string()); ss << string_bondorient_out << "_ave.l" << l << tag << ".dat"; fout.open(ss.str(), ios::app);
      }
      for(i=0;i<N;i++)
      {
        for(a=0;a<l_deg; a++) {
          re=real(qlm_ave[a][i]);
          im=imag(qlm_ave[a][i]);
          ql_ave[i] += re*re + im*im;
        }
        ql_ave[i] = ql_factor * sqrt(ql_ave[i]);
        if(verbose) fout << timestep << " " << i << " " << ql_ave[i] << endl;
        Q += ql_ave[i] / N;
        Q2 += ql_ave[i]*ql_ave[i] / N;
      }
      if(verbose) fout.close();
      if(verbose)
      {
        ss.str(std::string()); ss << string_bondorient_out << "_ave.l" << l << tag << ".ave"; fout.open(ss.str(), ios::app);
        fout << timestep << " " << Q << " " << sqrt((Q2-Q*Q)/N) << endl;
        fout.close();
      }
      if(debug) cout << " * Compute ql(i)_ave (averaged BOO) DONE\n";

      //---- Compute Cl_ij~qlm(i)*qlm(j) (I take the real part) ----//
      u=1; // C_ij runs over 2nd-sphere bonds.
      Cl_ij.resize( nb_list->bond_list[u].size() );
      for(b=0;b<nb_list->bond_list[u].size();b++){
        Cl_ij.set(b, 0.0);
        i = nb_list->int2i( nb_list->bond_list[u][b], N );
        j = nb_list->int2j( nb_list->bond_list[u][b], N );
        // t = nb_list->types2int(ps[i].label, ps[j].label); // useless...
        for(a=0;a<l_deg; a++) {
          Cl_ij[b] += real( qlm[a][i]*conj(qlm[a][j]) );
          // abs or real or complex Cl_ij?
        }
        if( ql[i]==0.0 || ql[j]==0.0) Cl_ij[b]=0.0; // safety condition for divergence (NOT JUSTIFIED!)
        else Cl_ij[b] *= (ql_factor*ql_factor)/(ql[i]*ql[j]); // normalization factors
        // Cl_ij[b] done.
      }
      if(debug) cout << " * Compute C_l(i,j) DONE\n";

      //---- Compute ql_dot(i) (BOC) ----//
      u=1; // q_l^dot is the average of C_l(ij) over the 2nd sphere
      for(i=0;i<N;i++){
        for(k=0;k<ps[i].neigh_list[u].size();k++){
          j = ps[i].neigh_list[u][k];
          if(i>j) continue; // avoid double counting
          rijSq = ps[i].rijSq_list[u][a]; // and use it to recover the radius
          fval = nb_list->fcut( rijSq/nb_list->rcutSq[u][t] ); // smooth cutoff
          b = nb_list->get_bond_index(i,j,N,u); // bond index in the bond list
          ql_dot[i] += fval*Cl_ij[b] / nb_list->neigh_anytype[u][i];
          ql_dot[j] += fval*Cl_ij[b] / nb_list->neigh_anytype[u][j];  // take into account i>j
        }
      }
      if(debug) cout << " * Compute ql_dot(i) (BOC) DONE\n";

      //---- Compute global average of ql_dot(i) ----//
      ss.str(std::string()); ss << string_bondcorr_out << ".l" << l << tag << ".dat"; fout.open(ss.str(), ios::app);
      Q = 0.0; //  <ql_dot>
      Q2 = 0.0; // <ql_dot^2>
      for(i=0;i<N;i++){
        fout << timestep << " " << i << " " << ql_dot[i] << endl;
        Q += ql_dot[i] / N;
        Q2 += ql_dot[i]*ql_dot[i] / N;
      }
      fout.close();
      ss.str(std::string()); ss << string_bondcorr_out << ".l" << l << tag << ".ave"; fout.open(ss.str(), ios::app);
      fout << timestep << " " << Q << " " << sqrt((Q2-Q*Q)/N) << endl;
      fout.close();
      if(debug) cout << " * Compute global average of BOC DONE\n";

      if(verbose) {
        //---- Compute a Locally Confined version of the global BOC ~ average <ql_dot(i)> around i ----//
        u=2; // 3rd sphere of neighbours
        ss.str(std::string()); ss << string_bondcorr_out << ".l" << l << tag << ".local_ave"; fout.open(ss.str(), ios::app);
        for(i=0;i<N;i++){
          for(k=0;k<ps[i].neigh_list[u].size();k++){
            j = ps[i].neigh_list[u][k];
            t = nb_list->types2int(ps[i].label, ps[j].label);
            rijSq = ps[i].rijSq_list[u][k];
            fval = nb_list->fcut( rijSq/nb_list->rcutSq[u][t] );
            Ql_dot[i] += fval*ql_dot[j] / nb_list->neigh_anytype[u][j];
          }
          fout << timestep << " " << i << " " << Ql_dot[i] << endl;
        }
        fout.close();
        if(debug) cout << " * Compute local average of BOC DONE\n";
      }

      //----- Compute Number of crystalline particles -----//
      ss.str(std::string()); ss << string_nxtal_out << ".l" << l << tag << ".indexes"; fout.open(ss.str(), ios::app);
      m=0; // number of crystalline particles
      fout << timestep;
      for(i=0;i<N;i++){
        if(ql_dot[i] > qldot_th){
          fout << " " << i; // indexes of crystalline particles
          m++;
        }
      }
      fout<<endl;
      fout.close();
      ss.str(std::string()); ss << string_nxtal_out << ".l" << l << tag << ".dat"; fout.open(ss.str(), ios::app);
      fout << timestep << " " << m << endl;
      fout.close();
      if(debug) cout << "*** "<<myName<<" computation for timestep " << timestep << " ENDED ***\n\n";
  }

};
#endif
