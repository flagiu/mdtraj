using namespace std;
//---------------------- Bond Orientational parameters: BOO (q_lm) and BOC (q_lm^dot) -------------------------------//
template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
init_bondorient()
    {
      /* // Useless: already done in init_neigh
      for(auto i=0;i<Nshells;i++){
        if(neigh[i].size()!=nTypePairs) neigh[i].resize(nTypePairs);
        for(auto t=0;t<nTypePairs;t++)
          if(neigh[i][t].length()!=N) neigh[i][t].resize(N);
      } */
      qlm.resize(l_deg);
      for(auto a=0;a<l_deg;a++) qlm[a].resize(N);
      ql.resize(N);
      ql_dot.resize(N);
      Ql_dot.resize(N);
      ss.str(std::string()); ss << s_bondorient << ".l" << l << tag << ".dat"; fout.open(ss.str(), ios::out);
      fout << "#Timestep, Particle index, Bond order orientation parameter q_l. # l = " << l << ", cutoffs = ";
      for(auto  i=0;i<Nshells;i++) fout << cutoff[i]<<", ";
      fout << endl;
      fout.close();
      ss.str(std::string()); ss << s_bondorient << ".l" << l << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "#Timestep, <q_l>, fluctuations. # l = " << l << ", cutoffs = ";
      for(auto  i=0;i<Nshells;i++) fout << cutoff[i]<<", ";
      fout << endl;
      fout.close();
      ss.str(std::string()); ss << s_bondcorr << ".l" << l << tag << ".dat"; fout.open(ss.str(), ios::out);
      fout << "#Timestep, Particle index, Bond order correlation parameter q_l_dot. # l = " << l << ", cutoffs = ";
      for(auto  i=0;i<Nshells;i++) fout << cutoff[i]<<", ";
      fout << endl;
      fout.close();
      ss.str(std::string()); ss << s_nxtal << ".l" << l << tag << ".dat"; fout.open(ss.str(), ios::out);
      fout << "#Timestep, Nc = number of crystalline particles # l = " << l << ", cutoffs = ";
      for(auto  i=0;i<Nshells;i++) fout << cutoff[i]<<", ";
      fout << endl;
      fout.close();
      ss.str(std::string()); ss << s_bondcorr << ".l" << l << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "#Timestep, <q_l_dot>, fluctuations. # l = " << l << ", cutoffs = ";
      for(auto  i=0;i<Nshells;i++) fout << cutoff[i]<<", ";
      fout << endl;
      fout.close();
      ss.str(std::string()); ss << s_bondcorr << ".l" << l << tag << ".xyz"; fout.open(ss.str(), ios::out);
      fout.close();
      ss.str(std::string()); ss << s_bondcorr << ".l" << l << tag << ".local.ave"; fout.open(ss.str(), ios::out);
      fout << "#Timestep, Particle index, local <q_l_dot(i)>. # l = " << l << ", cutoffs = ";
      for(auto  i=0;i<Nshells;i++) fout << cutoff[i]<<", ";
      fout << endl;
      fout.close();
}

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
compute_bondorient() {
    ntype fval, rijSq, Q,Q2, ql_factor = sqrt(4*M_PI/l_deg), re,im, tot_neigh;
    complex<ntype> Yval;
    vec rij;
    int i,j,k,a,m,u, t;
    if(debug) cout << "\n*** BOO,BOC computation STARTED ***\n";
    //---- Reset counters ----//
    for(i=0;i<N;i++){
        ql.set(i, 0.0);
        ql_dot.set(i, 0.0);
        Ql_dot.set(i, 0.0);
        for(a=0;a<l_deg; a++) qlm[a].set(i, 0.0);
        for(u=0;u<Nshells;u++)
          for(t=0;t<nTypePairs;t++)
            neigh[u][t][i] = 0.;
    }
    if(debug) cout << " * Reset counters DONE\n";
    //---- Compute Y_lm's and neighs ----//
    for(i=0;i<N;i++){
      for(u=0;u<Nshells;u++){
        for(k=0;k<ps[i].neigh_list[u].size();k++){
          j = ps[i].neigh_list[u][k];
          if(j>i) continue; // avoid double counting!
          rij = ps[i].rij_list[u][k];
          rijSq = ps[i].rijSq_list[u][k];
//          fval = fcut( rijSq/cutoffSq[u], p1half, p2half );
          if(u==0){ //Ider uses a step function for u==0 ?
            fval = ( rijSq <= cutoffSq[u] ? 1.0 : 0.0 );
          } else {
            fval = fcut( rijSq/cutoffSq[u], p1half, p2half );
          }
          t = types2int(ps[i].label, ps[j].label, nTypes);
          neigh[u][t][i] += fval;
          t = types2int(ps[j].label, ps[i].label, nTypes);
          neigh[u][t][j] += fval;
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
    ss.str(std::string()); ss << s_bondorient << ".l" << l << tag << ".dat"; fout.open(ss.str(), ios::app);
    for(i=0;i<N;i++){
      tot_neigh = 0.0;
      for(t=0;t<nTypePairs;t++) tot_neigh += neigh[0][t][i];
      for(a=0;a<l_deg; a++) {
        if(tot_neigh>0) qlm[a][i] /= tot_neigh; // qlm(i) = <Ylm>(i) completed
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
    ss.str(std::string()); ss << s_bondorient << ".l" << l << tag << ".ave"; fout.open(ss.str(), ios::app);
    fout << timestep << " " << Q << " " << sqrt((Q2-Q*Q)/N) << endl;
    fout.close();
    if(debug) cout << " * Compute ql(i) (BOO) DONE\n";
    //---- Compute Cl_ij~qlm(i)*qlm(j) (I take the real part) and ql_dot(i) (BOC) ----//
    Cl_ij.resize( bond_list[1].size() );
    for(k=0;k<bond_list[1].size();k++){ // Cij and BOC are computed on 2nd shell neighbours
      Cl_ij.set(k, 0.0);
      i = int2i( bond_list[1][k], N );
      j = int2j( bond_list[1][k], N );
      for(a=0;a<l_deg; a++) {
         Cl_ij[k] += ( real(qlm[a][i])*real(qlm[a][j]) + imag(qlm[a][i])*imag(qlm[a][j]) );
      }
      if( ql[i]==0.0 || ql[j]==0.0) Cl_ij[k]=0.0; // safety condition (NOT JUSTIFIED!)
      else Cl_ij[k] /= (ql[i]*ql[j]); // Cl_ij[k] done
      a = indexOf<int>( ps[i].neigh_list[1], j ); // find index of j in i's neighbour list
      rijSq = ps[i].rijSq_list[1][a]; // and use it to recover the radius
      fval = fcut( rijSq/cutoffSq[1], p1half, p2half );
      tot_neigh = 0.0;
      for(t=0;t<nTypePairs;t++) tot_neigh += neigh[1][t][i];
      ql_dot[i] += fval*Cl_ij[k] / tot_neigh;
      tot_neigh = 0.0;
      for(t=0;t<nTypePairs;t++) tot_neigh += neigh[1][t][j];
      ql_dot[j] += fval*Cl_ij[k] / tot_neigh;
    }
    if(debug) cout << " * Compute C_l(i,j) and ql_dot(i) (BOC) DONE\n";
    //---- Compute global average of ql_dot(i) ----//
    ss.str(std::string()); ss << s_bondcorr << ".l" << l << tag << ".dat"; fout.open(ss.str(), ios::app);
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
    ss.str(std::string()); ss << s_bondcorr << ".l" << l << tag << ".ave"; fout.open(ss.str(), ios::app);
    fout << timestep << " " << Q << " " << sqrt((Q2-Q*Q)/N) << endl;
    fout.close();
    if(debug) cout << " * Compute global average of BOC DONE\n";
    //---- Compute a Locally Confined version of the global BOC ~ average <ql_dot(i)> around i ----//
    ss.str(std::string()); ss << s_bondcorr << ".l" << l << tag << ".local.ave"; fout.open(ss.str(), ios::app);
    for(i=0;i<N;i++){
      for(k=0;k<ps[i].neigh_list[2].size();k++){
        j = ps[i].neigh_list[2][k];
        rijSq = ps[i].rijSq_list[2][k];
        fval = fcut( rijSq/cutoffSq[2], p1half, p2half );
        tot_neigh = 0.0;
        for(t=0;t<nTypePairs;t++) tot_neigh += neigh[2][t][i];
        Ql_dot[i] += fval*ql_dot[j] / tot_neigh;
      }
      fout << timestep << " " << i << " " << Ql_dot[i] << endl;
    }
    fout.close();
    if(debug) cout << " * Compute local average of BOC DONE\n";
    ss.str(std::string()); ss << s_bondcorr << ".l" << l << tag << ".xyz"; fout.open(ss.str(), ios::app);
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
    ss.str(std::string()); ss << s_nxtal << ".l" << l << tag << ".dat"; fout.open(ss.str(), ios::app);
    fout << timestep << " " << m << endl;
    fout.close();
    if(debug) cout << "*** BOO,BOC computation for timestep " << timestep << " ENDED ***\n\n";
}
