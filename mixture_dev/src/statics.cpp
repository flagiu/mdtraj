#include "lib/mdtraj.hpp"
using namespace std;

//---- "init" functions: prepare arrays, averages and headers for output files
//---- "compute" functions: compute quantities for each frame and total average

//------------------------------ Neighbour and Bond list --------------------------------------------//

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
build_neigh() {
  int u,i,j;
  vec rij, rij_mic;
  ntype rijSq, rijSq_mic;
  if(debug) cout << "\n*** Neighbours computation STARTED ***\n";
  //---- Reset counters and lists ----//
  for(u=0;u<maxshell;u++) {
    bond_list[u].clear();
  }
  for(i=0;i<N;i++){
      for(u=0;u<maxshell;u++){
        ps[i].neigh_list[u].clear();
        ps[i].rij_list[u].clear();
        ps[i].rijSq_list[u].clear();
        neigh[u].set(i, 0.0);
      }
  }
  if(debug) cout << " * Reset counters and lists DONE\n";
    //---- Build neighbour list (and save rij vectors) ----//
    for(i=0;i<N;i++){
      for(j=i+1;j<N;j++){
        rij = ps[j].r - ps[i].r; // real distance
        rijSq = rij.sq();
        rij_mic = mic(box, boxInv, rij); // first periodic image
        rijSq_mic = rij_mic.sq();
        if(rijSq_mic < rijSq){ // if closer, choose first periodic image
          rij = rij_mic;
          rijSq = rijSq_mic;
        }
        for(u=0;u<maxshell;u++){
          if(rijSq <= cutoffSq[u]){
            bond_list[u].push_back( ij2int(i,j,N) );

            ps[i].neigh_list[u].push_back(j);
            ps[j].neigh_list[u].push_back(i);

            ps[i].rij_list[u].push_back(   rij);
            ps[j].rij_list[u].push_back(-1*rij);

            ps[i].rijSq_list[u].push_back(rijSq);
            ps[j].rijSq_list[u].push_back(rijSq);
          }
        }
      }
    }
  if(debug) cout << " * Build neighbour and bond lists DONE\n";
  if(debug) cout << "*** Neighbour computation for timestep " << timestep << " ENDED ***\n\n";
  return;
}

//------------------------------------- Coordination Number ---------------------------------------//
template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
init_coordnum()
{
      if(neigh[0].length()!=N) neigh[0].resize(N);
      ss.str(std::string()); ss << s_coordnum << tag << ".dat"; fout.open(ss.str(), ios::out);
      fout << "#Timestep, Particle index, Coordination number. # cutoffs = ";
      for(auto  i=0;i<Nshells;i++) fout << cutoff[i]<<", ";
      fout << endl;
      fout.close();

      ss.str(std::string()); ss << s_coordnum << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "#Timestep, average coordination number, fluctuations. # cutoffs = ";
      for(auto  i=0;i<Nshells;i++) fout << cutoff[i]<<", ";
      fout << endl;
      fout.close();
}

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
compute_coordnum() {
  int i, j, k;
  ntype fval, rijSq;
  vec rij;
  for(i=0;i<N;i++) neigh[0].set(i,0.0); // start counter
  for(i=0;i<N;i++){
    for(k=0;k<ps[i].neigh_list[0].size();k++){ // search in 1st shell neighbour list
      j = ps[i].neigh_list[0][k];
      if(j>i) continue; // avoid double counting!
      rij = ps[i].rij_list[0][k];
      rijSq = ps[i].rijSq_list[0][k];
      fval = fcut( rijSq/cutoffSq[0], p1half, p2half );
      neigh[0][i] += fval;
      neigh[0][j] += fval;
    }
  }
  ss.str(std::string()); ss << s_coordnum << tag << ".dat"; fout.open(ss.str(), ios::app);
  for(i=0;i<N;i++) fout << timestep << " " << i << " " << neigh[0][i] << endl;
  fout.close();
  ss.str(std::string()); ss << s_coordnum << tag << ".ave"; fout.open(ss.str(), ios::app);
  fout << timestep << " " << neigh[0].mean() << " " << neigh[0].std()/sqrt(N) << endl;
  fout.close();
}

//---------------------- Bond Orientational parameters: BOO (q_lm) and BOC (q_lm^dot) -------------------------------//
template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
init_bondorient()
    {
      for(auto i=0;i<Nshells;i++) { if(neigh[i].length()!=N) neigh[i].resize(N); }
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
    ntype fval, rijSq, Q,Q2, ql_factor = sqrt(4*M_PI/l_deg), re,im;
    complex<ntype> Yval;
    vec rij;
    int i,j,k,a,m,u;
    if(debug) cout << "\n*** BOO,BOC computation STARTED ***\n";
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
      for(u=0;u<Nshells;u++){
        for(k=0;k<ps[i].neigh_list[u].size();k++){
          j = ps[i].neigh_list[u][k];
          if(j>i) continue; // avoid double counting!
          rij = ps[i].rij_list[u][k];
          rijSq = ps[i].rijSq_list[u][k];
//          fval = fcut( rijSq/cutoffSq[u], p1half, p2half );
          if(u==0){ //Ider uses a step function for u==0?
            fval = ( rijSq <= cutoffSq[u] ? 1.0 : 0.0 );
          } else {
            fval = fcut( rijSq/cutoffSq[u], p1half, p2half );
          }
          neigh[u][i] += fval;
          neigh[u][j] += fval;
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
        for(a=0;a<l_deg; a++) {
           if(neigh[0][i]>0) qlm[a][i] /= neigh[0][i]; // qlm(i) = <Ylm>(i) completed
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
      ql_dot[i] += fval*Cl_ij[k] / neigh[1][i];
      ql_dot[j] += fval*Cl_ij[k] / neigh[1][j];
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
        Ql_dot[i] += fval*ql_dot[j] / neigh[2][i];
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


//------------------------------- Radial Distribution Function ----------------------------------//
template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
init_rdf()
{
  int t1,t2,a;
  ntype r, shell1, shell2, normalization;
  rdf_binw = 0.5*(L.abs().min()) / (rdf_nbins-1);
  rdf_bins.resize(rdf_nbins);
  ss.str(std::string()); ss << s_rdf << tag << ".traj"; fout.open(ss.str(), ios::out);
  fout << "# First block: Radial distance; other blocks: g_00(r), g_01(r), ... for each frame.\n";
  for( auto i=0; i<rdf_nbins; i++)
  {
    rdf_bins[i] = (i+0.5)*rdf_binw; // take the center of the bin for histogram
    fout << rdf_bins[i] << endl;
  }
  rdf_norm.resize(nTypePairs);
  rdf.resize(nTypePairs);
  rdf_ave.resize(nTypePairs);
  rdf2_ave.resize(nTypePairs);
  for(t1=0;t1<nTypes;t1++)
  {
    for(t2=t1;t2<nTypes;t2++)
    {
      a = types2int(t1,t2,nTypes);
      if(debug) cout << " Inizializzo la coppia (" << t1 <<","<<t2<<")-->"<<  a << "\n";
      rdf_norm[a].resize(rdf_nbins);
      rdf[a].resize(rdf_nbins);
      rdf_ave[a].resize(rdf_nbins);
      rdf2_ave[a].resize(rdf_nbins);
      if(t1!=t2) normalization =     Nt[t2]/V * 4.0*M_PI/3.0 * Nt[t1]/2.0 * 2.0; // multiply by 2 because it's off-diagonal
      else       normalization = (Nt[t2]-1)/V * 4.0*M_PI/3.0 * Nt[t1]/2.0;
      shell1 = 0.0;
      for( auto i=0; i<rdf_nbins; i++)
      {
        r = (i+1)*rdf_binw; // take the upper side of the bin for shells
        shell2 = r*r*r;
        rdf_norm[a][i] = normalization*(shell2 - shell1); // normalize to ideal gas radial distribution
        shell1 = shell2;
        rdf_ave[a][i] = 0.0;
        rdf2_ave[a][i] = 0.0;
      }
    }
  }
  fout.close();
  ss.str(std::string()); ss << s_rdf << tag << ".ave"; fout.open(ss.str(), ios::out);
  fout << "# Radial distance | g_00(r), g_01(r), ... | error for each g(r).\n";
  fout.close();
}

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
compute_rdf(int frameidx)
{
    int i,j, k, tp, ti,tj;
    vec rij;
    ntype r,r_mic;
    for(tp=0;tp<nTypePairs;tp++){
      for(k=0; k<rdf_nbins; k++){
        rdf[tp][k] = 0.0;
      }
    }
    if(debug) cout << "*** RDF computation for timestep " << timestep << " STARTED ***\n";
    for(i=0;i<N-1;i++){
      for(j=i+1;j<N;j++){// i<j
         rij = ps[j].r - ps[i].r; // real distance
         r = rij.norm();
         rij = mic(box, boxInv, rij); // first periodic image
         r_mic = rij.norm();
         if( r_mic < r ) r = r_mic; // if closer, choose first periodic image
         k = int(floor( r/rdf_binw));
         if(k<rdf_nbins){
           ti = ps[i].label;
           tj = ps[j].label;
           tp = types2int(ti,tj,nTypes);
           rdf[tp][k] += 1.0;
         }
      }
    }
    ss.str(std::string()); ss << s_rdf << tag << ".traj"; fout.open(ss.str(), ios::app);
    fout << endl; // start new block
    for(k=0; k<rdf_nbins; k++)
    {
      for(tp=0;tp<nTypePairs;tp++)
      {
        rdf[tp][k] /= rdf_norm[tp][k];
        fout << rdf[tp][k];
        rdf_ave[tp][k] += rdf[tp][k];
        rdf2_ave[tp][k] += rdf[tp][k]*rdf[tp][k];
        if(tp<nTypePairs-1) fout << " ";
        else fout << endl;
      }
    }
    fout.close();
    if(frameidx == (nframes-1)){
        ss.str(std::string()); ss << s_rdf << tag << ".ave"; fout.open(ss.str(), ios::app);
        for(k=0; k<rdf_nbins; k++)
        {
          fout << rdf_bins[k] << " "; // print bin
          for(tp=0;tp<nTypePairs;tp++)
          {
            rdf_ave[tp][k] /= (frameidx+1);
            fout << rdf_ave[tp][k] << " "; // print average g(r)
          }
          for(tp=0;tp<nTypePairs;tp++)
          {
            rdf2_ave[tp][k] /= (frameidx+1);
            fout << sqrt( (rdf2_ave[tp][k]-rdf_ave[tp][k]*rdf_ave[tp][k])/(frameidx+1 -1) ); // print st.dev. of average g(r)
            if(tp<nTypePairs-1) fout << " ";
            else fout << endl;
          }
        }
        fout.close();
        if(debug) cout << "Average RDF printed to file\n";
    }
    if(debug) cout << "*** RDF computation for timestep " << timestep << " ENDED ***\n\n";
    return;
}

//--------------------- Angular Distribution Function (distribution of cos(theta_ijk) ----------------------------------//
// must be preceded by computation of nearest neighbours
template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
init_adf()
{
      adf_bins.resize(adf_nbins);
      adf.resize(adf_nbins);
      adf_ave.resize(adf_nbins);
      adf2_ave.resize(adf_nbins);
      adf_binw = 2.0 / (adf_nbins) ; // -1 < cos(angle) < 1
      ntype angle;
      ss.str(std::string()); ss << s_adf << tag << ".traj"; fout.open(ss.str(), ios::out);
      fout << "# First block: cos(angle); other blocks: ADF for each frame )\n";
      for( auto i=0; i<adf_nbins; i++){
        adf_bins[i] = -1.0 + (i+0.5)*adf_binw; // take the center of the bin for histogram
        fout << adf_bins[i] << endl;
        adf_ave[i] = 0.0;
        adf2_ave[i] = 0.0;
      }
      fout.close();
      ss.str(std::string()); ss << s_adf << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "# cos(angle), ADF, ADF error.\n";
      fout.close();
}

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
compute_adf(int frameidx)
{
    int i,j,k, a,b, bin, counts=0;
    vec rij, rik;
    ntype rijSq, rikSq, costheta;
    for(a=0; a<adf_nbins; a++) adf[a] = 0.0;
    if(debug) cout << "*** ADF computation for timestep " << timestep << " STARTED ***\n";
    for(i=0;i<N;i++){
      for(a=1;a<ps[i].neigh_list[0].size();a++){
          j = ps[i].neigh_list[0][a];
          if(j>i) continue; // avoid double counting!
          rij = ps[i].rij_list[0][a];
          rijSq = ps[i].rijSq_list[0][a];
          for(b=0;b<a;b++){
            k = ps[i].neigh_list[0][b]; //useless?
            rik = ps[i].rij_list[0][b];
            rikSq = ps[i].rijSq_list[0][b];
            costheta = (rij*rik) / sqrt(rijSq*rikSq);
            bin = int(floor( (costheta+1.0)/adf_binw)); // min value is -1.0 !
            if(bin<adf_nbins){
              adf[bin] += 1.0;
              counts++;
            }
          }
      }
    }
    ss.str(std::string()); ss << s_adf << tag << ".traj"; fout.open(ss.str(), ios::app);
    fout << endl; // start new block
    for(k=0; k<adf_nbins; k++){
      adf[k] /= (counts*adf_binw); // ci sta la larghezza del bin?
      fout << adf[k] << endl;
      adf_ave[k] += adf[k];
      adf2_ave[k] += adf[k]*adf[k];
    }
    fout.close();
    if(frameidx == (nframes-1)){
        ss.str(std::string()); ss << s_adf << tag << ".ave"; fout.open(ss.str(), ios::app);
        for(k=0; k<adf_nbins; k++){
          adf_ave[k] /= (frameidx+1);
          adf2_ave[k] /= (frameidx+1);
          fout << adf_bins[k] << " " << adf_ave[k] << " " << sqrt( (adf2_ave[k]-adf_ave[k]*adf_ave[k])/(frameidx+1 -1) ) << endl;
        }
        fout.close();
        if(debug) cout << "Average ADF printed to file\n";
    }
    if(debug) cout << "*** ADF computation for timestep " << timestep << " ENDED ***\n\n";
    return;
}

//------------ Minimum atomic distance ----------------//
template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
init_rmin()
{
  ss.str(std::string()); ss << s_rmin << tag << ".dat"; fout.open(ss.str(), ios::out);
  fout << "# Minimum atomic distance. # cutoff1 = " << cutoff[0]<< endl;
  fout.close();
}

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
compute_rmin()
{
  ntype rSq, rminSq = cutoffSq[0];
  int i,j, k;
  if(debug) cout << "*** RMIN computation for timestep " << timestep << " STARTED ***\n";
  for(i=0;i<N;i++){
    for(k=0;k<ps[i].neigh_list[0].size();k++){
      j = ps[i].neigh_list[0][k];
      if(j>i) continue; // avoid double counting!
      rSq = ps[i].rijSq_list[0][k];
      if(rSq<rminSq) rminSq=rSq;
    }
  }
  ss.str(std::string()); ss << s_rmin << tag << ".dat"; fout.open(ss.str(), ios::app);
  fout << sqrt(rminSq) << endl;
  fout.close();
  if(debug) cout << "*** RMIN computation for timestep " << timestep << " ENDED ***\n";
  return;
}

//--------------------- ALTBC ----------------------------------//
// must be preceded by computation of nearest neighbours (?)
template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
init_altbc()
{
      altbc_bins.resize(altbc_nbins);
      altbc.resize(altbc_nbins*altbc_nbins);
      altbc_ave.resize(altbc_nbins*altbc_nbins);
      altbc2_ave.resize(altbc_nbins*altbc_nbins);
      altbc_binw = (cutoff[0]-altbc_rmin) / (altbc_nbins-1);
      altbc_cos = cos( altbc_angle/180. * M_PI ); // should be close to 1.0
      if(debug) cout << "ALTBC cos threshold = " << altbc_cos<<endl;
      ss.str(std::string()); ss << s_altbc << tag << ".traj"; fout.open(ss.str(), ios::out);
      fout << "# 1st block: radial distance; other blocks: ALTBC 2D matrix for each frame; # deviation >= " << altbc_angle << " degrees\n";
      for( auto i=0; i<altbc_nbins; i++){
        altbc_bins[i] = altbc_rmin + (i+0.5)*altbc_binw; // take the center of the bin for histogram
        fout << altbc_bins[i] << endl;
      }
      fout << endl; // end of block
      for( auto i=0; i<altbc.length(); i++){
        altbc_ave[i] = 0.0;
        altbc2_ave[i] = 0.0;
      }
      fout.close();
      ss.str(std::string()); ss << s_altbc << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "# 1st block: radial distance; 2nd block: ALTBC 2D matrix; # deviation >= " << altbc_angle << " degrees\n";
      fout.close();
}

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
compute_altbc(int frameidx)
{
    int i,j,k,k0,k1, a,b, bin0, bin1, counts=0;
    vec rij, rik;
    ntype rijSq, rikSq, costheta, rijNorm, rikNorm, norm_factor= V*V/(N*(N-1)*(N-2));
    for(a=0; a<altbc.length(); a++) altbc[a] = 0.0;
    if(debug) cout << "*** altbc computation for timestep " << timestep << " STARTED ***\n";
    for(i=0;i<N;i++){
      for(a=1;a<ps[i].neigh_list[0].size();a++){
          j = ps[i].neigh_list[0][a];
          if(j>i) continue; // avoid double counting!
          rij = ps[i].rij_list[0][a];
          rijSq = ps[i].rijSq_list[0][a];
          rijNorm = sqrt(rijSq);
          if(rijNorm < altbc_rmin) continue;
          for(b=0;b<a;b++){
            k = ps[i].neigh_list[0][b]; //useless?
            rik = ps[i].rij_list[0][b];
            rikSq = ps[i].rijSq_list[0][b];
            rikNorm = sqrt(rikSq);
            if(rikNorm < altbc_rmin) continue;
            costheta = (rij*rik) / (rijNorm*rikNorm);
            bin0 = int(floor( (rijNorm-altbc_rmin)/altbc_binw ));
            bin1 = int(floor( (rikNorm-altbc_rmin)/altbc_binw ));
            if(bin0<altbc_nbins && bin1<altbc_nbins && fabs(costheta) >= altbc_cos){
              altbc[bin0 + altbc_nbins*bin1] += 1.0;
              counts++;
            }
          }
      }
    }
    ss.str(std::string()); ss << s_altbc << tag << ".traj"; fout.open(ss.str(), ios::app);
    fout << endl; // start new block
    for(k0=0; k0<altbc_nbins; k0++){
      for(k1=0; k1<altbc_nbins; k1++){
        k = k0 + altbc_nbins*k1;
        altbc[k] /= (counts*altbc_binw*altbc_binw); // ci sta la larghezza del bin?
        altbc[k] *= norm_factor;
        fout << altbc[k] << " "; // 2D matrix
        altbc_ave[k] += altbc[k];
        altbc2_ave[k] += altbc[k]*altbc[k];
      }
      fout << endl;
    }
    fout << endl; // end of block
    fout.close();
    if(frameidx == (nframes-1)){
        ss.str(std::string()); ss << s_altbc << tag << ".ave"; fout.open(ss.str(), ios::app);
        for(k0=0; k0<altbc_nbins; k0++){
          fout << altbc_bins[k0] << endl; // bins
        }
        fout << endl; // end of block
        for(k0=0; k0<altbc_nbins; k0++){
          for(k1=0; k1<altbc_nbins; k1++){
            k = k0 + altbc_nbins*k1;
            altbc_ave[k] /= (frameidx+1);
            altbc2_ave[k] /= (frameidx+1);
            fout << altbc_ave[k] << " "; // average
          }
          fout << endl;
        }
        /*
        fout << endl; // end of block
        for(k0=0; k0<altbc_nbins; k0++){
          for(k1=0; k1<altbc_nbins; k1++){
            k = k0 + altbc_nbins*k1;
            fout << sqrt( (altbc2_ave[k]-altbc_ave[k]*altbc_ave[k])/(frameidx+1 -1) ) << " "; // error
          }
          fout << endl;
        }
        */
        fout.close();
        if(debug) cout << "Average ALTBC printed to file\n";
    }
    if(debug) cout << "*** ALTBC computation for timestep " << timestep << " ENDED ***\n\n";
    return;
}