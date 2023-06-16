#include "lib/mdtraj.hpp"
using namespace std;

//---- "init" functions: prepare arrays, averages and headers for output files
//---- "compute" functions: compute quantities for each frame and total average

//------------------------------ Neighbour and Bond list --------------------------------------------//

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
build_neigh() {
  int u,i,j;
  vec rij;
  ntype rijSq;
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
        rij = ps[j].r - ps[i].r;
        //rij = ps[i].mic(rij, L); // minimum image convention
        rij = mic(box, boxInv, rij);
        rijSq = rij.sq();
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
    ntype fval;
    vec rij;
    for(auto i=0;i<N;i++) neigh[0].set(i,0.0); // start counter
    for(auto i=0;i<N-1;i++){
      for(auto j=i+1;j<N;j++){
        rij = ps[j].r - ps[i].r;
        //rij = ps[i].mic(rij, L); // minimum image convention
        rij = mic(box, boxInv, rij);
        fval = fcut( rij.sq()/cutoffSq[0], p1half, p2half );
        neigh[0][i] += fval;
        neigh[0][j] += fval;
//        if( ps[j].is_inside(ps[i], cutoff[0], L) ){
//          neigh[0][i]++;
//        }
      }
    }
    ss.str(std::string()); ss << s_coordnum << tag << ".dat"; fout.open(ss.str(), ios::app);
    for(auto i=0;i<N;i++) fout << timestep << " " << i << " " << neigh[0][i] << endl;
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
      Cl_ij[k] /= (ql[i]*ql[j]); // Cl_ij[k] done
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
      rdf_norm.resize(rdf_nbins);
      rdf_bins.resize(rdf_nbins);
      rdf.resize(rdf_nbins);
      rdf_ave.resize(rdf_nbins);
      rdf2_ave.resize(rdf_nbins);
      rdf_binw = 0.5*(L.abs().min()) / (rdf_nbins-1) ;
      ntype r, shell1, shell2, normalization = (N-1)/V * 4.0*M_PI/3.0 * N/2.0;
      shell1 = 0.0;
      ss.str(std::string()); ss << s_rdf << tag << ".traj"; fout.open(ss.str(), ios::out);
      fout << "# First block: radial distance; other blocks: g(r)\n";
      for( auto i=0; i<rdf_nbins; i++){
        rdf_bins[i] = (i+0.5)*rdf_binw; // take the center of the bin for histogram
        fout << rdf_bins[i] << endl;
        r = (i+1)*rdf_binw; // take the upper side of the bin for shells
        shell2 = r*r*r;
        rdf_norm[i] = normalization*(shell2 - shell1); // normalize to ideal gas radial distribution
        shell1 = shell2;
        rdf_ave[i] = 0.0;
        rdf2_ave[i] = 0.0;
      }
      fout.close();
      ss.str(std::string()); ss << s_rdf << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "# Radial distance, g(r), g(r) error.\n";
      fout.close();
}

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
compute_rdf(int frameidx)
{
    int i,j, k;
    vec rij;
    ntype r;
    for(k=0; k<rdf_nbins; k++)
      rdf[k] = 0.0;
    if(debug) cout << "*** RDF computation for timestep " << timestep << " STARTED ***\n";
    for(i=0;i<N;i++){
      for(j=i+1;j<N;j++){
         //rij = ( ps[j].r - ps[i].r ).mic(L);
         rij = mic(box, boxInv, ps[j].r - ps[i].r );
         r = rij.norm();
         k = int(floor( r/rdf_binw));
         if(k<rdf_nbins){
           rdf[k] += 1.0;
         }
      }
    }
    ss.str(std::string()); ss << s_rdf << tag << ".traj"; fout.open(ss.str(), ios::app);
    fout << endl; // start new block
    for(k=0; k<rdf_nbins; k++){
      rdf[k] /= rdf_norm[k];
      fout << rdf[k] << endl;
      rdf_ave[k] += rdf[k];
      rdf2_ave[k] += rdf[k]*rdf[k];
    }
    fout.close();
    if(frameidx == (nframes-1)){
        ss.str(std::string()); ss << s_rdf << tag << ".ave"; fout.open(ss.str(), ios::app);
        for(k=0; k<rdf_nbins; k++){
          rdf_ave[k] /= (frameidx+1);
          rdf2_ave[k] /= (frameidx+1);
          fout << rdf_bins[k] << " " << rdf_ave[k] << " " << sqrt( (rdf2_ave[k]-rdf_ave[k]*rdf_ave[k])/(frameidx+1 -1) ) << endl;
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

