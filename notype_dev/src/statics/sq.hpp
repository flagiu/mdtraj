using namespace std;
//------------------------------- Static Structure Factor ----------------------------------//
template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
init_sq()
{
  sq_nbins = int(floor( (qmodmax-qmodmin)/qmodstep )) + 1;
  sq_binw = M_PI / L[0]; // half mesh: it is in units of pi/L
  sq_norm.resize(sq_nbins);
  sq_bins.resize(sq_nbins);
  sq.resize(sq_nbins);
  sq2.resize(sq_nbins);
  sq_ave.resize(sq_nbins);
  sq2_ave.resize(sq_nbins);
  ss.str(std::string()); ss << s_sq << tag << ".traj"; fout.open(ss.str(), ios::out);
  fout << "# First block: q-wave-vectors; other blocks: S(q)\n";
  for( auto i=0; i<sq_nbins; i++){
    sq_bins[i] = (1+qmodmin+i*qmodstep)*sq_binw; // q-wave-vector values
    fout << sq_bins[i] << endl;
    sq_ave[i] = 0.0;
    sq2_ave[i] = 0.0;
  }
  fout.close();
  ss.str(std::string()); ss << s_sq << tag << ".ave"; fout.open(ss.str(), ios::out);
  fout << "# q-wave-vector, <S(q)>, <S(q)> error.\n";
  fout.close();
}

template <class ntype, class ptype>
void Trajectory<ntype, ptype>::
compute_sq(int frameidx)
{
  int i,k,qmod, qx,qy,qz, qcount;
  ntype arg, rho_real, rho_imag, sq_oriented;
  string line,a,b,c;
  if(debug) cout << "*** Sq computation for timestep " << timestep << " STARTED ***\n";
  for(k=0; k<sq_nbins; k++)
  {
    qmod = qmodmin + k*qmodstep;
    if(debug) cout << "  k="<<k<<", qmod="<<qmod<<endl;
    ss.str(std::string()); ss << qvectors_path << "/qvector." << setw(3) << setfill('0') << qmod; fout.open(ss.str(), ios::in);
    if(debug) cout << "  opening qvectors: " << ss.str() << endl;
    qcount=0;
    sq[k]=sq2[k]=0.0;
    while( getline(fout,line) )
    {
      istringstream(line) >> a >> b >> c;
      qx = stoi(a);
      qy = stoi(b);
      qz = stoi(c);
      if(debug) cout << "  read kx,ky,kz = " << qx<<", "<< qy<<", "<< qz<<endl;
      rho_real=rho_imag=0.0;
      for(i=0;i<N;i++)
      {
        arg = 2. * M_PI / L[0] * (qx*ps[i].r[0] + qy*ps[i].r[1] + qz*ps[i].r[2]);
        rho_real += cos(arg);
        rho_imag += sin(arg);
      }
      sq_oriented = ( rho_real*rho_real + rho_imag*rho_imag ) / N;
      if(debug) cout << "  S(qx,qy,qz) = "<<sq_oriented<<endl;
      sq[k] += sq_oriented;
      sq2[k] += sq_oriented*sq_oriented;
      qcount++;
    }
    sq[k] /= qcount;
    sq2[k] /= qcount;
    sq2[k] = sqrt( (sq2[k]-sq[k]*sq[k])/(qcount-1) ); // fluctuations of S(q;t0) over q
    if(debug) cout << "  Total "<<qcount<<" triplets of wavenumbers.\n";
    fout.close();
  }
  ss.str(std::string()); ss << s_sq << tag << ".traj"; fout.open(ss.str(), ios::app);
  fout << endl; // start new block
  for(k=0; k<sq_nbins; k++){
    fout << sq[k] << endl;
    sq_ave[k] += sq[k];
    sq2_ave[k] += sq[k]*sq[k];
  }
  fout.close();
  if(frameidx == (nframes-1)){
      ss.str(std::string()); ss << s_sq << tag << ".ave"; fout.open(ss.str(), ios::app);
      for(k=0; k<sq_nbins; k++){
        sq_ave[k] /= (frameidx+1);
        sq2_ave[k] /= (frameidx+1);
        fout << sq_bins[k] << " " << sq_ave[k] << " " << sqrt( (sq2_ave[k]-sq_ave[k]*sq_ave[k])/(frameidx+1 -1) ) << endl;
      }
      fout.close();
      if(debug) cout << "Average Sq printed to file\n";
  }
  if(debug) cout << "*** Sq computation for timestep " << timestep << " ENDED ***\n\n";
  return;
}
