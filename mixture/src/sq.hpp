#ifndef _SQ_H_
#define _SQ_H_
using namespace std;
//------------------------------- Static Structure Factor ----------------------------------//
const string qvectors_path="/home/flavio/programmi/mdtraj/QVECTORS";

template <class ntype, class ptype>
class SQ_Calculator
{
  using vec=myvec<ntype,3>;
  using mat=mymatrix<ntype,3,3>;
  private:
    ntype binw;
    int qm, qM, dq, nbins;
    vecflex<ntype> bins, norm, value, value2, ave, ave2;
    string string_out, myName, tag;
    fstream fout;
    stringstream ss;
    bool debug, verbose;
  public:
    SQ_Calculator(){
      myName = "SQ";
    }
    virtual ~SQ_Calculator(){}

    void init(int q_mod_min, int q_mod_max, int q_mod_step, vec box_diagonal,
      string string_out_, string tag_, bool debug_, bool verbose_)
    {
      qm = q_mod_min;
      qM = q_mod_max;
      dq = q_mod_step;
      string_out = string_out_;
      tag = tag_;
      debug = debug_;
      verbose = verbose_;

      nbins = int(floor( (qM-qm)/dq )) + 1;
      binw = M_PI / box_diagonal[0]; // half mesh: it is in units of pi/L
      bins.resize(nbins);
      value.resize(nbins);
      value2.resize(nbins);
      ave.resize(nbins);
      ave2.resize(nbins);

      if(verbose)
      {
        ss.str(std::string()); ss << string_out << tag << ".traj"; fout.open(ss.str(), ios::out);
        fout << "# First block: q-wave-vectors; other blocks: S(q)\n";
      }
      for( auto i=0; i<nbins; i++){
        bins[i] = (1+qm+i*dq)*binw; // q-wave-vector values
        if(verbose) fout << bins[i] << endl;
        ave[i] = 0.0;
        ave2[i] = 0.0;
      }
      if(verbose) fout.close();
      ss.str(std::string()); ss << string_out << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "# q-wave-vector, <S(q)>, <S(q)> error.\n";
      fout.close();
    }

    void compute(int frameidx, int nframes, int timestep, vector<ptype> ps)
    {
      const int N=ps.size();
      int i,k,qmod, qx,qy,qz, qcount, j,idx;
      ntype arg, rho_real, rho_imag, sq_oriented;
      string line,a,b,c;
      const int num_k_values = 1 + qM/2; // 1/2 is the mesh spacing!
      const int exp_storage_size = num_k_values * N * 3 * 2; // wavenumber, particle, cartesian, real/imag
      const ntype dk = 2.*binw; // 2*M_PI/L[0];
      ntype x0,x1, y0,y1, z0,z1;
      vecflex<ntype> exps;
      exps.resize(exp_storage_size);

      for(k=0;k<num_k_values;k++) // wavenumber
      {
        for(i=0;i<N;i++) // particle
        {
          for(j=0;j<3;j++) // x,y,z
          {
            idx = 2*3*N*k + 2*3*i + 2*j;
            // with real coordinates:
            //exps[ idx     ] = cos( dk*k*ps[i].r[j] );
            //exps[ idx + 1 ] = sin( dk*k*ps[i].r[j] );
            // with fractional coordinates (it is more generalizable):
            exps[ idx     ] = cos( 2*M_PI * k * ps[i].s[j] ); // dk*k*ps[i].r[j] );
            exps[ idx + 1 ] = sin( 2*M_PI * k * ps[i].s[j] ); // dk*k*ps[i].r[j] );
          }
        }
      }
      if(debug) cout << "*** "<<myName<<" computation for timestep " << timestep << " STARTED ***\n";
      for(k=0; k<nbins; k++)
      {
        qmod = qm + k*dq;
        if(debug) cout << "  k="<<k<<", qmod="<<qmod<<endl;
        ss.str(std::string()); ss << qvectors_path << "/qvector." << setw(3) << setfill('0') << qmod; fout.open(ss.str(), ios::in);
        if(debug) cout << "  opening qvectors: " << ss.str() << endl;

        qcount=0;
        value[k]=value2[k]=0.0;
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
            /* // Brute force method
            arg = dk * (qx*ps[i].r[0] + qy*ps[i].r[1] + qz*ps[i].r[2]);
            rho_real += cos(arg);
            rho_imag += sin(arg);
            */
            // Pre-computation method
            x0 = exps[ 2*3*N*abs(qx) + 2*3*i + 2*0     ]; // cos( kx * rx )
            x1 = exps[ 2*3*N*abs(qx) + 2*3*i + 2*0 + 1 ]; // sin(|kx|* rx )
            y0 = exps[ 2*3*N*abs(qy) + 2*3*i + 2*1     ]; // cos( ky * ry )
            y1 = exps[ 2*3*N*abs(qy) + 2*3*i + 2*1 + 1 ]; // sin(|ky|* ry )
            z0 = exps[ 2*3*N*abs(qz) + 2*3*i + 2*2     ]; // cos( kz * rz )
            z1 = exps[ 2*3*N*abs(qz) + 2*3*i + 2*2 + 1 ]; // sin(|kz|* rz )
            if(qx<0) x1 = -x1; // give correct sign to sin() !
            if(qy<0) y1 = -y1;
            if(qz<0) z1 = -z1;
            rho_real += ( x0*y0*z0 - x0*y1*z1 - x1*y0*z1 - x1*y1*z0 ); // cos( k*r )
            rho_imag += ( x1*y0*z0 + x0*y1*z0 + x0*y0*z1 - x1*y1*z1 ); // sin( k*r )
          }
          sq_oriented = ( rho_real*rho_real + rho_imag*rho_imag ) / N;
          if(debug) cout << "  S(qx,qy,qz) = "<<sq_oriented<<endl;
          value[k] += sq_oriented;
          value2[k] += sq_oriented*sq_oriented;
          qcount++;
        }
        value[k] /= qcount;
        value2[k] /= qcount;
        value2[k] = sqrt( (value2[k]-value[k]*value[k])/(qcount-1) ); // fluctuations of S(q;t0) over q
        if(debug) cout << "  Total "<<qcount<<" triplets of wavenumbers.\n";
        fout.close();
      }

      if(verbose)
      {
        ss.str(std::string()); ss << string_out << tag << ".traj"; fout.open(ss.str(), ios::app);
        fout << endl; // start new block
      }
      for(k=0; k<nbins; k++){
        if(verbose) fout << value[k] << endl;
        ave[k] += value[k];
        ave2[k] += value[k]*value[k];
      }
      if(verbose) fout.close();

      if(frameidx == (nframes-1)){
          ss.str(std::string()); ss << string_out << tag << ".ave"; fout.open(ss.str(), ios::app);
          for(k=0; k<nbins; k++){
            ave[k] /= (frameidx+1);
            ave2[k] /= (frameidx+1);
            fout << bins[k] << " " << ave[k] << " " << sqrt( (ave2[k]-ave[k]*ave[k])/(frameidx+1 -1) ) << endl;
          }
          fout.close();
          if(debug) cout << "Average "<<myName<<" printed to file\n";
      }
      if(debug) cout << "*** "<<myName<<" computation for timestep " << timestep << " ENDED ***\n\n";
      return;
    }
};

#endif
