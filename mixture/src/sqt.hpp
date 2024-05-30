#ifndef _SQT_H_
#define _SQT_H_
using namespace std;
//-------- Dynamic Structure Factor S(q,t), a.k.a. Intermediate Scattering Function F(q,t) -----------------------//
#ifndef _SQ_H_
const string qvectors_path="/home/flavio/programmi/mdtraj/QVECTORS";
#define MAX_QVECTORS 300
#endif

template <class ntype, class ptype>
class SQT_Calculator
{
  using vec=myvec<ntype,3>;
  using mat=mymatrix<ntype,3,3>;
  private:
    ntype binw;
    int qm, qM, dq, nbins;
    int period_in_dt_units, num_periods;
    vecflex<ntype> bins, rho_real, rho_imag, ave, ave2;
    string string_out, myName, tag, rho_tmpfile_prefix;
    fstream fout, f_rho;
    stringstream ss;
    bool debug, verbose;

    string get_rho_tmpfile(int idx) {
      ss.str(std::string());
      ss << rho_tmpfile_prefix << "_idx" << idx << ".dat";
      return ss.str();
    }
    void write_rho_tmpfile(int idx, ntype re, ntype im) {
      ofstream ftmp;
      ftmp.open(get_rho_tmpfile(idx),ios::out);
      ftmp<<re<<" "<<im<<endl;
      ftmp.close();
    }
    void append_rho_tmpfile(int idx, ntype re, ntype im) {
      ofstream ftmp;
      ftmp.open(get_rho_tmpfile(idx),ios::app);
      ftmp<<re<<" "<<im<<endl;
      ftmp.close();
    }
    void read_rho_tmpfile(int q, int idx, ntype *re, ntype *im) {
      ifstream ftmp;
      string line, a,b;
      int qq=-1;
      ftmp.open(get_rho_tmpfile(idx),ios::in);
      if (ftmp.is_open()){
        while ( getline(ftmp,line) && qq<q ){
          istringstream(line) >> a >> b;
          qq++;
        }
        ftmp.close();
        *re=stof(a);
        *im=stof(b);
      }
      else {
        cout << "ERROR: Unable to open file "<<get_rho_tmpfile(idx)<<endl;
        exit(1);
      }
    }

  public:
    SQT_Calculator(){
      myName = "SQT";
      rho_tmpfile_prefix = "_tmp_rho";
    }
    virtual ~SQT_Calculator(){}

    void init(int q_mod_min, int q_mod_max, int q_mod_step, int dtframe,
      int nframes, int period_in_real_units, vec box_diagonal,
      string string_out_, string tag_, bool debug_, bool verbose_)
    {
      int i,j,idx;
      string_out = string_out_;
      tag = tag_;
      debug = debug_;
      verbose = verbose_;
      qm = q_mod_min;
      qM = q_mod_max;
      dq = q_mod_step;
      if(period_in_real_units>0 && period_in_real_units<nframes*dtframe){
        period_in_dt_units = ceil(period_in_real_units/dtframe);
      } else{
        period_in_dt_units = nframes;
      }
      num_periods = nframes/period_in_dt_units;
      if(debug) cout << myName<<": Set period_in_dt_units = " << period_in_dt_units << ", num_periods = " << num_periods << endl;

      nbins = int(floor( (qM-qm)/dq )) + 1;
      binw = M_PI / box_diagonal[0]; // half mesh: it is in units of pi/L
      bins.resize(nbins);
      for(i=0; i<nbins; i++) bins[i] = (1+qm+i*dq)*binw; // q-wave-vector values
      // rho(q,t0+t) for each q,t
      rho_real.resize(nbins * period_in_dt_units);
      rho_imag.resize(nbins * period_in_dt_units);
      // S(q,t0+t) for each q,t
      ave.resize(nbins * period_in_dt_units);
      ave2.resize(nbins * period_in_dt_units);

      if(verbose)
      {
        ss.str(std::string()); ss << string_out << tag << ".traj"; fout.open(ss.str(), ios::out);
        fout << "# 1st block: q-wave-vectors; 2nd block: Delta timestep; Other blocks: S(q,t;t0) trajectory for each frame\n";
        for(i=0; i<nbins; i++) fout << bins[i] << endl;
        fout<<endl; // new block
        for(i=0;i<period_in_dt_units;i++)
          fout << i*dtframe << endl; // Delta timestep values
        fout.close();
      }

      // reset rho, <S(q,dt)>, <S(q,dt)^2> for each q=qmin,qmin+dq,...,qmax
      //    and for each dt=0,1,...,period_in_dt_units-1
      for(i=0; i<nbins; i++)
      {
        for(j=0;j<period_in_dt_units;j++)
        {
          idx = period_in_dt_units*i + j;
          ave[idx]=ave2[idx]=0.0;
        }
      }
      ss.str(std::string()); ss << string_out << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "# 1st block: q-wave-vector; 2nd block: Delta timestep; 3rd block: <S(q,t)> matrix; 4th block: <S(q,t)> error matrix.\n";
      for(i=0; i<nbins; i++) fout << bins[i] << endl;
      fout<<endl; // new block
      for(i=0;i<period_in_dt_units;i++) fout << i*dtframe << endl;
      fout.close();
    }

    void compute(int frameidx, int nframes, int timestep, vector<ptype> ps)
    {
      const int N=ps.size();
      const ntype invSqrtN=1.0/sqrt((ntype)N);
      int i,k,qmod, qx,qy,qz, q,num_qvectors, j,idx, idx_old, dframe, nperiod;
      int qvectors[3*MAX_QVECTORS];
      ntype rhoQvecRe,rhoQvecIm, rhoQvecRe_old,rhoQvecIm_old;
      ntype sqtQvecRe, sqtQvecIm;
      string line,a,b,c;
      const int num_k_values=1+qM/2; // 1/2 is the mesh spacing!
      const int exp_storage_size=num_k_values*N*3*2; // wavenumber, particle, cartesian, real/imag
      const ntype dk=2.*binw; // 2*M_PI/L[0];
      ntype x0,x1, y0,y1, z0,z1;
      vecflex<ntype> exps, ave_t0, ave2_t0;
      exps.resize(exp_storage_size);
      ave_t0.resize(nbins);
      ave2_t0.resize(nbins);

      nperiod = frameidx / period_in_dt_units;
      dframe = frameidx % period_in_dt_units;

      if(debug) cout << "*** "<<myName<<" computation for timestep " << timestep << " STARTED ***\n";
      if(dframe==0)
      { // reset rho(t0+t) for each q,t at the beginning of each cycle t=0
        for(i=0;i<nbins;i++)
        {
          for(j=0;j<period_in_dt_units;j++)
          {
            idx = period_in_dt_units*i + j;
            rho_real[idx]=rho_imag[idx]=0.0;
          }
        }
      }

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

      for(k=0; k<nbins; k++)
      {
        idx = period_in_dt_units*k + dframe;
        qmod = qm + k*dq; // in units of half mesh binw=pi/L
        if(debug) cout << "  k="<<k<<", qmod="<<qmod<<endl;
        // qvector.${qmod} file contains wavenumbers whose
        // modulus is in (qmod-0.5,qmod+0.5] in half mesh units
        ss.str(std::string()); ss << qvectors_path << "/qvector." << setw(3) << setfill('0') << qmod; fout.open(ss.str(), ios::in);
        if(debug) cout << "  opening qvectors: " << ss.str() << endl;
        num_qvectors=0;
        while( getline(fout,line) && num_qvectors<MAX_QVECTORS )
        {
          istringstream(line) >> a >> b >> c;
          qx = stoi(a); // in units of full mesh 2*pi/L
          qy = stoi(b);
          qz = stoi(c);
          if(debug) cout << "  read kx,ky,kz = " << qx<<", "<< qy<<", "<< qz<<endl;
          qvectors[3*num_qvectors  ]=qx;
          qvectors[3*num_qvectors+1]=qy;
          qvectors[3*num_qvectors+2]=qz;
          num_qvectors++;
        }
        fout.close();
        if(debug) cout << "  Total "<<num_qvectors<<" triplets of wavenumbers.\n";

        ave_t0[k]=ave2_t0[k]=0.0;
        for(q=0;q<num_qvectors;q++) {
          rhoQvecRe=rhoQvecIm=0.0;
          qx=qvectors[3*q  ];
          qy=qvectors[3*q+1];
          qz=qvectors[3*q+2];
          for(i=0;i<N;i++) {
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
            // Compute rho(t0+t)
            rhoQvecRe += ( x0*y0*z0 - x0*y1*z1 - x1*y0*z1 - x1*y1*z0 ); // cos( k*r )
            rhoQvecIm += ( x1*y0*z0 + x0*y1*z0 + x0*y0*z1 - x1*y1*z1 ); // sin( k*r )
          }
          rhoQvecRe*=invSqrtN;
          rhoQvecIm*=invSqrtN;
          if(q==0) write_rho_tmpfile(idx,rhoQvecRe,rhoQvecIm);
          else    append_rho_tmpfile(idx,rhoQvecRe,rhoQvecIm);
          // Project rho(t0+t) onto rho(t0)
          idx_old = period_in_dt_units*k + 0;
          read_rho_tmpfile(q,idx_old,&rhoQvecRe_old,&rhoQvecIm_old);
          sqtQvecRe = rhoQvecRe*rhoQvecRe_old + rhoQvecIm*rhoQvecIm_old;
          if(debug) cout << "  Re[ S(qx,qy,qz; dt="<<dframe<<") ] = "<<sqtQvecRe<<endl;
          sqtQvecIm = -rhoQvecRe*rhoQvecIm_old + rhoQvecIm*rhoQvecRe_old;
          if(debug) cout << "  Im[ S(qx,qy,qz; dt="<<dframe<<") ] = "<<sqtQvecIm<<endl;
          // Average over orientations at fixed q,t0,t
          ave_t0[k] += sqtQvecRe;
          ave2_t0[k] += sqtQvecRe*sqtQvecRe;
        }
        ave_t0[k] /= num_qvectors;
        ave2_t0[k] /= num_qvectors;
        if(num_qvectors>1){
          // fluctuations of S(q;t0) over q orientation
          ave2_t0[k] = sqrt( (ave2_t0[k]-ave_t0[k]*ave_t0[k])/(num_qvectors-1) );
        }
      }

      if(verbose)
      {
        ss.str(std::string()); ss << string_out << tag << ".traj"; fout.open(ss.str(), ios::app);
        fout << endl; // start new block
        for(k=0; k<nbins; k++) fout << ave_t0[k] << endl;
        fout.close();
      }

      for(k=0; k<nbins; k++){
        idx = period_in_dt_units*k + dframe;
        ave[idx] += ave_t0[k]; // average over t0 at fixed q,t
        ave2[idx] += ave_t0[k]*ave_t0[k];
      }

      if(frameidx == (nframes-1))
      {
          ss.str(std::string()); ss << string_out << tag << ".ave"; fout.open(ss.str(), ios::app);
          fout<<endl; // new block
          for(k=0; k<nbins; k++)
          {
            for(j=0;j<period_in_dt_units;j++)
            {
              idx = period_in_dt_units*k + j;
              ave[idx] /= num_periods; // average over t0 at fixed q,t
              fout << ave[idx] << " ";
            }
            fout << endl;
          }

          fout << endl; // new block
          for(k=0; k<nbins; k++){
            for(j=0;j<period_in_dt_units;j++){
              idx = period_in_dt_units*k + j;
              ave2[idx] /= num_periods;
              if(num_periods>1){
                ave2[idx] = sqrt( (ave2[idx]-ave[idx]*ave[idx])/(num_periods-1) );
              }
              fout << ave2[idx] << " ";
            }
            fout << endl;
          }
          fout.close();
          if(debug) cout << "Average "<<myName<<" printed to file\n";

          if(!debug){ // clear temporary files
            for(k=0; k<nbins; k++){
              for(j=0;j<period_in_dt_units;j++){
                idx = period_in_dt_units*k + j;
                remove(get_rho_tmpfile(idx).c_str());
              }
            }
            if(debug) cout << "Cleared temporary files of "<<myName<<"\n";
          }
      }
      if(debug) cout << "*** "<<myName<<" computation for timestep " << timestep << " ENDED ***\n\n";
      return;
    }
};

#endif
