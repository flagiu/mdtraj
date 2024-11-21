#ifndef _SQT_H_
#define _SQT_H_
using namespace std;
//-------- Dynamic Structure Factor S(q,t), a.k.a. Intermediate Scattering Function F(q,t) -----------------------//
//-------- Currently treats all particles as equal -----------------------//
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
    ntype binw, dk;
    int qm, qM, dq, nbins, nTypes, nTypePairs;
    vector<int> Nt;
    int period_in_dt_units, ntimes;
    int num_k_values, exp_storage_size;
    vecflex<ntype> bins, exps, ave_t0,ave2_t0, ave_t0_lin,ave2_t0_lin, ave,ave2;
    vecflex<ntype> num_avg, num_avg_predicted;
    string string_out, myName, tag, rho_tmpfile_prefix;
    fstream fout, f_rho;
    stringstream ss;
    bool debug, verbose, logtime;
    LogTimesteps logt;

    int types2int(int ti, int tj){ // map type pairs (ti,tj) in 0,1,...,nTypes-1 to integer index 0,1,...,nTypePairs
      if (ti>tj) return types2int(tj,ti); // map to ti<=tj
      if(ti==tj) return ti*nTypes - int(ti*(ti-1)/2);
      else return 1 + types2int(ti,tj-1);
    }

    string get_rho_tmpfile(int idx) {
      ss.str(std::string());
      ss<<rho_tmpfile_prefix<<"_idx"<<idx<<tag<<".dat";
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

    void initLogTime()
    {
      int i,j,idx;
      ntimes = logt.get_time_window_size();

      num_avg.resize(ntimes);
      num_avg_predicted.resize(ntimes);
      for(j=0;j<ntimes;j++) {
        num_avg[j]=0.0;
        num_avg_predicted[j]=logt.get_num_avg(j);
      }

      // S(q,t0+t) for each q,t
      ave.resize(nbins * ntimes);
      ave2.resize(nbins * ntimes);

      if(verbose)
      {
        ss.str(std::string()); ss << string_out << tag << ".traj"; fout.open(ss.str(), ios::out);
        fout << "# 1st block: q-wave-vectors; 2nd block: Delta timestep; Other blocks: S(q,t;t0) trajectory for each frame\n";
        for(i=0; i<nbins; i++) fout << bins[i] << endl;
        fout<<endl; // new block
        for(i=0;i<ntimes;i++)
          fout << logt.get_dt(i) << endl; // Delta timestep values
        fout.close();
      }

      // reset rho, <S(q,dt)>, <S(q,dt)^2> for each q=qmin,qmin+dq,...,qmax
      //    and for each dt=0,1,...,period_in_dt_units-1
      for(i=0; i<nbins; i++){
        for(j=0;j<ntimes;j++){
          idx = ntimes*i + j;
          ave[idx]=ave2[idx]=0.0;
        }
      }
      ss.str(std::string()); ss << string_out << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "# 1st block: q-wave-vector; 2nd block: Delta timestep; 3rd block: <S(q,t)> matrix; 4th block: <S(q,t)> error matrix; 5th block: <chi4(q,t)> suscettivity matrix.\n";
      for(i=0; i<nbins; i++) fout << bins[i] << endl;
      fout<<endl; // new block
      for(i=0;i<ntimes;i++) fout << logt.get_dt(i) << endl;
      fout.close();
    }

    void initLinearTime(int period_in_real_units, int dtframe, int nframes)
    {
      int i,j,idx;
      if(period_in_real_units>0 && period_in_real_units<nframes*dtframe){
        period_in_dt_units = ceil(period_in_real_units/dtframe);
      } else{
        period_in_dt_units = nframes;
      }

      ntimes = period_in_dt_units;
      num_avg.resize(ntimes);
      num_avg_predicted.resize(ntimes);
      i=nframes/period_in_dt_units;
      for(j=0;j<ntimes;j++) {
        num_avg[j]=0.0;
        num_avg_predicted[j]=i;
      }

      // S(q,t0+t) for each q,t
      ave.resize(nbins * ntimes);
      ave2.resize(nbins * ntimes);

      if(verbose)
      {
        ss.str(std::string()); ss << string_out << tag << ".traj"; fout.open(ss.str(), ios::out);
        fout << "# 1st block: q-wave-vectors; 2nd block: Delta timestep; Other blocks: S(q,t;t0) trajectory for each frame\n";
        for(i=0; i<nbins; i++) fout << bins[i] << endl;
        fout<<endl; // new block
        for(i=0;i<ntimes;i++)
          fout << i*dtframe << endl; // Delta timestep values
        fout.close();
      }

      // reset rho, <S(q,dt)>, <S(q,dt)^2> for each q=qmin,qmin+dq,...,qmax
      //    and for each dt=0,1,...,period_in_dt_units-1
      for(i=0; i<nbins; i++){
        for(j=0;j<ntimes;j++){
          idx = ntimes*i + j;
          ave[idx]=ave2[idx]=0.0;
        }
      }
      ss.str(std::string()); ss << string_out << tag << ".ave"; fout.open(ss.str(), ios::out);
      fout << "# 1st block: q-wave-vector; 2nd block: Delta timestep; 3rd block: <S(q,t)> matrix; 4th block: <S(q,t)> error matrix.\n";
      for(i=0; i<nbins; i++) fout << bins[i] << endl;
      fout<<endl; // new block
      for(i=0;i<ntimes;i++) fout << i*dtframe << endl;
      fout.close();
    }

    void init(int q_mod_min, int q_mod_max, int q_mod_step, int dtframe,
      int nframes, int period_in_real_units, vec box_diagonal,
      int nTypes_,int* Nt_, bool logtime_, LogTimesteps logt_,
      string string_out_, string tag_, bool debug_, bool verbose_)
    {
      int i,j,idx;
      string_out = string_out_;
      tag = tag_;
      debug = debug_;
      verbose = verbose_;
      logtime=logtime_;
      logt = logt_;
      nTypes = nTypes_;
      nTypePairs = nTypes*(nTypes+1)/2;
      Nt.resize(nTypes);
      int N=0;
      for(i=0;i<nTypes;i++) {
        Nt[i]=Nt_[i];
        N+=Nt[i];
      }
      qm = q_mod_min;
      qM = q_mod_max;
      dq = q_mod_step;

      nbins = int(floor( (qM-qm)/dq )) + 1;
      binw = M_PI / box_diagonal[0]; // half mesh: it is in units of pi/L
      bins.resize(nbins);
      for(i=0; i<nbins; i++) bins[i] = (qm+i*dq)*binw; // q-wave-vector values

      dk=2.*binw; // 2*M_PI/L[0];
      num_k_values=1+qM/2; // 1/2 is the mesh spacing!
      exp_storage_size=num_k_values*N*3*2; // wavenumber, particle, cartesian, real/imag
      exps.resize(exp_storage_size);
      ave_t0.resize(nbins);
      ave2_t0.resize(nbins);

      if(logtime){
        initLogTime();
      } else {
        initLinearTime(period_in_real_units, dtframe, nframes);
      }

      if(debug) cout << myName << " Initialization COMPLETED\n";
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
      ntype x0,x1, y0,y1, z0,z1;

      for(k=0;k<nbins;k++) { ave_t0[k]=ave2_t0[k]=0.0; }

      if(debug) cout << "*** "<<myName<<" computation for timestep " << timestep << " STARTED ***\n";

      if(logtime){
        nperiod = frameidx / logt.npc;
        dframe = frameidx % logt.npc; // log-spaced at short time
        if(dframe==0 && nperiod>0) dframe = logt.npc-1+nperiod; // linearly-spaced at long time
        /*
        frameidx 0 1 ... npc-1 npc npc+1 ... 2*npc-1 2*npc 2*npc+1 ...
        nperiod  0 0 ...     0  1    1   ...    1       2     2    ...
        dframe   0 1 ... npc-1 npc!! 1   ...  npc-1  npc+1!!  1    ...

        dframe   0 1 ... npc-1 npc npc+1 npc+2 ... npc-1+ncycles
        frameidx 0 1 ... npc-1 npc 2*npc 3*npc ... (ncycles-1+1)*npc
        nperiod  0 0 ...   0    1    2     3   ...     ncycles
        */
        if(dframe>=logt.npc){
          // Average over orientations at fixed q,t0, for each linear t
          ave_t0_lin.resize(nbins*nperiod);
          ave2_t0_lin.resize(nbins*nperiod);
          for(k=0;k<nbins*nperiod;k++) ave_t0_lin[k]=ave2_t0_lin[k]=0.0;
        }
      } else {
        nperiod = frameidx / ntimes;
        dframe = frameidx % ntimes;
      }
      if(debug) cout<<"  frameidx="<<frameidx<<" ==> nperiod="<<nperiod<<" , dframe="<<dframe<<endl;


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
        idx  = ntimes*k + dframe; // rho(q,t0+t)
        qmod = qm + k*dq; // in units of half mesh binw=pi/L
        if(debug) cout << "  k="<<k<<", qmod="<<qmod<<"; dframe="<<dframe<<", idx="<<idx<<endl;
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
          if(debug) cout << "   read kx,ky,kz = " << qx<<", "<< qy<<", "<< qz<<endl;
          qvectors[3*num_qvectors  ]=qx;
          qvectors[3*num_qvectors+1]=qy;
          qvectors[3*num_qvectors+2]=qz;
          num_qvectors++;
        }
        fout.close();
        if(debug) cout << "  Total "<<num_qvectors<<" triplets of wavenumbers.\n";

        for(q=0;q<num_qvectors;q++) {
          rhoQvecRe=rhoQvecIm=0.0;
          qx=qvectors[3*q  ];
          qy=qvectors[3*q+1];
          qz=qvectors[3*q+2];
          if(debug) cout << "  calculating on kx,ky,kz = " << qx<<", "<< qy<<", "<< qz<<endl;
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
            // Compute rho(qvec,t0+t)
            rhoQvecRe += ( x0*y0*z0 - x0*y1*z1 - x1*y0*z1 - x1*y1*z0 ); // cos( k*r )
            rhoQvecIm += ( x1*y0*z0 + x0*y1*z0 + x0*y0*z1 - x1*y1*z1 ); // sin( k*r )
          } // end of atoms
          rhoQvecRe*=invSqrtN;
          rhoQvecIm*=invSqrtN;
          // Save rho(qvec,t0+t) for each qvec in a file; for later projection
          if(!logtime || (dframe==0||dframe>=logt.npc))
          { // no need to save frames within the npc log cycle 0<dframe<npc
            if(q==0) write_rho_tmpfile(idx,rhoQvecRe,rhoQvecIm);
            else    append_rho_tmpfile(idx,rhoQvecRe,rhoQvecIm);
          }

          // Project rho(qvec,t0+t) onto rho(qvec,t0):
          // - linear sampling: t0 is the beginning of each cycle
          // - linlog sampling: depends
          if(logtime){
            if(dframe==0)    idx_old=ntimes*k+0; // if first absolute snapshot: t0 is itself
            else if(dframe<logt.npc){
              if(nperiod==0) idx_old=ntimes*k+0; // else if first period, t0 is the begining of the first npc log cycle
              else           idx_old=ntimes*k+logt.npc-1+nperiod; // else t0 is the begining of the current npc log cycle
            }
            else             idx_old=ntimes*k+0; // for linear subsampling, first compare with first absolute snapshot
          } else idx_old=ntimes*k+0;
          read_rho_tmpfile(q,idx_old,&rhoQvecRe_old,&rhoQvecIm_old);
          sqtQvecRe =  rhoQvecRe*rhoQvecRe_old + rhoQvecIm*rhoQvecIm_old;
          sqtQvecIm = -rhoQvecRe*rhoQvecIm_old + rhoQvecIm*rhoQvecRe_old;
          if(debug) {
            cout << "    rho(qxyz;    idx="<<idx<<") = "<<rhoQvecRe<<" + i* "<<rhoQvecIm<<endl;
            cout << "    rho(qxyz;idx_old="<<idx_old<<") = "<<rhoQvecRe_old<<" + i* "<<rhoQvecIm_old<<endl;
            cout << "   S(qx,qy,qz;dframe="<<dframe<<") = "<<sqtQvecRe<<" + i* "<<sqtQvecIm<<endl;
          }

          // Average over orientations at fixed q,t0,t
          ave_t0[k] += sqtQvecRe;
          ave2_t0[k] += sqtQvecRe*sqtQvecRe;
          //ave2_t0[k] += sqtQvecRe*sqtQvecRe + sqtQvecIm*sqtQvecIm;
          //ave_t0[k] += sqrt(sqtQvecRe*sqtQvecRe + sqtQvecIm*sqtQvecIm);

          // for log-lin sampling, when outside the log cycle,
          //   project also onto any linear sub-intervals:
          // dframe  : 0 1 ... npc-1 npc npc+1 npc+2 ... npc-1+ncycles
          // frameidx: 0 1 ... npc-1 npc 2*npc 3*npc ... (ncycles-1+1)*npc
          // nperiod : 0 0 ...   0    1    2     3   ...     ncycles
          // e.g. project dframe=npc+2 (nperiod=3), other than on dframe=0,
          //      on dframe=npc+1,npc (nperiod=2,1)
          if(logtime && dframe>=logt.npc) {
            if(debug) cout << "   Out-of-log-cycle linear subsampling:\n";
            for(i=0;i<nperiod;i++) {
              idx_old = idx-i;
              read_rho_tmpfile(q,idx_old,&rhoQvecRe_old,&rhoQvecIm_old);
              if(debug) {
                cout << "    rho(qxyz;    idx="<<idx<<") = "<<rhoQvecRe<<" + i* "<<rhoQvecIm<<endl;
                cout << "    rho(qxyz;idx_old="<<idx_old<<") = "<<rhoQvecRe_old<<" + i* "<<rhoQvecIm_old<<endl;
              }
              sqtQvecRe = rhoQvecRe*rhoQvecRe_old + rhoQvecIm*rhoQvecIm_old;
              sqtQvecIm = -rhoQvecRe*rhoQvecIm_old + rhoQvecIm*rhoQvecRe_old;
              if(debug) cout << "   S(qxyz;dt="<<i<<"*"<<logt.npc<<") = "<<sqtQvecRe<<" + i* "<<sqtQvecIm<<endl;
              // Average over orientations at fixed q,t0, for each linear t
              ave_t0_lin[k+nbins*i]  += sqtQvecRe;
              ave2_t0_lin[k+nbins*i] += sqtQvecRe*sqtQvecRe;
              //ave2_t0_lin[k+nbins*i] += sqtQvecRe*sqtQvecRe + sqtQvecIm*sqtQvecIm;
              //ave_t0_lin[k+nbins*i] += sqrt(sqtQvecRe*sqtQvecRe + sqtQvecIm*sqtQvecIm);
            }
          }
        } // end of qvectors
        if(debug) cout << "  done calculating qvectors\n";

        ave_t0[k] /= num_qvectors;
        ave2_t0[k] /= num_qvectors;
        if(debug){
          if(num_qvectors>1){ // fluctuations of S(q;t0) over q orientation
            cout << "   S(k="<<k<<";dframe="<<dframe<<") = "<<ave_t0[k]<<" +- "<<sqrt( (ave2_t0[k]-ave_t0[k]*ave_t0[k])/(num_qvectors-1) )<<endl;
          } else {
            cout << "   S(k="<<k<<";dframe="<<dframe<<") = "<<ave_t0[k]<<" +- 0.0\n";
          }
        }

        idx = ntimes*k + dframe;
        ave[idx] += ave_t0[k]; // average over t0 at fixed q,t
        ave2[idx] += ave_t0[k]*ave_t0[k];
        num_avg[dframe]+=1./(ntype)nbins;

        if(logtime && dframe>=logt.npc) {
          if(debug) cout << "  Out-of-log-cycle linear subsampling:\n";

          for(i=0;i<nperiod;i++) {
            j=k+nbins*i;
            ave_t0_lin[j]  /= num_qvectors;
            ave2_t0_lin[j] /= num_qvectors;
            if(debug){
              if(num_qvectors>1){ // fluctuations of S(q;t0) over q orientation
                cout << "   S(k="<<k<<";dt="<<i<<"*"<<logt.npc<<") = "<<ave_t0_lin[j]<<" +- "<<sqrt( (ave2_t0_lin[j]-ave_t0_lin[j]*ave_t0_lin[j])/(num_qvectors-1) )<<endl;
              } else {
                cout << "   S(k="<<k<<";dt="<<i<<"*"<<logt.npc<<") = "<<ave_t0_lin[j]<<" +- 0.0\n";
              }
            }
            if(i==0) {
              idx = ntimes*k + 0; // dt for the average among linear subsamples
              num_avg[0]+=1/(ntype)nbins;
            }else{
              idx = ntimes*k + logt.npc+i-1;
              num_avg[logt.npc+i-1]+=1/(ntype)nbins;
            }
            ave[idx] += ave_t0_lin[j]; // average over t0 at fixed q,t within the linear sub-tajectories
            ave2[idx] += ave_t0_lin[j]*ave_t0_lin[j];
          }
        }
      } // end of qvector modulus (index k)

      if(verbose)
      {
        ss.str(std::string()); ss << string_out << tag << ".traj"; fout.open(ss.str(), ios::app);
        fout << endl; // start new block
        for(k=0; k<nbins; k++) fout << ave_t0[k] << endl;
        fout.close();
      }

      if(frameidx == (nframes-1)) //last frame
      {
          if(debug){
            cout<<"Please check that num_avg is the same as the predicted one:\n";
            cout<<" index num_avg predicted\n";
            for(j=0;j<ntimes;j++)
              cout << " "<<j<<" "<<num_avg[j]<<" "<<num_avg_predicted[j]<<endl;
          }
          ss.str(std::string()); ss << string_out << tag << ".ave"; fout.open(ss.str(), ios::app);
          fout<<endl; // new block
          for(k=0;k<nbins;k++){
            for(j=0;j<ntimes;j++){
              idx = ntimes*k + j;
              ave[idx]/=num_avg[j]; // average over t0 at fixed q,t
              fout << ave[idx] << " ";
            }
            fout << endl;
          }

          fout << endl; // new block
          for(k=0;k<nbins;k++){
            for(j=0;j<ntimes;j++){
              idx = ntimes*k + j;
              ave2[idx]/=num_avg[j];
              if(num_avg[j]>1){
                fout << sqrt( (ave2[idx]-ave[idx]*ave[idx])/(num_avg[j]-1) ) << " ";
              }else{
                fout << "0.0 ";
              }
            }
            fout << endl;
          }
          /* // Susceptibility chi4
          fout << endl; // new block
          for(k=0;k<nbins;k++){
            for(j=0;j<ntimes;j++){
              idx = ntimes*k + j;
              //ave2[idx]/=num_avg[j];
              fout << ave2[idx]-ave[idx]*ave[idx] << " ";
            }
            fout << endl;
          }
          */
          fout.close();
          if(debug) cout << "Average "<<myName<<" printed to file\n";

          if(!debug){ // clear temporary files
            for(k=0; k<nbins; k++){
              for(j=0;j<ntimes;j++){
                idx = ntimes*k + j;
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
