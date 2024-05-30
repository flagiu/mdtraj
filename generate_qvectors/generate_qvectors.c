#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define DN 0.5           // mesh size, relative to 2*pi/L
#define MAXQARR 300      // max n. of qvectors with same modulus to be sampled
#define MAXCOUNT 1000    // max n. of extractions before quitting with error

double uniform(double a, double b);
void random_orient(double vec[3]);
double random_radius(double rmin, double rmax);
void random_wavenumbers(double nmin, double nmax, int n[3]);
int check_equivalence(int a[3], int b[3]);

int main(int argc, char** argv){
  int i,j, n, qvecs[MAXQARR][3], count, equivalent;
  double n1,n2;

  if(argc!=2) {
    fprintf(stderr,
	    "\n Usage: %s [n = shell number]\n"
	    "\n Generate random wavenumbers (max %d) nx,ny,nz in units of 2*pi/L"
      "\n within the spherical shell (n1, n2]."
	    "\n - n is an integer in mesh units; mesh is 2*pi/L*(DN=%.2f))"
	    "\n - n1,n2 = DN*(n -+ 1/2)"
      "\n Wavenumbers will be saved into QVECTORS/qvector.n."
      "\n\n", argv[0], MAXQARR, DN);
    exit(1);
  }
  n=atoi(argv[1]);
  n1=DN*((double)n-0.5);
  n2=DN*((double)n+0.5);
  srand48(time(NULL));

  fprintf(stdout,
    "\n n = %d : Random integer wavenumbers with norm in (%.2f,%.2f] ..."
    "\n\n",n,n1,n2
  );

  count=0;
  for(i=0;i<MAXQARR && count<MAXCOUNT;i++){
    count=0;
    do{
      random_wavenumbers(n1,n2,qvecs[i]);
      count++; //count the equivalece extractions at fixed i
      // check if an equivalent wavevector was already extracted
      equivalent=0;
      for(j=0;j<i && !equivalent;j++){
      	equivalent=check_equivalence(qvecs[i], qvecs[j]);
      	if(equivalent) {
          fprintf(stderr, "%d %d %d === %d %d %d\n",
      		  qvecs[i][0], qvecs[i][1], qvecs[i][2],
      			qvecs[j][0], qvecs[j][1], qvecs[j][2]
      		);
        }
      }
    } while(equivalent && count<MAXCOUNT);
    fprintf(stderr, "\n");

  }
  if(count>=MAXCOUNT) {
    i--; //don't print the last one: it's equivalent to some previous one
    fprintf(stdout, "Max n. of equivalence extractions exceeded in main()\n");
  }
  fprintf(stdout,"Terminating with %d wavenumbers\n\n", i);

  FILE *fp;
  char file[50]={0};
  sprintf(file, "QVECTORS/qvector.%03d", n);
  fprintf(stdout, "Saving into %s\n\n", file);
  fp=fopen(file, "w");
  for(j=0;j<i;j++) {
    fprintf(fp, "%d %d %d\n", qvecs[j][0], qvecs[j][1], qvecs[j][2]);
  }
  fclose(fp);

  return 0;
}

/////////////////////////////

double uniform(double a, double b) {
  if(a>b){
    fprintf(stderr, "Error: invalid sampling interval [%f,%f]\n",a,b);
    exit(1);
  }
  return a+drand48()*(b-a);
}

// Marsaglia algorithm
void random_orient(double vec[3]) {
  int i;
  double x[3], xi, norm2;
  do {
    norm2=0.0;
    for(i=0;i<3;i++) {
      x[i]=uniform(-1.,1.); //random in [-1,1]
      norm2+=(x[i]*x[i]);
    }
  } while(norm2>=1.0); //x norm must be <1
  xi=sqrt(1.-norm2);
  vec[0] = 2*x[0]*xi;
  vec[1] = 2*x[1]*xi;
  vec[2] = 1.-2*norm2;
}

double random_radius(double rmin, double rmax) {
  if(rmin>rmax || rmin<0. || rmax<0.){
    fprintf(stderr, "Error: invalid radial interval [%f,%f]\n",rmin,rmax);
    exit(1);
  }
  double umin,umax,u, r;
  umin=rmin*rmin*rmin;
  umax=rmax*rmax*rmax;
  u=uniform(umin,umax);
  r=cbrt(u); // r=u^1/3 gives P(r)~r^2
  return r;
}

void random_wavenumbers(double nmin, double nmax, int n[3]) {
  if(nmin>nmax || nmin<0 || nmax<0){
    fprintf(stderr, "Error: invalid wavenumber interval [%.2f,%.2f)\n",nmin,nmax);
    exit(1);
  }
  int i, m[3], norm2, err, count=0;
  double x[3], r, nmin2=nmin*nmin, nmax2=nmax*nmax;
  do{
    count++;
    random_orient(x); //estrai vettore sulla sfera unitaria
    r=random_radius(nmin,nmax); //estrai raggio ~r^2 in [rmin,rmax]
    norm2=0;
    for(i=0;i<3;i++){
      x[i]*=r; //moltiplica per il raggio
      m[i]=floor(x[i]+0.5); //arrotonda a intero
      norm2+=(m[i]*m[i]);
    }
    err=(norm2<=nmin2 || norm2>nmax2); //ripeti se la norma è fuori dall'intervallo
    fprintf(stderr, "( %d %d %d ) norm %f\n", m[0],m[1],m[2], sqrt(norm2));
  } while(err && count<MAXCOUNT);
  if(count>=MAXCOUNT){
    fprintf(stderr, "Max number of extractions exceeded in random_wavenumbers()\n\n");
    exit(1);
  }
  //termina algoritmo
  for(i=0;i<3;i++) n[i]=m[i];
}

// check if they are the same or the opposite
int check_equivalence(int a[3], int b[3]) {
  int equal=1, opp=1;
  //double anorm=0., bnorm=0.;
  for(int i=0; i<3; i++) {
    equal = (equal && (a[i]==b[i])); //same vector
    opp = (opp && (a[i]==-b[i])); //opposite vector
    //anorm+=(a[i]*a[i]);
    //bnorm+=(b[i]*b[i]);
  }
  //anorm=sqrt(anorm);
  //bnorm=sqrt(bnorm);

  return (equal || opp);
}
