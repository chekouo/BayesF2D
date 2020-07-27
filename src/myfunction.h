#include <gsl/gsl_rng.h>
#include <gsl/gsl_blas.h>
#include <stdlib.h>
#include <stdbool.h>
#define MIN(x, y) (((x) < (y)) ? (x) : (y))


 void MCMC(int* typeOutcome,int * n1,int *TT1,int *r1,int *nbrcov1,double * resp, double *X,int *KK, double* Base1, double* Xcov1, int *nbrsample1,int *burninsample1, double * Betasample,double *XsiMean1,double* rhoMean1,double * CovImageSample,double *LambdaMean1,double* seed1,double *hypsigm_rr,double *h1,int *covbool,int *rhosample1, double *yMean,double *sigma2xMean,double * sigma2yMean1,int * SamplingRho);


void Prediction(int* typeOutcome,int * n1,int* ntest1,int *TT1,int *r1,int *nbrcov1,double *yMean ,double *resp,double *X,double *Xtest, int *KK, double* Base1, double* Xcov1,double* Xcovtest1,int* Kf, int *nbrset,int * ListTest, int * ListTraining,int *nbrsample1,int* Predict1,double *XsiMean1,double *prob,double * CovImageMean1,int * rhosample1,double *LambdaMean1,double* seed1,double *hypsigm_rr,double *h1,int *covbool,double *sigma2xMean,double *sigma2yMean1);



void readDouble(char *filename, int nRows,  double * data );
void readInt(char *filename, int nRows,  int * data );
void SampleCovImage(gsl_rng * rr,int r,int *KK, double **CovImage,double ** Beta,double h, _Bool ** rho,double alpha0, double beta0);
double max(int n, double * x);
_Bool compar(int n,_Bool *u,_Bool *v);
double mean(int n,double * x);
void GetH (int l,int r, _Bool **rho,double *scal1,double *Muscal1,double **Beta, double h1,double **   Sigma, int *KK); 

//void MCMC(int *typeOutcome,int * n1,int* ntest1,int *TT1,int *r1,int *nbrcov1,double * resp, double *X,double *Xtest, int *KK,
   //        double* Base1, double* Xcov1,double* Xcovtest1,int* Kf, int *nbrset,int * ListTest, int * ListTraining, int *nbrsample1,  int *burninsample1,int* Predict1,
    //       double * Betasample,double *XsiMean1,double* rhoMean1,double *logpost,double *prob,double *  CovImageSample);
  static inline double r_lefttruncnorm(double a, double mean, double sd,const gsl_rng * r);

  void predictCrossVal(_Bool typeoutcome1,int PartNbr,int *nbrset,int ** Set, int **SetCompl,gsl_rng * rr,int n,int TT, int r, int *KK, int nbrsample, int P, int nbrcov, _Bool ** rhosample,   double *respy,double* y, double **Xcov, double ***X2D,double *** Base,double *** XsiMean, double * sigma2xMean,double ** LambdaMean,double ** Beta, double sigma2y,         double** CovImage,double h, double h0, double a0, double b0, int maxmodel, double *prob); 
  void predictNewdata(_Bool typeoutcome1,gsl_rng * rr,int n, int ntest,int TT, int r, int *KK, int nbrsample, int P, int nbrcov, _Bool ** rhosample, double *respy,double ***X2Dtest, double ** Xcov, double **Xcovtest,double *** Base,double *** XsiMean, double * sigma2xMean,double ** LambdaMean, double sigma2y,double **CovImage, double h, double    h0, double a0, double b0, int maxmodel, double *prob); 
  void SampleBetaCova(gsl_rng * rr,int n,int p, double ** Xcov, double h, double sigma2y,double *y, double *Beta);
double r_righttruncnorm(double b, double mean, double sd,const gsl_rng * r);
double ers_a_inf(double a,const gsl_rng * r);
double nrs_a_inf(double a,const gsl_rng * r);
double LogPost(int n, int r, int TT, int nbrcov, int *KK, double *** Base, double ***X2D, double **Xcov, 
               double *respy, double ***Xsi, double **Beta,double *BetaCova,_Bool ** rho,
               double sigma2y, double* sigma2x, double *w, double h,double h0,
               double** Lambda, double a0, double b0);
void free_bmatrix(_Bool **m, int nrl, int nrh, int ncl, int nch);
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch);
 void SampleRegressionCoeff(const gsl_rng * r,int n, int p, double *y, double ** X, double *beta,double *MuBetaPrior, double *Mu,
                            double* sigma2, double *SigmaPriorBeta);
int multivariate_gaussian (const gsl_rng * r,
                           const gsl_vector * mu,
                           const gsl_matrix * L,
                           gsl_vector * result);
double **dmatrix(int nrl, int nrh, int ncl, int nch); 
_Bool **bmatrix(int nrl, int nrh, int ncl, int nch);
void nrerror(char error_text[]);
void findc(int n,_Bool R[n],int a,int * IDX, int *nx);
int logmultigaussian(const gsl_vector * Ax,
                     const gsl_vector * x,
                     const gsl_matrix * L,
                     double sigma2,double logDetD,
                     double * result,double *quadF,
                     gsl_vector * work);
void proposal(int n,_Bool R[n],_Bool R1[n],float phi, gsl_rng * r);
void logLikelihood(int n,int p,double *y,double ** X,double sigma2,double Sigmall,double scal,double Muscal, double* quadForm,_Bool * Gam,           double * loggauss);
void SampleRho1(gsl_rng * rr,int n,int p,double *y,double ** X,double sigma2,double h,double scal,double Muscal,_Bool * Gam,double w);
void SampleRho(int l,int r,gsl_rng * rr,int n,int p,double *y,double ** X,double sigma2,double** CovImage,   double *scal,double *Muscal,_Bool ** rho,int *K,double w,double **Beta,double h);
void SampleSig2X(gsl_rng * rr,int n,int K, int TT,double a0,  double b0,double * sigma2, double **Base, double **resp, double **Xsi);
void SampleSig2(gsl_rng * rr,int n,int * KK, int r,double a0,  double b0,double * sigma2, double ***Xsi, double *resp, 
                double **Beta,double h,_Bool ** rho);
void SampleLambda(gsl_rng * rr,int n,int k,double a0, double b0,double * lambda, double **Xsi);
double Truncate(double mu,double sd, double lower, const gsl_rng * r);
void Sampley(gsl_rng * rr,int n,int * KK, int r, double ***Xsi, double *y, double *resp,
             double **Beta,double ** Xcov,double *       BetaCov, int nbrcov);
void SampleBeta(gsl_rng * rr,int n,int p,_Bool *rho, double ** Xsi, double h, double sigma2y,double *y, double *Beta,double scal,double Muscal);
void SampleXsi(gsl_rng * rr,_Bool includeOutcome,int n, int TT, int p, double ** Base, double *Beta,double **X2D, double *y, 
               double sigma2y, double sigma2x, double **Xsi,double **Mu, double* Lambda);
double scalar(int n,double A[n*n],double b_data[n], double *det);
double gsl_linalg_CD_lndet (gsl_matrix * CD);
_Bool ** UniqueModel(int nbrsample, int p, _Bool ** rhosample,int * modelidx,int *countd1);
void EstimateBeta(gsl_rng * rr,int n, int nbrcov, int P,
                   double *respy, double **XsiP, double *BetaP,_Bool * rhoP,
                   double sigma2y, double*InvSigmall,double h0);
void Combinefeatures(int n, int r, int nbrcov, int * KK, int P, double **Xcov, 
                     double ***Xsi,double ** XsiP);
double logpostRho1(int n, int nbrcov, int P,  
                  double *respy, double **XsiP, _Bool * rho,
                  double h,double h0, double a0, double b0);
 double logpostRho(int n, int nbrcov, int P,double sigma2y,double*InvSigmall,double logdetSigmall,
                    double *respy, double **XsiP, _Bool * rho,double h0, double a0, double b0);
void sort(int n,double *x,int *idx);
 double *GetInvSigmaBeta(_Bool *rhop,int r,int *KK, double **CovImage,double h, double *logdet);
