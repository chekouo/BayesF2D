#include <float.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <time.h>
#include <gsl/gsl_randist.h>
 #include <gsl/gsl_sf_erf.h>
#include <math.h>
 #include <gsl/gsl_cdf.h>
# include "myfunction.h"



void Prediction(int* typeOutcome,int * n1,int* ntest1,int *TT1,int *r1,int *nbrcov1,double * yMean,double *resp,double *X,double *Xtest, int *KK, double* Base1, double* Xcov1,double* Xcovtest1,int* Kf, int *nbrset,int * ListTest, int * ListTraining,int *nbrsample1,int* Predict1,double *XsiMean1,double *prob,double * CovImageMean1,int * rhosample1,double *LambdaMean1,double* seed1,double *hypsigm_rr,double *h1,int *covbool,double *sigma2xMean,double *sigma2yMean1){
_Bool typeOutcome1=typeOutcome[0];
int t,i,j,s,l,l1,k;int bs=0;
long seed=seed1[0];
 gsl_rng * rr = gsl_rng_alloc (gsl_rng_rand48);
 gsl_rng_set (rr, seed);
 int  Predict=Predict1[0];
int PartNbr=Kf[0];
int n=n1[0];
int ntest=ntest1[0];
/*
  if (Predict==1){
      printf("Cross-Validation\n");
         printf("Number of k-folds from cross-val is %d\n",Kf[0]);
  } else {
 printf("We have provided a new data set for test data of size %d\n",ntest);
  }
*/
int TT=TT1[0];
 int r=r1[0];//Number of imaging sequences;
int nbrsample=nbrsample1[0];
int nbrcov=nbrcov1[0];
int P=nbrcov;
 for (l=0;l<r;l++)
 for (k=0;k<KK[l];k++) P++;
_Bool ** rhosample=bmatrix(0,nbrsample-1,0,P-nbrcov-1);
bs=0;int bs1=0;
for (s=0;s<nbrsample;s++){
bs1=0;
for (l=0;l<r;l++){
  for (k=0;k<KK[l];k++){
   rhosample[s][bs1]=rhosample1[bs];
   bs++;bs1++;
}}}

double *** X2Dtest=malloc(r*sizeof(double **));
double ** Xcovtest=malloc(ntest*sizeof(double*));
  int j1=0;
   for (i=0;i<ntest;i++){
     Xcovtest[i]=malloc(nbrcov*sizeof(double));
     for (j=0;j<nbrcov;j++){
       Xcovtest[i][j]=Xcovtest1[j1];
       j1++;
     }
   }
    k=0;
    for (l=0;l<r;l++){
      X2Dtest[l]=dmatrix(0,ntest-1,0,TT-1);
      for (i=0;i<ntest;i++){
      for (t=0;t<TT;t++){
          X2Dtest[l][i][t]=Xtest[k];k++;
        }}
    }

j1=0;
double ** Xcov=malloc(n*sizeof(double*));
for (i=0;i<n;i++){
  Xcov[i]=malloc(nbrcov*sizeof(double));
  for (j=0;j<nbrcov;j++){
    Xcov[i][j]=Xcov1[j1];
    j1++;
  }
}

double *** Base=malloc(r*sizeof(double**));
 int k1=0;
 for (l=0;l<r;l++){
Base[l]=dmatrix(0,TT-1,0,KK[l]-1);
  for (t=0;t<TT;t++){
    for (k=0;k<KK[l];k++){
  Base[l][t][k]=Base1[k1];
   k1++;
 }
 }
 }
int rx=0;
 double *** XsiMean=malloc(r*sizeof(double **));
 for (l=0;l<r;l++){
  XsiMean[l]=dmatrix(0,n-1,0,KK[l]-1);
  for (k=0;k<KK[l];k++){
   for (i=0;i<n;i++){
      XsiMean[l][i][k]=XsiMean1[rx];
      rx++;
     }}
 }
double ** LambdaMean=malloc(r*sizeof(double*)); //variance for latent variables \xs
bs=0;
for (l=0;l<r;l++){
LambdaMean[l]=malloc(KK[l]*sizeof(double));
for (k=0;k<KK[l];k++){
LambdaMean[l][k]=LambdaMean1[bs];
 bs++;
}
}
double sigma2yMean=sigma2yMean1[0];

 double ** CovImageMean=malloc(r*sizeof(double*));
rx=0;
for (l=0;l<r;l++){
CovImageMean[l]=malloc(r*sizeof(double));
        for (l1=0;l1<r;l1++){
        CovImageMean[l][l1]=CovImageMean1[rx];
rx++;        
}
}
double h=h1[0]; //Shinkage hyper of regression effects
double h0=1000;
double a0=0.001; double b0=0.001;
int maxmodel=100;
if (Predict==0){
predictNewdata(typeOutcome1,rr, n, ntest,TT, r, KK, nbrsample,P, nbrcov, rhosample,yMean,X2Dtest, Xcov,Xcovtest, Base, XsiMean, sigma2xMean, LambdaMean, sigma2yMean, CovImageMean, h, h0,  a0,  b0, maxmodel, prob);
} else if (Predict==1){
int pp;int p1=0;int p2=0;
int ** Set=malloc(PartNbr*sizeof(int*));
 int ** SetCompl=malloc(PartNbr*sizeof(int*));
for  (pp=0;pp<PartNbr;pp++){
    Set[pp]=malloc(nbrset[pp]*sizeof(int));
    SetCompl[pp]=malloc((n-nbrset[pp])*sizeof(int));
    for (i=0;i<nbrset[pp];i++){
        Set[pp][i]=ListTest[p1]-1;p1++;
    //printf("%d ",Set[pp][i]);
    }
     for (i=0;i<n-nbrset[pp];i++){
         SetCompl[pp][i]=ListTraining[p2]-1;p2++;
     }
    //printf("\n");
}

double *** X2D=malloc(r*sizeof(double **));
k=0;
 for (l=0;l<r;l++){
   X2D[l]=dmatrix(0,n-1,0,TT-1);
   for (i=0;i<n;i++){
   for (t=0;t<TT;t++){
       X2D[l][i][t]=X[k];k++;
     }}
 }
double ** Beta=malloc(r*sizeof(double*));
 for (l=0;l<r;l++){
Beta[l]=malloc(KK[l]*sizeof(double));
for (k=0;k<KK[l];k++) Beta[l][k]=1;
}

predictCrossVal(typeOutcome1,PartNbr,nbrset,Set, SetCompl,rr,n,TT, r, KK,nbrsample,P, nbrcov, rhosample, yMean, resp, Xcov, X2D,Base,XsiMean,sigma2xMean,LambdaMean,Beta, sigma2yMean,CovImageMean, h, h0, a0, b0, maxmodel,prob);
 for  (pp=0;pp<PartNbr;pp++){
free(Set[pp]);free(SetCompl[pp]);
}
for (l=0;l<r;l++){
free_dmatrix(X2D[l],0,n-1,0,TT-1);
}
free(Set);free(SetCompl);free(X2D);
}

free_dmatrix(Xcov,0,n-1,0,nbrcov-1);free_dmatrix(Xcovtest,0,ntest-1,0,nbrcov-1);
free_bmatrix(rhosample,0,nbrsample-1,0,P-nbrcov-1);
for (l=0;l<r;l++){
free(LambdaMean[l]); 
free_dmatrix(Base[l],0,TT-1,0,KK[l]-1);
free_dmatrix(X2Dtest[l],0,ntest-1,0,TT-1);
free_dmatrix(XsiMean[l],0,n-1,0,KK[l]-1);
}
free(LambdaMean);   free(XsiMean);free(X2Dtest);free(Base);
free_dmatrix(CovImageMean,0,r-1,0,r-1);
}



void predictCrossVal(_Bool typeoutcome1,int PartNbr,int *nbrset,int ** Set, int **SetCompl,gsl_rng * rr,int n,int TT, int r, int *KK, int nbrsample, int P, int nbrcov, _Bool ** rhosample, double *respy,double* y, double **Xcov, double ***X2D,double *** Base,double *** XsiMean, double * sigma2xMean,double ** LambdaMean,double ** Beta, double sigma2y, double** CovImage,double h, double h0, double a0, double b0, int maxmodel, double *prob){

    int t,j,l,s,i,k;
    int countmod;  int *modelidx=malloc(nbrsample*sizeof(int));
    _Bool ** rhosampleUniq=UniqueModel(nbrsample, P-nbrcov, rhosample,modelidx,&countmod);
    double ** XsiP=dmatrix(0,n-1,0,P-1);
    Combinefeatures(n, r, nbrcov, KK,P, Xcov,XsiMean,XsiP);
int nbrpart=0;
for (nbrpart=0;nbrpart<PartNbr;nbrpart++){
int ntr=n-nbrset[nbrpart];
    double ** XsiPt=dmatrix(0,ntr-1,0,P-1);
    double ** XsiPte=dmatrix(0,n-ntr-1,0,P-1); 
    double * respy1=malloc(ntr*sizeof(double));
 
     for (i=0;i<ntr;i++){
        respy1[i]=respy[SetCompl[nbrpart][i]];
        for (j=0;j<P;j++){
            XsiPt[i][j]=XsiP[SetCompl[nbrpart][i]][j];
    }
    }

      for (i=0;i<n-ntr;i++){
         for (j=0;j<P;j++){
             XsiPte[i][j]=XsiP[Set[nbrpart][i]][j];
     }
     }

  //  if (strcmp(Predict,"CV")==1){
 double ** BetaS=dmatrix(0,countmod-1,0,P-1);
 double * Weight=malloc(countmod*sizeof(double));
 double * lpostRho=malloc(countmod*sizeof(double));
 //for (l=0;l<r;l++) for (s=0;s<l;s++) CovImage[s][l]=CovImage[l][s]=0;
 for (s=0;s<countmod;s++){
     double logdetSigma=0;
double *InvSigmall= GetInvSigmaBeta(rhosampleUniq[s],r,KK,CovImage,h,&logdetSigma); 
EstimateBeta(rr,ntr, nbrcov,P,respy1, XsiPt, BetaS[s],rhosampleUniq[s],sigma2y, InvSigmall, h0);
 Weight[s]=0;
 for (i=0;i<n-ntr;i++){
     double xb=0;
        for (j=0;j<P;j++){
      xb+=XsiPte[i][j]*BetaS[s][j];
        }
if (typeoutcome1==1){//binary response
        double logcdf=-log(2)+gsl_sf_log_erfc(-xb/sqrt(2));
        double logOneMinuscdf=-log(2)+gsl_sf_log_erfc(xb/sqrt(2));
        //double logcdf=log(gsl_cdf_ugaussian_P(xb));double logOneMinuscdf=log(gsl_cdf_ugaussian_Q(xb));

     Weight[s]+=-(respy[Set[nbrpart][i]]>0)*logcdf-(1-(respy[Set[nbrpart][i]]>0))*logOneMinuscdf;
 } else if (typeoutcome1==0){
Weight[s]+=-0.5*log(sigma2y)-pow(respy[Set[nbrpart][i]]-xb,2)/(2*sigma2y);
}

 }

 lpostRho[s]=logpostRho(n, nbrcov, P,sigma2y,InvSigmall,logdetSigma,respy, XsiP, rhosampleUniq[s],h0, a0, b0);
free(InvSigmall); 
 }
//int maxmodel=1000;
  int * highmodelidx=malloc(countmod*sizeof(int));
  sort(countmod,lpostRho,highmodelidx);
   
  double MaxWeight=max(countmod,Weight);
  int nbrmax=MIN(maxmodel,countmod);
   double * Weight1=malloc(nbrmax*sizeof(double));
  for (l=0;l<nbrmax;l++){
    Weight1[l]=exp(Weight[highmodelidx[l]]-MaxWeight);
   }
   double sumweight=nbrmax*mean(nbrmax,Weight1);

   for (l=0;l<nbrmax;l++){
     Weight1[l]=Weight1[l]/sumweight;
   }

   int ntest=nbrset[nbrpart];
  double *** Xsitest=malloc(r*sizeof(double **));
      for (l=0;l<r;l++){
       Xsitest[l]=dmatrix(0,ntest-1,0,KK[l]-1);
       for (i=0;i<ntest;i++){
         for (k=0;k<KK[l];k++){
           Xsitest[l][i][k]=0;
         }}}

       //Predict scores
 for (l=0;l<r;l++){
     double *yt=malloc(ntest*sizeof(double));
     double **X2Dtest1=dmatrix(0,ntest-1,0,TT-1);
     for (i=0;i<ntest;i++){ 
         yt[i]=0;
         for (t=0;t<TT;t++){
     X2Dtest1[i][t]=X2D[l][Set[nbrpart][i]][t]; 
 //printf("XXB=%lf ",X2Dtest1[i][t]);     
}
     }
  SampleXsi(rr,0,ntest, TT,KK[l], Base[l],Beta[l], X2Dtest1, yt,1, sigma2xMean[l], Xsitest[l],Xsitest[l],          LambdaMean[l]);
 free(yt);free_dmatrix(X2Dtest1,0,ntest-1,0,TT-1);
 }

 //Prediction of new samples using probit models
 //double ** XsiPtest=dmatrix(0,ntest-1,0,P-1);
 // Combinefeatures(ntest, r, nbrcov, K,P, Xcovtest,Xsitest,XsiPtest);

 int nm;
// printf("\nEstimated probabilities\n");
for (i=0;i<ntest;i++){
 prob[Set[nbrpart][i]]=0;
double * LogProb=malloc(nbrmax*sizeof(double));
double probb=0; 
// double sumw=0;    
 for (nm=0;nm<nbrmax;nm++){
LogProb[nm]=-DBL_MAX;
     double xb=0;
     for (j=0;j<nbrcov;j++){
     xb+=Xcov[Set[nbrpart][i]][j]*BetaS[highmodelidx[nm]][j];
     }
     j=0;
     for (l=0;l<r;l++){
         for (k=0;k<KK[l];k++){
     xb+=Xsitest[l][i][k]*BetaS[highmodelidx[nm]][j+nbrcov];j++;
  //printf("XXB=%lf ",Xsitest[l][i][k]);       
}
     }

//double phixb=gsl_cdf_ugaussian_P(xb);
if (typeoutcome1==1){
 double logphixb=-log(2)+gsl_sf_log_erfc(-xb/sqrt(2));
if (Weight1[nm]>0){
LogProb[nm]=logphixb+log(Weight1[nm]);
}
}
//prob[Set[nbrpart][i]]+=phixb*Weight1[nm];
if (typeoutcome1==0){
prob[Set[nbrpart][i]]+=xb*Weight1[nm];
}     
}
if (typeoutcome1==1){
  double logmax=max(nbrmax,LogProb);
 for (nm=0;nm<nbrmax;nm++){
probb+=exp(LogProb[nm]-logmax);
}
prob[Set[nbrpart][i]]=probb*exp(logmax);
}
free(LogProb);
// printf("Test=%d== %lf ",Set[nbrpart][i],probb);
// printf("sumw=%.3lf ",sumw);
// printf("%d== %lf ",Set[nbrpart][i],prob[Set[nbrpart][i]]);
 }
 //free(prob);
 //free_dmatrix(XsiPtest,0,ntest-1,0,P-1);

free_dmatrix(XsiPte,0,n-ntr-1,0,P-1);
free_dmatrix(XsiPt,0,ntr-1,0,P-1);
  for (l=0;l<r;l++){
       free_dmatrix(Xsitest[l],0,ntest-1,0,KK[l]-1); }
 free(Xsitest);
 
free(Weight);free(Weight1);
 free(highmodelidx);
 free_dmatrix(BetaS,0,countmod-1,0,P-1);
 free(lpostRho);


 
 free(respy1);
 
 } //end of nbrpart

free(modelidx);
 free_bmatrix(rhosampleUniq,0,countmod-1,0,P-nbrcov-1);
 free_dmatrix(XsiP,0,n-1,0,P-1);

 }

void predictNewdata(_Bool typeoutcome1,gsl_rng * rr,int n, int ntest,int TT, int r, int *KK, int nbrsample, int P, int nbrcov, _Bool ** rhosample, double *respy,double ***X2Dtest, double **Xcov, double **Xcovtest,double *** Base,double *** XsiMean, double * sigma2xMean,double ** LambdaMean, double sigma2y,double **CovImage, double h, double h0, double a0, double b0, int maxmodel, double *prob){
int j,l,s,i,k;
    int countmod;  int *modelidx=malloc(nbrsample*sizeof(int));
    _Bool ** rhosampleUniq=UniqueModel(nbrsample, P-nbrcov, rhosample,modelidx,&countmod);
     double ** XsiP=dmatrix(0,n-1,0,P-1);
    Combinefeatures(n, r, nbrcov, KK,P, Xcov,XsiMean,XsiP);
 double ** BetaS=dmatrix(0,countmod-1,0,P-1);
 double * lpostRho=malloc(countmod*sizeof(double));
 for (s=0;s<countmod;s++){
double logdetSigma=0; 
double *InvSigmall= GetInvSigmaBeta(rhosampleUniq[s],r,KK,CovImage,h,&logdetSigma);
EstimateBeta(rr,n, nbrcov,P,respy, XsiP, BetaS[s],rhosampleUniq[s],sigma2y, InvSigmall, h0);
lpostRho[s]=logpostRho(n, nbrcov, P,sigma2y,InvSigmall,logdetSigma,respy, XsiP, rhosampleUniq[s],h0, a0, b0);
free(InvSigmall);
 }
//int maxmodel=1000;
  int * highmodelidx=malloc(countmod*sizeof(int));
  sort(countmod,lpostRho,highmodelidx);
   double maxlogpost=lpostRho[0];

   int nbrmax=MIN(maxmodel,countmod);
   for (l=0;l<nbrmax;l++){
     lpostRho[l]=exp(lpostRho[l]-maxlogpost);
   }
   double sumpost=nbrmax*mean(nbrmax,lpostRho);

   for (l=0;l<nbrmax;l++){
     lpostRho[l]=lpostRho[l]/sumpost;
   }

  double *** Xsitest=malloc(r*sizeof(double **));
      for (l=0;l<r;l++){
       Xsitest[l]=dmatrix(0,ntest-1,0,KK[l]-1);
       for (i=0;i<ntest;i++){
         for (k=0;k<KK[l];k++){
           Xsitest[l][i][k]=0;
         }}}

       //Predict scores
 for (l=0;l<r;l++){
     double *yt=malloc(ntest*sizeof(double));
     for (i=0;i<ntest;i++) yt[i]=0;
double * Beta=malloc(KK[l]*sizeof(double));
for (k=0;k<KK[l];k++) Beta[k]=0;
  SampleXsi(rr,0,ntest, TT,KK[l], Base[l],Beta, X2Dtest[l], yt,1, sigma2xMean[l], Xsitest[l],Xsitest[l],LambdaMean[l]);
 free(yt);free(Beta);
 }
 //Prediction of new samples using probit models
 //printf("\nXsip\n");
 double ** XsiPtest=dmatrix(0,ntest-1,0,P-1);
  Combinefeatures(ntest, r, nbrcov, KK,P, Xcovtest,Xsitest,XsiPtest);
 int nm;
// printf("\nEstimated probabilities\n");
 for (i=0;i<ntest;i++){
 prob[i]=0;
     for (nm=0;nm<nbrmax;nm++){
     double xb=0;
     for (j=0;j<P;j++){
 xb+=XsiPtest[i][j]*BetaS[highmodelidx[nm]][j];
     }
if (typeoutcome1==1){ 
double phixb=gsl_cdf_ugaussian_P(xb);
 prob[i]+=phixb*lpostRho[nm];
 } else if (typeoutcome1==0){//continuous response
prob[i]+=xb*lpostRho[nm];
}
}
 //printf("%lf ",prob[i]);
 }
 //free(prob);
 free_dmatrix(XsiPtest,0,ntest-1,0,P-1);
 free(highmodelidx);
  for (l=0;l<r;l++){
       free_dmatrix(Xsitest[l],0,ntest-1,0,KK[l]-1); }
 free(Xsitest);
    free(modelidx);
 free_dmatrix(BetaS,0,countmod-1,0,P-1);
 free(lpostRho);
 free_bmatrix(rhosampleUniq,0,countmod-1,0,P-nbrcov-1);
 free_dmatrix(XsiP,0,n-1,0,P-1);
}


_Bool ** UniqueModel(int nbrsample, int p, _Bool ** rhosample,int * modelidx,int *countd1){
  //%%%% We remove the same models%%%%%%%%%%%%%%%%%%%
  //%%%% idxOfUnique  are the unique indexes%%%%%
  int i;int j=0;
  int countd=0;
  modelidx[countd]=0;
  for (i=0;i<nbrsample;i++){
    for (j=0;j<countd;j++){
    _Bool resultX =compar(p,rhosample[i], rhosample[modelidx[j]]);
      if (resultX==1)
            break;
    }
    if (j==countd){
      modelidx[countd]=i;
      countd++;
    }
    }
  *countd1=countd;


  _Bool ** UniqModel=malloc(countd*sizeof(_Bool*));
  for (i=0;i<countd;i++){
    UniqModel[i]=malloc(p*sizeof(_Bool));
    for (j=0;j<p;j++){
      UniqModel[i][j]=rhosample[modelidx[i]][j];
    }
  }
  printf("\n\n\nNumber of models is= %d",countd);
  //printf("\n");
  return UniqModel;
  }

void Combinefeatures(int n, int r, int nbrcov, int * KK, int P, double **Xcov, 
                      double ***Xsi,double ** XsiP){
int i,k,l;
  int b=0;
  for (k=0;k<nbrcov;k++){
    for (i=0;i<n;i++){
      XsiP[i][b]=Xcov[i][k];
    }
    b++;
  }
  for (l=0;l<r;l++){
    for (k=0;k<KK[l];k++){
      for (i=0;i<n;i++){
        XsiP[i][b]=Xsi[l][i][k];
      }
      b++;
    }}
}


void EstimateBeta(gsl_rng * rr,int n, int nbrcov, int P, 
                  double *respy, double **XsiP, double *BetaP,_Bool * rhoP,
                  double sigma2y, double*InvSigmall,double h0){
  int i,k,k1;
  for (k=0;k<P;k++) BetaP[k]=0;
  int kl; int * IDX=malloc((P-nbrcov)*sizeof(int));
  findc(P-nbrcov,rhoP,0,IDX, &kl);
kl+=nbrcov;
    double *SigmaPriorBeta1=malloc(kl*kl*sizeof(double));
    for (k=0;k<kl;k++){
        for (k1=0;k1<=k;k1++){
            if ((k<nbrcov) && (k1<nbrcov)){
                SigmaPriorBeta1[k*kl+k1]=SigmaPriorBeta1[k1*kl+k]=(1/(h0*sigma2y))*(k==k1);
            } 
        else if ((k>=nbrcov) && (k1>=nbrcov)){ 
            int kk=k-nbrcov;int kk1=k1-nbrcov;
            SigmaPriorBeta1[k*kl+k1]=SigmaPriorBeta1[k1*kl+k]=InvSigmall[kk*(kl-nbrcov)+kk1]/sigma2y;
        } else {
            SigmaPriorBeta1[k*kl+k1]=SigmaPriorBeta1[k1*kl+k]=0;
        }
    }}
    double *beta1=malloc(kl*sizeof(double));
    double ** XsiRho=dmatrix(0,n-1,0,kl-1);
    double *sigma21=malloc(n*sizeof(double));
    //printf("KL=%d ",kl);
    for (i=0;i<n;i++){
      sigma21[i]=sigma2y;
      for (k=0;k<(kl-nbrcov);k++){
        XsiRho[i][k+nbrcov]=XsiP[i][IDX[k]+nbrcov];
      }
      for (k=0;k<nbrcov;k++){
        XsiRho[i][k]=XsiP[i][k];
      }
      }
  
    double *MuPriorBeta=malloc(kl*sizeof(double));
    for (k=0;k<kl;k++) MuPriorBeta[k]=0;
  SampleRegressionCoeff(rr,n, kl, respy, XsiRho, beta1,MuPriorBeta, beta1,
                        sigma21,SigmaPriorBeta1);
  free(sigma21);free(MuPriorBeta); free_dmatrix(XsiRho,0,n-1,0,kl-1);

  for (k=0;k<(kl-nbrcov);k++){
    BetaP[IDX[k]+nbrcov]= beta1[k+nbrcov];
  }
  for (k=0;k<nbrcov;k++){
    BetaP[k]= beta1[k];
  }
  
  free(beta1); free(SigmaPriorBeta1); free(IDX);
}
double LogPostXsi(int n, int r, int nbrcov, int *KK, double **Xcov, 
                  double *respy, double ***Xsi, double **Beta,double *BetaCova,_Bool ** rho,
                  double sigma2y, double *w, double h,double h0, double a0, double b0){
  int i,k,l;
  double logpost=0;
  for (i=0;i<n;i++){
    double xb=0; 
    for (l=0;l<r;l++){
      for (k=0;k<KK[l];k++){
        xb+=Xsi[l][i][k]*Beta[l][k];
      }   
    }
    for (k=0;k<nbrcov;k++){
      xb+=Xcov[i][k]*BetaCova[k];
    }
    
    logpost+= -0.5*pow(respy[i]-xb,2)/sigma2y;
    
  }
  logpost+=-n*log(sigma2y)/2;
  for (l=0;l<r;l++){
    for (k=0;k<KK[l];k++){
      logpost+=-0.5*pow(Beta[l][k],2)/(h*sigma2y)-0.5*log(h*sigma2y);
      logpost+= rho[l][k]*log(w[l])+(1-rho[l][k])*log(1-w[l]);
    }
  }
  ///Prior for BetaCova
  for (k=0;k<nbrcov;k++){
    logpost+= -0.5*pow(BetaCova[k],2)/(h0)-0.5*log(h0*sigma2y);
  }
  //Prior for sigma2y
  logpost+=(a0-1)*log(sigma2y)-b0/sigma2y;
  return logpost;
}

/////
//
 double logpostRho(int n, int nbrcov, int P,double sigma2y,double*InvSigmall,double logdetSigmall,
                   double *respy, double **XsiP, _Bool * rho,double h0, double a0, double b0){  
int kl;int k,k1,i; int * IDX=malloc((P-nbrcov)*sizeof(int));
   findc(P-nbrcov,rho,0,IDX, &kl);
 kl+=nbrcov;
  double ** XsiRho=dmatrix(0,n-1,0,kl-1);
     double *sigma21=malloc(n*sizeof(double));
     //printf("KL=%d ",kl);
     for (i=0;i<n;i++){
       sigma21[i]=sigma2y;
       for (k=0;k<(kl-nbrcov);k++){
         XsiRho[i][k+nbrcov]=XsiP[i][IDX[k]+nbrcov];
       }
       for (k=0;k<nbrcov;k++){
         XsiRho[i][k]=XsiP[i][k];
       }
       }    
 double * InvSigma=malloc(kl*kl*sizeof(double));
 double *SigmaPriorBeta1=malloc(kl*kl*sizeof(double));
     for (k=0;k<kl;k++){
         for (k1=0;k1<=k;k1++){
             if ((k<nbrcov) && (k1<nbrcov)){
                 SigmaPriorBeta1[k*kl+k1]=SigmaPriorBeta1[k1*kl+k]=(1/(h0*sigma2y))*(k==k1);
             }
         else if ((k>=nbrcov) && (k1>=nbrcov)){
             int kk=k-nbrcov;int kk1=k1-nbrcov;
             SigmaPriorBeta1[k*kl+k1]=SigmaPriorBeta1[k1*kl+k]=InvSigmall[kk*(kl-nbrcov)+kk1]/sigma2y;
         } else {
             SigmaPriorBeta1[k*kl+k1]=SigmaPriorBeta1[k1*kl+k]=0;
         }
         double a=0;
         for (i=0;i<n;i++){
a+=XsiRho[i][k]*XsiRho[i][k1]/sigma2y;
         }

InvSigma[k*kl+k1]= InvSigma[k1*kl+k]=a+SigmaPriorBeta1[k*kl+k1];

     }}

  gsl_matrix_view m  = gsl_matrix_view_array (InvSigma,kl,kl);
   gsl_linalg_cholesky_decomp (&m.matrix);
   gsl_vector *mu =gsl_vector_alloc (kl);
   int s;
   double xy[kl];
   for (s=0; s<kl;s++){
     xy[s]=0;
     for (i=0; i<n;i++){
       xy[s]+=XsiRho[i][s]*respy[i]/sigma2y;
       }
   }
   double logdet=gsl_linalg_CD_lndet(&m.matrix);
   gsl_vector_view b= gsl_vector_view_array (xy, kl);
   gsl_linalg_cholesky_solve (&m.matrix, &b.vector, mu);
   double scal=0;
   for (s=0; s<kl;s++){
scal+=xy[s]*gsl_vector_get (mu, s);
 }
   gsl_vector_free (mu);free(InvSigma);free(SigmaPriorBeta1);
 logdet+=nbrcov*log(h0)+logdetSigmall;
  free_dmatrix(XsiRho,0,n-1,0,kl-1);
  free(sigma21); free(IDX);
 return -0.5*logdet+0.5*scal;
 }



double logpostRho1(int n, int nbrcov, int P,  
                  double *respy, double **XsiP, _Bool * rho,
                   double h,double h0, double a0, double b0){
 
  double logdetnewi;double scalnewi=0;
  int i,j,s;
    double a;
    double Sigma[n*n];
    for (i=0;i<n;i++){
      for (j=0;j<=i;j++){
        a=0;
        for (s=0;s<P;s++){
          if (s<nbrcov)
          a+=h0*XsiP[i][s]*XsiP[j][s];
          else a+=h*rho[s-nbrcov]*XsiP[i][s]*XsiP[j][s];
          }
        Sigma[i*n+j]=Sigma[j*n+i]=a;
      }
      Sigma[i*n+i]+=1;
    }
    
    scalnewi=scalar(n,Sigma,respy, &logdetnewi);
    return -0.5*logdetnewi-((n/2.0)+a0)*log(1+scalnewi/(2*b0));
  
}

double scalar(int n,double A[n*n],double b_data[n], double *det)
{
  gsl_matrix_view m
  = gsl_matrix_view_array (A, n,n);
  
  gsl_vector_view b
    = gsl_vector_view_array (b_data, n);
  
  gsl_vector *x = gsl_vector_alloc (n);
  
  //  int s;
  
  //gsl_permutation * p = gsl_permutation_alloc (4);
  
  gsl_linalg_cholesky_decomp (&m.matrix);
  
  gsl_linalg_cholesky_solve (&m.matrix, &b.vector, x);
  *det=gsl_linalg_CD_lndet(&m.matrix);
  double scal=0;
  int q;    
  for (q=0;q<n;q++){
    scal+=b_data[q]*gsl_vector_get (x, q);
  }
  // gsl_permutation_free (p);
  gsl_vector_free (x);
  return  scal;
}

double gsl_linalg_CD_lndet (gsl_matrix * CD)
{
  size_t i, n = CD->size1;
  
  double lndet = 0.0;
  
  for (i = 0; i < n; i++)
  {
    lndet += log (gsl_matrix_get (CD, i, i));
  }
  //gsl_matrix_free(CD);
  return 2*lndet;
}



_Bool compar(int n,_Bool *u,_Bool *v){
  int i;
  for (i=0;i<n;i++){
    if (u[i]!=v[i])
      return 0;
  }
  return 1;
}

