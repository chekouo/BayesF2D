#include <gsl/gsl_rng.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <time.h>
#include <gsl/gsl_randist.h>
#include<stdio.h>
#include <math.h>
#include <string.h>
#include <gsl/gsl_cdf.h>
# include "myfunction.h"


//void MCMC(char **typeOutcome,int * n1,int* ntest1,int *TT1,int *r1,int *nbrcov1,double * resp, double *X,double *Xtest, int *KK, double* Base1, double* Xcov1,double* Xcovtest1,int* Kf, int *nbrset,int * ListTest, int * ListTraining, int *nbrsample1,int *burninsample1,int* Predict1, 
       //   double * Betasample,double *XsiMean1,double* rhoMean1,double *logpost,double *prob,double * CovImageSample){
  
  void MCMC(int* typeOutcome,int * n1,int *TT1,int *r1,int *nbrcov1,double * resp, double *X,int *KK, double* Base1, double* Xcov1, int *nbrsample1,int *burninsample1, double * Betasample,double *XsiMean1,double* rhoMean1,double * CovImageSample,double *LambdaMean1,double* seed1,double *hypsigm_rr,double *h1,int *covbool,int *rhosample1, double *yMean,double *sigma2xMean,double * sigma2yMean1){
  int i,l,k,t;
  int bs=0;
  clock_t t1 = clock();
  //char * typeOutcome1=typeOutcome[0];
  int typeOutcome1=typeOutcome[0];
  if (typeOutcome1==1){
  //printf("The Type of outcome is %s\n",typeOutcome1);
  printf("The Type of outcome is %s\n","binary"); 
  } else if (typeOutcome1==2){
    printf("The Type of outcome is %s\n","continuous"); 
  }
  int nbrsample=nbrsample1[0];
  printf("Number of MCMC samples after burn-in is %d\n",nbrsample);
  int burninsample=burninsample1[0];
  printf("Number of burn-in is %d\n",burninsample);
  int n=n1[0];
  printf("Number of samples is %d\n",n);
  int TT=TT1[0];
  printf("Number of intensities squared is %d\n",TT);
  int r=r1[0];
  printf("Number of imaging sequences is %d\n",r);
  int nbrcov=nbrcov1[0];
  printf("Number of covariates is %d\n",nbrcov);
  int j;
  
double *** X2D=malloc(r*sizeof(double **));
double ** Xcov=malloc(n*sizeof(double*));
int j1=0;
for (i=0;i<n;i++){
  Xcov[i]=malloc(nbrcov*sizeof(double));
  for (j=0;j<nbrcov;j++){
    Xcov[i][j]=Xcov1[j1];
    j1++;
  }
}
 k=0;
 for (l=0;l<r;l++){
   X2D[l]=dmatrix(0,n-1,0,TT-1);
   for (i=0;i<n;i++){
   for (t=0;t<TT;t++){
       X2D[l][i][t]=X[k];k++;
     }}
 }

//get test data//

 double *** Base=malloc(r*sizeof(double**));
 int k1=0;
 for (l=0;l<r;l++){
Base[l]=dmatrix(0,TT-1,0,KK[l]-1);
  for (t=0;t<TT;t++){
    for (k=0;k<KK[l];k++){
  Base[l][t][k]=Base1[k1];
   //printf("%lf ", Base1[k1]);
   k1++;
 }
 }
 }

 /* We set hyperparameters */
 //double h=1; //Shinkage hyper of regression effects
 double h=h1[0]; //Shinkage hyper of regression effects
 double h0=1000;
 double a0=0.001; double b0=0.001;
// double alpha0= 0.5 ;double beta0=0.5; 
/*Hyperparameter of sigma_rr'*/ 
double alpha0=hypsigm_rr[0] ;double beta0=hypsigm_rr[1]; 
long seed=seed1[0];
 gsl_rng * rr = gsl_rng_alloc (gsl_rng_rand48);
 gsl_rng_set (rr, seed);
 double* w=malloc(r*sizeof(double));
 /* Parameters*/
 int l1;  
 double ** CovImage=malloc(r*sizeof(double*));
 double ** CovImageMean=malloc(r*sizeof(double*)); 
 for (l=0;l<r;l++) {
      CovImage[l]=malloc(r*sizeof(double));
      CovImage[l][l]=0;
       CovImageMean[l]=malloc(r*sizeof(double));
       CovImageMean[l][l]=0;
  for (l1=0;l1<l;l1++) {
      CovImage[l][l1]=CovImage[l1][l]=0.5;
  CovImageMean[l][l1]=CovImageMean[l1][l]=0;
 CovImage[l][l1]=CovImage[l1][l]=0;
  }
  }
for (l=0;l<r;l++) {
 for (l1=0;l1<r;l1++){
       if (l1!=l) CovImage[l][l]+=CovImage[l][l1];
   }
   CovImage[l][l]+=h;
   //printf("h=\n");
   }
 
  //printf("%lf ",CovImage[l][l1]);
 double sigma2y=0.1;//variance of the response
 double sigma2yMean=0;//Mean of the response variable
  double *sigma2x=malloc(r*sizeof(double));//variance of Xsi
 // double *sigma2xMean=malloc(r*sizeof(double));//variance of Xsi
  double *sigma21=malloc(n*sizeof(double));
  for (i=0;i<n;i++){sigma21[i]=sigma2y;}
double ** Lambda=malloc(r*sizeof(double*)); //variance for latent variables \xsi
 double ** LambdaMean=malloc(r*sizeof(double*)); //variance for latent variables \xsi 
double ** Beta=malloc(r*sizeof(double*));
 double * BetaCova=malloc(nbrcov*sizeof(double));
 for (j=0;j<nbrcov;j++){
  BetaCova[j]=0;
 }
 _Bool ** rho=malloc(r*sizeof(_Bool*));
 int P=nbrcov;
 for (l=0;l<r;l++)
 for (k=0;k<KK[l];k++) P++;
 _Bool ** rhosample=bmatrix(0,nbrsample-1,0,P-nbrcov-1);
 double ** rhoMean=malloc(r*sizeof(double*));
 for (l=0;l<r;l++){
   sigma2x[l]=0.1;sigma2xMean[l]=0;
   rho[l]=malloc(KK[l]*sizeof(_Bool));
   rhoMean[l]=malloc(KK[l]*sizeof(double));
   Beta[l]=malloc(KK[l]*sizeof(double));
   Lambda[l]=malloc(KK[l]*sizeof(double));
    LambdaMean[l]=malloc(KK[l]*sizeof(double));
   w[l]=0.5;
   for (k=0;k<KK[l];k++){
        LambdaMean[l][k]=0;
     Lambda[l][k]=0.1; rhoMean[l][k]=0;
   double  uni=gsl_ran_flat (rr, 0, 1);
   if (uni<=0.5) { rho[l][k]=1; Beta[l][k]=1;  
   } else {rho[l][k]=0;Beta[l][k]=0;
   }
   }
 }
 //Initialization of Xsi's
 double *** Xsi=malloc(r*sizeof(double **));
  double *** XsiDelete=malloc(r*sizeof(double **));
  double *** XsiMean=malloc(r*sizeof(double **));
 for (l=0;l<r;l++){
  Xsi[l]=dmatrix(0,n-1,0,KK[l]-1);
  XsiMean[l]=dmatrix(0,n-1,0,KK[l]-1);
  XsiDelete[l]=dmatrix(0,n-1,0,KK[l]-1);
   for (i=0;i<n;i++){
     for (k=0;k<KK[l];k++){
       Xsi[l][i][k]=gsl_ran_ugaussian(rr);
      XsiMean[l][i][k]=0; XsiDelete[l][i][k]=0;
     }}
 }



double *respy=malloc(n*sizeof(double));
// double *yMean=malloc(n*sizeof(double));
 for (i=0;i<n;i++) yMean[i]=0;
//if (strcmp(typeOutcome1,"continuous")==0){
   if (typeOutcome1==2){
   for (i=0;i<n;i++){
     respy[i]=resp[i];
   }
// } else if (strcmp(typeOutcome1,"binary")==0){
   } else if (typeOutcome1==1){
   for (i=0;i<n;i++){
   if (resp[i]==1) respy[i]=Truncate(0,1,0,rr);
     else respy[i]=-Truncate(0,1,0,rr);
 }
 }
 

int s;int cvs=0;
for (s=0;s<burninsample+nbrsample;s++){
   int rs=0; 
    double *y=malloc(n*sizeof(double));
  //Sample sigma2x
  //if (strcmp(typeOutcome1,"continuous")==0){
    if (typeOutcome1==2){
  SampleSig2(rr,n,KK, r,a0,b0,&sigma2y,Xsi, respy,Beta,h,rho);
    } else if (typeOutcome1==1){
    sigma2y=1;
    Sampley(rr,n,KK,r,Xsi, respy,resp,Beta,Xcov,BetaCova,nbrcov);
  }
  if (s>=burninsample) sigma2yMean+=sigma2y/nbrsample;
//Sample BetaCovariate
  for (i=0;i<n;i++){
    double xb=0;
    for (l=0;l<r;l++){
        for (k=0;k<KK[l];k++){
      xb+=Xsi[l][i][k]*Beta[l][k];
        }}
    y[i]=respy[i]-xb;
  }
  SampleBetaCova(rr,n,nbrcov,Xcov, h0, sigma2y,y, BetaCova);
  
  
  /// sample Beta's
  int l1;
  for (l=0;l<r;l++){
      double scal=0;
     double Muscal=0;
    for (i=0;i<n;i++){
      double xb=0; 
      for (l1=0;l1<r;l1++){
        if (l1!=l){
          for (k=0;k<KK[l1];k++){
            xb+=Beta[l1][k]*Xsi[l1][i][k];
          }}}
      for (j=0;j<nbrcov;j++){
        xb+=Xcov[i][j]*BetaCova[j];
      }
      y[i]=respy[i]-xb;
   
    }
  SampleRho(l,r,rr,n,KK[l],y,Xsi[l],sigma2y,CovImage,   &scal,&Muscal,rho,KK,w[l],Beta,h);   
    if (s>=burninsample){
      for (k=0;k<KK[l];k++){
          rhosample[s-burninsample][rs]=rho[l][k];
          rs++;
        rhoMean[l][k]+= (double) rho[l][k]/nbrsample;
      } 
    }
  SampleBeta(rr,n,KK[l],rho[l], Xsi[l], CovImage[l][l], sigma2y,y, Beta[l],scal,Muscal);
   
   // We sample Xsi
  SampleXsi(rr,1,n, TT, KK[l], Base[l],Beta[l], X2D[l], y, 
             sigma2y, sigma2x[l], Xsi[l],XsiDelete[l],Lambda[l]);

      for (k=0;k<KK[l];k++){
        SampleLambda(rr,n,k,a0, b0,&Lambda[l][k], Xsi[l]) ;
      }
     SampleSig2X(rr,n,KK[l], TT, a0,  b0,&sigma2x[l], Base[l], X2D[l], Xsi[l]);
   
  if (s>=burninsample){
  sigma2xMean[l]+=sigma2x[l]/nbrsample;
      for (k=0;k<KK[l];k++){
   Betasample[bs]=Beta[l][k];
   bs++;
   LambdaMean[l][k]+=Lambda[l][k]/nbrsample;
   for (i=0;i<n;i++){
XsiMean[l][i][k]+=Xsi[l][i][k]/nbrsample;
//Xsisample[xs]=Xsi[l][i][k];
//xs++;
   }
   }
   }
  }// End for l
  free(y);
//logpost[s]=LogPost(n, r, TT, nbrcov, KK, Base, X2D, Xcov,respy,Xsi,Beta,BetaCova,rho,sigma2y, sigma2x,w, h, h0,Lambda, a0,b0);
if (covbool[0]==1){
SampleCovImage(rr,r,KK, CovImage,Beta,h, rho,alpha0,beta0);
}
if (s>=burninsample){
      for (i=0;i<n;i++){
      yMean[i]+=respy[i]/nbrsample;}
    for (l=0;l<r;l++){
        for (l1=0;l1<=l;l1++){
        CovImageMean[l][l1]+=CovImage[l][l1]/nbrsample;
        }
        for (l1=0;l1<r;l1++){
           //CovImageMean[l][l1]=CovImageMean[l1][l];
           CovImageSample[cvs]=CovImage[l][l1];cvs++;
        }
}
}

//Print some values
if (s%((nbrsample+burninsample+5)/5)==1){
printf("The number of MCMC iterations is %d\n",s);
   for (l=0;l<r;l++){
    printf("variance of x for %d is %lf\n",l,sigma2x[l]);
   }
   for (l=0;l<r;l++){
printf("The %dth lambda\n",l);
       for (k=0;k<KK[l];k++)
       printf("%lf ",Lambda[l][k]);
   printf("\n");
   }
    printf("variance of response y is %lf\n",sigma2y);
   for (l=0;l<r;l++){
    printf("Rho for sequence %d\n",l);
     for (k=0;k<KK[l];k++)
      printf("%d ",rho[l][k]);
    printf("\n");
    }
    printf("\n");
for (l=0;l<r;l++){
 printf("\n Beta values for sequence %d\n",l);
  for (k=0;k<KK[l];k++)
printf("%lf ",Beta[l][k]);
printf("\n");
}
    printf("\n");
  printf("\n Beta for covariates\n");
  for (j=0;j<nbrcov;j++){
    printf("%lf ",BetaCova[j]);
  }
printf("\n");
//R_CheckUserInterrupt();
}

}
//get some outcomes

bs=0;int bs1=0;
for (s=0;s<nbrsample;s++){
bs1=0;
for (l=0;l<r;l++){
  for (k=0;k<KK[l];k++){
   rhosample1[bs]=rhosample[s][bs1];
   bs++;bs1++;
}
}
}

bs=0;int rx=0;
for (l=0;l<r;l++){
for (k=0;k<KK[l];k++){
 rhoMean1[bs]=rhoMean[l][k];
LambdaMean1[bs]=LambdaMean[l][k];
 bs++;
 for (i=0;i<n;i++){
     XsiMean1[rx]=XsiMean[l][i][k];
     rx++;
 }
    printf("%lf ",rhoMean[l][k]);
}
printf("\n");
}
 printf("Cov Estimate\n");
for (l=0;l<r;l++){
  for (l1=0;l1<=l;l1++){
    CovImageMean[l1][l]=CovImageMean[l][l1];
  }
}

for (l=0;l<r;l++){
    for (l1=0;l1<r;l1++){
        printf("%lf, ",CovImageMean[l][l1]);
}
 printf("\n");
}

sigma2yMean1[0]=sigma2yMean;




/*** Prediction performance*/

//List unique models


for (l=0;l<r;l++){
//    free_dmatrix(X2Dtest[l],0,ntest-1,0,TT-1);
     free(rho[l]);free(rhoMean[l]);free(Beta[l]);free(Lambda[l]);free(LambdaMean[l]);
   free_dmatrix(X2D[l],0,n-1,0,TT-1);
   free_dmatrix(Xsi[l],0,n-1,0,KK[l]-1); 
   free_dmatrix(XsiMean[l],0,n-1,0,KK[l]-1);
   free_dmatrix(XsiDelete[l],0,n-1,0,KK[l]-1);
  free_dmatrix(Base[l],0,TT-1,0,KK[l]-1); 
   }
   free(X2D);
   free_bmatrix(rhosample,0,nbrsample-1,0,P-nbrcov-1);free(Beta);free(rho);
   free(Lambda);free(LambdaMean);free_dmatrix(Xcov,0,n-1,0,nbrcov-1);
   gsl_rng_free (rr);free(rhoMean);free(sigma2x);
   free(Base);
   free(respy);//free(sigma2xMean);//free(yMean);
   free(Xsi);
   free(BetaCova);
   free(XsiMean);free(XsiDelete);
   free(w);free(sigma21);free_dmatrix(CovImage,0,r-1,0,r-1);free_dmatrix(CovImageMean,0,r-1,0,r-1);
   t1 = clock() - t1;
   double  time_taken = ((double)t1)/CLOCKS_PER_SEC; // in seconds
   printf("\n\nTime taken in seconds is %f\n",time_taken);
   printf("\nTime taken in minutes is %f\n",time_taken/60);
   // printf("\nTime taken in hours is %f\n",time_taken/3600);
   }
