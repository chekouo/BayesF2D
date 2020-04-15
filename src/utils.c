#include <gsl/gsl_rng.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include <time.h>
#include <gsl/gsl_randist.h>
#include <math.h>
# include "myfunction.h"
static const double t4 = 0.45;

void SampleCovImage(gsl_rng * rr,int r,int *KK1, double **CovImage,double ** Beta,double h, _Bool ** rho,double alpha0, double beta0){
int k,k1,l1,l2,l3;
int *nx=malloc((r+1)*sizeof(int));
 nx[0]=0;
   int *nxCum=malloc((r+1)*sizeof(int));
 //for (l1=0;l1<r+1;l1++) nxCum[l1]=0;
   int T=0;
for (l1=0;l1<r;l1++){
    T+=KK1[l1];
}

_Bool rhop[T];
int t=0;
for (l1=0;l1<r;l1++){
    for (k=0;k<KK1[l1];k++){
        rhop[t]=rho[l1][k];t++;
    }
}

//for (t=0;t<T;t++) printf("Rhop=%d ",rhop[t]);

  int Kl=0;
   int **IDX=malloc(r*sizeof(int*));
   for (l1=0;l1<r;l1++){
       int KK=KK1[l1];
   IDX[l1]=malloc(KK*sizeof(int));
   findc(KK,rho[l1],0,IDX[l1], &nx[l1+1]);
   Kl+=nx[l1+1];
   }
for (l1=0;l1<r+1;l1++){
       nxCum[l1]=0;
   for (l2=0;l2<=l1;l2++)
        nxCum[l1]+=nx[l2];}

 for (l1=0;l1<r;l1++){
  CovImage[l1][l1]=h;
       for (l2=0;l2<r;l2++){
   if (l1!=l2) CovImage[l1][l1]+=CovImage[l1][l2]*nx[l2+1];
    }
    }
/*
printf("\n");
 for (l1=0;l1<r;l1++){
        for (l2=0;l2<r;l2++){
    printf("%lf ",CovImage[l1][l2]);
     }
       printf("\n");
     }
*/
     double BetaList[Kl];
 for (k1=0;k1<Kl;k1++){
    for (l3=0;l3<r;l3++){
   if ((nxCum[l3]<=k1)&&(k1<nxCum[l3+1]))
     BetaList[k1]=Beta[l3][IDX[l3][k1-nxCum[l3]]];
     }}
double alphap=1; double betap=1;
for (l1=0;l1<r;l1++){
    for (l2=0;l2<l1;l2++){
        double oldcov=CovImage[l1][l2];double oldll=CovImage[l1][l1]; double oldll2=CovImage[l2][l2];
        double logdet=0;
     double *InvSigmaBeta = GetInvSigmaBeta(rhop,r,KK1, CovImage,h, &logdet); 
         double logpostold=0;
         double scal=0;
         for (k=0;k<Kl;k++){
             double xx=0;
             for (k1=0;k1<Kl;k1++){
                 xx+=InvSigmaBeta[k*Kl+k1]*BetaList[k1];}
    scal+=BetaList[k]*xx;
         }
         logpostold=-0.5*logdet-0.5*scal+log(gsl_ran_gamma_pdf(oldcov,alpha0,beta0));
    double sigmaprop=gsl_ran_gamma(rr, alphap, betap);
   CovImage[l1][l2]=CovImage[l2][l1]=sigmaprop;
 //printf("SCAL1=%.2lf ", scal);
    double * InvSigmaBeta1 = GetInvSigmaBeta(rhop,r,KK1, CovImage,h, &logdet);     
   scal=0;
    for (k=0;k<Kl;k++){
              double xx=0;
              for (k1=0;k1<Kl;k1++){
                  xx+=InvSigmaBeta1[k*Kl+k1]*BetaList[k1];}
     scal+=BetaList[k]*xx;
          }
//printf("SCAL2=%.2lf ", scal);

    double logpostnew=-0.5*logdet-0.5*scal+log(gsl_ran_gamma_pdf(CovImage[l1][l2],alpha0,beta0));;
 double rat=logpostnew-logpostold+log(gsl_ran_gamma_pdf(oldcov,alphap,betap))-log(gsl_ran_gamma_pdf(CovImage[l1][l2],alphap,betap));
double u=gsl_ran_flat (rr, 0, 1);
if (log(u)>rat){
CovImage[l1][l2]=CovImage[l2][l1]=oldcov;
CovImage[l1][l1]=oldll;CovImage[l2][l2]=oldll2;
}
free(InvSigmaBeta); free(InvSigmaBeta1);
}
}
for (l1=0;l1<r;l1++) free(IDX[l1]);
free(IDX);free(nx);free(nxCum);
}

double *GetInvSigmaBeta(_Bool *rhop,int r,int *KK, double **CovImage,double h, double *logdet){
int *nx=malloc((r+1)*sizeof(int));
nx[0]=0;
   int *nxCum=malloc((r+1)*sizeof(int));
  int l1,l2,k,k1;
   int Kl=0;
  int **IDX=malloc(r*sizeof(int*));
      int Kk=0;
  for (l1=0;l1<r;l1++){
      int Kp=KK[l1]; _Bool rho[Kp];
      //if (l1!=0){Kk=KK[l1-1];}
    for (k=0;k<Kp;k++){
     rho[k]=rhop[k+Kk];
   }
    Kk=Kk+Kp;
  IDX[l1]=malloc(Kp*sizeof(int));
  findc(Kp,rho,0,IDX[l1], &nx[l1+1]);
  Kl+=nx[l1+1];
  }
  for (l1=0;l1<r;l1++){
 CovImage[l1][l1]=h;
      for (l2=0;l2<r;l2++){
  if (l1!=l2) CovImage[l1][l1]+=CovImage[l1][l2]*nx[l2+1];
  //printf("%lf ",CovImage[l1][l2]);
   }
    //  printf("\n");
   }

  // printf("\n");
/*
   for (l1=0;l1<r;l1++){
       for (l2=0;l2<r;l2++){
   printf("%lf ",CovImage[l1][l2]);
    }
      printf("\n");
    }
  */
  for (l1=0;l1<r+1;l1++){
      nxCum[l1]=0;
  for (l2=0;l2<=l1;l2++)
       nxCum[l1]+=nx[l2];}
 double * Sigmal=malloc(Kl*Kl*sizeof(double));
  for (k=0;k<Kl;k++){
  for (k1=0;k1<=k;k1++){
  for (l1=0;l1<r;l1++){
     for (l2=0;l2<=l1;l2++){
  if ((nxCum[l1]<=k)&&(k<nxCum[l1+1])&&(nxCum[l2]<=k1)&&(k1<nxCum[l2+1])){
  if ((l1==l2)&&(k1!=k)){
  Sigmal[k*Kl+k1]=Sigmal[k1*Kl+k]=0;
  } else {
          Sigmal[k*Kl+k1]=Sigmal[k1*Kl+k]=CovImage[l1][l2];
      }
  }
  }}
  }
  }
/*
   for (k=0;k<Kl;k++){
   for (k1=0;k1<Kl;k1++){
       printf("%lf, ",Sigmal[k*Kl+k1]);
   }
   printf("\n");
   }
*/
   gsl_matrix_view m  = gsl_matrix_view_array (Sigmal, Kl,Kl);
    gsl_linalg_cholesky_decomp (&m.matrix); 
 *logdet=gsl_linalg_CD_lndet(&m.matrix);
    gsl_linalg_cholesky_invert (&m.matrix);

for (k=0;k<Kl;k++){
     for (k1=0;k1<=k;k1++){
Sigmal[k*Kl+k1]=Sigmal[k1*Kl+k]=gsl_matrix_get(&m.matrix,k,k1);
     }
}
free(nx);free(nxCum);
for (l1=0;l1<r;l1++) {free(IDX[l1]);}
 free(IDX);
return Sigmal;
}

void GetH (int l,int r, _Bool **rho,double *scal1,double *Muscal1,double **Beta, double h1,double ** Sigma, int *KK){
 double scal=0;double Muscal=0;
 Sigma[l][l]=h1;
 int l1,l2,k1;
 int Kl=0;
 int k=0;
 int * ID=malloc((r-1)*sizeof(int));
 int *nx=malloc(r*sizeof(int));
  int *nxCum=malloc(r*sizeof(int));
 nx[0]=0;
 int ** IDX=malloc((r-1)*sizeof(int*));
 for (l1=0;l1<r;l1++) nxCum[l1]=0;
 for (l1=0;l1<r;l1++){
 if (l1!=l){
 ID[k]=l1;
 IDX[k]=malloc(KK[l1]*sizeof(int));
 findc(KK[l1],rho[l1],0,IDX[k], &nx[k+1]);
 Kl+=nx[k+1];
 Sigma[l][l]+=Sigma[l][l1]*nx[k+1];
 k++;
     }
 }

 for (l1=0;l1<r;l1++){
     Sigma[l1][l1]=h1;
for (l2=0;l2<r;l2++){
    if (l1!=l2){
    for (k=0;k<KK[l2];k++)
    Sigma[l1][l1]+=rho[l2][k]*Sigma[l1][l2];

 }}}
/*
 printf("\nCovImage\n");
  for (l1=0;l1<r;l1++){
      for (l2=0;l2<r;l2++){
printf("%lf ",Sigma[l1][l2]);
      }
   printf("\n");
  }
*/
  for (l1=0;l1<r;l1++)
 for (l2=0;l2<=l1;l2++)
      nxCum[l1]+=nx[l2];
if (Kl!=0){
 double * Sigmal=malloc(Kl*Kl*sizeof(double));
 for (k=0;k<Kl;k++){
 for (k1=0;k1<=k;k1++){
     for (l1=0;l1<r-1;l1++){
    for (l2=0;l2<=l1;l2++){
 if ((nxCum[l1]<=k)&&(k<nxCum[l1+1])&&(nxCum[l2]<=k1)&&(k1<nxCum[l2+1])){
 if ((l1==l2)&&(k1!=k)){
 Sigmal[k*Kl+k1]=Sigmal[k1*Kl+k]=0;
 } else {
         Sigmal[k*Kl+k1]=Sigmal[k1*Kl+k]=Sigma[ID[l1]][ID[l2]];
 //printf("xx=%lf " ,Sigma[ID[l1]][ID[l2]]);
 }}

 }}
 //printf("%lf ",Sigmal[k*Kl+k1]);
 }
 //printf("\nHH=%d\n",l);
 }
/*
 printf("\nHH=%d\n",l);
 for (k=0;k<Kl;k++){
     for (k1=0;k1<Kl;k1++){
printf("%lf, ",Sigmal[k*Kl+k1]);
     }
 printf("\n");
 }
 */
 
 gsl_matrix_view m  = gsl_matrix_view_array (Sigmal, Kl,Kl);
 gsl_linalg_cholesky_decomp (&m.matrix);
 gsl_vector *mu =gsl_vector_alloc (Kl);
    double xy[Kl];double bet[Kl];
    for (k=0; k<Kl;k++){
        for (l1=0;l1<r-1;l1++){
            if ((nxCum[l1]<=k)&&(k<nxCum[l1+1])){
                xy[k]=Sigma[l][ID[l1]];
 //           printf("\nID=%d kk=%d\n",ID[l1],IDX[l1][k-nxCum[l1]]);
                bet[k]=Beta[ID[l1]][IDX[l1][k-nxCum[l1]]];
            }
        }
  //        printf("XY=%lf ",bet[k]);
    }
  //   printf("\n\n");
    gsl_vector_view b= gsl_vector_view_array (xy, Kl);
    gsl_linalg_cholesky_solve (&m.matrix, &b.vector, mu);
 for (k=0; k<Kl;k++){
     scal+=xy[k]*gsl_vector_get(mu,k);
 }

gsl_vector_view b1= gsl_vector_view_array (bet, Kl);
   gsl_linalg_cholesky_solve (&m.matrix, &b1.vector, mu);
    for (k=0; k<Kl;k++){
       Muscal+=xy[k]*gsl_vector_get(mu,k);
    }

 //printf("Scal=%lf\n  HH= %d" ,scal,l);
free(Sigmal); 
 gsl_vector_free (mu);//free(IDx);
 }
  *scal1=scal;*Muscal1=Muscal;
 // *scalinv=scal/(Sigma[l][l]*(Sigma[l][l]-nxl*scal));
 //printf("Logdet=%lf",*logdethl);
free(ID);free(nx);free(nxCum);
for (l1=0;l1<r-1;l1++) {free(IDX[l1]);} free(IDX);

}


double LogPost(int n, int r, int TT, int nbrcov, int *KK, double *** Base, double ***X2D, double **Xcov, 
               double *respy, double ***Xsi, double **Beta,double *BetaCova,_Bool ** rho,
               double sigma2y, double* sigma2x, double *w, double h,double h0,
               double** Lambda, double a0, double b0){
  int i,k,t,l;
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
    for (t=0;t<TT;t++){
      for (l=0;l<r;l++){
      double xsb=0;
      for (k=0;k<KK[l];k++){
        xsb+=Xsi[l][i][k]*Base[l][t][k];
      }
      logpost+=-0.5*pow(X2D[l][i][t]-xsb,2)/sigma2x[l]-0.5*log(sigma2x[l]);
      }
    }
    }
  logpost+=-n*log(sigma2y)/2;
  for (l=0;l<r;l++){
    for (k=0;k<KK[l];k++){
  logpost+=rho[l][k]*(-0.5*pow(Beta[l][k],2)/(h*sigma2y)-0.5*log(h*sigma2y));
  logpost+= rho[l][k]*log(w[l])+(1-rho[l][k])*log(1-w[l]);
  //Prior for lambda
  logpost+=(a0-1)*log(Lambda[l][k])-b0/Lambda[l][k];
    }
    //Prior for sigma2x
    logpost+=(a0-1)*log(sigma2x[l])-b0/sigma2x[l];
    }
  ///Prior for BetaCova
  for (k=0;k<nbrcov;k++){
    logpost+= -0.5*pow(BetaCova[k],2)/(h0)-0.5*log(h0*sigma2y);
  }
  //Prior for sigma2y
  logpost+=(a0-1)*log(sigma2y)-b0/sigma2y;
  return logpost;
}


void SampleXsi(gsl_rng * rr,_Bool includeOutcome,int n, int TT, int p, double ** Base, double *Beta,double **X2D, double *y,  double sigma2y, double sigma2x, double **Xsi,double **Mu, double* Lambda){
  int i,k,t;
  double **  Xx=dmatrix(0,TT,0,p-1);
  for (k=0;k<p;k++){
    for (t=0;t<TT;t++){
      Xx[t][k]=Base[t][k];
    }
    Xx[TT][k]=Beta[k];
  }
  double * sigma2vx=malloc((TT+1)*sizeof(double));
  for (t=0;t<TT;t++){
    sigma2vx[t]=sigma2x;
  }
  sigma2vx[TT]=sigma2y;
  for (i=0;i<n;i++){
    double* yx=malloc((TT+1)*sizeof(double));
    for (t=0;t<TT;t++){
      yx[t]=X2D[i][t];
    }
    yx[TT]=y[i];
int k1;
    //double *Mu=malloc(p*sizeof(double));
double *SigmaPriorBeta1=malloc(p*p*sizeof(double));
for (k=0;k<p;k++)
         for (k1=0;k1<=k;k1++)
         SigmaPriorBeta1[k*p+k1]=SigmaPriorBeta1[k1*p+k]=(1/Lambda[k])*(k==k1);    
double *MuPrior=malloc(p*sizeof(double));
for (k=0;k<p;k++){
       MuPrior[k]=0;
   }

if (includeOutcome==1){
    SampleRegressionCoeff(rr,TT+1,p, yx, Xx, Xsi[i],MuPrior,Mu[i], sigma2vx,SigmaPriorBeta1);
    } else {
      SampleRegressionCoeff(rr,TT,p, yx, Xx, Xsi[i],MuPrior, Mu[i], sigma2vx,SigmaPriorBeta1);
    }
    //free(Mu);
    free(MuPrior);free(yx); free(SigmaPriorBeta1);
  }
  free_dmatrix(Xx,0,TT,0,p-1);free(sigma2vx);
}
void SampleBetaCova(gsl_rng * rr,int n,int p, double ** Xcov, double h, double sigma2y,double *y, double *Beta){
  int i,k,k1;
  int kl=p;
    double *SigmaPriorBeta1=malloc(kl*kl*sizeof(double));
    for (k=0;k<kl;k++)
        for (k1=0;k1<=k;k1++)
        SigmaPriorBeta1[k*kl+k1]=SigmaPriorBeta1[k1*kl+k]=(1/(h*sigma2y))*(k==k1);
    double *sigma21=malloc(n*sizeof(double));
    for (i=0;i<n;i++){
      sigma21[i]=sigma2y;
    }
    double *Mu=malloc(kl*sizeof(double));
    double *MuPrior=malloc(p*sizeof(double));
 for (k=0;k<kl;k++){
        MuPrior[k]=0;
    }
    SampleRegressionCoeff(rr,n, kl, y, Xcov, Beta,MuPrior,Mu, sigma21,SigmaPriorBeta1);
    free(Mu);free(MuPrior);
   free(SigmaPriorBeta1);free(sigma21);
}

void SampleBeta(gsl_rng * rr,int n,int p,_Bool *rho, double ** Xsi, double Sigmall, double sigma2y,double *y, double *Beta,double scal,double Muscal){
  int i,k,k1;
for (k=0;k<p;k++) Beta[k]=0;
int kl; int * IDX=malloc(p*sizeof(int));
findc(p,rho,0,IDX, &kl);

double scalinv=scal/(Sigmall*(Sigmall-kl*scal));
if (kl>0){
  double *SigmaPriorBeta1=malloc(kl*kl*sizeof(double));
  //for (k=0;k<kl;k++){SigmaPriorBeta1[k]=Sigmall*sigma2y;}
 for (k=0;k<kl;k++)
          for (k1=0;k1<=k;k1++)
          SigmaPriorBeta1[k*kl+k1]=SigmaPriorBeta1[k1*kl+k]=(1/sigma2y)*(scalinv+(1/Sigmall)*(k==k1));
  
  double *beta1=malloc(kl*sizeof(double));
  double ** XsiRho=dmatrix(0,n-1,0,kl-1);
  double *sigma21=malloc(n*sizeof(double));
  for (i=0;i<n;i++){
    sigma21[i]=sigma2y;
    for (k=0;k<kl;k++){
      XsiRho[i][k]=Xsi[i][IDX[k]];
    }}
  double *Mu=malloc(kl*sizeof(double));
  double *MuPrior=malloc(kl*sizeof(double));
  for (k=0;k<kl;k++){
      MuPrior[k]=Muscal;
  }
  SampleRegressionCoeff(rr,n, kl, y, XsiRho, beta1,MuPrior,Mu,sigma21,SigmaPriorBeta1);
  for (k=0;k<kl;k++){
  Beta[IDX[k]]= beta1[k];
// printf("Bet=%lf ",beta1[k]);
  }
  free(beta1);free(SigmaPriorBeta1);free(sigma21);
  free_dmatrix(XsiRho,0,n-1,0,kl-1);free(Mu);free(MuPrior);
}
free(IDX);
}

void Sampley(gsl_rng * rr,int n,int * KK, int r, double ***Xsi, double *y, double *resp, double **Beta,double ** Xcov,double *BetaCov, int nbrcov){
  int l,i,j,k;

  for (i=0;i<n;i++){
    double xb=0;
  for (l=0;l<r;l++){
    for (k=0;k<KK[l];k++){
      xb+=Xsi[l][i][k]*Beta[l][k];
    }   
    }
  for (j=0;j<nbrcov;j++){
         xb+=Xcov[i][j]*BetaCov[j];
       }
  if (resp[i]==1){
    //y[i]= gsl_ran_gaussian_tail (rr, double a, 1)
  //y[i]=Truncate(xb,1,0,rr);
  y[i]=r_lefttruncnorm(0, xb, 1,rr);
  } else {
    //y[i]=-Truncate(-xb,1,0,rr);
    y[i]=r_righttruncnorm(0, xb, 1,rr);
  }
  //printf("%lf ",y[i]);
  //printf("xb=%lf ",xb);
  } 
  
}


void SampleLambda(gsl_rng * rr,int n,int k,double a0, double b0,double * lambda, double **Xsi){
  int i;
  double sumsq=0;
  for (i=0;i<n;i++){
    sumsq+=pow(Xsi[i][k],2);
  }
  double inv=1/(b0+0.5*(sumsq));
  double  sig2=1/gsl_ran_gamma (rr, a0+n/2.0, inv);
  *lambda=sig2;
}


void SampleSig2(gsl_rng * rr,int n,int * KK, int r,double a0, double b0,double * sigma2, double ***Xsi, double *resp, 
                double **Beta,double h,_Bool ** rho){
  int i,l,k;
double sumsq=0;

  double sig2;
  for (i=0;i<n;i++){
    double xb=0;
    for (l=0;l<r;l++){
        for (k=0;k<KK[l];k++){
          xb+=Xsi[l][i][k]*Beta[l][k];
          
        }}
    sumsq+=pow(resp[i]-xb,2);
  }
  int nr=0;double sumb2=0;
  for (l=0;l<r;l++){
    for (k=0;k<KK[l];k++){
    sumb2+=pow(Beta[l][k],2)/h;
      nr+=rho[l][k];
    }}
 
    double inv=1/(b0+0.5*(sumsq+sumb2));
  sig2=1/gsl_ran_gamma (rr, a0+(n+nr)/2.0, inv);
  *sigma2=sig2;
}

void SampleSig2X(gsl_rng * rr,int n,int K, int TT,double a0,  double b0,double * sigma2, double **Base, double **resp, double **Xsi){
  int i,t,k;
  double sumsq=0;
  double sig2;
  for (i=0;i<n;i++){
    for (t=0;t<TT;t++){
    double xb=0;
      for (k=0;k<K;k++){
        xb+=Xsi[i][k]*Base[t][k];
      }
    sumsq+=pow(resp[i][t]-xb,2);
  }
}
  
  double inv=1/(b0+0.5*sumsq);
  sig2=1/gsl_ran_gamma (rr, a0+n*TT/2.0, inv);
  *sigma2=sig2;
}

void SampleRho(int l,int r,gsl_rng * rr,int n,int p,double *y,double ** X,double sigma2,double** CovImage,double *scal,double *Muscal,_Bool ** rho,int *KK,double w,double **Beta,double h){
  int j;

//void logLikelihood(int n,int p,double *y,double ** X,double sigma2,double Sigmall,double scal, double* quadForm,_Bool * Gam,           double * loggauss)  
  
  double logpostold;double quadForm;
 double scal1=0;
      double Muscal1=0;
  GetH (l,r, rho, &scal1,&Muscal1,Beta, h,CovImage, KK);  
  _Bool Gam[p];
   for (j=0;j<p;j++){
      Gam[j]=rho[l][j];}
  logLikelihood(n,p,y,X,sigma2,CovImage[l][l],scal1,Muscal1,&quadForm,rho[l],&logpostold);
  _Bool Gamnew[p];
  float phi=0.5;
  proposal(p,rho[l],Gamnew,phi, rr);
 for (j=0;j<p;j++){
     rho[l][j]=Gamnew[j];
 }
double scal2=0;
       double Muscal2=0;
   GetH (l,r, rho, &scal2,&Muscal2,Beta, h,CovImage, KK);
  //printf("\nSCAL2=%lf\n",scal2);
   double logpostnew;
  logLikelihood(n,p,y,X,sigma2,CovImage[l][l],scal2,Muscal2,&quadForm,rho[l],&logpostnew);
  for (j=0;j<p;j++){
    logpostold+=Gam[j]*log(w)+(1-Gam[j])*log(1-w);
    logpostnew+=Gamnew[j]*log(w)+(1-Gamnew[j])*log(1-w);
  }
  double u=gsl_ran_flat (rr, 0, 1);
  if ((log(u)>logpostnew-logpostold)){
    for (j=0;j<p;j++){
      rho[l][j]=Gam[j];
      *scal=scal1;*Muscal=Muscal1;  }
  } else {
        *scal=scal2;*Muscal=Muscal2;
    }
  }

void SampleRho1(gsl_rng * rr,int n,int p,double *y,double ** X,double sigma2,double Sigmall,double scal,double Muscal,_Bool * Gam,double w){
  int j;

//void logLikelihood(int n,int p,double *y,double ** X,double sigma2,double Sigmall,double scal, double* quadForm,_Bool * Gam,           double * loggauss)  
  
  double logpostold;double quadForm;
  
  //double scal=0;
  logLikelihood(n,p,y,X,sigma2,Sigmall,scal,Muscal,&quadForm,Gam,&logpostold);
  _Bool Gamnew[p];
  float phi=0.5;
  proposal(p,Gam,Gamnew,phi, rr);
  double logpostnew;
  logLikelihood(n,p,y,X,sigma2,Sigmall,scal,Muscal,&quadForm,Gamnew,&logpostnew);
  for (j=0;j<p;j++){
    logpostold+=Gam[j]*log(w)+(1-Gam[j])*log(1-w);
    logpostnew+=Gamnew[j]*log(w)+(1-Gamnew[j])*log(1-w);
  }
  double u=gsl_ran_flat (rr, 0, 1);
  if ((log(u)<logpostnew-logpostold)){
    for (j=0;j<p;j++){
      Gam[j]=Gamnew[j];
    }
  }
}

  void logLikelihood(int n,int p,double *y,double ** X,double sigma2,double Sigmall,double scal,double Muscal, double* quadForm,_Bool * Gam,double * loggauss){
      int i,s,s1;
  int NZ1[p];
  int nx1=0;
  findc(p,Gam,0,NZ1, &nx1);
  double scalinv=scal/(Sigmall*(Sigmall-nx1*scal));

 // printf("\nSCAL=%lf\n",scal);
 // printf("\nSCALINV=%lf\n",scalinv);  
  // printf("\nSigmall=%lf\n\n",Sigmall);
  double result;double quadF=0;
  if (nx1>0){
    gsl_vector *work1 =gsl_vector_alloc (nx1);
    double * Sigma1=malloc(nx1*nx1*sizeof(double));
    for (s=0;s<nx1;s++){
      for (s1=0;s1<=s;s1++){
        double a=0;
        for (i=0;i<n;i++){
          a+=X[i][NZ1[s]]*X[i][NZ1[s1]];
        }
        Sigma1[s*nx1+s1]=Sigma1[s1*nx1+s]=a+scalinv;
        
      }
      Sigma1[s*nx1+s]+=(1.0/Sigmall);
      //printf("q1=%lf\n",Tau[NZ1[s]][j]);
    }
 gsl_matrix_view m1  = gsl_matrix_view_array (Sigma1, nx1,nx1);

    gsl_linalg_cholesky_decomp (&m1.matrix);
    gsl_vector *yi =gsl_vector_alloc (n);
    
   double muu[n];
    for (i = 0; i < n; i++)
    {
muu[i]=0;
for (s=0;s<nx1;s++){
    muu[i]+=X[i][NZ1[s]]*Muscal;
}
      gsl_vector_set (yi, i, y[i]-muu[i]);
    }
    
    gsl_vector *yy =gsl_vector_alloc (nx1);
    double sumT=0;
    sumT=nx1*log(Sigmall)+log(1-nx1*scal/Sigmall);;
    for (s=0;s<nx1;s++){
      double a=0;
      //if (Gam[NZ1[s]]==1){
      //sumT+=log(h);
      for (i=0;i<n;i++){
        //a+=X[i][NZ1[s]]*y[i]/sqrt(sigma2);
        a+=X[i][NZ1[s]]*(y[i]-muu[i]);
      }
      //}
      gsl_vector_set (yy, s, a);
    
    }
    
    logmultigaussian(yy,yi,&m1.matrix,sigma2,sumT,&result,&quadF,work1);
    //logmultigaussianT(yi, yy,  &m1.matrix,&result,&quadF, work1);
    gsl_vector_free (yy);gsl_vector_free (work1);gsl_vector_free (yi);
    free(Sigma1);
  } else {
    for (i = 0; i < n; i++){
      quadF+=pow(y[i],2);
    }
    quadF=quadF/sigma2;
    result=-(n/2.0)*log(sigma2)- 0.5*n*log(2.0*M_PI)-0.5*quadF;
  }
  *quadForm=quadF;
  *loggauss=result;
  
}
 
 
void SampleRegressionCoeff( const gsl_rng * r,int n, int p, double *y, double ** X, double *beta,double *MuBetaPrior, double *Mu,
                           double* sigma2, double *SigmaPriorBeta){
  int s,s1,i;
  double * InvSigma=malloc(p*p*sizeof(double));
  for (s=0;s<p;s++){
    for (s1=0;s1<=s;s1++){
      double a=0;
      for (i=0;i<n;i++){
        a+=X[i][s]*X[i][s1]/sigma2[i];
      }
      //a+=SigmaPriorBeta[s][s1];
      InvSigma[s*p+s1]=InvSigma[s1*p+s]=a+SigmaPriorBeta[s1*p+s];
    }
   // InvSigma[s*p+s]+=1/(SigmaPriorBeta[s]);
  }

  //printf("\n");
  //for (s=0;s<p*p;s++) printf("%lf ",SigmaPriorBeta[s]);
  //printf("\n\n");

  gsl_matrix_view m  = gsl_matrix_view_array (InvSigma,p,p);
  gsl_linalg_cholesky_decomp (&m.matrix);
  gsl_vector *mu =gsl_vector_alloc (p);
  double xy[p];
  for (s=0; s<p;s++){
    xy[s]=0;
    for (i=0; i<n;i++){
      xy[s]+=X[i][s]*y[i]/sigma2[i];
      }
    //printf("y=%lf ",xy[s]); 
    for (s1=0;s1<p;s1++){
xy[s]+=SigmaPriorBeta[s1*p+s]*MuBetaPrior[s1];
  }
    
    }
 // for (i=0; i<n;i++){
  //  printf("y=%lf ",y[i]); 
 // }
 // printf("\n");
  gsl_vector_view b= gsl_vector_view_array (xy, p);
  gsl_linalg_cholesky_solve (&m.matrix, &b.vector, mu);
  //printf("Mu\n");
  gsl_vector *result=gsl_vector_alloc (p);
  multivariate_gaussian (r, mu, &m.matrix,result);
  for (s=0; s<p;s++){
    beta[s]=gsl_vector_get(result,s);
    Mu[s]=gsl_vector_get(mu,s);
  }
  free(InvSigma); gsl_vector_free (mu);gsl_vector_free (result);
}




int multivariate_gaussian (const gsl_rng * r,
                           const gsl_vector * mu,
                           const gsl_matrix * L,
                           gsl_vector * result)
{
  // This function generates a multivariate distribution with input the 
  //cholesky decompostion of the inverce cova matrix
  /*
  * L     matrix resulting from the Cholesky decomposition of
  *  *      the inverse of   variance-covariance matrix Sigma^-1 = L L^T (dimension d x d)
  */
const size_t M = L->size1;
const size_t N = L->size2;

if (M != N)
{
  GSL_ERROR("requires square matrix", GSL_ENOTSQR);
}
else if (mu->size != M)
{
  GSL_ERROR("incompatible dimension of mean vector with variance-covariance matrix", GSL_EBADLEN);
}
else if (result->size != M)
{
  GSL_ERROR("incompatible dimension of result vector", GSL_EBADLEN);
}
else
{
  size_t i;
  
  for (i = 0; i < M; ++i)
    gsl_vector_set(result, i, gsl_ran_ugaussian(r));
  
  gsl_blas_dtrsv(CblasLower, CblasTrans, CblasNonUnit, L, result);
  //gsl_blas_dtrmv(CblasLower, CblasNoTrans, CblasNonUnit, L, result);
  gsl_vector_add(result, mu);
  
  return GSL_SUCCESS;
}
}


void nrerror(char error_text[])
{
  printf("Utils run-time error...\n");
  printf("%s\n",error_text);
  printf("...now exiting to system...\n");
  exit(1);
}
 double **dmatrix(int nrl, int nrh, int ncl, int nch)
 {
   int i;
   double **m;
   
   m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
   if (!m) nrerror("allocation failure 1 in dmatrix()");
   m -= nrl;
   
   for(i=nrl;i<=nrh;i++)
   {
     m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
     if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
     m[i] -= ncl;
   }
   return m;
 }

void free_bmatrix(_Bool **m, int nrl, int nrh, int ncl, int nch)
{
  int i;
  
  for(i=nrh;i>=nrl;i--) free((_Bool*) (m[i]+ncl));
  
  free((_Bool*) (m+nrl));
}
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch)
{
  int i;
  
  for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
  
  free((char*) (m+nrl));
}

void findc(int n,_Bool R[n],int a,int * IDX, int *nx)
{
  int ii_data[n];
  int idx = 0;
  int  ii = 0;
  _Bool  exitg2 = 0;
  _Bool guard2=0;
  while ((exitg2 == 0) && (ii < n)) {
    guard2 = 0;
    if (R[ii]!= a) {
      idx++;
      ii_data[idx - 1] = ii;
      if (idx >= n) {
        exitg2 = 1;
      } else {
        guard2 = 1;
      }
    } else {
      guard2 = 1;
    }
    
    if (guard2 == 1) {
      ii++;
    }
  }
  
  int loop_ub=idx;
  for (idx = 0; idx < loop_ub; idx++) {
    IDX[idx] = ii_data[idx];
  }
  *nx=loop_ub;
  
}

int logmultigaussian(const gsl_vector * Ax,
                     const gsl_vector * x,
                     const gsl_matrix * L,
                     double sigma2,double logDetD,
                     double * result,double *quadF,
                     gsl_vector * work){
  /* 
   * I_n+A*D*t(A) is the covariance matrix 
   * inv(I_n+A*D*t(A))=I_n-A*inv(inv(D)+t(A)*A)*t(A)
   * The covariance matrix is actually sigma2*(I_n+A*D*t(A))
   * and L is the lower triangular matrix obtained from the cholesky decomposition of inv(D)+t(A)*A=LL^t.
   * 
   * In particular D=Diag(h)
   */
  
  const size_t M = L->size1;
  const size_t N = L->size2;
  if (M != N)
  {
    GSL_ERROR("requires square matrix", GSL_ENOTSQR);
  }
  else if (Ax->size != M)
  {
    GSL_ERROR("incompatible dimension of quantile vector", GSL_EBADLEN);
  }
  else if (work->size != M)
  {
    GSL_ERROR("incompatible dimension of work vector", GSL_EBADLEN);

}
else
{
  size_t i;
  double quadForm;        /* (x - mu)' Sigma^{-1} (x - mu)  with mu=0 but without sigma2*/
  double logSqrtDetSigma; /* log [ sqrt(|Sigma|) ] */
  
  /* compute: work = x - mu with mu=0*/
  for (i = 0; i < M; ++i)
  {
    double xi = gsl_vector_get(Ax, i);
   // double mui = gsl_vector_get(mu, i);
    gsl_vector_set(work, i, xi);
  }
  
  /* compute: work = L^{-1} * (Ax) */
  gsl_blas_dtrsv(CblasLower, CblasNoTrans, CblasNonUnit, L, work);
  
  /* compute: quadForm = (Ax)' Sigma^{-1} (Ax) */
  gsl_blas_ddot(work, work, &quadForm);
  /* compute: quadForm = (x)' (x) */
  
  double quadForm2;gsl_blas_ddot(x, x, &quadForm2);
  quadForm=quadForm2-quadForm; 
    
  /* compute: log [ sqrt(|Sigma|) ] = sum_i log L_{ii} Note that the 
   * |I_n+A*D*t(A)|=|inv(D)+t(A)*A|*|D|
   */
  logSqrtDetSigma = 0.0;
  for (i = 0; i < M; ++i)
  {
    double Lii = gsl_matrix_get(L, i, i);
    logSqrtDetSigma += log(Lii);
  }
  logSqrtDetSigma+=0.5*logDetD;
  *quadF=quadForm;
  *result = -(M/2.0)*log(sigma2)-0.5*quadForm/sigma2 - logSqrtDetSigma - 0.5*M*log(2.0*M_PI);
  return GSL_SUCCESS;
}
}
void proposal(int n,_Bool R[n],_Bool R1[n],float phi, gsl_rng * r)
{
  int i;
  for (i=0; i<n;i++)
    R1[i]=R[i];
  int n1=0;//number of 1
  int n0=0;//number of zeros
  int v1[n];
  findc(n,R,1,v1,&n0);//find indices different from 1 i.e ==0;
  int v2[n-n0];
  findc(n,R,0,v2,&n1);//find indices different of zeros i.e ==1
  double u=gsl_ran_flat (r, 0, 1);
  //printf("No==%d, N1==%d\n",n0,n1);
  if ((u < phi) || (n0 == 0) || (n1 == 0)) {
    int l= gsl_rng_uniform_int(r,n);
    
    R1[l] = 1 - R[l];
  } else {
    int l1=gsl_rng_uniform_int(r,n0);
    int l2=gsl_rng_uniform_int(r,n1);
    
    R1[v1[l1]] = R[v2[l2]];
    R1[v2[l2]] = R[v1[l1]];
  }
}


double Truncate(double mu,double sd, double lower, const gsl_rng * r){
  double lowern=(lower-mu)/sd;
  double alphaopt=(lowern+sqrt(pow(lowern,2)+4))/2;
  double z=lowern+gsl_ran_exponential (r, 1/alphaopt);
  double qz=exp(-pow(z-alphaopt,2)/2);
  double u=gsl_ran_flat (r, 0,1);
  int i=0;
  while ((u>qz)&& (i<10)){
    z=lowern+gsl_ran_exponential (r, 1/alphaopt);
    qz=exp(-pow(z-alphaopt,2)/2);
    u=gsl_ran_flat (r, 0,1);
    i++;
  }
  return z*sd+mu;
}



/* Exponential rejection sampling (a,inf) 
double ers_a_inf(double a,const gsl_rng * r) {
  //SAMPLER_DEBUG("ers_a_inf", a, R_PosInf);
  const double ainv = 1.0 / a;
  double x, rho;
  do {
    x=gsl_ran_exponential (r, ainv)+a;
    rho = exp(-0.5 * pow((x - a), 2));
  } while (gsl_ran_flat (r, 0,1) > rho);
  return x;
}
*/
/* Normal rejection sampling (a,inf) */
double nrs_a_inf(double a,const gsl_rng * r) {
  //SAMPLER_DEBUG("nrs_a_inf", a, R_PosInf);
  //double x = -DBL_AX;
  double x = gsl_ran_ugaussian(r);
  while (x < a) {
    x = gsl_ran_ugaussian(r);
  }
  return x;
}


double r_lefttruncnorm(double a, double mean, double sd,const gsl_rng * r) {
  const double alpha = (a - mean) / sd;
  if (alpha < t4) {
    return mean + sd * nrs_a_inf(alpha,r);
  } else {
    return mean + sd * ers_a_inf(alpha,r);
  }
}
double r_righttruncnorm(double b, double mean, double sd,const gsl_rng * r) {
  const double beta = (b - mean) / sd;
  /* Exploit symmetry: */
  return mean - sd * r_lefttruncnorm(-beta, 0.0, 1.0,r);
}

_Bool **bmatrix(int nrl, int nrh, int ncl, int nch)
{
  int i;
  _Bool **m;
  
  m=(_Bool **) malloc( (nrh-nrl+1)*sizeof(_Bool*));
  if (!m) nrerror("allocation failure 1 in dmatrix()");
  m -= nrl;
  
  for(i=nrl;i<=nrh;i++)
  {
    m[i]=(_Bool *) malloc((nch-ncl+1)*sizeof(_Bool));
    if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
    m[i] -= ncl;
  }
  return m;
}

void sort(int n,double *x,int *idx)
{
  int i,j;
  double a;
  int id;
  for (i = 0; i < n; i++)
    idx[i]=i;
  for (i = 0; i < n; ++i)
  {
    for (j = i + 1; j < n; ++j)
    {
      if (x[i] <= x[j])
      {
        a =  x[i];
        id=idx[i];
        idx[i]=idx[j];
        x[i] = x[j];
        idx[j]=id;
        x[j] = a;
      }
    }
  }
  
}

double mean(int n,double * x){
  int i;
  double me=0;
  for (i=0;i<n;i++)
    me+=x[i];
  return me/n;
}

double max(int n, double * x){
double xmax=x[0];
int i;
for (i=0;i<n;i++){
if (x[i]>xmax)
xmax=x[i];
}
return xmax;
}
  
  

