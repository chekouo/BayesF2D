BayesF2DMCMC<-function(covbCov=TRUE,SamplingGamma=TRUE,y=y,typeoutcome="binary", TwoDX=TwoDX,Xcov=NULL,nbriter=1000,nbrburnin=100,chainNber=1,hypsigm_rr=c(0.5,0.5),h=1){
  if (length(hypsigm_rr)!=2){
stop("hypsigm_rr must be a vector of dimension 2, hyperparameters of sigma_rr")
}

  if (typeoutcome=="binary"){typeoutcome=1}
  if (typeoutcome=="continuous"){typeoutcome=2}


n=length(y);
Dim=dim(TwoDX);
  r=Dim[1];
Ssym=0;
for (l in 1:r){
    Ssym=Ssym+sum(apply(TwoDX[l,,,],1,isSymmetric))
}
if (Ssym==0){
nsym=1;
print("Our 2D-functional features (e.g. GLCM) are non-symmetric")
} else { 
nsym=0;
print("Our 2D-functional features (e.g. GLCM) are symmetric")
}

TwoDX11=Centered(X=TwoDX);
TwoDX=TwoDX11$X;
Bas=Base(X=TwoDX,nsym)
NbrComp=Bas$NbrComp;Base2D=Bas$Base;

#print(paste("The number of eigen values functions is", NbrComp))

if (nsym==1){
  TT=Dim[3]*Dim[4];
  X=as.vector(aperm(TwoDX));
  TNbrComp=NbrComp[,1]*NbrComp[,2];
} else {
XX=aperm(apply(TwoDX,c(2,1),function(y) as.vector(y[lower.tri(y,diag = T)])))
TT=dim(XX)[3];
X=as.vector(aperm(XX));
TNbrComp=NbrComp[,1]*(NbrComp[,2]+1)/2
}
 # KK=as.vector(t(NbrComp));
KK=as.vector(TNbrComp);
dimbeta=nbriter*sum(TNbrComp)
  dimXsi=n*sum(TNbrComp)
  MeanCov=NULL;
  if (!is.null(Xcov)){
  MeanCov=apply(Xcov,2,mean);
  Xcov=scale(Xcov,center = TRUE,scale=FALSE)
  }
  nbrcov=ncol(cbind(rep(1,n),Xcov));
  Xcov1=as.vector(t(cbind(rep(1,n),Xcov)));
  Mean=TwoDX11$Mean;

result<-.C("MCMC",typeOutcome=as.integer(typeoutcome),n1=as.integer(n),
            TT1=as.integer(TT),r1=as.integer(r),nbrcov1=as.integer(nbrcov),resp=as.double(y),X=as.double(X),KK=as.integer(KK),Base1=as.double(Base2D),Xcov1=as.double(Xcov1),
         nbrsample1=as.integer(nbriter),
           burninsample1=as.integer(nbrburnin),
           Betasample=as.double(rep(0,dimbeta)),XsiMean1=as.double(rep(0,dimXsi)),
           rhoMean1=as.double(rep(0,sum(TNbrComp))),
           CovImageSample=as.double(rep(0,nbriter*r*r)),LambdaMean1=as.double(rep(0,sum(TNbrComp))),seed1=as.double(chainNber),hypsigm_rr=as.double(hypsigm_rr),h1=as.double(h),covbool=as.integer(1*covbCov),rhosample1=as.integer(rep(0,dimbeta)),yMean=as.double(rep(0,n)),sigma2xMean=as.double(rep(0,r)),sigma2yMean1=as.double(0),as.integer(1*SamplingGamma))
  BetaSample=list()
  BetaMedian=list()
  Beta2.5=list();
Beta97.5=list()
  XsiSampleMean=list()
  GamMean=list()
  LambdaMean=list()
  bs=0;xs=0;
  for (l in 1:r){
  BetaSample[[l]]=matrix(0,nbriter,TNbrComp[l])
  BetaMedian[[l]]=rep(0,TNbrComp[l])
   Beta2.5[[l]]=rep(0,TNbrComp[l])
  Beta97.5[[l]]=rep(0,TNbrComp[l])
  GamMean[[l]]=rep(0,TNbrComp[l])
  }
  for (s in 1:nbriter){
  for (l in 1:r){
    for (k in 1:(TNbrComp[l])){
    bs=bs+1;
    BetaSample[[l]][s,k]=result$Betasample[bs]
    for (i in 1:n){
      xs=xs+1
    }
    }
  }
  }

  XsiMean=result$XsiMean1
  xx=0;xr=0;
  BetaFunc=list();kkk=1;
  T1=Dim[3]
  for (l in 1:r){
    BetaFunc[[l]]=matrix(0,T1,T1)
    b=n*TNbrComp[l];
    XsiSampleMean[[l]]=matrix(XsiMean[(1+xx):(xx+b)],TNbrComp[l],byrow = T)
    xx=xx+b;
    BetaMedian[[l]]=apply(BetaSample[[l]],2,median)
     Beta2.5[[l]]=apply(BetaSample[[l]],2,function (x) quantile(x, probs = 0.025))
     Beta97.5[[l]]=apply(BetaSample[[l]],2,function (x) quantile(x, probs = 0.975))
    GamMean[[l]]=result$rhoMean1[(1+xr):(xr+TNbrComp[l])]
    LambdaMean[[l]]=result$LambdaMean1[(1+xr):(xr+TNbrComp[l])]
  xr=xr+TNbrComp[l];
  for (s in 1:T1){
    if (nsym==1){xy=1} else {xy=s}
    for (t in xy:T1){
      kk=1;
      for (k1 in 1:NbrComp[l,1]){
        if (nsym==1){xx=1} else {xx=k1}
        for (k2 in xx:NbrComp[l,2]){
      BetaFunc[[l]][s,t]= BetaFunc[[l]][s,t]+BetaMedian[[l]][kk]*Base2D[kkk];
      if  (nsym==0){BetaFunc[[l]][t,s]=BetaFunc[[l]][s,t]}
      kkk=kkk+1;
      kk=kk+1;
      }
    }
  }
  }
  }
CovSample1=aperm(array(result$CovImageSample,dim=c(r,r,nbriter)))
CovImageMean=apply(CovSample1,c(2,3),mean);
CovImageMean1=as.vector(t(apply(CovSample1,c(2,3),mean)));
  return (list(NbrComp=NbrComp,BetaMedian=BetaMedian,BetaFunc=BetaFunc,Beta2.5=Beta2.5,Beta97.5=Beta97.5,XsiSampleMean=XsiSampleMean,GamMean=GamMean,LambdaMean=LambdaMean,CovImageMean=CovImageMean,CovSample=CovSample1,
 ForPredict=list(nsym=nsym,MeanCov=MeanCov,n=n,Mean=Mean,typeoutcome=typeoutcome,TT=TT,r=r,nbrcov=nbrcov,y=y,KK=KK,Base2D=Base2D,Xcov1=Xcov1,nbriter=nbriter,XsiMean1=result$XsiMean1,CovImageMean1=CovImageMean1,rhosample1=result$rhosample1,LambdaMean1=result$LambdaMean1,hypsigm_rr=hypsigm_rr,h=h,covboolean=covbCov*1,yMean=result$yMean,X=X,sigma2xMean=result$sigma2xMean,sigma2yMean1=result$sigma2yMean1)))
}

