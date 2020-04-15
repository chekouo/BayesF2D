Centered<-function(X=X){
Dim=dim(X);
R=Dim[1];N=Dim[2];
Mean=array(0,dim=c(R,Dim[3],Dim[4]))
for (r in 1:R){
  for (n in 1:N){
    X[r,n,,]=log(X[r,n,,]+0.5)
  }
  Mean[r,,]=apply(X[r,,,],c(2,3),mean)
  #Sd=apply(X[r,,,],c(2,3),sd)
  for (n in 1:N){
    X[r,n,,]=(X[r,n,,]-Mean[r,,])
  }
}
list(X=X,Mean=Mean)
}
 CenteredTest<-function(Xtest=Xtest,Mean=Mean){
 Dim=dim(Xtest);
 R=Dim[1];N=Dim[2];
 for (r in 1:R){
   for (n in 1:N){
     Xtest[r,n,,]=(log(Xtest[r,n,,]+0.5)-Mean[r,,])
   }
 }
 Xtest
 }


library(fdapace)
require("mgcv")

Xst=function(s,t,mod){
  newd=data.frame(x=s,y=t)
  predict.gam(mod,newd)
}
Xs=Vectorize(function(s=s,mod){
  integrate(function(t) Xst(s,t,mod),lower = 0,upper = TT+1)$value
},"s")
Xt=Vectorize(function(t=t,mod){
  integrate(function(s) Xst(s,t,mod),lower = 0,upper = TT+1)$value
},"t")

#s=1:TT
EigenFunction=function(s=s,XB=XB,nsym=1){
set.seed(1)
  K=rep(0,2)
  M=length(s)
  Dim=dim(XB);N=Dim[1];
  #for (r in 1:R){
  yPred=matrix(0,N,M);yPredT=matrix(0,N,M)
  for (n in 1:N){
    m=XB[n,,]
    df <- data.frame(x = rep(seq_len(ncol(m)), each = nrow(m)),
                     y = rep(seq_len(nrow(m)), times = ncol(m)),
                     z = c(m))
    mod <- gam(z ~ te(x, y), data = df)
    yPred[n,]=Xs(s=s,mod);
    if (nsym==1) {
        yPredT[n,]=Xt(t=s,mod)
    } else {
        yPredT[n,]=yPred[n,];
    }
  }
  L3 <- MakeFPCAInputs(IDs = rep(1:N, each=M), tVec=rep(s,N), t(yPred));
  FPCAdense <- FPCA(L3$Ly, L3$Lt,optns = list(FVEthreshold=0.95))
  EigenFXs=FPCAdense$phi
  K[1]=FPCAdense$selectK;
  if (nsym==1) {
  L3 <- MakeFPCAInputs(IDs = rep(1:N, each=M), tVec=rep(s,N), t(yPredT));
  FPCAdense <- FPCA(L3$Ly, L3$Lt,optns = list(FVEthreshold=0.95))
  EigenFXt=FPCAdense$phi;K[2]=FPCAdense$selectK;
  } else {
      EigenFXt=EigenFXs;K[2]=K[1];
  }
  list(EigenFXs=EigenFXs,EigenFXt=EigenFXt,K=K)
}
Base<-function(X=X,nsym=1){
    Dim=dim(X)
 R=Dim[1];N=Dim[2];T1=Dim[3]
EigenFunctionsXs=list();EigenFunctionsXt=list(); K=matrix(0,R,2)
for (r in 1:R){
  Eig=EigenFunction(s=1:T1,XB=X[r,,,],nsym)
  EigenFunctionsXs[[r]]=Eig$EigenFXs; EigenFunctionsXt[[r]]=Eig$EigenFXt;
  K[r,]=Eig$K
}
Base2D=NULL;
kk=0;
for (r in 1:R){
  for (s in 1:T1){
 if (nsym==1){xy=1} else {xy=s}
      for (t in xy:T1){
   for (k1 in 1:K[r,1]){
   if (nsym==1){xx=1} else {xx=k1}
    for (k2 in xx:K[r,2]){
      kk=kk+1;
      Base2D[kk]=EigenFunctionsXs[[r]][s,k1]*EigenFunctionsXt[[r]][t,k2]
        }
      }
    }
    }
}
list(NbrComp=K,Base=Base2D)
}

