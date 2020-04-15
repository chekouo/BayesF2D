Prediction<-function(Bayes2DObject,prediction="CV",K=NULL,TwoDXtest=NULL, Xcovtest=NULL){
Bayes2DObject=Bayes2DObject$ForPredict;
n=Bayes2DObject$n;
r=Bayes2DObject$r;

 if (prediction=="CV")
  {
      predict=1;
      if (is.null(K)){
          stop("You must specify K, the number of folds for cross validation")
      }
##library(caret);
TestList <- createFolds(1:n, k = K, list = TRUE, returnTrain = FALSE);
nbrset=unlist(lapply(TestList,length))
ListTest=unlist(TestList);
ListTraining=unlist(lapply(TestList,function(x) seq(1,n)[-x]))
  } else {
  K=1;
  predict=0;
  nbrset=1;
  ListTraining=1:n;
  }
if (predict==0){
    if (is.null(TwoDXtest)){
        stop("You must specify the new testing data for imaging sequence for cross validation")
    } else {
if (length(Bayes2DObject$ForPredict$Xcov1)>n){       
        if (is.null(Xcovtest)){
 stop("You must specify the new testing data for the set of covarites (e.g age, gender) for cross validation")
        }
      }
    }

}

nsym=Bayes2DObject$nsym
MeanCov=Bayes2DObject$MeanCov; 
Xcov1=Bayes2DObject$Xcov1; 
X=Bayes2DObject$X;
if (predict==1){
 Xtest=X;
  Xcovtest1=Xcov1;
  ntest=0;
  nr=n;
  } else if (!is.null(TwoDXtest)){
  ntest=dim(TwoDXtest)[2]
 ListTest=1:ntest;
  Ssym1=0;
 for (l in 1:r){
     Ssym1=Ssym1+sum(apply(TwoDXtest[l,,,],1,isSymmetric))
 }
 if ((nsym==0 & Ssym1==0)|| (nsym==1 & Ssym1>0)){
      stop("Both the training data and the test data should be (non)symmetric")
 }
Mean=Bayes2DObject$Mean;
TwoDXtest=CenteredTest(Xtest=TwoDXtest,Mean=Mean);
    if (!is.null(Xcovtest)){
      Xcovtest=t(t(Xcovtest)-MeanCov)
        if (nrow(Xcovtest)!=ntest){
        stop("The number of test samples must be the same both in the clinical and the imaging data")
        }
    }
    if (nsym==1){
    Xtest=as.vector(aperm(TwoDXtest));
    } else {
    XXtest=aperm(apply(TwoDXtest,c(2,1),function(y) as.vector(y[lower.tri(y,diag = T)])))
    Xtest=as.vector(aperm(XXtest));
    }
    Xcovtest1=as.vector(t(cbind(rep(1,ntest),Xcovtest)));
    nr=ntest;
  }



seed1=1;
typeoutcome=Bayes2DObject$typeOutcome;TT=Bayes2DObject$TT;
nbrcov=Bayes2DObject$nbrcov;y=Bayes2DObject$y;
KK=Bayes2DObject$KK;Base2D=Bayes2DObject$Base2D;
nbriter=Bayes2DObject$nbriter;
XsiMean1=Bayes2DObject$XsiMean1;CovImageMean1=Bayes2DObject$CovImageMean1;
rhosample1=Bayes2DObject$rhosample1;LambdaMean1=Bayes2DObject$LambdaMean1;
hypsigm_rr=Bayes2DObject$hypsigm_rr;h=Bayes2DObject$h;covboolean=Bayes2DObject$covboolean
yMean=Bayes2DObject$yMean; sigma2xMean=Bayes2DObject$sigma2xMean;
sigma2yMean1=Bayes2DObject$sigma2yMean1;
result<-.C("Prediction",typeOutcome=as.integer(typeoutcome),n1=as.integer(n),
           ntest1=as.integer(nr), TT1=as.integer(TT),r1=as.integer(r),
           nbrcov1=as.integer(nbrcov),yMean=as.double(yMean),resp=as.double(y),
           X=as.double(X),
           Xtest=as.double(Xtest),KK=as.integer(KK),Base1=as.double(Base2D),
           Xcov1=as.double(Xcov1),Xcovtest1=as.double(Xcovtest1),Kf=as.integer(K),
           nbrset=as.integer(nbrset), ListTest=as.integer(ListTest),
           ListTraining=as.integer(ListTraining),nbrsample1=as.integer(nbriter),
           Predict1=as.integer(predict),XsiMean1=as.double(XsiMean1),
           prob=as.double(rep(0,nr)),CovImageMean1=as.double(CovImageMean1),
           rhosample1=as.integer(rhosample1), LambdaMean1=as.double(LambdaMean1),
           seed1=as.double(seed1),hypsigm_rr=as.double(hypsigm_rr),h1=as.double(h),
           covbool=as.integer(covboolean),sigma2xMean=as.double(sigma2xMean),sigma2yMean1=as.double(sigma2yMean1))

return (list(predprob=result$prob))
}

