\name{Prediction}
\alias{Prediction}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Prediction performance
}
\description{
This function evaluates the prediction performance using cross-validation techniques or by computing the prediction of new observations. 
}
\usage{
Prediction(Bayes2DObject,prediction="CV",K=NULL,TestList=NULL,TwoDXtest=NULL, Xcovtest=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Bayes2DObject}{
This is Bayes2D object obtained by fitting an MCMC algorithm on the training dataset.   
}
  \item{prediction}{
%%     ~~Describe \code{Ind.Var} here~~
It takes  two values:  "CV" for cross-valiadation technique  or "NoCV" if we predict  new observations. 
}
  \item{K}{
The number of folds if prediction="CV". 
}
  \item{TestList}{
Sample partition as a list if prediction="CV". For instance, if n=5 and K=2, Teslist[[1]]=1,3,5 and Teslist[[1]]=2,4. However, when  TesList=NULL, the partition is done ramdomly.
}
\item{TwoDXtest}{
%%     ~~Describe \code{sample} here~~
New dataset (e.g log-GLCM data) for prediction of dimension \eqn{r\times n\times T\times T} where T is the number of gray levels, r is the number of imaging sequences, and n is the number of individuals. 
}
  \item{Xcovtest }{ 
Matrix of covariates of the new indiduals. Dimension \eqn{n\times nbrcov} where nbrcov is the number of covariates (e.g. clinical covariates).
}
}
\details{
The function will return predictive probabilities if the response is binary or scalar values if the response is continuous. 
}
\references{
Thierry Chekouo et al (2020), \emph{A Bayesian 2D Functional Linear Model for Gray Level Co-occurrence Matrices in Texture Analysis of LowerGrade Gliomas}, \emph{submitted}.
}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{BayesF2DMCMC}}
}
\examples{
library(BayesF2D)
#### Generate data
data(nonsymCor)
#Cor = BayesF2D::nonsymCor
Cor = nonsymCor
s.cov0 = 5*matrix(c(1,-0.1,-0.1,1), nrow=2)
S0 = list(s.cov0, s.cov0, s.cov0, s.cov0)

s.cov1 = 10*matrix(c(1,-0.1,-0.1,1), nrow=2)
S1 = list(s.cov1, s.cov1, s.cov1, s.cov1)
set.seed(1)
glcm = BayesF2D::simulateGLCM_binary(n0=30,n1=20,GLCM.type = "nonsym",
                                        m0 = list(c(1,7),c(3,5),c(5,3),c(7,1)),
                                        S0 = S0,
                                        m1 = list(c(2,6),c(2,2),c(5,2),c(7,1)),
                                        S1 = S1,
                                        Gam = Cor,
                                        noise = 3)

meth1="NoCorBeta"
y=glcm$y ## binary response
N=length(y)
R=4# number of imaging sequences
TT=8 # number of image intensities
X=array(0,dim=c(R,N,TT,TT))
for (i in 1:N){
X[,i,,]=glcm$GLCM[[i]];
}

Res<-BayesF2DMCMC(covbCov=(meth1=="CorBeta"),y=y,typeoutcome="binary",TwoDX=X,nbriter=5000,nbrburnin=1000,chainNber=1,hypsigm_rr=c(0.5,0.5),h=1)
library(AUC)
### Cross-validation prediction
Pred<-Prediction(Bayes2DObject=Res,prediction="CV",K=10,TwoDXtest=NULL, Xcovtest=NULL)
Pred$predprob
au1=auc(roc(Pred$predprob,as.factor(y)));
print(au1)

### Prediction of new individuals
## Generate test data
glcmTest = BayesF2D::simulateGLCM_binary(n0=30,n1=20,GLCM.type = "nonsym",
                                        m0 = list(c(1,7),c(3,5),c(5,3),c(7,1)),
                                        S0 = S0,
                                        m1 = list(c(2,6),c(2,2),c(5,2),c(7,1)),
                                        S1 = S1,
                                        Gam = Cor,
                                        noise = 3)
Xtest=array(0,dim=c(R,N,TT,TT))
for (i in 1:N){
Xtest[,i,,]=glcmTest$GLCM[[i]];
}

Pred<-Prediction(Bayes2DObject=Res,prediction="NoCV",TwoDXtest=Xtest)
Pred$predprob
au1=auc(roc(Pred$predprob,as.factor(glcmTest$y)));
print(au1)


}
% add one or more standard keywords, see file 'keywords' in the
% r documentation directory.
%\keyword{ ~kwd1 }% use one of  rshowdoc("keywords")
%\keyword{ ~kwd2 }% __only one__ keyword per line
