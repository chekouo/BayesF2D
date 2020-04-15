\name{BayesF2DMCMC}
\alias{BayesF2DMCMC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
An MCMC algorithm to perform a Bayesian  two-dimensional functional linear models for scalar responses.  
}
\description{
The algorithm computes (i) the marginal posterior probabilities of inclusions for each feature (i.e combinations of eigenvectors obtained from each dimension using FPCA), (ii) provide estimates of functional regression effects, and (ii) evaluate prediction performance of the method for both continuous and binary responses.  
}
\usage{
BayesF2DMCMC(covbCov=TRUE,y=y,typeoutcome="binary", TwoDX=TwoDX,Xcov=NULL,nbriter=1000,
nbrburnin=100,chainNber=1,hypsigm_rr=c(0.5,0.5),h=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{covbCov}{
If TRUE, then we assume a positive correlation between the functional regression effects if the feature types (e.g if TRUE, the correlation is assume between the 4 imaging sequences mentioned in the referenced paper). Otherwise, no correlation. 
}
\item{y}{
The response variable which can be binary or continuous. If it's binary, it should be  enter as 1 or  0.
}
\item{typeoutcome}{
Type of the response variable. If the response is binary, then typeoutcome="binary". Otherwise, it's typeoutcome="continuous". 
}
\item{TwoDX}{
Array of features values with dimension \eqn{r\times n\times T\times T} where \eqn{r} is the number of feature types (e.g. the number of imaging sequences), \eqn{n} is the number of individuals, \eqn{T} is the one-dimensional (functional) size (e.g. the number of imaging intensities in the  GLCM matrix). 
}
\item{Xcov}{
Other nonfunctional covariates such as sex, age and other clinical covariates
}
\item{nbriter}{
%%     ~~Describe \code{sample} here~~
Total number of MCMC draws. It must be larger than burnin.
}
\item{nbrburnin }{ 
Number of draws to discard for burn-in
}
\item{chainNber }{
This defines the seed to generate random numbers. It may help to generate multiple MCMC chains. 
}
\item{hypsigm_rr }{
Shape and rate hyperparameters of the prior of the covariance \eqn{\sigma_{rr'}}, covariance of  the regression effects between feature types
}
\item{h}{
Shrinkage parameter for regression coefficients defined  the variance of coefficiens effects.
}
%  \item{TwoDXtest}{
%%     ~~Describe \code{b.tau} here~~
%}
}
\details{
The function will return several R objects, which can be assigned to a variable. To see the results, use the ``$" operator. 
}
\value{
\item{NbrComp}{
Matrix of the number of FPCA (functional PCA) components used for each feature type. It was chosen using a fraction of variance explained (FVE) of at least 0.95. 
}
 \item{BetaMedian}{ This is the median of the effect of each basis function (pair of eigunfucntions) obtained from FPCA}
\item{BetaFunc}{ Estimated functional effects using median estimaes of the effect of each base vector}
\item{XsiMean}{ Posterior mean estimmates of \eqn{\xi^{(r)}_{ikl}}'s,  principal component scores associed with pair of eigenfunctions \eqn{(\phi_k, \psi_l)}'s.
}
\item{GamMean}{ Estimated  marginal  posterior probabilities, \eqn{Prob(\gamma^{(r)}_{kl}=1)} of component selection.} 
%% ...
\item{LambdaMean}{
Posterior mean \eqn{\lambda^{(r)}_{kl}}'s, eigen values associated with each pair  of eigen functions.
}
\item{CovImageMean}{
Posterior means of the covariance matrix  of the effects between the feature types. It's dimensio is \eqn{r\times r}, where \eqn{r} is the number of feature types (e.g. \eqn{r=4}imaging sequences).
}
\item{ForPredict}{
List of objects to perform prediction, using Bayesian model averaging
}
}
\references{
Thierry Chekouo et al (2020), \emph{A Bayesian 2D Functional Linear Model for Gray Level Co-occurrence Matrices in Texture Analysis of Lower Grade Gliomas}, \emph{submitted}
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
\code{\link{simulateGLCM_binary}},\code{\link{Prediction}}
}
\examples{
library(BayesF2D)
data(nonsymCor)
#Cor = BayesF2D::nonsymCor
Cor = nonsymCor
s.cov0 = 5*matrix(c(1,-0.1,-0.1,1), nrow=2)
S0 = list(s.cov0, s.cov0, s.cov0, s.cov0)

s.cov1 = 10*matrix(c(1,-0.1,-0.1,1), nrow=2)
S1 = list(s.cov1, s.cov1, s.cov1, s.cov1)

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
}
% add one or more standard keywords, see file 'keywords' in the
% r documentation directory.
%\keyword{ ~kwd1 }% use one of  rshowdoc("keywords")
%\keyword{ ~kwd2 }% __only one__ keyword per line
