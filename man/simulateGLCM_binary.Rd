\name{simulateGLCM_binary}
\alias{simulateGLCM_binary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate GLCMs for binary response
}
\description{
Simulate symmetric or non-symmetric GLCMs for binary response data.
}
\usage{
simulateGLCM_binary(n0 = 10,n1 = 10, m0, S0, m1, S1, GLCM.type = 'sym', 
Gam = NULL, noise = 1, counts = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{n0}{
Number of samples for y=0 for binary data; defaults to 10.
}
\item{n1}{
Number of samples for y=1 for binary data; defaults to 10.
}
\item{m0}{
List of 4 2x1 means bivariate indicating the highest intensity coordinates in GLCMs for y=0..
}
\item{m1}{
List of 4 2x1 means bivariate indicating the highest intensity coordinates in GLCMs for y=1.
}
\item{S0}{
List of 4 2x2 covariance matrices of the bivariate normal for individuals for y=0
}
\item{S1}{
List of 4 2x2 covariance matrices of the bivariate normal for individuals for y=1
}
\item{GLCM.type}{
It can take either "sym" or "nonsym" for symmateric and non-symmetric GLCMs; defaults to symmetric GLCM
}
\item{Gam}{
List of 4x4 correlation matrices - 64 non 'nonsym' case and 36 for 'sym' case; defaults to covariance matrices provided as data in the package
}
\item{noise}{
 Noise defaults to 1.
}
  \item{counts }{
counts a 'n0+n1' vector for the number of pixels assumed in the image; First n0 belong to y=0 and last n1 belong to y=1
}
}
\details{
Return the GLCMs and the response.
}
\references{
Thierry Chekouo et al (2020), \emph{A Bayesian 2D Functional Linear Model for Gray LevelCo-occurrence Matrices in Texture Analysis of LowerGrade Gliomas}, \emph{submitted}.
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
# Example to generate "non-symmetric" GLCMs (4 correlated sequences for each subject)
# for 20 subjects in two categories (10 in each category). The mean for the four
# imaging sequences for y=0 is assumed to be (1,7), (3,5) (5,3) and (7,1), and for
# y=1 is assumed to be (2,6), (2,2), (5,2) and (7,1). The covariance matrix for the
# bivariate normal is assumed to be [5 -0.5; -0.5 5] for y=0 and [10 -1; -1 10] for y=1.
# The correlation between the imaging sequences is included from the data provided in
# the package (nonsymCor.RData) which is loaded into the variable 'Cor'.

# Similar appraoch could be used for symmetric GLCM. Only changes would be loading the
# appropriate correlation matrix in the first step and specifycing the 'GLCM.type' as 'sym'.
data(nonsymCor)
#Cor = BayesF2D::nonsymCor
Cor = nonsymCor
s.cov0 = 5*matrix(c(1,-0.1,-0.1,1), nrow=2)
S0 = list(s.cov0, s.cov0, s.cov0, s.cov0)

s.cov1 = 10*matrix(c(1,-0.1,-0.1,1), nrow=2)
S1 = list(s.cov1, s.cov1, s.cov1, s.cov1)

glcm = BayesF2D::simulateGLCM_binary(GLCM.type = "nonsym",
                                        m0 = list(c(1,7),c(3,5),c(5,3),c(7,1)),
                                        S0 = S0,
                                        m1 = list(c(2,6),c(2,2),c(5,2),c(7,1)),
                                        S1 = S1,
                                        Gam = Cor,
                                        noise = 3)


glcm


}
% add one or more standard keywords, see file 'keywords' in the
% r documentation directory.
%\keyword{ ~kwd1 }% use one of  rshowdoc("keywords")
%\keyword{ ~kwd2 }% __only one__ keyword per line
