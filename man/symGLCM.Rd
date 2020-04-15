\name{symGLCM}
\alias{symGLCM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simulate symmetric GLCMs.
}
\description{
Simulate symmetric GLCMs which are correlated across imaging modalities for one subject
}
\usage{
symGLCM(m, S, Gam, noise, g_size, count)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{m}{
List of 4 2x1 means bivariate indicating the highest intensity coordinates in GLCMs
}
  \item{S}{
List of 4 2x2 covariance matrices of the bivariate normal
}
  \item{Gam}{
Gam list of 64 4x4 correlation matrices
}
  \item{noise}{
 noise.
}
  \item{g_size }{
Number of gray levels.
}
  \item{count }{
Number of pixels in the image.
}
}
\details{
Return the GLCMs as a list. The  process to generate data are explained in the reference paper. 
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
\code{\link{simulateGLCM_binary}}
}
\examples{
library(BayesF2D)
}
% add one or more standard keywords, see file 'keywords' in the
% r documentation directory.
%\keyword{ ~kwd1 }% use one of  rshowdoc("keywords")
%\keyword{ ~kwd2 }% __only one__ keyword per line
