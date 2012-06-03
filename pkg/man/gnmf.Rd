\name{gnmf}
\alias{gnmf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generalized non-negative matrix factorization
}
\description{
Performs generalized non-negative matrix factorization
}
\usage{
gnmf(V, scheme, nsteps = 2000, repeats = 20, ranks = 2, cltarget = "PATTERN",
clscheme = "Binary", reffile = "", scaling = "F", normalizing = "F", alphas = 1,
runtype = "simulation", cstepsize = 20, idealization = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{V}{
Input data matrix
}
  \item{scheme}{
KL, Renyi, or ED
}
  \item{nsteps}{
Update steps, default 2000
}
  \item{repeats}{
Repeats, default 20
}
  \item{ranks}{
The number of components into which matrix \code{V} is to be factored, default 2 (a scalar)
}
  \item{cltarget}{
Clustering target, default 'PATTERN' (\code{H} matrix)
  either PATTERN or ALTERNATE
}
  \item{clscheme}{
Clustering scheme, default 'Binary', could be 'PearsonHC'
}
  \item{reffile}{
Default none
}
  \item{scaling}{
Boolean, default F
}
  \item{normalizing}{
Boolean, H matrix normalization, default 'F'
}
  \item{alphas}{
Renyi parameter, default 1.0 (a scalar), ignored if \code{scheme} is not Renyi
}
  \item{runtype}{
simulation (default) or evaluation or whole
}
  \item{cstepsize}{
Convergence test step size, default 20
}
  \item{idealization}{
Default 1
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
Devarajan K. Nonnegative matrix factorization: an analytical and interpretive 
tool in computational biology. PLoS Comput Biol. 2008 Jul 25;4(7):e1000029.

\code{\url{http://devarajan.fccc.edu}}
}
\author{
Jose M. Maisog, Guoli Wang, Karthik Devarajan
}
\note{
Further notes...
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\link{gnmf-package},
\link{V}
}
\examples{
# Load sample data.
data(V)

# Compute NMF with 20 repeats.
result <- gnmf(V,scheme="KL")

# Extract H and W from the result.
# H and W are lists, each containing the result of 20 repeats.
H <- result$H
W <- result$W

# Get the H and W matrices of the first repeat.
H1 <- H[[1]]
W1 <- W[[1]]
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ Clustering }
\keyword{ Multivariate Techniques }% __ONLY ONE__ keyword per line
