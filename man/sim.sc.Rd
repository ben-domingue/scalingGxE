\name{sim.sc}
\alias{sim.sc}
\title{Simulate data consistent with scaling model}
\description{
}
\usage{
sim.sc(i=1,E,b0=.8,b1=.2,h=sqrt(.6),a=.5)
}
\arguments{
}
\references{
}
\author{Ben Domingue \email{ben.domingue@gmail.com}} 
\examples{
set.seed(8675309)
E<-function() rnorm(5000)
df<-sim.sc(E=E)
}

