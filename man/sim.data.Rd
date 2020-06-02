\name{sim.data}
\alias{sim.data}
\title{Simulate GxE data}
\description{
    Simulates data corresponding to either the scaling or vanilla GxE model. 
    Assumes G ~ N(0,1)
}
\usage{
sim.data(E,scaling=TRUE,b0=.8,b1=.2,b2=0,b3=.05,h=sqrt(.6),a=.5,sigma=1)
}
\arguments{
    \item{E}{User-supplied E values}
    \item{scaling}{Controls whether data is generated from scaling model or vanilla GxE}
    \item{b0}{Only used in scaling}
    \item{b1}{Used in both scaling and vanilla}
    \item{b2}{Only used in vanilla}
    \item{b3}{Only used in vanilla}
    \item{h}{Only used in scaling}
    \item{a}{Only used in scaling}
    \item{sigma}{Error variance, used in both scaling and vanilla}
}
\references{
}
\author{Klint Kanopka \email{klint.kanopka@gmail.com}} 
\examples{
set.seed(8675309)
df_s <- sim.data(E=rnorm(5000))
df_v <- sim.data(E=rnorm(5000), scaling=FALSE)
}

