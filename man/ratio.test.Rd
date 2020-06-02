\name{ratio.test}
\alias{ratio.test}
\title{Estimate the heteroscedastic GxE model}
\description{
forthcoming
}
\usage{
ratio.test<-function(df,hess=FALSE,ctl.list=list(maxit=1000,reltol=sqrt(.Machine$double.eps)/1000))
}
\arguments{
\item{df}{data frame.}
\item{hess}{Analytic hessian (FALSE) or empirical hessian (TRUE).}
\item{ctl.list}{Options to pass to optim}
}
\references{
forthcoming
}
\author{Ben Domingue \email{ben.domingue@gmail.com}} 
\examples{
set.seed(8675309)
E<-function() rnorm(5000)
df<-sim.sc(E=E)
est<-ratio.test(df,hess=FALSE)
}

