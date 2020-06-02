\name{xi.test}
\alias{xi.test}
\title{Test xi}
\description{
}
\usage{
xi.test(est,p.value=TRUE)
}
\arguments{
\item{est}{Object returned by ratio.test}
\item{p.value}{If true, will return probability associated with test
  statisic. If false, will return CI.}
}
\references{
}
\author{Ben Domingue \email{ben.domingue@gmail.com}} 
\examples{
##example wherein the scaling model is the true model
set.seed(8675309)
E<-function() rnorm(5000)
df<-sim.sc(E=E)
est<-ratio.test(df,hess=FALSE)
xi.test(est)
##note we fail to reject the null
##############################################################
##example wherein vanilla GxE is the true model
simdata<-function(i=1,E,b1=.1,b2=0,b3=.05,sigma=1) {
    E<-E()
    N<-length(E)
    G=rnorm(N,0,1)
    y=b1*G+b2*E+b3*G*E+rnorm(N,sd=sigma)
    df<-data.frame(E=E,y=y,g=G)#,e=eps)
    df
}
E<-function() rnorm(8000)
df<-simdata(E=E)
est<-ratio.test(df,hess=TRUE)
xi.test(est)
##note we reject the null
}

