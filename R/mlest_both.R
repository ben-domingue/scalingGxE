##note this has lambda0/1/2




mlest_2lambda<-function(df,
                start.pars=list(
                    m=rnorm(1,sd=.05),
                    pi0=rnorm(1,sd=.05),
                    pi1=rnorm(1,sd=.05),
                    lambda0=1+runif(1),
                    lambda1=runif(1,max=.2),
                    lambda2=runif(1,max=.2)
                ),
                #hess=FALSE, #if FALSE, compute empirical hessian
                ctl.list=list(maxit=1000,reltol=sqrt(.Machine$double.eps)/1000)
                ) {
    ll4<-function(pars,df,pr) { #this is the likelihood function to be maximized 
        for (nm in names(pars)) assign(nm,pars[[nm]])
        mu<-m*df$E+
            pi0*df$g+
            pi1*df$g*df$E
        sig<-sqrt((lambda0+df$E*lambda1+df$g*lambda2)^2)
        ##
        d<-dnorm(df$y,mean=mu,sd=sig,log=TRUE) 
        tr<- -1*sum(d)
        #print(pars)
        #print(tr)
        tr
    }
    gr4<-function(pars,df) { #this is the gradient
        for (nm in names(pars)) assign(nm,pars[[nm]])
        mu<-m*df$E+
            pi0*df$g+
            pi1*df$g*df$E
        sig<-sqrt((lambda0+df$E*lambda1+df$g*lambda2)^2)
        ##
        dfdmu<-(df$y-mu)*sig^(-2)
        dfdm<-sum(dfdmu*df$E)
        dfdpi0<-sum(dfdmu*df$g)
        dfdpi1<-sum(dfdmu*df$g*df$E)
        dfdsigma <- (df$y - mu)^2*sig^(-3) - sig^(-1)
        dfdlambda0 <- sum(dfdsigma)
        dfdlambda1 <- sum(dfdsigma*df$E)
        dfdlambda2 <- sum(dfdsigma*df$g)
        -1*c(dfdm,dfdpi0,dfdpi1,dfdlambda0,dfdlambda1,dfdlambda2)
    }
    ##
    fit<-optim(par=start.pars,ll4,
               gr=gr4,
               df=df,
               control=ctl.list,
               hessian=TRUE
               )
    est<-c(fit$par)
    #if (hess) { #see https://stats.stackexchange.com/questions/27033/in-r-given-an-output-from-optim-with-a-hessian-matrix-how-to-calculate-paramet
    vc<-solve(fit$hessian) #note, no -1 given that i am doing that directly in the above
    return(list(est=est,var=vc))
}
