
mlest<-function(df,
                hess=FALSE, #if FALSE, compute empirical hessian
                covars=NULL, #names of covariates
                ctl.list=list(maxit=1000,reltol=sqrt(.Machine$double.eps)/1000)
                ) {
    ##
    anal_vcov <- function(pars, df){
        neg_hess <- function(pars,df) { #computes the negative hessian
            for (nm in names(pars)) assign(nm,pars[[nm]])
            mu <- m*df$E + pi0*df$g + pi1*df$g*df$E
            sig <- sqrt((lambda0 + df$E*lambda1)^2)
            ##
            H <- matrix(data=NA, nrow=5, ncol=5)
            # second partials wrt beta_0
            H[1,1] <- sum(df$E^2*sig^(-2))
            H[1,2] <- H[2,1] <- sum(df$E*df$g*sig^(-2))
            H[1,3] <- H[3,1] <- sum(df$E^2*df$g*sig^(-2))
            H[1,4] <- H[4,1] <- 2*sum(df$E*(df$y - mu)*sig^(-3))
            H[1,5] <- H[5,1] <- 2*sum(df$E^2*(df$y - mu)*sig^(-3))
            # second partials wrt pi_0
            H[2,2] <- sum(df$g^2*sig^(-2))
            H[2,3] <- H[3,2] <- sum(df$E*df$g^2*sig^(-2))
            H[2,4] <- H[4,2] <- 2*sum(df$g*(df$y - mu)*sig^(-3))
            H[2,5] <- H[5,2] <- 2*sum(df$E*df$g*(df$y - mu)*sig^(-3))
            # second partials wrt pi_1
            H[3,3] <- sum(df$E^2*df$g^2*sig^(-2))
            H[3,4] <- H[4,3] <- 2*sum(df$E*df$g*(df$y - mu)*sig^(-3))
            H[3,5] <- H[5,3] <- 2*sum(df$E^2*df$g*(df$y - mu)*sig^(-3))
            # second partials wrt lambda_0
            H[4,4] <- sum(3*(df$y - mu)^2*sig^(-4) - sig^(-2))
            H[4,5] <- H[5,4] <- sum(3*df$E*(df$y - mu)^2*sig^(-4) - df$E*sig^(-2))
            # second partials wrt lambda_1
            H[5,5] <- sum(3*df$E^2*(df$y - mu)^2*sig^(-4) - df$E^2*sig^(-2))
            return(H)
        }
        H <- neg_hess(pars, df)
        vcov <- solve(H)
        rownames(vcov) <- colnames(vcov) <- names(pars)
        return(vcov)
    }
    ##
    ll<-function(pars,df,xx) { #this is the likelihood function to be maximized 
        for (nm in names(pars)) assign(nm,pars[[nm]])
        beta<-pars[-(1:5)]
        mu<-m*df$E+
            pi0*df$g+
            pi1*df$g*df$E+
            as.numeric(as.matrix(xx) %*% matrix(beta,ncol=1))
        sig<-sqrt((lambda0+df$E*lambda1)^2)
        ##
        d<-dnorm(df$y,mean=mu,sd=sig,log=TRUE) #sqrt(lambda0^2+df$E^2*lambda1^2+2*lambda0*lambda1*df$E))
        tr<- -1*sum(d)
        #print(pars)
        #print(tr)
        tr
    }
    gr<-function(pars,df,xx) { #this is the gradient
        for (nm in names(pars)) assign(nm,pars[[nm]])
        beta<-pars[-(1:5)]
        mu<-m*df$E+
            pi0*df$g+
            pi1*df$g*df$E+
            as.numeric(as.matrix(xx) %*% matrix(beta,ncol=1))
        sig<-sqrt((lambda0+df$E*lambda1)^2)
        ##
        dfdmu<-(df$y-mu)*sig^(-2)
        dfdm<-sum(dfdmu*df$E)
        dfdpi0<-sum(dfdmu*df$g)
        dfdpi1<-sum(dfdmu*df$g*df$E)
        dfdsigma <- (df$y - mu)^2*sig^(-3) - sig^(-1)
        dfdlambda0 <- sum(dfdsigma)
        dfdlambda1 <- sum(dfdsigma*df$E)
        dfdbeta<-rep(0,ncol(xx))
        for (i in 1:ncol(xx)) dfdbeta[i]<-sum(dfdmu*xx[,i])
        ##
        -1*c(dfdm,dfdpi0,dfdpi1,dfdlambda0,dfdlambda1,dfdbeta)
    }
    ##working out covariates
    if (length(covars)>0) {
        xx<-df[,covars]
        xx<-cbind(rep(1,nrow(xx)),xx) #for intercept
    } else xx<-matrix(rep(1,nrow(df)),ncol=1)
    nb<-ncol(xx)
    ##
    pars<-c(m=rnorm(1,sd=.05),pi0=rnorm(1,sd=.05),pi1=rnorm(1,sd=.05),lambda0=1+runif(1),lambda1=runif(1,max=.2),beta=rep(0,nb))
    fit<-optim(pars,ll,
               gr=gr,
               df=df,
               xx=xx,
               control=ctl.list,
               hessian=hess
               )
    est<-c(fit$par)
    if (hess) { #see https://stats.stackexchange.com/questions/27033/in-r-given-an-output-from-optim-with-a-hessian-matrix-how-to-calculate-paramet
        vc<-solve(fit$hessian) #note, no -1 given that i am doing that directly in the above
    } else {
        test<-try(vc<-anal_vcov(est,df))
        test<-class(test)
        count<-1
        while (test=='try-error' & count<50) {
            pars<-list(m=rnorm(1,sd=.05),pi0=rnorm(1,sd=.05),pi1=rnorm(1,sd=.05),lambda0=1+runif(1),lambda1=runif(1,max=.2))
            fit<-optim(pars,ll,
                       gr=gr,
                       df=df,
                       control=ctl.list,
                       hessian=TRUE
                       )
            test<-try(vc<-solve(fit$hessian))
            test<-class(test)
            count<-count+1
        }
        if (test=='try-error') return(list(est=est))
    }
    return(list(est=est,var=vc))
}



