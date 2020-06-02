sim.data <- function(E,i=1,b0=.8,b1=.2,b2=0,b3=.05,h=sqrt(.6),a=.5,sigma=1,scaling=TRUE) {
    N <- length(E)
    G <- rnorm(N,0,1)
    eps <- rnorm(N,0,sigma)
    if (scaling){
        e <- sqrt(1-h^2)
        ystar <- h*G+e*eps
        y <- a*E+(b0 + b1*E)*ystar
    } else {
        y <- b1*G+b2*E+b3*G*E+eps
    }
    df <- data.frame(E=E,y=y,g=G)
    df
}
