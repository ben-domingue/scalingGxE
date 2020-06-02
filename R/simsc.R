sim.sc<-function(i=1,E,b0=.8,b1=.2,h=sqrt(.6),a=.5) {
    E<-E()
    N<-length(E)
    e<-sqrt(1-h^2)
    G<-rnorm(N,0,1)
    eps<-rnorm(N,0,1)
    ystar<-h*G+e*eps
    y<-a*E+(b0+b1*E)*ystar 
    df<-data.frame(E=E,y=y,g=G)#,e=eps)
    df
}
