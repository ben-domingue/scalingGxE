##this does the inference on xi based on output of ratiotest above
xi.test<-function(est,p.value=TRUE) {
    testfun<-"b[2]*b[5]-b[3]*b[4]=0"
    library(nlWaldTest)
    if (p.value) {
        nlWaldtest(coeff=est$est,Vcov=est$var,texts=testfun) 
    } else {
        nlConfint(coeff=est$est,Vcov=est$var,texts=sub("=0","",testfun,fixed=TRUE))
    }
}

