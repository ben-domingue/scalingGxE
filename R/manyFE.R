#this performs the stuff explained in section 2.1 of guimaraes and portugal
manyFE<-function(x,f,fe.name="fe1",id.name="tch_id1",tol=NULL,verb=TRUE,model.keep=FALSE) { 
  #making sure that the fe column is present and in the formula.
  if (!(fe.name %in% all.vars(f))) update.formula(f,paste(".~.+",fe.name))->f
  if (!(fe.name %in% names(x))) 1->x[[fe.name]]
  #
  c(id.name,all.vars(f))->nam
  x[,names(x) %in% nam]->x2
  x2[rowSums(is.na(x2))==0,]->x2
  rss1<-0
  dif<-100000
  counter<-1
  if (is.null(tol)) {
    test<-TRUE #convergence based upon difference of fe coefficient and unity
    tol<-.00005
  } else test<-FALSE #convergence based on change in RSS
  while (abs(dif)>tol) {
    grep(id.name,names(x2))->id.index
    grep(fe.name,names(x2))->fe.index
    grep(paste("^",all.vars(f)[1],"$",sep=""),names(x2))->outcome.index
    lm(f,data=x2)->m
    rss2<-rss1
    sum(resid(m)^2)->rss1
    if (!test) dif<-rss2-rss1 else {
      grep(fe.name,names(m$coefficients))->tmp.index
      if (!is.na(m$coef[tmp.index])) dif<-abs(m$coef[tmp.index]-1)
    }
    predict(m)->yhat
    coef(m)->foo
    grep(fe.name,names(foo))->index
    coef(m)[index]->hold
    if (is.na(hold)) hold<-0
    x2[,fe.index]->foo
    x2[,fe.index]<-x2[,outcome.index]-yhat+hold*foo
    unique(x2[,id.index])->id
    aggregate(x2[fe.index],list(x2[,id.index]),mean)->agg
    names(agg)<-c(id.name,fe.name)
    NULL->x2[,fe.index]
    merge(x2,agg,by=id.name)->x2
    if (verb) print(c(counter, dif)); counter<-counter+1
  }
  aggregate(x2[,fe.index],list(x2[,id.index]),mean)->fe
  c(id.name,fe.name)->names(fe)
  lm(f,data=x2,model=model.keep)->m
  x2[,id.index]->id
  list(m,fe,id)->to.ret
  names(to.ret)<-c("est","fe","id")
  to.ret
}

