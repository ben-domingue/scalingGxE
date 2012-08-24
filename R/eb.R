
#this produces empirical bayes estimates with an object returned by many.factors()
eb<-function(obj) {
  obj[[1]]->mod
  obj[[2]]->fe
  obj[[3]]->id
  mod$model->dat
  mod$residuals->dat$res
  id->dat$id
  aggregate(dat$res,list(dat$id),sd)->out
  names(out)<-c("id","sd")
  data.frame(table(dat$id))->tab
  names(tab)<-c("id","n")
  merge(out,tab,by="id")->out
  out$se<-out$sd^2/out$n
  merge(fe,out,by=1)->out
  var(out$fe1)-mean(out$se,na.rm=TRUE)->v2
  out$shrink<-v2/(v2+out$se)
  out$fe1*out$shrink->out$eb
  out
}

