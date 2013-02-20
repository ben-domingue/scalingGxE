a package to compute large numbers of fixed effects. 

some example code. suppose you have a data frame "dat" with columns "outcome", "prior", and "teacher" a teacher id. you can run a model where teacher fixed effects for predicted "outcome" scores conditioned on "prior" scores as follows:
 
###begin
library(manyFE)
fm<-formula('outcome~prior+fe1')
#traditional VAM
1->dat$fe1
manyFE(dat,f=fm,id.name="teacher",model.keep=TRUE)->VA.full
eb(VA.full)->VA
###end

you have to include the "fe1" column in the data frame. the "manyFE" function computes fixed effects and then "eb" shrinks them.

