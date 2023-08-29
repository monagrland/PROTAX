#############################################
# NOTE: this example is for Fungi, if normal taxonomy, include also LEVEL 1
#############################################

# in this example there are 4 predictors, if creating new predicotors, do the changes in 3 Perl scripts:
# create_xdata4Q.pl, trainclassify4Q.pl, classify4Q.pl
# difference between trainclassify.pl and classify.pl is that trainclassify automatically removes the sequence being classifyid from references based on its seqid

# cd into model directory, if protaxscripts not located as shown below, change accordingly

PROTAX=../protaxscripts

# map training sequences against reference sequences

perl $PROTAX/fastagrep.pl train.ids refs.fa > train.fa
~/Work/prog/usearch10.0.240_i86linux32 -usearch_global train.fa -db refs.udb -id 0.75 -maxaccepts 1000 -strand both -userfields query+target+id -userout train.m8

for LEVEL in 2 3 4 5 6 7
do
 echo $LEVEL
 perl $PROTAX/create_xdata4Q.pl 0.05 train$LEVEL tax$LEVEL ref.tax$LEVEL rseqs$LEVEL train.m8 train${LEVEL}.xdat 1
done

### MCMC can be run separately (in parallel) for each taxonomy level
R
source("../protaxscripts/amcmc.rcode.txt")
library(compiler)
logprior=cmpfun(logprior)
loglikelihood=cmpfun(loglikelihood)
adaptiveMCMC=cmpfun(adaptiveMCMC)

# if changing number of predictors, remember to update num.params
num.params=1+4
ind=1001:2000

# here doing re-initialization of parameter covariance matrix after first run of MCMC for each chain
# change this and number of interations based on the look of traceplots
# from each MCMC, save MAP estimate

# pdf("mcmc.pdf",width=10,height=5)
dev.new(width=10,height=5)

dat=read.xdata("train2.xdat")
pp1=adaptiveMCMC(dat,num.params,10000,2000,1000,rseed=1,info=1)
traceplot.all(pp1,ind,num.levels=1, title="L2a")

k=which.max(pp1$postli[ind])
write.postparams(pp1,"mcmc2a",ind[k])

initstate=initialize.adaptation(pp1$params[2000,])
pp1=adaptiveMCMC(dat,num.params,10000,2000,1000,rseed=1,info=1,prev.state=initstate)
traceplot.all(pp1,ind,num.levels=1, title="L2b")

k=which.max(pp1$postli[ind])
write.postparams(pp1,"mcmc2",ind[k])

dat=read.xdata("train3.xdat")
pp1=adaptiveMCMC(dat,num.params,10000,2000,1000,rseed=1,info=1)
traceplot.all(pp1,ind,num.levels=1, title="L3a")

k=which.max(pp1$postli[ind])
write.postparams(pp1,"mcmc3a",ind[k])

initstate=initialize.adaptation(pp1$params[2000,])
pp1=adaptiveMCMC(dat,num.params,10000,2000,1000,rseed=1,info=1,prev.state=initstate)
traceplot.all(pp1,ind,num.levels=1, title="L3b")

k=which.max(pp1$postli[ind])
write.postparams(pp1,"mcmc3",ind[k])

dat=read.xdata("train4.xdat")
pp1=adaptiveMCMC(dat,num.params,10000,2000,1000,rseed=1,info=1)
traceplot.all(pp1,ind,num.levels=1, title="L4a")

k=which.max(pp1$postli[ind])
write.postparams(pp1,"mcmc4a",ind[k])

initstate=initialize.adaptation(pp1$params[2000,])
pp1=adaptiveMCMC(dat,num.params,10000,2000,1000,rseed=1,info=1,prev.state=initstate)
traceplot.all(pp1,ind,num.levels=1, title="L4b")

k=which.max(pp1$postli[ind])
write.postparams(pp1,"mcmc4",ind[k])

dat=read.xdata("train5.xdat")
pp1=adaptiveMCMC(dat,num.params,10000,2000,1000,rseed=1,info=1)
traceplot.all(pp1,ind,num.levels=1, title="L5a")

k=which.max(pp1$postli[ind])
write.postparams(pp1,"mcmc5a",ind[k])

initstate=initialize.adaptation(pp1$params[2000,])
pp1=adaptiveMCMC(dat,num.params,10000,2000,1000,rseed=1,info=1,prev.state=initstate)
traceplot.all(pp1,ind,num.levels=1, title="L5b")

k=which.max(pp1$postli[ind])
write.postparams(pp1,"mcmc5",ind[k])

dat=read.xdata("train6.xdat")
pp1=adaptiveMCMC(dat,num.params,10000,2000,1000,rseed=1,info=1)
traceplot.all(pp1,ind,num.levels=1, title="L6a")

k=which.max(pp1$postli[ind])
write.postparams(pp1,"mcmc6a",ind[k])

initstate=initialize.adaptation(pp1$params[2000,])
pp1=adaptiveMCMC(dat,num.params,10000,2000,1000,rseed=1,info=1,prev.state=initstate)
traceplot.all(pp1,ind,num.levels=1, title="L6b")

k=which.max(pp1$postli[ind])
write.postparams(pp1,"mcmc6",ind[k])

dat=read.xdata("train7.xdat")
pp1=adaptiveMCMC(dat,num.params,10000,2000,1000,rseed=1,info=1)
traceplot.all(pp1,ind,num.levels=1, title="L7a")

k=which.max(pp1$postli[ind])
write.postparams(pp1,"mcmc7a",ind[k])

initstate=initialize.adaptation(pp1$params[2000,])
pp1=adaptiveMCMC(dat,num.params,10000,2000,1000,rseed=1,info=1,prev.state=initstate)
traceplot.all(pp1,ind,num.levels=1, title="L7b")

k=which.max(pp1$postli[ind])
write.postparams(pp1,"mcmc7",ind[k])

#dev.off()

