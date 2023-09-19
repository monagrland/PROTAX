read.xdata = function (file) {
 # file format, each line contains 3 columns: weight node.priprob xdat1;xdat2;...;xdatL
 # each xdati contains: level,node.category,m,n,m-by-1 weightvector,elements of m-by-n matrix row by row separated by ','
 # white space must be only used for separating 3 columns in each line, otherwise ; and ,

 a=read.table(file, colClasses=c("numeric","numeric","character"))
 dat=rep(list(0),nrow(a))

 for (i in 1:nrow(a)) {
  xdats=strsplit(a[i,3],';')[[1]]  
  foo.levels = rep(0, length(xdats)) 
  foo.ncs = rep(0, length(xdats)) 
  foo.xs = rep(list(0), length(xdats)) 
  foo.xw = rep(list(0), length(xdats)) 
  for (j in 1:length(xdats)) {
   xdat=as.numeric(strsplit(xdats[j],',')[[1]])
   nr=xdat[3]
   nc=xdat[4]
   if (length(xdat) != (4+nr+nr*nc)) {
    cat(sprintf("ERROR (read.xdata): file '%s' line %d/%d, xdat %d/%d matrix dimensions (%d-by-%d) don't agree with data (%d elements). (Note: file should include %d-by-1 weight vector)\n",file,i,nrow(a),j,length(xdats),nr,nc,length(xdat)-nr-4,nr))
    return(NA)
   }
   foo.levels[j]=xdat[1]
   foo.ncs[j]=xdat[2]
   foo.xw[[j]] = xdat[5:(4+nr)]
   foo.xs[[j]] = matrix(xdat[(5+nr):length(xdat)],nrow=nr,ncol=nc,byrow=T)
  }
  dat[[i]] = list(weight=a[i,1],node.priprob=a[i,2],level=foo.levels, node.category=foo.ncs, Xweight=foo.xw,X=foo.xs)
 }
 dat
}

merge.xdata = function (dat1,col1,dat2,col2) {
 # col1,col2 are column indices of dat1,dat2 to be merged to dat

 dat=dat1
 for (i in 1:length(dat)) {
  for (j in 1:length(dat[[i]]$level)) {
   dat[[i]]$X[[j]] = cbind(dat1[[i]]$X[[j]][,col1,drop=F], dat2[[i]]$X[[j]][,col2,drop=F])
  }
 }
 dat
}

write.postparams = function(pp, outfile, ind) {
 if (length(ind)==1) {
  # write row vector
  write.table(t(c(pp$postli[ind],pp$params[ind,])),outfile,col.names=F,row.names=F)
 }
 else {
  write.table(cbind(pp$postli[ind],pp$params[ind,]),outfile,col.names=F,row.names=F)
 }
}

initialize.adaptation = function(startparams) {
 # initval is initial values of parameters
 num.params = length(startparams)
 k=rep(1,num.params)
 ac=rep(0,num.params) 
 cc=rep(0,num.params)
 eig.sqrtval = matrix(1,nrow=num.params,ncol=1)
 eig.vec = diag(1,num.params,num.params)

 prevstate=list(k=k,ac=ac,cc=cc,eig.sqrtval=eig.sqrtval,eig.vec=eig.vec,params=matrix(startparams,nrow=num.params,ncol=1))
}

loglikelihood = function (xlist, params) {
 # xlist: weight node.priprob, level-identifiers, node.category-identifiers, X-matrix row weight vector list, X-matrix list
 # params: q level1params level2params ... levelLparams
 # d: parameter vector dimension in each level

 d=ncol(xlist$X[[1]])
 li=1

 for (i in 1:length(xlist$level)) {
  j=2+(xlist$level[[i]]-1)*d
  z = xlist$X[[i]] %*% params[j:(j+d-1),1,drop=F]
  # scale z so that maximum equals to zero
  ez = xlist$Xweight[[i]] * exp(z - max(z)) + 1e-10
  li = li * (sum(ez[xlist$node.category[i]])/sum(ez))
 }
 # p(mislabeling): inverse logit of param[1]
 eq=exp(min(c(params[1],600)))
 prob.q = eq/(1+eq)

 xlist$weight * log( prob.q * xlist$node.priprob + (1-prob.q) * li )
}

logprior = function(params, s2) {
 # N(mean=0,var=s2), ignore normalizing constant
 sum(-params*params/s2)
}

adaptiveMCMC = function (dat, num.params, s2, num.iterations, num.burnin, rseed=1, info=0, prev.state=NA) {
 
 # seed for random number generator

 set.seed(rseed)

 # dat is list, each element contains one or more X matrices with node information

 # save all samples (also from burnin iterations) for diagnostic purposes

 sample.params = matrix(0, nrow=num.iterations, ncol=num.params)
 sample.postli = rep(0, num.iterations)

 S1=matrix(0,nrow=num.params,ncol=1)
 S2=matrix(0,nrow=num.params,ncol=num.params)
 ns1s2=0

 if (length(prev.state)==1) {
  # start from scratch
  k=rep(1,num.params)
  ac=rep(0,num.params) 
  cc=rep(0,num.params)
  eig.sqrtval = matrix(1,nrow=num.params,ncol=1)
  eig.vec = diag(1,num.params,num.params)
  params = matrix(0,nrow=num.params, ncol=1)
 }
 else {
  # continue MCMC from given state
  k = prev.state$k
  ac = prev.state$ac
  cc = prev.state$cc
  eig.sqrtval = prev.state$eig.sqrtval
  eig.vec = prev.state$eig.vec
  params = prev.state$params
 }

 old.logli=logprior(params, s2)
 for (di in 1:length(dat)) {
  old.logli = old.logli + loglikelihood(dat[[di]], params)
 }

 for (iter in 1:num.iterations) {

  if (info) {
   cat(".")
   if (iter %% 50 == 0) {cat(sprintf(" %d\n",iter))}
  }

  # rotate SCAM 
  for (i in 1:num.params) {

  if (info>1) {cat(sprintf("  i %d\n",i))}

   # propose new parameter values

   new.params = params + rnorm(1)*k[i]*eig.sqrtval[i]*eig.vec[,i,drop=F]       

   old.logprior = logprior(params, s2)
   new.logprior = logprior(new.params, s2)
   new.logli=0
   for (di in 1:length(dat)) {
    new.logli = new.logli + loglikelihood(dat[[di]], new.params)
   }

   # Metropolis acceptance (priors already included in new and old logli)

   cc[i] = cc[i]+1

   if (info>1) {
    print(new.params)
    cat(sprintf("iter %d new.logli %g old.logli %g\n",iter, new.logli, old.logli))
   }

   logratio = new.logli + new.logprior - old.logli - old.logprior
   if (is.finite(logratio) & (runif(1) < exp(logratio))) {
    ac[i] = ac[i]+1
    params = new.params
    old.logli = new.logli
   }
  } # end of rotate SCAM: all eigenvectors utilized
  
  # save sample values from all iterations
  
  sample.params[iter,] = params
  sample.postli[iter] = old.logli

  if (iter <= num.burnin) {
   q = 1 + exp(-iter/500)
   w = 1 - 0.1 * exp(-iter/500)
   S1 = S1 + w * params
#   S2 = S2 + w * params %*% t(params)
   S2 = S2 + w * tcrossprod(params)
   ns1s2 = ns1s2 + w
   acr = ac/cc
   k = k * q ^ (acr-0.44)
   k[k>1e5]=1e5
   k[k<1e-5]=1e-5
   ac = w * ac
   cc = w * cc

   if (iter > 50) {
#    v = eigen( (S2 - S1 %*% t(S1) /ns1s2)/(ns1s2-1) + 1e-5*diag(num.params), symmetric=T )
    v = eigen( (S2 - tcrossprod(S1) /ns1s2)/(ns1s2-1) + 1e-5*diag(num.params), symmetric=T )
    eig.sqrtval = sqrt(v$values)
    eig.vec = v$vectors
   }
  }
 }

 list(params=sample.params, postli=sample.postli, ac=ac,cc=cc,k=k, eig.sqrtval=eig.sqrtval, eig.vec=eig.vec)
}

traceplot.one = function(x,y,mi=c(),num.histbars=20,name="") {
 yrange=range(y)
 # bottom, left, top, right
 par(mar=c(3,2,2,1))
 plot(x,y,pch=".", main=name, xlab="",ylab="", ylim=yrange)
 points(x[mi],y[mi],col="red",lwd=3,cex=2)

 # bottom, left, top, right
 par(mar=c(3,0,2,1))

 h=hist(y, breaks=seq(yrange[1],yrange[2],length.out=num.histbars), plot=F)
 barplot(h$density, axes=F, horiz=T, space=0, col="white")

 par(mar=c(3,2,2,1))
}

traceplot.all = function(pp,ind=0,num.levels=1,num.histbars=20,title="") {
 # utilizes pp$postli and pp$params
 betadim = (ncol(pp$params)-1)/num.levels

 num.xpanels = max(betadim, 4)
# par(mfrow=c(num.levels+1, 2*num.xpanels))
 # bottom, left, top, right
 par(mar=c(3,2,2,1))

 n=(num.levels+1)*2*num.xpanels
 layout(matrix(1:n,num.levels+1,2*num.xpanels,byrow=T), widths=rep(c(3,1),num.xpanels))

 # title panel
 plot(0,pch="",bty="n",xlab="",ylab="",xaxt="n",yaxt="n")
 text(0,title,cex=2)
 plot.new()

 if (length(ind)==1) {ind=1:length(pp$postli)}
 mi = which.max(pp$postli[ind])

 traceplot.one(ind,pp$postli[ind],mi,num.histbars,"log posterior")
 traceplot.one(ind,pp$params[ind,1],mi,num.histbars,"logit(mislabeling prob)")
 hist(exp(pp$params[ind,1])/(1+exp(pp$params[ind,1])), xlim=c(0,1), freq=F, main="mislabeling prob")
 plot.new()

 if (num.xpanels > 4) {
  for (j in 5:num.xpanels) {plot.new(); plot.new()}
 }

 i=1
 for (level in 1:num.levels) {
  for (beta in 1:betadim) {
   i=i+1
   name=sprintf("LEVEL%d beta%d",level,beta) 
   traceplot.one(ind,pp$params[ind,i],mi,num.histbars,name)
  }
  if (num.xpanels > betadim) {  
   for (j in (betadim+1):num.xpanels) {plot.new(); plot.new()}
  }
 }
}

amcmc.diagnostic.plot = function(pp) {
 par(mfrow=c(2,2))
 # c(bottom, left, top, right), default 5,4,4,2
 par(mar=c(3,3,2,2))
 d=length(pp$ac)
 barplot(pp$ac/pp$cc, names=1:d, main="acceptance ratios",las=1)
 barplot(pp$k, names=1:d, main="k-step sizes",las=1)
 barplot(pp$eig.sqrtval, names=1:d, main="sqrt(eigenvalues)",las=1)
 co=gray(seq(0,1,by=1/12))
 image(1:d,1:d,pp$eig.vec, zlim=c(-1,1),col=co, main="eigenvectors in rows\nwhite:+1 black:-1",las=1,xlab="",ylab="")
}

accuracy.plot = function(prob,correct,name="") {
  n=length(prob)
  plot(0,type="n",xlab="cumulative prob %",ylab="cumulative correct %",xlim=c(0,100),ylim=c(0,100),las=1)
  grid()
  abline(0,1,col="gray")
  s=sort(prob,dec=F,index.return=T)
  x=cumsum(prob[s$ix])/n*100
  y=cumsum(correct[s$ix])/n*100
  lines(x,y,col="black")
  points(x[n],y[n],pch=19,cex=1,col="black")
  title(sprintf("%s\n%d items, correct %.1f %%",name,n,sum(correct)/n*100))
}
