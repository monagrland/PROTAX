#############################################
# NOTE: this example is for Fungi, if normal taxonomy, include also LEVEL 1
#############################################

#### check classification with training samples (here using those which were used for training level 7)
# cd into model directory, if protaxscripts not located as shown below, change accordingly

mkdir checktrain
cd checktrain

# MDIR is directory where the taxonomy and refseq related parts of model are, PDIR contains parameters (although here set to be the same directory)
PROTAX=../../protaxscripts
MDIR=../
PDIR=$MDIR
SIMFILE=$PDIR/train.m8

# this is for Fungi taxonomy, so classification starts from node 1 (Fungi)
perl $PROTAX/trainsample2init.pl 1 $MDIR/train7 > query1.logprob

# parent probs from previous level classification
for LEVEL in 2 3 4 5 6 7
do
 echo "LEVEL $LEVEL" 
 PREVLEVEL=$((LEVEL-1))
 IFILE=query${PREVLEVEL}.logprob
 OFILE=query${LEVEL}.logprob
 perl $PROTAX/trainclassify4Q.pl 0.05 $IFILE $MDIR/tax$LEVEL $MDIR/ref.tax$LEVEL $MDIR/rseqs$LEVEL $PDIR/mcmc${LEVEL} map $SIMFILE 0 .1 $OFILE 1
done

##### calculate correctness (note: correct label is the one which training sample mimicked, it can be e.g. unknown species)

for LEVEL in 2 3 4 5 6 7
do
 perl $PROTAX/trainsample2correct.pl $LEVEL $MDIR/taxonomy $MDIR/train7 > query${LEVEL}.tax
 perl $PROTAX/trainsample2addcor.pl query${LEVEL}.logprob query${LEVEL}.tax $MDIR/tax$LEVEL > query${LEVEL}.cor
done

# check nseq2 samples
paste query7.cor query7.tax | grep nseq2 | cut -f2,3 > foo7

##### bias accuracy plots
R

aplot = function(prob,correct,add=F,col="black",name="") {
  n=length(prob)
  if (!add) {
   plot(0,type="n",xlab="cumulative prob %",ylab="cumulative correct %",xlim=c(0,100),ylim=c(0,100),las=1)
   grid()
   abline(0,1,col="gray")
   title(sprintf("%s\n%d items, correct %.1f %%",name,n,sum(correct)/n*100))
  }
  s=sort(prob,dec=F,index.return=T)
  x=cumsum(prob[s$ix])/n*100
  y=cumsum(correct[s$ix])/n*100
  lines(x,y,col=col)
  points(x[n],y[n],pch=19,cex=1,col=col)
}

nimi=c("kingdom","phylum","class","order","family","genus","species")

pdf("train_biasaccuracy.pdf",width=10,height=6)

par(mfrow=c(2,4))
for (i in 2:7) {
 file=sprintf("query%d.cor",i)
 a=read.table(file,header=F)
 aplot(a[,2],a[,3],name=sprintf("ITS2: %s",nimi[i]))
}

# check nseq2 samples
#paste query7.cor query7.tax | grep nseq2 | cut -f2,3 > foo7
a=read.table("foo7",header=F)
aplot(a[,1],a[,2], name="species nseq2")

dev.off()
