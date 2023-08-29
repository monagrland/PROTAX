PROTAX=protaxscripts
MDIR=model

ODIR=test
mkdir $ODIR

SIMFILE1=test.m8

~/Work/prog/usearch10.0.240_i86linux32 -usearch_global test.fa -db $MDIR/refs.udb -id 0.75 -maxaccepts 1000 -strand both -userfields query+target+id -userout $SIMFILE1

grep '^>' test.fa | cut -c2- > $ODIR/query.ids

perl $PROTAX/testsample2init.pl 1 $ODIR/query.ids > $ODIR/query1.logprob

# parent probs from previous level classification
for LEVEL in 2 3 4 5 6 7
do
 echo "LEVEL $LEVEL" 
 PREVLEVEL=$((LEVEL-1))
 IFILE=$ODIR/query${PREVLEVEL}.logprob
 OFILE=$ODIR/query${LEVEL}.logprob
 perl $PROTAX/classify4Q.pl 0.05 $IFILE $MDIR/tax$LEVEL $MDIR/ref.tax$LEVEL $MDIR/rseqs$LEVEL $MDIR/mcmc${LEVEL} map $SIMFILE1 0 .1 $OFILE 1
done

###################
# .logprob files contain list of node ids and logprobs, convert node ids to taxonomic names and logprobs to probs
# NOTE: for final output, we can add information from best matching reference sequences etc., this is just an example

for LEVEL in 2 3 4 5 6 7
do
 perl $PROTAX/nameprob.pl $MDIR/taxonomy $ODIR/query${LEVEL}.logprob > $ODIR/query${LEVEL}.nameprob
done
