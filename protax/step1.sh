#############################################
# Preliminaries: prepare 3 files
#         1) taxonomy (file 'taxonomy')
#         2) reference sequences (file 'refs.fa')
#         3) mapping from reference sequences to taxonomy (file 'seqid2tax')
#############################################

# define MODEL directory freely and set PROTAX to directory which contains the scripts

PROTAX=`pwd`/protaxscripts
MODEL=`pwd`/model

mkdir $MODEL
cp taxonomy $MODEL
cp refs.fa $MODEL
cp seqid2tax $MODEL
cd $MODEL

#############################################
# PART I: Reference sequence database
#############################################

# create USEARCH index:
~/Work/prog/usearch10.0.240_i86linux32 -makeudb_usearch refs.fa -output refs.udb

#############################################
# PART II: Taxonomy (choose A or B depending on your case)
#############################################

# create separate taxonomy file for each level# 
# in this example there are 7 levels in the taxonomy, change it if needed

#######################################
# A) Special treatment for Fungi taxonomy: skip nonFungi node in level 1
#######################################

perl $PROTAX/cptaxonomy.pl taxonomy > taxonomy.fungi
perl $PROTAX/thintaxonomy.pl 1 taxonomy.fungi > tax1
# remove nonFungi from taxonomy after kingdom level (1) since it is not divided further (this is quick and dirty hack so that no changes need to be done in script  generate_training_data.pl)
for LEVEL in 2 3 4 5 6 7
do
 perl $PROTAX/thintaxonomy.pl $LEVEL taxonomy.fungi | grep -v -w nonFungi > tax$LEVEL
done

#######################################
# B) Normal case where no need to skip any nodes
#######################################

for LEVEL in 1 2 3 4 5 6 7
do
 perl $PROTAX/thintaxonomy.pl $LEVEL taxonomy.priors > tax$LEVEL
done

#############################################
# PART III: Generating training data lists
#############################################

# here 5000 training samples out of which 5% represent unknown taxa
# NOTE: include also LEVEL 1 if normal taxonomy (this is for Fungi example)

for LEVEL in 2 3 4 5 6 7
do
 echo "LEVEL $LEVEL"
 perl $PROTAX/seqid2taxlevel.pl $LEVEL seqid2tax > ref.tax$LEVEL
 perl $PROTAX/get_all_reference_sequences.pl $LEVEL tax$LEVEL ref.tax$LEVEL rseqs$LEVEL
 perl $PROTAX/generate_training_data.pl tax$LEVEL ref.tax$LEVEL rseqs$LEVEL 4750 2 no train.level$LEVEL 
 perl $PROTAX/generate_unk_training_data.pl $LEVEL tax$LEVEL ref.tax$LEVEL rseqs$LEVEL 250 2 no train.unk$LEVEL 
 cat train.level$LEVEL train.unk$LEVEL > train$LEVEL
 cut -f6 -d" " train$LEVEL | sort | uniq > train${LEVEL}.id
done

cat train[2,3,4,5,6,7].id | sort | uniq > train.ids
