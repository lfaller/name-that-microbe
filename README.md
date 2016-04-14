# name-that-microbe

The idea is to make a version of probeseq that works with not only oral bugs but any kind of reference database to quickly identify 16S sequences from a miseq experiment.

# V regions based on E. coli annotation

| Region | Start | End |
|--------|-------|-----|
| V1 | 69 | 99 |
| V2 | 137 | 242 |
| V3 | 338 | 533 |
| V4 | 576 | 682 |
| V5 | 822 | 879 |
| V6 | 967 | 1046 |
| V7 | 1117 | 1173 |
| V8 | 1243 | 1294 |
| V9 | 1435 | 1465 |

# References used

* gg_otus-13_8-release/rep_set/97_otus.fasta
* NCBI's 16S database (to be added)

# Run a small test dataset

## Generate a reference database

```{bash}
cd ~/repos/name-that-microbe

IN=tests/small.fasta
OUT=tests/small.out
LOG=tests/small.log
AMBIG=tests/small_ambig.fasta

#V1-3 region
#START=69
#END=533

#V3-4 region
START=338
END=682

python utils/generate_ref_set.py -i $IN -o $OUT -s $START -e $END -l $LOG -a $AMBIG
```

## Annotate some fake sequence data

```{bash}
cd ~/repos/name-that-microbe

REF_FILE=tests/small.out
FASTA=tests/small_dataset.fasta
OUTFILE=tests/small_dataset_annotated.txt
LOG=tests/small_dataset.log

python annotate.py -r $REF_FILE -i $FASTA -o $OUTFILE -l $LOG
```
