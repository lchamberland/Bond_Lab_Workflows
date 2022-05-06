# CCGP_aptostichus
10x genomics assembly


```
ragtag.py scaffold -f 1200 -t 10 ref/YBT-1520.fasta output_dir/ZZQ-130.fasta.
```
DOI : https://doi.org/10.1128/mra.00887-21


cats :https://doi.org/10.1093/jhered/esaa057

Universal Adapters - these do not change 

```
[adapters]
i7:GATCGGAAGAGCACACGTCTGAACTCCAGTCAC*ATCTCGTATGCCGTCTTCTGCTTG
i5:AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT*GTGTAGATCTCGGTGGTCGCCGTATCATT
```
Tag sequence are your unique barcodes (i7 and i5) for each well
```
[tag sequences]
```
tag map
```

```

```
illumiprocessor \
    --input raw-reads/ \
    --output clean-reads \
    --config illumiprocessor.conf \
    --cores 20
```
