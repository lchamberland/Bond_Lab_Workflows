
You need to edit your partition file (.charsets file) to get it into RAXML format format. 

```
vim yourpartitionfile.charsets
```

once you're in editing mode in the file
1. delete the first two lines (blank line and begin characters) and last two lines (charsets and END) using the command dd to delete a single line 

```
dd
```

2. Find and replace using the three commands below one at a time
delete all apostrophes
```
:%s/'//g
```
delete all semicolins
```
:%s/;//g
```
find and replace charset with DNA,
```
:%s/charset/DNA,/g
```

3. write (save) file and quit
```
:wq
```


Once you've edited your partition file, I usually rename my partition file to a shorter name. Make a slurm script using the command below. Make sure you add your personalized BATCH command block at the top. 
```
[your batch command block]

module load iqtree

iqtree2 -s YOUR-ALIGNMENT.phylip -p YOUR-EDITED-PARTITION-FILE.charsets --prefix your-taxa-IQTREE -m MFP -B 1000 -T AUTO
```
