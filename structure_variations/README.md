### detection of structure variation
A series of custom scripts (https://github.com/comery/marmoset) identified and sorted our SNPs, indels in the alignments. We employed svmu (v0.4-alpha), Assemblytics (v1.2) and SyRi (v1.0), to detect SVs from Mummer alignment. After several test rounds, we found that svmu reported more accurate large indels, Assemblytics detected CNVs, particularly tandem repeats, whereas SyRi detected other SVs well. Hence, we employed these three methods and combined the results as confident SVs. For Mummer alignment, we used the settings of “nucmer --maxmatch -c 500 -b 500 -l 100” and following filtering of “delta-filter -m -i 90 -l 100”, which was recommended by SyRi (https://schneebergerlab.github.io/syri/), while we used default parameters for svmu and Assemblytics. 

#### run syri

```shell
sh syri.sh
```


### run svmu

```shell
sh svmu.sh
```

### run Assemblytics

```shell
sh Assemblytics.sh
```