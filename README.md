# HaplogroupEstimation
Script used for the haplogroup estimation across smaller variant windows.

Usage: HaplogroupEstimation.py [-h] --sam [SAM file] --vcf [vcf file] --window [window size]

The SAM file can be read from the stdin by using `--sam -`. The default window size is 35.

The output format is tab seperated and contains the following:

* reference or contig
* Positions in the SNP window (dash seperated)
* haplotype 1
* count of haplotype 1
* haplotype 2
* count of haplotype 2
* etc.

For example:
```
CP	10571-10577	AT	11	TA	32	
CP	10719-10725	TT	14	AG	16	
```

