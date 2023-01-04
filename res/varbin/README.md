# File description

## Gene level index at 5k, 20k, and 50k bin resolution
_*grch37.5k.gene.index.txt.gz, grch37.20k.gene.index.txt.gz, grch37.50k.gene.index.txt.gz:*_
Table containg gene coordinates, and bin.id, which maps back to chromInfo.  The gene table 
itself is constant across all three, however, the bin.id mapping changes.  Positions were
mapped with _`mapBins2Genes.R`_.

e.g.

```bash
zcat grch37.20k.gene.index.txt.gz | head | column -t
```
| seqnames |  start |    end | width | strand | ensembl_gene_id | hgnc.symbol | gene.type      | bin.id | arm | band    | stain |
|----------+--------+--------+-------+--------+-----------------+-------------+----------------+--------+-----+---------+-------|
|        1 |  69091 |  70008 |   918 | *      | ENSG00000186092 | OR4F5       | protein_coding |      1 | 1p  | 1p36.33 | gneg  |
|        1 | 134901 | 139379 |  4479 | *      | ENSG00000237683 |             | protein_coding |      1 | 1p  | 1p36.33 | gneg  |
|        1 | 367640 | 368634 |   995 | *      | ENSG00000235249 | OR4F29      | protein_coding |      1 | 1p  | 1p36.33 | gneg  |
|        1 | 621059 | 622053 |   995 | *      | ENSG00000185097 | OR4F16      | protein_coding |      1 | 1p  | 1p36.33 | gneg  |
|        1 | 738532 | 739137 |   606 | *      | ENSG00000269831 |             | protein_coding |      1 | 1p  | 1p36.33 | gneg  |
|        1 | 818043 | 819983 |  1941 | *      | ENSG00000269308 |             | protein_coding |      1 | 1p  | 1p36.33 | gneg  |
|        1 | 860260 | 879955 | 19696 | *      | ENSG00000187634 | SAMD11      | protein_coding |      1 | 1p  | 1p36.33 | gneg  |
|        1 | 861264 | 866445 |  5182 | *      | ENSG00000268179 |             | protein_coding |      1 | 1p  | 1p36.33 | gneg  |
|        1 | 879584 | 894689 | 15106 | *      | ENSG00000188976 | NOC2L       | protein_coding |      1 | 1p  | 1p36.33 | gneg  |


## Varbin.py internal files

### varbin bin boundaries
Tables used by `varbin.py` to estimate read count matrices, and perform integer copy number 
calculations.  File contains no header information. Columns represent chromosome, bin start, 
bin start absoulte position (cummulative sum of start), bin end, bin end absolute position, 
and the number of mappable positions in the bin

_*grch37.bin.boundaries.5k.bowtie.k50.sorted.txt.gz, 
  grch37.bin.boundaries.20k.bowtie.k50.sorted.txt.gz,
  grch37.bin.boundaries.50k.bowtie.k50.sorted.txt.gz:*_

e.g.

```bash
zcat varbin.gc.content.20k.bowtie.k50.grch37.txt | head 
```
| 1 |       0 |       0 |  892036 | 892036 | 124672 |
| 1 |  892036 |  892036 | 1033528 | 141492 | 124672 |
| 1 | 1033528 | 1033528 | 1172212 | 138684 | 124672 |
| 1 | 1172212 | 1172212 | 1311902 | 139690 | 124672 |
| 1 | 1311902 | 1311902 | 1487403 | 175501 | 124672 |
| 1 | 1487403 | 1487403 | 1699006 | 211603 | 124672 |
| 1 | 1699006 | 1699006 | 1847901 | 148895 | 124672 |
| 1 | 1847901 | 1847901 | 1989199 | 141298 | 124672 |
| 1 | 1989199 | 1989199 | 2127306 | 138107 | 124672 |
| 1 | 2127306 | 2127306 | 2264445 | 137139 | 124673 |


### grch37 chromosome sizes
size in basepairs of each chromosme and cumulative start position
_*grch37.chrom.sizes.txt.gz:*_


### GC content
Extended bin boundary file containing the gc content of each bin, number of genes, etc. 
Files contain headers.

_*varbin.gc.content.20k.bowtie.k50.grch37.txt.gz:*_
_*varbin.gc.content.50k.bowtie.k50.grch37.txt.gz:*_
_*varbin.gc.content.5k.bowtie.k50.grch37.txt.gz:*_


e.g.

```bash
zcat varbin.gc.content.20k.bowtie.k50.grch37.txt.gz | head
```
| bin.chrom | bin.start | bin.end | bin.length | gene.count | cgi.count | dist.telomere |  gc.content | arm | band    | stain |
|-----------+-----------+---------+------------+------------+-----------+---------------+-------------+-----+---------+-------|
|         1 |         0 |  892035 |     892036 |         17 |        19 |             0 | 0.445098328 | 1p  | 1p36.33 | gneg  |
|         1 |    892036 | 1033527 |     141492 |          7 |        27 |        892036 | 0.625886976 | 1p  | 1p36.33 | gneg  |
|         1 |   1033528 | 1172211 |     138684 |          8 |        20 |       1033528 |  0.60129503 | 1p  | 1p36.33 | gneg  |
|         1 |   1172212 | 1311901 |     139690 |         11 |        21 |       1172212 | 0.624797766 | 1p  | 1p36.33 | gneg  |
|         1 |   1311902 | 1487402 |     175501 |         11 |        25 |       1311902 | 0.573102148 | 1p  | 1p36.33 | gneg  |
|         1 |   1487403 | 1699005 |     211603 |          8 |        25 |       1487403 | 0.553673625 | 1p  | 1p36.33 | gneg  |
|         1 |   1699006 | 1847900 |     148895 |          2 |         4 |       1699006 | 0.491836529 | 1p  | 1p36.33 | gneg  |
|         1 |   1847901 | 1989198 |     141298 |          4 |        17 |       1847901 | 0.570015145 | 1p  | 1p36.33 | gneg  |
|         1 |   1989199 | 2127305 |     138107 |          2 |        13 |       1989199 | 0.565662856 | 1p  | 1p36.33 | gneg  |
