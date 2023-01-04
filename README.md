# ANEVA-h

A toolkit for quantifying genetic variation in gene dosage from allelic specific expression (ASE) data. ANEVA-h implements the ANEVA method for haplotype-based AE data generated from phASER. For more details about ANEVA, please see Mohammadi P et al (2019) doi:10.1126/science.aay0256

### Installation Instructions
***

* The following are dependencies of the package. To install, run the following R code in an interactive Rstudio session or terminal
``` R
install.packages(c("pracma", "progress"))
```

* Clone this [repository](https://github.com/ejsong37/ANEVA-h/) and navigate to the ANEVA-h folder.  

### Input Files
***

Both refcounts and altcounts are tab-separated (TSV) files with an index column and a column `name` containing the names of genes in addition to the columns for each sample. Ref counts will have the reference allele counts, and Alt counts will have the alternative allele counts from ASE Data.

The format is outlined below:

|   | name   | sample 1 | sample 2 | ... | sample n |
| - | ------ | -------- | -------- | --- | -------- |
| 0 | Gene 1 | Count 1  | Count 2  | ... | Count n  |
| 1 | Gene 2 | Count 3  | Count 4  | ... | Count m  |
| ... | ... | ... | ... | ... | ... |



