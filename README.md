# Companion code for "AKT mutant allele-specific activation dictates pharmacologic sensitivities"

![alt text](https://github.com/agorelick/akt1/blob/master/data/wiki.png "Fig 3. Structural and signaling impact of AKT in-frame indels")

Heatmaps comparing the displacement of AKT1 residues in the (time-averaged) E17K and P68-77dup protein structures (Fig 3a) can be regenerated from precomputed residue-residue distance maps as follows:

Clone the github repo:
```shell
git clone https://github.com/agorelick/akt1
cd akt1
```

Install required R packages from CRAN:
```r
## check for missing required packages, install them.
required.packages <- c('data.table','ggplot2','cowplot','here')
new.packages <- required.packages[!(required.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)>0) install.packages(new.packages)
```

Run the fig_3a.R script:
```r
## create a PDF with figure 3a: fig_3a.pdf
Rscript fig_3a.R
```

