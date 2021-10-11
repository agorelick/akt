# AKT mutant allele-specific activation dictates pharmacologic sensitivities

![alt text](https://github.com/agorelick/akt1/blob/master/data/wiki.png "Fig 3. Structural and signaling impact of AKT in-frame indels")

# Instructions

## Fig 3a, Supplementary Figure 3a-b:

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

Run the do.R script:
```r
## create PDFs with figures 3a, 3a-b:
Rscript do.R
```

## Fig 3b-d

AKT1 structural diagrams in Figures 3b-d were generated using UCSF Chimera. A Chimera scene (.py format) containing an overlay of AKT1 WT and P68-77dup structures (3o96 template) is contained in this repo (`data/chimera_wt_dup_overlay.py`), and can be opened with UCSF Chimera to regenerate Figures 3b-d. UCSF Chimera can be freely downloaded here: https://www.cgl.ucsf.edu/chimera/download.html.

With UCSF Chimera, open  the scene using: `File > Open > [path-to-repo]/data/chimera_wt_dup_overlay.py`. By default, only the wild-type AKT1 protein is shown (as in Fig 3b). To see the P68-77dup protein overlay, first select the P68-77dup protein: `Select > Chain > A > average_dupe.pdb`. Then you can show the ribbon for the P68-77dup structure: `Actions > Ribbon > Show`.



### Citation
- URL: PENDING
- DOI: PENDING

### Contact
E-mail any questions to [tbhattarai@loxooncology.com](mailto:tbhattarai@loxooncology.com?subject=[GitHub]%20AKT1%20paper).

