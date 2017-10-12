## Bias, robustness and scalability in differential expression analysis of single-cell RNA-seq data

This app allows interactive browsing and exploration of the results presented in the paper

* C Soneson & MD Robinson: [Bias, robustness and scalability in differential expression analysis of single-cell RNA-seq data](http://biorxiv.org/content/early/2017/05/28/143289). bioRxiv doi:10.1101/143289 (2017).

For more information about the study design and discussion of the results, please see the full paper. All the code used to perform the analyses can be found on GitHub: [https://github.com/csoneson/conquer_comparison](https://github.com/csoneson/conquer_comparison). 

### Getting started

To get started, select the methods you want to include in the visualization in the panel on the left side (this panel will be accessible from the individual visualization tabs as well). Then go to the **"Select input"** tab and select the data sets, sample sizes and filterings to include in the visualizations. After that, the individual evaluations can be browsed in the respective tabs. 

### Notes

Most of the tabs are responsive to the selections made in the "Select input" tab. However, some of the evaluations are only performed on a subset of the data set instances:

- *Fraction NA adj.p*: only the 17 real scRNA-seq data sets are considered. 
- *Number of DEGs*: only the 9 real "signal" scRNA-seq data sets are considered.
- *Type 1 error*: only the 8 real "null" scRNA-seq data sets are considered.
- *FDP at adj.p = 0.05 cutoff*, *TPR at adj.p = 0.05 cutoff*, *AUROC*: only the three simulated "signal" scRNA-seq data sets are considered.
- *Time requirement*: all 23 scRNA-seq data sets (real and simulated) are considered.
- *Concordance scores*: only the 8 real scRNA-seq data sets for which both "signal" and "null" instances are available are considered. 
- *DE genes characteristics*: only the 8 real "null" scRNA-seq data set instances with no expression filtering are considered. 
- *Cross-method similarity*: only the 17 real scRNA-seq data sets are considered.

In addition, the results in the "Performance summary" tab are not affected by any filtering in the "Select input" tab. The performance scores are pre-computed for each method as described in the original paper. However, the results can be filtered by excluding some evaluation criteria and/or some methods. 