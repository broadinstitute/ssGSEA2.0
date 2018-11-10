## ssGSEA2.0 / PTM-SEA / PTMsigDB

Tools to perform gene-centric **single sample Gene Set Enrichment Analysis (ssGSEA)** and site-centric **PTM Signature Enrichment Analysis (PTM-SEA)**.


#### PTM signatures database (PTMsigDB)
PTMsigDB is a database comprised of modification site-specific signatures of perturbations, kinase activities and signaling pathways curated from more than 2,500 publications. 

[v1.8.1](https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/db/c1.all.v6.1.symbols.gmt)


#### PTM-Signature Enrichment analysis (PTM-SEA)
PMT-SEA is an adaption of ssGSEA to perform site-specific siganture analysis by scoring PTMsigDB's bi-derectional signature-sets.



https://github.com/broadinstitute/ssGSEA2.0/blob/master/db/c1.all.v6.1.symbols.gmt


### Single-sample GSEA 2.0



### Instructions to run ssGSEA2.0 / PTM-SEA

The tools can be executed in three different XX

#### Command line interface
```ssgsea-cli.R``` can be executed from the command line, please see ```ssgsea-cli.R --help``` for instructions. 

#### R


The R-scripts ```ssgsea-gui.R``` and ```ssgsea-cli.R``` provide wrapper to the ssGSEA/PTM-SEA R-program which performs the actual enrichment analysis. 



#### ```ssgsea-gui.R``` 
- is aimed towards users with little or no knowledge of Linux command line.


### Instructions to run ```ssgsea-gui.R```

##### **Windows OS**
To run the script source it into a running R-session.

- RStudio: open the file and press 'Source' in the upper right part of the editor window 
- R-GUI: drag and drop this file into an R-GUI window

Input files and databases can be specified via Windows file dialogs that will be automatically invoked. The first dialog lets you choose a folder containing input files in [GCT v1.2](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GCT:_Gene_Cluster_Text_file_format_.28.2A.gct.29) or [GCT v1.3](https://clue.io/connectopedia/gct_format) format. The script loops over all GCT files in this directory and runs ssGSEA on each file separately. The second dialog window lets the user choose one or multiple gene set databases in [GMT](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29) format such as [MSigDB](http://software.broadinstitute.org/gsea/msigdb/). A current version of MSigDB databases can be found in the ```db``` subfolder. 

##### **iOS/MAC** 
In order to invoke file dialogs as decribed above, the [XQuartz](https://www.xquartz.org) X Window System is required. Once installed ```ssgsea-gui.R``` can be sourced into an R session.   

##### **ssGSEA/PTM-SEA parameters**

Other paramaters for ssGSEA/PTM-SEA can be altered inside the parameters section in ```ssgsea-gui.R```. The default parameters have been choosen carefully and should provide reliable results for most use-case scenerios. 

### Gene Set Enrichment Analysis (GSEA)
For information about the method please visit http://software.broadinstitute.org/gsea/.

### Changes to the original ssGSEA implementation
Original code written by Pablo Tamayo. Adapted with additional modifications by D. R. Mani and Karsten Krug. Adaptions include:

- support of multiple CPU cores (```doParallel``` R-package)
- support of GCT v1.3 format using functions from [cmapR](https://github.com/cmap/cmapR)
- improved handling of missing values
- scoring of directional gene sets (PTMsigDB)
- basic error handling
- general performance improvements
- additional output files like rank plots and parameter files


### References

1. Abazeed, M. E., Adams, D. J., Hurov, K. E., Tamayo, P., Creighton, C. J., Sonkin, D., et al. (2013).
       Integrative Radiogenomic Profiling of Squamous Cell Lung Cancer. Cancer Research, 73(20), 6289-6298.
       http://doi.org/10.1158/0008-5472.CAN-13-1616

2.  Barbie, D. A., Tamayo, P., Boehm, J. S., Kim, S. Y., Susan, E., Dunn, I. F., . Hahn, W. C. (2010). Systematic RNA interference reveals that oncogenic KRAS- driven cancers require TBK1, 462(7269), 108-112. https://doi.org/10.1038/nature08460
       
3. Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., et al. (2005).
   Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles.
  Proceedings of the National Academy of Sciences of the United States of America, 102(43), 15545-15550. http://doi.org/10.1073/pnas.0506580102

***
