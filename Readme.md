## Single sample GSEA
The R-script ```1_run_ssGSEA.r``` provides a wrapper to the ssGSEA R-program which performs the actual single sample Gene Set Enrichment analysis (GSEA).

### Instructions
To execute the script source it into a running R-session.

- RStudio: open the file and press 'Source' in the upper right part of the editor window 
- R-GUI: drag and drop this file into an R-GUI window

Input files and databases can be specified via Windows file dialogs that will be automatically invoked. The first dialog lets you choose a folder containing input files in [GTC v1.2](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GCT:_Gene_Cluster_Text_file_format_.28.2A.gct.29) format. The script will loop over all gct files in this directory and run ssGSEA on each file separately. The second dialog window lets the user choose a gene set database such as [MSigDB](http://software.broadinstitute.org/gsea/msigdb/). A current version of MSigDB databases can be found in the 'db' subfolder. 


Other paramaters for ssGSEA can be altered inside the parameters section in ```1_run_ssGSEA.r```. The default parameters have been choosen carfully and should provide reliable results for most use-case scenerios. 

### Gene Set Enrichment Analysis
For information about the method please visit http://software.broadinstitute.org/gsea/.

### Changes to the original ssGSEA
Original code written by Pablo Tamayo. Adapted with additional modifications by D. R. Mani and Karsten Krug. Adaptions include:

- support of multiple CPU cores (```doParallel``` R-package)
- improved handling of missing values
- scoring of directional gene sets
- basic error handling
- general performance improvements
- additional output files like rank plots and parameters


### Reference

1. Abazeed, M. E., Adams, D. J., Hurov, K. E., Tamayo, P., Creighton, C. J., Sonkin, D., et al. (2013).
       Integrative Radiogenomic Profiling of Squamous Cell Lung Cancer. Cancer Research, 73(20), 6289-6298.
       http://doi.org/10.1158/0008-5472.CAN-13-1616
       
       
2. Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., et al. (2005).
   Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles.
  Proceedings of the National Academy of Sciences of the United States of America, 102(43), 15545-15550.
