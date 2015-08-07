# Consistency in large pharmacogenomic studies

This repository contains code to reproduce the results that we reported in our Brief Communication Arising article "Consistency in large pharmacogenomic studies" (Nature). All of this code is based on the code of Haibe-Kains et al. (www.nature.com/articles/nature12831), which was released as a supplementary file that accompanied their paper (http://www.nature.com/nature/journal/v504/n7480/extref/nature12831-s2.zip).

We offer two methods to reproduce our results; the first contained in the "pipeline" folder, offers absolute transparency, as it is a version of the original Haibe-Kains et al. code with the fewest possible modification and allows our results to be reproduced exactly using the modified full Haibe-Kains pipeline. Because this pipeline can be difficult to run and requires very large amounts of data to be downloaded, we also offer a simplified method to reproduce our results, contained in the folder "Simplified". This folder contain the .RData files (downloaded and created using the original pipeline) that are required to reproduce only our results, hence, using this pipeline, our results can be reproduced very easily.

Specific instructions on how to run pipelines:

# Simplified Pipeline

First, clone this repository (using git) or download the files. Next, to run the "simplified" pipeline, simply navigate to the "Simplified" folder, open an R terminal and run the commands in the "simplePipeline.R" script.

Note, the R libraries ggplot2 and survcomp are required by this script. You can install these using the following R code:

> source("http://bioconductor.org/biocLite.R")
> biocLite(c("ggplot2", "survcomp"))

Note: This pipeline has been tested in R version 3.2.0.

# Full pipeline

We have modified two of the files provided by Haibe-Kains et. al., namely  "CDRUG_correlation_data.R" and "CDRUG_analysis.R". In "CDRUG_correlation_data.R", our code can be found between lines 98 and 123. In  "CDRUG_analysis.R" our code can be found between lines 289 and 487. All other files and content are unedited and provided exactly as they were by Haibe-Kains et al.

To run these files and reproduce our results in their entirety, the steps taken are identical to those described in the original Haibe-Kains et al paper. Detailed instructions on how to do this are also included in the Supplementary Information of Haibe-Kains et al. (http://www.nature.com/nature/journal/v504/n7480/extref/nature12831-s1.pdf).

Once the full Haibe-Kains et al. pipeline has been run once in its entirety (i.e. all of the data files have been downloaded and prepared), it is possible to step through our code in an interactive session in the R-terminal: First, open the R terminal and set the working directory to the directory that contains all of the scripts (i.e. the "pipeline" folder). Then, run lines 1-78 of "CDRUG_pipeline.R" (to initialize some variables required by the other scripts) by pasting them into the R terminal. Next, paste line 1-98 of "CDRUG_correlation_data.R"; the next lines (98-123) of "CDRUG_correlation_data.R" will now allow you to reproduce the result for correlations in expression data "between" and "across" samples, so go ahead and paste them into the R-terminal. Next, open "CDRUG_analysis.R" and run all code from lines 1-289. The code from lines 289-487 of "CDRUG_analysis.R" will now allow you to reproduce the remainder of our results in an interactive R session.

Note that some of the files to be downloaded for the initial run of the Haibe-Kains et al. pipeline are very large (total 30GB), thus running the pipeline the first time can be slow. We also recommend using a machine with at least 16GB of RAM. Also, if the software environment is not set up exactly as described in the Haibe-Kains et al. Supplementary Information, errors may occur.





