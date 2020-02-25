To run the code contained in this repository, we recommend the following:

1. Download supplementary dataset 6 at https://www.nature.com/articles/nbt.3981, and ensure that the 
downloaded file is named "nbt.3981-S9.zip" in your working directory.

2. Download the integrated OTU table available at http://downloads.ihmpdcc.org/data/MBQC/mbqc_integrated_otus.tsv.gz" and
ensure that the downloaded fileis named "mbqc_integrated_otus.tsv.gz" in your working directory

3. Ensure that you have installed the R libraries this code requires:

	xgboost
	magrittr
	Matrix
	tidyr
	dplyr
	RDS
	data.table
	batchtools
	ggplot2
	knitr
	gridExtra
	RColorBrewer
	
4. Create or download appropriate configuration and template files for batchtools, a cluster computing interface for R. We have provided files 
appropriate for a Sun Grid Engine environment (note that two template files are necessary, as some cluster jobs recruit multiple
nodes for faster computing). While it is possible to run this code outside of a parallel environment, we expect computation to be 
prohibitively time consuming in this case.
