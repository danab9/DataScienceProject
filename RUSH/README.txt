RNA-seq analysis - RUSH 

there are three R scripts: 
- prepare_rushdata.R: 
	this script is dealing with importing, creating overivew files and ID conversion ensembl to hgnc symbols
- direct_comparsion_RUSH_main.R:
	attempt of a combined analysis of rush and nativio et al. data including also batch correction methods 
- RUSH_RNAseq_analysis.Rmd:
	Data exploration (FC comparison,PCA,heatmap,...)  is done here as well as DESeq2 analysis 


Furthermore, the folder "plots" contains all plots, "results" contains the diff. expressed genes and "overviews" contains all sample infos.

