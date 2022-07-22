ChIPseq identifies DNA binding sites of transcription factors and other proteins. 
We had data for binding sites of four histone modifications (H3K9ac, H3K27ac, H3K122ac, H3K4me1).
In the following, the code in the GitHub is described which we used for analyzing the ChIPseq data, which unfortunately had bad quality.

ChIPseq Analyses Code and files:
The code is divided into three main folders.
Since the DiffBind analysis was a bit bigger we have an extra README/ notebook file in this folder.


DiffBindanalysis: differential binding analysis - see own README file


ML_matrix_creation: creation of iClusterBayes input matrix
- patient_merge.sh: bash script for creating the big unfiltered matrix for the ML technique
- iCluster_matrix_creation.R: R script for filtering matrix for ML approach
- bedfiles_analysis_iCluster_tryout: for analysis of matrix - testing around
- iCluster_mat_input.txt: final input matrix for iClusterBayes
- correlationmatrix_iCluster_input.R: script for creation of correlation matrix of input matrix
- correlationmatrix_iCluster_input.png: picture of correlation matrix of input matrix


additional_analysis: other analyses
- MIT.AD_H3K9ac_merged.bed, MIT.Y_H3K9ac_merged.bed, MIT.O_H3K9ac_merged.bed: required for creation of Venn diagram showing overlap of regions between different patients
- Overlap_HistMarks.R: script for creating Venn diagram for howing range overlap between three groups
- Overlap_HistMarks.png: Venn Diagram showing result from R script (see one line above)
- MIT_file_creation.sh: shell script for creating MIT files which are input for different tools
- combine_files_to_MTLs_script.sh: shell script for combining MIT files into MTL files - required for different things (basically patient files merged based on histone mark)
- script_correlationmatrix.sh: for creation of correlation matrix of the 4 MTL files - for getting overview over similarity (using MTL files)
- heatmap_PearsonCorr_readCounts.png: from script one line above - Pearson correlation between MTL files
- heatmap_SpearmanCorr_readCounts.png: same as one above but using Spearman Correlation coefficient
- peak_sized_ChIPs.R: script for getting overview over length of ranges in MTL files
- peak_sized_ChIPs.pdf: resulting figure from script above
[- we did some additional statistical tests and overviews such as histograms but did not write them in a script]