RNA-seq analysis identifies upregulated and downregulated genes. For this purpose,
DESeq2 is used in this analysis and the results are compared with RUSH and the results
of the paper. For further analysis, RNA was normalized (per gene length).  
The results were further visualized and in addition a linear model was performed 
with the nomalized RNA and the clinical data. 

For the plots created in Rstudio, you can find a description under the code lines for the plots.

RNA seq analysis code:

Clinical distribution in the tree groups (young, old and AD)
- gender
- age
- cause of death
-> Clinicaldistribution.R

Differential gene expression analysis:
- read RNAseq dataset
- grouping
- differential genexpression AD vs old with DESeq2
- up- and downregulated genes 
- quick comparison to results paper (no visualization)
- visualization result differential expression:
	- MA plot
	- PCA
	- different Heatmaps
- differential expression with other groups

	Output: 
	- upregulated.txt, downreguated.txt
	- Only gene names: upregulatedgenes.txt, downregulatedgenes.txt
	- PCARNAseq.png, MAplot.png

-> differentialexpression.R and differentialexpression_function.R

- Comparision own and paper/rush results differential expression:
	 	
	Output: VennUpregualted.png, VennDownregulated.png
-> Comp_OWNvsPAPER.R and Comp_MAINvsRUSH


RNA normalization:
- ID convertion: hngc_symbols to ensemble ID
- used genelength for normalization (function getGeneLengthAndGCContent)	
- actual function: fpkm
- comparison RNA and normalized RNA
- boxplot of normalized RNA in different groups
	- with kruskal-wallis or wilcox test

	Output: normalizedRNA.txt
 
-> normalizedRNA.R


Visualization of GO term (David) "regulation of transcription" with different heatmaps:
	
	Output: 
	- RT_heatmap.png
	- RT_heatmapNew.png
	- GOenrichment.png

-> GODavid_differentialexpression.R

STRING analysis of upregulated regulation of transcription genes:
- visualization in cytoscape:
	- calculation of differential gene expression changes (AD vs old) for color intensity of nodes
	- calculation of mean gene expression in AD for nodesize
- visualization of the number of interactions from genes with the most interactions

	Output: String_interactions.png 

- Other figures: 
	- STRING network.png: Interaction network of STRING analysis visualized in cytoscape. 
	Reveald an interaction network of 75 genes and SIRT1, EP300 and CREBBP together with RXRA
	are located in the center of the network. 52 genes don't have an interaction with other genes.
	Nodesize: mean expression in AD
	Color intensity: log2-fold (AD vs old)

Relating groups (AD vs old) to gene expression variation:
- normalized RNA
- use of ExpressionSet
- linear model (lmFit)
	Output: lmFit.png
->lmfit.R