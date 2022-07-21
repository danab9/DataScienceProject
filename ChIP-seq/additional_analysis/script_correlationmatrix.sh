# for creation of correlation overview between the different classes
# starting from the merged bed files:


#transform bed file into bam
bedtools2/bin/bedToBam -i MTL.H3K9ac.bed -g bedtools2/genomes/human.hg19.genome > MTL.H3K9ac.bam
bedtools2/bin/bedToBam -i MTL.H3K4me1.bed -g bedtools2/genomes/human.hg19.genome > MTL.H3K4me1.bam
bedtools2/bin/bedToBam -i MTL.H3K27ac.bed -g bedtools2/genomes/human.hg19.genome > MTL.H3K27ac.bam
bedtools2/bin/bedToBam -i MTL.H3K122ac.bed -g bedtools2/genomes/human.hg19.genome > MTL.H3K122ac.bam

#index bam file (required for next step
samtools index MTL.H3K9ac.bam
samtools index MTL.H3K4me1.bam
samtools index MTL.H3K27ac.bam
samtools index MTL.H3K122ac.bam

# run multiBamSummary
multiBamSummary bins --bamfiles MTL.H3K9ac.bam MTL.H3K4me1.bam MTL.H3K27ac.bam MTL.H3K122ac.bam -o results_MTLs.npz

#create plot with plotCorrelation
plotCorrelation --corData results_MTLs.npz -c pearson -p heatmap --outFileCorMatrix SpearmanCorr_readCounts.tab -o heatmap_SpearmanCorr_readCounts.png --plotTitle "Spearman Correlation of Read Counts" --colorMap RdYlBu --plotNumbers
