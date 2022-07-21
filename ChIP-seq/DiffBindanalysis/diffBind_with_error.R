library(DiffBind)

file <- "samplesheet_diffbind.txt"

#read in samples
samples <- read.table(file, sep="\t")
colnames(samples) <- c("SampleID", "Peaks", "bamReads", "PeakCaller", "Factor")


H3K9ac <- dba(sampleSheet=samples)

# correlation plots of diffBind between three classes Y, O and AD for H3K9ac
plot(H3K9ac)
# Y and O more similar than to AD, also O and AD more correlated than Y to AD which makes sense


H3K9ac_blacklist <- dba.blacklist(H3K9ac)

#count reads
H3K9ac_count <- dba.count(H3K9ac_blacklist)


#normalize the data
H3K9ac_norm <- dba.normalize(H3K9ac_count)


#calculate contrast so DiffBind knows how to model data
H3K9ac_contrast <- dba.contrast(H3K9ac_norm, reorderMeta = list(Factor="AD"), minMembers=1)
#not working! more groups are required/ with more replicates!


#perform differential analysis - DID NOT WORK - tried out different ideas on how to solve but did not get it running -> performed own diffBind analysis
# error from DESeq2, which is called within dba.analyze, saying that their only zero counts for every gene. possible reason might be the bam files, that are recovered from the bed files and thus lost a lot of information
H3K9ac_analyzed <- dba.analyze(H3K9ac_norm)


