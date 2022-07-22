# 
library('org.Hs.eg.db')
library(dplyr)
library(stringr)

# load Nativio et. al. RNA-seq STAR results
rnaseqdata <- read.table("../RNA-seq/GSE153873_summary_count.star.txt", sep = "\t",header = TRUE,row.names = 1)
colnames(rnaseqdata) <- paste0(data.frame(str_split(colnames(rnaseqdata),'[.]'))[3,],data.frame(str_split(colnames(rnaseqdata),'[.]'))[1,])

# get all hgnc symbols from old reference to new reference name, to search again
# for ensembl ids.
# (if there is no ensemble ID, one can not normalize the gene by it's length)
library(EDASeq)
translateSymbolToID <- function(symbols_list) {
  # get entrez id
  ensembl <- mapIds(org.Hs.eg.db, symbols_list, 'ENSEMBL', 'SYMBOL')
  
  ids<- cbind(as.data.frame(symbols_list), as.data.frame(ensembl))#####
  no_ids<- ids[which(is.na(ids$ensembl)), ]  # check with G list
  
  mart <- useMart(dataset="hsapiens_gene_ensembl","ensembl")
  G_list <- getBM(filters= "hgnc_symbol", 
                  attributes=c('hgnc_symbol','ensembl_gene_id'),
                  values= no_ids[,1],
                  mart= mart)
  
  G_list_foundids <- G_list[which(G_list$ensembl_gene_id != ''),]
 
  # get rid of duplicates  
  # there are not much in our case. Should be checked specifically for other cases.
  G_list_foundids_nodups <- G_list_foundids %>% group_by(hgnc_symbol) %>% filter(row_number()==1)
  
  colnames(ids) <- colnames(G_list_foundids_nodups)
  result <- rbind(G_list_foundids_nodups, ids)
  return(result)
}

# get all possible new symbols
general_symbols <- translateSymbolToID(rownames(rnaseqdata))
general_symbols_notfound <- general_symbols[which(general_symbols$ensembl_gene_id %in% c('',NA)),] # not found ensmbl ids

# get approved symbols for non-found ids 
symbolscheck <- read.csv('data/hgnc-symbol-check.csv', skip = 1, header = TRUE)
symbolscheck <- symbolscheck[,c(1,3)]  # to do: get from online source?
has_approved_symbol <- symbolscheck %>% filter(!Approved.symbol %in% c(NA, '')) 
appSymbols_ids <- translateSymbolToID(has_approved_symbol[,2])

has_approved_symbol <- has_approved_symbol %>% rename(hgnc_symbol=Approved.symbol) # change col name 
has_approved_symbol <- has_approved_symbol %>% left_join(appSymbols_ids)
dim(has_approved_symbol)
has_approved_symbol <- has_approved_symbol %>% filter(!ensembl_gene_id %in% c(NA, ''))
dim(has_approved_symbol) # 129 symbols; cols: Input(old symbol), hgnc_symbol(new), ensembl_gene_id


general_symbols_notfound_new <- data.frame(general_symbols_notfound[which(general_symbols_notfound$hgnc_symbol %in% has_approved_symbol$Input),])
general_symbols_notfound_new <- rename(general_symbols_notfound_new, Input=hgnc_symbol)
has_approved_symbol <- has_approved_symbol[!duplicated(has_approved_symbol$Input),] # remove double inputs

general_symbols_notfound_new <- general_symbols_notfound_new %>% left_join(has_approved_symbol, by="Input")
general_symbols_notfound_new <- general_symbols_notfound_new %>% select (c("Input","ensembl_gene_id.y"))
general_symbols_notfound_new <- general_symbols_notfound_new[,c("Input", "ensembl_gene_id.y")]
colnames(general_symbols_notfound_new) <- c("hgnc_symbol", "ensembl_gene_id")



#### RNA-seq Normalization with all symbols ####
library(EDASeq)

# get lengths of all samples 
general_symbols_copy <- general_symbols[which(!general_symbols$ensembl_gene_id %in% c('',NA)),] 
general_symbols_found <- rbind(general_symbols_copy, general_symbols_notfound_new)
dim(general_symbols_found)
general_symbols_found <- general_symbols_found[!duplicated(general_symbols_found$hgnc_symbol),] # remove duplicates

# get entrez id
ensembl <- mapIds(org.Hs.eg.db, rownames(rnaseqdata), 'ENSEMBL', 'SYMBOL')
ids<- cbind(as.data.frame(rownames(rnaseqdata)), as.data.frame(ensembl))#####
ids<- ids[which(is.na(ids$ensembl)), ]
ensemblnew <- ensembl[!is.na(ensembl)]
genelength <- getGeneLengthAndGCContent(id=ensemblnew,'hsa')


general_symbols_notfound_vector <- c(general_symbols_notfound_new$ensembl_gene_id)
genekength_notfound_list <- getGeneLengthAndGCContent(general_symbols_notfound_vector, "hsa")

# combine both:
all_lengths <- data.frame(rbind(genelength, genekength_notfound_list))
all_lengths$ensembl_gene_id <- rownames(all_lengths)
all_lengths <- all_lengths %>% left_join(as.data.frame(general_symbols_found))
found_ids_rnaseq <- rnaseqdata[which(rownames(rnaseqdata) %in% all_lengths$hgnc_symbol),]
all_lengths <- all_lengths[which(all_lengths$hgnc_symbol %in% rownames(rnaseqdata)),]

# normalize again
condition <- factor(colnames(rnaseqdata))
deseqrna <- DESeqDataSetFromMatrix(found_ids_rnaseq, DataFrame(condition), ~ condition)
mcols(deseqrna)$basepairs <- all_lengths$length

# Counts normalized per kilobase
normalizedRNAseq <- log2(fpkm(deseqrna, robust = TRUE)+1)
write.csv(normalizedRNAseq, 'results/normalizedRNA_new.tsv', dec = ".",
          row.names = TRUE, col.names = TRUE)



