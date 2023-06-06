library(tidyfst)


taxid <- function(tax) {
  rank <- tax$V1 %>% strsplit(';') %>% lapply(rev) %>% sapply('[', 1)
  class_rank <- paste(substr(tax$V3, 1, 1), rank, sep = '__')
  rank_class <- data.table(rank = rank, class_rank = class_rank, class = tax$V3, depth =  sapply(strsplit(tax$V1, ';'), length), tax = tax$V1)
}

tax <- fread('G:/Database/silva/tax/tax_slv_ssu_138.txt.gz')

tax_Bac <- tax[which(tax$V1 %>% strsplit(';') %>% sapply('[[', 1) == 'Bacteria'), ]
rank_class_bac <- taxid(tax_Bac)
fwrite(rank_class_bac, 'tax_138_Bacteria.txt', sep = '\t')

tax_Arc <- tax[which(tax$V1 %>% strsplit(';') %>% sapply('[[', 1) == 'Archaea'), ]
rank_class_arc <- taxid(tax = tax_Arc)
fwrite(rank_class_arc, 'tax_138_Archaea.txt', sep = '\t')

tax_Euk <- tax[which(tax$V1 %>% strsplit(';') %>% sapply('[[', 1) == 'Eukaryota'), ]
rank_class_euk <- taxid(tax = tax_Euk)
fwrite(rank_class_euk, 'tax_138_Eukaryota.txt', sep = '\t')
