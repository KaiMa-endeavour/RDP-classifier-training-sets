library(MicroEcoTk)
library(tidyfst)


fa <- read_fasta('G:/Database/silva/Silva_138.fa')
domain <- fa$Header %>% strsplit(';') %>% sapply('[[', 1) %>% strsplit(' ') %>% sapply('[[', 2)

fa_bac <- fa[which(domain == 'Bacteria'), ]
write_fasta(fa_bac, 'SILVA_138_Bacteria.fa')

fa_arc <- fa[which(domain == 'Archaea'), ]
write_fasta(fa_arc, 'SILVA_138_Archaea.fa')

fa_euk <- fa[which(domain == 'Eukaryota'), ]
write_fasta(fa_euk, 'SILVA_138_Eukaryota.fa')
