library(tidyfst)


taxid <- function(tax) {
  
  rank <- tax$V1 %>% strsplit(';') %>% lapply(rev) %>% sapply('[', 1)
  
  class_rank <- paste(substr(tax$V3, 1, 1), rank, sep = '__')
  # father <- tax$V1 %>% strsplit(';') %>% lapply(rev) %>% sapply('[', 2)
  # for (p in 1:length(table(son))) {
  #   if (length(table(father[which(son == names(table(son))[p])])) != 1 & !names(table(son))[p] %in% c('Archaea', 'Bacteria', 'Eukaryota')) {
  #     cat(names(table(son))[p], '\t')
  #     for (r in id <- which(son == names(table(son))[p])) {
  #       son[r] <- paste(son[r], father[r], sep = '_')
  #     }
  #   }
  # }
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

# taxid <- function(tax) {
#   file <- matrix(NA, nrow(tax), 1)
#   # file[1, ] <- '0*Root*-1*0*rootrank'
#   for (r in 1:(nrow(tax))) {
#     cat(r, '\n')
#     tag <- unlist(tax[r, 1]) %>% strsplit(split = ';') %>% unlist() %>% rev()
#     if(!tag[1] %in% c('Archaea', 'Bacteria', 'Eukaryota')) id_p <- id_parent(parent_name = tag[2], file) else id_p <- -1
#     file[r, ] <- paste(r-1, tag[1], rev(id_p)[1], length(tag)-1, tax[r, 3], sep = '*')
#   }
#   dup_row_id <- which(duplicated(do.call('rbind', strsplit(file, '\\*'))[, -c(1, 3, 4)]))
#   if (length(dup_row_id) != 0) {
#     for (r in dup_row_id) {
#       cat(r, '\n')
#       tag <- unlist(tax[r, 1]) %>% strsplit(split = ';') %>% unlist() %>% rev()
#       tag[1] <- paste(tag[1], tag[2], sep = '_')
#       if(!tag[1] %in% c('Archaea', 'Bacteria', 'Eukaryota')) id_p <- id_parent(parent_name = tag[2], file) else id_p <- -1
#       file[r, ] <- paste(r-1, tag[1], rev(id_p)[1], length(tag)-1, tax[r, 3], sep = '*')
#     }
#   }
#   
#   no_genus_id <- which(do.call('rbind', strsplit(file, '\\*'))[, 5] != 'genus')
# 
# 
#   # str_sub(unlist(tax[1, 1]), 1, -2) %in% word(unlist(tax[, 1]), 1, -2, sep = ';')
#   # 'Archaea;Nanoarchaeota;Nanoarchaeia;Woesearchaeales' %in% word(unlist(tax[, 1]), 1, -2, sep = ';')
#   no_genus_noson_id <- no_genus_id[which(!do.call('rbind', strsplit(file, '\\*'))[no_genus_id, 1] %in% do.call('rbind', strsplit(file, '\\*'))[no_genus_id, 3])]
#   for (r in no_genus_noson_id) {
#     cat(r, '\n')
#     tag <- unlist(tax[r, 1]) %>% strsplit(split = ';') %>% unlist() %>% rev()
#     if(!tag[1] %in% c('Archaea', 'Bacteria', 'Eukaryota')) id_p <- id_parent(parent_name = tag[2], file) else id_p <- -1
#     file[r, ] <- paste(r-1, tag[1], rev(id_p)[1], length(tag)-1, 'genus', sep = '*')
#   }
#   # dup_row_id <- which(duplicated(do.call('rbind', strsplit(file, '\\*'))[, -c(1, 3, 4)]))
#   # if (length(dup_row_id) != 0) {
#   #   for (r in dup_row_id) {
#   #     cat(r, '\n')
#   #     tag <- unlist(tax[r, 1]) %>% strsplit(split = ';') %>% unlist() %>% rev()
#   #     tag[1] <- paste(tag[1], tag[2], sep = '_')
#   #     if(!tag[1] %in% c('Archaea', 'Bacteria', 'Eukaryota')) id_p <- id_parent(parent_name = tag[2], file) else id_p <- -1
#   #     file[r, ] <- paste(r-1, tag[1], rev(id_p)[1], length(tag)-1, 'genus', sep = '*')
#   #   }
#   # }
#   # file <- as.data.table(file[which(!duplicated(do.call('rbind', strsplit(file, '\\*'))[, -c(1, 3, 4)])), ])
#   file
# }
# 
# taxid_bac <- taxid(tax = tax_Bac)
# fwrite(taxid_bac, 'trainset_taxid_Bacteria.txt', col.names = F)
# taxid_arc <- taxid(tax = tax_Arc)
# fwrite(taxid_arc, 'trainset_taxid_Archaea.txt', col.names = F)
# taxid_euk <- taxid(tax = tax_Euk)
# fwrite(taxid_euk, 'trainset_taxid_Eukaryota.txt', col.names = F)
