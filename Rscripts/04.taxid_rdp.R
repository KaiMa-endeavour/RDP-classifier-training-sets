id_parent <- function(parent_name, file, depth) {
  file <- file[file$V1 != 'xx']
  if (any(grepl(paste0('\\*', parent_name, '\\*'), file$V1))) {
    id_p <- unlist(strsplit(file[grepl(paste0('\\*', parent_name, '\\*'), file$V1) & grepl(depth, file$V1)]$V1, split = '\\*'))[1]
  } else {
    cat(parent_name, 'was not found.', '\n')
    break
  }
  # tag_split <- as.data.table(do.call('rbind', strsplit(file[!is.na(file)], split = '\\*')))
  # if(parent_name %in% tag_split$V2) {
  #   id_p <- tag_split[which(parent_name == tag_split$V2), ]$V1
  # } else {
  #   cat(parent_name, 'was not found.', '\n')
  #   break
  # }
  id_p
}
# fasta <- './xx/mod_UNITE_public_29.11.2022.fa'
# output_path = './xx/trainset_taxid_Fungi.txt'
# ref_object <- 'Fungi'
# nworker = 6
taxid_rdp <- function(fasta, ref_object, nworker, output_path) {
  require(MicroEcoTk)
  require(tidyfst)
  require(stringr)
  require(parallel)
  options(warn = 2)
  fa <- read_fasta(fasta)
  # rank_class <- data.table(rank = tax$V1 %>% strsplit(';') %>% lapply(rev) %>% sapply('[', 1), calss = tax$V3, depth = tax$V1 %>% strsplit(';') %>% sapply(length))
  header <- word(fa$Header, 2, -1, sep = ' ')
  items <- strsplit(header, ';')
  
  if (ref_object %in% c('Bacteria', 'Archaea')) {
    ranks <- c('domain', 'phylum', 'class', 'order', 'family', 'genus')
  } else if (ref_object %in% 'Eukaryota') {
    ranks <- c('domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus')
    # ranks <- c('domain', 'superkingdom', 'kingdom', 'subkingdom', 'superphylum', 'phylum', 'subphylum', 'infraphylum', 'superclass', 'class', 'subclass', 'infraclass', 'superorder', 'order', 'suborder', 'superfamily', 'family', 'subfamily', 'genus')
  } else if (ref_object %in% 'Fungi') {
    ranks <- c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
  }
  ranks_depth <- data.table(rank = ranks, depth = 1:length(ranks)-1)
  
  cl <- makeCluster(nworker, type = "PSOCK") # PSOCK (windows); FORK (linux);
  items_table <- parLapply(cl, parLapply(cl, items, t), as.data.table) %>% rbindlist(fill = T)
  stopCluster(cl)
  colnames(items_table) <- ranks
  # items_table <- items_table[1:3145]
  
  i = 1
  file <- matrix('xx', length(unique(unlist(items)))*2, 1) %>% as.data.table()
  for (r in 1:length(ranks)) {
    # r = 7
    sons <- items_table[, c(get('r'))] %>% table() %>% names()
    rank <- ranks[r]
    
    for (n in 1:length(sons)) {
      # n = 805
      son <- sons[n]
      if (!is.na(son)) {
        son_l <- ifelse(grepl('\\(', son), gsub('\\(', '\\\\(', son), son)
        son_l <- ifelse(grepl('\\)', son), gsub('\\)', '\\\\)', son_l), son)
        son_l <- ifelse(grepl('\\[', son), gsub('\\[', '\\\\[', son), son)
        son_l <- ifelse(grepl('\\]', son), gsub('\\]', '\\\\]', son_l), son)
        
        if (r>1) father <- items_table[which(son == unlist(items_table[, c(get('r'))])), c(get('r'))-1] %>% table() %>% names() else father <- NULL
        if(!son %in% c('d__Archaea', 'd__Bacteria', 'd__Eukaryota', 'k__Fungi')) id_p <- id_parent(parent_name = father, file, depth =  ranks_depth$depth[ranks_depth$rank==rank]-1) else id_p <- -1
        if (r < length(ranks)) {
          son_after <- unlist(items_table[which(son == unlist(items_table[, c(get('r'))])), c(get('r')+1)])
          if (any(is.na(son_after))) {
            if (!any(grepl(paste0('\\*', son_l, '\\*'), file) & grepl(paste0('\\*', rev(ranks)[1], '$'), file))) {
              file[i, ] <- paste(i-1, son, min(id_p), ranks_depth$depth[which(ranks_depth$rank == rank)], rev(ranks)[1], sep = '*')
              cat(file[i, ]$V1, '\n')
              i <- i + 1
            }
          }
        }
        if ((length(son_after) > 0 & all(!is.na(son_after))) | (r == length(ranks))) {
          if (!any(grepl(paste0('\\*', son_l, '\\*'), file) & grepl(paste0('\\*', rank, '$'), file))) {
            file[i, ] <- paste(i-1, son, min(id_p), ranks_depth$depth[which(ranks_depth$rank == rank)], rank, sep = '*')
            cat(file[i, ]$V1, '\n')
            i <- i + 1
          }
        }
      }
    }
  }
  file <- file[file$V1 != 'xx']
  file_non_redundant <- data.table(file[!duplicated(do.call('rbind', strsplit(file$V1, '\\*'))[, -c(1, 3, 4)])])
  if(length(file$V1) == nrow(file_non_redundant)) cat('\n', 'Taxid non redundancy.', '\n')
  # fwrite(data.table(file), output_path, col.names = F)
  fwrite(file_non_redundant, output_path, col.names = F)
}

