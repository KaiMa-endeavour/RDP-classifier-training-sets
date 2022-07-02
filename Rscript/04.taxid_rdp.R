id_parent <- function(parent_name, file) {
  tag_split <- as.data.table(do.call('rbind', strsplit(file[!is.na(file)], split = '\\*')))
  if(parent_name %in% tag_split$V2) {
    id_p <- tag_split[which(parent_name == tag_split$V2), ]$V1
  } else {
    cat(parent_name, 'was not found.', '\n')
    break
  }
  id_p
}
taxid_rdp <- function(fasta, ref_object, nworker, output_path) {
  require(MicroEcoTk)
  require(tidyfst)
  require(stringr)
  require(parallel)
  options(warn = 2)
  fa <- read_fasta(fasta)
  header <- word(fa$Header, 2, -1, sep = ' ')
  items <- strsplit(header, ';')
  
  if (ref_object %in% c('Bacteria', 'Archaea', 'Fungi')) {
    ranks <- c('domain', 'phylum', 'class', 'order', 'family', 'genus')
  } else if (ref_object %in% 'Eukaryota') {
    ranks <- c('domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus')
  }
  ranks_depth <- data.table(rank = ranks, depth = 1:length(ranks)-1)
  
  cl <- makeCluster(nworker, type = "PSOCK") # PSOCK (windows); FORK (linux);
  items_table <- parLapply(cl, parLapply(cl, items, t), as.data.table) %>% rbindlist(fill = T)
  stopCluster(cl)
  colnames(items_table) <- ranks
  
  i = 1
  file <- matrix(NA, sum(sapply(items, length)), 1)
  for (r in 1:length(ranks)) {
    sons <- items_table[, c(get('r'))] %>% table() %>% names()
    rank <- ranks[r]
    
    for (n in 1:length(sons)) {
      son <- sons[n]
      if (!is.na(son)) {
        son_l <- ifelse(grepl('\\(', son), gsub('\\(', '\\\\(', son), son)
        son_l <- ifelse(grepl('\\)', son), gsub('\\)', '\\\\)', son_l), son)
        son_l <- ifelse(grepl('\\[', son), gsub('\\[', '\\\\[', son), son)
        son_l <- ifelse(grepl('\\]', son), gsub('\\]', '\\\\]', son_l), son)
        
        if (r>1) father <- items_table[which(son == unlist(items_table[, c(get('r'))])), c(get('r'))-1] %>% table() %>% names() else father <- NULL
        if(!son %in% c('d__Archaea', 'd__Bacteria', 'd__Eukaryota', 'd__Fungi')) id_p <- id_parent(parent_name = father, file) else id_p <- -1
        if (r < length(ranks)) {
          son_after <- unlist(items_table[which(son == unlist(items_table[, c(get('r'))])), c(get('r')+1)])
          if (any(is.na(son_after))) {
            if (!any(grepl(paste0('\\*', son_l, '\\*'), file) & grepl('\\*genus$', file))) {
              file[i, ] <- paste(i-1, son, min(id_p), ranks_depth$depth[which(ranks_depth$rank == rank)], 'genus', sep = '*')
              cat(file[i, ], '\n')
              i <- i + 1
            }
          }
        }
        if (length(son_after) > 0 & all(!is.na(son_after)) | (r == length(ranks))) {
          if (!any(grepl(paste0('\\*', son_l, '\\*'), file) & grepl(paste0('\\*', rank, '$'), file))) {
            file[i, ] <- paste(i-1, son, min(id_p), ranks_depth$depth[which(ranks_depth$rank == rank)], rank, sep = '*')
            cat(file[i, ], '\n')
            i <- i + 1
          }
        }
      }
    }
  }
  file <- file[!is.na(file)]
  file_non_redundant <- data.table(file[!duplicated(do.call('rbind', strsplit(file, '\\*'))[, -c(1, 3, 4)])])
  if(length(file) == nrow(file_non_redundant)) cat('\n', 'Taxid non redundancy.', '\n')
  fwrite(file_non_redundant, output_path, col.names = F)
}

library(optparse)
option_list <- list(
  make_option(opt_str = c('--fasta'), type = 'character', default = F), 
  make_option(opt_str = c('--ref_object'), type = 'character', default = F), 
  make_option(opt_str = c('--nworker'), type = 'integer', default = F), 
  make_option(opt_str = c('--output_path'), type = 'character', default = F)
)
args <- parse_args(OptionParser(option_list=option_list, usage = 'RDP taxid'))

taxid_rdp(fasta = args$fasta, ref_object = args$ref_object, nworker = args$nworker, output_path = args$output_path)
