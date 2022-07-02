fa_mod <- function(fasta, nworker, output_path) {
  require(MicroEcoTk)
  require(tidyfst)
  require(stringr)
  require(doParallel)
  
  fa <- read_fasta(fasta)
  fa$Sequence <- toupper(fa$Sequence)
  fa <- fa[!duplicated(fa$Sequence), ]
  
  fa$Header <- gsub('\\"', '', fa$Header)
  lineage <- strsplit(fa$Header, 'Lineage=') %>% sapply('[', 2) %>% word(3, -1, sep = ';')
  header_l <- strsplit(lineage, ';')
  ranks_num <- max(sapply(header_l, length))
  ranks <- c('domain', 'phylum', 'class', 'order', 'family', 'genus')
  
  famod_item <- function(seq) {
    xx <- rep(12344321, length(ranks))
    for (n in 1:(ranks_num/2)) {
      rank <- sapply(header_l[seq], '[', 2*n)
      if (rank %in% ranks) {
        item <- sapply(header_l[seq], '[', 2*n-1)
        xx[which(rank == ranks)] <- paste(substr(rank, 1, 1), item, sep = '__')
      }
    }
    xx
  }
  cl <- makeCluster(nworker, type = "PSOCK") # PSOCK (windows); FORK (linux);
  registerDoParallel(cl)
  header_table <- foreach(seq = 1:nrow(fa), .packages = 'stringr', .combine = 'rbind', .verbose = T) %dopar% famod_item(seq)
  stopCluster(cl)
  header_table[header_table == 12344321] <- ''
  colnames(header_table) <- ranks
  
  cat('\n')
  if(!dir.exists(output_path)) dir.create(output_path)
  
  for (i in 2:ncol(header_table)) {
    item_c <- names(table(header_table[, c(get('i'))]))[names(table(header_table[, c(get('i'))])) != '']
    for (r in 1:length(item_c)) {
      f_num <- table(header_table[which(header_table[, c(get('i'))] == item_c[r]), c(get('i')-1)])
      if (length(f_num) != 1) {
        header_table[, c(get('i'))] <- gsub(pattern = paste0('^', item_c[r], '$'), replacement = '', unlist(header_table[, c(get('i'))]))
        cat('Delete', item_c[r], '\n')
      }
    }
    names(table(header_table[, c(get('i'))]))
    yes_no <- header_table[, c(get('i'))] == ''
    
    if(any(yes_no)) {
      sons_id <- which(yes_no)
      fathers <- header_table[, c(get('i'))-1]
      table(fathers)
      for (n in sons_id) {
        father_id <- which(fathers == header_table[n, c(get('i'))-1])
        if (all(header_table[father_id, c(get('i'))] == '') & all(header_table[n, c(get('i')):ncol(header_table)] == '')) {
          cat('Row ID:', n, 'ok! ', header_table[n, c(get('i'))-1], '\n')
        } else {
          xx <- paste(substr(ranks[i], 1, 1), paste('NA', sapply(strsplit(header_table[n, c(get('i'))-1], '__'), '[', 2), sep = '.'), sep = '__')
          header_table[father_id, c(get('i'))] [header_table[father_id, c(get('i'))] == ''] <- xx
          cat('\n', 'Row ID:', n, 'add:', xx, '\n', '\n')
        }
      }
    }
  }
  
  header <- sapply(1:nrow(fa), function(n) paste(unlist(header_table[n, ]), collapse = ';'))
  header <- gsub('\\;+$', '', header)
  header <- gsub('\\;+', '\\;', header)
  
  fa$Header <- paste(word(fa$Header, 1, end = 1, sep = ' '), header, sep = ' ')
  fwrite(data.table(`Seq ID` = word(fa$Header, 1, end = 1), header_table), paste0(output_path, 'header_table.csv'))
  fasta <- strsplit(fasta, '\\/') %>% sapply(nth, -1)
  write_fasta(fa, paste0(output_path, fasta))
}

library(optparse)
option_list <- list(
  make_option(opt_str = c('--fasta'), type = 'character', default = F), 
  make_option(opt_str = c('--nworker'), type = 'integer', default = F), 
  make_option(opt_str = c('--output_path'), type = 'character', default = F)
)
args <- parse_args(OptionParser(option_list=option_list, usage = 'RDP fa mod'))

fa_mod(fasta = args$fasta, nworker = args$nworker, output_path = args$output_path)
