fa_mod <- function(fasta, rank, ref_object, nworker, output_path) {
  require(MicroEcoTk)
  require(tidyfst)
  require(stringr)
  require(doParallel)
  
  fa <- read_fasta(fasta)
  tax <- fread(rank, na.strings = '') %>% as.data.frame
  fa$Header <- word(fa$Header, 1, end = -2, sep = ';')
  fa$Sequence <- gsub('U', 'T', fa$Sequence)
  fa <- fa[!duplicated(fa$Sequence), ]
  header <- word(fa$Header, 2, -1, sep = ' ')
  
  if (ref_object %in% c('Bacteria', 'Archaea')) {
    ranks <- c('domain', 'phylum', 'class', 'order', 'family', 'genus')
  } else if (ref_object %in% 'Eukaryota') {
    ranks <- c('domain', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus')
    # ranks <- c('domain', 'superkingdom', 'kingdom', 'subkingdom', 'superphylum', 'phylum', 'subphylum', 'infraphylum', 'superclass', 'class', 'subclass', 'infraclass', 'superorder', 'order', 'suborder', 'superfamily', 'family', 'subfamily', 'genus')
    tax <- tax[which(tax$class %in% ranks), ]
  }
  
  header_l <- strsplit(header, ';')
  famod_item <- function(seq) {
    # seq=1
    s <- header_l[[seq]]
    xx <- rep(12344321, length(ranks))
    for (n in 1:length(s)) {
      # n=7
      # s[which(sapply(1:length(s), function(i) tax$class[which(tax$rank == s[i])]) %in% s[n])]
      # print(ranks)
      # & tax$class == ranks[n]
      item_rank <- tax[tax$rank == s[n], ]$class
      item_id <- min(as.numeric(word(match(item_rank, ranks), 1, -1, sep = '')))
      if (n < length(s)) {
        item_after_rank <- word(tax[tax$rank == s[n+1], ]$class, 1, -1, sep = '')
        item_after_id <- match(item_after_rank, ranks)
        if (length(item_after_id) == 1) {
          if (is.na(item_after_id)) item_after_id <- item_id + 1
        } else {
          item_after_id <- min(item_after_id)
        }
      } else {
        item_after_id <- item_id + 1
      }
      if (!is.na(item_id)) {
        if (item_id < item_after_id) {
          item_table <- tax[tax$rank == s[n], ]
          item <- item_table[item_table$depth == min(item_table$depth), ]$class_rank[1]
          xx[c(get('item_id'))] <- item
        }
      }
    }
    # if (any(which(xx == '12344321') < which(xx != '12344321'))) xx[min(which(which(xx == '12344321') < which(xx != '12344321'))):length(xx)] <- '12344321'
    xx
  }
  cl <- makeCluster(nworker, type = "PSOCK") # PSOCK (windows); FORK (linux);
  registerDoParallel(cl)
  header_table <- foreach(seq = 1:nrow(fa), .packages = 'stringr', .combine = 'rbind', .verbose = T) %dopar% famod_item(seq)
  stopCluster(cl)
  header_table[header_table == 12344321] <- ''
  # header_table <- matrix('', nrow(fa), length(ranks))
  colnames(header_table) <- ranks
  # for (seq in 1:nrow(fa)) {
  #   # seq = 1710
  #   s <- header_l[[seq]]
  #   for (n in 1:length(s)) {
  #     # n=7
  #     # s[which(sapply(1:length(s), function(i) tax$class[which(tax$rank == s[i])]) %in% s[n])]
  #     # print(ranks)
  #     # & tax$class == ranks[n]
  #     item_rank <- tax[tax$rank == s[n], ]$class
  #     item_id <- min(as.numeric(word(match(item_rank, ranks), 1, -1, sep = '')))
  #     if (n < length(s)) {
  #       item_after_rank <- word(tax[tax$rank == s[n+1], ]$class, 1, -1, sep = '')
  #       item_after_id <- match(item_after_rank, ranks)
  #       if (length(item_after_id) == 1) {
  #         if (is.na(item_after_id)) item_after_id <- item_id + 1
  #       } else {
  #         item_after_id <- min(item_after_id)
  #       }
  #     } else {
  #       item_after_id <- item_id + 1
  #     }
  #     if (!is.na(item_id)) {
  #       if (item_id < item_after_id) {
  #         item_table <- tax[tax$rank == s[n], ]
  #         item <- item_table[item_table$depth == min(item_table$depth), ]$class_rank[1]
  #         header_table[seq, c(get('item_id'))] <- item
  #       }
  #     }
  #   }
  #   cat('sequence', paste0(seq, '/', nrow(fa)), '----', header_table[seq, ], '\n')
  # }
  cat('\n')
  if(!dir.exists(output_path)) dir.create(output_path)
  # header_table <- strsplit(header, ';') %>% sapply(rbind) %>% sapply(as.data.table) %>% rbindlist(fill = T)
  # header_table <- as.data.frame(rbindlist(sapply(sapply(strsplit(header, ';'), rbind), as.data.table), fill = T))
  tax_add <- matrix(NA, nrow(header_table)*ncol(header_table), 5)
  tax_id <- 1
  for (i in 2:ncol(header_table)) {
    # i=4
    # for (r in 1:nrow(sub_tax)) {
    #   header_table[, c(get('i'))] <- gsub(pattern = sub_tax$rank[r], replacement = sub_tax$class_rank[r], unlist(header_table[, c(get('i'))]))
    # }
    item_c <- names(table(header_table[, c(get('i'))]))[names(table(header_table[, c(get('i'))])) != '']
    for (r in 1:length(item_c)) {
      # r=26
      f_num <- table(header_table[which(header_table[, c(get('i'))] == item_c[r]), c(get('i')-1)])
      if (length(f_num) != 1) {
        header_table[, c(get('i'))] <- gsub(pattern = paste0('^', item_c[r], '$'), replacement = '', unlist(header_table[, c(get('i'))]))
        cat('Delete', item_c[r], '\n')
      }
    }
    names(table(header_table[, c(get('i'))]))
    yes_no <- header_table[, c(get('i'))] == ''
    
    if(any(yes_no)) {
      # sons_id <- setdiff(which(yes_no), sons_id)
      sons_id <- which(yes_no)
      fathers <- header_table[, c(get('i'))-1]
      table(fathers)
      for (n in sons_id) {
        # n = 6
        # & (header_table[n, c(get('i'))-1] != ''
        father_id <- which(fathers == header_table[n, c(get('i'))-1])
        if (all(header_table[father_id, c(get('i'))] == '') & all(header_table[n, c(get('i')):ncol(header_table)] == '')) {
          cat('Row ID:', n, 'ok! ', header_table[n, c(get('i'))-1], '\n')
        } else {
          xx <- paste(substr(ranks[i], 1, 1), paste('NA', sapply(strsplit(header_table[n, c(get('i'))-1], '__'), '[', 2), sep = '.'), sep = '__')
          header_table[father_id, c(get('i'))] [header_table[father_id, c(get('i'))] == ''] <- xx
          tax_add[tax_id, ] <- c('NA', xx, ranks[i], i, paste(header[n], xx, sep = ';'))
          tax_id <- tax_id + 1
          cat('\n', 'Row ID:', n, 'add:', xx, '\n', '\n')
        }
      }
    }
    # header_table[, c(get('i'))][which(is.na(header_table[, c(get('i'))]))] <- paste(substr(tax$class[tax$rank == e[c(get('i'))]], 1, 1), 'NA', sep = '__')
  }
  colnames(tax_add) <- c('rank', 'class_rank', 'class', 'depth', 'tax')
  tax_add <- tax_add[which(rowSums(!is.na(tax_add)) > 0), ]
  tax_add <- tax_add[!duplicated(tax_add), ]
  tax <- rbind(tax, tax_add)
  rank <- strsplit(rank, '\\/') %>% sapply(nth, -1)
  fwrite(tax, paste0(output_path, rank), sep = '\t')
  
  header <- sapply(1:length(header), function(n) paste(unlist(header_table[n, ]), collapse = ';'))
  header <- gsub('\\;+$', '', header)
  header <- gsub('\\;+', '\\;', header)
  
  fa$Header <- paste(word(fa$Header, 1, end = 1), header, sep = ' ')
  fwrite(data.table(`Seq ID` = word(fa$Header, 1, end = 1), header_table), paste0(output_path, 'header_table.csv'))
  fasta <- strsplit(fasta, '\\/') %>% sapply(nth, -1)
  write_fasta(fa, paste0(output_path, fasta))
}

library(optparse)
option_list <- list(
  make_option(opt_str = c('--fasta'), type = 'character', default = F), 
  make_option(opt_str = c('--rank'), type = 'character', default = F), 
  make_option(opt_str = c('--ref_object'), type = 'character', default = F), 
  make_option(opt_str = c('--nworker'), type = 'integer', default = F), 
  make_option(opt_str = c('--output_path'), type = 'character', default = F)
)
args <- parse_args(OptionParser(option_list=option_list, usage = 'fa mod'))

fa_mod(fasta = args$fasta, rank = args$rank, ref_object = args$ref_object, nworker = args$nworker, output_path = args$output_path)
