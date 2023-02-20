library(MicroEcoTk)
library(tidyfst)
library(stringr)
library(doParallel)

# output_path <- 'C:/Users/KaiMa/Desktop/xx'
output_path <- '/clusterfs/node1/makai/Database/Unite2/modFun/'
fa <- read_fasta('/clusterfs/node1/makai/Database/Unite2/UNITE_public_29.11.2022.fasta')
# fa <- read_fasta('G:/UNITE_public_29.11.2022.fasta')
nworker <- 62

fa_fungi <- fa[which(grepl('|k__Fungi;', fa$Header)), ]

fa_fungi <- fa_fungi[!duplicated(fa_fungi$Sequence), ]
header <- fa_fungi$Header
lineage <- strsplit(header, '\\|') %>% sapply('[[', 2)

ranks <- c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
w <- c('k', 'p', 'c', 'o', 'f', 'g', 's')
# header_table <- matrix('', nrow(fa_fungi), length(ranks))

# tax <- matrix(NA, nrow(fa_fungi)*length(ranks), 5)

# tax_id <- 1
famod_item <- function(seq) {
  xx <- rep(12344321, length(ranks))
  s <- unlist(strsplit(lineage[seq], ';'))
  tax <- matrix(NA, length(s), 5)
  for (n in 1:length(s)) {
    # n = 1
    if (strsplit(s[n], '__') %>% sapply('[', 1)  == w[n]) {
      xx[n] <- paste(w[n], (strsplit(s, '__') %>% sapply('[', 2))[n], sep = '__')
      tax[n, ] <- c(sapply(strsplit(s[n], '__'), '[', 2), xx[n], ranks[n], match(ranks[n], ranks), lineage[seq])
    }
  }
  list(xx, as.data.table(tax))
}

# for (seq in 1:nrow(fa_fungi)) {
#   # seq = 1
#   
#   for (n in 1:length(s)) {
#     # n = 1
#     if (strsplit(s[n], '__') %>% sapply('[', 1)  == w[n]) {
#       header_table[seq, n] <- paste(ifelse(n == 1, 'd', w[n]), (strsplit(s, '__') %>% sapply('[', 2))[n], sep = '__')
#       if(all(!grepl(paste0('^', header_table[seq, c(get('n'))], '$'), tax[, 2]))) {
#         tax[tax_id, ] <- c(sapply(strsplit(s[n], '__'), '[', 2), header_table[seq, c(get('n'))], ranks[n], match(ranks[n], ranks), lineage[seq])
#         tax_id <- tax_id + 1
#       }
#     }
#   }
#   cat('sequence', paste0(seq, '/', nrow(fa_fungi)), '----', header_table[seq, ], '\n')
# }
cl <- makeCluster(nworker, type = "PSOCK") # PSOCK (windows); FORK (linux);
registerDoParallel(cl)
header_out <- foreach(seq = 1:nrow(fa_fungi), .packages = c('stringr', 'tidyfst')) %dopar% famod_item(seq)

header_table <- t(sapply(1:length(header_out), function(i) {
  header_out[[i]][[1]]
}, simplify = T))
header_table[header_table == 12344321] <- ''
colnames(header_table) <- ranks
tax <- lapply(1:length(header_out), function(i) {
  header_out[[i]][[2]]
}) %>% rbindlist()
tax <- tax[!duplicated(tax[, -5]), ]
colnames(tax) <- c('rank', 'class_rank', 'class', 'depth', 'tax')

tax <- as.data.frame(tax)
tax <- tax[which(rowSums(!is.na(tax)) > 0), ]
cat('\n')
if(!dir.exists(output_path)) dir.create(output_path)

tax_add <- matrix(NA, nrow(header_table)*ncol(header_table), 5)
tax_id <- 1
for (i in 2:ncol(header_table)) {
  # i=2
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
    sons_id <- which(yes_no)
    fathers <- header_table[, c(get('i'))-1]
    table(fathers)
    for (n in sons_id) {
      # n = 69
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
}
colnames(tax_add) <- c('rank', 'class_rank', 'class', 'depth', 'tax')
tax_add <- tax_add[which(rowSums(!is.na(tax_add)) > 0), ]
tax_add <- tax_add[!duplicated(tax_add), ]
tax <- rbind(tax, as.data.table(tax_add))
fwrite(tax, paste0(output_path, '/tax_UNITE_public_29.11.2022.txt'), sep = '\t')

header <- sapply(1:length(header), function(n) paste(unlist(header_table[n, ]), collapse = ';'))

fa_fungi$Header <- paste(word(fa_fungi$Header, 1, end = 1, sep = '\\|'), header, sep = ' ')
fwrite(data.table(`Seq ID` = word(fa_fungi$Header, 1, end = 1), header_table), paste0(output_path, '/header_table.csv'))
write_fasta(fa_fungi, paste0(output_path, '/mod_UNITE_public_29.11.2022.fa'))
