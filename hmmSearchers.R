library(dplyr)

# system('~everett/software/hmmer-3.1b2/bin/hmmsearch --tblout query.vfam.tbl --domtblout query.vfam.dom.tbl vFam-B_2014.hmm query.ff > /dev/null')
# system('~everett/software/hmmer-3.1b2/bin/hmmsearch --tblout query.tbl --domtblout query.dom.tbl Pfam-A.hmm query.ff > /dev/null')

parseHMMSearch <- function(file){
  r <- scan(file, what = 'character', sep = '\n')
  r <- r[-which(grepl('^#', r))]
  colNames <- c('target_name', 'accession', 'tlen', 'query_name',  'accession',   'qlen',   'full.E-value',  'full.score', 'full.bias',   
                'domain.num',  'domain.total', 'domain.c-Evalue',  'domain.i-Evalue',  'domain.score',  'domain.bias',  'domain.hmm.from',
                'domain.hhm.to',  'domain.ali.from',    'domain.ali.to',  'domain.env.from',    'domain.env.to')
    
  bind_rows(lapply(base::strsplit(r, '\\s+'), function(x){
                     x <- x[1:length(colNames)]
                     names(x) <- colNames
                     data.frame(t(x))
                   }))
}

v <- parseHMMSearch('query.vfam.dom.tbl') %>%
     group_by(target_name) %>%
     mutate(full.score = as.numeric(full.score)) %>%
     mutate(domain.score = as.numeric(domain.score)) %>%
     summarise(score = ifelse(max(full.score) > max(domain.score), max(full.score), max(domain.score))) %>%
     ungroup()

p <- parseHMMSearch('query.dom.tbl') %>%
     group_by(target_name) %>%
     mutate(full.score = as.numeric(full.score)) %>%
     mutate(domain.score = as.numeric(domain.score)) %>%
     summarise(score = ifelse(max(full.score) > max(domain.score), max(full.score), max(domain.score))) %>%
     ungroup()

v
p

    