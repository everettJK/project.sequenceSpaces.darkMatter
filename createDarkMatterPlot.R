library(Biostrings)
library(tidyverse)

a <- readAAStringSet('data/all_orf_unannotated_aaseq.fa')
seqAttributes <- data.frame(id = names(a), length = width(a))

blast <- read.table('data/all_orf_unannotated_aaseq.fa.chunk-1.blastp', sep = '\t', header = FALSE)

system('scp microb120:/home/everett/dm/all_orf_unannotated_aaseq.blastp data/swprot.blastp')
blast <- read.table('data/swprot.blastp', sep = '\t', header = FALSE)


# # Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
names(blast) <- c('qseqid','sseqid','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore')
blast$qAlignLen <- blast$qend - blast$qstart + 1
blast$qLength <- seqAttributes[match(blast$qseqid, seqAttributes$id),]$length
blast$percentCoverage <- (blast$qAlignLen / blast$qLength)*100
blast <- select(blast, qseqid, pident, percentCoverage)

percentSeqsFound <- round((n_distinct(blast$qseqid) / n_distinct(seqAttributes$id))*100, 1)

ppNum <- function(n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)

ggplot(blast, aes(percentCoverage, pident)) +
  theme_bw() +
  geom_point(alpha = 0.1) +
  ylim(c(0,100)) +
  scale_x_continuous(limits = c(0, 100)) +
  labs(x = '% sequence coverage', y = '% identity') +
  geom_hline(yintercept=50, color = 'blue') +
  geom_vline(xintercept=50, color = 'blue') +
  ggtitle(paste0('Dark matter alignments to SwissProt - ', ppNum(n_distinct(blast$qseqid)), ' ORFs'))


