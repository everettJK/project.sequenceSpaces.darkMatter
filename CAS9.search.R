library(dplyr)
library(Biostrings)
options(stringsAsFactors = FALSE)

mafftBin    <- '~/software/mafft-7.397/bin'
hmmerBin    <- '~/software/hmmer-3.1b2/bin'

system(paste0(hmmerBin, "/hmmsearch --cpu 25 -Z 1000000  --domZ 1000000 --tblout data/CAS9_swprot.hmmsearch ",
              "data/hmms/CAS9_swprot.hmm  data/all_orf_unannotated_aaseq.fa"), ignore.stdout = TRUE, ignore.stderr = TRUE)


CAS9_swprot.ids <- scan('data/CAS9_swprot.ids', what = 'character', sep = '\n')

RunC_seqs <- bind_rows(lapply(CAS9_swprot.ids, function(x){
  p <- readLines(url(paste0('https://www.uniprot.org/uniprot/', x, '.txt')))
  if(length(grep('REGION\\s+\\d+\\s+\\d+\\s+RuvC\\-', p) == 3)){
    sequence <- readLines(url(paste0('https://www.uniprot.org/uniprot/', x, '.fasta')))
    sequence <- paste0(sequence[2:length(sequence)], collapse = '')
    o <- stringr::str_match_all(p[grep('REGION\\s+\\d+\\s+\\d+\\s+RuvC\\-', p)], 'REGION\\s+(\\d+)\\s+(\\d+)\\s+RuvC\\-')
    if(is.unsorted(as.integer(unlist(lapply(o, '[[', 2))))) browser() #stop('The RuvC records are not ordered.')
    o <- unlist(lapply(o, function(x2){ base::substr(sequence, as.numeric(as.vector(x2)[2]), as.numeric(as.vector(x2)[3])) }))
    return(data.frame(id = x, RuvC1 = o[1], RuvC2 = o[2], RuvC3 = o[3]))
  } else{
    return(data.frame())
  }
}))

HNH_seqs <- bind_rows(lapply(CAS9_swprot.ids, function(x){
  p <- readLines(url(paste0('https://www.uniprot.org/uniprot/', x, '.txt')))
  if(length(grep('DOMAIN\\s+\\d+\\s+\\d+\\s+HNH', p) == 1)){
    sequence <- readLines(url(paste0('https://www.uniprot.org/uniprot/', x, '.fasta')))
    sequence <- paste0(sequence[2:length(sequence)], collapse = '')
    o <- stringr::str_match_all(p[grep('DOMAIN\\s+\\d+\\s+\\d+\\s+HNH', p)], 'DOMAIN\\s+(\\d+)\\s+(\\d+)\\s+HNH')
    return(data.frame(id = x, HNH = unlist(lapply(o, function(x2){ base::substr(sequence, as.numeric(as.vector(x2)[2]), as.numeric(as.vector(x2)[3])) }))))
  } else{
    return(data.frame())
  }
}))

createHMM <- function(name, x, n, dir = 'data'){
  writeXStringSet(AAStringSet(setNames(x, n)), filepath = file.path(dir, paste0(name, '.ff')))
  system(paste0(mafftBin, '/mafft --maxiterate 5000 --globalpair ', file.path(dir, paste0(name, '.ff')), ' > ', file.path(dir, paste0(name, '.mafft'))))
  system(paste0(hmmerBin, '/hmmbuild -n ', paste0(name, '.hmm'), ' --amino ', file.path(dir, paste0(name, '.hmm')), ' ', file.path(dir, paste0(name, '.mafft'))))
}

createHMM('HNH', HNH_seqs$HNH, HNH_seqs$id)
createHMM('RuvC1', RunC_seqs$RuvC1, RunC_seqs$id)
createHMM('RuvC2', RunC_seqs$RuvC2, RunC_seqs$id)
createHMM('RuvC3', RunC_seqs$RuvC3, RunC_seqs$id)
          
system(paste0(hmmerBin, "/hmmsearch --cpu 25 -Z 1000000  --domZ 1000000 --tblout data/CAS9_swprot.hmmsearch data/CAS9_swprot.hmm data/all_orf_unannotated_aaseq.fa"))
system(paste0(hmmerBin, "/hmmsearch --cpu 25 -Z 1000000  --domZ 1000000 --tblout data/HNH.hmmsearch data/HNH.hmm data/all_orf_unannotated_aaseq.fa"))
system(paste0(hmmerBin, "/hmmsearch --cpu 25 -Z 1000000  --domZ 1000000 --tblout data/RuvC1.hmmsearch data/RuvC1.hmm data/all_orf_unannotated_aaseq.fa"))
system(paste0(hmmerBin, "/hmmsearch --cpu 25 -Z 1000000  --domZ 1000000 --tblout data/RuvC2.hmmsearch data/RuvC2.hmm data/all_orf_unannotated_aaseq.fa"))
system(paste0(hmmerBin, "/hmmsearch --cpu 25 -Z 1000000  --domZ 1000000 --tblout data/RuvC3.hmmsearch data/RuvC3.hmm data/all_orf_unannotated_aaseq.fa"))


parse_hmmsearch_output <- function(file){
   h <- scan(file, what = 'character', sep = '\n')
   h <- dplyr::bind_rows(lapply(h[-grep('^#', h)], function(x){
          x <- unlist(stringr::str_split(x, '\\s+'))
          data.frame(targetName          = x[1],
                     fullSeq.eval        = as.numeric(x[5]), 
                     fullSeq.bitScore    = as.numeric(x[6]), 
                     fullSeq.bias        = as.numeric(x[7]), 
                     bestDomain.eval     = as.numeric(x[8]), 
                     bestDomain.bitScore = as.numeric(x[9]), 
                     bestDomain.bias     = as.numeric(x[10]), 
                     desc                = paste(x[19:length(x)], collapse = ' '))
  })) 
  dplyr::arrange(h, fullSeq.eval)
}
   
h <- parse_hmmsearch_output('data/CAS9_swprot.hmmsearch') 
HNH <- parse_hmmsearch_output('data/HNH.hmmsearch') 
RuvC1 <- parse_hmmsearch_output('data/RuvC1.hmmsearch') 
RuvC2 <- parse_hmmsearch_output('data/RuvC2.hmmsearch') 
RuvC3 <- parse_hmmsearch_output('data/RuvC3.hmmsearch') 

bind_rows(lapply(1:nrow(h), function(x){
  d <- h[x,]
  d <- select(d, -fullSeq.bitScore, -desc, -fullSeq.bias, -bestDomain.eval, -bestDomain.bitScore, -bestDomain.bias)
  d$HNH.eval <- RuvC1[match(d$targetName, HNH$targetName),]$fullSeq.eval
  d$RuvC1.eval <- RuvC1[match(d$targetName, RuvC1$targetName),]$fullSeq.eval
  d$RuvC2.eval <- RuvC1[match(d$targetName, RuvC2$targetName),]$fullSeq.eval
  d$RuvC3.eval <- RuvC1[match(d$targetName, RuvC3$targetName),]$fullSeq.eval
  d
})) %>%
  filter(fullSeq.eval <= 1e-10)
