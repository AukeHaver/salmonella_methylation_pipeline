### Function for reading .fasta files
ReadFasta <- function(file){
  fasta <- readLines(file)
  index <- grep(">",fasta)
  s <- data.frame(index=index,from=index+1, to=c((index-1)[-1],length(fasta)))
  seqs <- rep(NA,length(index))
  for(i in 1:length(index)){
    seqs[i]<-paste(fasta[s$from[i]:s$to[i]],collapse="")
  }
  return(data.frame(name=gsub(">","",fasta[index]),sequence=seqs))
}

### Function for downloading Uniprot fasta files
Read_Uniprot_Fasta <- function(uniprot_id){
  fasta <- readLines(paste0("https://www.uniprot.org/uniprot/",uniprot_id,".fasta"))
  description = fasta[1]
  sequence = ""
  for(i in 2:length(fasta)){
    sequence <- paste0(sequence,fasta[i])
  }
  return(c(description, sequence))
}

## Function for calculating GC ratios
calculate_gc_ratio <- function(sequence){
  cg_count <- sum(str_count(sequence,c("C","G")))
  sequence_length <- str_length(sequence)
  ratio <- cg_count/sequence_length*100
  return(ratio %>%
           round(1) %>%
           format(nsmall=1) %>%
           as.character())
}

## Function for generating genome descriptions based on contigs and GC ratios
generate_genome_description <- function(directory,genome_name){
  return(
    paste0(directory,genome_name) %>%
      ReadFasta() %>%
      rowwise() %>%
      mutate(genome = genome_name,
             contig_length = str_length(sequence),
             GC_ratio = calculate_gc_ratio(sequence)) %>%
      ungroup()%>%
      arrange(desc(contig_length))%>%
      mutate(contig = row_number())%>%
      select(genome,contig,contig_length,GC_ratio)
  )
}