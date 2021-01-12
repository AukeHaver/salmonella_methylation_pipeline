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