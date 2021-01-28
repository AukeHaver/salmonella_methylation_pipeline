### AUKE HAVER
### 10-12-2020
### Scripts for PacBio basemods/motif analysis


# Function to read GFF files without getting warnings
readgff<-function(gff_file){
  suppressWarnings(gff_file <- readGFF(gff_file) %>% as_tibble())
  return(gff_file)
}

# Function for cleaning gff file 
cleangff <- function(dataframe){
  return(
    dataframe %>% 
      rowwise()%>%
      # obtain protein id for specific protein sequence
      mutate(protein_id = as.character(inference),
             protein_id = gsub(".*.:","",protein_id),
             protein_id = paste0(">",gsub("[\")]","",protein_id)),
             # clean prokka name format
             Name = gsub("XXX","",Name),
             # clean prokka gene format
             gene = gsub("XXX","",gene),
             product = ifelse(str_length(product)<3,NA,product),
             # Extract contig number
             contig = gsub(".*_","",seqid)%>% as.integer())%>%
      select(-c(source,
                score,
                phase,
                ID,
                locus_tag,
                inference,
                note,
                rpt_family,
                rpt_type,
                rpt_unit_seq))
  )
}
# Function below appears highly similar to function above, perhaps remove one or the other.
# Function to clean gene annotations
clean_gene_annotation <- function(gene_dataframe){
  return(
    gene_dataframe %>%
      rowwise()%>%
      mutate(protein_id = as.character(inference),
             protein_id = gsub(".*.:","",protein_id),
             protein_id = paste0(">",gsub("[\")]","",protein_id)),
             product = ifelse(str_length(product)<3,NA,product),
             contig = gsub(".*_","",seqid)%>% as.integer()) %>%
      select(-c(source,
                score,
                phase,
                ID,
                locus_tag,
                inference,
                note,
                rpt_family,
                rpt_type,
                rpt_unit_seq))
  )
}

# Function to determine the annotation type of a CSD
determine_annotation_category <- function(annotation_type,annotation_product,annotation_source){
  if(annotation_source == "UniProt_VF" | annotation_source=="VFDB"){
    return("Virulence Factor")}
  if(annotation_source == "Uniprot_Achtmann"|annotation_source =="cgMLST"){
    return("MLST")
  }
  if(annotation_type=="CDS" & is.na(annotation_product)){
    return(annotation_source)}
  else if(annotation_type=="CDS" & annotation_product=="hypothetical protein"){
    return("hypothetical protein")}
  else if (annotation_type=="CDS"){
    return(annotation_source)}
  else {
    return(annotation_type)
  }}

# Function to filter motif dataframe for basemods linked to motifs and annotated to a gene or locus. 
# If a basemod is linked to multiple motifs, that basemod is split into multiple rows.
motif_df_filter <- function(motif_dataframe){
  return(
    motif_dataframe %>%
      filter(!is.na(motif)) %>%
      filter(motif!="character(0)") %>%
      filter(Name !="unmapped") %>%
      filter(Name !="NA") %>%
      separate_rows(motif,sep=",") %>% 
      mutate(motif = gsub("[\"()c ]","",motif),
             ID = row_number()) #removes formatting from MotifMaker (regex)
  )
}

# Function to return the reverse complement of a dna sequence
reverse_complement <- function(dna_string){
  ans <- ""
  base_list <- str_split(dna_string,"")[[1]]
  for(i in 1:length(base_list)){
    ans <- paste0(if(base_list[i]=="A"){"T"}
                  else if(base_list[i]=="T"){"A"}
                  else if(base_list[i]=="G"){"C"}
                  else if(base_list[i]=="C"){"G"}
                  else if(base_list[i]=="Y"){"R"}#added 17-12-2020 for motif mtase match
                  else if(base_list[i]=="R"){"Y"}#added 17-12-2020 for motif mtase match
                  else {"N"},ans)
  }
  return(ans)
}
# Function to check for any two bases wheter they are complementary (according to IUPAC nucleiotide codes)
is_complementary_base <- function(base1,base2){
  return(
    if(base1=="N"& base2=="N"){TRUE}
  else if(base1=="C" & base2=="G"){TRUE}
  else if(base1=="C" & base2=="R"){TRUE} 
  else if(base1=="G" & base2=="C"){TRUE}
  else if(base1=="G" & base2=="Y"){TRUE}
  else if(base1=="A" & base2=="T"){TRUE}
  else if(base1=="A" & base2=="Y"){TRUE}
  else if(base1=="T" & base2=="A"){TRUE}
  else if(base1=="T" & base2=="R"){TRUE}
  else if(base1=="R" & base2=="Y"){TRUE}
  else if(base1=="R" & base2=="C"){TRUE}
  else if(base1=="R" & base2=="T"){TRUE}
  else if(base1=="Y" & base2=="R"){TRUE}
  else if(base1=="Y" & base2=="A"){TRUE}
  else if(base1=="Y" & base2=="G"){TRUE}
  else{FALSE})
}

# Function to check if a sequence is the reverse complement of another
check_reverse_complement <- function(sequence1,sequence2){
  # Split the first sequence into a list of characters
  sequence1 <- str_split(sequence1,"",simplify=TRUE)
  # Reverse the second sequence and split it into a list of characters
  sequence2 <- reverse(sequence2) %>% str_split("",simplify=TRUE)
  # Count the number of matches between the sequences
  matches <- list()
  for(i in 1:length(sequence1)){
    matches[i] <- is_complementary_base(sequence1[i],sequence2[i])}
  return(length(sequence1)==sum(unlist(matches)))
  }


# Function to identify gene matches
is_match <- function(query_index,match_index,motif_list,ID_list,strand_list){
  # Check if basemods are on opposite strands
      query_strand <- strand_list[query_index]
      match_strand <- strand_list[match_index]
      return(
      if(query_strand!=match_strand){
        query_motif <- motif_list[query_index]
        match_motif <- motif_list[match_index]
        if(str_length(query_motif)!=str_length(match_motif)){FALSE}
        else if(query_motif==match_motif | check_reverse_complement(query_motif,match_motif))
        {TRUE} 
        else{FALSE}}
      else{FALSE}
  )}

# Function to return a motif. If one string is NA, pick the other. If both are of equal size, sort them alphabetically and pick the first.
decide_motif <- function(string_1,string_2){
  ans <- 
    if(is.na(string_1))
    {string_2}
  else if (is.na(string_2))
  {string_1}
  else
  {string_1}
  return(sort(c(string_1,string_2))[1])
}


# Function to identify basemods in each occurence of a motif
match_motifs <- function(input_df){
  return_list <- list()
  # Filter the dataframe
  input_df<- input_df %>% 
    motif_df_filter()
  # Make a list of all basemod positions, of all IDs, of all motifs, and of all strand attributes
  bm_position_list <- input_df$start
  bm_ID_list <- input_df$ID
  bm_motif_list <- input_df$motif
  bm_strand_list <- input_df$sense_antisense
  # Initialize lists for all sense and antisense basemods belonging to 1 instance of a motif
  sense_list = vector(mode="list")
  antisense_list = vector(mode="list")
  # Loop over all basemods to identify matching basemods based solely on location
  for(i in 1:nrow(input_df)){
    query_id <- input_df$ID[i]
    query_pos <- input_df$start[i]
    motif_size <- bm_motif_list[query_id] %>% str_length()
    match_list <- c(query_id,input_df$ID[(input_df$start > query_pos - motif_size & input_df$start < query_pos + motif_size)]) %>%
      unique()
    # Filter matches for (1) location on opposite strands and identical motif
    filtered_match_list <- match_list[1]
    if(length(match_list)>1){
      for(j in 2:length(match_list)){
        if(is_match(match_list[1],match_list[j],bm_motif_list,bm_ID_list,bm_strand_list)){
          filtered_match_list <- c(filtered_match_list,match_list[j])
        }
      }
    }
    # If the query has 1 match, add it and the match to the appropriate strand list
    if(length(filtered_match_list)==2){
      match_id <- filtered_match_list[2]
      if(bm_strand_list[query_id]=="sense"){
        sense_list <- c(sense_list,query_id)
        antisense_list <- c(antisense_list,match_id)
      }
      else {
        sense_list <- c(sense_list,match_id)
        antisense_list <- c(antisense_list,query_id)
      }
    }
    # Else, if the basemod has no match, add it to the correct list and add NA to the other
    else if (length(filtered_match_list)==1){
      if(bm_strand_list[query_id]=="sense"){
        sense_list <- c(sense_list,query_id)
        antisense_list <- c(antisense_list,NA)
      }
    # Else do the opposite
    else {
      sense_list <- c(sense_list,NA)
      antisense_list <- c(antisense_list,query_id)}}}
  result_df <- tibble(s_index = unlist(sense_list), a_index=unlist(antisense_list)) %>% 
    distinct() %>% 
    rowwise %>% 
    mutate(seqid = input_df$seqid[ifelse(is.na(s_index),a_index,s_index)],
           # The antisense strand is read into mRNA.
           # therefore, if the motifs slightly differ (Y or R instead of C,A,T,G), the motif call on the antisense strand is leading.
           motif = decide_motif(input_df$motif[a_index][1],input_df$motif[s_index][1]),
           s_basemod = input_df$start[s_index],
           a_basemod = input_df$start[a_index],
           s_type = input_df$type[s_index],
           a_type = input_df$type[a_index],
           s_IPD = input_df$IPDRatio[s_index],
           a_IPD = input_df$IPDRatio[a_index],
           Name = input_df$Name[ifelse(is.na(s_index),a_index,s_index)],
           hemi_methylation = ifelse(is.na(s_basemod)|is.na(a_basemod),TRUE,FALSE))%>% 
    select(-c("a_index","s_index")) %>% 
    group_by(s_basemod) %>% 
    mutate(s_duplicate=n()) %>%
    ungroup %>% 
    group_by(a_basemod) %>% 
    mutate(a_duplicate=n()) %>% ungroup() %>%
    filter(!(hemi_methylation & s_duplicate>1 & a_duplicate>1)) %>%
    select(-c("a_duplicate","s_duplicate"))
  return(result_df)}

# Function to map basemods back to the location on the gene.
map_on_gene <- function(basemod,gene,annotation_df){
  # if either the basemod or gene is NA, return NA
  if(is.na(basemod)|is.na(gene)){return(NA)}
  # Else, subtract the start position of the gene from the basemod position, and add 1 so the 1st position is not index 0
  else {
    return(basemod - annotation_df$start[annotation_df$Name == gene] +1)
  }
}

# Function for reading annotation file into dataframe
read_annot_df <- function(reference_file,reference_folder){
  return(
    readgff(paste0(reference_folder,"/",reference_file)) %>% 
      mutate(Name = gsub("XXX","",Name),
             gene = gsub("XXX","",gene),
             gene = gsub("_.$","",gene)) %>%
      filter(!is.na(Name)))}

# Function for identifying correctly annotated motifs.
# Both one of the present values needs to be positive
correctly_annotated <- function(basemod1,basemod2,gene,annotation_df){
  gene_start <- annotation_df$start[annotation_df$Name==gene]
  gene_end <- annotation_df$end[annotation_df$Name==gene]
  basemod1_outside <- ifelse(is.na(basemod1),TRUE,
                             if(basemod1<gene_start|basemod1>gene_end){TRUE}
                             else{FALSE})
  basemod2_outside <- ifelse(is.na(basemod1),TRUE,
                             if(basemod1<gene_start|basemod1>gene_end){TRUE}
                             else{FALSE})
  return(ifelse(basemod1_outside&basemod2_outside,FALSE,TRUE))
}
 

# Map to Gene
# The function takes 3 inputs: (1) the contig of the basemod, (2) the location of that basemod on the specified contig, (3) a dataframe containing the annotations for that genome.
map_to_gene <- function(contigs,basemod,annotation_df){
  # First filter the annotation dataframe for genes on the specified contig
  annotation_df <- annotation_df %>% filter(contig == contigs)
  # Make a list of the indices of the final base of each gene
  end_list <- annotation_df$end
  # Make a list of which strand those genes are located on
  strand_list <- annotation_df$strand
  # Make a list of the indicies of the first base of each gene
  start_list <- annotation_df$start
  # Make a list of the names of those genes
  Name_list <- annotation_df$Name
  # Select the first occurance of: the index for the final base of a gene which is located further in the strand than the basemod.
  index <- which(basemod-20<=end_list)[1]
  # Return statement
  return(
    # if the location of the basemod is not compatible than, or higher than any final base index of any gene, return "unmapped"
    if(is.na(index)){"unmapped"} 
    else{
      # if the location of the basemod is not located between the START and END of the first gene with an END higher than or equal to the index of the basemod, then surely the basemod is not located on a gene: "unmapped"
      if(!start_list[index]<=basemod+20){"unmapped"}
      # return: the Name of the gene belonging to the gene we found the basemod on and the strand on which the gene is situated (for sense/antisense comparison, separated by a tab)
      else{paste0(Name_list[index],"\t",strand_list[index])}
    }
  )
  
}

# Annotate basemods
# The function takes two inputs, the sample name and the corresponding reference (split, as pangenomes use the same reference).
annotate_basemods <- function(sample_name,reference){
  # First, make a dataframe from the Prokka annotation .gff file.
  annot_df <- 
    # The file is located in the REF directory
    readgff(paste0("../REF/",reference))%>%
    # A leading and trailing "XXX" was added before annotation with Prokka to preserve Locustag names.
    mutate(Name = gsub("XXX","",Name),
           # The sequence ID contains the contig identifier.
           # It is vital that any matching gene is checked for being located on the same contig. 
           # This is saved in an additional column.
           contig = gsub("^.*_","",seqid)) %>%
    filter(!is.na(Name))%>%
    # As little information as possible is stored in the environment.
    select(contig,start,end,Name,strand)
  # Then, make a dataframe from the motifmaker basemods and match basemods to genes. 
  # This dataframe is not stored in the environment but immediately exported to a .gff file.
  readgff(paste0("../smrtlink/motifs/motifs_",sample_name,".gff")) %>%
    # Again the contig number is extracted from the sequence ID.
    mutate(contig = gsub("^.*_","",seqid)) %>% 
    # All gene mapping needs to be performed in a rowwise fashion. 
    rowwise() %>%
    # Execute the function in the chunk above.
    mutate(Name = map_to_gene(contig,start,annot_df)) %>% 
    # The function in the chunk above returns the name and strand of the corresponding gene in one character string.
    # This needs to be split.
    separate(Name,into=c("Name","sense_antisense"),sep="\t",fill="right") %>%
    mutate(sense_antisense = ifelse(sense_antisense==strand,"sense","antisense")) %>%
    select(-contig) %>%   
    export(paste0("../smrtlink/annotated_motifs/motifs_",sample_name,".gff"),format="GFF3")
  return(paste0("finished ", sample_name))
}

# Function to filter motifs previously incorrectly atributed as being located on a specific gene
filter_motifs <- function(motif_df, annotation_df){
  # list the genes, start and end positions from the annotation dataframe.
  return(
    motif_df %>% 
      rowwise()%>%
      mutate(gene_start = annotation_df$start[annotation_df$Name==Name],
             gene_end = annotation_df$end[annotation_df$Name==Name],
             # In case of presumed hemi methylation, replace NA values with the location of the remaining basemod.
             basemod1 = ifelse(is.na(a_basemod),s_basemod,a_basemod),
             basemod2 = ifelse(is.na(s_basemod),a_basemod,s_basemod),
             # For each basemod in a motif, return whether that basemod is located on the gene
             basemod1_in_gene = ifelse(basemod1>=gene_start & basemod1 <=gene_end,
                                       TRUE,
                                       FALSE),
             basemod2_in_gene = ifelse(basemod2>=gene_start & basemod2 <=gene_end,
                                       TRUE,
                                       FALSE),
             # If both basemods of a motif are not located on the attributed gene, the annotation is incorrect.
             correct_annotation = basemod1_in_gene | basemod2_in_gene,
             gene_length=gene_end-gene_start+1) %>% 
      ungroup()%>%
      # Filter for correct annotations
      filter(correct_annotation) %>%
      # Select only desired columns
      select(-c(basemod1,
                basemod2,
                basemod1_in_gene,
                basemod2_in_gene,
                correct_annotation,
                gene_end)))
}

# Function to change the location of a basemod into the location of the basemod on its associated gene
return_gene_index <- function(motif_df){
  return(
    motif_df %>% 
      rowwise() %>%
      mutate(s_basemod = ifelse(is.na(s_basemod),NA,s_basemod-gene_start),
             a_basemod = ifelse(is.na(a_basemod),NA,a_basemod-gene_start))
  )
}

make_motif_matches<-function(info_file,destination_folder){
  for(i in 1:nrow(info_file)){
    annot <- read_annot_df(info_file$Reference[i],"../REF/")
    overview_df <- readGFF(paste0("../smrtlink/annotated_motifs/motifs_",info_file$Sample[i],".gff")) %>%
      match_motifs %>%
      filter_motifs(annot) %>% 
      return_gene_index() %>%
      mutate(contig = as.integer(gsub(".*_","",seqid)),
             sample = info_file$Sample[i],
             serovar= info_file$Serovar[i]) %>%
      write_tsv(paste0(destination_folder,"/matched_motifs_",info_file$Sample[i],".tsv"))}}

### Function to read annotated CDS into dataframe and combine with protein database 
# If the sample is a pangenome, only change sample name to "pangenome"
read_annotation_overview <-function(design_file,annotation_directory,protein_database){
  for(i in 1:nrow(design_file)){
    if(i ==1){
      annotated_genes_overview <- 
        paste0(annotation_directory,design_file$Reference[i]) %>% 
        readgff() %>% 
        # If the sample is a pangenome, only change sample name to "pangenome"
        mutate(sample = ifelse(design_file$Sequence_type[i]=="pangenome",
                               "pangenome", 
                               design_file$Sample[i])) %>%
        cleangff()
    } else {
      annotated_genes_overview <- 
        paste0("../REF/",design_file$Reference[i]) %>% 
        readgff() %>%
        mutate(sample = ifelse(design_file$Sequence_type[i]=="pangenome",
                               "pangenome", 
                               design_file$Sample[i])) %>%
        cleangff() %>%
        full_join(annotated_genes_overview,
                  by=colnames(annotated_genes_overview))}
  }
  # Add metadata info from non_redundant_protein_db to the dataframe of annotated genes
  annotated_genes_overview <- protein_database %>%
    rowwise()%>%
    mutate(protein_id = SeqID,
           product = as.character(metadata),
           Gene = as.character(Gene),
           Locus = as.character(Locus),
           Source=as.character(Source)) %>%
    select(-c(sequence,SeqID,duplicates,metadata)) %>%
    right_join(annotated_genes_overview,by="protein_id") %>%
    rowwise() %>%
    mutate(Source=ifelse(is.na(Source),"UniProtKB",Source),
           gene = ifelse(is.na(Gene),gene,Gene),
           product = ifelse(is.na(product.x),product.y,product.x),
           sequence_length=end-start+1)%>%
    select(-c(seqid,Gene,product.x,product.y)) %>%
    relocate(sample,contig,type,start,end,strand,gene,Locus,Name,Source,sequence_length,product) %>%
    distinct()
  return(annotated_genes_overview)
  
}

# Function to write motif is short form
short_motif_form <- function(motif_string){
  return(
    gsub(pattern="N+",
         replacement=paste0("N",str_count(motif_string,"N")),
         motif_string))
}
### EXTRA
# Function to read mapping counts 
read_mapping_counts <- function(index){
  return(readLines(paste0("../output/mapping_counts/sample_",samples[index],".txt"))%>% as.integer())
}

# Function to read a motifs file
read_motifs_file <- function(directory,sample_name){
  return(readGFF(
    paste0(directory,"motifs_",sample_name,".gff")) %>%
      as_tibble()) %>%
    rowwise()%>%
    mutate(in_motif = is.na(motif) == FALSE & motif !="character(0)",
           Sample = sample_name) %>%
    select(-motif,-id)
}
