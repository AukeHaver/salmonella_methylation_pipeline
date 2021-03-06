### Function to make a modified design file based on a setup file
# Currently only available for pangeome of all included samples (no filtering)
# This can later be expanded for addition of multiple pagenomes
make_analysis_file <- function(design_dataframe,setup_dataframe){
  # Read design dataframe into another dataframe
  output_design <- design_dataframe
  # Remove References if required
  if(!setup_dataframe$Include_Reference){
    output_design <- output_design%>% filter(!Sequence_type=="reference")
  }
  # Add Hybrid assembly if required
  if(setup_dataframe$Include_Hybrid){
    output_design <-
      # Read design
      design_dataframe %>%
      rowwise() %>%
      # Change varables to correct type
      mutate(Sample = gsub(".*_","hybrid_",Sample),
             Sequence_type = "hybrid",
             Reference = paste0(
               gsub(".*_","hybrid_",Sample),
               ".gff")) %>%
      full_join(output_design,by=colnames(output_design))}
  # Add pangenome assembly if required
  if(setup_dataframe$Include_Pangenome){
    output_design <-
      # Read design
      design_dataframe %>%
      rowwise() %>%
      # Change varables to correct type
      mutate(Sample = gsub(".*_","pangenome_",Sample),
             Sequence_type = "pangenome",
             Reference = paste0(
               gsub(".*_","pangenome_",Sample),
               ".gff")) %>%
      full_join(output_design,by=colnames(output_design))}
  # Add external assembly if required
  if(setup_dataframe$Include_External){
    output_design <-
      # Read design
      design_dataframe %>%
      rowwise() %>%
      # Change varables to correct type
      mutate(Sample = gsub(".*_","external_",Sample),
             Sequence_type = "external",
             Reference = paste0(
               gsub(".*_","external_",Sample),
               ".gff")) %>%
      full_join(output_design,by=colnames(output_design))}
  return(output_design %>%
           mutate(Reference = paste0(Sample,".gff")))}


### Function for generating unicycler assembly commands
generate_unicycler_commands <- function(design_dataframe,directory){
  commands <- list()
  # Parse over all samples
  for(i in 1:nrow(design_dataframe)){
    # copy sample number
    sample_number <- design_dataframe$Sample[i]%>% gsub(pattern="reference_",replacement="")
    # Append information
    commands <- c(commands,paste0(
      "bsub -q bio 'unicycler -1 ../sequencing_data/Illumina_reads/Sample-", 
      sample_number,
      "_R1.fastq.gz ",
      "-2 ../sequencing_data/Illumina_reads/Sample-", 
      sample_number,
      "_R2.fastq.gz ",
      "-l ../sequencing_data/PacBio_reads/Sample-", 
      sample_number,
      ".fastq ",
      "-o hybrid_", 
      sample_number,"'"
    ))
  }
  # Make a list of all lines to write to the script file
  c(
    "### Script for Unicycler assembly of hybrid sequences (Wick et al., 2017)",
    paste0("### Script generated by the Salmonella Methylation Pipeline (SMP) on ",date()),
    "### GitHub url: github.com/AukeHaver/Salmonella_Methylation_Pipeline",
    "",
    "# Activate environment:",
    "conda activate unicycler_env",
    "",
    "# Move to correct working directory",
    paste0("cd ",
           file.path(directory,"unicycler",fsep="/")
    ),
    "",
    "# make output directory",
    "mkdir output",
    "",
    "# Generate unicycler Assemblies",
    commands,
    "",
    "# Deactivate environment",
    "conda deactivate") %>%
    write_lines("../scripts/3_unicycler_bash.txt")
}

### Function for generating pangenome assembly commands
# Can be modified for pangenome generation from other sources
generate_seq_seq_pan_commands <- function(design_dataframe,directory){
  # Create a assembly list file from unicycler assemblies
  paste0("../hybrid/hybrid",
         gsub("reference","",design_dataframe$Sample),
         "/assembly.fasta") %>%
    write_lines("../seq_seq_pan/assembly_list.txt")
  c(# list of commands
    "### Script for pangenomistation of hybrid sequences (Jandrasits et al., 2018)",
    paste0("### Script generated by the Salmonella Methylation Pipeline (SMP) on ",date()),
    "### GitHub url: github.com/AukeHaver/Salmonella_Methylation_Pipeline",
    "",
    "# Activate environment:",
    "conda activate pangenome_env",
    "",
    "# Move to correct working directory",
    paste0("cd ",
           file.path(directory,"seq_seq_pan",fsep="/")
    ),
    "",
    "# Assemble from assembly list",
    "bsub -q bio 'seq-seq-pan-wga --config genomefile=assembly_list.txt outfilename=output/pangenome'",
    "",
    "# Deactivate environment",
    "conda deactivate")%>%
    write_lines("../scripts/4_pangenome_bash.txt")
}

### Function for generating individual annotatino command
generate_prokka_command <- function(sample,sequence_type,reference){
  sample_number <- gsub(".*_","",sample)
  return(
    paste0(
      "bsub -q bio 'prokka --outdir ../prokka/", sample, " \\\n",
      "--force --usegenus --rawproduct --centre X \\\n",
      "--prefix ", sample, " \\\n",
      "--genus Salmonella \\\n",
      "--species enterica \\\n",
      "--proteins ../protein_database/database/senterica_pdb.fasta \\\n",
      "--locustag SE",
      str_split(sequence_type,"")[[1]][1] %>% toupper(),
      sample_number,
      format(Sys.Date(),format= "%d%m%y")," \\\n",
      # In case the sample is a reference genome
      ifelse(sequence_type == "reference", paste0("../REF/",gsub(pattern=".gff",replacement =".fasta",reference)),
             # In case the sample is a hybrid
             ifelse(sequence_type == "hybrid",paste0("../unicycler/",sample,"/assembly.fasta"),
                    # In case the sample is from external
                    ifelse(sequence_type== "external",paste0("../external/",sample,".fasta"),
                           # In case the sample is the pangenome
                           ifelse(sequence_type=="pangenome",paste0("../seq_seq_pan/output/pangenome_consensus.fasta"),"error")))),
      " && \\\n",
      "cp ",sample,"/",sample,".fna ../REF/",sample,".fasta && \\\n",
      "cp ",sample,"/",sample,".gff ../REF/",sample,".gff'\n", ""))}

### Function for generating assembly annotation commands
generate_prokka_commands <- function(design_dataframe,directory){
  for(i in 1:nrow(design_dataframe)){
    if(i == 1){
      commands <- generate_prokka_command(design_dataframe$Sample[i],
                                          design_dataframe$Sequence_type[i],
                                          design_dataframe$Reference[i])
    } else {
      commands <- c(commands,
                    generate_prokka_command(design_dataframe$Sample[i],
                                            design_dataframe$Sequence_type[i],
                                            design_dataframe$Reference[i]))
    }}
  c("### Script for annotation of bacterial assemblies (Seemann. T, 2014)",
    paste0("### Script generated by the Salmonella Methylation Pipeline (SMP) on ",date()),
    "### GitHub url: github.com/AukeHaver/Salmonella_Methylation_Pipeline",
    "",
    "# Activate environment:",
    "conda activate prokka_env",
    "",
    "# Move to correct working directory:",
    paste0("cd ",file.path(directory,"prokka",fsep="/")),
    "",
    "# Annotate Genomes:",
    unique(commands),
    "",
    "# Deactivate environment",
    "conda deactivate") %>%
    write_lines("../scripts/5_prokka_bash.txt")
}

### Function for generating demultiplexing commands
# Currently poorly implemented demultiplex script
generate_smrtlink_commands_part_1 <- function(read_files,directory){
  c("### Script for demultiplexing of reads with Lima",
    paste0("### Script generated by the Salmonella Methylation Pipeline (SMP) on ",date()),
    "### GitHub url: github.com/AukeHaver/Salmonella_Methylation_Pipeline",
    "",
    "# Activate environment:",
    "conda activate smrtlink",
    "",
    "# Move to correct working directory:",
    paste0("cd ",file.path(directory,"sequencing_data",fsep="/")),
    "",
    "# Make neccessary directory:",
    "mkdir demultiplex",
    "",
    "# Demultiplex reads:",
    "bsub -q bio lima \\",
    "same \\",
    "--split-bam-named \\",
    "--min-score 26 \\",
    "--min-length 5000 \\",
    "--num-threads 12  \\",
    paste0("rawdata/",read_files[1]," \\"),
    paste0("rawdata/",read_files[2]," \\"),
    paste0("demultiplex/",gsub(".subreads.bam","demux_strict.subreadset.xml",read_files[1])," \\"),
    "",
    "# Rename output files:",
    "cd demultiplex",
    "for file in * ; do mv \"${file#*--}\"; done",
    "",
    "# Deactivate environment",
    "conda deactivate")%>%
    write_lines("../scripts/6_smrtlink_bash_part_1.txt")
}

#### Function for generating alignment commands
# Sub function for index command
generate_index_command <- function(string){
  return(paste0("bsub -q bio pbmm2 index ../REF/",string,".fasta ../REF/",string,".mmi"))}
# Sub function for align command
generate_align_command <- function(string){
  sample_number <- gsub(".*_","",string)
  sample_type <- gsub(sample_number,"",string)
  return(paste0("bsub -q bio pbmm2 align --sort ../REF/",
                string,
                ".mmi ../sequencing_data/demultiplex/sample_",
                sample_number,".bam alignment/",sample_number,"_",sample_type,"aligned_pbmm2.bam"))}
# Sub function for align command
generate_bpindex_command <- function(string){
  sample_number <- gsub(".*_","",string)
  sample_type <- gsub(sample_number,"",string)
  return(paste0("bsub -q bio pbindex alignment/",sample_number,"_",sample_type,"aligned_pbmm2.bam"))
}
# Main function
generate_smrtlink_commands_part_2 <- function(design_dataframe,directory){
  c("### Script for read alignment with Samtools",
    paste0("### Script generated by the Salmonella Methylation Pipeline (SMP) on ",date()),
    "### GitHub url: github.com/AukeHaver/Salmonella_Methylation_Pipeline",
    "",
    "# Activate environment:",
    "conda activate smrtlink",
    "",
    "# Move to correct working directory:",
    paste0("cd ",file.path(directory,"smrtlink",fsep="/")),
    "",
    "# Make an alignment directory:",
    "mkdir alignment",
    "",
    "# Index reference genomes:",
    lapply(design_dataframe$Sample,FUN=generate_index_command) %>% unlist(),
    "",
    "# Align and sort files:",
    lapply(design_dataframe$Sample,FUN=generate_align_command)%>% unlist(),
    "",
    "# pbindex bam files:",
    lapply(design_dataframe$Sample,FUN=generate_bpindex_command)%>% unlist(),
    "",
    "# Deactivate environment",
    "conda deactivate")%>%
    write_lines("../scripts/6_smrtlink_bash_part_2.txt")
}

#### Function for basemod calling commands
# Sub function for fai file creation
generate_fai_index_command <- function(string){
  return(paste0("bsub -q bio samtools faidx ../REF/",string,".fasta"))
}
# Sub function for basemod calling
generate_basemod_calling_command <- function(string){
  sample_number <- gsub(".*_","",string)
  sample_type <- gsub(sample_number,"",string)
  return(paste0(
    "bsub -q bio 'ipdSummary alignment/",sample_number,"_",sample_type,"aligned_pbmm2.bam \\\n",
    "--reference ../REF/", string,".fasta \\\n",
    "--numWorkers 6 \\\n",
    "--identify m6A,m4C \\\n",
    "--minCoverage 75 \\\n",
    "--mapQvThreshold 60 \\\n",
    "--methylFraction \\\n",
    "--pvalue 0.001 \\\n",
    "--outfile basemods/basemod_",string," \\\n",
    "--gff basemods/basemod_",string,".gff' \n"
  ))
}
# Main function
generate_smrtlink_commands_part_3 <- function(design_dataframe,directory){
  c("### Script for basemod calling with IPDSummary",
    paste0("### Script generated by the Salmonella Methylation Pipeline (SMP) on ",date()),
    "### GitHub url: github.com/AukeHaver/Salmonella_Methylation_Pipeline",
    "",
    "# Activate environment:",
    "conda activate smrtlink",
    "",
    "# Move to correct working directory:",
    paste0("cd ",file.path(directory,"smrtlink",fsep="/")),
    "",
    "# Make a basemods directory:",
    "mkdir basemods",
    "",
    "# create fai files for reference genomes:",
    lapply(design_dataframe$Sample,FUN=generate_fai_index_command) %>% unlist(),
    "",
    "# call basemods:",
    lapply(design_dataframe$Sample,FUN=generate_basemod_calling_command)%>% unlist(),
    "",
    "# Deactivate environment",
    "conda deactivate")%>%
    write_lines("../scripts/6_smrtlink_bash_part_3.txt")
}

#### Function for generating motif commands
# Sub function for find command
generate_find_command <- function(string){
  return(paste0("bsub -q bio motifMaker find -f ../REF/",string,".fasta -g basemods/basemod_",string,".gff -o motifs/motifs_",string,".csv"))
}
# Sub function for reprocess command
generate_reprocess_command <- function(string){
  return(paste0("bsub -q bio motifMaker reprocess -f ../REF/",string,".fasta -g basemods/basemod_",string,".gff -m motifs/motifs_",string,".csv -o motifs/motifs_",string,".gff"))
}
# Main function
generate_smrtlink_commands_part_4 <- function(design_dataframe,directory){
  c("### Script for motif analysis with Motifmaker",
    paste0("### Script generated by the Salmonella Methylation Pipeline (SMP) on ",date()),
    "### GitHub url: github.com/AukeHaver/Salmonella_Methylation_Pipeline",
    "",
    "# Activate environment:",
    "conda activate smrtlink",
    "",
    "# Move to correct working directory:",
    paste0("cd ",file.path(directory,"smrtlink",fsep="/")),
    "",
    "# Make a motifs directory",
    "mkdir motifs",
    "",
    "# Find motif and pass reference:",
    lapply(design_dataframe$Sample,FUN=generate_find_command)%>% unlist(),
    "",
    "# Reprocess gff files with motif information:",
    lapply(design_dataframe$Sample,FUN=generate_reprocess_command)%>% unlist(),
    "",
    "# Deactivate environment",
    "conda deactivate")%>%
    write_lines("../scripts/6_smrtlink_bash_part_4.txt")
}