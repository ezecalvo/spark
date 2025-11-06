#!/usr/bin/env Rscript
# Author: Shaimae Elhajjajy
#Modified by Ezequiel Calvo-Roitberg

# Purpose: This R script will analyze the GENCODE annotations for the most highly expressed genes in a cell line the goal is to curate a robust, diverse, and representative gene 
# set to be used in simulations.

### SET UP ###

# Import libraries
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(scales))

#Parse arguments
#####
option_list = list(
  make_option(c("-g", "--gtf"), type="character", default=NULL, 
              help="genome annotation", metavar="character"),
  make_option(c("-f", "--genome_fasta"), type="character", default=NULL, 
              help="genome fasta", metavar="character"),
  make_option(c("-o", "--dir_out"), type="character", default=".", 
              help="dir out for the bed and metric files", metavar="character"),
  make_option(c("--protein_coding_only"), action="store_true", default=FALSE,
              help="only use protein-coding genes"),
  make_option(c("-t", "--threads"), type="numeric", default="1", 
              help="Number of threads", metavar="numeric"),
  make_option(c("-c", "--n_gene_clusters"), type="numeric", default="6", 
              help="Number of gene clusters to get", metavar="numeric"),
  make_option(c("-n", "--number_of_genes"), type="numeric", default="120", 
              help="Number of total genes to simulate", metavar="numeric"),
  make_option(c("--seed"), type="numeric", default=NULL,
              help="seed for randomization")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$gtf)){
  stop("GTF file is missing", call.=FALSE)
}
if (is.null(opt$genome_fasta)){
  stop("Genome fasta file is missing", call.=FALSE)
}

if (!is.null(opt$seed)){
  set.seed(opt$seed)
}

#Set number of threads for data table
data.table::setDTthreads(threads = opt$threads)

#Make directory to temp files
temp_dir <- paste(opt$dir_out,'/temp',sep = '')
system(paste("mkdir ", temp_dir,sep = ''))

# Read in the gencode annotation data
gencode_data = fread(opt$gtf, sep = "\t",header = F)

### FUNCTIONS ###

# Extract all GENCODE gene IDs so the RNA-seq gene IDs can be searched for in the GENCODE annotation
getGencodeGeneIDs <- function(gencode_data){
  gene_ids_gencode <- sub(".*gene_id \\\"(ENSG[0-9]+)\\\".*", "\\1", gencode_data$V9)
  return(gene_ids_gencode)
}

getGencodeTranscriptIDs <- function(gencode_data){
  transcript_id <- sub('.*transcript_id "([^"]+)".*', '\\1', gencode_data$V9)
  return(transcript_id)
}

# Find the type of gene that the longest transcript is from (as opposed to the transcript type)
getGeneType <- function(gene_annotation){
  gene_type <- sub('.*gene_biotype "([^"]+)".*', '\\1', gencode_data$V9)  
  return(gene_type)
  
}

# Find the longest transcript for the gene
getLongestTranscript <- function(gencode_data){
  longest_txpt <- gencode_data[gencode_data$V3=='transcript',] %>% dplyr::group_by(gene_id) %>% dplyr::mutate(transcript_length=V5-V4) %>% dplyr::filter(transcript_length==max(transcript_length))
  return(longest_txpt)
  
}

getfeature_metrics <- function(gencode_data,longest_txpt){
  all_features <- gencode_data[gencode_data$transcript_id%in%longest_txpt$transcript_id,] #Filter to make it faster
  all_exons <- all_features[all_features$V3=='exon',]
  all_exons$exon_length <- all_exons$V5-all_exons$V4
  #Split by strand and sort (doesn't have to be per gene)
  exons_pos_strand <- all_exons[all_exons$V7=='+',]
  exons_neg_strand <- all_exons[all_exons$V7=='-',]
  #Sort
  exons_pos_strand <- exons_pos_strand[order(exons_pos_strand$V4,decreasing = F),]
  exons_neg_strand <- exons_neg_strand[order(exons_neg_strand$V5,decreasing = T),]
  
  #Assign genomic order
  exons_pos_strand <- exons_pos_strand %>% 
    dplyr::group_by(transcript_id) %>% 
    dplyr::mutate(genomic_order=dplyr::row_number())
  
  exons_neg_strand <- exons_neg_strand %>% 
    dplyr::group_by(transcript_id) %>% 
    dplyr::mutate(genomic_order=dplyr::row_number())
  
  all_exons <- rbind(exons_pos_strand,exons_neg_strand)
  
  exon_coordinates <- data.frame(chr=all_exons$V1,start=all_exons$V4,end=all_exons$V5,id=all_exons$transcript_id,score=all_exons$genomic_order,strand=all_exons$V7)
  
  #exon_coordinates$start <- exon_coordinates$start-1
  #exon_coordinates$end <- exon_coordinates$end-1
  
  #Get exons and introns metrics
  exon_intron_metrics <- all_exons %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::mutate(intron_length = dplyr::case_when(
      V7 == "+" ~ dplyr::lead(V4) - V5,
      V7 == "-" ~ V4 - dplyr::lead(V5),
      TRUE ~ NA_real_
    )) %>%
    dplyr::reframe(
      number_of_exons=dplyr::n(),
      number_of_introns=number_of_exons-1,
      mean_exon_length = mean(exon_length, na.rm = TRUE),
      mean_intron_length = mean(intron_length, na.rm = TRUE),
      first_exon_length = dplyr::first(exon_length),
      first_intron_length = dplyr::first(intron_length)
    )
  
  #Get intron coordinates
  intron_coordinates_raw <- all_exons[,c(1,4,5,14,11,7)] %>%
    dplyr::group_by(transcript_id) %>% 
    dplyr::mutate(
      V4_new = dplyr::case_when(
        V7 == "+" ~ V5 + 1,
        V7 == "-" ~ dplyr::lead(V5) + 1,
        TRUE ~ NA_real_
      ),
      V5_new = dplyr::case_when(
        V7 == "+" ~ dplyr::lead(V4) - 1,
        V7 == "-" ~ V4 - 1,
        TRUE ~ NA_real_
      )) %>% 
    dplyr::filter(!is.na(V4_new) & !is.na(V5_new))
  
  intron_coordinates <- data.frame(chr=intron_coordinates_raw$V1,start=intron_coordinates_raw$V4_new,end=intron_coordinates_raw$V5_new,id=intron_coordinates_raw$transcript_id,score=intron_coordinates_raw$genomic_order,strand=intron_coordinates_raw$V7)
  #intron_coordinates$start <- intron_coordinates$start-1
  #intron_coordinates$end <- intron_coordinates$end-1
  
  #Get UTR metrics
  utrs <- all_features[all_features$V3%in%c("five_prime_utr","three_prime_utr"),]
  utrs$utr_length <- utrs$V5-utrs$V4
  #GTFs normally might have different features (rows) for utrs if they across multiple exons
  utrs <- utrs %>% dplyr::group_by(transcript_id,V3) %>% dplyr::reframe(utr_length=sum(utr_length))
  utrs_pivot <- tidyr::pivot_wider(utrs[,c('transcript_id','utr_length','V3')],id_cols = transcript_id,values_from = utr_length,names_from = V3)
  
  return(list(exon_intron_metrics=exon_intron_metrics,utrs_pivot=utrs_pivot,exon_coordinates=exon_coordinates,intron_coordinates=intron_coordinates))
}

getLongestTxptSequence <- function(longest_txpt,exon_coordinates,intron_coordinates,genome_fasta,temp_dir){
  
  coordinates_transcripts_for_bedtools <- data.frame(chr=longest_txpt$V1,start=longest_txpt$V4,end=longest_txpt$V5,id=paste(longest_txpt$transcript_id,'transcript',sep = '_'),score=longest_txpt$V6,strand=longest_txpt$V7)
  
  #coordinates_transcripts_for_bedtools$start <- coordinates_transcripts_for_bedtools$start-1
  #coordinates_transcripts_for_bedtools$end <- coordinates_transcripts_for_bedtools$end-1
  
  exon_coordinates$id <- paste(exon_coordinates$id,'exon',exon_coordinates$score,sep = '_')
  intron_coordinates$id <- paste(intron_coordinates$id,'intron',intron_coordinates$score,sep = '_')
  
  #Merge all coordinates to run bedtools only once (parsing happens downstream)
  all_coordinates <- rbind(coordinates_transcripts_for_bedtools,exon_coordinates,intron_coordinates)
  all_coordinates$start <- as.integer(all_coordinates$start)
  all_coordinates$end <- as.integer(all_coordinates$end)
  
  all_coordinates$start <- all_coordinates$start-1
  
  path_for_coordinates <- paste(temp_dir,"/longest_transcript_and_features.bed",sep = '')
  fwrite(all_coordinates,path_for_coordinates, sep = "\t", append = FALSE, col.names = FALSE,row.names = F)
  
  #Run nucleotide metrics
  get_script_path <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- grep("^--file=", args, value = TRUE)
    
    if (length(file_arg) > 0) {
      # If running via Rscript
      return(normalizePath(sub("^--file=", "", file_arg)))
    } else if (interactive()) {
      # If running interactively (e.g., RStudio)
      if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
        return(rstudioapi::getSourceEditorContext()$path)
      } else {
        stop("Cannot determine script path in interactive mode without rstudioapi.")
      }
    } else {
      stop("Script path could not be determined.")
    }
  }
  
  script_path <-dirname(get_script_path())
  
  system(paste("python ", script_path,"/nucleotide_metrics.py ", path_for_coordinates,' ',genome_fasta,' ',temp_dir,sep = ''))
  
  
  #Read nucleotide metrics
  nucleotide_metrics_path <- paste(temp_dir,'/nucleotide_metrics.tsv',sep='')
  nucleotide_metrics <- fread(nucleotide_metrics_path,sep = '\t',header = T)
  #system(paste('rm',nucleotide_metrics_path,sep=' '))
  
  #Parse and pivot
  nucleotide_metrics$Feature <- gsub("(_\\d+\\(.+\\))$", "", nucleotide_metrics$Feature)
  nucleotide_metrics$transcript_id <- sub("(ENST\\d+)_([a-zA-Z]+).*", "\\1", nucleotide_metrics$Feature)
  nucleotide_metrics$feature <- sub("(ENST\\d+)_([a-zA-Z]+).*", "\\2", nucleotide_metrics$Feature)

  nucleotide_metrics_pivot <- nucleotide_metrics %>% 
    dplyr::select(-Feature) %>% 
    tidyr::pivot_longer(cols = A:G, names_to = "base", values_to = "value") %>%
    dplyr::group_by(transcript_id,feature,base) %>% 
    dplyr::reframe(new_value=mean(value)) %>% 
    tidyr::pivot_wider(names_from = c(feature, base), values_from = new_value)
  
  return(nucleotide_metrics_pivot)
  
}

#Wrangle to make parsing more efficient
gencode_data$gene_id <- getGencodeGeneIDs(gencode_data)
gencode_data$transcript_id <- getGencodeTranscriptIDs(gencode_data)
gencode_data$gene_type  <- getGeneType(gencode_data)

#Optional protein coding filter
if (opt$protein_coding_only) {
  gencode_data <- gencode_data[gencode_data$gene_type=='protein_coding',]
  
}

#These are modular metrics instead of a unique function to make it more customizable
#Get gene length
gene_length_ref <- gencode_data[gencode_data$V3=='gene',] %>% dplyr::group_by(gene_id) %>% dplyr::reframe(gene_length=V5-V4)

#Find the longest transcript for the gene
longest_txpt <- getLongestTranscript(gencode_data)

#Add longest transcript info to gene ids to generate the final df with metrics
db_gene_transcript <- gencode_data[gencode_data$transcript_id%in%longest_txpt$transcript_id,c('gene_id','transcript_id','V3','gene_type')]
db_gene_transcript <- db_gene_transcript[db_gene_transcript$V3=='transcript',]
gene_length_ref <- merge(gene_length_ref,db_gene_transcript[,-'V3'],by='gene_id')

#Get the number of exons and introns per transcript
feature_metrics <- getfeature_metrics(gencode_data,longest_txpt)
exon_introns_metrics <- feature_metrics[["exon_intron_metrics"]]
#Get length of UTRs
utr_metrics <- feature_metrics[["utrs_pivot"]]

utr_gene_has <- colnames(utr_metrics)[-1]

exon_coordinates <- feature_metrics[["exon_coordinates"]]
intron_coordinates <- feature_metrics[["intron_coordinates"]]

#Get longest transcript sequence and metrics then parse
nucleotide_metrics_pivot <- getLongestTxptSequence(longest_txpt,exon_coordinates,intron_coordinates,opt$genome_fasta,temp_dir)
#nucleotide_metrics_pivot <- getLongestTxptSequence(longest_txpt,exon_coordinates,intron_coordinates,'/pi/athma.pai-umw/genomes/hg38/Homo_sapiens.GRCh38.dna.primary_assembly.fa','/pi/athma.pai-umw/analyses/ezequiel/kinetic_barcoding/simulations/drb_4sU_simulation/AEBP2/min4/temp/')

#Build the output
metric_df <- merge(merge(merge(gene_length_ref,merge(longest_txpt[,c(1,4,5,11,6,7,13)],exon_introns_metrics,by='transcript_id'),by='transcript_id'),utr_metrics,by='transcript_id',all.x=T),nucleotide_metrics_pivot,by='transcript_id')

colnames(metric_df) = c('transcript_id','gene_id','gene_length','gene_type','chr','start','end','score','strand','transcript_length','number_of_exons','number_of_introns','mean_exon_length','mean_intron_length','first_exon_length', 'first_intron_length', utr_gene_has, 'transcript_A', 'transcript_C', 'transcript_T','transcript_G', 'exon_A', 'exon_C', 'exon_T', 'exon_G', 'intron_A', 'intron_C', 'intron_T','intron_G')

#Write metrics
path_metrics <- paste(opt$dir_out,'/metrics_per_transcript_df.tsv',sep='')
fwrite(metric_df, path_metrics, sep = "\t")

######
#Hierarchical clustering

data_scaling <- function(metric_df){
  #Create a separate dataframe with only the variables that will be used in the hiearchical clustering
  safe_cols <- c("transcript_id", "gene_length", "transcript_length", "number_of_exons", "number_of_introns",
                    "mean_exon_length", "mean_intron_length", "first_exon_length", "first_intron_length",
                    "five_prime_utr", "three_prime_utr",
                    "transcript_A", "transcript_C", "transcript_T", "transcript_G",
                    "exon_A", "exon_C", "exon_T", "exon_G",
                    "intron_A", "intron_C", "intron_T", "intron_G")
  safe_cols_present <- intersect(safe_cols, colnames(metric_df))
  data_to_scale <- metric_df[, ..safe_cols_present]

  
  #Replace NAs by 0 for genes that don't have a given feature (for example intronless genes)
  data_to_scale <- replace(data_to_scale, is.na(data_to_scale), 0)
  # Standardize the data by subtracting the mean and dividing by the standard deviation
  data_scaled <- scale(data_to_scale[, 2:ncol(data_to_scale)], center = TRUE, scale = TRUE)
  
  #For positive scaled values, take the log; for negative scaled values, take the absolute value and then the negative log
  log10_transform <- function(x) {
    sign(x) * log10(pmax(abs(x), .Machine$double.eps))
  }
  
  data_scaled <- log10_transform(data_scaled)
  data_scaled <- cbind(data_to_scale[,c("transcript_id")],data_scaled)
  
  return(data_scaled) 
}

cluster_making <- function(data_scaled,gene_clusters){
  #Run HC
  rownames_data_matrix <- data_scaled$transcript_id
  data_matrix <- as.matrix(data_scaled[, -c(1)])
  rownames(data_matrix) <- rownames_data_matrix
  data.clust <- hclust(d = dist(x = data_matrix), method = "complete")
  #Get a given number of clusters
  gene_clusters_df = cutree(data.clust, k = gene_clusters)
  #Wrangle
  gene_clusters_df <- as.data.frame(gene_clusters_df)
  gene_clusters_df$transcript_id <- rownames(gene_clusters_df)
  rownames(gene_clusters_df) <- NULL
  return(gene_clusters_df) 
}

gene_sampling <- function(gene_clusters_df, number_of_genes, n_gene_clusters) {
  # Calculate the initial target per cluster
  n_genes_per_cluster_to_sample <- floor(number_of_genes / n_gene_clusters)
  n_genes_per_cluster <- table(gene_clusters_df$gene_clusters) %>% as.data.frame()
  colnames(n_genes_per_cluster) <- c('gene_clusters','total_genes')
  
  # Identify clusters with fewer genes than the target and calculate leftovers
  n_genes_per_cluster <- n_genes_per_cluster %>%
    dplyr::mutate(
      to_sample = pmin(total_genes, n_genes_per_cluster_to_sample),  # Limit to available genes
      leftover = pmax(0, n_genes_per_cluster_to_sample - total_genes) # Shortfall per cluster
    )
  
  # Calculate total leftover genes
  total_leftover <- sum(n_genes_per_cluster$leftover)
  
  # Redistribute leftover genes to clusters that can accommodate more
  n_genes_per_cluster <- n_genes_per_cluster %>%
    dplyr::mutate(
      can_take_more = total_genes - to_sample,
      additional = pmin(total_leftover, can_take_more))%>% 
    dplyr::mutate(additional=floor(pmin(total_leftover, can_take_more)/length(additional>0)))
  
  # Adjust to_sample based on redistribution
  n_genes_per_cluster <- n_genes_per_cluster %>%
    dplyr::mutate(
      to_sample = to_sample + additional
    )
  
  # Perform sampling based on adjusted to_sample per cluster
  sampled_genes <- merge(gene_clusters_df, n_genes_per_cluster, by = "gene_clusters") %>%
    dplyr::group_by(gene_clusters) %>%
    dplyr::group_map(~ dplyr::slice_sample(.x, n = .x$to_sample[1]), .keep = TRUE) %>%
    dplyr::bind_rows() %>%
    dplyr::select(gene_clusters, transcript_id)
  
  return(sampled_genes)
}

gtf_generator <- function(path_to_fastatsv,sampled_genes,exon_coordinates,intron_coordinates,db_gene_transcript,path_for_gtf){
  #Paese the tsv to only have features
  path_to_parsed_file <- paste(path_to_fastatsv,'longest_transcript_only_features.tsv',sep='')
  command_for_parsingtsv <- paste("awk -F'\t' '$1 !~ /transcript/' ", path_to_fastatsv,'longest_transcript_and_features.tsv ','>',path_to_parsed_file,sep = '')
  system(command_for_parsingtsv)
  
  #Import and parse
  seq_per_feature <- fread(path_to_parsed_file,header = F)
  seq_per_feature <- seq_per_feature %>% 
    tidyr::separate(V1, into = c("transcript_id", "feature", "position"), sep = "_") %>%
    dplyr::mutate(position = gsub("\\([+-]\\)", "", position))
  
  colnames(seq_per_feature)[4] <- 'sequence'
  #Filter
  seq_per_feature <- seq_per_feature[seq_per_feature$transcript_id%in%sampled_genes$transcript_id,]
  
  #Add coordinate info per feature
  exon_coordinates$feature <- 'exon'
  intron_coordinates$feature <- 'intron'
  all_features <- rbind(exon_coordinates,intron_coordinates)
  colnames(all_features)[c(4,5)] <- c('transcript_id','position')
  #Add gene id for export
  all_features <- merge(all_features,db_gene_transcript[,c('transcript_id','gene_id')],by='transcript_id')
  
  seq_per_feature_full <- merge(seq_per_feature,all_features,by=c('transcript_id','feature','position'))
  seq_per_feature_full_for_export <- seq_per_feature_full[,c('chr','start','end','gene_id','feature','position','strand','sequence')]
  
  #Generate 0 based coordinates
  #seq_per_feature_full_for_export$start <- seq_per_feature_full_for_export$start-1
  #seq_per_feature_full_for_export$end <- seq_per_feature_full_for_export$end-1
  
  #Sort
  seq_per_feature_full_for_export_pos <- seq_per_feature_full_for_export[seq_per_feature_full_for_export$strand=='+',]
  seq_per_feature_full_for_export_neg <- seq_per_feature_full_for_export[seq_per_feature_full_for_export$strand=='-',]
  seq_per_feature_full_for_export_pos <- seq_per_feature_full_for_export_pos[order(seq_per_feature_full_for_export_pos$start,decreasing = F),]
  seq_per_feature_full_for_export_neg <- seq_per_feature_full_for_export_neg[order(seq_per_feature_full_for_export_neg$end,decreasing = T),]
  seq_per_feature_full_for_export <- rbind(seq_per_feature_full_for_export_pos,seq_per_feature_full_for_export_neg)
  #Generate one file per transcript id
  path_for_gtf <- paste(path_for_gtf,'gtf/',sep='')
  system(paste('mkdir ',path_for_gtf,sep=''))
  
  seq_per_feature_full_for_export %>% 
    dplyr::group_by(gene_id) %>%
    dplyr::group_walk(~ write.table(.x,
                                    file = paste(path_for_gtf, .y$gene_id, ".tsv", sep=''), 
                                    sep = "\t", 
                                    row.names = FALSE, 
                                    quote = FALSE,col.names=T),.keep = TRUE)
  
  system(paste('gzip ',path_for_gtf,'/*.tsv',sep='')) #gzip everything
  system(paste('rm ',path_to_fastatsv,'longest_transcript_and_features.tsv',sep='')) #Remove temp file
  
}

# Standardize the data by subtracting the mean and dividing by the standard deviation
metric_df <- fread(path_metrics,sep='\t',header = T)
data_scaled <- data_scaling(metric_df)
#Run HC and get a given number of clusters
if (nrow(data_scaled)>=opt$n_gene_clusters & nrow(data_scaled)>1){
  gene_clusters_df <- cluster_making(data_scaled,opt$n_gene_clusters)
  colnames(gene_clusters_df)[1] <- 'gene_clusters'
#Sample clusters
  sampled_genes <- gene_sampling(gene_clusters_df, opt$number_of_genes, opt$n_gene_clusters)
}else{
  sampled_genes <- dplyr::sample_n(data_scaled,opt$number_of_genes)
  
}

#############
#GTF for mRNA generation

gtf_generator(paste(temp_dir,'/',sep=''),sampled_genes,exon_coordinates,intron_coordinates,db_gene_transcript,opt$dir_out)
















