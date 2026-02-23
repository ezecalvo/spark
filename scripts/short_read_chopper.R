suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(R.utils))

option_list = list(
  make_option(c("--tsv"), type="character", default=NULL, 
              help="path to tsv with simulated mRNAs", metavar="character"),
  make_option(c("--insert_size"), type="character", default="200,300", 
              help="library insert size range, e.g. 200,300", metavar="character"),
  make_option(c("--read_length"), type="numeric", default=100, 
              help="length of each individual read", metavar="numeric"),
  make_option(c("-t", "--threads"), type="numeric", default=1, 
              help="Number of threads", metavar="numeric"),
  make_option(c("--seq_depth"), type="numeric", default=20000000, 
              help="sequencing depth", metavar="numeric"),
  make_option(c("--tpm_lower_limit"), type="numeric", default=5, 
              help="TPM lower limit", metavar="numeric"),
  make_option(c("--tpm_upper_limit"), type="numeric", default=200, 
              help="TPM upper limit", metavar="numeric"),
  make_option(c("-o", "--dir_out"), type="character", default=".", 
              help="dir out to create temp files", metavar="character"),
  make_option(c("--seed"), type="numeric", default=NULL,
              help="seed for randomization"),
  make_option(c("--fragments"), type="character", default='no', 
              help="to export ground truth", metavar="character"),
  make_option(c("--no_sizeselection"), action="store_true", default=FALSE,
              help="If specified, do not filter fragments by size (keep everything)"),
  make_option(c("--no_fragmentation"), action="store_true", default=FALSE,
              help="If specified, do not fragment; keep full molecule length")
)

#Parse input
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

if (is.null(opt$tsv)){
  stop("path to tsv with simulated mRNAs is missing", call.=FALSE)
}

if (!is.null(opt$seed)){
  set.seed(opt$seed)
}

setDTthreads(threads = opt$threads)
insert_size <- opt$insert_size
insert_size_split <- strsplit(insert_size, ",")[[1]]
insert_size <- c(as.numeric(insert_size_split[1]), as.numeric(insert_size_split[2]))

file <- opt$tsv

#File importing and formating
file_importer <- function(file){
  full_reads <- data.table::fread(file=file, sep = '\t', header = T, stringsAsFactors = FALSE)
  # Ensure transcript_id matches the row index logic used later
  full_reads$transcript_id <- as.numeric(1:nrow(full_reads))
  full_reads$sequence_length <- nchar(full_reads$full_molecule_sequence)
  return(full_reads)
}

#Chopper function
get_reads <- function(lengths, eta_val = 200, insert_size, transcript_id, no_fragmentation=FALSE, no_sizeselection=FALSE) {
  
  if (no_fragmentation) {
    # return a single fragment from 1 to the end of the molecule
    fragments <- data.frame(
      transcript_id = transcript_id,
      read_start = 1,
      read_end = lengths,
      length = lengths, # length is the full length
      size_selection = 'passed_size_selection' 
    )
    
    if (!no_sizeselection) {
      if (fragments$length < insert_size[1] || fragments$length > insert_size[2]) {
        fragments$size_selection <- 'filtered_out'
      }
    }
    
    
    fragments$length <- NULL 
    
    return(fragments)
  }
  
  eta_val <- mean(c(insert_size[1],insert_size[2]))
  
  deltas <- log10(lengths)
  ns_minus_1 <- pmax(round(lengths / eta_val / gamma(1 / deltas + 1)) - 1, 0)
  xis <- lapply(ns_minus_1, function(n) { diff(sort(c(runif(n), 0, 1))) })
  xis_transformed <- mapply(function(x, d) { x^(1 / d) }, xis, deltas, SIMPLIFY = FALSE)
  delta_is <- mapply(function(len, x_t) {
    round(len * x_t / sum(x_t))
  }, lengths, xis_transformed, SIMPLIFY = FALSE)
  
  # Get all the start and end points of the fragments
  starts <- lapply(delta_is, function(d) {
    if (length(d) > 1) {
      c(sample(min(insert_size[1], d[1]), 1), cumsum(d[1:(length(d) - 1)]))
    } else {
      sample(min(insert_size[1], d), 1)
    }
  })
  
  ends <- lapply(delta_is, function(d) {
    if (length(d) > 1) {
      c(cumsum(d[1:(length(d) - 1)]), sum(d) - sample(min(insert_size[1], sum(d) - d[length(d)]), 1))
    } else {
      d
    }
  })
  
  fragments <- data.frame(
    transcript_id = rep(transcript_id, length(delta_is[[1]])),
    read_start = unlist(starts),
    read_end = unlist(ends)
  )
  
  # Create first row from 1 to the first read_start
  first_row <- data.frame(
    transcript_id = fragments$transcript_id[1],
    read_start = 1,
    read_end = fragments$read_start[1]
  )
  
  # Create last row from the last read_end to the value of lengths
  last_row <- data.frame(
    transcript_id = fragments$transcript_id[nrow(fragments)],
    read_start = fragments$read_end[nrow(fragments)],
    read_end = lengths
  )
  
  # Combine all together
  fragments <- rbind(first_row, fragments, last_row)
  fragments$length <- fragments$read_end - fragments$read_start
  
  # Filter by insert size
  fragments$size_selection <- 'filtered_out'
  
  
  if (no_sizeselection) {
    # ff the flag is set, everything passes the selection
  } else {
    fragments$size_selection <- 'passed_size_selection'
    fragments[fragments$length >= insert_size[1] & fragments$length <= insert_size[2], 'size_selection'] <- 'passed_size_selection'
  }
  
  fragments$length <- NULL  # Optional: drop if not needed
  
  return(fragments)
}


random_string_gen <- function(n = 5000) {
  a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
  paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}

full_reads <- file_importer(file)
#Get gene that was simulated
gene_id <- sub(".*/(ENSG[0-9]+(?:_background)?)(?:\\.tsv\\.gz)$", "\\1", file)

list_with_all_fragments <- lapply(seq_len(nrow(full_reads)), function(i) {
  get_reads(lengths = full_reads$sequence_length[i],
            insert_size = insert_size,
            transcript_id = full_reads$transcript_id[i],
            # Pass the new options from the opt object
            no_fragmentation = opt$no_fragmentation,
            no_sizeselection = opt$no_sizeselection)
})

reads_list <- dplyr::bind_rows(list_with_all_fragments,.id = 'transcript_id')

#Keep all fragments before filtering by their length
reads_list_all_fragments <- reads_list
#Keep fragments that are long enough
reads_list <- reads_list[reads_list$size_selection=='passed_size_selection',]


if (nrow(reads_list)>10){ #Some genes have really short mRNAs and will get filtered out during size selection
  reads_list$transcript_id <- as.numeric(reads_list$transcript_id)
  #Assign a TPM and sample the number of reads (rows) accordingly
  gene_length <- mean(full_reads$sequence_length)/1000
  gene_tpm <- runif(n=1, min=as.integer(opt$tpm_lower_limit), max=as.integer(opt$tpm_upper_limit))
  seq_depth <- opt$seq_depth/10^6
  
  #Get the number of fragments for the gene. 
  reads_to_get <- gene_length*gene_tpm*seq_depth
  
  reads_list$read_start <- as.integer(reads_list$read_start)
  reads_list$read_end <- as.integer(reads_list$read_end)
  
  if (opt$fragments=='with_ground_truth'){
    #Get a ground truth that belongs to fragmented and size selected fragments (but not mapped)
    all_positions_selected <- unlist(mapply(seq, reads_list$read_start, reads_list$read_end))
    freq_table_selected <- as.data.table(table(all_positions_selected))
    setnames(freq_table_selected, c("position", "frequency"))
    freq_table_selected[, position := as.integer(as.character(position))]
    
    dir_ground_truth <- paste(opt$dir_out,'/ground_truth_after_size_selection',sep = '')
    dir.create(dir_ground_truth, showWarnings = FALSE, recursive = TRUE)
    path_for_file <- paste(dir_ground_truth,'/',gene_id,'.tsv.gz',sep = '')
    fwrite(freq_table_selected,path_for_file,sep = '\t',col.names = T,row.names = F,quote = F,compress = 'gzip')
  }
  
  
  reads_list_collapsed <- reads_list %>%
    dplyr::mutate(read_coordinate = paste0(read_start, "-", read_end)) %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::summarise(read_coordinates = paste(read_coordinate, collapse = ","), .groups = "drop")
  
  #### Get ground truth for all the fragments before fragmentation and size selection
  if (opt$fragments=='with_ground_truth'){
    
    # Avoid division by zero if all fragments pass size selection (which happens with no_sizeselection)
    n_passed <- nrow(reads_list_all_fragments[reads_list_all_fragments$size_selection=='passed_size_selection',])
    if(n_passed == 0) {
      ratio_total_passed <- 0
    } else {
      ratio_total_passed <- nrow(reads_list_all_fragments)/n_passed
    }
    
    #Sample the total reads based on the ratio and the assigned gene TPM
    all_fragments_to_get <- reads_to_get
    
    #Sample multiple times to smooth the distribution
    all_sampling <- data.frame()
    for (j in 1:100){
      if (nrow(reads_list_all_fragments)>=all_fragments_to_get){
        reads_list_all_fragments_sampled <- dplyr::sample_n(reads_list_all_fragments,size=all_fragments_to_get,replace = F) #No replacement
      }else{
        sampled_reads <- dplyr::sample_n(reads_list_all_fragments,size=all_fragments_to_get-nrow(reads_list_all_fragments),replace = T) #With replacement for the reads that are missing
        reads_list_all_fragments_sampled <- rbind(reads_list_all_fragments,sampled_reads)
      }
      
      #Get ground truth for the number of mRNAs that are going into this -> used for edge effect correction
      all_positions <- unlist(mapply(seq, reads_list_all_fragments_sampled$read_start, reads_list_all_fragments_sampled$read_end))
      freq_table <- as.data.table(table(all_positions))
      setnames(freq_table, c("position", "frequency"))
      freq_table[, position := as.integer(as.character(position))]
      all_sampling <- rbind(all_sampling,freq_table)
    }
    final_freqs <- all_sampling %>% dplyr::group_by(position) %>% dplyr::reframe(frequency=mean(frequency))
    final_freqs <- final_freqs[-nrow(final_freqs),]#Remove last nucleotide which will always be noisy
    
    dir_ground_truth <- paste(opt$dir_out,'/ground_truth_pre_size_selection',sep = '')
    dir.create(dir_ground_truth, showWarnings = FALSE, recursive = TRUE)
    path_for_file <- paste(dir_ground_truth,'/', gene_id,'.tsv.gz',sep = '')
    fwrite(final_freqs,path_for_file,sep = '\t',col.names = T,row.names = F,quote = F,compress = 'gzip')
  }
  
  temp_dir <- paste(opt$dir_out,'/temp/mRNAs_with_fragments',sep = '')
  dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
  path_for_file <- paste(temp_dir,'/',gene_id,'_fragments.tsv',sep = '')
  
  fwrite(reads_list_collapsed, path_for_file, sep = '\t', col.names = T, row.names = F, quote = F)
  
}


