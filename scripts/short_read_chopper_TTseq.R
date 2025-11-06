

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
              help="sequencing depth", metavar="numeric"),
  make_option(c("--tpm_upper_limit"), type="numeric", default=200, 
              help="sequencing depth", metavar="numeric"),
  make_option(c("-o", "--dir_out"), type="character", default=".", 
              help="dir out to create temp files", metavar="character"),
  make_option(c("--seed"), type="numeric", default=NULL,
              help="seed for randomization")
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
  full_reads$transcript_id <- as.numeric(1:nrow(full_reads))
  full_reads$sequence_length <- nchar(full_reads$full_molecule_sequence)
  return(full_reads)
}

#Chopper function
get_reads <- function(lengths, eta_val = 200, insert_size, transcript_id) {
  
  eta_val <- mean(c(insert_size[1],insert_size[2]))
  # Select a fragment from each transcript, size select, and return the start and end positions
  # Inputs: 
  #   lengths - the length of the transcript
  #   eta_val - eta parameter for Weibull distribution
  #   insert_size - vector of length two with the lower and upper bounds for valid fragments
  #   transcript_id - the transcript ID to assign to each resulting fragment
  # Output:
  #   A data frame with transcript_id, read_start, and read_end
  
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
  fragments$size_selection <- 'passed_size_selection'
  #fragments[fragments$length >= insert_size[1] & fragments$length <= insert_size[2], 'size_selection'] <- 'passed_size_selection'
  fragments$length <- NULL  # Optional: drop if not needed
  
  return(fragments)
}


random_string_gen <- function(n = 5000) {
  a <- do.call(paste0, replicate(5, sample(LETTERS, n, TRUE), FALSE))
  paste0(a, sprintf("%04d", sample(9999, n, TRUE)), sample(LETTERS, n, TRUE))
}

full_reads <- file_importer(file)
#Get gene that was simulated
gene_id <- sub(".*\\/([A-Z]+[0-9]+(?:_[a-z]+)?)\\.tsv\\.gz$", "\\1", file)


list_with_all_fragments <- lapply(seq_len(nrow(full_reads)), function(i) {
  get_reads(lengths = full_reads$sequence_length[i],
            insert_size = insert_size,
            transcript_id = full_reads$transcript_id[i])
})

reads_list <- dplyr::bind_rows(list_with_all_fragments,.id = 'transcript_id')

#Keep all fragments before filtering by their length
reads_list_all_fragments <- reads_list
#Keep fragments that are long enough
reads_list <- reads_list[reads_list$size_selection=='passed_size_selection',]


if (nrow(reads_list)>10){ #Some genes have really short mRNAs and will get filtered out during size selection
  reads_list$transcript_id <- as.numeric(reads_list$transcript_id)
  
  #TPM assignment for TT-seq happens after fragmentation and 4sU PD
  reads_list$read_start <- as.integer(reads_list$read_start)
  reads_list$read_end <- as.integer(reads_list$read_end)
  
  #Keep fragments longer than 10nt
  reads_list <- reads_list %>% dplyr::filter(read_end-read_start>10)
  
  #Get a ground truth that belongs to fragmented and size selected fragments (but not mapped)
  # all_positions_selected <- unlist(mapply(seq, reads_list$read_start, reads_list$read_end))
  # freq_table_selected <- as.data.table(table(all_positions_selected))
  # setnames(freq_table_selected, c("position", "frequency"))
  # freq_table_selected[, position := as.integer(as.character(position))]
  # 
  # dir_ground_truth <- paste(opt$dir_out,'/ground_truth_after_size_selection',sep = '')
  # dir.create(dir_ground_truth, showWarnings = FALSE, recursive = TRUE)
  # path_for_file <- paste(dir_ground_truth,'/',gene_id,'.tsv.gz',sep = '')
  # fwrite(freq_table_selected,path_for_file,sep = '\t',col.names = T,row.names = F,quote = F,compress = 'gzip')
  
  reads_list <- reads_list %>%
    dplyr::mutate(read_coordinate = paste0(read_start, "-", read_end)) %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::summarise(read_coordinates = paste(read_coordinate, collapse = ","), .groups = "drop")

  
  full_reads$transcript_id <- as.numeric(full_reads$transcript_id)
  reads_list$transcript_id <- as.numeric(reads_list$transcript_id)
  
  full_reads <- merge(full_reads,reads_list,by='transcript_id')
  
  #Export
  temp_dir <- paste(opt$dir_out,'/temp/mRNAs_with_fragments',sep = '')
  dir.create(temp_dir, showWarnings = FALSE, recursive = TRUE)
  path_for_file <- paste(temp_dir,'/',gene_id,'_fragments.tsv',sep = '')
  fwrite(full_reads,path_for_file,sep = '\t',col.names = T,row.names = F,quote = F)
}


















