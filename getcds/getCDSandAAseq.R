# getCDSandAAseq
# Main function get_sequences() is used to translate and transcribe transcript sequences from a GTF/GFF
# dependencies are instaeed through a docker image: docker pull murphydaviducl/getorf

options(warn = -1)

# Load libraries
# packages <- c("dplyr", 
#               "magrittr", 
#               "readr", 
#               "GenomicRanges", 
#               "BSgenome.Hsapiens.UCSC.hg38", 
#               "Biostrings", 
#               "ORFik", 
#               "optparse")
# 
# for (package in packages) {
#   if (!requireNamespace(package, quietly = TRUE)) {
#     install.packages(package)
#   }
# }

suppressPackageStartupMessages({
  library(dplyr)
  library(magrittr)
  library(readr)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(Biostrings)
  library(ORFik)
  library(optparse)
})

# Command options
option_list <- list(
  make_option("--gtf", type = "character", default = NULL,
              help = "Input GTF name"),
  make_option("--transcript", type = "character", default = NULL,
              help = "Input transcript file as TSV in a single column"),
  make_option("--output_csv", type = "character", default = NULL,
              help = "Output CSV file name"),
  make_option("--output_fasta", type = "character", default = NULL,
              help = "Output FASTA file name")
)

# Create option parser
parser <- OptionParser(option_list = option_list)

# Parse command-line arguments
args <- parse_args(parser)

# Read GTF file
if (!file.exists(args$gtf)) {
  stop(paste("Input GTF file", args$gtf, "does not exist."))
}

if (file.access(args$gtf, 4) == -1) {
  stop(paste("Input GTF file", args$gtf, "is not readable."))
}

gtf <- rtracklayer::import(args$gtf) # GTF input

# Read transcript file
if (!file.exists(args$transcript)) {
  stop(paste("Input transcript file", args$transcript, "does not exist."))
}

if (file.access(args$transcript, 4) == -1) {
  stop(paste("Input transcript file", args$transcript, "is not readable."))
}

# Read input data
# Transcript IDs
transcripts <- readr::read_lines(args$transcript) %>% as.character() # transcripts of interest

# GTF/GFF with transcript structures
gtf <- rtracklayer::import(args$gtf) # GTF file
gtf <- gtf[gtf$transcript_id %in% transcripts, ]

# Need to make sure chromosomes present in the genome reference match the input GTF/GFF. Have found that alt contigs does have non-matching IDs.
seqnames_not_in_ref <- seqnames(gtf) %>% 
                       unique() %>% 
                       setdiff(seqnames(BSgenome.Hsapiens.UCSC.hg38) %>% unique())

# Function to get nucleotide (NT), open reading frame (ORF) and amino acid (AA) sequences
get_sequences <- function(gtf){
  
  # Get sequences
  NT_seq <-
    subset(gtf, type == "exon") %>%
    split(., mcols(.)$transcript_id) %>%
    lapply(., function(x) {
      if (unique(strand(x)) == "-") {
        x %>%
          sort(., decreasing = TRUE) %>%  # Sort in decreasing order
          getSeq(Hsapiens, .) %>% paste(., collapse = "")
      } else {
        x %>%
          sort(.) %>%  # Sort in increasing order
          getSeq(Hsapiens, .) %>% paste(., collapse = "")
      }
    }) %>%
    utils::stack() %>%
    setNames(c("NT_seq", "Transcript")) %>%
    subset(., select=c("Transcript", "NT_seq")) 
  
  # ORF prediction, get ORF sequences and translate to AA sequences. 
  # I have used "ATG" as start codon and only use the longest predicted ORF
  
  NT_seq$CDS_seq <- NA
  NT_seq$AA_seq <- NA
  
  for (i in 1:NROW(NT_seq)) {
    ORF_range <- ORFik::findORFs(
      NT_seq$NT_seq[i], longestORF = TRUE, startCodon = "ATG") %>%
      unlist() %>%
      data.frame()
    
    # Check if the resulting data frame is empty and adjust accordingly
    if (nrow(ORF_range) == 0) {
      ORF_range <- data.frame(width = NA)
    } else {
      ORF_range <- ORF_range %>%
        dplyr::filter(width == max(width))
    }
    
    # Use ifelse to handle the condition and assign values accordingly
    NT_seq$CDS_seq[i] <- ifelse(is.na(ORF_range$width),
                                NA,
                                substr(as.character(NT_seq$NT_seq[i]),
                                       start = ORF_range$start, stop = ORF_range$end))
    
    # Get AA
    NT_seq$AA_seq[i] <- ifelse(is.na(ORF_range$width),
                               NA,
                               NT_seq$CDS_seq[i] %>% 
                                 DNAString() %>% 
                                 Biostrings::translate() %>% 
                                 as.character())
  }
   
  return(NT_seq)
}

# Run main function
output_data <-
  get_sequences(gtf = gtf)

# write output as CSV
write_csv(output_data, args$output_csv)

# write output as fasta
aastring <- AAStringSet(output_data$AA_seq)
names(aastring) <- output_data$Transcript

writeXStringSet(x = aastring, filepath = args$output_fasta, format = "fasta")