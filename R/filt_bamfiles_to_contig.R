library(Rsamtools)
library(GenomicRanges)
library(BiocParallel)
library(tools)

# Function to filter BAM file to a specific contig
filt_bamfiles_to_contig <- function(bamfile, contig_name, output_dir) {
  # Check if BAM file exists
  if (!file.exists(bamfile)) {
    stop(paste("BAM file not found:", bamfile))
  }

  # Generate output file name
  bamfile_base <- basename(bamfile)
  output_bam <- file.path(output_dir, paste0(file_path_sans_ext(bamfile_base), "_", contig_name, ".bam"))

  # Get sequence lengths from BAM file header
  seq_lengths <- scanBamHeader(bamfile)[[1]]$targets

  # Check if the contig exists in the BAM file
  if (!contig_name %in% names(seq_lengths)) {
    warning(paste("Contig", contig_name, "not found in BAM file:", bamfile))
    return(NULL)  # Skip this BAM file
  }

  # Get the length of the contig
  contig_length <- seq_lengths[contig_name]

  # Create a GRanges object specifying the entire contig
  which_range <- GRanges(seqnames = contig_name, ranges = IRanges(1, contig_length))

  # Create a ScanBamParam object with the specified contig
  param <- ScanBamParam(which = which_range)

  # Ensure the input BAM file is indexed
  if (!file.exists(paste0(bamfile, ".bai"))) {
    indexBam(bamfile)
  }

  # Filter the BAM file
  filterBam(file = bamfile, destination = output_bam, param = param)

  # Index the output BAM file
  indexBam(output_bam)

  # Return the path to the filtered BAM file
  return(output_bam)
}
