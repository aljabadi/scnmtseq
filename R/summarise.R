########################################################################################### -
##  Script to quantify DNA methylation and chromatin accessibility over genomic features ##
########################################################################################### -

# This script overlaps the individual CpG sites with genomic features such as promoters, gene bodies, etc.

## Input ##

# (1) For every cell, a long data.table with (at least ) columns ["chr", "pos", "rate"]
# Example:
#   chr     pos      rate
#   19    3152031     1
#   19    3152424     0

# IMPORTANT: For every CpG/GpC site, the the rate must be 0 or 1

# (2) genomic feature annotation files in BED6 format
#   chr start end strand id anno
#   1	3531624	3531843	*	CGI_1	CGI
#   1	3670619	3671074	*	CGI_2	CGI
#   1	3671654	3672156	*	CGI_3	CGI

#' Quantify DNA methylation and chromatin accessibility over genomic features
#'
#' These functions overlap the individual CpG sites with genomic features such as
#' promoters, gene bodies, etc.
#'
#' @param filepath path to cov.gz file
#' @param annotation_file path to the annotation file that has \code{c("chr","start","end","strand","id")} columns
#' @param anno_name Character, annotation name to use and add to the summarised table, e.g. 'Enhancer'
#' @param context One of c('CG', 'GC')
#' @param valid_chromosomes Character vector of chromosomes to keep. Must match the chromosome notation in the data
#' @param a,b Beta prior parameters. See SN 1 \url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4117646/}
#' @param out.dir Output directory to save the tsv file. The output file will be of form: outdir/CpG_or_GpC/anno_name/cov_file_name.tsv
#' @param force Logical, whether to overwrite even if the output file exists
#' @name summarise
#' @return Writes file to out.dir and returns NULL
NULL

## suppress R CMD check complain for NSE
utils::globalVariables( c("N", "Nmet", "anno", "chr", "end", "id", "rate", "set_names", "start", "."))

#' @details
#' \itemize{
#' \item{\code{summarise_sample} summarises the calls for a given cov.gz file and annotation file}
#' \item{\code{summarise_anno} summarises the calls for a given annotation file for all cells}
#' \item{\code{summarise_context} summarises the calls for a all annotation file for all cells for a given context (CG or GC)}
#' \item{\code{summarise} summarises the calls for a all annotation file for all cells for both contexts (CG and GC)}

#' }

#' @import data.table
#' @import stringr
#' @import BiocParallel
#' @import parallel
#' @rdname summarise
#' @export
summarise_sample <- function(filepath,
                             annotation_file,
                             anno_name, # name to use
                             context = c('CG', 'GC') ,
                             valid_chromosomes = as.character(c(1:22, 'MT', 'X', 'Y')),
                             a = 0,
                             b = 0,
                             out.dir,
                             force=FALSE)
{
  # fname.out = create_output_filename(filepath, anno_name)
  # filepath <- 'data/CpG/10A_DAC_060320_cytosine.NOMe.CpG.cov.gz'
  out.dir <- gsub("/$", "", out.dir) # remove trailing slash, if any
  # sample name as well_condition_batch
  sample.name <- str_split(basename(filepath), pattern = '_cytosine')[[1]][1]
  # output file same as cov file with tsv extension
  fname.out.base <- gsub('.cov.gz$', '.tsv', basename(filepath))
  fname.out <- file.path(out.dir, ifelse(context == 'CG', 'CpG', 'GpC'), anno_name, fname.out.base)
  if (!dir.exists(dirname(fname.out)))
    dir.create(dirname(fname.out), recursive = TRUE, showWarnings = FALSE)

  if (!file.exists(fname.out) | force)
  {
    context <- toupper(context)
    context <- match.arg(context)
    sample_data <- fread(filepath)
    sample_data <- setnames(sample_data, c('chr', 'start', 'end', 'rate',
                                           'ns', # n_support
                                           'nn'  #n_notSupport
    ))
    sample_data <- sample_data[chr %in% valid_chromosomes]
    # rate should be 0 or 1 at cytosine level
    sample_data <- sample_data[, rate :=rate/100]
    sample_data <- sample_data[, rate := round(rate)]
    # create a factor from chromosomes
    sample_data <- sample_data[,chr := factor(chr)]
    # sort (key) for faster search when overlapping
    sample_data <- setkey(x = sample_data, chr, start, end)
    # head(file)
    # unique(file$chr)
    # length(unique(file$chr))

    ## ---- annotation
    ## bigBed
    # annotation_file <- 'data/Annotations/HL60_H3K27ac_ENCFF828ATQ.bigBed'
    # annotation_file <- 'data/Annotations/CGI_SeqMonk.txt'
    # annotation_file <- 'data/Annotations/bed/Exons.bed'
    anno_data <- fread(annotation_file)
    colnames(anno_data) <- c("chr","start","end","strand","id")
    anno_data <- anno_data[id != ''] # remove blank IDs
    anno_data[,anno := anno_name]
    anno_data <- setkey(x = anno_data, chr, start, end)
    ov <- foverlaps(sample_data, anno_data, nomatch=0)
    ov <- ov[,"i.end":=NULL] %>% setnames("i.start","pos")
    # .[,c("sample","anno") := list(sample,anno)] %>%
    # Compute number of methylated CpGs and the corresponding methylation rates
    ov <- ov[,.(Nmet=sum(rate==1), N=.N), keyby=.(id)]
    # a and b correspond to the beta prior - see Smallwood et al 2014 SN1
    ov <- ov[,rate:=(Nmet + a)/(N + b)]
    ov <- ov[, `:=`(sample = sample.name, anno = anno_name)]
    ov <- ov[, c("anno","sample","id","Nmet","N","rate")]
    # Store and save results
    if (any(is.na(ov$rate)))
      stop("NA found in rates: ", filepath, " with anno: ", annotation_file)
    fwrite(ov, fname.out, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE)

  } else
  {
    file_exists(fname.out)
  }

  return(TRUE)
}

#' @param samplefiles.path path to all cov.gz files for a given context
#' @param cores Integer, number of CPUS to use
#' @rdname summarise
#' @export
summarise_anno <-
  function(samplefiles.path,
           annotation_file,
           anno_name,
           context = c('CG', 'GC'),
           valid_chromosomes = as.character(c(1:22, 'MT', 'X', 'Y')),
           a = 0,
           b = 0,
           out.dir,
           cores = parallel::detectCores(),
           force = FALSE)
  {

    filepaths <- list.files(samplefiles.path, pattern = '.cov.gz$', full.names = TRUE)
    mget(names(formals()), sys.frame(sys.nframe()))
    if (cores > 1)
    {
      ## parallel backend
      BPPARAM <- if (.Platform$OS.type == "unix") MulticoreParam(workers = cores) else SnowParam(workers = cores)
    } else {
      BPPARAM <- SerialParam()
    }
    context = match.arg(context)
    BiocParallel::bplapply(filepaths, function(filepath)
        {
        summarise_sample(
            filepath = filepath,
            anno_name = anno_name,
            annotation_file = annotation_file,
            context = context,
            valid_chromosomes = valid_chromosomes,
            a = a,
            b = b,
            out.dir = out.dir,
            force = force
        )
        }
    , BPPARAM = BPPARAM)
    BiocParallel::bpstop(BPPARAM)

    return(invisible(NULL))
}

# summarise_anno(samplefiles.path = samplefiles.path, annotation_file = annotation_file)

#' @rdname summarise
#' @param annotation_files Named character vector (name = 'path') where names will be used as annotation names
#' @export
summarise_context <- function(samplefiles.path,
                              annotation_files, # named character (name = 'path') where names will be used as assay names
                              context = c('CG', 'GC'),
                              valid_chromosomes = as.character(c(1:22, 'MT', 'X', 'Y')),
                              a = 0,
                              b = 0,
                              out.dir,
                              cores = parallel::detectCores(),
                              force = FALSE)
{
    for (i in seq_along(annotation_files))
    {
        annotation_file <- annotation_files[i]
        anno_name <- names(annotation_files)[i]
        cat2("\nPROCESSING CONTEXT: %s | ANNOTATION: %s (%s / %s)", context, anno_name, i, length(annotation_files))
        summarise_anno(
            samplefiles.path = samplefiles.path,
            annotation_file = annotation_file,
            anno_name = anno_name,
            context = match.arg(context),
            valid_chromosomes = valid_chromosomes,
            a = a,
            b = b,
            out.dir = out.dir,
            cores = cores,
            force = force
        )
    }
  return(invisible(NULL))
}

#' @param datapath Path to the data. It includes datapath/CpG and datapath/GpC subdirectories which contain the cov files for all cells
#' @rdname summarise
#' @export
summarise <-
  function(datapath = 'data',
           annotation_files,
           valid_chromosomes = as.character(c(1:22, 'MT', 'X', 'Y')),
           a = 0,
           b = 0,
           force,
           cores = parallel::detectCores(),
           out.dir = 'output/scnomeseq') {
    # This function wraps \code{\link{summarise_context}} for both CG and GC contexts
    datapath <- rm_trailing(datapath)

    # GpC
    GpC.path <- file.path(datapath, 'GpC')
    gc()
    summarise_context(
      samplefiles.path = GpC.path,
      annotation_files = annotation_files,
      context = 'GC',
      valid_chromosomes = valid_chromosomes,
      a = a,
      b = b,
      out.dir = out.dir,
      force = force,
      cores = cores
    )

    # CpG
    CpG.path <- file.path(datapath, 'CpG')
    gc()
    summarise_context(
      samplefiles.path = CpG.path,
      annotation_files = annotation_files,
      context = 'CG',
      valid_chromosomes = valid_chromosomes,
      a = a,
      b = b,
      out.dir = out.dir,
      force = force,
      cores = cores
    )

    return(invisible(NULL))
  }

