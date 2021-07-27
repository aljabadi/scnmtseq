#' Merge tsv outputs from \code{\link{summarise}} functions into one file.
#'
#' These functions merge the tsv files in a given path into a single tsv file.
#'
#' @param anno_path path to the directory containing the summarised tsv files for each cell for a given annotation
#' @param context_path path to the directory containing the summarised tsv files for each annotation for a given context
#' @param force Logical, whether to overwrite even if the output file exists
#' @name merge
#' @return Writes merged tsv files to the same directory as the annotation/context files and returns TRUE
NULL
## suppress check warnings for undeclared symbols used in the function
utils::globalVariables("id", package = 'scnmtseq', add = FALSE)

#' @import data.table
#' @rdname merge
#' @export
merge_anno <- function(anno_path, force = FALSE) {
    ## for an annotation folder path containing the summarised samples files,
    ## create '{anno_path}.tsv' and merge all sample files

  # anno_path <- 'output/CpG/H3K27ac'
  anno_path <- gsub('/$', '', anno_path)
  merged_anno_path <- paste0(anno_path, '.tsv')
  files <- list.files(anno_path, full.names = TRUE)

  if (force | !file.exists(merged_anno_path))
  {
    for (file in files)
    {
      tsv <- fread(file) %>%  setnames(c("anno","sample","id","Nmet","N","rate"))
      tsv <- tsv[id != ''] # remove blank IDs, if any
      fwrite(tsv, merged_anno_path, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE, append = TRUE)
    }
  } else
  {
    file_exists(merged_anno_path)
  }

  TRUE
}

# merge_anno(anno_path = 'output/CpG/H3K27ac/')
#' @import data.table
#' @rdname merge
#' @export
merge_context <- function(context_path, force = FALSE, cores = 1) {
    context_path <- gsub('/$', '', context_path)
    ## for a context folder path containing the summarised anno files (.tsv) from merge_anno,
    ## create '{CpG/GpC}.tsv' and merge all annotation files

    anno_paths <- list.dirs(context_path, full.names = TRUE)[-1]

    if (cores > 1)
    {
      ## parallel backend
      BPPARAM <- if (.Platform$OS.type == "unix") MulticoreParam(workers = cores) else SnowParam(workers = cores)
      ## announce before a long parallel run with no serial outputs
      cat2("\nPATH: %s | merging annotations\n%s", context_path, paste0(anno_paths, collapse = "\n"))
    } else {
      BPPARAM <- SerialParam()
    }



    bplapply(seq_along(anno_paths), FUN = function(i, cores){
      anno_path <- anno_paths[i]
      if (cores == 1)
        cat2("\nPATH: %s | merging annotation %s / %s", context_path, i, length(anno_paths))
      merge_anno(anno_path = anno_path, force = force)
    }, BPPARAM = BPPARAM, cores = cores)

    merged_context_path <- paste0(context_path, '.tsv')

    files <- list.files(context_path, pattern = '.tsv$', full.names = TRUE)

    if (force | !file.exists(merged_context_path))
    {
      cat2("\nMerging all the annotations into", merged_context_path, "...\n")
      for (file in files)
      {
        tsv <- fread(file)
        fwrite(tsv, merged_context_path, quote=FALSE, sep="\t", col.names=FALSE, row.names=FALSE, append = TRUE)
      }
    } else
    {
      file_exists(merged_context_path)
    }

    TRUE
}
