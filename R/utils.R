#' Utilities
#'
#' Utility functions
#'
#' @param fmt,... passed to \code{\link[base]{sprintf}}
#' @param filepath,path path to file
#' @name utils
#' @return None
NULL
## ---- colourised 'cat's
#' @examples
#' cat2("this is a '%s' text", 'yellow')
#' @rdname utils
#' @export
cat2 <- function(fmt, ...)
{ # a yellow cat
    cat("\033[38;5;214m", sprintf(fmt, ...), "\033[39m")
}

#' @examples
#' cat3("this is a '%s' text", 'purple')
#' @rdname utils
#' @export
cat3 <- function(fmt, ...)
{ # a purple cat
    cat("\033[35m", sprintf(fmt, ...), "\033[39m")
}

#' @examples
#' cat4("this is a '%s' text", 'blue')
#' @rdname utils
#' @export
cat4 <- function(fmt, ...)
{ # a blue cat
    cat("\033[34m", sprintf(fmt, ...), "\033[39m")
}
## ---- END colourised 'cat's
#' @examples
#' rm_trailing("foo/bar/")
#' # "foo/bar"
#' @rdname utils
#' @export
## remove possible trailing slash
rm_trailing <- function(path)
{
    gsub('/$', '', path)
}

#' @rdname utils
#' @export
file_exists <- function(filepath)
{
    # cat that the file exists and run was skipped
    # file_exists(LETTERS[1:3])
    cat2("\nFile(s) exist(s):\n")
    cat3(paste(filepath, collapse = '\n'), sep = '')
    cat2("\nUse --force 'TRUE' to overwrite.\n\n")
}

#' @rdname utils
#' @export
check_file_exists <- function(filepath)
{
    # if filepath does not exist, throw error
    if (!file.exists(filepath))
        stop(filepath, " does not exist\n")
    invisible(NULL)
}
