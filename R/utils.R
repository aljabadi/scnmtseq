## ---- colourised 'cat's
cat2 <- function(fmt, ...)
{ # a yellow cat
    cat("\033[38;5;214m", sprintf(fmt, ...), "\033[39m")
}

# cat2("this is a '%s' text", 'yellow')

cat3 <- function(fmt, ...)
{ # a purple cat
    cat("\033[35m", sprintf(fmt, ...), "\033[39m")
}

# cat3("this is a '%s' text", 'purple')

cat4 <- function(fmt, ...)
{ # a blue cat
    cat("\033[34m", sprintf(fmt, ...), "\033[39m")
}

# cat4("this is a '%s' text", 'blue')
## ---- END colourised 'cat's

## remove possible trailing slash
rm_trailing <- function(path)
{
    gsub('/$', '', path)
}

file_exists <- function(filepath)
{
    # cat that the file exists and run was skipped
    # file_exists(LETTERS[1:3])
    cat2("\nFile(s) exist(s):\n")
    cat3(paste(filepath, collapse = '\n'), sep = '')
    cat2("\nUse --force 'TRUE' to overwrite.\n\n")
}

check_file_exists <- function(filepath)
{
    if (!file.exists(filepath))
        stop(filepath, " does not exist\n")
    invisible(NULL)
}
