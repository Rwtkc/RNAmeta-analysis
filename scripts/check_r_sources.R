args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[[1]] else "."
files <- list.files(root, pattern = "\\.R$", full.names = TRUE, recursive = TRUE)
failures <- character()
for (file in files) {
  tryCatch(parse(file = file), error = function(e) {
    failures <<- c(failures, paste(file, conditionMessage(e)))
  })
}
if (length(failures)) stop(paste(failures, collapse = "\n"))
cat("Parsed", length(files), "R source files\n")
