load_stored_mcmc <- function(path = NULL) {
  candidate_paths <- c(
    path,
    "data/storedMCMCiter",
    "app/data/storedMCMCiter"
  )
  
  candidate_paths <- candidate_paths[!is.na(candidate_paths)]
  
  existing_path <- candidate_paths[dir.exists(candidate_paths)][1]
  
  if (is.na(existing_path)) {
    stop(
      "No stored MCMC directory found. Tried:\n",
      paste(candidate_paths, collapse = "\n")
    )
  }
  
  files <- list.files(existing_path, pattern = "\\.rds$", full.names = TRUE)
  
  if (length(files) == 0) {
    stop("No stored MCMC files found in: ", existing_path)
  }
  
  out <- lapply(files, readRDS)
  names(out) <- sub("\\.rds$", "", basename(files))
  
  out
}

result.CAPPMx.base <- load_stored_mcmc()