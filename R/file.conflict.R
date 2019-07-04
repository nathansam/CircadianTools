#' FileConflict
#' @description Checks if a file which will be created already exists and, if necessary asks the user if this file should be overwritten.
#' @param filename The name of a file (with extension included)
#' @export

FileConflict <- function(filename){

  if (file.exists(filename) == TRUE) { # If the file already exists, should it be overwritten?
    prompt <- cat(paste("The file, ", filename, ", already exists. Do you wish to overwrite the file? [Y/n]  \n"))
    continue <- readline(prompt = prompt) # Ask the user

    if (continue == "y" | continue == "Y" | continue == "yes" | continue == "Yes") {
      cat("Overwriting File")
      file.remove(filename) # Overwrite the file if yes
    } else {
      cat("Command Terminated")
      opt <- options(show.error.messages = FALSE)
      on.exit(options(opt)) # Stop the function
      stop()
    }
  }
}
