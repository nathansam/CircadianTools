
LoadingGen <- function(genenumber) {
    loading_values <- rep(0, 20)
    for (i in 1:20) {
        loading_values[i] <- round(genenumber/20 * i)
    }
    return(loading_values)
}

LoadingPrint <- function(iteration, loading_values) {
    if (iteration == 1) {
        cat(crayon::red(noquote("This may take a while if using a large dataset! \n")))
        cat(noquote("Progress: \n"))
        cat(noquote(paste(replicate(20, "\u25A1"), collapse = "")))
        cat(noquote("\n"))
    }
    if (iteration %in% loading_values) {
        position <- match(iteration, loading_values)
        hashes <- paste(replicate(position, "\u25A0"), collapse = "")
        dashes <- paste(replicate(20 - position, "\u25A1"), collapse = "")
        cat(noquote(paste("\r", hashes, dashes, "\n", sep = "")))
        utils::flush.console()
    }
}
