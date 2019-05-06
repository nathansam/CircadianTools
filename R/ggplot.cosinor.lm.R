

#' Plot a cosinor model
#'
#' Adapted from the Cosinor package by Michael Sachs. Given a cosinor.lm model fit, generate a plot of the data with the fitted values.
#' Optionally allows for plotting by covariate levels 0 and 1. Unlike the original version, the function will plot for the full timecourse rather just one full period.
#'
#'
#' @param object An object of class \code{cosinor.lm}
#' @param x_str Character vector naming the covariate(s) to be plotted. May be NULL to plot overall curve
#' @param endtime The last time value for the time course
#'
#'
#'
#'
#'
#'
ggplot.cosinor.lm <- function(object, x_str = NULL, endtime) {
    
    timeax <- seq(0, endtime, length.out = 200)
    covars <- grep("(rrr|sss)", attr(object$fit$terms, "term.labels"), invert = TRUE, 
        value = TRUE)
    
    newdata <- data.frame(time = timeax, rrr = cos(2 * pi * timeax/object$period), 
        sss = sin(2 * pi * timeax/object$period))
    for (j in covars) {
        newdata[, j] <- 0
    }
    if (!is.null(x_str)) {
        
        for (d in x_str) {
            
            tdat <- newdata
            tdat[, d] <- 1
            newdata <- rbind(newdata, tdat)
            
        }
        newdata$levels <- ""
        for (d in x_str) {
            
            newdata$levels <- paste(newdata$levels, paste(d, "=", newdata[, 
                d]))
            
        }
        
        
    }
    
    
    newdata$Y.hat <- predict(object$fit, newdata = newdata)
    
    if (missing(x_str) || is.null(x_str)) {
        
        ggplot2::ggplot(newdata, ggplot2::aes_string(x = "time", y = "Y.hat")) + 
            ggplot2::geom_line()
        
    } else {
        
        ggplot2:ggplot(newdata, aes_string(x = "time", y = "Y.hat", col = "levels")) + 
            ggplot2::geom_line()
        
    }
}
