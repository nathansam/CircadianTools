#' CosinorAnalysis:
#' @description Fits cosinor models to transcriptomics data and plots the best-fitting models using ggplot2.
#'
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @param period The period of rhythmicity which is being tested for. Defaults to 24 (circadian).
#' @param timelag The value of time for which observations begin. Defaults to 6.
#' @param threshold Set significance value for which genes are considered significant and thus plotted. Defaults to 6e-07
#' @param save Logical. If TRUE, saves plots to working directory. Defaults to FALSE.
#' @param print Logical. If TRUE renders significant genes in the plot viewer. Defaults to TRUE
#' @param df Logical. If TRUE a dataframe containing the results of the cosinor analysis will be returned. Defaults to TRUE.
#' @param adj String. p-value adjustment. Defaults to 'bonferonni'. 'none is also supported'
#' @return Prints or saves ggplot2 object(s). Optionally returns dataframe containing gene name and p values from F-test ranking of cosinor models
#' @examples
#' cosinor_results <- CosinorAnalysis(Laurasmappings)
#'
#' @export


CosinorAnalysis <- function(dataset, period = 24, timelag = 6, threshold = 0.05, adj = "bonferroni", save = FALSE,
                            print = TRUE, df = TRUE) {
  expected_adj <- c("bonferroni", "Bonferroni", "none")

  if (adj %in% expected_adj == FALSE) {
    stop(paste("The adjustment method ", adj, " is not recognized"))
  }


  #dataset <- CircadianTools::GeneClean(dataset)
  genenumber <- nrow(dataset)  #number of genes in the dataset
  if (adj == "bonferroni" | adj == "Bonferroni") {
    threshold <- threshold/genenumber
  }
  pvalues <- rep(0, genenumber)  #init list of pvalues
  cosinor.pvalue.df <- data.frame(sample = dplyr::select(dataset, 1), pVal = pvalues)  #first column gene name, second:pvalue
  timevector <- CircadianTools::MakeTimevector(dataset)
  loading_values <- CircadianTools::LoadingGen(genenumber)

  for (i in 1:genenumber) {
    CircadianTools::LoadingPrint(i, loading_values)

    genematrix <- dplyr::filter(dataset, dplyr::row_number() == i)

    genename <- genematrix[1, 1]
    genematrix <- genematrix[-1]
    genematrix <- t(genematrix)
    geneexpression <- data.frame(timevector - timelag, genematrix)
    names(geneexpression) <- c("timevector", "activity")
    cosinormodel <- cosinor::cosinor.lm(activity ~ time(timevector), period = period, data = geneexpression)
    cosinor.pvalue.df[i, 2] <- cosinor2::cosinor.detect(cosinormodel)[4]
    if (cosinor2::cosinor.detect(cosinormodel)[4] < threshold) {
      if (adj == "bonferroni" | adj == "Bonferroni") {
        plot_title <- paste("Gene=", genename, ", P-value=", round(cosinor2::cosinor.detect(cosinormodel)[4] *
                                                                     genenumber, 10))
      }

      if (adj == "none") {
        plot_title <- paste("Gene=", genename, ", P-value=", round(cosinor2::cosinor.detect(cosinormodel)[4],
                                                                   10))
      }

      cosinorplot <- CircadianTools::ggplot.cosinor.lm(cosinormodel, endtime = tail(timevector, n = 1) - timelag) +
        ggplot2::geom_point(ggplot2::aes(y = activity, x = timevector), data = geneexpression, size = 3,
                            alpha = 0.5, color = "#39A5AE") + ggplot2::ggtitle(plot_title) + ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 1)) +
        ggplot2::theme(text = ggplot2::element_text(size = 12)) + ggplot2::xlab("Time (hours)") + ggplot2::ylab("Trancripts Per Million (TPM)")


      if (save == TRUE) {
        ggplot2::ggsave(paste("Cosinor_", genename, ".png"), cosinorplot, width = 10, height = 4.5, units = "in")
      }
      if (print == TRUE) {
        print(cosinorplot)
      }
    }
  }
  if (df == TRUE) {
    return(cosinor.pvalue.df)
  }
}


#' CosinorAnalysisPar:
#' @description Parallel Implementation of \link{CosinorAnalysis}. Fits cosinor models to transcriptomics data and plots the best-fitting models using ggplot2.
#'
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @param period The period of rhythmicity which is being tested for. Defaults to 24 (circadian).
#' @param timelag The value of time for which observations begin. Defaults to 6.
#' @param threshold Set significance value for which genes are considered significant and thus plotted. Defaults to 6e-07
#' @param save Logical. If TRUE, saves plots to working directory. Defaults to FALSE.
#' @param print Logical. If TRUE renders significant genes in the plot viewer. Defaults to TRUE
#' @param df Logical. If TRUE a dataframe containing the results of the cosinor analysis will be returned. Defaults to TRUE.
#' @return Prints or saves ggplot2 object(s). Optionally returns dataframe containing gene name and p values from F-test ranking of cosinor models
#' @examples
#' cosinor.results <- CosinorAnalysisPar(Laurasmappings)
#'
#' @export


CosinorAnalysisPar <- function(dataset, period = 24, nthreads = NULL, timelag = 6) {
  if (is.null(nthreads) == TRUE) {
    nthreads <- parallel::detectCores()
  }
  `%dopar%` <- foreach::`%dopar%` # Load the dopar binary operator from foreach package
  if (is.null(nthreads) == TRUE) {
    nthreads <- parallel::detectCores()
  }


  dataset <- CircadianTools::GeneClean(dataset)
  genenumber <- nrow(dataset)  #number of genes in the dataset
  pvalues <- rep(0, genenumber)  #init list of pvalues
  cosinor.pvalue.df <- data.frame(sample = dplyr::select(dataset, 1), pVal = pvalues)  #first column gene name, second:pvalue
  timevector <- CircadianTools::MakeTimevector(dataset)
  loading_values <- CircadianTools::LoadingGen(genenumber)


  cl <- parallel::makeForkCluster(nthreads)
  doParallel::registerDoParallel(cl)
  cosinor.pvalue.df <- foreach::foreach(i = 1:genenumber, .combine = rbind) %dopar% {
    CircadianTools::LoadingPrint(i, loading_values)

    genematrix <- dplyr::filter(dataset, dplyr::row_number() == i)
    sample <- genematrix[1, 1]
    genematrix <- genematrix[-1]
    genematrix <- t(genematrix)
    geneexpression <- data.frame(timevector - timelag, genematrix)
    names(geneexpression) <- c("timevector", "activity")
    cosinormodel <- cosinor::cosinor.lm(activity ~ time(timevector), period = period, data = geneexpression)
    pVal <- cosinor2::cosinor.detect(cosinormodel)[4]
    data.frame(sample, pVal)
  }

  return(cosinor.pvalue.df)
}


#' CosinorPlot:
#' @description Fits a cosinor model to a given gene in a given dataset and plots the model
#'
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @param genename The name of a gene intended for plotting. Must be a string.
#' @param timelag The value of time for which observations begin. Defaults to 6.
#' @param save Logical. If TRUE, saves plot to working directory. Defaults to FALSE.
#' @param print Logical. If TRUE renders plot in the plot viewer. Defaults to TRUE
#' @return Prints or saves ggplot2 object(s)
#' @examples
#' CosinorPlot('comp100000_c0_seq2', Laurasmappings )
#'
#' @export


CosinorPlot <- function(genename, dataset, timelag = 6, period = 24, print = TRUE, save = FALSE) {
  genematrix <- subset(dataset, dataset[1] == genename)
  timevector <- CircadianTools::MakeTimevector(genematrix)  # makes timevector

  genematrix <- genematrix[-1]  #remove gene name from subset
  genematrix <- t(genematrix)
  geneexpression <- data.frame(timevector - timelag, genematrix)
  names(geneexpression) <- c("timevector", "activity")
  cosinormodel <- cosinor::cosinor.lm(activity ~ time(timevector), period = period, data = geneexpression)
  cosinor.plot.object <- CircadianTools::ggplot.cosinor.lm(cosinormodel, endtime = 21) + ggplot2::geom_point(ggplot2::aes(y = activity,
                                                                                                                          x = timevector), data = geneexpression, size = 3, alpha = 0.5, color = "#008dd5") + ggplot2::ggtitle(paste("Gene=",
                                                                                                                                                                                                                                     genename, ", P-value=", round(cosinor2::cosinor.detect(cosinormodel)[4], 10))) + ggplot2::theme_bw() + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 1)) +
    ggplot2::theme(text = ggplot2::element_text(size = 12)) + ggplot2::xlab("Time (hours)") + ggplot2::ylab("Trancripts Per Million (TPM)")

  if (save == TRUE) {
    ggplot2::ggsave(paste("Cosinor_", genename, ".png"), cosinor.plot.object, width = 10, height = 4.5, units = "in")
  }
  if (print == TRUE) {
    return(cosinor.plot.object)
  }
}


#' CosinorSignificantPlot :
#' @description Prints or saves the genes found to be most significant by \link{CosinorAnalysis}
#'
#' @param results A dataframe generated by \code{CosinorAnalysis}
#' @param dataset A transcriptomics dataset which was used with \code{CosinorAnalysis} to generate the results dataframe. First columns should be gene names. All other columns should be expression levels.
#' @param number The number of most significant genes printed or saved
#' @param period The period of rhythmicity which is being tested for. Defaults to 24 (circadian).
#' @param save Logical. If TRUE, saves plots to working directory. Defaults to FALSE.
#' @param print Logical. If TRUE renders significant genes in the plot viewer. Defaults to TRUE
#' @return Prints or saves ggplot2 object(s).
#' @examples
#' cosinor.results <- CosinorAnalysis(Laurasmappings)
#' CosinorSignificantPlot(cosinor.results, Laurasmappings, number = 15, period=24 ,save=TRUE)
#'
#' @export
CosinorSignificantPlot <- function(results, dataset, number = 10, period = 24, print = TRUE, save = FALSE) {
  results <- results[order(results$pVal), ]  # Order by most significant p-value

  for (i in 1:number) {
    myplot <- CircadianTools::CosinorPlot(as.character(results[i, 1]), dataset, period = period)
    myplot <- myplot + ggplot2::ggtitle(paste("Gene = ", as.character(results[i, 1]), " P-Value = ", as.character(results[i,
                                                                                                                          2])))

    if (print == TRUE) {
      print(myplot)
    }

    if (save == TRUE) {
      ggplot2::ggsave(paste("rank=", i, "Cosinor_", as.character(results[i, 1]), ".png"), myplot, width = 10,
                      height = 4.5, units = "in")
    }
  }
}

