#' tanalysis:
#' @description Experimental! A t.test is carried out on gene activity levels between time points and the number of significant increases & decreases is returned.
#' @param row.no The row number of the gene of interest
#' @param dataset A transcriptomics dataset. First columns should be gene names. All other columns should be expression levels.
#' @param psignificance The maximum p-value for which a result of a t-test is classed as significant
#' @return The number of signficant increases and decreases as a vector
#' @examples
#' ups.downs <- tanalysis(row.no = 1, dataset = Laurasmappings,  psignificance = 0.01)
#'
#' @export

tanalysis <- function(row.no, dataset, psignificance = 0.05) {
    timevector <- CircadianTools::maketimevector(dataset)  # vector of timevalues
    unique.timevector <- unique(timevector)  #vector of unique timevalues (ignores repetitions)
    gene <- dplyr::filter(dataset, dplyr::row_number() == row.no)  # select gene by row number
    gene <- gene[-1]  #remove sample column
    gene <- as.matrix(t(gene))  # transpose gene to column
    genedf <- data.frame(gene, timevector)  # create dataframe for gene activity and timevector
    names(genedf) <- c("activity", "timevector")
    n.ups <- 0  # Set variable to list the number of significant positive changes found via the t-tests.
    n.downs <- 0  # Set variable to list the number of significant negative changes found via the t-tests.

    for (i in 1:(length(unique.timevector) - 1)) {
        group1 <- subset(genedf, timevector == unique.timevector[i])
        group2 <- subset(genedf, timevector == unique.timevector[i + 1])

        zerogroup <- FALSE  # Flag for if both groups entirely consists of zero values. Assume false
        if (all(group1[, 1] == 0)) {
            # Is group 1 entirely zero values?  Is group 2 entirely zero values?
            if (all(group2[, 1] == 0)) {
                zerogroup <- TRUE  # Both groups are entirely zero values. Don't attempt t-test
            }
        }

        if (zerogroup == FALSE) {
            t.test.obj <- t.test(group1[, 1], group2[, 1], var.equal = TRUE)  # Carry out t-test assuming equal variance
            if (t.test.obj$p.value < psignificance) {
                # If p-value from t-test is below significance threshold Is the t statistic positive?
                if (t.test.obj$statistic > 0) {
                  n.ups <- n.ups + 1
                } else {
                  # Must be negative then
                  n.downs <- n.downs + 1
                }
            }
        }
    }

    return(c(n.ups, n.downs))  # Return the number of signficant increases and decreases
}
