#' This function calculates the quantile matrix for a given quantile q, using
#' the data_rain dataset.
#' The count_min parameter specifies the minimum number of observations
#' required to calculate a quantile.
#' The step parameter determines the resolution of the quantile matrix.
#'
#' @param q The quantile value.
#' @param data_rain The dataset containing the rainfall data.
#' @param count_min The minimum number of observations required to calculate a
#' quantile. Default is 80.
#' @param step The resolution of the quantile matrix. Default is 0.005.
#' @param zeros A logical value indicating whether to include zero values in the
#'              dataset. Default is TRUE.
#' @param qlim A logical value indicating a quantile limit, more precised some
#'             times. Default is TRUE.
#'
#' @return A matrix containing the quantile values for each pair of variables
#' in the dataset.
#'
#' @import evd
#' @import dplyr
#' @importFrom graphics par
#' @importFrom graphics abline
#'
#' @examples
#' quantile_matrix(0.5, data_rain)
#' quantile_matrix(0.9, data_rain, count_min = 100, step = 0.01)
#'
#' @export
quantile_matrix <- function(q, data_rain, count_min = 80, step = 0.005,
                            zeros = TRUE, qlim = TRUE) {
    N <- ncol(data_rain)
    quant_mat <- matrix(NA, nrow = N, ncol = N)
    count_excess <- matrix(NA, nrow = N, ncol = N)
    par(mfrow = c(3, 3))
    for (i in 1:(N)){
      for (j in (i):N) {
        # without NA values
        # if (i != j) {
        data <- drop_na(data_rain[c(i, j)]) # for all pairs
        if (!zeros) {
          data <- data[rowSums(data != 0, na.rm = TRUE) > 0, ]
        }
        n <- nrow(data)
        data <- cbind(rank(data[, 1]) / (n + 1), rank(data[, 2]) / (n + 1))
        colnames(data) <- c("s1", "s2")
        q_modif <- q
        # check excess above a threshold q
        cp_cond <- as.data.frame(data[data[, 2] > q, ])
        # nb of simultaneous excesses
        count <- sum(cp_cond[, 1] > q)
        # if there are not enough data above the quantile
        while (count < count_min) {
            q_modif <- q_modif - step # we reduce it
            cp_cond <- as.data.frame(data[data[, 2] > q_modif, ])
            count <- sum(cp_cond[, 1] > q_modif)
        }
        quant_mat[i, j] <- q_modif
        if (q != q_modif) {
            if (qlim) {
              chiplot(data, qlim = c(0.99, 0.999), xlim = c(0.99, 0.999),
                  which = 1) # to check on plot
            } else {
              chiplot(data, which = 1)
            }
            abline(v = q_modif)
        }
        count_excess[i, j] <- count
        # }
        }
    }
    return(list(quant_mat, count_excess))
}