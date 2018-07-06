#' @name PanelSurv
#' @rdname PanelSurv
#' @title Create a PanelSurv Object
#'
#' @description Create a panel count survival object,
#' usually used as a response variable in a model formula.
#'
#' @param ID Observation subject's ID.
#' @param time Observation time.
#' @param count Observation subject's ID.
#' @param x An \code{PanelSurv} object.
#'
#' @return An object of S3 class \code{"PanelSurv"}.
#' \describe{
#'   \item{psDF}{a data frame, part of original input data frame with variable
#' "ID", "time" and "count".}
#'   \item{timeGrid}{ordered distinct observation times in the set of all
#'     observation times.}
#'   \item{panelMatrix}{a matrix representation of panel count data, one
#'     row per subject, one column per time point in \code{"timeGrid"}.}
#' }
#'
#'   In the case of \code{is.PanelSurv}, a logical value \code{TRUE} if
#'   \code{x} inherits from class \code{"PanelSurv"}, otherwise an \code{FALSE}.
#'
#'   In the case of \code{plot.PanelSurv}, a tile plot of
#'   \code{panelMatrix} produced by package \code{ggplot2} with color
#' indicating number of counts since last observation time.
#'
#' @seealso \code{\link{panelReg}}
#' @examples
#' data(blaTum)
#' response <- with(blaTum, PanelSurv(id, time, count))
#' is.PanelSurv(response)
#' plot(response)
NULL

#' @rdname PanelSurv
#' @export
PanelSurv <- function(ID, time, count) {
    if (sum(time <= 0) > 0)
        stop("Observation time must be positive.")
    index <- which(!duplicated(ID))
    N <- length(index)
    uniqueID <- ID[index]
    timeGrid <- sort(unique(time))
    panelMatrix <- matrix(NA, N, length(timeGrid))
    for (i in 1:N) {
        rowSet <- which(ID == uniqueID[i])
        panelMatrix[i, which(timeGrid %in% time[rowSet])] <- count[rowSet]
    }
    ps <- list(psDF=data.frame(ID = ID, time=time, count=count),
               timeGrid=timeGrid, panelMatrix=panelMatrix)
    class(ps) <- "PanelSurv"
    ps
}

#' @rdname PanelSurv
#' @export
is.PanelSurv <- function(x) inherits(x, "PanelSurv")

#' @export
plot.PanelSurv <- function(x, ...) {
    with(x$psDF, ggplot(data = x$psDF, aes(time, ID, height = 2, width = 15)) +
                 geom_tile(aes(fill = count)) + theme_bw() +
                 scale_fill_gradient(low = "grey", high = "black"))
}