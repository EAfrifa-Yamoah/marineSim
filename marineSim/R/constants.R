#' Simulation constants (Section 3.1 of the article)
#'
#' Fixed design quantities of the ground truth generator and the forest backbone.
#' They are exported so that analysis scripts can refer to them by name; changing
#' them changes the random stream and therefore every reported number.
#'
#' \describe{
#'   \item{GRID}{Number of cells along each side of the square domain (50).}
#'   \item{EXTENT_KM}{Side length of the domain in kilometres (400).}
#'   \item{CELL_KM}{Cell size in kilometres (8).}
#'   \item{N_TEST}{Number of held out test cells per dataset (1,500).}
#'   \item{RF_TREES}{Trees per forest in the benchmark (150).}
#'   \item{N_ANCHORS}{Number of k-means anchors for the distance field block (15).}
#'   \item{RADII_KM}{Radii of the multi scale spatial lags in kilometres (25, 75, 150).}
#'   \item{BASE_SEED}{Base seed from which world and dataset seeds are derived.}
#'   \item{C0_COUNT, ETA_SCALE}{Intercept (log 3) and slope (0.5) of the abundance
#'     link, \eqn{\mu = \mathrm{mask} \cdot \exp(C_0 + \eta / 2)}.}
#'   \item{BLOCKS}{Names of the four engineered feature blocks: COORD, EDF, LAG, TEMP.}
#' }
#' @name constants
#' @aliases GRID EXTENT_KM CELL_KM N_TEST RF_TREES N_ANCHORS RADII_KM BASE_SEED C0_COUNT ETA_SCALE BLOCKS
NULL

#' @rdname constants
#' @export
GRID <- 50L
#' @rdname constants
#' @export
EXTENT_KM <- 400
#' @rdname constants
#' @export
CELL_KM <- EXTENT_KM / GRID
#' @rdname constants
#' @export
N_TEST <- 1500L
#' @rdname constants
#' @export
RF_TREES <- 150L
#' @rdname constants
#' @export
N_ANCHORS <- 15L
#' @rdname constants
#' @export
RADII_KM <- c(25, 75, 150)
#' @rdname constants
#' @export
BASE_SEED <- 20260601L
#' @rdname constants
#' @export
C0_COUNT <- log(3)
#' @rdname constants
#' @export
ETA_SCALE <- 0.5
#' @rdname constants
#' @export
BLOCKS <- c("COORD", "EDF", "LAG", "TEMP")
