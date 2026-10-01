#' marineSim: the stRF simulation engine
#'
#' Functions are grouped as follows.
#' \itemize{
#'   \item World: \code{\link{matern_grf}}, \code{\link{make_world}},
#'     \code{\link{world_grid}}, \code{\link{cell_xy}}, \code{\link{covariates}}.
#'   \item Designs and datasets: \code{\link{sample_sites}},
#'     \code{\link{build_dataset}}, \code{\link{dataset_seed}}.
#'   \item Feature blocks: \code{\link{spatial_lags}}, \code{\link{edf}},
#'     \code{\link{nearest_train}}, \code{\link{dist_km}}.
#'   \item Models and scoring: \code{\link{score}}, \code{\link{glm_score}},
#'     \code{\link{auc}}.
#'   \item Attribution: \code{\link{required_subsets}}, \code{\link{skey}},
#'     \code{\link{shapley_row}}, \code{\link{add_shapley}}, \code{\link{boot_ci}},
#'     \code{\link{decomp_table}}.
#'   \item Extrapolation diagnostic: \code{\link{area_of_applicability}}.
#'   \item Constants: \code{\link{constants}}.
#' }
#' The analysis scripts in the repository's \code{analysis/} directory drive the
#' full benchmark, the auxiliary studies and the seagrass case study.
#'
#' @keywords internal
"_PACKAGE"
