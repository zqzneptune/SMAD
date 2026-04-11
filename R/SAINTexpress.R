#' SAINTexpress-spc
#' SAINTexpress scoring using spectral count data (spc mode).
#' Implements the SAINT-MRF model with Generalized Poisson distribution
#' for spectral count-based AP-MS data.
#'
#' @title SAINTexpress-spc
#' @param inter_df A data frame with columns: ip_id, bait_id, prey_id, quant.
#'   Each row represents one interaction observation (spectral count).
#' @param prey_df A data frame with columns: prey_id.
#'   Optional columns: prey_length, gene_name.
#' @param bait_df A data frame with columns: ip_id, bait_id, test_ctrl.
#'   test_ctrl must be "T" (test) or "C" (control).
#' @param GO_df Optional data frame with columns: go_term, genes.
#'   genes should be space-separated gene names.
#' @param f Frequency parameter (0 < f < 1). Default 0.5.
#' @param R Maximum number of replicates for scoring. Default 100.
#' @param L Maximum number of controls per prey. Default 100.
#'
#' @return A data frame with columns: Bait, Prey, PreyGene, Spec, SpecSum,
#'   AvgSpec, NumReplicates, ctrlCounts, AvgP, MaxP, TopoAvgP, TopoMaxP,
#'   SaintScore, logOddsScore, FoldChange, BFDR, boosted_by.
#'
#' @author Hyungwon Choi, Damian Fermin (original SAINTexpress),
#'   Qingzhou Zhang (Rcpp integration)
#' @references Teo G, et al. SAINTexpress: improvements and additional
#'   features in Significance Analysis of INTeractome software.
#'   J Proteomics. 2014;100:37-43.
#'
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib SMAD
#' @export
#' @examples
#' bait_path <- system.file("exdata", "TIP49", "bait.dat", package = "SMAD")
#' prey_path <- system.file("exdata", "TIP49", "prey.dat", package = "SMAD")
#' inter_path <- system.file("exdata", "TIP49", "inter.dat", package = "SMAD")
#' 
#' bait <- read.table(bait_path, sep = "\t", header = FALSE, 
#'                    col.names = c("ip_id", "bait_id", "test_ctrl"))
#' prey <- read.table(prey_path, sep = "\t", header = FALSE, 
#'                    col.names = c("prey_id", "prey_length"))
#' inter <- read.table(inter_path, sep = "\t", header = FALSE, 
#'                     col.names = c("ip_id", "bait_id", "prey_id", "quant"))
#' 
#' result <- SAINTexpress_spc(inter, prey, bait)
#' head(result)

SAINTexpress_spc <- function(inter_df, prey_df, bait_df,
                              GO_df = NULL,
                              f = 0, R = 100L, L = 100L) {
    # Input validation
    if(!is.data.frame(inter_df))
        stop("inter_df must be a data.frame")
    if(!is.data.frame(prey_df))
        stop("prey_df must be a data.frame")
    if(!is.data.frame(bait_df))
        stop("bait_df must be a data.frame")

    inter_cols <- c("ip_id", "bait_id", "prey_id", "quant")
    if(!all(inter_cols %in% colnames(inter_df))) {
        missing <- setdiff(inter_cols, colnames(inter_df))
        stop("inter_df missing columns: ", paste(missing, collapse = ", "))
    }

    if(!("prey_id" %in% colnames(prey_df)))
        stop("prey_df must have column: prey_id")

    bait_cols <- c("ip_id", "bait_id", "test_ctrl")
    if(!all(bait_cols %in% colnames(bait_df))) {
        missing <- setdiff(bait_cols, colnames(bait_df))
        stop("bait_df missing columns: ", paste(missing, collapse = ", "))
    }

    if(f < 0 || f >= 1) stop("f must be between 0 (inclusive) and 1 (exclusive)")
    # Map f to internal format (1-f) to match SAINTexpress CLI logic
    f_internal <- 1.0 - as.double(f)
    if(R < 1) stop("R must be a positive integer")
    if(L < 1) stop("L must be a positive integer")

    # Ensure correct types
    inter_df$ip_id <- as.character(inter_df$ip_id)
    inter_df$bait_id <- as.character(inter_df$bait_id)
    inter_df$prey_id <- as.character(inter_df$prey_id)
    inter_df$quant <- as.numeric(inter_df$quant)

    prey_df$prey_id <- as.character(prey_df$prey_id)
    if("gene_name" %in% colnames(prey_df))
        prey_df$gene_name <- as.character(prey_df$gene_name)

    bait_df$ip_id <- as.character(bait_df$ip_id)
    bait_df$bait_id <- as.character(bait_df$bait_id)
    bait_df$test_ctrl <- as.character(bait_df$test_ctrl)

    .SAINTexpress_spc_impl(inter_df, prey_df, bait_df, GO_df,
                            f_internal, as.integer(R), as.integer(L))
}


#' SAINTexpress-int
#' SAINTexpress scoring using intensity data (int mode).
#' Implements the SAINT-MRF model with Normal distribution
#' for intensity-based AP-MS data. Intensities are log-transformed
#' and z-score normalized internally.
#'
#' @title SAINTexpress-int
#' @param inter_df A data frame with columns: ip_id, bait_id, prey_id, quant.
#'   Each row represents one interaction observation (raw intensity, NOT logged).
#' @param prey_df A data frame with columns: prey_id.
#'   Optional columns: prey_length, gene_name.
#' @param bait_df A data frame with columns: ip_id, bait_id, test_ctrl.
#'   test_ctrl must be "T" (test) or "C" (control).
#' @param GO_df Optional data frame with columns: go_term, genes.
#'   genes should be space-separated gene names.
#' @param f Frequency parameter (0 <= f < 1). Default 0.
#' @param R Maximum number of replicates for scoring. Default 100.
#' @param L Maximum number of controls per prey. Default 100.
#'
#' @return A data frame with columns: Bait, Prey, PreyGene, Intensity,
#'   IntensitySum, AvgIntensity, NumReplicates, ctrlIntensity, AvgP, MaxP,
#'   TopoAvgP, TopoMaxP, SaintScore, OddsScore, FoldChange, BFDR, boosted_by.
#'
#' @author Hyungwon Choi, Damian Fermin (original SAINTexpress),
#'   Qingzhou Zhang (Rcpp integration)
#' @references Teo G, et al. SAINTexpress: improvements and additional
#'   features in Significance Analysis of INTeractome software.
#'   J Proteomics. 2014;100:37-43.
#'
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib SMAD
#' @export
#' @examples
#' bait_path <- system.file("exdata", "TIP49", "bait.dat", package = "SMAD")
#' prey_path <- system.file("exdata", "TIP49", "prey.dat", package = "SMAD")
#' inter_path <- system.file("exdata", "TIP49", "inter.dat", package = "SMAD")
#' 
#' bait <- read.table(bait_path, sep = "\t", header = FALSE, 
#'                    col.names = c("ip_id", "bait_id", "test_ctrl"))
#' prey <- read.table(prey_path, sep = "\t", header = FALSE, 
#'                    col.names = c("prey_id", "prey_length"))
#' inter <- read.table(inter_path, sep = "\t", header = FALSE, 
#'                     col.names = c("ip_id", "bait_id", "prey_id", "quant"))
#' 
#' result <- SAINTexpress_int(inter, prey, bait)
#' head(result)

SAINTexpress_int <- function(inter_df, prey_df, bait_df,
                              GO_df = NULL,
                              f = 0, R = 100L, L = 100L) {
    # Input validation
    if(!is.data.frame(inter_df))
        stop("inter_df must be a data.frame")
    if(!is.data.frame(prey_df))
        stop("prey_df must be a data.frame")
    if(!is.data.frame(bait_df))
        stop("bait_df must be a data.frame")

    inter_cols <- c("ip_id", "bait_id", "prey_id", "quant")
    if(!all(inter_cols %in% colnames(inter_df))) {
        missing <- setdiff(inter_cols, colnames(inter_df))
        stop("inter_df missing columns: ", paste(missing, collapse = ", "))
    }

    if(!("prey_id" %in% colnames(prey_df)))
        stop("prey_df must have column: prey_id")

    bait_cols <- c("ip_id", "bait_id", "test_ctrl")
    if(!all(bait_cols %in% colnames(bait_df))) {
        missing <- setdiff(bait_cols, colnames(bait_df))
        stop("bait_df missing columns: ", paste(missing, collapse = ", "))
    }

    if(f < 0 || f >= 1) stop("f must be between 0 (inclusive) and 1 (exclusive)")
    # Map f to internal format (1-f) to match SAINTexpress CLI logic
    f_internal <- 1.0 - as.double(f)
    if(R < 1) stop("R must be a positive integer")
    if(L < 1) stop("L must be a positive integer")

    # Ensure correct types
    inter_df$ip_id <- as.character(inter_df$ip_id)
    inter_df$bait_id <- as.character(inter_df$bait_id)
    inter_df$prey_id <- as.character(inter_df$prey_id)
    inter_df$quant <- as.numeric(inter_df$quant)

    prey_df$prey_id <- as.character(prey_df$prey_id)
    if("gene_name" %in% colnames(prey_df))
        prey_df$gene_name <- as.character(prey_df$gene_name)

    bait_df$ip_id <- as.character(bait_df$ip_id)
    bait_df$bait_id <- as.character(bait_df$bait_id)
    bait_df$test_ctrl <- as.character(bait_df$test_ctrl)

    .SAINTexpress_int_impl(inter_df, prey_df, bait_df, GO_df,
                            f_internal, as.integer(R), as.integer(L))
}
