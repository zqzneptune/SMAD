#' CompPASS
#' Comparative Proteomic Analysis Software Suite (CompPASS) is based on spoke
#' model. This algorithm was developed by Dr. Mathew Sowa for defining the
#' human deubiquitinating enzyme interaction landscape
#' (Sowa, Mathew E., et al., 2009). The implementation of this algorithm was
#' inspired by Dr. Sowa's online tutorial
#' (\url{http://besra.hms.harvard.edu/ipmsmsdbs/cgi-bin/tutorial.cgi}).
#' The output includes Z-score, S-score, D-score and WD-score.
#' This function also computes entropy and normalized 
#' WD-score. The source  code for this function was based on the source code. 
#'  \url{https://github.com/dnusinow/cRomppass}
#'
#'
#' @title CompPASS
#' @param datInput A dataframe with column names: idRun, idBait, idPrey, 
#' countPrey. Each row represent one unique protein captured in one 
#' pull-down experiment
#' @return A data frame consists of unique bait-prey pairs with 
#' Z-score, S-score,D-score and WD-score indicating interacting probabilities.

#' @author Qingzhou Zhang, \email{zqzneptune@hotmail.com}
#' 
#' @references Sowa, Mathew E., et al. "Defining the human deubiquitinating
#' enzyme interaction landscape." Cell 138.2 (2009): 389-403.
#' \url{https://doi.org/10.1016/j.cell.2009.04.042}
#' @references Huttlin, Edward L., et al. "The BioPlex network: a systematic 
#' exploration of the human interactome." Cell 162.2 (2015): 425-440. 
#' \url{https://doi.org/10.1016/j.cell.2015.06.043}
#' @references Huttlin, Edward L., et al. "Architecture of the human 
#' interactome defines protein communities and disease networks." 
#' Nature 545.7655 (2017): 505. 
#' \url{https://www.nature.com/articles/nature22366}

#' @importFrom dplyr group_by
#' @importFrom dplyr summarise
#' @importFrom dplyr mutate
#' @importFrom dplyr left_join
#' @importFrom dplyr filter
#' @importFrom dplyr n
#' @importFrom dplyr ungroup
#' @importFrom tidyr spread
#' @importFrom magrittr %>%
#' @importFrom stats quantile
#' @export
#' @examples
#' data(TestDatInput)
#' datScore <- CompPASS(TestDatInput)
#' head(datScore)
CompPASS <- function(datInput){
    .CompPASS_validate_input(datInput)
    
    . <- NULL
    idPrey <- NULL
    idBait <- NULL
    idRun <- NULL
    countPrey <- NULL
    MaxTSC <- NULL
    nBait <- NULL
    f_sum <- NULL
    AvePSM <- NULL
    MeanDiff <- NULL
    Mean <- NULL
    SD <- NULL
    WD_inner <- NULL
    WD_raw <- NULL
    WD_raw.factor <- NULL
    
    # Pre-calculations
    k <- length(unique(datInput$idBait))
    f_stats <- .CompPASS_f_stats(datInput)
    p_stats <- .CompPASS_p_stats(datInput)
    entropy_stats <- .CompPASS_entropy_stats(datInput)
    ave_psm_stats <- .CompPASS_ave_psm_stats(datInput)
    
    # Combine stats
    dat <- ave_psm_stats %>%
        left_join(., unique(entropy_stats[, c("BP", "Entropy")]), by = "BP") %>%
        left_join(., f_stats, by = "idPrey") %>%
        left_join(., p_stats[, c("BP", "p")], by = "BP") %>%
        mutate(`k` = k)
    
    # Prey stats (Mean, SD)
    prey_stats <- .CompPASS_prey_stats(datInput, k, f_stats)
    
    # Scoring
    output <- dat %>%
        left_join(., prey_stats[, c("idPrey", "Mean", "SD")], by = "idPrey") %>%
        mutate(`scoreZ` = (`AvePSM` - `Mean`) / (`SD`)) %>%
        mutate(`scoreS` = sqrt((`AvePSM`) * (`k`) / (`f_sum`))) %>%
        mutate(`scoreD` = sqrt((`AvePSM`) * (((`k`) / (`f_sum`))^`p`))) %>%
        mutate(`WD_inner` = (`k` / `f_sum`) * (`SD` / `Mean`)) %>%
        mutate(`WD_raw` = sqrt(`AvePSM` * (`WD_inner`^`p`)))
    
    # Weighted WD score
    i <- !(is.na(output$WD_raw) | is.nan(output$WD_raw))
    WD_raw_val <- output$WD_raw[i]
    WD_raw_factor <- unname(quantile(WD_raw_val, 0.98)[1])
    output[i, "scoreWD"] <- WD_raw_val / WD_raw_factor
    
    return(as.data.frame(output[, c("idBait", "idPrey", "AvePSM", "scoreZ", 
                                "scoreS", "scoreD", "Entropy", "scoreWD")]))
}

.CompPASS_validate_input <- function(datInput) {
    colInput <- c("idRun", "idBait", "idPrey", "countPrey")
    if(!is.data.frame(datInput)){
        stop("Input data should be data.frame")
    }
    if(!all(colInput %in% colnames(datInput))){
        missingCol <- setdiff(colInput, colnames(datInput))
        stop("Input data missing: ", paste(missingCol, collapse = ", "))
    }
}

.CompPASS_f_stats <- function(datInput) {
    unique(datInput[, c("idBait", "idPrey")]) %>%
        group_by(`idPrey`) %>%
        summarise(`f_sum` = n(), .groups = "drop")
}

.CompPASS_p_stats <- function(datInput) {
    datInput[, c("idBait", "idPrey")] %>%
        group_by(`idBait`, `idPrey`) %>%
        summarise(`p` = n(), .groups = "drop") %>%
        mutate(`BP` = paste(`idBait`, `idPrey`, sep = "~"))
}

.CompPASS_entropy_stats <- function(datInput) {
    datInput %>%
        group_by(`idBait`, `idPrey`, `idRun`) %>%
        summarise(`MaxTSC` = max(`countPrey`), .groups = "drop") %>% 
        group_by(`idBait`, `idPrey`) %>%
        mutate(`p` = (`MaxTSC` + 1/n())/(sum(`MaxTSC`) + 1)) %>% 
        mutate(`Entropy` = sum(-1*`p`*log(`p`, 2))) %>%
        mutate(`BP` = paste(`idBait`, `idPrey`, sep = "~")) %>%
        ungroup()
}

.CompPASS_ave_psm_stats <- function(datInput) {
    datInput %>%
        group_by(`idBait`, `idPrey`) %>%
        summarise(`AvePSM` = mean(`countPrey`), .groups = "drop") %>%
        mutate(`BP` = paste(`idBait`, `idPrey`, sep = "~"))
}

.CompPASS_prey_stats <- function(datInput, k, f_stats) {
    stats <- .CompPASS_ave_psm_stats(datInput)
    preyWithBaitOnly <- stats %>%
        group_by(`idPrey`) %>%
        summarise(`nBait` = n(), .groups = "drop") %>%
        filter(`nBait` == 1) %>%
        .$idPrey
    
    statsM <- stats %>%
        mutate(`AvePSM` = ifelse(`idBait` != `idPrey`, `AvePSM`,
                                    ifelse(`idPrey` %in% preyWithBaitOnly, 
                                        `AvePSM`, NA)))
    
    statsMatrix <- spread(`statsM`[, c("idBait", "idPrey", "AvePSM")], `idBait`, `AvePSM`)
    m <- as.matrix(statsMatrix[, -1])
    rownames(m) <- statsMatrix$idPrey 
    prey.mean <- rowSums(m, na.rm = TRUE) / k
    
    f_num_no_prey_map <- ifelse(names(prey.mean) %in% preyWithBaitOnly & 
                                    names(prey.mean) %in% unique(datInput$idBait),
                                k, k - f_stats$f_sum[match(names(prey.mean), f_stats$idPrey)])
    
    prey.stats <- data.frame(`idPrey` = names(prey.mean),
                                `Mean` = prey.mean,
                                `MeanDiff` = rowSums((m - prey.mean)^2, na.rm = TRUE),
                                `f_num_no_prey` = f_num_no_prey_map,
                                stringsAsFactors = FALSE) %>%
        mutate(`SD` = sqrt((`MeanDiff` + ((`Mean`^2) * (`f_num_no_prey`)))/(k - 1)))
    return(prey.stats)
}
