#' Purification Enrichment (PE) Score Calculation
#'
#' Computes PE scores for protein-protein interactions from AP-MS data.
#' Based on Collins et al. (2007) Bayes classifier framework.
#'
#' @title PE
#' @param datInput  data.frame with columns: idRun, idBait, idPrey
#' @param rBait     probability a true interaction is detected (default 0.37)
#' @param cntPseudo pseudocount for prey frequency estimation (default 1)
#' @return data.frame with columns: PPI, InteractorA, InteractorB, BP, PB,
#'         spokeBP, spokePB, matrixPP, PE
#'
#' @author Qingzhou Zhang, \email{zqzneptune@hotmail.com}
#' @references Collins, Sean R., et al. "Toward a comprehensive atlas of the 
#' physical interactome of Saccharomyces cerevisiae." Molecular & Cellular 
#' Proteomics 6.3 (2007): 439-450.
#'
#' @importFrom data.table as.data.table
#' @importFrom data.table setDT
#' @importFrom data.table setnames
#' @importFrom data.table uniqueN
#' @importFrom data.table CJ
#' @importFrom RcppAlgos comboGeneral
#' @export
#' @examples
#' data(TestDatInput)
#' datScore <- PE(TestDatInput)
#' head(datScore)
PE <- function(datInput, rBait = 0.37, cntPseudo = 1) {
  .PE_validate_input(datInput)

  # Remove self-interactions
  datInput <- datInput[datInput$idBait != datInput$idPrey, ]

  # Convert to data.table for fast joins/aggregations
  dt <- data.table::as.data.table(datInput)

  # ── Prey frequency (Eq. 9) ────────────────────────────────────────────────
  fPrey <- .PE_prey_freq(dt, cntPseudo)

  # ── PPI space ─────────────────────────────────────────────────────────────
  ppi_space <- .PE_generate_space(dt)

  # ── Spoke scores (Eq. 5-7) ────────────────────────────────────────────────
  spoke <- .PE_calculate_spoke(dt, fPrey, rBait)

  # ── Matrix scores (Eq. 10-12) ─────────────────────────────────────────────
  mat <- .PE_calculate_matrix(fPrey, ppi_space$pkPPS, rBait)

  # ── Assemble output (Eq. 4) ───────────────────────────────────────────────
  out <- ppi_space$ppsOut
  out <- merge(out, spoke$spokeBP, by = "BP",  all.x = TRUE)
  out <- merge(out, spoke$spokePB, by = "PB",  all.x = TRUE)
  out <- merge(out, mat,           by = "PPI", all.x = TRUE)

  # PE = S_ij + S_ji + M_ij  (NA → 0 when component absent)
  cols <- c("spokeBP", "spokePB", "matrixPP")
  for (col in cols) {
    if (!col %in% names(out)) out[[col]] <- NA_real_
  }
  out$PE <- rowSums(out[, cols], na.rm = TRUE)
  out
}

# ── Input validation ──────────────────────────────────────────────────────────
.PE_validate_input <- function(datInput) {
  required <- c("idRun", "idBait", "idPrey")
  if (!is.data.frame(datInput))
    stop("Input data should be a data.frame")
  missing <- setdiff(required, colnames(datInput))
  if (length(missing))
    stop("Input data missing columns: ", paste(missing, collapse = ", "))
}

# ── Prey frequencies (Eq. 9) ─────────────────────────────────────────────────
.PE_prey_freq <- function(dt, cntPseudo) {
  totObs   <- nrow(dt)
  nDistPry <- data.table::uniqueN(dt$idPrey)

  # Count occurrences of each prey
  obs <- dt[, .N, by = idPrey]
  data.table::setnames(obs, "N", "obsPry")

  obs[, f := (obsPry + cntPseudo) / (totObs + nDistPry * cntPseudo)]
  obs
}

# ── PPI space (bait-prey + prey-prey pairs) ───────────────────────────────────
.PE_generate_space <- function(dt) {
  # Bait-prey pairs (sorted so A < B lexicographically)
  bp_mat  <- dt[, {
    s <- pmin(idBait, idPrey)
    b <- pmax(idBait, idPrey)
    list(A = s, B = b)
  }]
  pkBP <- unique(paste(bp_mat$A, bp_mat$B, sep = "~"))

  # Prey-prey pairs per run
  runs   <- split(dt$idPrey, dt$idRun)
  pkPPS  <- lapply(runs, function(preys) {
    preys <- unique(preys)          # avoid duplicates within a run
    if (length(preys) < 2) return(character(0))
    combos <- RcppAlgos::comboGeneral(sort(preys), 2)
    paste(combos[, 1], combos[, 2], sep = "~")
  })

  pkPP   <- unique(unlist(pkPPS, use.names = FALSE))

  # Union of both pair types
  allPPI <- unique(c(pkBP, pkPP))

  # Build output data.frame
  parts <- strsplit(allPPI, "~", fixed = TRUE)
  ppsOut <- data.frame(
    PPI         = allPPI,
    InteractorA = vapply(parts, `[[`, character(1), 1L),
    InteractorB = vapply(parts, `[[`, character(1), 2L),
    stringsAsFactors = FALSE
  )
  ppsOut$BP <- paste(ppsOut$InteractorA, ppsOut$InteractorB, sep = ":")
  ppsOut$PB <- paste(ppsOut$InteractorB, ppsOut$InteractorA, sep = ":")

  list(ppsOut = ppsOut, pkPPS = pkPPS)
}

# ── Spoke (bait-prey) scores (Eq. 5-7) ───────────────────────────────────────
#
# For each bait i, prey j, purification k:
#   s_ijk = log10( (r + (1-r)*p_ijk) / p_ijk )   if j seen in k   (Eq. 6)
#   s_ijk = log10(1 - r)                            if j not seen    (Eq. 7)
#
# p_ijk = 1 - exp( -f_j * n_prey_k * n_bait_i )                    (Eq. 8)
#
# S_ij  = sum_k s_ijk  over all runs where i is bait
#
.PE_calculate_spoke <- function(dt, fPrey, rBait) {
  data.table::setDT(dt)
  log1mr <- log10(1 - rBait)

  # Number of times each bait was used
  nBait <- dt[, .N, by = idBait]            # N = number of runs per bait
  data.table::setnames(nBait, "N", "cntBait")

  # Actually: cntBait = number of distinct runs per bait
  # (datInput rows already correspond to individual bait-prey observations)
  cntBait_dt <- unique(dt[, .(idRun, idBait)])[, .N, by = idBait]
  data.table::setnames(cntBait_dt, "N", "cntBait")

  # Prey count per run (n_prey_k)
  cntPry_dt <- dt[, .(cntPry = .N), by = idRun]

  # Merge everything onto observed bait-prey rows
  dS <- merge(dt,      cntBait_dt, by = "idBait", all.x = TRUE)
  dS <- merge(dS,      fPrey[, .(idPrey, f)], by = "idPrey", all.x = TRUE)
  dS <- merge(dS,      cntPry_dt,  by = "idRun",  all.x = TRUE)

  # Eq. 8 & 6: p and s for observed pairs
  dS[, p  := 1 - exp(-f * cntPry * cntBait)]
  dS[, s  := log10((rBait + (1 - rBait) * p) / p)]
  dS[, BP := paste(idBait, idPrey, sep = ":")]

  # Aggregate observed s over runs (cntBP = number of observations, sSum = sum)
  pkSU <- dS[, .(cntBP = .N, sSum = sum(s)), by = BP]

  # Full bait × prey space (only bait proteins × all unique preys)
  allBaits <- unique(dt$idBait)
  allPreys <- unique(dt$idPrey)
  spacePk  <- data.table::CJ(idBait = allBaits, idPrey = allPreys)[idBait != idPrey]
  spacePk[, BP := paste(idBait, idPrey, sep = ":")]

  spacePk <- merge(spacePk, pkSU,       by = "BP",     all.x = TRUE)
  spacePk <- merge(spacePk, cntBait_dt, by = "idBait", all.x = TRUE)
  spacePk[is.na(cntBP), cntBP := 0L]
  spacePk[is.na(sSum),  sSum  := 0]

  # nonBP runs = total bait runs − runs where prey was observed (Eq. 7 contribution)
  spacePk[, finalS := sSum + (cntBait - cntBP) * log1mr]

  # Return as BP and PB keyed tables
  spokeBP <- spacePk[, .(BP, spokeBP = finalS)]
  spokePB <- spacePk[, .(PB = BP, spokePB = finalS)]
  # Rename PB column: swap direction
  spacePk[, PB := paste(idPrey, idBait, sep = ":")]
  spokePB <- spacePk[, .(PB, spokePB = finalS)]

  list(spokeBP = spokeBP, spokePB = spokePB)
}

# ── Matrix (prey-prey) scores (Eq. 10-12) ────────────────────────────────────
#
# m_ijk = log10( (r + (1-r)*p_ijk) / p_ijk )                       (Eq. 11)
# p_ijk = 1 - exp( -f_i * f_j * n_tot_prey_prey )                  (Eq. 12)
# M_ij  = sum_k m_ijk  (sum over co-occurrence runs)
#
.PE_calculate_matrix <- function(fPrey, pkPPS, rBait) {
  # Flatten prey-prey pairs to a frequency table
  allPP  <- unlist(pkPPS, use.names = FALSE)
  if (length(allPP) == 0L)
    return(data.frame(PPI = character(0), matrixPP = numeric(0)))

  totPPS <- length(allPP)
  freq   <- data.table::as.data.table(table(allPP))
  data.table::setnames(freq, c("allPP", "N"), c("PPI", "Freq"))

  # Split PPI key into individual prey IDs
  parts <- strsplit(freq$PPI, "~", fixed = TRUE)
  freq[, pry1 := vapply(parts, `[[`, character(1), 1L)]
  freq[, pry2 := vapply(parts, `[[`, character(1), 2L)]

  # Attach f values for each prey
  f_dt <- fPrey[, .(idPrey, f)]
  freq <- merge(freq, f_dt, by.x = "pry1", by.y = "idPrey", all.x = TRUE)
  data.table::setnames(freq, "f", "f1")
  freq <- merge(freq, f_dt, by.x = "pry2", by.y = "idPrey", all.x = TRUE)
  data.table::setnames(freq, "f", "f2")

  # Eq. 12 & 11
  freq[, p      := 1 - exp(-f1 * f2 * totPPS)]
  freq[, finalM := log10((rBait + (1 - rBait) * p) / p) * Freq]

  freq[, .(PPI, matrixPP = finalM)]
}