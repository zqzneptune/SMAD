## Suppress R CMD check NOTEs for variables used in data.table / dplyr NSE
## (variables that look like "no visible binding for global variable" to the
##  static analyser but are perfectly valid in the respective DSLs)

utils::globalVariables(c(
  ## data.table special symbols
  ".", ".N", ":=",
  ## Shared column names (dplyr / data.table NSE)
  "idBait", "idPrey", "idRun", "countPrey", "lenPrey",
  ## CompPASS helpers
  "AvePSM", "MaxTSC", "Mean", "MeanDiff", "f_num_no_prey", "nBait", "p",
  ## HG helpers
  "NormalSpec", "SumNS", "NSAF", "NormalNSAF",
  "UniprotID", "ppiTN", "tnA", "NMinTn", "tnB",
  ## Hart helpers
  "InteractorA", "InteractorB", "Freq", "TnA", "TnB", "totTn",
  ## PE helpers – .PE_prey_freq
  "obsPry", "f",
  ## PE helpers – .PE_calculate_spoke
  "cntBait", "cntPry", "cntBP", "sSum", "s", "BP", "finalS", "PB",
  ## PE helpers – .PE_calculate_matrix
  "pry1", "pry2", "f1", "f2", "finalM", "PPI"
))
