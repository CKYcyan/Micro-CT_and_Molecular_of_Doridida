args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
root <- if (length(script_arg)) {
  normalizePath(file.path(dirname(sub("^--file=", "", script_arg[1])), ".."), winslash = "/", mustWork = TRUE)
} else normalizePath(getwd(), winslash = "/", mustWork = TRUE)
source(file.path(root, "R", "workflow_helpers.R"))

require_packages("ape")
args <- commandArgs(trailingOnly = TRUE)
repo <- if (length(args) >= 1) normalizePath(args[1], winslash = "/", mustWork = TRUE) else dirname(root)
output <- if (length(args) >= 2) args[2] else file.path(root, "results", "r_phylogeny_validation_20261008")
iq <- "Molecular_Phylogeny/2024_09_18-12_36_39 (IQtree)"
bi <- "Molecular_Phylogeny/MrBayes6133_1000to4000"
alignment <- file.path(repo, iq, "merge_align_L-ins-i__6133__602_625_1972_2934__218.phy")
lines <- readLines(alignment, warn = FALSE)
dims <- scan(text = lines[1], quiet = TRUE)
sequence_lines <- lines[-1][nzchar(trimws(lines[-1]))]
parts <- strsplit(trimws(sequence_lines), "\\s+")
taxa <- vapply(parts, "[", character(1), 1L)
sequences <- vapply(parts, function(x) paste(x[-1], collapse = ""), character(1))
stopifnot(identical(as.integer(dims), c(218L, 6133L)), length(taxa) == dims[1],
          !anyDuplicated(taxa), all(nchar(sequences) == dims[2]),
          all(grepl("^[ACGTRYSWKMBDHVN?acgtryswkmbdhvn-]+$", sequences)))
files <- c(file.path(iq, c("IQ_partition.nex.treefile", "IQ_partition.nex.contree")),
           file.path(bi, c("infile.nex.con.tre", "infile.nex.run1.t consensus.newick",
                          "infile.nex.run1.t consensus_1000to4000.nex",
                          "infile.nex.run2.t consensus_1000to4000.nex", "mcc_tree.nex", "mcc_tree.nwk")),
           file.path("Trace_History", c("pruned_tree_154.tre.nex", "pruned_tree_154_IQtree.tre.nex",
                      "pruned_tree_154_MrBayes.tre.nex", "pruned_tree_row_154.tre",
                      "pruned_tree_154_IQtree(ForTaxaRow).tre")))
rows <- lapply(files, function(relative) {
  path <- file.path(repo, relative)
  first <- readLines(path, n = 1, warn = FALSE)
  nexus <- grepl("#NEXUS", first, ignore.case = TRUE)
  warnings <- character()
  trees <- tryCatch(withCallingHandlers(
    if (nexus) ape::read.nexus(path) else ape::read.tree(path),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }),
    error = identity)
  is_trace <- startsWith(relative, "Trace_History")
  if (inherits(trees, "error")) {
    if (!is_trace) stop(relative, ": ", conditionMessage(trees))
    return(data.frame(File = relative, Tree = NA_integer_, Tips = NA_integer_,
      UniqueTips = NA, AlignmentTaxaMatch = NA, BranchLengthsPresent = NA,
      FiniteNonnegativeLengths = NA, MaxNumericNodeLabel = NA_real_, Status = "REVIEW",
      Note = paste("ape could not parse this Mesquite asset:", conditionMessage(trees))))
  }
  if (inherits(trees, "phylo")) trees <- list(trees)
  do.call(rbind, lapply(seq_along(trees), function(i) {
    tree <- trees[[i]]
    expected <- if (is_trace) 154L else 218L
    unique_tips <- !anyDuplicated(tree$tip.label)
    set_ok <- if (is_trace) all(tree$tip.label %in% taxa) else setequal(tree$tip.label, taxa)
    branches_ok <- is.null(tree$edge.length) || all(is.finite(tree$edge.length) & tree$edge.length >= 0)
    structural <- length(tree$tip.label) == expected && unique_tips && branches_ok
    if (!is_trace && (!structural || !set_ok)) stop("Primary tree structure mismatch: ", relative)
    node_numbers <- suppressWarnings(as.numeric(tree$node.label))
    valid <- node_numbers[is.finite(node_numbers)]
    maximum <- if (length(valid)) max(valid) else NA_real_
    mcc <- grepl("mcc_tree", relative)
    note <- if (mcc) "Node-label scale/provenance unresolved; do not interpret as posterior probabilities."
            else if (is_trace && !set_ok) "Renamed/numeric Mesquite tips need a verified mapping to alignment identifiers."
            else "Structure checked; support interpretation not validated."
    if (length(warnings)) note <- paste(note, "ape warning:", paste(unique(warnings), collapse = "; "))
    data.frame(File = relative, Tree = i, Tips = length(tree$tip.label), UniqueTips = unique_tips,
       AlignmentTaxaMatch = set_ok, BranchLengthsPresent = !is.null(tree$edge.length),
       FiniteNonnegativeLengths = branches_ok, MaxNumericNodeLabel = maximum,
       Status = if (structural && set_ok && !mcc && !length(warnings)) "STRUCTURE_PASS" else "REVIEW",
       Note = note)
  }))
})
dir.create(output, recursive = TRUE, showWarnings = FALSE)
checks <- do.call(rbind, rows)
write.csv(checks, file.path(output, "tree_structure_checks.csv"), row.names = FALSE)
write.csv(data.frame(Alignment = basename(alignment), Taxa = dims[1], Sites = dims[2],
                     MD5 = unname(tools::md5sum(alignment))),
          file.path(output, "alignment_check.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(output, "sessionInfo.txt"))
cat("PASS: 218 x 6133 alignment and eight primary tree structures/tip sets.\n")
print(table(checks$Status))
cat("REVIEW rows retain unresolved provenance/format limits; exit zero means the audit completed.\n")
cat("LIMIT: no inference rerun, convergence certification, MCC support repair or ancestral-state validation.\n")
