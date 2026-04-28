#!/usr/bin/env Rscript
# =============================================================================
# predictPASC_reduced.R - Predict PASC using the parsimonious reduced RF
#
#
# Usage:
#    Rscript predictPASC_reduced.R <bsseqObject.rda> [output.csv] [--offset N]
#
# Uses the top 25 CpGs selected by MeanDecreaseGini importance.
# --offset N  Add N to CpG positions before matching (default: 0).
# =============================================================================
suppressPackageStartupMessages({
  library(bsseq); library(GenomicRanges); library(randomForest)
})

args <- commandArgs(trailingOnly = TRUE)
coord.offset <- 0L
offset.idx <- which(args == "--offset")
if (length(offset.idx) > 0) {
    coord.offset <- as.integer(args[offset.idx[1] + 1])
    args <- args[-c(offset.idx[1], offset.idx[1] + 1)]
}
if (length(args) < 1) {
    cat("Usage: Rscript predictPASC_reduced.R <bsseq.rda> [output.csv] [--offset N]\n")
    quit(status = 1)
}
input.rda <- args[1]
output.csv <- if (length(args) >= 2) args[2] else {
    baseName <- sub("\\.(rda|RData)$", "", basename(input.rda), ignore.case = TRUE)
    file.path(dirname(normalizePath(input.rda)), paste0(baseName, "_pasc_predictions.csv"))
}

modelsRDA <- "rf_reduced_models.rda"
load(modelsRDA)
cat("Loaded reduced RF model (top", parsimonious_k, "CpGs)\n")

model.use <- parsimonious_model$model
cpg.ids <- parsimonious_model$cpg_ids
thr <- parsimonious_model$threshold
n.features <- parsimonious_k

envNew <- new.env()
load(input.rda, envir = envNew)
bsNew <- NULL
for (obj in ls(envNew)) {
    if (inherits(get(obj, envir = envNew), "BSseq")) { bsNew <- get(obj, envir = envNew); break }
}
if (is.null(bsNew)) stop("No BSseq object found in: ", input.rda)
cat("Samples:", ncol(bsNew), "| CpGs:", nrow(bsNew), "\n")

sample.ids <- sampleNames(bsNew)
if (!length(sample.ids)) sample.ids <- paste0("Sample_", seq_len(ncol(bsNew)))

new_pos <- start(bsNew)
if (coord.offset != 0L) new_pos <- new_pos + coord.offset
new.ids <- paste(as.character(seqnames(bsNew)), new_pos, sep = ":")
match.idx <- match(cpg.ids, new.ids)
n.found <- sum(!is.na(match.idx))
cat("Features matched:", n.found, "/", n.features, "\n")
if (n.found < n.features * 0.8) warning("Ruh-roh! Less than 80% of features found.")

bs.sel <- bsNew[match.idx[!is.na(match.idx)], ]
cov.sel <- getCoverage(bs.sel, type = "Cov")
M.sel <- getCoverage(bs.sel, type = "M")
beta.sel <- M.sel / cov.sel
beta.sel[is.nan(beta.sel) | cov.sel < 5L] <- NA
beta.sel <- pmax(pmin(beta.sel, 0.999), 0.001)
m.sel <- log2(beta.sel / (1 - beta.sel))

X_new <- matrix(NA_real_, nrow = ncol(bsNew), ncol = n.features)
X_new[,!is.na(match.idx)] <- t(m.sel)

feat_means_red <- if (!is.null(parsimonious_model$feat_means)) {
    parsimonious_model$feat_means
} else {
    apply(X_new, 2, function(x) { m <- mean(x, na.rm=TRUE); if(is.nan(m)||is.na(m)) 0 else m })
}
for (j in seq_len(n.features)) {
    na_j <- which(is.na(X_new[, j]) | is.nan(X_new[, j]) | is.infinite(X_new[, j]))
    if (length(na_j) > 0) X_new[na_j, j] <- feat_means_red[j]
}

probs <- predict(model.use, X_new, type = "prob")[, "PASC"]
labels <- ifelse(probs >= thr, "PASC", "Control")
pct <- round(rank(probs) / length(probs) * 100, 1)

out <- data.frame(SampleID = sample.ids, Probability_PASC = round(probs, 6),
    Predicted_Label = labels, Score_Percentile = pct,
    stringsAsFactors = FALSE)

pd <- tryCatch(as.data.frame(pData(bsNew), stringsAsFactors = FALSE), error = function(e) NULL)
if (!is.null(pd) && ncol(pd) > 0) {
    out <- cbind(out, pd)
    cat("pData columns appended:", paste(colnames(pd), collapse = ", "), "\n")
}

write.csv(out, output.csv, row.names = FALSE)
cat("Predictions saved to:", output.csv, "\n")
cat("Predicted PASC:", sum(labels == "PASC"), "| Control:", sum(labels == "Control"), "\n")

