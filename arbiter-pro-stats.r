# Compute Arbiter Pro Statistics
# Unlicense (Ø) 2024 CLP, IJmuiden
# https://www.fmjd.org/downloads/FMJD_Annexes_2024_4_2-sig.pdf
# https://www.fmjd.org/downloads/DA/Dra_sort_criteria.pdf
# a. Solkoff or Buchholz:     the sum of scores of all opponents. Opponents Bye score is 0.
# b. Solkoff median:          the sum of scores of all opponents minus the highest score minus the lowest score.
# c. Short Truncated Solkoff: the sum of scores of all opponents minus the lowest score.
# d. Full Truncated Solkoff:  the sum of scores of all opponents minus the lowest score;
#                             if this is equal the sum of all opponents minus the 2 lowest scores …. Etc.
# e. Solkoff plus:            the sum of the opponents Solkoffs.
# Run gf.r first.

stopifnot("Fed." %in% colnames(rdtable))  

if ("Total" %in% colnames(rdtable))
  Total <- as.integer(rdtable[, "Total"]) else
  if ("Pt" %in% colnames(rdtable))
    Total <- as.integer(rdtable[, "Pt"]) else
    if ("Pts." %in% colnames(rdtable))
      Total <- as.integer(rdtable[, "Pts."]) else
      Total <- points
dim(Total) <- c(npls, 1)                                              # Score points indexed by player.

Lows <- apply(matrix(Total[as.vector(opponents)], nrow(Total)), 1, FUN = min, na.rm = TRUE) # Lowest opponents value.
Solk  <- rowSums(matrix(Total[as.vector(opponents)], nrow(Total)), na.rm = TRUE)            # Buchholz, Weerstand, Solkoff.
SSolk <- Solk - Lows                                                                        # wpl, Short Truncated Solkoff.
Solkp <-  rowSums(matrix(Solk[as.vector(opponents)], nrow(Total)), na.rm = TRUE)            # Solkoff plus.
SB   <- rowSums(gfile * Total[as.vector(opponents)], na.rm = TRUE)    # Sonneborg-Berger, Neustadtl.

mmm    <- matrix(Total[as.vector(opponents)], nrow(Total))            # Opponents scores.
mmm[is.na(mmm)] <- 0                                                  # Opponents Bye score is 0.
FTSolk <-
  t(
    apply(mmm, 1,
      function(x) head(sum(x) - cumsum(sort(x, decreasing = FALSE)), -1) # Full Truncated Solkoff.
    )
  )                                                                # Truncated Solkoffs, indexed by player.
colnames(FTSolk) <- paste0("Ts", seq_len(ncol(FTSolk)))               # Ts1, Ts2, ...

# Sort by Federation, Total, Full Truncated Solkoff (decr)
# Previous FMJD standard
dd <- data.frame(rdtable[, gfrmcols], Tot. = Total, Solk, SB, SSolk, Solkp, FTSolk)

id2 <- with(dd, order(Fed., -Tot., -Ts1, -Ts2, -Ts3, -Ts4, -Ts5, -Ts6, -Ts7, -Ts8)) # Naive sort.


if ("Fed." %in% colnames(dd)) {
  # Number of Full Truncated Solkoff columns: number of rounds - 1 (Swiss)
  bycols <- dd[, c("Fed.", "Tot.", tail(names(dd), nrds - 1))]       # Data frame containing By keys.
  bysort <- c(FALSE, TRUE, rep(TRUE, nrds - 1))                      # Decending, Ascending ...
  id3     <- do.call("order",
                     c(bycols,
                       decreasing = list(bysort),
                       method = "radix"))
}
all(id2 == id3)
d3 <- dd[id3, ]
subset(d3, Fed. == "LTU")

# ------------------------------------------------------------------
library(bench)
bench::mark(
  "order"   = with(dd, order(Fed., -Tot., -Ts1, -Ts2, -Ts3, -Ts4, -Ts5, -Ts6, -Ts7, -Ts8)),
  "do.call" = { # Number of Full Truncated Solkoff columns: number of rounds - 1 (Swiss)
                bycols <- dd[, c("Fed.", "Tot.", tail(names(dd), nrds - 1))]  # Data frame containing By keys.
                bysort <- c(FALSE, TRUE, rep(TRUE, nrds - 1))                 # Decending, Ascending ...
                do.call("order",
                        c(bycols,
                          decreasing = list(bysort),
                          method = "radix"))
  }
)
