# Compute Arbiter Pro Statistics
# Unlicense (Ø) 2024 CLP, IJmuiden
# https://www.fmjd.org/downloads/FMJD_Annexes_2024_4_2-sig.pdf
# a. Solkoff or Buchholz:     the sum of scores of all opponents.
# b. Solkoff median:          the sum of scores of all opponents minus the highest score minus the lowest score.
# c. Short Truncated Solkoff: the sum of scores of all opponents minus the lowest score.
# d. Full Truncated Solkoff:  the sum of scores of all opponents minus the lowest score;
#                             if this is equal the sum of all opponents minus the 2 lowest scores …. Etc.
# e. Solkoff plus:            the sum of the opponents Solkoff.

if ("Total" %in% colnames(rdtable) )
  Total <- as.integer(rdtable[, "Total"]) else
  Total <- points
dim(Total) <- c(npls, 1)

Lows <- apply(matrix(Total[as.vector(opponents)], nrow(Total)), 1, FUN = min, na.rm = TRUE) # Lowest opponents value
Solk  <- rowSums(matrix(Total[as.vector(opponents)], nrow(Total)), na.rm=TRUE)              # Buchholz, Weerstand, Solkoff
SSolk <- Solk - Lows                                                                        # wpl, Short Truncated Solkoff
Solkp <-  rowSums(matrix(Solk[as.vector(opponents)], nrow(Total)), na.rm=TRUE)              # Solkoff plus

SB   <- rowSums(gfile * Total[as.vector(opponents)], na.rm=TRUE)  # Sonneborg-Berger, Neustadtl

head(cbind(rdtable[,gfrmcols], Total=Total, Solk, SB, SSolk, Solkp))
