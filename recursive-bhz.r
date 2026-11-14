# Run gf.r first to process input
# compute Recursive Buchholz from FMJD matrix
# Single round robin tournament - rbhz = 0.5 x (W-L) / N, add draw against self, then all opponents are equal
# Paired comparisons: an axiomatic approach to ranking methods, p6, 2013
# ronde dossier: one game by round
# Gonzalez-Diaz: aggregated scores in A, Mmat = A + t(A)
#
# The rb-iteration oscilates when all cycles are even (bipartite graph).
# This might also occur when the domain is divisible.
# One solution is to introduce a tie against yourself in the crosstable.
# This does not effect the Laplacian and therefore also not the outcome of the iteration.
#
# 2023-1-31 - redefine maxit
#----------------------------------------------------------------------------------------

cat('\n', paste("----------------------------- recursive-bhz.r", sep=", "), '\n')

most_common_SCC <- which.max(table(SCC$membership))
rbhz <- as.matrix(ifelse(SCC$membership %in% most_common_SCC, 0, NA)) # base is 0, NA = SCCs, <<<<<<<<<<<<<<<<<<<<<<<
mask_SCC  <- rbhz[as.vector(opponents)]             # Set excluded opponents to NA.
mask_SCC  <- mask_SCC + as.vector(rbhz)             # Set opponents of excluded players to NA.
dim(mask_SCC) <- dim(opponents)                     # Restore matrix (n x r).

oppSCC <- opponents + mask_SCC                      # Opponents maximal SCC.
resSCC <- results + mask_SCC                        # Results maximal SCC.

Mmat   <- rowSums(!is.na(resSCC))                   # Number of games in max. SCC.
s     <- matrix(rowSums(resSCC, na.rm = TRUE))   
rbopp <- matrix(0, nrow(opponents), 1)              # average rbhz-ratings opponents, indexed by player.

iterations <- rbhz                                  # store iterations, one column for each iteration step.
convergence <- c(NA)                                # concatenate missing value ( n x it. step).

# iterate until rbhz is stable ---------------------
npls <- nrow(opponents)
stopifnot(npls > 1)                                 # it needs two to tango.

maxit <- max(npls * log(npls), 100L)                # max number of iterations (guess).
steptol <- 1e-4 / 3                                 # change in diff after <<<e-four>>> decimals
                                                    # change in diff after <<<three>>> decimals
it <- 0;

system.time(
while ((it <- it+1) <= maxit) {                     # a loop might occur when bipartite

  rboppm <- matrix(rbhz[as.vector(oppSCC)], nrow = npls) # rbhz of opponents  (n x r.
  rbopp = rowSums(rboppm, na.rm=TRUE) / Mmat        # average rbhz of opponents, indexed by player, ignore missing values
  
  rb.new <- rbopp  + s / Mmat                       # Next step in Banach fixpoint iteration. Gonzalex Diaz 2014, p.6 rbhz       == lsq

  try( iterations  <- cbind(iterations, rb.new), silent=TRUE) # column bind, append next iteration
  diff <- rb.new - rbhz                             # change in iteration
  xn_diff <- max(abs(diff), na.rm=TRUE)             # max norm
  convergence <- cbind(convergence, xn_diff )       # append to convergence

  rbhz <- rb.new                                    # next rbhz approximation
  if ( all(rbhz == 0L, na.rm = TRUE) ) break        # All zeros, this is going to repeat itself 

  if (xn_diff < steptol) break                      # change in rating below limit
  if (it %% 100 == 0) cat(sprintf("Progress = %d, sd(diff) = %7.4g\n", it, xn_diff))
}
)

{
cat(sprintf("It / Maxit        : %d, %g\n", it, maxit))
cat(sprintf("XN x-Diff / x_tol : %5.2g, %5.2g\n", xn_diff, steptol) )

if (!exists("report_tpr")) {
    cat("rbhz", sep="\n"); print(rbhz, digits = 4)
    cat("iterations ", sep="\n"); print(iterations, digits = 4)
}
cat("Convergence rbhz", sep="\n"); cat(convergence, '\n')
cat("#Iterations", sep="\n"); print(it)
if (it >= maxit) print("Maximum number of iteration reached.")
}

