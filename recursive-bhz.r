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

Mmat <- rowSums(!is.na(gfile))                      # matches matrix, number of games
pctpts = points / (Par * Mmat) - 0.5                # average plus/min points score
message("pctpts: ", all(abs(pctpts - s / (wins + draws + losses) / Par) < .Machine$double.eps, na.rm=TRUE)) # test pctpt, % plus score

rbhz <-  matrix(0, nrow(opponents), 1)              # rbhz-rating player
rbopp <- matrix(0, nrow(opponents), 1)              # average rbhz-ratings opponents, indexed by player
iterations <- rbhz                                  # store iterations, one column for each iteration step
convergence <- c(NA)                                # concatenate missing value ( n x it. step)

# iterate until rbhz is stable ---------------------
npls <- nrow(opponents)
stopifnot(npls > 1)                                 # it needs two to tango

maxit <- max(npls * npls, 100L)                     # max number of iterations (guess)

steptol <- 1e-4 / 3                                 # change in diff after <<<e-four>>> decimals                                              # change in diff after <<<three>>> decimals
it <- 0;

system.time(
while ((it <- it+1) <= maxit) {                     # a loop might occur when bipartite
  rboppm = rbhz[as.vector(opponents)]               # rbhz of opponents  (n x r) for debugging
  dim(rboppm) = dim(opponents)                      # restore matrix
  rbopp = rowSums(rboppm, na.rm=TRUE) / Mmat        # average rbhz of opponents, indexed by player, ignore missing values
  
# rb.new <- rbopp  + pctpts                         # Next step in Banach fixpoint iteration. Gonzalex Diaz 2014, p.6 rbhz / Par == lsq
  rb.new <- rbopp  + s / Mmat                       # Next step in Banach fixpoint iteration. Gonzalex Diaz 2014, p.6 rbhz       == lsq
  rb.new <- rb.new - mean(rb.new, na.rm = TRUE)     # reset sum to zero, ignore missing values
  
  try( iterations  <- cbind(iterations, rb.new), silent=TRUE) # column bind, append next iteration
  diff <- rb.new - rbhz                             # change in iteration
  xn_diff <- max(abs(diff), na.rm=TRUE)             # max norm
  convergence <- cbind(convergence, xn_diff )       # append to convergence

  if ( all(rb.new == 0L, na.rm = TRUE) ) break      # All zeros, this is going to repeat itself 
  rbhz <- rb.new                                    # next rbhz approximation

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

