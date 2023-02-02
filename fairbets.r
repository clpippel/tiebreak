# Fairbets.r
# Run gf.r first to process input
# Principle eigenvector and fairbets ranking.
# Volij, Slutzki, Marco Slikker · Peter Borm · René van den Brink
# Scoring of web pages and tournaments—axiomatizations, Slutzki Volij
# http://volij.co.il/publications/papers/tourna3.pdf
#
# Paired comparisons analysis:  an axiomatic approach to ranking methods
# http://eio.usc.es/pub/julio/papers/15._Paired_comparisons_analysis-_an_axiomatic_approach_to_rankings_in_tournaments_web.pdf
#
# Limiting behaviour of Markov Chains - Cross Validated.
# MIT Open Courseware, Discrete Stochastic Processes, chapetr 3
# https://ocw.mit.edu/courses/6-262-discrete-stochastic-processes-spring-2011/pages/course-notes/
#
# Calculation by poweriteration of the game matrix
# s1 = A.1, s2 = A.s1, sk = A.sk-1, ...
# It folows from the Perron-Frobenius theorem that in case of indivisible profiles,
# the normalized sequence (sk) converges to the principal eigenvector.
# If the fb matrix is reducible and a non-periodic then
# the result of the calculation is the probability of being in an absorbing state.
# If the fb matrix is periodic, the iteration does not converge.
#--------------------------------------------------------------------

cat('\n', "----------------------------- fairbets.r", "pefbev", exists("pefbev"), '\n')
cat("\n",ifelse(exists("pefbev"),"Perron Frobenius eigenvector","Fairbets") ,"\n")

# pefbev <- 1L                                      # test Perron eigenvector

# select most common SCC / league
largest_SCC <- which.max(table(SCC$membership))
# column vector of relative ratings
fb <- matrix(ifelse(SCC$membership %in% largest_SCC, 0, 0)) # set non-SCC players to NA (or zero when one group) <<<<<<
fb[apply(is.na(opponents), 1, all)] <- NA           # remove isolates, no opponent <<<<<<

mask_SCC  <- fb[c(opponents)]                       # set excluded opponents to NA
mask_SCC  <- mask_SCC + fb[,1]                      # set opponents of excluded players to NA
dim(mask_SCC) <- dim(opponents)                     # restore matrix (n x r)
fb <- fb + 1

# Fairbets ranking is invariant under adding win and loss points equally to a player
pts_ag <- rowSums(t.sparse(gfile + mask_SCC, opponents), na.rm=TRUE) + 1 # points against, sum of col(gfile), incl game against self

if (exists("pefbev") && pefbev == 1) pts_ag[] <- 1  # Eigen vector 
                       
fbfile <- (gfile + mask_SCC) / pts_ag               # A / losses, losses = sum col(A), indexed by player(n x r) 

iterations <- fb                                    # store iterations, one column for each iteration step
convergence <- c(NA)                                # concatenate missing value ( n x it. steps)

# iterate until fb is stable ---------------------
npls <- nrow(opponents)

maxit <- max(npls * npls, 100L)
steptol <- (1e-4 / 3) / npls                        # change in diff after <<<four>>> decimals
it <- 0L

# fixpoint, eigenvalue of: fb = ((Gf+I)/Losses).fb, Gf = Game file
stime <- system.time( 
while ((it <- it+1L) < maxit) {

  fbm = fb[c(opponents)]                            # rating of opponents  (n x r) for debugging
  dim(fbm) = dim(opponents)                         # restore matrix

  # ((Gf + I)/L).fb, sparse matrix/vector multiplication
  # add I for speed up, avoid division by zero  
  fbn <- rowSums(fbfile * fb[as.vector(opponents)], na.rm = TRUE) + (1/pts_ag)*fb # matrix multiplication + point by self
  fbn <- fbn / sum(fbn, na.rm=TRUE)                 # normalize by sum
  
  if (min(fbn == 0, na.rm = TRUE) == 1) break       # All zeros, this is going to repeat itself
  
  iterations  <- cbind(iterations, fbn)             # column bind, append next iteration
  diff <- fbn - fb                                  # change in iteration
  sd_diff <- sd(diff, na.rm=TRUE)

  convergence <- cbind(convergence, sd_diff)		# append to convergence

  fb <- fbn                                         # next approximation
 
  if (sd_diff < steptol ) break        # change in rating below limit
  if (!(it %% 100)) message(paste0("Fb progress: ", it, ", x-tol = ", sd_diff)) 
}
)

# validate FB solution
if (!exists("pefbev") ) {
  Wfb <- rowSums((gfile + mask_SCC) * fb[as.vector(opponents)], na.rm = TRUE) # gewonnen weddenschappen
  Lfb <- rowSums(t.sparse(gfile + mask_SCC, opponents), na.rm = TRUE) * fb    # verloren weddenschappen
  res <- Wfb - Lfb
  cat(sprintf("Validate FB, L2 Residual   = %5.2g\n" , sqrt(sum(res*res, na.rm=TRUE)) )) # Euclidean norm L2
}

cat("\n")
cat(sprintf("It / Maxit        : %d, %g\n", it, maxit) )
cat(sprintf("SD x-Diff / x_tol : %5.2g, %5.2g\n", sd(diff, na.rm=TRUE), steptol) )

if (!exists("report_tpr")) {
cat("", sep="\n"); print("Fair bets / Pev"); print(zapsmall(fb, digits = 4))
cat("", sep="\n"); print("Fb / Pev iterations"); print(iterations, digits = 4)
}
cat("Convergence FB", sep="\n"); print(c(convergence), digits = 5)
if (it <  maxit) cat("Convergence reached\n")
if (it >= maxit) cat(sprintf(">>> No convergence, steptol = %5.2g, difference = %5.2g, it = %d, maxit = %d\n" , steptol, sd_diff, it, maxit )) 
print(stime)

# plot fair bet main SCC
try( {
  dev.new()
  lytfb     <- layout.sugiyama(g, hgap = 1L, vgap = 1L)$layout
  lytfb[,2] <- fb
  lytfb[is.na(fb),2] <- max(fb[!is.na(fb)])           # replace NA by max
  plot(g, layout=lytfb, main = paste0(ifelse(exists("pefbev"),"Perron Frobenius eigenvector", "Fairbets")
                                       , ifelse(any(is.na(fb)), " (non-SCC, isolates removed)", " (all games, no SCCs)")
))
})
