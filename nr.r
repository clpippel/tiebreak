# ---------------------------------------------------------------------
# nr.nr
# calculate relative ratings (x) as roots of We(x) - W
# root finding with Newton Raphson
# compute ratings of most frequent SCC
# re-calibrate ratings such that minimum rating = 1 (Sevilla)
# http://www.math.wichita.edu/~cma/stat774/ch2
# (2.6 Newton's method for a system of nonlineair equations)
# test tolerance
# L2 sqrt(sum(x.x)), euclidean distance
# BB,fnsane, L2 / sqrt(n) < tol
# Pracma, broyden, L2 < tol
# ---------------------------------------------------------------------

cat('\n', paste("----------------------------- nr.r", sep=", "), '\n')

celo <- 400/log(10)                                 # Elo constant in Logistic distribution
Pr  <- function(x) {return(1./(1. + exp(-x)))}      # logistic probability function
dPr <- function(x) {e <- Pr(x); return(e * (1-e))}  # derivative

#-------------------------------------------------------
# W_expected, expected score
# f(x) = ( f1(x), f2(x) ... fn(x) )
# fi(x) = Sum(We(xi - xj) - W), j = 1..rounds
# Par = max(∀ plusscore - minscore)
# ------------------------------------------------------
W_expected <- function(x) {
  dpopp <- x[,1] - x[as.vector(opponents)];         # difference between own rating, opponent rating
  dim(dpopp) = dim(opponents)                       # indexed by player, round         

  Jf     <<- -dPr(dpopp / celo) * Par               # Jacobian opponents: Jf(x) = δfi/δxk, k != i (global assign)
  Jfdiag <<- -rowSums(Jf, na.rm = TRUE)             # Gradient (diag):    Jf(x) = δfi/δxi, Sum(Jf) == 0

  Wem <- (Pr(dpopp / celo) - 0.5)                   # expected score above average(n x r)
  stopifnot(sum(Wem, na.rm=TRUE) < .Machine$double.eps * length(Wem)) # sum of rating differences
  We <- as.matrix(rowSums(Wem, na.rm=TRUE)) * Par   # expected score (n)
  
  return(rowSums(Wem, na.rm=TRUE) * Par)            # expected score (indexed by player)
}            

#----------------------------------------------------------------
# Conjugate gradients, solve A.x = b
# Sparse representation. (A, opponents) indexed by player x round
# A must be Positive Definite or
# semi PD and b in column space of A (Kaasschieter)
#-----------------------------------------------------------------
cg.solve  <- function(A, Adiag, b, reltol=NULL) {
  # Conjugate gradient method
  # A                                               # matrix
  # Adiag                                           # A[i,i]
  # opponents                                       # sparse matrix representation
  # reltol =										# relative error tolerance 

  if (is.null(reltol)) reltol <- sqrt(.Machine$double.eps)
  
  if (length(b) == 0L) return(b)
  x <- matrix(0, length(b), 1)                      # we assume x0 = 0
  r <- b                                            # Residual, r = b - A.x; x = 0
  p <- r                                            # search direction

  rsold <- sum(r * r)                               # r.r
  rs0   <- rsold                                    # initial error

  for (it in seq_along(b)) {
    if (rsold < rs0 * reltol * reltol) break        # found solution
    # as.vector(), otherwise it won't work for two dimensions
    Ap <- rowSums(A * p[as.vector(opponents)], na.rm = TRUE) # Ap <- A.p
    cgcum <<- cgcum + 1                             # count A.p iterations (global assignment)
    Ap <- Ap + Adiag * p                            # sparse file multiplication
    pAp <- sum(p * Ap)                              # A-orthogonal 
    if (pAp < .Machine$double.eps) break            # search directions exhausted

    alpha <- rsold / pAp
    x <- x + alpha * p
    r <- r - alpha * Ap                             # cumulative update (b - Ax)
    rsnew <- sum(r * r)                             # r.r next
    p = r + (rsnew / rsold) * p                     # next search direction
    rsold = rsnew
  }
  # browser()                                       # open browser environment to debug or inspect locals
  return(x)
}

#----------------------------------------------------
# choose largest indivisible profile domain
# set rrtg of excluded players to NA, otherwise 0
# recalculate W within main domain
# test sum wins and losses (otherwise nr will fail)
# ---------------------------------------------------

# column vector of relative ratings
rrtg <- mask_rtg <- as.matrix(ifelse(SCC$membership %in% largest_SCC, 0, NA)) # base is 0, exclude non largest SCC
diff <- rrtg

W <- as.matrix(rowSums(results + mask_rtg[as.vector(opponents)] + mask_rtg[,1] , na.rm=TRUE)) # "Copeland" score points (above/below average), within SCC largest.
stopifnot(sum(W) < .Machine$double.eps)             # zero sum

iterations <- rbind(rrtg, 0, 0, 0)                  # store iterations, one column for each iteration step
convergence <- c(NA)                                # concatenate missing value ( n x it. step)

cgcum <- 0                                          # count cg iterations cumulatively

# sqrt(kappa)/2 * log(2/cgeps) > Npls
cgeps <- exp(-2)                                    # cg residual error relative to initial residual, not too small.

steptol <- (1.0e-3 / 3)                             # stdev(diff) < steptol, 3 sigma change in diff after <<<e-three>>> decimals
maxit <- npls + 5L                                  # maximum number of nr iterations, 5 for small npls
maxit <- 100
We    <- W_expected(rrtg)                           # expected score
res0 <- We - W                                      # residual at x=0 

L2step <- NA
#----------------------------------------------------
# find roots of f(x) = We(x) - W 
# ---------------------------------------------------
stime <- system.time(   
for (it in seq_len(maxit)){

  We  <- W_expected(rrtg)                           # Expected score
  res <- We - W                                     # compute difference We, W, + Jacobian
  if (crossprod(res) < .Machine$double.eps) break   # solution found 
  kappa = max(Jfdiag) / min(Jfdiag[Jfdiag[]>0])     # condition number, λmax /  λmin 
  cg_th <- round(sqrt(kappa) * log(2 / cgeps) / 2 +0.5)  # rate of convergence CG method
  
# --------------------------------------------------------------------
# calculate next step by multivariate Newton Raphson update
# dx = Inverse([Jf(x)]).f(x), or solve unknown dx in Jf(x).dx = -f(x)
# --------------------------------------------------------------------
  nr_step <- cg.solve(Jf, Jfdiag, res, cgeps)
  rrtg.n <- rrtg - (nr_step * celo)
  stopifnot( abs(mean(rrtg.n, na.rm = TRUE)) < sqrt(.Machine$double.eps)) # average rating is invariant  
#---------------------------------------------------------------------
                                                   
  iterations  <- cbind(iterations, rbind(rrtg.n, kappa, cg_th, cgcum)) # column bind, append next iteration
  diff <- rrtg.n - rrtg                             # change in iteration
  L2step <- sqrt(sum(diff*diff, na.rm=TRUE) )       # L2 norm, sqrt(dot product diff.diff)
  convergence <- cbind(convergence, L2step)         # append to convergence
                                                    
  rrtg <- rrtg.n                                    # next rrtg approximation
                                                    
  if (sd(diff, na.rm=TRUE) < steptol ) break        # change in rating below limit
}

)

# rrtg <- rrtg - mean(rrtg, na.rm = TRUE)           # reset average rating to zero    
# if (min(rrtg) < 1) rrtg <- rrtg - min(rrtg) + 1   # assure ratings are strictly positive (compatible with Sevilla)
# rrtg <- round(rrtg)                               # round to integer
#---------------------------------------------------------------------

cat("\n")
if (!exists("report_tpr")) {
  cat("rrtg", sep="\n"); print(zapsmall(rrtg))
  cat("Iterations NR", sep="\n"); print(zapsmall(iterations))
  if (it >= maxit) print("Maximum number of iteration reached.")
}
cat("Convergence NR", sep="\n"); print(c(convergence), digits = 4)

# import extra player data in Sevilla (match=full name, data=Rating)
# write.csv2(cbind(rdtable[,2], W, rrtg), quote=FALSE)

# excluded players by SCC level
SCC$membership[SCC$membership != largest_SCC]
{
cat("\n")
cat(sprintf("# spelers     = %5d      , # rondes = %d\n", npls,  nrds)      )
cat(sprintf("CG iterations = %5d (cum), reltol   = %g\n", cgcum,  cgeps) )
cat(sprintf("NR iterations = %5d      , maxit    = %#.1f\n", it, maxit) )
cat(sprintf("SD x-Diff     = %7.2g    , x_tol    = %5.2g\n", sd(diff, na.rm=TRUE), steptol) )
cat(sprintf("L2 Residual   = %5.2g\n" , sqrt(sum(res*res)) )) 
cat(sprintf("SD Residual   = %5.2g\n" , sd(res, na.rm=TRUE) ))
cat(sprintf("XN Residual   = %5.2g\n" , max(abs(res), na.rm=TRUE) )) # max norm 
cat("\n")                                                                                           
}
print(stime)

uz   <- exp(rrtg / 400 * log(10))                   # Spielstärken, Zermelo, p. 451
uz[] <- 100* uz / sum(uz, na.rm=TRUE)               # genormeerd op som = 100
colnames(uz) <- "Zermelo, Berechnung der Turnier-ergebnisse, p.437"; uz[1,1]
PrElo <- 1 / (1 + exp(-rrtg/400*log(10)))           # probability Elo

try( {
lytnr     <- layout_with_fr(g)
lytnr[,2][which(!is.na(rrtg))] <- rrtg[which(!is.na(rrtg))]
lytnr <- norm_coords(lytnr, ymin = -1, ymax = 1)
} )
    
if (npls < 500L) {
 dev.new()
 plot(g, layout=lytnr, main = "Relative Elo Rating" )
}
