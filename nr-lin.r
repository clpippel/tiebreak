# ---------------------------------------------------------------------
# nr-lin.nr
# calculate relative ratings (x) as roots of We(x) - W = 0
# root finding with Newton Raphson
# compute ratings of most frequent SCC
# re-calibrate ratings such that minimum rating = 1
# ---------------------------------------------------------------------

cat('\n', paste("----------------------------- nr-lin.r", sep=", "), '\n')

slope <- 4                                          # inverse of slope
celo <-  200                                        # Elo class interval
Pr  <- function(DP) {return(DP/slope + 0.5)}        # lineair with slope 1/4
dPr <- function(DP) {return( (!is.na(DP)) / slope) }# derivative

# Conjugate gradients, solve linear equation  A.x = b -----------------
# A must be Positive Definite or
# semi PD and b in column space of A, e.g. Σb = 0
conjgrad  <- function(A, Adiag, b, Eps) {

# Conjugate gradient method
# matrix A in gamefile format (sparse)
# Solve A.x = b
# A     = Off diagonal elements matrix A (indexed through opponents)
# Adiag = diagonal elements matrix A
# first approximation x = 0

# set breakpoint for debugging
# browser()
  if (length(b) == 0) return(b)
  
# we assume x0 = 0  
  x <- matrix(0, length(b), 1)

# Residual, r = b - A * x; x = 0
  r <- b

# search direction
  p <- r

# r.r
  rsold <- drop(crossprod(r))
  rs0   <- rsold
  
  for (it in seq_along(b)) {
  # solution found within error limit
    if (rsold < rs0 * Eps * Eps) break
    
  # update global counter
    cgcum <<- cgcum + 1
  # sparse matrix multiplication A.p
  # force vector when two columns
    Ap <- rowSums(A * p[as.vector(opponents)], na.rm = T)
    Ap <- Ap + Adiag * p

  # A-orthogonal
    pAp <- drop(crossprod(p, Ap))
    
  # search directions exhausted
    if (pAp < .Machine$double.eps) break

    alpha <- rsold / pAp
    x <- x + alpha * p
    
  # cumulative update (b - Ax)
    r <- r - alpha * Ap
    
  # r.r next
    rsnew <- drop(crossprod(r))
    
  # next search direction
    p = r + (rsnew / rsold) * p
    rsold = rsnew 
  }
# set breakpoint for debugging, inspecting local
# browser()
  return(x)
}

#----------------------------------------------------
# choose largest indivisible profile domain
# set rrtg of excluded players to NA, otherwise 0
# calculate mask_SCC (n x r), NA for excluded players
# recalculate wins, losses, draws within main domain
# test sum wins and losses (otherwise nr will fail)
# ---------------------------------------------------

# column vector of relative ratings
most_common_SCC <- which.max(table(SCC$membership))
Par <- max(results, na.rm=TRUE) - min(results, na.rm=TRUE) # maximal score
rrtg <- matrix(ifelse(SCC$membership %in% most_common_SCC, 0, 0)) # base is 0, no SCCs, <<<<<<<<<<<<<<<<<<<<<<<
diff <- rrtg - rrtg                                 # set to zero
mask_SCC  <- rrtg[as.vector(opponents)]             # set excluded opponents to NA
mask_SCC  <- mask_SCC + rrtg[,1]                    # set opponents of excluded players to NA
dim(mask_SCC) <- dim(opponents)                     # restore matrix (n x r)

W <- matrix(rowSums(results + mask_SCC, na.rm=TRUE)) # actual score points (above/below average)
stopifnot(sum(W) < .Machine$double)                 # zero sum

iterations <- rbind(rrtg, 0, 0, 0)                  # store iterations, one column for each iteration step
convergence <- c(NA)                                # concatenate missing value ( n x it. step)

cgcum <- 0                                          # count cg iterations cumulatively
cgeps <- 1E-5                                       # cg residual error relative to initial residual  <<<<<<<<<<<<<
steptol <- (1.0e-3 / 3)^2 * npls                    # sigma(diff) / 3 < 0.00001, change in diff after <<<e-five>>> decimals
maxit <- 2 * log(npls) + 5                          # maximum number of nr iterations, 5 for small npls

NRstart <- Sys.time()
#------------------------------------------
# find roots of f(x) = We(x) - W 
# f(x) = ( f1(x), f2(x) ... fn(x) )
# fi(x) = Sum(We(xi - xj) - W), j = 1..rounds
# -----------------------------------------
for (it in seq_len(maxit)){

  rtopp     = rrtg[as.vector(opponents)]            # rating opponents
  dim(rtopp) = dim(opponents)                       # restore matrix (n x r)
  dpopp <- (rrtg[,1] - rtopp)                       # rating difference (dp) with opponents
  dpopp <- dpopp + mask_SCC                         # remove excluded opponents

  Wem <- (Pr(dpopp / celo) - 0.5)                    # expected score above average(n x r)
  stopifnot(sum(Wem, na.rm=TRUE) < .Machine$double.eps * length(Wem)) # sum of rating differences
  We <- matrix(rowSums(Wem, na.rm=TRUE)) * Par      # expected score (n)

  res = We - W                                      # residual = f(x)
  dd <- drop(crossprod(res))                        # initialize dd
  
  if (drop(crossprod(res)) < .Machine$double.eps) break # solution found
  Jf <- -dPr(dpopp / celo) * Par                    # Jacobian opponents:   Jf(x) = δfi/δxk, k != i
  
  Jfdiag <- -rowSums(Jf, na.rm = TRUE)              # Jacobian self (diag): Jf(x) = δfi/δxi, Sum(Jf) == 0
  kappa = max(Jfdiag) / min(Jfdiag[Jfdiag[]>0])     #  λmax /  λmin 
  cg_th <- round(sqrt(kappa) * log(2 / cgeps) / 2 +0.5)  # rate of convergence CG method
  
# --------------------------------------------------------------------
# calculate next step by multivariate Newton Raphson update
# dx = Inverse([Jf(x)]).f(x), or solve unknown dx in Jf(x).dx = -f(x)
# --------------------------------------------------------------------
  nr_step <- conjgrad(Jf, Jfdiag, res, cgeps)
  t <- rrtg - (nr_step * celo)
#---------------------------------------------------------------------
    
  tx = mean(t, na.rm = TRUE)                        # average rating, ignore missing values
  t <- t - tx                                       # reset average rating to zero
                                                    
  iterations  <- cbind(iterations, rbind(t, kappa, cg_th, cgcum)) # column bind, append next iteration
  diff <- t - rrtg                                  # change in iteration
  dd <- sum(diff*diff, na.rm=TRUE)                  # dot product diff.diff
  convergence <- cbind(convergence, dd)             # append to convergence
                                                    
  rrtg <- t                                         # next rrtg approximation
                                                    
  if (dd < steptol) break                           # change in rating below limit
}
NRend <- Sys.time()
# if (min(rrtg) < 1) rrtg <- rrtg - min(rrtg) + 1   # assure ratings are strictly positive (compatible with Sevilla)

if (!exists("report_tpr")) {
cat("rrtg lin", sep="\n"); print(rrtg, digits = 4)
cat("Iterations NR-lin", sep="\n"); print(iterations, digits = 4)
}
cat("Convergence NR-lin", sep="\n"); print(c(convergence), digits = 4)

# import extra player data in Sevilla (match=full name, data=Rating)
# write.csv2(cbind(rdtable[,2], W, rrtg), quote=FALSE)

# excluded players by SCC
c(which(SCC$membership != most_common_SCC))

{
cat("\n")
cat(sprintf("# spelers     = %5d      , # rondes = %d\n", npls,  nrds)      )
cat(sprintf("CG iterations = %5d (cum), reltol   = %g\n", cgcum,  cgeps) )
cat(sprintf("NR iterations = %5d      , maxit    = %#.1f\n", it, maxit) )
cat(sprintf("SD x-Diff     = %5.2g      , x_tol    = %5.2g\n", sd(diff, na.rm=TRUE), steptol) )
cat(sprintf("L2 Residual   = %5.2g\n" , sqrt(sum(res*res)) )) 
cat(sprintf("SD Residual   = %5.2g\n" , sd(res, na.rm=TRUE) ))
cat(sprintf("XN Residual   = %5.2g\n" , max(abs(res), na.rm=TRUE) )) # max norm 
cat("\n")                                                                                           
}

NRend - NRstart                                     # CPU NR iteration

uz   <- exp(rrtg / 400 * log(10))                   # Spielstärken, Zermelo, p. 451
uz[] <- uz / sum(uz)                                # genormeerd op som = 1

