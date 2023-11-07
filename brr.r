## brr.r
## http://www.jellyjuke.com/mathematical-explanation-of-the-bayesian-resume-rating.html
## Bayesian Resume Rating (Matt).
##--------------------------------------------------------------------------------
## https://github.com/clpippel/tiebreak (repository).
## run gf.r first to read game file input.
## Ctrl L to clear window.
## -------------------------------------------------------------------------------

cat('\n', paste("----------------------------- brr.r", sep=", "), '\n')

## brr parameters, default Jelly Juke
parity           <- 1                               #
ownSigmaIncluded <- 0                               # 0 = No 
oppSigmaIncluded <- 1                               # 1 = Yes
priorIncluded    <- TRUE                            # Bayesian prior
fixedS           <- FALSE                           # Set S = 1

                                   # run gf.r first
if ( !exists("results") )   q()    # indexed by player, round (npls x nrds), skew symmetric scores
if ( !exists("opponents") ) q()    # indexed by player, round (npls x nrds), opponents
if ( !exists("SCC") )       q()    # stronly connected components
npls <- nrow(opponents)
nrds <- ncol(opponents)
## if ( SCC$no != 1) q()           # Test results are strongly connected
sprintf("# Players = %d, # Rounds = %d, # SCCs = %d", npls, nrds, SCC$no)

##-----------------------------------------------------
## Calculate mask_SCC (n x r).
## Universal domain: set all ratings to 0.
##   or
## Choose largest indivisible profile domain and
## set rating of non main players to NA.
## ----------------------------------------------------
most_common_SCC <- which.max(table(SCC$membership))
brrB  <- matrix(ifelse(SCC$membership %in% most_common_SCC, 0, 0)) # column vector of BR Ratings (B), single domain
colnames(brrB) <- "B"

brrS <- cbind(rep(1, npls)); colnames(brrS) <- "S"
brrS[which(is.na(brrB))] <- NA                       # column vector of BRR Sigma's

iterations <- brrB                                   # store brrB, brrS
iterations <- cbind(iterations, brrS)                # indexed by player, iteration (2x)
convergence <- c(Inf)                                # concatenate missing value

limit <- (1.0e-4 / 2)                                # convergence after <<<four>>> decimals < .00005
maxit <- max(npls * nrds, 100L)                      # maximum number of iterations

##------------------------------------------
## calculate brrB, brrS.
##-----------------------------------------
for (it in seq(maxit)){

  oppS <- brrS[c(opponents)]                         # opponents sigma's, indexed by player x round
  ownS <- matrix(rep(brrS,each=nrds), ncol=nrds, byrow=TRUE) # own S, indexed by player x round
  gfS  <- sqrt(ownSigmaIncluded*ownS^2 + parity + oppSigmaIncluded*oppS^2 + parity) # S per game, indexed by player, round

  ## For each player p
  ## find distribution fp R --> R
  brrB_next <- brrB
  brrS_next <- brrS

  for (pp in seq(npls) ) {
    if (is.na(brrB[pp])) break
    
    ## Try different options for B and calculate most probable B.
    ## Find f R --> R.
    ## Discrete approximation of integral.
    ## Note: f is not a probability distribution,
    ## because the area under the curve ∫f does not add up to 1.
    ## Assume: f / ∫f(x)dx is normally distributed.
    ## Approximate integral by discrete steps.
    step <- 0.5
    approx <- seq(from = -5, to = 5, by = step)
    brrB_try <- brrB

    gfPr <- matrix(NA, 1, nrds)                     # probabilities win, loss, draw (1 x r)
    fofx <- rep(0, length(approx))                  # store resulting probability here   
    for (x in seq_along(approx) ) {
      ## Calculate joint probability of results of player pp,
      ## given brr = x.
      ## Pr(win) : N(dp, 0, sd)
      ## Pr(loss): 1 - Pr(win)
      ## Pr(draw): Pr(Win) * Pr(loss)
      brrB_try[pp] <- approx[x]
      
      oppopp     <- opponents[pp,]      
      rtopp      <- rbind(brrB_try[oppopp])         # rating opponents
      dpopp      <- approx[x] - rtopp               # rating difference (dp) with opponents
      
      ## Assuming brrB_try, calculate probability of actual results.
      ## Select rating differences when win, loss, draw.
      ## draw is modelled as a win followed by a loss.
      prior <- ifelse(priorIncluded, dnorm(approx[x], sd=1), 1) # Bayesian prior 
      gfPrwin <- pnorm(dpopp, 0, gfS[pp,])          # keep for efficiency
      resultspp                  <- results[pp,]
      gfPr[which(resultspp >  0)] <- gfPrwin[which(resultspp > 0)]
      gfPr[which(resultspp <  0)] <- 1 -  gfPrwin[which(resultspp < 0)]
      gfPr[which(resultspp == 0)] <- sqrt(gfPrwin[which(resultspp == 0)] * (1 - gfPrwin[which(resultspp == 0)]))
      fofx[x] <- exp(rowSums(log(gfPr), na.rm = TRUE) + log(prior) )
    }
    ## calculate next B, S
    sumf <- sum(fofx, na.rm = TRUE)
    expf <- sum(approx * fofx, na.rm = TRUE)        # expectation fofx: Σx.f(x).
    stopifnot(sumf * step > 0)                      # test area under the curve.
    brrBpp        <- expf /sumf 
    brrB_next[pp] <- brrBpp
    m2B <- sum((approx - brrBpp)^2 * fofx, na.rm = TRUE) # second moment location = B
    stopifnot(m2B > 0)                              # S is strictly positive
    brrS_next[pp] <- sqrt(m2B/sumf)                 # S
    if (fixedS) brrS_next[pp] <- 1                  # S is constant
    if (m2B == 0) break                             # for debugging
    ## assuming a normal distribution, derive sigma from max (alternative)
    ## autc <- sumf * step / max(fofx)              # area under the curve
    ## brrS_next[pp] <- autc / sqrt(2 * pi)         # sigma.
 } 
 iterations  <- cbind(iterations, brrB_next, brrS_next) # append next iteration
 dd <- 3 * sd(brrB_next - brrB, na.rm=TRUE)         # 3σ, 99.7%
 xdd <- brrB - brrB_next
 if (is.na(dd)) dd <- -1                            # single player
 convergence <- cbind(convergence, dd)              # append to convergence
                                                    
 brrB <- brrB_next                                  # next brrB approximation
 brrS <- brrS_next                                  # Sigma of f

 if (sum(convergence == dd, na.rm=TRUE) > 2) break  # cycle, no convergence
 if (dd < limit) break                              # change in rating below limit
}

brrN <- brrB/brrS; colnames(brrN) <- "B-norm"
rk  <- cbind(rank(-brrN, ties.method= "min" )); colnames(rk) <- "Rank"
iterations <- cbind(iterations, brrN, rk)
rm(rk)

{
cat(sprintf("Parity = %.1f, With Own σ = %d, With Opp σ = %d, With Prior = %d\n", parity, ownSigmaIncluded, oppSigmaIncluded, priorIncluded) )
print("Iterations, (base = 0)"); print(round(iterations, 4) )
print(convergence)
}

{
cat(gfheader, sep='\n')
cat(sprintf("players    = %d, rounds      = %d\n", npls, nrds) )
cat(sprintf("dd         = %g, conv. limit = %g\n", dd, limit) )
cat(sprintf("iterations = %d, maxit       = %g\n", it, maxit) )
if(dd >= limit) warning("No convergence <<<<<")
}
