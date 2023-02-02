# zermelo.r
# compute relative rating from gamefile (produced by Sevilla)
# relative rating, The Rating of Chessplayers, Past&Present, Arpad E. Elo,
# p54, The Method of Successive Approximations

if ( !exists("opponents") ) q()
if ( !exists("SCC") ) q()
# if ( SCC$no != 1) q()
#--------------------------------------------------------------------------------
# run gf.r to read input
# -------------------------------------------------------------------------------

#----------------------------------------------------
# choose largest indivisible profile domain
# set rrtg of excluded players to NA, otherwise 0
# calculate mask_SCC (n x r), NA for non main player
# recalculate wins, losses, draws within main domain
# test sum wins and losses (otherwise nr will fail)
# ---------------------------------------------------
most_common_SCC <- which.max(table(SCC$membership))
rrtg <- matrix(ifelse(SCC$membership %in% most_common_SCC, 0, NA)) # column vector of relative ratings
mask_SCC  <- rrtg[c(opponents)]                     # remove non-main opponents
mask_SCC  <- mask_SCC + rrtg[,1]                    # remove opponents of non-main players
dim(mask_SCC) <- dim(opponents)                     # restore matrix (n x r)

nr_wins     <- matrix(rowSums((results >  0) + mask_SCC , na.rm = TRUE)) # number of wins
nr_draws    <- matrix(rowSums((results == 0) + mask_SCC , na.rm = TRUE)) # number of draws
nr_losses   <- matrix(rowSums((results <  0) + mask_SCC , na.rm = TRUE)) # number of losses

W  = nr_wins + nr_draws / 2                         # Game points, actual result
games = nr_wins + nr_draws + nr_losses              # Number of games
We = games / 2                                      # expected score
iterations <- rrtg                                  # store iterations, one column for each iteration step
convergence <- c(NA)                                # concatenate missing value ( n x it. step)

npls <- nrow(games)
epsR  <- 1e-5                                       # test residual
limit <- (3 * 1.0e-3)^2 * npls                      # sigma(diff) / 3 < 1, change in diff after <<<three>>> decimals
maxit <- max(npls * npls, 100L)                     # maximum number of iterations
celo <- 400/log(10)                                 # Elo constant in Logistic distribution

#------------------------------------------
# find roots of f(x) = We(x) - W 
# f(x) = ( f1(x), f2(x) ... fn(x) )
# fi(x) = Σ(We(xi - xj) - W), j = 1..rounds
# -----------------------------------------
for (it in 1:maxit){

  rtopp     = rrtg[c(opponents)]                    # rating opponents
  dim(rtopp) = dim(opponents)                       # restore matrix (n x r)
  dpopp <- (rrtg[,1] - rtopp)                       # rating difference (dp) with opponents 
  dpopp < - dpopp + mask_SCC                        # remove non-main opponents

  Wem = 1 / (1 + exp(-dpopp/celo))                  # expected score = 1 / (1 + 10^(-dp/400)) (n x r) 
  We = matrix(rowSums(Wem, na.rm=TRUE))             # expected score (by player)
  res = We - W                                      # residual = f(x)
  if (sum(res*res) < .Machine$double.eps) break     # solution found

#---------------------------------------------------------
# Zermelo iteration, slow, stable (shutter weir, klepstuw)
# Tangent update, Newton-Raphson in one dimension
# --------------------------------------------------------
  it_step <- log(We/W)                              # Zermelo update
# it_step <- ((We-W)/games) * 4                     # Tangent update, 1/4 = tangent at E(dp) = 0
  t <- rrtg - (it_step * celo)
#----------------------------------------------------
                                                    
  tx = mean(t, na.rm = TRUE)                        # average rating, ignore missing values
  t <- t - tx                                       # reset average rating to zero
                                                    
  iterations  <- cbind(iterations, t)               # column bind, append next iteration
  diff <- t - rrtg                                  # change in iteration
  dd <- sum(diff*diff, na.rm=TRUE)                  # dot product diff.diff
  convergence <- cbind(convergence, dd)             # append to convergence
                                                    
  rrtg <- t                                         # next r approximation
                                                    
  if (dd < limit) break                             # change in rating below limit
}

if (min(rrtg) < 1) rrtg <- rrtg - min(rrtg) + 1     # assure ratings are strictly positive
rrtg <- round(rrtg)                                 # round to integer
rrtg <- round(rrtg)                                 # round to integer

print("Iterations, (base = 0)"); print(round(iterations, 2) )
print(convergence)

# import extra player data in Sevilla (match=full name, data=Rating)
write.csv2(cbind(rdtable[,2], rrtg), quote=FALSE)

sprintf("dd = %g, convergence limit = %g", dd, limit)
sprintf("iterations = %d, maxit = %g", it, maxit)


