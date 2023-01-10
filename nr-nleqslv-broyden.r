# ---------------------------------------------------------------------
# nr-nleqslv.r
# Solving systems of nonlinear equations with Broyden
# https://cran.r-project.org/web/packages/nleqslv/nleqslv.pdf
# calculate relative ratings (x) as roots of We(x) - W = 0
# ---------------------------------------------------------------------
# timing packages
# function      user  system elapsed  Tournament  
# nr-nleqslv.r  4.77    0.00    4.77  JBFhuge 
# nr            0.03    0.00    0.03  JBFhuge
#
# install.packages("nleqslv")
# Run gf.r first
# list.files(getwd(), pattern = "*.r")
# source('nr-nleqslv.r')

require(nleqslv)
ls("package:nleqslv")
require(igraph)

celo <- 400/log(10)                                 # Elo constant in Logistic distribution
E  <- function(DP) {return(1./(1. + exp(-DP)))}     # logistic function
dE <- function(DP) {e <- E(DP); return(e * (1-e))}  # derivative

#--------------------------------------------
# find roots of f(x) = We(x) - W 
# f(x) = ( f1(x), f2(x) ... fn(x), n = number of players
# fi(x) = Sum(We(xi - xj) - W), j = 1..rounds
# fn+1 =  ∑x
# Par = max(∀ plusscore - minscore)
# -------------------------------------------
f <- function(x) {
   x <- head(x, -1) + mask_rtg                                # remove non_SCC and set dimension

  dpopp <- x[,1] + mask_rtg[,1] - x[as.vector(opponents)]  
  dim(dpopp) = dim(opponents)    

  Wem <- (E(dpopp / celo) - 0.5)                    # expected score above average(n x r)
  stopifnot(sum(Wem, na.rm=TRUE) < .Machine$double.eps * length(Wem)) # sum of rating differences
  We <- matrix(rowSums(Wem, na.rm=TRUE)) * Par      # expected score (n)
  return(c(We-W, sum(x)))                           # residual = f(x)
}

#----------------------------------------------------
# choose largest indivisible profile domain
# calculate mask_SCC (n x r), NA for excluded players
# recalculate wins, losses, draws within main domain
# test sum wins and losses (otherwise nr will fail)
# ---------------------------------------------------

# column vector of relative ratings
# stopifnot(SCC$no==1)                              # OK iff strongly connected 

mask_rtg <- matrix(ifelse(SCC$membership %in% largest_SCC, 0, NA)) # base is 0, exclude largest SCC
W <- matrix(rowSums(results + mask_rtg[opponents] + mask_rtg[,1] , na.rm=TRUE)) # "Copeland" score points (above/below average), within SCC largest.
stopifnot(sum(W) < .Machine$double.eps)             # zero sum
maxit <- npls * log(npls) + 5                       # maximum number of nr iterations, 5 for small npls

p0 <- rep(0, nrow(W) + 1)                           # ratings, ∑ratings

stime <- system.time(                                                
sol <- nleqslv (
                , x  = p0                           # start at x=0
                , fn = f                            # root of W - We(x), x = rating difference
                , method = "Broyden"                # Newton requires jacobian
                , control=list(allowSingular = TRUE, ftol=1e-6)
                )
)
#---------------------------------------------------------------------
rrtg_sol <- head(cbind(sol$x), -1) + mask_rtg        # remove sum and mask main SCC
print(zapsmall(rrtg_sol, digits=4))                  # relative rating broyden

uz   <- exp(rrtg_sol / 400 * log(10))                # Spielstärken, Zermelo, p. 451
uz[] <- uz / sum(uz)                                 # genormeerd op som = 1

print(sol)
print(rrtg_sol)
cat(sprintf("Stdev f$vec = %5.2g\n\n"  , sd(sol$fvec) ) )
print(stime)

