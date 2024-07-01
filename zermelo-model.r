# zermelo-model.r
# Unlicense (Ø) 2021 CLP, IJmuiden
# Zermelo, E. "Die Berechnung der Turnier-Ergebnisse als ein Maximumproblem der Wahrscheinlichkeitsrechnung.."
# Mathematische Zeitschrift 29 (1929): 436-460. <http://eudml.org/doc/168081>.
# https://gdz.sub.uni-goettingen.de/id/PPN266833020_0029
# https://www.digizeitschriften.de/download/pdf/266833020_0029/log31.pdf
#
# (p.437) Unser verfahren kommt darauf hinaus, daß die ralativen Spielststärken
# als Wahscheinlichkeiten aufgefaßt und so bestimmt werden, daß die
# Wahrscheinligkeit für das Eintreten des beobachteten Turnier-Ergebnisses
# eine möchlichst große wird. 
#
# Test example
# rElo 104.008415, 104.008415, 4.855994, -212.872824
# npls      <- 4
# opponents <- matrix(c(4,3,2,1,3,4,1,2,2,1,4,3), nrow=npls) # indexed by player, round.
# gfile     <- matrix(c(2,1,1,0,1,2,1,0,1,1,1,1), nrow=npls) # indexed by player, round.
# Par       <- 2
# ------------------------------------------------------------
                                            
if (!exists("SCC") || SCC$no > 1 ) q()          # Single SCC.
if (!exists("npls")) q()                        # Run gf.r first

# U is vector of ratings indexed by player.
# Probability(Ar beats As),
#   Ur,s = Ur / (Ur + Us), p. 437. Eq. 1.
# Combined probability of all games,
#   W    =  ∏ Wr,s ∀r,s and r ≠ s. Eq. 3.
# Combined probabity of Ar beats As with Gr,s against Gs,r,
#   Wr,s = (Ur,s ^ Gr,s) * (Us,r ^ Gs,r). Eq 2, binomial model.
# Gr,s is number of half wins player Ar against As.
# Kr,s   = Gr,s + Gs,r, no draws.
#
# Function Φ, eq. 3, p.438.
# Returns log of combined probability (W).
fie <- function(u) {
  u <- matrix(u, nrow = npls)
  uropp <- u[, 1] / (u[, 1] + u[as.vector(opponents)]) # as.vector for two rounds.
  dim(uropp) <- dim(opponents)
  return(sum(gfile * log(uropp), na.rm = TRUE)) # gfile = number of halfwins,
}                                               # upto a multiplicative constant.
u0 <- rep(1, npls)                              # All players have equal strength.

# # p. 438, eq. (3).
# Maximize combined probability (wz) of all outcomes.
wz <- optim(par = u0
            , fn = fie
            , gr = NULL
            , method = "L-BFGS-B"               # Low memory, bounded, L-BFGS-B.
            , lower = 1E-7, upper = 20          # Boxed. (Mannheim)
            , control = list(maxit = 1000, pgtol = 0, fnscale = -1))
c(wz$value, fie(u0))
wz$message
if ( abs(fie(wz$par) - wz$value) > .Machine$double.eps) warning("Something wrong in Φ, optim")
if (fie(wz$par) < fie(u0)) warning("Solution less likely then all draws")

# Convert to Elo domain.
rElo <- log(wz$par) * (400 / log(10))
rElo <- rElo - mean(rElo)                       # Normalize at sum is zero.
round(rElo)
                                   
# ---------------------------------------------------
# Expected result Ar against As.
#   ΣKrt * Urt, ∀t, r = 1,2,...,n. Eq 8a, p. 443.
# Krt is number of games between Ar and At.
# In the gamefile representation Krt = 1 for all players. One game per round.
# The crux: validate actual score (Waz) equals expected score Wez(u).
Waz <- rowSums(gfile, na.rm=TRUE)               # Actual score from game file.

uz <- matrix(wz$par, nrow = npls)               # Most likely ratings.
uzopp <- uz[, 1] / (uz[, 1] + uz[as.vector(opponents)]) # Indexed by player, round.
dim(uzopp) <- dim(opponents)                    # Likely result.
Wez <- rowSums(uzopp, na.rm = TRUE) * Par       # Expected score.
stopifnot(sd(Waz - Wez) < 1E-3)                 # Actual score equals expected score.                                         
sd(Waz - Wez)

# --------------------------------------------------------------------------------------
# Calculate rating by Zermelo iteration
# p. 453, §5.
maxit <- max(npls * log(npls), 1000)             # Maximum number of iterations.
tol   <- 1E-3
uz <- matrix(u0, nrow=npls)

# Iterate until Waz == Wez.
for (it in seq(maxit)) {
  uzopp <- uz[, 1] / (uz[, 1] + uz[as.vector(opponents)]) # Indexed by player, round.
  dim(uzopp) <- dim(opponents)
  Wez <- rowSums(uzopp, na.rm = TRUE) * Par     # Expected score.
  if (any(Wez < .Machine$double.eps)) break     # Strictly positive.
  if (sd(Waz - Wez) < tol) break
        
  uz.n <- (Waz / Wez) * uz                      # Next rating.
  uz.n <- uz.n / exp(sum(log(uz.n)) / length(uz.n)) # Renormalize to product is one.
  if (sd(uz - uz.n) < 1E-5) {
    print("No progression in iteration step")
    print(sd(uz - uz.n))
    break                                       # No convergence
  }
  uz <- uz.n
}

if (it >= maxit) print("Maximum number of iteration reached.")
sd(Waz - Wez)
c(fie(uz), fie(wz$par))

# Convert to Elo domain.
rEloT <- log(uz) * (400 / log(10))
rEloT <- rEloT - mean(rEloT)                    # Normalize at sum is zero.
# zapsmall(rEloT)                               # Calculated by iteration.
