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

# Run gf.r first
if (exists("SCC")) stopifnot(SCC$no == 1)       # Single SCC.

# U is vector of ratings indexed by player.
# Probability Ar beats As: Ur,s.
#   Ur,s = Ur / (Ur + Us), p. 437. Eq. 1.
# Combined probability of all games.
#   W    =  ∏ Wr,s ∀r,s and r ≠ s. Eq. 3.
# Combined probabity of Ar beats As with Gr,s against Gs,r.
#   Wr,s = (Ur,s ^ Gr,s) * (Us,r ^ Gs,r). Eq 2, binomial model.
# Gr,s is number of half wins player Ar against As.
# Kr,s   = Gr,s + Gs,r, no draws.
#
# Returns log of combined probability (W).
fie <- function(u) {
  u <- matrix(u, nrow = npls)
  uropp <- u[, 1] / (u[, 1] + u[as.vector(opponents)]) # as.vector for two rounds.
  dim(uropp) <- dim(opponents)
  return(sum(gfile * log(uropp), na.rm = TRUE)) # gfile = number of halfwins,
}                                               # upto a multiplicative constant.
u0 <- rep(1, npls)                              # All players have equal strength.
fie(u0)

# # p. 438, eq. (3).
# Minimize combined probability (wz) of all outcomes.
wz <- optim(par = u0
            , fn = fie
            , gr = NULL
            , method = "L-BFGS-B"                # Low memory, bounded, L-BFGS-B.
            , lower = 1E-7, upper = 20           # Boxed.
            , control = list(maxit = 1000, pgtol = 1E-10, fnscale = -1))
wz$message

# Convert to Elo domain.
rElo <- log(wz$par) * (400 / log(10))
rElo <- round(rElo - mean(rElo))                 # Normalize at sum is zero.
rElo

# ---------------------------------------------------
# Expected result Ar against As.
#   ΣKrt * Urt, ∀t, r = 1,2,...,n. Eq 8a, p. 443.
# Krt is number of games between Ar and At.
# In the gamefile representation Krt = 1 for all players. One game per round.
# The crux: validate actual score (Waz) equals expected score Wez(u).
Waz <- rowSums(gfile, na.rm=TRUE)                 # Actual score from game file.

uz <- matrix(wz$par, nrow = npls)                 # Most likely ratings.
uzopp <- uz[, 1] / (uz[, 1] + uz[as.vector(opponents)]) # Indexed by player, round.
dim(uzopp) <- dim(opponents)                      # Likely result.
Wez <- rowSums(uzopp, na.rm = TRUE) * Par
stopifnot(sd(Waz - Wez) < 1E-3)                   # Actual score equals expected score.
sd(Waz - Wez)
