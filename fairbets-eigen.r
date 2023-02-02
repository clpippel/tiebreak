# compute fair bets tiebreak 
# R eigen function, check fairbets.r
# Paired comparisons analysis:  an axiomatic approach to ranking methods
# http://eio.usc.es/pub/julio/papers/15._Paired_comparisons_analysis-_an_axiomatic_approach_to_rankings_in_tournaments_web.pdf 
#---------------------------------------------------

# Crosstable results
# Ranking Participants in Generalized Tournaments
# Giora Slutzki ∗ and Oscar Volij
# September 8, 2004
# Ctab =
# matrix(
# c(0, 1, 1, 0,
#  0, 0, 1, 1,
#  0, 0, 0, 1,
#  1, 0, 0, 0
# ), nrow = 4, byrow=TRUE)

Ctab <- as.matrix(g[])                              # Cross table results (indexed by player)

# Fairbets matrix
Ftab <- Ctab; diag(Ftab) <- 1                       # against self, avoid division by zero
Ftab <- Ftab / colSums(Ftab)                        # Fairbets tab

Xtab <- Ctab                                        # Ctab <<<<<<<<<<<<<<<<<<<< solve eigen values/vectors

# compute principle eigenvector
system.time(eigenv <- eigen(Xtab))
mev <- which(eigenv$value == max(Re(eigenv$value[Im(eigenv$value)==0L])) ) # index of maximum real eigenvalue
pevalR <- setNames(Re(eigenv$value[mev]), "Principle real eigen value")
pevecR <- (Re(eigenv$vectors[, mev]) )              # principal real eigenvector
pevecR <- pevecR / pevecR[tail(which(round(pevecR, digits=7)!=0),1)] # normalize minimum (>0) = 1
pevecR <- matrix(pevecR, dimnames= list(c(), "Principle real eigenvector"))
print(zapsmall(pevalR))
print(zapsmall(pevecR))

# check Av = λv
Av <- Xtab %*% pevecR                               # matrix time eigenvector
print(zapsmall(unique(matrix(round(Av / pevecR, 10), dimnames=list(c(), "A.v / λ.v") ))))  # show eigen values
print(intToUtf8(0x03BB))

