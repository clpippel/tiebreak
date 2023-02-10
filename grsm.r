# Run gf.r first to process input
# Generalized Row Sum Method (Pavel Chebotarev)
# define LSM when alternative model is required (LSM = 1,10,800 for LSM, ECF, AVG models)
#
# On the ranking of a Swiss system chess team tournament, 2015, https://arxiv.org/abs/1507.05045v5
# A graph interpretation of the least squares ranking method (https://arxiv.org/abs/1508.06778)
# Preference fusion when the number of alternatives exceeds two: indirect scoring procedures, 1999), https://arxiv.org/abs/math/0602171v3
# Generalization of the Row Sum Method for Incomplete Paired Comparisons,” 1989, Automation and Remote Control, 50, 1103-1113,
# https://www.researchgate.net/publication/258514016_Generalization_of_the_Row_Sum_Method_for_Incomplete_Paired_Comparisons_Automation_and_Remote_Control_50_1103-1113
# Generalized Row Sums (x) are a solution of:
#   (I + εL)x(ε)(N, R, M) = (1 + εmn)s(N, R, M), where ε > 0, ε ≤ 1 / (m(n - 2)) or
#   (e'.I + L)x = y.s, y = e' + m.n, e' > m.n - 2, e' = 1 /e (epsilon in above pubs)
#
# ECF  grading (by game) : Opponent's grade - 50  + 100n, n = (1, ½ or 0)
# USCF initial rating    : Opponent's grade - 400 + 800n, n = (1, ½ or 0)
# CA   Rp = Rc + 400 (W - L) / N (https://chess.ca/node/444)
#
# Conjugate gradients, solve A.x = b
# Sparse representation. A, opponents indexed by player x round
# A must be Positive Definite or
# semi PD and b in column space of A (Kaasschieter)
# terminology in:
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.newton_krylov.html#scipy.optimize.newton_krylov
# 19-01-2023, LSM naar elon
# 10-02-2023, ys

maxlines = 500                                      # max players to print

cat( '\n', paste("----------------------------- grsm.r", "LSM", exists("LSM"), sep=", "), '\n' )

cg.solve  <- function(A, Adiag, b, reltol=NULL) {
  # Conjugate gradient method
  # A                                               # matrix
  # Adiag                                           # A[i,i]
  # opponents                                       # sparse matrix representation
  # reltol = 

  if (is.null(reltol)) reltol <- sqrt(.Machine$double.eps)
  
  if (length(b) == 0) return(b)
  x <- matrix(0, length(b), 1)                      # we assume x0 = 0
  r <- b                                            # Residual, r = b - A.x; x = 0
  p <- r                                            # search direction

  rsold <- sum(r * r)                               # r.r
  rs0   <- rsold                                    # initial error

  for (it in seq_along(b)) {
    if (rsold < rs0 * reltol * reltol) break        # found solution
    # as.vector(), otherwise it won't work for two dimensions
    Ap <- rowSums(A * p[as.vector(opponents)], na.rm = TRUE) # Ap <- A.p
    Ap <- Ap + Adiag * p                            # sparse file multiplication
    pAp <- sum(p * Ap)                              # A-orthogonal 
    if (pAp < .Machine$double.eps) break            # search directions exhausted

    alpha <- rsold / pAp
    x <- x + alpha * p
    r <- r - alpha * Ap                             # cumulative update (b - Ax)
    rsnew <- sum(r * r)                             # r.r next
    cgcum <<- cgcum + 1                             # count iterations (global assignment)
    p = r + (rsnew / rsold) * p                     # next search direction
    rsold = rsnew
  }
  # browser()                                       # open browser environment to debug or inspect locals
  return(x)
}

stopifnot(npls > 1)                                 # two players
#----------------------------------------------------
nrrr = max(unlist(apply(opponents,1, table)))       # (m) number of round robin rounds (max same opponent)
s    = matrix(rowSums(results, na.rm=TRUE))         # saldo wins and losses, skew symmetric score
growsums <- matrix(0, nrow(opponents), 1)           # Generalized Row Sums

if (exists("LSM")) {                                # choose LSM or GRS
    elon <- LSM; y <- 1 ;                           # ε = 0, γ = scale factor
} else {
    elon <- (npls-2) * nrrr;                        # elon = round(log(npls)), approx Elo
	# elon <- (npls-1) * nrrr;                      # cheb 1989 example 1
    y    <- elon + npls * nrrr;                     # GRS: y = ε' + m.n,  well chosen: ε' ≥ n.m - 2 (Gonzalez-Diaz) 
}

# ---------------------------------------------------
# Original equation Cheboratev
# solve x | (e'.I + L)x = y.s
# e'is parameter, y = e'+ m.n, n = players, m = max number 
# L is Laplacian of the adjacency game graph
# s is skew symmetric score
# Lmat = (e'.I + L)
# ---------------------------------------------------
Lmat     <- (-!is.na(opponents)) * 1                # Gamefile: Match matrix M = 1
Ldiag <- -rowSums(Lmat) + elon                      # diag(L) + e'.I
kappa = max(Ldiag) / min(Ldiag[Ldiag[] >0 ])        # λmax /  λmin 
cg_th <- round(sqrt(kappa) * log(2 / sqrt(.Machine$double.eps) ) / 2 +0.5)  # rate of convergence CG method

ys <- y * s                                         # intercept
cgcum <- 0                                          # statistics
growsums <- cg.solve(Lmat, Ldiag, ys)               # solve lineair equation
# ---------------------------------------------------

{
message("Solve: x | (e'.I + L)x = y.s where y = e'+ n.m, L = Laplacian of game graph")
message(sprintf("e' = %a, y = %a, n = %d, RR-rounds = %d", elon, y, npls, nrrr) )  

# validate solution
x <- growsums
Ax <- rowSums(Lmat * x[as.vector(opponents)], na.rm = TRUE) # A.x
Ax <- Ax + Ldiag * x								# sparse matrix multiplication
res <- Ax - ys 										# A.x - y.s

cat("\nValidate grsm solution\n")
cat(sprintf("L2 Residual   = %5.2g\n" , sqrt(sum(res*res)) )) # Euclidean norm L2
cat(sprintf("SD Residual   = %5.2g\n" , sd(res, na.rm=TRUE) )) # standard deviation
cat(sprintf("XN Residual   = %5.2g\n" , max(abs(res), na.rm=TRUE) )) # max norm 
cat("\n")
}

# e' is well choosen:
# Generalization_of_the_Row_Sum_Method_for_Incomplete_Paired_Comparisons, 1989
# 4. Generalized Row Sums and Paired Comparisons with bounded outcomes (eq.21)
# Expected outcome between i, k in interval [Rmax, Rmin]
gsopp     = growsums[as.vector(opponents)]          # rating opponents
dim(gsopp) = dim(opponents)                         # restore matrix (n x r)
dgrs <- ( growsums[,1]) - gsopp                     # difference growsums player and opponents (xi - xj)
rmax  <- max(results, na.rm=TRUE)                   # maximal win

mask_rmax  <- ifelse(results == rmax, 0, NA)        # wins, maximal win, win / loss are symmetric

# fij = results + (xj - xi + results.m.n) / elon    # 1 / elon = Eps (1989, eq 20)
# fij = (results + (xj - xi)/y) * y/elon
# payoff function f must be non-negative/positive for rmax/rmin (1989, theorem 3, eq 20)
{
  gain <- results + mask_rmax - dgrs / y            # gain of max result must be positive (loon naar werken)
  if ( TRUE %in% (min(gain, na.rm=TRUE) < 0) ) { 
	cat(sprintf("Negative gain in rmax: %g\n", rmax))
    if (!exists("report_tpr")) print(ifelse(gain > 0, NA, gain), na.print = "." , quote = FALSE, digits = 2)
  }
  else {
	cat(sprintf("Monotone, max gain : %g\n", max(gain, na.rm=TRUE)))
    sc_monotone <- TRUE
  }
}

# show rank, ratings
rkgrs <- cbind(rank(-growsums, ties.method= "min"), growsums)
colnames(rkgrs) <- c("Rk", ifelse(exists("LSM"), "LSSums", "GRSums"))
{ 
cat(sprintf("GRS parameters: e=%.2f, y=%.2f\n", elon, y))
cat(sprintf("CG convergence: kappa=%g, cg_th=%d, # of iterations=%g\n", kappa, cg_th, cgcum))
if (!exists("report_tpr")) print(rkgrs[1:(min(nrow(rkgrs),maxlines)),], digits=3, noquotes=TRUE)
}
 
dev.new()
try({
lytgrs     <- layout_with_fr(g)
# lytgrs     <- layout.circle(g)
lytgrs[,2] <- growsums
plot(g, layout=lytgrs, main = paste0("SCCs with ", ifelse(exists("LSM"),"LSM", "GRS"), " layout") )
})

growsums[which(apply(is.na(opponents), 1, all))] <- NA    # remove isolates <<<<<<

# in [-Par, Par]
# print(matrix(outer(growsums, growsums, FUN="-")/max(elon,1), nrow=length(growsums)), digits=2)
