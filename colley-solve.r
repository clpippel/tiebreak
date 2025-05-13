# Colley's Bias Free College Football Ranking Method
# create crossref table and Laplacian matrix  indexed by player
# Colley Matrix is equivalent to the ''Generalized row sum'' method, a parametric family of ranking methods
# developed by P.Yu. Cheboratev (1989).
# <ref>
# Csató, L. On the ranking of a Swiss system chess team tournament,
# Annals of Operations Research, 254, 17-36 (2017).
# https://arxiv.org/abs/1507.05045v5.
# https://www.colleyrankings.com/matrate.pdf
#
# 12-05-2024, Check with equation 24
#             L400: ranking in Elo domain [-400, 400]          

rm(crossref, Lmat)
idx      <- data.frame(x = c(row(opponents)), y = c(opponents), z = c(res_sym))
crossref <- tapply(idx$z, idx[1:2], FUN = length)
cbind(rdtable["Name"], crossref)

# Create Laplacian
N <- rowSums(crossref, na.rm=TRUE)
Lmat <- -as.numeric(as.matrix(crossref))
dim(Lmat) <- dim(crossref)
Lmat[is.na(Lmat)] <- 0 
diag(Lmat) <- N
stopifnot(sum(Lmat)==0L)

# Colley matrix
# Mean ranking is .5
elon       <- 2;
b          <- (wins - losses) /2 + 1
Amat       <- Lmat
diag(Amat) <- diag(Amat) + elon
system.time(r.colley <- matrix(solve(Amat, b)))     # A.x = b

## Colley model
print(gfheader)
cbind(rdtable[,gfrmcols]
              , N, wins, draws, losses
              , Rk   = rank(r.colley, ties.method= "min")
              , R    = r.colley
              , L400 = round((r.colley - .5) * 800)
) -> colley
print(colley, digits = 4)

## source("scm_gain.r")

## https://www.arndt-bruenner.de/mathe/scripts/gleichungssysteme.htm
## https://arxiv.org/abs/2005.02280
##  5  0 -1 -1 -1   .5
##  0  4 -1  0 -1   1
## -1 -1  6 -1 -1   1
## -1  0 -1  4  0   1
## -1 -1 -1  0  5  1.5
##  
## c(19/46, 12/23, 1/2, 11/23, 27/46)  
## c(0.4130 0.5217 0.5000 0.4783 0.5870)

