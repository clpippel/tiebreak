# create report of tiebreak methods
# Run gf.r first to process input
# 19-01-2023, Colley matrix, Laplace rule of succession

report_tpr <- TRUE									# Signal reporting

if (exists("LSM")) rm(LSM)

source('nr-lin.r', echo=FALSE, print=TRUE)          # relative AVG ratings, newton iteration
avg400  <- rrtg                                     # P(x) = ½ + avg400 / 800   ; AVG400 probability

source('nr.r', echo=FALSE)                          # relative Elo ratings, newton iteration

source('grsm.r', echo=FALSE)                        # Generalised Row Sum Ratings
grs     <- growsums
colnames(grs) <- paste0("grs", elon)
elongrs <- elon
ygrs    <- y                                        # P(x) = ½ + grs / ygrs / par * 2; grs probability

source('recursive-bhz.r', echo=FALSE, print=TRUE)   # Recursive Buchholz

LSM <- 0                                            # Least Square Ratings, add LSM to diagonal Laplacian matrix
source('grsm.r', echo=FALSE)
lsq     <- growsums                                 # P(x) = ½ + lsq / par / 2  ; lsq probability
elonlsq <- elon
ylsq    <- y

# https://en.wikipedia.org/wiki/Rule_of_succession  # Laplace: (s + n) / N + 2)
# https://www.colleyrankings.com/matrate.pdf
LSM <- 2                                            # Colleys matrix, ε = 2
source('grsm.r', echo=FALSE)
Lapl <- growsums + Par / 2                          # recalibrate to draw (Colley matrix)
elonLapl <- elon
yLapl    <- y

colnames(Lapl)   <- "Lapl"
colnames(avg400) <- "A400"
colnames(lsq)    <- "lsq"

stopifnot( sd(avg400 - (lsq/Par)*800, na.rm = TRUE) < 1E5) # lsq ~ AVG400

pefbev <- 1                                          # Perron Frobenius (principle) eigen vector
source('fairbets.r', echo=FALSE, print=TRUE)
pev <- fb * sum(points)                             # https://www.arndt-bruenner.de/mathe/scripts/engl_eigenwert2.htm

rm(pefbev)                                          # fair bets
source('fairbets.r', echo=FALSE, print=TRUE)
fbpts <- fb * sum(points)                           # normalize by all points (ping-pong model (Volij)

frk <- seq_len(npls)                                # final rank
rank1 <- matrix(rank(-lsq   , ties.method= "min")); colnames(rank1) <- "  Rk"
rank2 <- matrix(rank(-rbhz,   ties.method= "min")); colnames(rank2) <- "  Rk"
rank3 <- matrix(rank(-avg400, ties.method= "min")); colnames(rank3) <- "  Rk"
rank4 <- matrix(rank(-rrtg,   ties.method= "min")); colnames(rank4) <- "  Rk"
rank5 <- matrix(rank(-Lapl,   ties.method= "min")); colnames(rank5) <- "  Rk"
rank6 <- matrix(rank(-grs,    ties.method= "min")); colnames(rank6) <- "  Rk"
rank7 <- matrix(rank(-pev,    ties.method= "min")); colnames(rank7) <- "  Rk"
rank8 <- matrix(rank(-fbpts,  ties.method= "min")); colnames(rank8) <- "  Rk"

rrtgf <- rrtg
Roffset <- NA
if (length(raticol) == 1L) Roffset <- mean(as.numeric(sub("[ABC] ", "", rdtable[,raticol]) ), na.rm =TRUE)
if (is.finite(Roffset)) rrtgf <- rrtgf + Roffset

N                <- wins + draws + losses
pctW             <- W / (N * Par)

colnames(points) <- "Pts"
colnames(N)      <- "N"
colnames(wins)   <- "+"
colnames(draws)  <- "="
colnames(losses) <- "-"
colnames(s)      <- "t/2"                           # t = Copeland-index, t / 2 == s == rowSums(R)
colnames(W)      <- paste0(intToUtf8(177), "W")     # plus-minus sign compatible with Ansi
colnames(pctW)   <- paste0(intToUtf8(177), "%W")    # plus-minus sign compatible with Ansi
report1 <- matrix(round(lsq   ,3)); colnames(report1) <- "lsq"
report2 <- matrix(round(rbhz  ,3)); colnames(report2) <- "RBhz"
report3 <- matrix(round(avg400,0)); colnames(report3) <- "A400"
report4 <- matrix(round(rrtgf))   ; colnames(report4) <- "rElo"
report5 <- matrix(round(Lapl  ,3)); colnames(report5) <- "Lapl"
report6 <- matrix(round(grs   ,3)); colnames(report6) <- paste0("grs",elongrs)
report7 <- matrix(round(pev   ,3)); colnames(report7) <- "Pev"
report8 <- matrix(round(fbpts ,3)); colnames(report8) <- "f-bets"
score   <- apply(as.character(as.hexmode(cbind(rank1, rank2, rank3, rank4, rank5, rank6, rank7, rank8)), width=1), 1, paste, collapse="") #drop 3,4 
score   <- ifelse(nchar(score) > min(nchar(score)), "", score) # single digit ranks
report  <- cbind(rdtable[,gfrmcols]
                , N, wins, draws, points, s, W, round(pctW, 2)
                , rank1, report1, rank2, report2, rank3, report3, rank4, report4, rank5, report5, rank6, report6, rank7, report7, rank8, report8, score)

# report[is.na(report)] <- Inf

# sort report by Points, rElo, AVG400
# report <- report[order(report[,"Pts"],report[,"rElo"],report[,"A400"],decreasing=c(TRUE, FALSE, FALSE)),]

crlf <- '\n'
cat(crlf, paste("----------------------------- report.r", sep=", ") )
{report1
cat(sprintf(crlf))
cat(sprintf("%s", fcsv), crlf)
cat(gfheader, sep='\n')
cat(sprintf("\nGamefile: %dx%d(nrrr = %d), Par = %g", npls, nrds, nrrr, Par), crlf)
cat(sprintf("lsq  : e'= %2d, y = %d", elonlsq, ylsq), crlf)
cat(sprintf("A400 : e'= %2d, y = %g", 0, 800 / Par), crlf)
cat(sprintf("Lapl : e'= %2d, y = %d", elonLapl, yLapl), crlf)
cat(sprintf("grs  : e'= %2d, y = %d (= e' + n.m)", elongrs, ygrs), crlf, crlf)

print(report, quote=FALSE, row.names = FALSE, , na.print = ".", max=200*ncol(report))
}
rm(LSM)
rm(report_tpr)

cat(crlf)
cat(sprintf("A400/800 -/- lsq/Par : %e (sdev)\n", sd(avg400 / 800 - lsq / Par, na.rm=TRUE)) )
cat(sprintf("Rbhz     -/- lsq     : %e (sdev)\n", sd(rbhz         - lsq      , na.rm=TRUE)) )

colnames(rrtg) <- "rElo"

# when round robin print probabilities, avg400, lsq, grs, s / N+nrrr, rElo,
# 
# probs <- cbind( avg400 / 800
#               , lsq    / Par
#               , grs    / ( (N + 1) * Par)
#               , s      / ( (N +1 ) * Par)
#               , 1      / ( 1 + exp(-rrtg /(400 / log(10))) ) - 0.5
#               )
# 
# colnames(probs) <- c("A400", "lsq", "grs",  "% t/2", "% rElo")
# writeLines("\nProbabilities")
# print(zapsmall(probs, digits = 4))

# {
# x_as <- seq_len(npls)
# dev.new()
# plot (x_as, probs[,1], main="Probs", type="l", lwd=3, col=rainbow(3)[1], xlab="Players")
# lines(x_as, probs[,3], type="l", lwd=3, col=rainbow(4)[2])
# lines(x_as, probs[,5], type="l", lwd=3, col=rainbow(4)[3])
# legend("topright", legend=c("avg400", "Grs", "rElo"), col=rainbow(4), lty=c(1,1,1), lwd=3)
# }

