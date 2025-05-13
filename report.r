# create report of tiebreak methods
# Run gf.r first to process input
# 19-01-2023, Colley matrix, Laplace rule of succession
#  6-02-2023, Recalibrate Lapl
# 10-02-2023, Colley Matrix, recalibrate
# 15-05-2023, score calculation corrected for R3.0.0
# 11-05-2025, Colley rankings in Elo domain [-400, +400]

report_tpr <- TRUE									# Signal reporting

source('nr-lin.r', echo=FALSE)                      # relative AVG ratings, newton iteration
avg400  <- rrtg                                     # P(x) = ½ + avg400 / 800   ; AVG400 probability

source('nr-elo.r', echo=FALSE)                      # relative Elo ratings, Log distribution, Newton Raphson iteration

grspar <- c()                            
source('grsm.r', echo=FALSE)                        # Generalised Row Sum Ratings
grs     <- growsums
colnames(grs) <- paste0("grs", elon)
elongrs <- elon
ygrs    <- y                                        # P(x) = ½ + grs / ygrs / par * 2; grs probability

source('recursive-bhz.r', echo=FALSE)               # Recursive Buchholz

grspar <- 0                                         # Least Square Ratings, add LSM to diagonal Laplacian matrix
source('grsm.r', echo=FALSE)
lsq     <- growsums                                 # P(x) = ½ + lsq / par / 2  ; lsq probability
elonlsq <- elon
ylsq    <- y

# https://en.wikipedia.org/wiki/Rule_of_succession  # Laplace: (1 + successes) / (N + 2)
# https://www.colleyrankings.com/matrate.pdf
# Colleys matrix, ε = 3, three outcomes: 0, 1, 2
# or 2 because number of games in swiss is small
# Colley ranking with mean = 0
grspar = c(2, 1 / Par); 
source('grsm.r', echo=FALSE)
L400     <- growsums * 800                          # Convert to [0% - 100%], Elo domain
elonLap  <- elon
yLap     <- y 

colnames(L400)   <- "L400"
colnames(avg400) <- "A400"
colnames(lsq)    <- "lsq"

stopifnot( sd(avg400 - (lsq/Par)*800, na.rm = TRUE) < 1E5) # lsq ~ AVG400

pefbev <- 1                                          # Perron Frobenius (principle) eigen vector
source('fairbets.r', echo=FALSE)
pev <- fb * sum(points)                              # https://www.arndt-bruenner.de/mathe/scripts/engl_eigenwert2.htm

rm(pefbev)                                           # fair bets
source('fairbets.r', echo=FALSE)
fbpts <- fb * sum(points)                            # normalize by all points (ping-pong model (Volij)

frk <- seq_len(npls)                                 # final rank
rank1 <- matrix(rank(-lsq   , ties.method= "min")); colnames(rank1) <- "  Rk"
rank2 <- matrix(rank(-rbhz,   ties.method= "min")); colnames(rank2) <- "  Rk"
rank3 <- matrix(rank(-avg400, ties.method= "min")); colnames(rank3) <- "  Rk"
rank4 <- matrix(rank(-rrtg,   ties.method= "min")); colnames(rank4) <- "  Rk"
rank5 <- matrix(rank(-L400,   ties.method= "min")); colnames(rank5) <- "  Rk"
rank6 <- matrix(rank(-grs,    ties.method= "min")); colnames(rank6) <- "  Rk"
rank7 <- matrix(rank(-pev,    ties.method= "min")); colnames(rank7) <- "  Rk"
rank8 <- matrix(rank(-fbpts,  ties.method= "min")); colnames(rank8) <- "  Rk"

rrtgf <- rrtg
# Roffset <- NA
# if (length(raticol) == 1L) Roffset <- mean(as.numeric(sub("[ABC] ", "", rdtable[,raticol]) ), na.rm =TRUE)
# if (is.finite(Roffset)) rrtgf <- rrtgf + Roffset

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
report5 <- matrix(round(L400,0))  ; colnames(report5) <- "L400"
report6 <- matrix(round(grs   ,3)); colnames(report6) <- paste0("grs",elongrs)
report7 <- matrix(round(pev   ,3)); colnames(report7) <- "Pev"
report8 <- matrix(round(fbpts ,3)); colnames(report8) <- "f-bets"

# calculate score summary.
cbind(rank1, rank2, rank3, rank4, rank5, rank6, rank7, rank8)      -> t1    # drop 3,4
{t1[t1 > 15] <- 0  ; t1} |> as.character.hexmode(keepStr=TRUE)     -> t2    # single hex digit 
{t2[t2=="0"] <- "."; t2} |> apply(1, paste, collapse = "")         -> score
score[lengths(regmatches(score, gregexpr(".", score, fixed=TRUE))) == ncol(t1)] <- "" # remove all blanks
rm(t1, t2)

report  <- cbind(rdtable[,gfrmcols]
                , N, wins, losses, points, s, W, round(pctW, 2)
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
cat(sprintf("L400 : e'= %2d, y = %.2f", elonLap, yLap), crlf)
cat(sprintf("grs  : e'= %2d, y = %.2f (= e' + n.m)", elongrs, ygrs), crlf, crlf)

print(report, quote=FALSE, row.names = FALSE, , na.print = ".", max=200*ncol(report))
}
rm(grspar)
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

