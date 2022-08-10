# convert crosstable to gamefile 
# select round robin round (rrrounds)
# Example from:
#   https://chess-results.com/tnr437642.aspx?lan=1&art=4
# input: 
#   ;Final Ranking crosstable after 5 Rounds
#   Rk.;;Name;Rtg;FED;1;2;3;4;5;6;Pts.;TB1;TB2 
#   1;GM;Zaragatski Ilja;2482;GER;*;1;1;1;1;1;5;10,00;0
#   2;;Kornitzky Tino;2055;GER;0;*;½;1;1;1;3,5;5,00;0
#   3;;Wüest Andrin;2225;SUI;0;½;*;½;1;1;3;4,00;0
#   4;;Aeschbacher Johann;2053;FRA;0;0;½;*;1;1;2,5;2,50;0
#   5;;Vögtlin André;1719;SUI;0;0;0;0;*;1;1;0,00;0
#   6;;Züllig Flavian;1379;SUI;0;0;0;0;0;*;0;0,00;0
#   # Annotation:
#   # Tie Break1: Sonneborn-Berger-Tie-Break variable
#   # Tie Break2: Direct Encounter (The results of the players in the same point group)
# output:
#     Rk. Var.2               Name  Rtg FED Pts.   TB1 TB2   1   2   3   4   5
#   1   1    GM    Zaragatski Ilja 2482 GER    5 10,00   0 6w1 2w1 3b1 4w1 5b1
#   2   2           Kornitzky Tino 2055 GER  3,5  5,00   0 5w1 1b0 6b1 3w½ 4b1
#   3   3             Wüest Andrin 2225 SUI    3  4,00   0 4w½ 5b1 1w0 2b½ 6b1
#   4   4       Aeschbacher Johann 2053 FRA  2,5  2,50   0 3b½ 6b1 5w1 1b0 2w0
#   5   5            Vögtlin André 1719 SUI    1  0,00   0 2b0 3w0 4b0 6b1 1w0
#   6   6           Züllig Flavian 1379 SUI    0  0,00   0 1b0 4w0 2b0 5w0 3b0

cat(rm(list = ls()))                                # remove all objects
if (!require(igraph)) q()                           # install.packages("igraph") or use menu
options(error=traceback)
cat("\014")                                         # clear console
sessionInfo()

rrrounds <- 1

fcsv <- ifelse(exists("choose.files")
              ,choose.files(default="*.csv", multi = FALSE)
              ,"stdin")                             # online ideone.com (no graphics)
if (is.na(fcsv)) q()                                # quit when no input
fcon <- file(fcsv, "r")

rdtable <- read.csv(fcon                            # read game file (n x r)
                   ,header = FALSE
                   ,sep=";", dec=",", strip.white=TRUE, blank.lines.skip = TRUE
                   ,comment.char = "#"
                   ,stringsAsFactors=FALSE
                   ,encoding = "UTF-8"              # Sevilla = UTF-8, W10 = ANSI
                   )
close(fcon)

trow <- min(which(!is.na(suppressWarnings(as.numeric(rdtable[,1]))))) - 1 # find first title row
stopifnot( trow > 0)

# concatenate nonblank columns by row
gfheader <- cbind(apply(rdtable[seq_len(trow-1),], 1, function(row) {paste(row[row != ""], collapse = "; ")} ) )


names(rdtable) <- rdtable[trow,]                    # set column names from input
rdtable <- rdtable[-c(seq_len(trow)),]              # remove header rows, if any
rownames(rdtable) <- c()                            # remove rownames

npls = nrow(rdtable)                                # number of players
stopifnot(npls > 1)                                 # non trivial                  

# Keep: S9+, X9+, V9+, 9+, + = one or more occurrences
# gfrmcols <- grep("^([SVX.]+.*\\d+.*|\\d+.*)$", names(rdtable), ignore.case=TRUE, invert=TRUE)
# gfrmcols <- grep("^[SXV]?\\d+$", names(rdtable), ignore.case=TRUE, invert=TRUE)


# columns to remove. Keep: R*10*, .*10*  | 10*..., where * is any sequence
gfrmcols <- grep("^([R.]+.*\\d+.*|\\d+.*)$", names(rdtable), ignore.case=TRUE, invert=TRUE)

rdm <- trimws(as.matrix(rdtable[, -gfrmcols]), which="both") # remove  # remove Tiebreak, Pts columns, keep results
rdm[grep('R$', rdm, ignore.case=TRUE)] <- ""        # remove reglementaire uitslagen, xR

res_sym <- substr(rdm,rrrounds, rrrounds)[,1:npls]  # results symbols
nrds <- npls - (1 - npls%%2)                        # of rounds = nr. players minus 1 when even

{ print(sprintf("Players: %d", npls), quote=FALSE);
  print(sprintf("Rounds: %d" , nrds), quote=FALSE);
  print(unlist(dimnames(rdm)[2]), quote=FALSE)
}                                                 # show players, rounds

# Bergertable with N = 11,12 and step N/2
# Ronde                         
# 1     1-12    2-11    3-10    4-9     5-8     6-7
# 2     12-7    8-6     9-5     10-4    11-3    1-2
# 3     2-12    3-1     4-11    5-10    6-9     7-8
# 4     12-8    9-7     10-6    11-5    1-4     2-3
# 5     3-12    4-2     5-1     6-11    7-10    8-9
# 6     12-9    10-8    11-7    1-6     2-5     3-4
# 7     4-12    5-3     6-2     7-1     8-11    9-10
# 8     12-10   11-9    1-8     2-7     3-6     4-5
# 9     5-12    6-4     7-3     8-2     9-1     10-11
# 10    12-11   1-10    2-9     3-8     4-7     5-6
# 11    6-12    7-5     8-4     9-3     10-2    11-1 
# ---------------------------------------------------
rrN <- ifelse(npls%%2, npls + 1 , npls )            # make schedule for even players
opp <- matrix(0, rrN, npls)                         # opponents
col <- matrix("", rrN, npls)                        # colors

btab <- matrix(c(seq_len(rrN/2), seq(from=rrN, to=rrN/2+1)), rrN/2, 2) # berger table round 1
if (rrN > npls) btab[1,2] <- 1                      # oneven, create game against self
rrstep <- rrN / 2                                   # step for next berger schedule

# create opponents, colors with berger table
for (r in seq_len(nrds)) {
  rrself <- btab[,1]
  rropp  <- btab[,2]
  opp[matrix(rrself),r] <- rropp
  opp[matrix(rropp),r]  <- rrself
  opp[which(opp[, r] == seq_len(rrN)), r] <- (if (rrN > npls) NA else rrN) # by, or N
  col[matrix(rrself),r] <- "w"
  col[matrix(rropp) ,r] <- "b"
  if (rrN == npls & r > 1) {
    opp[rrN, r] <- (opp[rrN, r-1] - 1 + rrN/2) %% (rrN-1) + 1 # seperate schedule for N = even
    col[rrN, r] <- ifelse(col[rrN, r-1] == "w", "b", "w")
  }
  btab <- ((btab + rrstep - 1) %% (rrN-1) ) + 1     # next berger round 
 }
if (rrN > npls) opp <- opp[-rrN,]                   # remove N when number of players is oneven

# combine gamefile (N x R) ---------------------------
# format: opp-col-res (14b1)
gfile <- matrix("", npls, nrds)
for (r in seq_len(nrds) ) {
  for (i in seq_len(nrow(opp)) ) {
    val <- opp[i, r]
     if (!is.na(val)) {
        if (res_sym[i, val] != "") gfile[i, r] <- paste0(opp[i,r], col[i,r], res_sym[i, val])
     }
  }
}

# gfile <- cbind(gfile1, gfile2)                    # combine round robin rounds (multiple round robin)
# add headers
gamefile <- cbind(rdtable[,gfrmcols], gfile)

{
writeLines(gfheader)
write.table(gamefile, quote=FALSE, row.names = FALSE, sep=";")
}
# rm(btab, col, gfile, i, opp, r, res_sym, rrN, rropp, rrrounds, rrself, rrstep, val)
# copy / paste from console to Notepad++, check encoding =UTF-8
