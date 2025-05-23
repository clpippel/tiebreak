# gf.r
# Read game file, make tables to calculate rankings.
# Unlicense (Ø) 2021 CLP, IJmuiden
# Hieke Hylkema (Lin Yutan):
# - Een goede reiziger weet niet waar hij naar toe gaat.
# - Een uitstekende reiziger weet niet waar hij vandaan komt.
#
#   2022-aug-14, SCC, make_cluster.
#   2022-aug-14, handling of nonexisting opponents improved.
#   2023-jan-10, calculation of Par corrected.
#   2023-jan-23, tidy source.
#   2023-jan-31, avoid 2-column matrix when indexing matrix.
#   2023-feb-07, correct laplacian calculation.
#   2023-feb-21, SCCs*Players indexed by level simplified by data.frame.
#                Set working directory using R_USER
#   2024-jul-08  European style csv, force automatic rownames.
#   2024-aug-10  Crosstable, allow for no games.
#   2025-may-06  Vectorized version of t.sparse
#                https://stackoverflow.com/questions/79608748/
#                how-to-remove-and-speed-up-a-for-loop-over-matrix-columns-by-vectorisation
#   2025-may-09  Crosstable by tapply
#                https://stackoverflow.com/questions/68970385/
#
#   CSV, Comma Separated Values (RFC 4180)
#   remove manually "=" from ="0" (excel)
# - read and validate gamefile from Chess-Results, Sevilla / FMJD format into:
#   rdtable, rd = rondedossier (Sevilla)
# - "FMJD" in second line to recognize FMJD style input
# - header-lines before title row
# - title row before first numerical value in first column
# - comment lines out in csv with #
# Format gamefile (rondedossier)
# - Sevilla           Pos Name 1 2 3 4 5 6
#   Results           4b2;1w0;4b1;1w1;
# - Chess results:    Rk. Naam Rtg 1.Rd 2.Rd 3.Rd Pts. TB1 TB2 Rp
#   Results           9w1;1b0;9w½;1b½
# - FMJD Arbiter Pro: Place S.No. Title Name Fed. FMJD R1 R2 R3 R4 R5 R6 R7 R8 R9 Total SSolk
#   Results           9w2;1b0;9w1;1b1;  +2;
#   Tournooibase:     Pl FMJD-id Title Name Cn FMJD-\nrating R1 R2 R3 R4 R5 R6 R7 R8 R9 Pt
#   Results:          2/15w 2/14b 2/4w 1/3w 2/6b 2/8w 1/2b 1/5w 1/9w
#   website:          https://toernooibase.kndb.nl/opvraag/standen.php, edit into
#                     https://toernooibase.kndb.nl/opvraag/matrixnewfmjd.php
#
#   Program flow:
#   - Copy, paste this code into R.
#   - read inputfile into rdtable
#   - delete columns not containing games, output in rdm (ronde dossier matrix)
#   - create opponent and results matrix (n x r)
#   - find weakly and strongly connected components (components in R)
#   - create matrix crosstab for excel / FMJD forum (npls < 100, nrrr = 1)
#   - start report
#
# Tables created:
# rdtable   - csv file
# gfheader  - header lines
# gfrmcols  - columns removed from rdtable
#
# indexed by player x round:
# rdm       - results extracted from rdtable, indexed by player x round
# res_sym   - game results symbols, numeric expressions (.8)
# opponents - opponents
# results   - skew symmetric game results (numeric)
# gfile     - game results (numeric)
#
# crosstab  - indexed by player x player
# SCC       - Strongly connected components, $membership, $csize, $no, $groups
# graph g and quotient graph qg
#
# ls() to list all objects
# edit(table) to inspect
#
# useful R:
# typeof(), mode(), storage.mode()
# str(), structure(), dim(), attributes(), summary()
# class(), methods(class="igraph")
# dput(), length(), object.size(), tracemem(), nchar()
# vector("character", 10), numeric(5), logical(5), list()
# names, dimnames, dim, setNames
# Integer constant, 1L, 2L
# head, tail, nrow(), ncol(), complete.cases()
# sapply(dataframe, class)
# serialize(), readRDS(), saveRDS(), sink("r-output.txt"), sink()
# as_edgelist, get.edgelist
# ls("package:igraph"), lsf.str("package:igraph"), methods("print")
# .libPaths(), sessionInfo()
# https://cran.r-project.org/doc/manuals/r-devel/R-lang.html
#
# edit D:\Program Files\R\R-3.6.1\etc\Rprofile.site to set default directory to \Work
# https://stackoverflow.com/questions/22432344/how-do-you-change-the-default-directory-in-rstudio-or-r
# setwd("d:\\Users\\username\\Documents\\Work")
# clear workspace (ctrl L), attach library and set working directory
# remove.packages('igraph')
# install.packages('igraph')
# lsf.str("package:base")        # Apply lsf.str function
#            # Apply lsf.str function
# installed.packages()           # List ...
# remove.packages(pkgs, lib)     # remove package
#

cat(rm(list = ls()))                                # remove all objects <-----------------+

setwd(paste0(Sys.getenv("R_USER"), "/Work"))        # set working directory
if (!require(igraph)) install.packages("igraph")    # install.packages("igraph") or use menu
igraph_version()                                    # igraph version
sessionInfo()

# ----------------------------------------------------
# A graph interpretation of the least squares ranking method
# https://arxiv.org/abs/1508.06778v1
# Matches matrix Mij: number of comparisons between i, j
# Laplacian matrix L, lij = -mij (i<>k), lii sum of mij
# Gamefile format: mij = 1
#
# From Incomplete Preferences to Ranking via Optimization, 1997, p5
# https://arxiv.org/abs/math/0602552v1
# Score system: Chess, FMJD
# Teams: 20-0, 19-1 .. 0 - 20
# Delftse telling: 12-0, 9-1, 8-2, 7-3 en 6-4
# Skew symmetric score Rik (p18, formula 25) in:
# Rik = (Own score -/- Opp score) / 2, in:
# Rij is unweighted: zero or one game
# Preference fusion when the number of alternatives exceeds two
# https://arxiv.org/abs/math/0602171v3
# ----------------------------------------------------

# Transpose sparse matrix, no for loop.
# Reverse the outcomes between two opponents.
t.sparse <- function(m, opp) {
  matrix( m[ matrix(c(opp, col(opp)), ncol = 2)],
          ncol = ncol(opp)
        )
}

# transpose directed weighted graph, including isolates
# Note that unlike graph_from_edgelist, add_edges needs a vertex sequence (=transposed edgelist).
t.graph <- function(g) {
  tg <-
    add_edges(delete_edges(g, edges = E(g)), matrix(get.edgelist(g, names = FALSE)[, 2:1], nrow = 2, byrow = TRUE)) # revert edges
  edge_attr(tg) <- edge_attr(g)                                      # save edge attributes
  return(tg)
}

# contract vertices, one by equivalence class (=SCC)
# See quotient group / graph, condensation
make_quotient_graph <- function(g, membership, title) {
  g %>%
    contract(membership, vertex.attr.comb = list(
      category = "first",
      color = "first",
      name = function(x) paste(x[seq_len(min(2, length(x)))], collapse = "/")
    )) %>%
    simplify(edge.attr.comb = "sum", remove.loops = FALSE) %>% # combine multiple edges
    set_graph_attr("name", title) %>%                          # add title to graph
    return
}
#-----------------------------------------------------------------------------------

fcsv <-
  ifelse(exists("choose.files")                     # file name game file
         , choose.files(default = "*.csv", multi = FALSE)
         , "stdin")                                 # online ideone.com (no graphics)
if (is.na(fcsv)) q(); print(fcsv)                   # quit when no input
fcon <- file(fcsv, "r")                             # open connection to read inut from

# first five lines of input must contain the maximum number of columns <<<
# no quoting assumed
# edit(rdtable) to inspect
# "ronde dossier" or game file (npls x nrds),
rdtable <- read.csv2(fcon                           # read game file (n x r)
                     , header = FALSE
                     , row.names = NULL             # automatic numbering
                     , colClasses = "character"
                     , strip.white = TRUE
                     , blank.lines.skip = TRUE
                     , comment.char = "#"
                     , stringsAsFactors = FALSE
                     , encoding = "UTF-8")          # Sevilla = UTF-8, W10 = ANSI
close(fcon)
stopifnot(length(names(rdtable)) != 1)              # read.csv did not recognize table (add ;;;; in first csv line)

trow <- min(which(!is.na(suppressWarnings(as.numeric(rdtable[, 1]))))) - 1L # title row: first row with numeric in first column
stopifnot(trow > 0L)
stopifnot(trow <= 5L)                               # header row must be smaller then 5, to guess max number of collumns

# concatenate nonblank columns by row
gfheader <- cbind(apply(rdtable[seq_len(trow - 1), ], 1, function(row) {paste(row[row != ""], collapse = "; ")}))

names(rdtable) <- trimws(rdtable[trow, ], whitespace = "[\\h\\v]") # set column names from input, remove UTF ws, nbsp
rdtable <- rdtable[-as.vector(seq_len(trow)), ]     # remove header rows, if any
rownames(rdtable) <- NULL                           # recalculate rownames

npls <- nrow(rdtable); names(npls) <- "Number of players"; npls # number of players
stopifnot(isTRUE(npls > 0))                         # non trivial

FMJD <- max(regexpr("*FMJD-report", gfheader)) > 0L # FMJD report style

# columns to remove. Keep: R*10*, .*10*  | 10*..., where * is any sequence
# drop columns with empty header
gfrmcols <- grep("^([R.]+.*\\d+.*|\\d+.*)$", names(rdtable), ignore.case = TRUE, invert = TRUE)

#-----------------------------------------------------------------------------------
rdm <- as.matrix(rdtable[, -gfrmcols])              # ronde dossier matrix remove all columns excepts rounds with games

# Q&D, assume latin1 when non UTF-8 symbols in rdm
if (!isTRUE(all.equal(rdm, iconv(rdm, from = "UTF-8", to = "UTF-8")))) {
  cat("Repair encoding", "\n", all.equal(rdm, iconv(rdm, from = "UTF-8", to = "UTF-8")), "\n")
  rdm <- iconv(rdm, from = "latin1", to = "UTF-8")
}

nrds <- ncol(rdm); names(nrds) <- "Number of rounds"; nrds # number of rounds

halfch <- intToUtf8(0x00BD)                         # ½, indifferent to latin1, ANSI encodings in source program
# create opponent and results matrix from game file, rows: number of players, columns: number of rounds
# regex pattern: lookback = digit, split = color, lookahead = result character
if (isTRUE(FMJD)) {                                 # FMJD, 10x10 draughts, format 1/14b
  pt1 <- "/"
  sbs <- regexpr(pt1, rdm, perl = TRUE)             # index of the result character
  opponents <- structure(as.integer(substr(rdm, sbs + 1, nchar(rdm) - 1)), dim = dim(rdm), dimnames = dimnames(rdm)) # force matrix
  res_sym   <-        substr(rdm, 1, sbs - 1)       # symmetric representation
  col_sym   <- substr(rdm, nchar(rdm), nchar(rdm))  # color symbol
} else {                                            # chess format 14b1
  pt1 <- paste("(?<=[0-9])[wbzshax][", halfch, "0-9]", sep = "") # regex pattern to find result and opponent
  sbs <- regexpr(pt1, rdm, perl = TRUE)             # start, base of the result character
  opponents <- structure(as.integer(substr(rdm, 1, sbs - 1)), dim = dim(rdm), dimnames = dimnames(rdm)) # force matrix
  res_sym   <- substr(rdm, sbs + 1, nchar(rdm))     # symmetric representation
  col_sym   <- substr(rdm, sbs, sbs)                # color symbol
}
oppNOK <- !is.na(opponents) & (opponents > npls | opponents < 1) # quick test: negative or to big and not NA
stopifnot(all(apply(col_sym, 2, function(x) {sum(x == "w")}) == apply(col_sym, 2, function(x) {sum(x == "b")})))

if (any(oppNOK)) {
  print(ifelse(!oppNOK, ".", opponents), quote = FALSE)
  stop("Opponent(s) out of range")
}

opponents[res_sym %in% c("+", "-", "=")] <- NA      # clear opponents with forfeit results
res_sym[is.na(opponents)] <- NA                     # clear results with no opponent

rm(fcon, pt1, sbs, trow)                            # tidy up intermediate vars

# https://stackoverflow.com/questions/58693905/how-to-remove-this-for-loop
# https://www.youtube.com/watch?v=7ha78yWRDlE
# if player opponent exists then player.opponent.opponent == player
if (!all(apply(opponents, 2, function(v) all(seq_along(v) == v[v], na.rm = TRUE)))) {
  print("Columns NOK")
  print(apply(opponents, 2, function(v) all(seq_along(v) == v[v], na.rm = TRUE)))
  # v <- opponents                                  # check player == opponent of opponent
  # cbind(v[,1], v[,1][v[,1]], seq_along(v[,1]))    # check column 1
  stop("\nGamefile not symmetric")
}

# compute results, skew symmetric results
gfile <- gsub(halfch, ".5", res_sym)                # replace ½
gfile <- gsub("\\+$", "",   gfile)                  # replace +, draughts plus score
gfile <- gsub("-$",   "",   gfile)                  # replace -, draughts min  score
gfile <- as.numeric(gfile)                          # game file
stopifnot(!all(is.na(gfile)))                       # decimal symbol must be point (as in R)

dim(gfile) <- dim(opponents)
dimnames(gfile) <- dimnames(opponents)

Par     <- max(gfile + t.sparse(gfile, opponents), na.rm = TRUE); # max score range
names(Par) <- "Maximal win score (1=chess, 2=draughts, 3 = football)"; Par
points <- matrix(rowSums(gfile, na.rm = TRUE))      # points, ex bye results
bhlz   <- rowSums(matrix(points[as.vector(opponents)], nrow(points)), na.rm = TRUE) # Buchholz, Weerstand
bhzp   <- rowSums(matrix(bhlz[as.vector(opponents)], nrow(points)), na.rm = TRUE)   # Solkoff plus
nstl   <- rowSums(gfile * points[as.vector(opponents)], na.rm = TRUE) # Sonneborg-Berger, Neustadtl
FB     <- nstl / rowSums(t.sparse(gfile, opponents), na.rm = TRUE)    # Fairbets first iteration

raticol <- which(regexpr("^([Rr]a?ti?n?g$|FMJD|APRO)$", trimws(names(rdtable))) > 0) # Rating column
aor <- rowMeans(matrix(as.numeric(rdtable[, raticol][opponents]), nrow(rdtable)), na.rm = TRUE) #average opponent rtg

results <- (gfile - t.sparse(gfile, opponents)) / 2 # skew symmetric results
stopifnot(all(abs(colSums(results, na.rm = TRUE) - 0) < .Machine$double.eps ^ 0.5)) # zero sum game

wins   <- matrix(rowSums(results >  0, na.rm = TRUE)) # number of wins
draws  <- matrix(rowSums(results == 0, na.rm = TRUE)) # number of draws
losses <- matrix(rowSums(results <  0, na.rm = TRUE)) # number of losses

stopifnot(sum(draws) %% 2 == 0)                     # number of draws must be even

s    <- matrix(rowSums(results, na.rm = TRUE))      # saldo wins and losses, skew symmetric score
nrrr <- max(unlist(apply(opponents, 1, table))); names(nrrr) <- "Round robin rounds, max same opponent"; nrrr
Rmax <- max(results, na.rm = TRUE); names(Rmax) <- "Maximal skew symmetric win score"; Rmax
Rmin <- min(results, na.rm = TRUE)                  # min win score
Parx <- Rmax - Rmin
Parx; Par
warning(Parx != Par)                                # All draws

# data.entry(opponents)                             # inspect for debugging
# test opponents and results, if not OK check encoding input file != ANSI (½)

if (sum(!is.na(opponents)) != sum(!is.na(results))) {
  sum(!is.na(opponents))
  sum(!is.na(results))
  stop("Input error: Opponents <> Results, check invalid result characters (2, ½, -, +)")
}
stopifnot(all.equal(!is.na(results), !is.na(opponents))) # orphaned result or opponent

# --------------------------------------------------#
# Plot results                                      #
# Make edge list of opponents                       #
# --------------------------------------------------#
# Create edge triple list: (player, opponent, result)
# Subsetting drop=FALSE: keep original dimensions, when elist has a single row.
# For example, if the competition consists of one win.
elist <- na.omit(cbind(seq_len(npls), as.vector(opponents), as.vector(gfile)))
elist <- elist[elist[, 3] != 0, , drop = FALSE]     # keep draws, wins
if (length(elist) == 0L) {                          # test nonempty graph
  print(results)
  stop("graph is empty")
}

# --------------------------------------------------#
# Create directed, named, weighted graph.           #
# Vertices are players. Edges are games.            #
# Draw is half-win and half-loss.                   #
# Weights are game scores.                          #
# Find SCC's in results                             #
# --------------------------------------------------#
cat(sprintf("Igraph version: %s\n", packageVersion("igraph"))) # 0.7.1 or above
g <- make_graph(edges = c(), n = npls)              # include one vertex per player
g <- set_vertex_attr(g, name = "name", value = paste(seq_len(vcount(g)), sep = ""))
g <- add_edges(g, as.vector(t(elist[, 1:2])), weight = elist[, 3]) # add results (draw, wins)

# partition graph into strongly connected components -------
SCC    <- make_clusters(g, components(g, mode = "strong")$membership)
SCC$no <- length(SCC)                               # Number of SCC's
largest_SCC <- which.max(table(SCC$membership))     # Index most populated SCC (main)

g <- set_graph_attr(g, name = "main"
                    , value = paste("Players = ", npls, ", Games = ", sum(wins, draws, losses) / 2, ", SCCs = ", SCC$no))
g <- set_graph_attr(g, name = "layout", value = ifelse(SCC$no == 1, layout_in_circle, layout_with_fr))
g <- set_vertex_attr(g, name = "color", value = rainbow(SCC$no, alpha = 0.7)[SCC$membership])

# plot game file
# Add colours and use the mark.group argument
if (ecount(g) < 1E6) {
  dev.new()
  system.time(plot(g, mark.groups = groups(SCC)))
}
qg <- make_quotient_graph(g, SCC$membership, title = "SCCs in games")

# layered layout, sugiyama -------------------------
dev.new()                                           # comment out for rextester.com
wts <- edge_attr(qg)$weight                         # weight = points by edge
wts[wts == 1] <- NA                                 # remove one game, improve readability

qg <- set_graph_attr(qg, name = "main",
                     value = ifelse(vcount(qg) == 1
                                    , "Graph is strongly connected"
                                    , paste(qg$name, "Sugiyama layers:", sep = ", ")))
qg <- set_graph_attr(qg, name = "layout", value = layout_with_sugiyama)
qg <- set_edge_attr(qg, name = "label", value = wts)    # show numner of games
qg <- set_vertex_attr(qg, name = "V(qg)$label.cex", value = .8) # font size
plot(qg)

# ---------------------------------------------------------------------------------------
# compute layers, shortest paths with negative weights (bellman-ford)
qt <- simplify(qg, remove.loops = TRUE)             # remove edges to self
stopifnot(is_dag(qt))                               # no cycles, acyclic

qt <- set_edge_attr(qt, name = "weight", value = - 1)  # longest path = shortest negative
dis_SCC <- (-distances(qt, v = (V(qt)), mode = "in"))  # matrix of VxV distances of shortest paths
layers <- apply(dis_SCC, 1, max)                       # max per row

# SCCs*Players indexed by level
if (SCC$no > 1L)
  print(data.frame(layers, SCC = groups(SCC)))
print(qg)

# Analyse eigenvalues, condition number
if (npls < 64) try({
  laplacian <- -table(rep(seq(npls), nrds), c(opponents),  useNA="no") # frequency table of games, including games against self
  diag(laplacian) <- -rowSums(laplacian, na.rm = TRUE)    # Laplacian matrix, indexed by player
  lev <- Mod(eigen(laplacian)$values)                     # eigen vector Laplacian
  kappa <- max(lev) / min(lev[lev > 1e-5]); names(kappa) <- "Kappa matches matrix";
  print(kappa)
}, silent = FALSE)

# create tournament cross matrix from opponents, res_sym
# Stackoverflow, questions/79612407/
idx      <- data.frame(x = c(row(opponents)), y = c(opponents), z = c(res_sym))
crosstab <- tapply(idx$z, idx[1:2], FUN = paste0, collapse="")
diag(crosstab) <- "x"
crosstab <- cbind(rdtable[, gfrmcols], crosstab)    # Add header columns
print("Cross table")
print(head(crosstab, 20), quote = FALSE, na.print=".") # print first lines of tournament

print(sprintf("# spelers = %d, aantal rondes = %d ", npls,  nrds), quote = FALSE)

edge_density(as_undirected(g, mode = c("collapse")))
gfheader
rm(raticol, wts, elist, oppNOK, Parx, idx)          # tidy up intermediate vars

# --------------------------------------------------#
# source('report.r', echo=FALSE, print=TRUE)        #
