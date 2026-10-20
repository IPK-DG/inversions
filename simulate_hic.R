library(data.table)

#functions to generate simulated Hi-C links for a single chromosome matching real barley data
#chrlen ... chromosome length in bins
#start/end ... start and end coordinates of inversions, start and end need to have the same length
#noise ... noise parameter
make_sim <- function(chrlen=500, start = c(100,400), end = c(120, 410), noise=0){

 # simulate Hi-C link chromosome without inversion
 b=0:chrlen
 cenpos <- chrlen / 2
 setnames(data.table(expand.grid(b, b, stringsAsFactors=F)), c("bin1", "bin2"))[]->x
 x[, dist := abs(bin1 - bin2)]
 x[dist > 0]->x
 x[, cen := cenpos]
 x[, arm1 := ifelse(bin1 <= cenpos, 0, 1)]
 x[, arm2 := ifelse(bin2 <= cenpos, 0, 1)]
 x[arm1 == 0, cdist1 := (cen - bin1)]
 x[arm2 == 0, cdist2 := (cen - bin2)]
 x[arm1 == 1, cdist1 := (bin1 - cen)]
 x[arm2 == 1, cdist2 := (bin2 - cen)]
 x[, cdist := abs(cdist1 - cdist2)]
 x[, nlinks1 := 10**(ceiling(max(log10(0.5+dist))) - log10(0.5+dist))]
 x[, nlinks2 := 10**(ceiling(max(log10(0.5+cdist))) - log10(0.5+1/8*cdist))]
 x[, nlinks := 0.05*(nlinks1 + 0.03*nlinks2)]
 x[, nlinks := nlinks + noise * runif(.N)]

 # generate inverted chromosome
 x[, .(bin1, bin2)] -> y
 y[, ref_bin1 := bin1]
 y[, ref_bin2 := bin2]
 for(i in 1:length(start)){
  s = start[i]
  e = end[i]
  y[bin1 %between% c(s, e), bin1 := as.integer(s + e - bin1)]
  y[bin2 %between% c(s, e), bin2 := as.integer(s + e - bin2)]
 }

 # simulate Hi-C links for inverted chromosome
 y[, dist := abs(bin1 - bin2)]
 y[dist > 0]->y
 y[, cen := cenpos]
 y[, arm1 := ifelse(bin1 <= cenpos, 0, 1)]
 y[, arm2 := ifelse(bin2 <= cenpos, 0, 1)]
 y[arm1 == 0, cdist1 := (cen - bin1)]
 y[arm2 == 0, cdist2 := (cen - bin2)]
 y[arm1 == 1, cdist1 := (bin1 - cen)]
 y[arm2 == 1, cdist2 := (bin2 - cen)]
 y[, cdist := abs(cdist1 - cdist2)]
 y[, nlinks1 := 10**(ceiling(max(log10(0.5+dist))) - log10(0.5+dist))]
 y[, nlinks2 := 10**(ceiling(max(log10(0.5+cdist))) - log10(0.5+1/8*cdist))]
 y[, nlinks := 0.05*(nlinks1 + 0.03*nlinks2)]
 y[, nlinks := nlinks + noise * runif(.N)]

 # create heterozygous inversion
 x[, .(bin1, bin2, nlinks_x=nlinks)][y[, .(bin1=ref_bin1, bin2=ref_bin2, nlinks_y=nlinks)], on=c("bin1", "bin2")]->xy
 xy[, nlinks_xy := 1/2 * (nlinks_x + nlinks_y)]

 # return list with Hi-C matrices for chromosomes without inversions and with inversions occuring in homozygous and heterozygous state
 list(no=x, hom=y, het=xy)
}
