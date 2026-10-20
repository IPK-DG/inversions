library(igraph)
library(data.table)

findInversions<-function(links, genotypes, chromosome_lengths,
			 chromosomes=NULL, merge_dist=2e6, maxdist=1e8, multiplier=1e8, 
			 cab=c(1,1.8), coefficient=-2, rsquare=0.4, dcab=1.6){
 if(length(genotypes) > 2){
  stop("More than two genotypes specified.")
 }

 p1 <- genotypes[1]
 p2 <- genotypes[2]

 setnames(links[chr1 == chr2][, chr2 := NULL], "chr1", "chr")[]->zz

 if(is.null(chromosomes)){
  chrs = unique(zz$chr)
 } else {
  zz[chr %in% chromosomes]->zz
 }

 zz[, dist := abs(bin2 - bin1)]

 #calculate directionality biases
 zz[dist <= maxdist & dist > 0]->zz
 zz[, l := log(multiplier) + log(nlinks) - log(sum(nlinks)), by=genotype][l > 0]->zz
 chromosome_lengths[, .(chr, length)][zz, on="chr"]->zz
 zz[dist <= bin1 & dist <= bin2 & (length - bin1) >= dist & (length - bin2) >= dist]->zz
 dcast(zz, chr + bin1 + bin2 ~ genotype, value.var='l', fill=0)->z
 z[ , .(r = sum(sign(bin1 - bin2) * (get(p1) - get(p2)))), key=.(chr, bin=bin1)]->bias

 #pick candidate inversions
 zz[, r := l - mean(l), by=.(dist, genotype)]
 dcast(zz[genotype %in% c(p1, p2)], chr + dist + bin1 + bin2 ~ genotype, value.var="r")->w
 w[!is.na(get(p1)) & !is.na(get(p2))]->w
 w[, rr := get(p1) - get(p2)]
 w[get(p1) < cab[1] & get(p2) > cab[2]]->ww
 w0 <- w

 if(nrow(ww) > 0){
  ww[bin1 < bin2] -> ww
  setorder(ww, chr, bin1, bin2)
  ww[, id := 1:.N]

  #fit linear models for each candidate
  rbindlist(lapply(ww$id, function(i){
   w<-bias[chr == ww[id == i, chr] & bin >= ww[id == i, bin1] & bin <= ww[id == i, bin2]]
   w[, x := bin/1e6]
   w[, lm(r ~ x)] -> l
   data.table(id=i, r2=summary(l)$adj.r.squared, lm=as.numeric(l$coefficients[2]))
  }))[ww, on="id"]->ww
  ww[lm <= coefficient & r2 >= rsquare & abs(rr) >= dcab]->ww
  ww -> ww0

  #merge overlapping inversions
  if(nrow(ww) > 0){
   ww[, .(id1=id, chr, len1=dist, start=bin1, end=bin2)]->a
   ww[, .(id2=id, chr, len2=dist, start=bin1, end=bin2)]->b
   setkey(a, chr, start, end)
   setkey(b, chr, start, end)

   binsize <- merge_dist
   foverlaps(a,b)[id1 != id2] -> ff
   ff[start < i.start | (start == i.start & end < i.end)]->ff
   ff[, o := end - i.start]
   ff[abs(i.start - start) <= binsize | abs(i.end - end) <= binsize | o == pmin(len1, len2)] -> ff

   graph.edgelist(directed=F, as.matrix(ff[, .(paste(id1), paste(id2))])) -> g
   data.table(id=as.integer(V(g)$name), cl=clusters(g)$mem)[ww, on="id"]->ww
   if(any(!is.na(ww$cl))){
    ww[is.na(cl), cl := max(na.omit(ww$cl)) + 1:.N]
   } else {
    ww[, cl := 1:.N]
   }
   ww[order(rr)][!duplicated(cl)][, cl := NULL][]->ww2

   setnames(ww2, c("bin1", "bin2"), c("start", "end"))
   ww2[, .(chr, start, end, size=dist, p1=get(p1), p2=get(p2), coefficient=lm, rsquare=r2, dcab=abs(rr))]->res
   setnames(res, c("p1", "p2"), c(p1, p2))
   setorder(res, chr, start, end)
  } else {
   res <- NULL
  }
 } else {
  res <- NULL
 }

 setnames(bias, "r", "directionality_bias")
 list(putative_inversions=res, directionality_bias=bias)

}
