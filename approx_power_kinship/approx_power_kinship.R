##
## This example shows how to compute the power to discriminate between sibs and halfsib
## using autosomal loci
##

require(DNAprofiles)
data(freqsNLngm)
fr <- freqsNLngm

fs.vs.hs.fs <- ki.dist(hyp.1="FS",hyp.2="HS",hyp.true="FS",freqs.ki=fr)
fs.vs.hs.hs <- ki.dist(hyp.1="FS",hyp.2="HS",hyp.true="HS",freqs.ki=fr)
prod(sapply(fs.vs.hs.fs,function(y) length(y$x))) # we have no hope of fitting this in memory

# `round' to lower values
fs.cdf.l <- dist.pair.cdf(dists.product.pair(fs.vs.hs.fs,appr=TRUE,appr.method=1))
# `round' to higher values
fs.cdf.u <- dist.pair.cdf(dists.product.pair(fs.vs.hs.fs,appr=TRUE,appr.method=2))

hs.cdf.l <- dist.pair.cdf(dists.product.pair(fs.vs.hs.hs,appr=TRUE,appr.method=1))
hs.cdf.u <- dist.pair.cdf(dists.product.pair(fs.vs.hs.hs,appr=TRUE,appr.method=2))

t0 <- 1e3
# prob. for true full sib to have a `large' lr
1-fs.cdf.l(t0)
1-fs.cdf.u(t0)

# prob. for half sib to have a `large' lr
1-hs.cdf.l(t0)
1-hs.cdf.u(t0)