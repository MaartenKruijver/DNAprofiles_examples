require(DNAprofiles)
data(freqsNLngm)
fr <- freqsNLngm

# distribution of siblind index for pairs of true sibs at all loci
si.fs.dist <- ki.dist(hyp.1="FS",hyp.2="UN",hyp.true="UN",freqs.ki=fr)

prod(sapply(si.fs.dist,function(y) length(y$x))) # we have no hope of fitting this in memory

# but we can get exact bounds
si.fs.cdf.l <- dist.duo.cdf(dists.product.duo.appr(si.fs.dist,appr.method=1)) # `round' to lower values
si.fs.cdf.u <- dist.duo.cdf(dists.product.duo.appr(si.fs.dist,appr.method=2)) # `round' to higher values

si.fs.cdf.l(1e6,exc.prob=TRUE)
si.fs.cdf.u(1e6,exc.prob=TRUE)