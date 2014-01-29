##
## This example shows how to obtain the cdf of the conditional ki as a function
##

require(DNAprofiles)
data(freqsNLngm)

fr <- freqsNLngm # allelic freqs on 15 NGM loci

x <- sample.profiles(N=1,fr) # sample a single profile

# obtain the dist of the sibling index (lr comparing full siblings versus unrelated)
# under H_{true} being full siblings
fs.si.dist <- cond.ki.dist(x,hyp.1="FS",hyp.2="UN",hyp.true="FS",freqs.ki=fr)

# assuming independence between the marker loci, the dist of the SI  is the product of the marginals
prod(sapply(fs.si.dist,function(y) length(y$x))) # quite some events!

# obtain the dist of the product for two subsets of the loci
fs.si.dist.duo <- dists.product.duo(fs.si.dist)
#str(fs.si.dist.duo) # contains the cumulative dist on one subset of the loci; dist on rest of the loci
# obtain the cdf of the sibling index for true sibs as a function
fs.si.cdf <- dist.duo.cdf(fs.si.dist.duo)

## we could also obtain the cdf directly
## fs.si.cdf <- cond.ki.cdf(x,hyp.1="FS",hyp.2="UN",hyp.true="FS",freqs=fr)

fs.si.cdf(1) # pr. that SI<=1 for true sibs
fs.si.cdf(1,exc.prob=TRUE) # pr. that SI>1 for true sibs

# plot the cdf for true sibs
x0 <- seq(from=-10,to=20,length=50)
plot(x0,fs.si.cdf(10^x0),type="l")

# add the cdf for unrelated profiles
unr.si.cdf <- cond.ki.cdf(x,hyp.1="FS",hyp.true="UN",freqs=fr)
lines(x0,unr.si.cdf(10^x0),lty=2)