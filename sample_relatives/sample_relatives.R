##
## This example shows how to sample a profile and relatives (sibs, offspring)
##

require(DNAprofiles)
data(freqsNLngm)

fr <- freqsNLngm # allelic freqs on 15 NGM loci

x <- sample.profiles(N=1,fr) # sample a single profile
fs <- sample.relatives(x,N=1e5,type="FS",freqs=fr) # sample sibs of this profile
## alternatively specify IBD probs: fs <- sample.relatives(x,N=1e5,type=c(1/4,1/2,1/4),freqs=fr)

# compute sibling index (likelihood ratio of full sibs vs. unrelated) with all the true sibs
fs.si <- ki.db(x=x,db=fs,type="FS",freqs=fr) 

# plot an estimate of the conditional ki-distribution of full sibs of x
plot(density(log10(fs.si))) 