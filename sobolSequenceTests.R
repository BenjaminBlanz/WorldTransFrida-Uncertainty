# test sobolSequence

require(SobolSequence)
# Usage
# sobolSequence.points(dimR, dimF2 = 10, count, digitalShift = FALSE)
# Arguments
# dimR	        dimention.
# dimF2         F2-dimention of each element.                       No idea what this is
# count	        number of points.                                   
# digitalShift	use digital shift or not.                           No idea what this is
# 
# Value
# matrix of points where every row contains dimR dimensional point.
# 

# two dimensional ####
sobolPoints2d <- sobolSequence.points(2,30,10)
plot(sobolPoints2d,xlim=c(0,1),ylim=c(0,1))
abline(v=0.5,h=0.5,col='red')

# three dimensional ####
require(plot3Drgl)
sobolPoints3d <- sobolSequence.points(3,30,20)
plot3d(sobolPoints3d,xlim=c(0,1),ylim=c(0,1),zlim=c(0,1))
surface3d(c(0,1),c(0,1),matrix(rep(0.5,4),nrow=2),col='white',shininess=1,alpha=0.9)
surface3d(c(0,1),matrix(rep(0.5,4),nrow=2),c(0,1),col='white',shininess=1,alpha=0.9)
surface3d(matrix(rep(0.5,4),nrow=2),c(0,1),c(0,1),col='white',shininess=1,alpha=0.9)

# n dimensional ####
N <- 100
N.sample <- 1e6
sobolPointsNd <- sobolSequence.points(N,10,N.sample)
numInTopNdrant <- sum(colSums(t(sobolPointsNd)>rep(0.5,N))==N)
numInBottomNdrant <- sum(colSums(t(sobolPointsNd)<rep(0.5,N))==N)

cat(sprintf('Expected Value of samples per Ndrant %f\nActual samples in top and bottom Ndrant %i and %i\n',
						N.sample/2^N,
						numInTopNdrant,numInBottomNdrant))
