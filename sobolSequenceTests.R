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

# two dimensional
sobolPoints2d <- sobolSequence.points(2,10,200)
plot(sobolPoints2d)

# three dimensional
require(plot3Drgl)
sobolPoints3d <- sobolSequence.points(3,10,200)
plot3d(sobolPoints3d)

