
import sys, os, itk
from optparse import OptionParser

parser = OptionParser(usage = '''"Usage: agrobi.py image"
  image: the input image''')

parser.add_option("-v", "--verbose", action="store_true", dest="verbose", help="display infos about progress")
parser.add_option("-t", "--threads", type="int", dest="threads", default=0, help="number of threads to use. Defaults to the number of procs")

opts, args = parser.parse_args()

# check the arguments
if len(args) != 1:
	parser.error("incorrect number of arguments")
inputFileName = args[0]
inputFile = open( inputFileName )


if opts.verbose:
	itk.auto_progress()

if opts.threads > 0:
	itk.MultiThreader.SetGlobalDefaultNumberOfThreads( opts.threads )




##########################
# nuclei
##########################

labelRobustNuclei = itk.ImageFileReader.IUC3.New()

# select a single nucleus - see the loop at the end
singleMaskRobustNuclei = itk.BinaryThresholdImageFilter.IUC3IUC3.New(labelRobustNuclei, UpperThreshold=1, LowerThreshold=1)
# and compute the distance map. It is used too get the distance of the spot from the nuclear envelop
# -- the mask is inverted to avoid the 0 distance pixels on the border of the object. We prefer to have them outside of the object.
# roiNucleus = itk.RegionOfInterestImageFilter.IUC3IUC3.New(singleMaskRobustNuclei)
# invertedSingleMaskRobustNuclei = itk.InvertIntensityImageFilter.IUC3IUC3.New(roiNucleus)
invertedSingleMaskRobustNuclei = itk.InvertIntensityImageFilter.IUC3IUC3.New(singleMaskRobustNuclei)
maurerSingleNuclei = itk.SignedMaurerDistanceMapImageFilter.IUC3IF3.New(invertedSingleMaskRobustNuclei, UseImageSpacing=True, SquaredDistance=False) #, InsideIsPositive=True)
ciSingleRobustNuclei = itk.CentralIndexMapImageFilter.IF3IF3.New(maurerSingleNuclei)


# print the header
print"stimulation", "img", "nucleus", "x.nucleus", "y.nucleus", "z.nucleus", "threshold", "gene", "x.gene", "y.gene", "z.gene", "px", "py", "pz", "dist", "ci"


# skip the first (header) line
inputFile.readline()

for line in inputFile :
	# read the parameters
	paramsStr = line.split()
	
	# convert numbers to numbers
	params = []
	for p in paramsStr :
		if p.isdigit():
			params.append( int(p) )
		else :
			params.append( p )
	
	stimulation, img, nucleus, xNucleus, yNucleus, zNucleus, threshold, gene, xGene, yGene, zGene = params

	labelRobustNuclei.SetFileName( img[1:-1] + "-nuclei-segmentation.nrrd" ) # [1:-1] to remove the "" around the file name
	labelRobustNuclei.UpdateLargestPossibleRegion()
	
	# get the label of the nucleus
	l = labelRobustNuclei.GetOutput().GetPixel( [xNucleus, yNucleus, zNucleus] )
	
	# and use it to select only that nucleus in the image
	singleMaskRobustNuclei.SetUpperThreshold( l )
	singleMaskRobustNuclei.SetLowerThreshold( l )
	
	# then we can update the ci calculator
	ciSingleRobustNuclei.UpdateLargestPossibleRegion()

	# and read the CI and get the distance from the map
	ci = ciSingleRobustNuclei.GetOutput().GetPixel( [xGene, yGene, zGene] )
	dist = maurerSingleNuclei.GetOutput().GetPixel( [xGene, yGene, zGene] )
	
	# finally compute the physical position 
	pp = itk.index_to_physical_point( labelRobustNuclei, [xGene, yGene, zGene] )
	
	# and print the results
	print stimulation, img, nucleus, xNucleus, yNucleus, zNucleus, threshold, gene, xGene, yGene, zGene, pp[0], pp[1], pp[2], dist, ci
	
