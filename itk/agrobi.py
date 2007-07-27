
import sys, os, itk
from optparse import OptionParser

parser = OptionParser(usage = '''"Usage: agrobi.py image"
  image: the input image''')

parser.add_option("-v", "--verbose", action="store_true", dest="verbose", help="display infos about progress")
parser.add_option("-V", "--visual-validation", action="store_true", dest="visualValidation", help="write 3 images with the segmented zones overlayed on top of the input channels")
parser.add_option("-t", "--threads", type="int", dest="threads", default=0, help="number of threads to use. Defaults to the number of procs")
parser.add_option("-s", "--stimulation", dest="stimulation", default="?", help="the stimulation text in the output")
parser.add_option("-n", "--nuclei", dest="nuclei", default="/dev/null", help="the output file for the nuclei data. Defaults to /dev/null")
parser.add_option("-g", "--genes", dest="genes", default="-", help="the output file for the genes data. Defaults to standard output")
parser.add_option("-S", "--save-segmentation", action="store_true", dest="saveSegmentation", help="write the segmentation results, so they can be reused later")

opts, args = parser.parse_args()

# check the arguments
if len(args) != 1:
	parser.error("incorrect number of arguments")

inputImageName = args[0]

if opts.verbose:
	itk.auto_progress()
	print >> sys.stderr, "Processing", inputImageName

if opts.threads > 0:
	itk.MultiThreader.SetGlobalDefaultNumberOfThreads( opts.threads )


# init the output files, if needed
if opts.nuclei == "-":
	nucleiFile = sys.stdout
else:
	nucleiFile = file( opts.nuclei, "a" )
	if os.path.getsize( opts.nuclei ) == 0:
		nucleiFile.write( '"stimulation" "img" "nucleus" "size" "elongation" "x" "y" "z" "mean" "sigma" "threshold" "otsu"\n' )

if opts.genes == "-":
	genesFile = sys.stdout
else:
	genesFile = file( opts.genes, "a" )
	if os.path.getsize( opts.genes ) == 0:
		genesFile.write( '"stimulation" "img" "nucleus" "gene" "x" "y" "z" "px" "py" "pz" "dist" "ci" "mean" "max" "median" "nMean" "nMinimum" "nMaximum" "nMedian" "nSigma" "nKurtosis" "nSkewness"\n' )



# itk.auto_not_in_place()


##########################
# nuclei
##########################

readerNuclei = itk.lsm(channel=1)
# fill the holes - nucleoles etc - in the nuclei
fillHoles2D = itk.GrayscaleFillholeImageFilter.IUC2IUC2.New(auto_progress=False)
fillHolesNuclei = itk.SliceBySliceImageFilter.IUC3IUC3.New(readerNuclei, Filter=fillHoles2D.GetPointer())
# remove some noise and smooth the image
medianNuclei = itk.MedianImageFilter.IUC3IUC3.New(readerNuclei)
gaussianNuclei = itk.SmoothingRecursiveGaussianImageFilter.IUC3IUC3.New(medianNuclei, Sigma=0.2)
inputNuclei = gaussianNuclei

# to compare with different thresholding methods
otsuNuclei = itk.OtsuThresholdImageCalculator.IUC3.New(gaussianNuclei)

# now we have 2 things: a large shape used to find the spots, and a smaller one used to precisely find the
# border of the nucleus.

# use the rats algorithm to find the nucleus with precision
gradientNuclei = itk.GradientMagnitudeRecursiveGaussianImageFilter.IUC3IUC3.New(inputNuclei, Sigma=0.36)
robustNuclei = itk.RobustAutomaticThresholdImageFilter.IUC3IUC3IUC3.New(inputNuclei, gradientNuclei)
# remove the objects too small to be a nucleus
binarySizeOpeningRobustNuclei = itk.BinaryShapeOpeningImageFilter.IUC3.New(robustNuclei, Attribute="PhysicalSize", Lambda=150)
fillHolesRobustNuclei = itk.SliceBySliceImageFilter.IUC3IUC3.New(binarySizeOpeningRobustNuclei, Filter=fillHoles2D.GetPointer())
# we have the mask of our nuclei
maskRobustNuclei = fillHolesRobustNuclei
# split and labelize the nuclei
maurerRobustNuclei = itk.SignedMaurerDistanceMapImageFilter.IUC3IF3.New(maskRobustNuclei, UseImageSpacing=True)
watershedRobustNuclei = itk.MorphologicalWatershedImageFilter.IF3IUC3.New(maurerRobustNuclei, Level=1.5, MarkWatershedLine=False) #, FullyConnected=True)
maskWatershedRobustNuclei = itk.MaskImageFilter.IUC3IUC3IUC3.New(watershedRobustNuclei, maskRobustNuclei)
# remove (again) the objects too small to be a nucleus
labelSizeOpeningRobustNuclei = itk.LabelShapeOpeningImageFilter.IUC3.New(maskWatershedRobustNuclei, Attribute="PhysicalSize", Lambda=150)
# remove the nucleus on the border - note that they can touch "a little" the border
labelSizeOnBorderOpeningRobustNuclei = itk.LabelShapeOpeningImageFilter.IUC3.New(labelSizeOpeningRobustNuclei, Attribute="SizeOnBorder", Lambda=1500, ReverseOrdering=True)
# # relabel the objects so we are sure to get consecutive labels after the opening
# relabelRobustNuclei = itk.ShapeRelabelImageFilter.IUC3.New(labelSizeOnBorderOpeningRobustNuclei)
# labelRobustNuclei = relabelNRobustuclei
labelRobustNuclei = labelSizeOnBorderOpeningRobustNuclei

# dilate the mask so it can be used to search the spots a little outside the nuclei
thresholdNuclei = itk.BinaryThresholdImageFilter.IF3IUC3.New(maurerRobustNuclei, UpperThreshold=0.3)
maskNuclei = thresholdNuclei
# split and labelize the nuclei
watershedNuclei = itk.MorphologicalWatershedFromMarkersImageFilter.IF3IUC3.New(maurerRobustNuclei, maskWatershedRobustNuclei, MarkWatershedLine=False) #, FullyConnected=True)
maskWatershedNuclei = itk.MaskImageFilter.IUC3IUC3IUC3.New(watershedNuclei, maskNuclei)
#
# the nuclei on the border are already removed in the robust procedure
#
# # remove the nucleus on the border - note that they can touch "a little" the border
# labelSizeOnBorderOpeningNuclei = itk.LabelShapeOpeningImageFilter.IUC3.New(maskWatershedNuclei, Attribute="SizeOnBorder", Lambda=1500, ReverseOrdering=True)
# # relabel the objects so we are sure to get consecutive labels after the opening
# relabelNuclei = itk.ShapeRelabelImageFilter.IUC3.New(labelSizeOnBorderOpeningNuclei)
# labelNuclei = relabelNuclei
labelNuclei = maskWatershedNuclei


# create the label collection to get some data about the nucleus
labelCollectionRobustNuclei = itk.LabelImageToLabelCollectionImageFilter.IUC3LI3.New(labelRobustNuclei, UseBackground=True)
statisticsLabelCollectionRobustNuclei = itk.StatisticsLabelCollectionImageFilter.LI3IUC3.New(labelCollectionRobustNuclei, inputNuclei)

labelCollectionNuclei = itk.LabelImageToLabelCollectionImageFilter.IUC3LI3.New(labelNuclei, UseBackground=True)
shapeLabelCollectionNuclei = itk.ShapeLabelCollectionImageFilter.LI3.New(labelCollectionRobustNuclei)

# select a single nucleus - see the loop at the end
singleMaskNuclei = itk.BinaryThresholdImageFilter.IUC3IUC3.New(labelNuclei, UpperThreshold=1, LowerThreshold=1)
singleMaskRobustNuclei = itk.BinaryThresholdImageFilter.IUC3IUC3.New(labelRobustNuclei, UpperThreshold=1, LowerThreshold=1)
# and compute the distance map. It is used too get the distance of the spot from the nuclear envelop
# -- the mask is inverted to avoid the 0 distance pixels on the border of the object. We prefer to have them outside of the object.
# roiNucleus = itk.RegionOfInterestImageFilter.IUC3IUC3.New(singleMaskRobustNuclei)
# invertedSingleMaskRobustNuclei = itk.InvertIntensityImageFilter.IUC3IUC3.New(roiNucleus)
invertedSingleMaskRobustNuclei = itk.InvertIntensityImageFilter.IUC3IUC3.New(singleMaskRobustNuclei)
maurerSingleNuclei = itk.SignedMaurerDistanceMapImageFilter.IUC3IF3.New(invertedSingleMaskRobustNuclei, UseImageSpacing=True, SquaredDistance=False) #, InsideIsPositive=True)
ciSingleRobustNuclei = itk.ErodedVolumeFractionMapImageFilter.IF3IF3.New(maurerSingleNuclei)
# use an interpolator to get the distance at the exact center of gravity position
maurerInterpolator = itk.LinearInterpolateImageFunction.IF3D.New(maurerSingleNuclei)
# and another one for the CI
ciInterpolator = itk.LinearInterpolateImageFunction.IF3D.New(ciSingleRobustNuclei)


##########################
# wap
##########################

readerWap = itk.lsm(channel=2)
# mask the cytoplasm: there is too much noise
maskNWap = itk.MaskImageFilter.IUC3IUC3IUC3.New(readerWap, maskNuclei)
# roiWap = itk.RegionOfInterestImageFilter.IUC3IUC3.New(maskNWap)
# again, remove some noise
# medianWap = itk.MedianImageFilter.IUC3IUC3.New(roiWap)
medianWap = itk.MedianImageFilter.IUC3IUC3.New(maskNWap)
gaussianWap = itk.SmoothingRecursiveGaussianImageFilter.IUC3IUC3.New(medianWap, Sigma=0.1)
inputWap = gaussianWap
# lets select a single nucleus in the cleaned input image. 2 versions are used:
# - a first one using the robust nuclei procedure, to detect the spot
# - a second one using the enlarged mask to extend the spots outside the nucleus if needed
maskRobustNucleiWap = itk.MaskImageFilter.IUC3IUC3IUC3.New(inputWap, singleMaskRobustNuclei)
maskNucleiWap = itk.MaskImageFilter.IUC3IUC3IUC3.New(inputWap, singleMaskNuclei)
# wap images are quite difficult to segment, because their is lot of noise which look like... wap spots.
# so keep the 4 more visible spots in the image, after have removed the region in the componenents too
# big to be a spot
maxtreeWap = itk.ImageToMaximumTreeFilter.IUC3CTUC3D.New(maskRobustNucleiWap)
sizeMaxtreeWap = itk.PhysicalSizeComponentTreeFilter.CTUC3D.New(maxtreeWap)
filteredSizeMaxtreeWap = itk.AttributeFilteringComponentTreeFilter.CTUC3D.New(sizeMaxtreeWap, Lambda=0.7, ReverseOrdering=True, FilteringType="Direct")
intensityMaxtreeWap = itk.VolumeLevellingComponentTreeFilter.CTUC3D.New(filteredSizeMaxtreeWap)
keepMaxtreeWap = itk.KeepNLobesComponentTreeFilter.CTUC3D.New(intensityMaxtreeWap, NumberOfLobes=4)
# leavesWap = itk.ComponentTreeLeavesToBinaryImageFilter.CTUC3DIUC3.New(keepMaxtreeWap)
leavesWap = itk.ComponentTreeToImageFilter.CTUC3DIUC3.New(keepMaxtreeWap)
reconsWap = itk.ReconstructionByDilationImageFilter.IUC3IUC3.New(leavesWap, maskNucleiWap)
maximaWap = itk.RegionalMaximaImageFilter.IUC3IUC3.New(reconsWap)
maskWap = maximaWap

labelCollectionWap = itk.BinaryImageToLabelCollectionImageFilter.IUC3LI3.New(maskWap)
statsLabelCollectionWap = itk.StatisticsLabelCollectionImageFilter.LI3IUC3.New(labelCollectionWap, inputWap, InPlace=False)
statisticsNucleiLabelCollectionWap = itk.StatisticsLabelCollectionImageFilter.LI3IUC3.New(labelCollectionWap, readerNuclei)


##########################
# casein
##########################

readerCas = itk.lsm(channel=0)
maskNCas = itk.MaskImageFilter.IUC3IUC3IUC3.New(readerCas, maskNuclei)
# again, remove some noise
medianCas = itk.MedianImageFilter.IUC3IUC3.New(maskNCas)
gaussianCas = itk.SmoothingRecursiveGaussianImageFilter.IUC3IUC3.New(medianCas, Sigma=0.1)
# remove the bigger objects in the background with a white tophat by attribute
# Will do nearly nothing in most of the images, but will drop the nucleus image in some of
# them, for a small cost
sizeOpeningCas = itk.PhysicalSizeOpeningImageFilter.IUC3IUC3.New(gaussianCas, Lambda=2.0)
subtractCas = itk.SubtractImageFilter.IUC3IUC3IUC3.New(gaussianCas, sizeOpeningCas)
inputCas = subtractCas
maskNucleiCas = itk.MaskImageFilter.IUC3IUC3IUC3.New(inputCas, singleMaskNuclei)
# select the spot with a simple threshold
thresholdCas = itk.BinaryThresholdImageFilter.IUC3IUC3.New(maskNucleiCas, LowerThreshold=59)
# keep only the ones at least partially in the nucleus
reconsCas = itk.BinaryReconstructionByDilationImageFilter.IUC3.New(singleMaskRobustNuclei, thresholdCas)
binarySizeOpeningCas = itk.BinaryShapeOpeningImageFilter.IUC3.New(reconsCas, Attribute="PhysicalSize", Lambda=0.02)
maskCas = binarySizeOpeningCas

labelCollectionCas = itk.BinaryImageToLabelCollectionImageFilter.IUC3LI3.New(maskCas)
statisticsLabelCollectionCas = itk.StatisticsLabelCollectionImageFilter.LI3IUC3.New(labelCollectionCas, inputCas, InPlace=False)
statisticsNucleiLabelCollectionCas = itk.StatisticsLabelCollectionImageFilter.LI3IUC3.New(labelCollectionCas, readerNuclei)


##########################
# the input file
##########################

def set_file_name( name ):
	readerNuclei.SetFileName( name )
	readerWap.SetFileName( name )
	readerCas.SetFileName( name )


# set_file_name( "wap_cas_20070320_10m2.lsm" )
# set_file_name( "wap_cas_20070320_12m2.lsm" )
# set_file_name( "wap_cas_20070320_3m1.lsm" )
# set_file_name( "wap_cas_20070320_9m2.lsm" )
# set_file_name( "wap_cas_20070320_8m2.lsm" )
# set_file_name( "wap_cas_20070320_6m1.lsm" )

set_file_name( inputImageName )


##########################
# visual validation
##########################

if opts.visualValidation:
	subLabelNuclei = itk.SubtractImageFilter.IUC3IUC3IUC3.New(labelNuclei, labelRobustNuclei)
	overlayNuclei = itk.LabelOverlayImageFilter.IUC3IUC3IRGBUC3.New(readerNuclei, subLabelNuclei, UseBackground=True)
	labels = itk.NaryRelabelImageFilter.IUC3IUC3.New(subLabelNuclei)
	overlays = itk.LabelOverlayImageFilter.IUC3IUC3IRGBUC3.New(readerWap, labels, UseBackground=True)
	
	# create a cross structuring element - why I haven't implemented that in the FlatStructuringElement class, as well as a Ball()
	# method which support the spacing ?
	crossImage = itk.Image.UC3.New(Regions=[11,11,7])
	crossImage.Allocate()
	crossImage.FillBuffer(0)
	for i in range(11):
		crossImage.SetPixel((5,i,3), 255)
		crossImage.SetPixel((i,5,3), 255)
	for i in range(7):
		crossImage.SetPixel((5,5,i), 255)
	crossKernel = itk.FlatStructuringElement._3.FromImageUC(crossImage)
	
	# and create a dilate filter to draw the crosses on the cas and wap images
	crossDilate = itk.BinaryDilateImageFilter.IUC3IUC3SE3.New(Kernel=crossKernel, ForegroundValue=200)
	
	# a function to copy the output of a filter in a new image
	def copyImage( f ):
		f.Update()
		i = f.GetOutput()
		imgDuplicator = itk.ImageDuplicator[i].New(i)
		imgDuplicator.Update()
		return imgDuplicator.GetOutput()
	
	# the list to store the images of the wap and caseins spots
	caseins = []
	waps = []

	# v = itk.show(labels, MaxOpacity=0.05)
	itk.write(overlayNuclei, readerNuclei.GetFileName()+"-nuclei.tif") #, True)




# let start, really
statisticsLabelCollectionRobustNuclei.Update()
shapeLabelCollectionNuclei.Update()
otsuNuclei.Compute()

if opts.saveSegmentation:
	itk.write( labelRobustNuclei, readerNuclei.GetFileName()+"-nuclei-segmentation.nrrd", True)

# to be reused later
spacing = itk.spacing(readerNuclei)

# find the labels used - we are not sure to have all the labels in the range because of the attribute openongs
ls = [l+1 for l in range(*itk.range(labelRobustNuclei)) if statisticsLabelCollectionRobustNuclei.GetOutput().HasLabel(l+1)]

for l in ls :
	if opts.verbose:
		print >> sys.stderr, "  nuclei", l
	
	# set the label
	singleMaskNuclei.SetUpperThreshold( l )
	singleMaskNuclei.SetLowerThreshold( l )
	singleMaskRobustNuclei.SetUpperThreshold( l )
	singleMaskRobustNuclei.SetLowerThreshold( l )
	
	nucleusObject = statisticsLabelCollectionRobustNuclei.GetOutput().GetLabelObject(l)
	nucleusSize = nucleusObject.GetPhysicalSize()
	nucleusElongation = nucleusObject.GetRegionElongation()
	nucleusIdx = [int(round(v)) for v in nucleusObject.GetCentroid()]
	nucleusMean = nucleusObject.GetMean()
	nucleusSigma = nucleusObject.GetSigma()
	
	print >> nucleiFile, '"%s"' % opts.stimulation, '"%s"' % readerNuclei.GetFileName(), l, nucleusElongation, nucleusSize, nucleusIdx[0], nucleusIdx[1], nucleusIdx[2], nucleusMean, nucleusSigma, robustNuclei.GetThreshold(), otsuNuclei.GetThreshold()
	
# 	roi = shapeLabelCollectionNuclei.GetOutput().GetLabelObject(l).GetRegion()
# 	roiNucleus.SetRegionOfInterest( roi )
# 	roiWap.SetRegionOfInterest( roi )
	
	# update the distance map
	ciSingleRobustNuclei.UpdateLargestPossibleRegion()
	
	if opts.visualValidation:
		# put the segmented waps and caseins in a new image
		tempLabelWap = copyImage( maskWap )
		tempLabelCas = copyImage( maskCas )
		
	# get info for wap
	statsLabelCollectionWap.UpdateLargestPossibleRegion()
	statisticsNucleiLabelCollectionWap.UpdateLargestPossibleRegion()
	wapObjects = statsLabelCollectionWap.GetOutput()
	for wl in range(1, wapObjects.GetNumberOfObjects()+1) :
		labelObject = wapObjects.GetLabelObject(wl)
		cog = labelObject.GetCenterOfGravity()
		centerContinuousIdx = [v/s for v, s in zip(cog, spacing)]
		centerIdx = [int(round(v)) for v in centerContinuousIdx]
		dist = maurerInterpolator.EvaluateAtContinuousIndex( centerContinuousIdx )
		ci = ciInterpolator.EvaluateAtContinuousIndex( centerContinuousIdx )
		  
		labelObjectNuclei = statisticsNucleiLabelCollectionWap.GetOutput().GetLabelObject(wl)
		
		print >> genesFile, '"%s"' % opts.stimulation, '"%s"' % readerNuclei.GetFileName(), l, '"wap"', centerIdx[0], centerIdx[1], centerIdx[2], cog[0], cog[1], cog[2], dist, ci, labelObject.GetMean(), labelObject.GetMaximum(), labelObject.GetMedian(), labelObjectNuclei.GetMean(), labelObjectNuclei.GetMinimum(), labelObjectNuclei.GetMaximum(), labelObjectNuclei.GetMedian(), labelObjectNuclei.GetSigma(), labelObjectNuclei.GetKurtosis(), labelObjectNuclei.GetSkewness()
		
		if opts.visualValidation:
			# write a single pixel in the output image to mark the center of the spot
			tempLabelWap.SetPixel( centerIdx, 200 )
	
	# get info for cas
	statisticsLabelCollectionCas.Update()
	statisticsNucleiLabelCollectionCas.UpdateLargestPossibleRegion()
	casObjects = statisticsLabelCollectionCas.GetOutput()
	for wl in range(1, casObjects.GetNumberOfObjects()+1) :
		labelObject = casObjects.GetLabelObject(wl)
		cog = labelObject.GetCenterOfGravity()
		centerContinuousIdx = [v/s for v, s in zip(cog, spacing)]
		centerIdx = [int(round(v)) for v in centerContinuousIdx]
		dist = maurerInterpolator.EvaluateAtContinuousIndex( centerContinuousIdx )
		ci = ciInterpolator.EvaluateAtContinuousIndex( centerContinuousIdx )
		
		labelObjectNuclei = statisticsNucleiLabelCollectionCas.GetOutput().GetLabelObject(wl)
		
		print >> genesFile, '"%s"' % opts.stimulation, '"%s"' % readerNuclei.GetFileName(), l, '"cas"', centerIdx[0], centerIdx[1], centerIdx[2], cog[0], cog[1], cog[2], dist, ci, labelObject.GetMean(), labelObject.GetMaximum(), labelObject.GetMedian(), labelObjectNuclei.GetMean(), labelObjectNuclei.GetMinimum(), labelObjectNuclei.GetMaximum(), labelObjectNuclei.GetMedian(), labelObjectNuclei.GetSigma(), labelObjectNuclei.GetKurtosis(), labelObjectNuclei.GetSkewness()
		
		if opts.visualValidation:
			# write a single pixel in the output image to mark the center of the spot
			tempLabelCas.SetPixel( centerIdx, 200 )
		
	if opts.visualValidation:
		# draw the crosses on the image and copy the result to the image list
		crossDilate.SetInput( tempLabelWap )
		waps.append( copyImage( crossDilate ) )
		crossDilate.SetInput( tempLabelCas )
		caseins.append( copyImage( crossDilate ) )
	

if opts.visualValidation:
	overlays.SetInput( readerWap.GetOutput() )
	for i, wap in enumerate( waps ):
		labels.SetInput( i+1, wap)
	itk.write(overlays, readerWap.GetFileName()+"-wap.tif") #, True)

	overlays.SetInput( readerCas.GetOutput() )
	for i, cas in enumerate( caseins ):
		labels.SetInput( i+1, cas)
	itk.write(overlays, readerWap.GetFileName()+"-cas.tif") #, True)


