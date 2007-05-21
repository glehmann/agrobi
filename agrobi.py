
import sys, os, itk
itk.auto_progress()
itk.auto_not_in_place()

# to be used with hyper threading procs only!
itk.MultiThreader.SetGlobalDefaultNumberOfThreads(1)

# check the arguments
if len(sys.argv) != 5:
	print >> sys.stderr, "Usage:", sys.argv[0], "stimulation nuclei.txt genes.txt image"
	print >> sys.stderr, '  stimulation: usually "yes" or "no"'
	print >> sys.stderr, '  nuclei.txt: the file where the nuclei data will be put'
	print >> sys.stderr, '  genes.txt: the file where the genes data will be put'
	print >> sys.stderr, '  image: the input image'
	sys.exit(1)

stimulation = sys.argv[1]
inputImageName = sys.argv[4]


print "Processing", sys.argv[4]

# init the output files, if needed
nucleiFileName = sys.argv[2]
nucleiFile = file( nucleiFileName,"a" )
if os.path.getsize( nucleiFileName ) == 0:
  nucleiFile.write( '"stimulation" "img" "nucleus" "size" "elongation" "x" "y" "z" "mean" "sigma" "threshold"\n' )
  
genesFileName = sys.argv[3]
genesFile = file( genesFileName, "a" )
if os.path.getsize( genesFileName ) == 0:
  genesFile.write( '"stimulation" "img" "nucleus" "gene" "x" "y" "z" "px" "py" "pz" "dist" "ci"\n' )
  

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

# now we have 2 things: a large shape used to find the spots, and a smaller one used to precisely find the
# border of the nucleus.

# use the rats algorithm to find the nucleus with precision
gradientNuclei = itk.GradientMagnitudeRecursiveGaussianImageFilter.IUC3IUC3.New(inputNuclei, Sigma=0.36)
robustNuclei = itk.RobustAutomaticThresholdImageFilter.IUC3IUC3IUC3.New(inputNuclei, gradientNuclei)
# remove the objects too small to be a nucleus
binarySizeOpeningRobustNuclei = itk.BinaryShapeOpeningImageFilter.IUC3.New(robustNuclei, Lambda=100000)
fillHolesRobustNuclei = itk.SliceBySliceImageFilter.IUC3IUC3.New(binarySizeOpeningRobustNuclei, Filter=fillHoles2D.GetPointer())
# we have the mask of our nuclei
maskRobustNuclei = fillHolesRobustNuclei
# split and labelize the nuclei
maurerRobustNuclei = itk.SignedMaurerDistanceMapImageFilter.IUC3IF3.New(maskRobustNuclei, UseImageSpacing=True)
watershedRobustNuclei = itk.MorphologicalWatershedImageFilter.IF3IUC3.New(maurerRobustNuclei, Level=1.5, MarkWatershedLine=False) #, FullyConnected=True)
maskWatershedRobustNuclei = itk.MaskImageFilter.IUC3IUC3IUC3.New(watershedRobustNuclei, maskRobustNuclei)
# remove the nucleus on the border - note that they can touch "a little" the border
labelSizeOnBorderOpeningRobustNuclei = itk.LabelShapeOpeningImageFilter.IUC3.New(maskWatershedRobustNuclei, Attribute="SizeOnBorder", Lambda=1500, ReverseOrdering=True)
# # relabel the objects so we are sure to get consecutive labels after the opening
# relabelRobustNuclei = itk.ShapeRelabelImageFilter.IUC3.New(labelSizeOnBorderOpeningRobustNuclei)
# labelRobustNuclei = relabelNRobustuclei
labelRobustNuclei = labelSizeOnBorderOpeningRobustNuclei

# dilate the mask so it can be used to search the spots a little outside the nuclei
thresholdNuclei = itk.BinaryThresholdImageFilter.IF3IUC3.New(maurerRobustNuclei, UpperThreshold=0.3)
maskNuclei = thresholdNuclei
# split and labelize the nuclei
maurerNuclei = itk.SignedMaurerDistanceMapImageFilter.IUC3IF3.New(maskNuclei, UseImageSpacing=True)
watershedNuclei = itk.MorphologicalWatershedFromMarkersImageFilter.IF3IUC3.New(maurerNuclei, maskWatershedRobustNuclei, MarkWatershedLine=False) #, FullyConnected=True)
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
labelCollectionNuclei = itk.LabelImageToLabelCollectionImageFilter.IUC3LI3.New(labelRobustNuclei, UseBackground=True)
statisticsLabelCollectionNuclei = itk.StatisticsLabelCollectionImageFilter.LI3IUC3.New(labelCollectionNuclei, inputNuclei)

subLabelNuclei = itk.SubtractImageFilter.IUC3IUC3IUC3.New(labelNuclei, labelRobustNuclei)
overlayNuclei = itk.LabelOverlayImageFilter.IUC3IUC3IRGBUC3.New(readerNuclei, subLabelNuclei, UseBackground=True)

# select a single nucleus - see the loop at the end
singleMaskNuclei = itk.BinaryThresholdImageFilter.IUC3IUC3.New(labelNuclei, UpperThreshold=1, LowerThreshold=1)
singleMaskRobustNuclei = itk.BinaryThresholdImageFilter.IUC3IUC3.New(labelRobustNuclei, UpperThreshold=1, LowerThreshold=1)
# and compute the distance map. It is used too get the distance of the spot from the nuclear envelop
# -- the mask is inverted to avoid the 0 distance pixels on the border of the object. We prefer to have them outside of the object.
invertedSingleMaskRobustNuclei = itk.InvertIntensityImageFilter.IUC3IUC3.New(singleMaskRobustNuclei)
maurerSingleNuclei = itk.SignedMaurerDistanceMapImageFilter.IUC3IF3.New(invertedSingleMaskRobustNuclei, UseImageSpacing=True, SquaredDistance=False) #, InsideIsPositive=True)
ciSingleRobustNuclei = itk.CentralIndexMapImageFilter.IF3IF3.New(maurerSingleNuclei)
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
# again, remove some noise
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
connectedWap = itk.ConnectedComponentImageFilter.IUC3IUC3.New(leavesWap, FullyConnected=True)
labelWap = connectedWap

labelWapNuclei = itk.NaryRelabelImageFilter.IUC3IUC3.New(subLabelNuclei)
overlayWap = itk.LabelOverlayImageFilter.IUC3IUC3IRGBUC3.New(readerWap, labelWapNuclei, UseBackground=True)

labelCollectionWap = itk.LabelImageToLabelCollectionImageFilter.IUC3LI3.New(labelWap, UseBackground=True)
statsLabelCollectionWap = itk.StatisticsLabelCollectionImageFilter.LI3IUC3.New(labelCollectionWap, inputWap)

##########################
# casein
##########################

readerCas = itk.lsm(channel=0)
maskNCas = itk.MaskImageFilter.IUC3IUC3IUC3.New(readerCas, maskNuclei)
# again, remove some noise
medianCas = itk.MedianImageFilter.IUC3IUC3.New(maskNCas)
gaussianCas = itk.SmoothingRecursiveGaussianImageFilter.IUC3IUC3.New(medianCas, Sigma=0.1)
inputCas = gaussianCas
maskNucleiCas = itk.MaskImageFilter.IUC3IUC3IUC3.New(inputCas, singleMaskNuclei)
# select the spot with a simple threshold
thresholdCas = itk.BinaryThresholdImageFilter.IUC3IUC3.New(maskNucleiCas, LowerThreshold=59)
# keep only the ones at least partially in the nucleus
reconsCas = itk.BinaryReconstructionByDilationImageFilter.IUC3.New(singleMaskRobustNuclei, thresholdCas)
binarySizeOpeningCas = itk.BinaryShapeOpeningImageFilter.IUC3.New(reconsCas, Attribute="PhysicalSize", Lambda=0.02)
maskCas = binarySizeOpeningCas
connectedCas = itk.ConnectedComponentImageFilter.IUC3IUC3.New(maskCas, FullyConnected=True)
labelCas = connectedCas

labelCasNuclei = itk.NaryRelabelImageFilter.IUC3IUC3.New(subLabelNuclei)
overlayCas = itk.LabelOverlayImageFilter.IUC3IUC3IRGBUC3.New(readerCas, labelCasNuclei, UseBackground=True)

labelCollectionCas = itk.LabelImageToLabelCollectionImageFilter.IUC3LI3.New(labelCas, UseBackground=True)
statisticsLabelCollectionCas = itk.StatisticsLabelCollectionImageFilter.LI3IUC3.New(labelCollectionCas, inputCas)


labels = itk.NaryBinaryToLabelImageFilter.IUC3IUC3.New(maskRobustNuclei, maskWap, maskCas)


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

# v = itk.show(labels, MaxOpacity=0.05)
itk.write(overlayNuclei, readerNuclei.GetFileName()+"-nuclei.tif") #, True)

statisticsLabelCollectionNuclei.Update()


caseins = []
waps = []

# to be reused later
spacing = itk.spacing(readerNuclei)

# find the labels used - we are not sure to have all the labels in the range because of the attribute openongs
ls = [l+1 for l in range(*itk.range(labelRobustNuclei)) if statisticsLabelCollectionNuclei.GetOutput().HasLabel(l+1)]

for l in ls :
	# set the label
	singleMaskNuclei.SetUpperThreshold( l )
	singleMaskNuclei.SetLowerThreshold( l )
	singleMaskRobustNuclei.SetUpperThreshold( l )
	singleMaskRobustNuclei.SetLowerThreshold( l )
	
	# update the distance map
	ciSingleRobustNuclei.Update()
	
	nucleusObject = statisticsLabelCollectionNuclei.GetOutput().GetLabelObject(l)
	nucleusSize = nucleusObject.GetPhysicalSize()
	nucleusElongation = nucleusObject.GetRegionElongation()
	nucleusIdx = [int(round(v)) for v in nucleusObject.GetCentroid()]
	nucleusMean = nucleusObject.GetMean()
	nucleusSigma = nucleusObject.GetSigma()
	
	print >> nucleiFile, '"%s"' % stimulation, '"%s"' % readerNuclei.GetFileName(), l, nucleusElongation, nucleusSize, nucleusIdx[0], nucleusIdx[1], nucleusIdx[2], nucleusMean, nucleusSigma, robustNuclei.GetThreshold()
	
	
	# put the segmented wap in a new image
	tempLabelWap = copyImage( maskWap )
		
	# get info for wap
	statsLabelCollectionWap.Update()
	wapObjects = statsLabelCollectionWap.GetOutput()
	for wl in range(1, wapObjects.GetNumberOfObjects()+1) :
		cog = wapObjects.GetLabelObject(wl).GetCenterOfGravity()
		centerContinuousIdx = [v/s for v, s in zip(cog, spacing)]
		centerIdx = [int(round(v)) for v in centerContinuousIdx]
		dist = maurerInterpolator.EvaluateAtContinuousIndex( centerContinuousIdx )
		ci = ciInterpolator.EvaluateAtContinuousIndex( centerContinuousIdx )
		  
		print >> genesFile, '"%s"' % stimulation, '"%s"' % readerNuclei.GetFileName(), l, '"wap"', centerIdx[0], centerIdx[1], centerIdx[2], cog[0], cog[1], cog[2], dist, ci
		
		# write a single pixel in the output image to mark the center of the spot
		tempLabelWap.SetPixel( centerIdx, 200 )
	
	# draw the crosses on the image
	crossDilate.SetInput( tempLabelWap )
	# and copy the result to the image list
	waps.append( copyImage( crossDilate ) )
	
	
	# put the segmented cas in a new image
	tempLabelCas = copyImage( maskCas )
	
	# get info for cas
	statisticsLabelCollectionCas.Update()
	casObjects = statisticsLabelCollectionCas.GetOutput()
	for wl in range(1, casObjects.GetNumberOfObjects()+1) :
		cog = casObjects.GetLabelObject(wl).GetCenterOfGravity()
		centerContinuousIdx = [v/s for v, s in zip(cog, spacing)]
		centerIdx = [int(round(v)) for v in centerContinuousIdx]
		dist = maurerInterpolator.EvaluateAtContinuousIndex( centerContinuousIdx )
		ci = ciInterpolator.EvaluateAtContinuousIndex( centerContinuousIdx )
		
		print >> genesFile, '"%s"' % stimulation, '"%s"' % readerNuclei.GetFileName(), l, '"cas"', centerIdx[0], centerIdx[1], centerIdx[2], cog[0], cog[1], cog[2], dist, ci
		
		# write a single pixel in the output image to mark the center of the spot
		tempLabelCas.SetPixel( centerIdx, 200 )
		
	# draw the crosses on the image
	crossDilate.SetInput( tempLabelCas )
	# and copy the result to the image list
	caseins.append( copyImage( crossDilate ) )
	

for i, (cas, wap) in enumerate( zip( caseins, waps ) ):
	labelCasNuclei.SetInput( i+1, cas)
	labelWapNuclei.SetInput( i+1, wap)
	
itk.write(overlayWap, readerWap.GetFileName()+"-wap.tif") #, True)
itk.write(overlayCas, readerCas.GetFileName()+"-cas.tif") #, True)

