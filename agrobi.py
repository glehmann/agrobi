
import sys, os, itk
itk.auto_progress()
itk.auto_not_in_place()

print "Processing", sys.argv[4]

# init the output files, if needed
nucleiFileName = sys.argv[2]
nucleiFile = file( nucleiFileName,"a" )
if os.path.getsize( nucleiFileName ) == 0:
  nucleiFile.write( '"stimulation" "img" "nucleus" "elongation" "size" "x" "y" "z" "mean" "sigma"\n' )
  
genesFileName = sys.argv[3]
genesFile = file( genesFileName, "a" )
if os.path.getsize( genesFileName ) == 0:
  genesFile.write( '"stimulation" "img" "nucleus" "gene" "x" "y" "z" "dist" "ci"\n' )
  


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

# an automatic threshold to find the spots even when they are a little outside the nucleus
kappaNuclei = itk.KappaSigmaThresholdImageFilter.IUC3IUC3.New(inputNuclei)
# remove the objects too small to be a nucleus
binarySizeOpeningNuclei = itk.BinaryShapeOpeningImageFilter.IUC3.New(kappaNuclei, Lambda=100000)
# we have the mask of our nuclei
maskNuclei = binarySizeOpeningNuclei
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

overlayNuclei = itk.LabelOverlayImageFilter.IUC3IUC3IRGBUC3.New(readerNuclei, labelRobustNuclei, UseBackground=True)

# select a single nucleus - see the loop at the end
singleMaskNuclei = itk.BinaryThresholdImageFilter.IUC3IUC3.New(labelNuclei, UpperThreshold=1, LowerThreshold=1)
singleMaskRobustNuclei = itk.BinaryThresholdImageFilter.IUC3IUC3.New(labelRobustNuclei, UpperThreshold=1, LowerThreshold=1)
# and compute the distance map. It is used too get the distance of the spot from the nuclear envelop
maurerSingleNuclei = itk.SignedMaurerDistanceMapImageFilter.IUC3IF3.New(singleMaskRobustNuclei, UseImageSpacing=True, SquaredDistance=False, InsideIsPositive=True)
# the thresholded distance map, to compute the CI
thresholdMaurerSingleNuclei = itk.BinaryThresholdImageFilter.IF3IUC3.New(maurerSingleNuclei)
labelCollectionSingleNuclei = itk.LabelImageToLabelCollectionImageFilter.IUC3LI3.New(thresholdMaurerSingleNuclei, UseBackground=True)
shapeLabelCollectionSingleNuclei = itk.ShapeLabelCollectionImageFilter.LI3.New(labelCollectionSingleNuclei)


##########################
# wap
##########################

readerWap = itk.lsm(channel=2)
# mask the cytoplasm: there is too much noise
maskNWap = itk.MaskImageFilter.IUC3IUC3IUC3.New(readerWap, singleMaskNuclei)
# again, remove some noise
medianWap = itk.MedianImageFilter.IUC3IUC3.New(maskNWap)
gaussianWap = itk.SmoothingRecursiveGaussianImageFilter.IUC3IUC3.New(medianWap, Sigma=0.1)
inputWap = gaussianWap
# wap images are quite difficult to segment, because their is lot of noise which look like... wap spots.
# so keep the 4 more visible spots in the image, after have removed the region in the componenents too
# big to be a spot
maxtreeWap = itk.ImageToMaximumTreeFilter.IUC3CTUC3D.New(inputWap)
sizeMaxtreeWap = itk.PhysicalSizeComponentTreeFilter.CTUC3D.New(maxtreeWap)
filteredSizeMaxtreeWap = itk.AttributeFilteringComponentTreeFilter.CTUC3D.New(sizeMaxtreeWap, Lambda=0.8, ReverseOrdering=True, FilteringType="Subtract")
intensityMaxtreeWap = itk.LocalIntensityComponentTreeFilter.CTUC3D.New(filteredSizeMaxtreeWap)
keepMaxtreeWap = itk.KeepNLobesComponentTreeFilter.CTUC3D.New(intensityMaxtreeWap, NumberOfLobes=4)
leavesWap = itk.ComponentTreeLeavesToBinaryImageFilter.CTUC3DIUC3.New(keepMaxtreeWap)
maskWap = leavesWap
connectedWap = itk.ConnectedComponentImageFilter.IUC3IUC3.New(leavesWap, FullyConnected=True)
labelWap = connectedWap

labelWapNuclei = itk.NaryRelabelImageFilter.IUC3IUC3.New(labelRobustNuclei)
overlayWap = itk.LabelOverlayImageFilter.IUC3IUC3IRGBUC3.New(readerWap, labelWapNuclei, UseBackground=True)

labelCollectionWap = itk.LabelImageToLabelCollectionImageFilter.IUC3LI3.New(labelWap, UseBackground=True)
statsLabelCollectionWap = itk.StatisticsLabelCollectionImageFilter.LI3IUC3.New(labelCollectionWap, inputWap)

##########################
# casein
##########################

readerCas = itk.lsm(channel=0)
maskNCas = itk.MaskImageFilter.IUC3IUC3IUC3.New(readerCas, singleMaskNuclei)
# again, remove some noise
medianCas = itk.MedianImageFilter.IUC3IUC3.New(maskNCas)
gaussianCas = itk.SmoothingRecursiveGaussianImageFilter.IUC3IUC3.New(medianCas, Sigma=0.1)
inputCas = gaussianCas
# keep the 4 more visible spots
thresholdCas = itk.BinaryThresholdImageFilter.IUC3IUC3.New(inputCas, LowerThreshold=59)
binarySizeOpeningCas = itk.BinaryShapeOpeningImageFilter.IUC3.New(thresholdCas, Attribute="PhysicalSize", Lambda=0.02)
maskCas = binarySizeOpeningCas
connectedCas = itk.ConnectedComponentImageFilter.IUC3IUC3.New(maskCas, FullyConnected=True)
labelCas = connectedCas

labelCasNuclei = itk.NaryRelabelImageFilter.IUC3IUC3.New(labelRobustNuclei)
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

set_file_name( sys.argv[4] )

# v = itk.show(labels, MaxOpacity=0.05)
itk.write(overlayNuclei, readerNuclei.GetFileName()+"-nuclei.tif") #, True)

statisticsLabelCollectionNuclei.Update()


caseins = []
waps = []
imgDuplicator = itk.ImageDuplicator.IUC3.New()

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
	maurerSingleNuclei.Update()
	distMap = maurerSingleNuclei.GetOutput()
	
	nucleusObject = statisticsLabelCollectionNuclei.GetOutput().GetLabelObject(l)
	nucleusSize = nucleusObject.GetPhysicalSize()
	nucleusElongation = nucleusObject.GetRegionElongation()
	nucleusIdx = [int(round(v)) for v in nucleusObject.GetCentroid()]
	nucleusMean = nucleusObject.GetMean()
	nucleusSigma = nucleusObject.GetSigma()
	
	print >> nucleiFile, '"%s"' % sys.argv[1], '"%s"' % readerNuclei.GetFileName(), l, nucleusSize, nucleusElongation, nucleusIdx[0], nucleusIdx[1], nucleusIdx[2], nucleusMean, nucleusSigma
	
	# get info for wap
	statsLabelCollectionWap.Update()
	wapObjects = statsLabelCollectionWap.GetOutput()
	for wl in range(1, wapObjects.GetNumberOfObjects()+1) :
		centroid = wapObjects.GetLabelObject(wl).GetCenterOfGravity()
		centerIdx = [int(round(v/s)) for v, s in zip(centroid, spacing)]
		dist = distMap.GetPixel( centerIdx )
		
		if dist > 0:
		  thresholdMaurerSingleNuclei.SetLowerThreshold( dist )
		  shapeLabelCollectionSingleNuclei.Update()
		  innerSize = shapeLabelCollectionSingleNuclei.GetOutput().GetLabelObject(255).GetPhysicalSize()
		  ci = ( nucleusSize - innerSize ) / nucleusSize
	        else:
		  ci = 0
		  
		print >> genesFile, '"%s"' % sys.argv[1], '"%s"' % readerNuclei.GetFileName(), l, '"wap"', centerIdx[0], centerIdx[1], centerIdx[2], dist, ci
		
		# put the segmented wap in a new image
		imgDuplicator.SetInputImage( maskWap.GetOutput() )
		imgDuplicator.Update()
		waps.append( imgDuplicator.GetOutput() )
		
	# get info for cas
	statisticsLabelCollectionCas.Update()
	casObjects = statisticsLabelCollectionCas.GetOutput()
	for wl in range(1, casObjects.GetNumberOfObjects()+1) :
		centroid = casObjects.GetLabelObject(wl).GetCenterOfGravity()
		centerIdx = [int(round(v/s)) for v, s in zip(centroid, spacing)]
		dist = distMap.GetPixel( centerIdx )
		
		if dist > 0:
		  thresholdMaurerSingleNuclei.SetLowerThreshold( dist )
		  shapeLabelCollectionSingleNuclei.Update()
		  innerSize = shapeLabelCollectionSingleNuclei.GetOutput().GetLabelObject(255).GetPhysicalSize()
		  ci = ( nucleusSize - innerSize ) / nucleusSize
	        else:
		  ci = 0
		
		print >> genesFile, '"%s"' % sys.argv[1], '"%s"' % readerNuclei.GetFileName(), l, '"cas"', centerIdx[0], centerIdx[1], centerIdx[2], dist, ci
		
		# put the segmented cas in a new image
		imgDuplicator.SetInputImage( maskCas.GetOutput() )
		imgDuplicator.Update()
		caseins.append( imgDuplicator.GetOutput() )
		

for i, (cas, wap) in enumerate( zip( caseins, waps ) ):
	labelCasNuclei.SetInput( i+1, cas)
	labelWapNuclei.SetInput( i+1, wap)
	
itk.write(overlayWap, readerWap.GetFileName()+"-wap.tif") #, True)
itk.write(overlayCas, readerCas.GetFileName()+"-cas.tif") #, True)

