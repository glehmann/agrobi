
import sys
import itk
itk.auto_progress()
itk.auto_not_in_place()

##########################
# nuclei
##########################

readerNuclei = itk.lsm(channel=1)
# remove some noise and smooth the image
medianNuclei = itk.MedianImageFilter.IUC3IUC3.New(readerNuclei)
gaussianNuclei = itk.SmoothingRecursiveGaussianImageFilter.IUC3IUC3.New(medianNuclei, Sigma=0.36)
inputNuclei = gaussianNuclei
# an automatic threshold
kappaNuclei = itk.KappaSigmaThresholdImageFilter.IUC3IUC3.New(inputNuclei, Kappa=2.3)
# fill the holes
fillHoles2D = itk.GrayscaleFillholeImageFilter.IUC2IUC2.New(auto_progress=False)
fillHolesNuclei = itk.SliceBySliceImageFilter.IUC3IUC3.New(kappaNuclei, Filter=fillHoles2D.GetPointer())
# remove the objects too small to be a nucleus
binarySizeOpeningNuclei = itk.BinaryShapeOpeningImageFilter.IUC3.New(fillHolesNuclei, Lambda=100000)
# we have the mask of our nuclei
maskNuclei = binarySizeOpeningNuclei
# split and labelize the nuclei
maurerNuclei = itk.SignedMaurerDistanceMapImageFilter.IUC3IF3.New(maskNuclei, UseImageSpacing=True)
watershedNuclei = itk.MorphologicalWatershedImageFilter.IF3IUC3.New(maurerNuclei, Level=1.5, MarkWatershedLine=False) #, FullyConnected=True)
maskWatershedNuclei = itk.MaskImageFilter.IUC3IUC3IUC3.New(watershedNuclei, maskNuclei)
labelSizeOnBorderOpeningNuclei = itk.LabelShapeOpeningImageFilter.IUC3.New(maskWatershedNuclei, Attribute="SizeOnBorder", Lambda=1000, ReverseOrdering=True)
relabelNuclei = itk.ShapeRelabelImageFilter.IUC3.New(labelSizeOnBorderOpeningNuclei)
labelNuclei = relabelNuclei

labelCollectionNuclei = itk.LabelImageToLabelCollectionImageFilter.IUC3LI3.New(labelNuclei, UseBackground=True)
shapeLabelCollectionNuclei = itk.ShapeLabelCollectionImageFilter.LI3.New(labelCollectionNuclei)

overlayNuclei = itk.LabelOverlayImageFilter.IUC3IUC3IRGBUC3.New(readerNuclei, labelNuclei, UseBackground=True)

singleMaskNuclei = itk.BinaryThresholdImageFilter.IUC3IUC3.New(labelNuclei, UpperThreshold=1, LowerThreshold=1)
maurerSingleNuclei = itk.SignedMaurerDistanceMapImageFilter.IUC3IF3.New(singleMaskNuclei, UseImageSpacing=True, SquaredDistance=False, InsideIsPositive=True)
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
# keep the 4 more visible spots
maxtreeWap = itk.ImageToMaximumTreeFilter.IUC3CTUC3D.New(gaussianWap)
intensityMaxtreeWap = itk.LocalIntensityComponentTreeFilter.CTUC3D.New(maxtreeWap)
keepMaxtreeWap = itk.KeepNLobesComponentTreeFilter.CTUC3D.New(intensityMaxtreeWap, NumberOfLobes=4)
leavesWap = itk.ComponentTreeLeavesToBinaryImageFilter.CTUC3DIUC3.New(keepMaxtreeWap)
maskWap = leavesWap
connectedWap = itk.ConnectedComponentImageFilter.IUC3IUC3.New(leavesWap, FullyConnected=True)
labelWap = connectedWap

labelWapNuclei = itk.NaryBinaryToLabelImageFilter.IUC3IUC3.New(singleMaskNuclei, maskWap)
overlayWap = itk.LabelOverlayImageFilter.IUC3IUC3IRGBUC3.New(readerWap, labelWapNuclei, UseBackground=True)

labelCollectionWap = itk.LabelImageToLabelCollectionImageFilter.IUC3LI3.New(labelWap, UseBackground=True)
shapeLabelCollectionWap = itk.ShapeLabelCollectionImageFilter.LI3.New(labelCollectionWap)

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
maskCas = thresholdCas
connectedCas = itk.ConnectedComponentImageFilter.IUC3IUC3.New(maskCas, FullyConnected=True)
labelCas = connectedCas

labelCasNuclei = itk.NaryBinaryToLabelImageFilter.IUC3IUC3.New(singleMaskNuclei, maskCas)
overlayCas = itk.LabelOverlayImageFilter.IUC3IUC3IRGBUC3.New(readerCas, labelCasNuclei, UseBackground=True)

labelCollectionCas = itk.LabelImageToLabelCollectionImageFilter.IUC3LI3.New(labelCas, UseBackground=True)
shapeLabelCollectionCas = itk.ShapeLabelCollectionImageFilter.LI3.New(labelCollectionCas)


labels = itk.NaryBinaryToLabelImageFilter.IUC3IUC3.New(singleMaskNuclei, maskWap, maskCas)



def set_file_name( name ):
	readerNuclei.SetFileName( name )
	readerWap.SetFileName( name )
	readerCas.SetFileName( name )


# set_file_name( "wap_cas_20070320_12m2.lsm" )
# set_file_name( "wap_cas_20070320_3m1.lsm" )
# set_file_name( "wap_cas_20070320_9m2.lsm" )
# set_file_name( "wap_cas_20070320_8m2.lsm" )
# set_file_name( "wap_cas_20070320_6m1.lsm" )

set_file_name( sys.argv[1] )

# v = itk.show(labels, MaxOpacity=0.05)
itk.write(overlayNuclei, readerNuclei.GetFileName()+"-nuclei.tif", True)
itk.write(overlayWap, readerWap.GetFileName()+"-wap.tif", True)
itk.write(overlayCas, readerCas.GetFileName()+"-cas.tif", True)

print '"img"', '"nucleus"', '"gene"', '"x"', '"y"', '"z"', '"dist"', '"ci"'

shapeLabelCollectionNuclei.Update()

ls = [l+1 for l in range(*itk.range(labelNuclei))]

for l in ls :
	# set the label
	singleMaskNuclei.SetUpperThreshold( l )
	singleMaskNuclei.SetLowerThreshold( l )
	
	# update the distance map
	maurerSingleNuclei.Update()
	distMap = maurerSingleNuclei.GetOutput()
	
	nucleusSize = shapeLabelCollectionNuclei.GetOutput().GetLabelObject(l).GetPhysicalSize()
	
	# get info for wap
	shapeLabelCollectionWap.Update()
	wapObjects = shapeLabelCollectionWap.GetOutput()
	for wl in range(1, wapObjects.GetNumberOfObjects()+1) :
		centroid = wapObjects.GetLabelObject(wl).GetCentroid()
		centerIdx = [int(round(v)) for v in centroid]
		dist = distMap.GetPixel( centerIdx )
		
		thresholdMaurerSingleNuclei.SetLowerThreshold( dist )
		shapeLabelCollectionSingleNuclei.Update()
		innerSize = shapeLabelCollectionSingleNuclei.GetOutput().GetLabelObject(255).GetPhysicalSize()
		ci = ( nucleusSize - innerSize ) / nucleusSize
		
		print '"%s"' % readerNuclei.GetFileName(), l, '"wap"', centerIdx[0], centerIdx[1], centerIdx[2], dist, ci
		
	# get info for cas
	shapeLabelCollectionCas.Update()
	casObjects = shapeLabelCollectionCas.GetOutput()
	for wl in range(1, casObjects.GetNumberOfObjects()+1) :
		centroid = casObjects.GetLabelObject(wl).GetCentroid()
		centerIdx = [int(round(v)) for v in centroid]
		dist = distMap.GetPixel( centerIdx )
		
		thresholdMaurerSingleNuclei.SetLowerThreshold( dist )
		shapeLabelCollectionSingleNuclei.Update()
		innerSize = shapeLabelCollectionSingleNuclei.GetOutput().GetLabelObject(255).GetPhysicalSize()
		ci = ( nucleusSize - innerSize ) / nucleusSize
		
		print '"%s"' % readerNuclei.GetFileName(), l, '"cas"', centerIdx[0], centerIdx[1], centerIdx[2], dist, ci
		

