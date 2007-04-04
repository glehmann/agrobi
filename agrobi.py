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
# an automatic threshold
kappaNuclei = itk.KappaSigmaThresholdImageFilter.IUC3IUC3.New(gaussianNuclei)
# fill the holes
fillHoles2D = itk.GrayscaleFillholeImageFilter.IUC2IUC2.New(auto_progress=False)
fillHolesNuclei = itk.SliceBySliceImageFilter.IUC3IUC3.New(kappaNuclei, Filter=fillHoles2D.GetPointer())
# remove the objects too small to be a nucleus
binarySizeOpeningNuclei = itk.BinaryShapeOpeningImageFilter.IUC3.New(fillHolesNuclei, Lambda=100000)
# we have the mask of our nuclei
maskNuclei = binarySizeOpeningNuclei
# labelize the nuclei
connectedNuclei = itk.ConnectedComponentImageFilter.IUC3IUC3.New(maskNuclei)
labelNuclei = connectedNuclei

overlayNuclei = itk.LabelOverlayImageFilter.IUC3IUC3IRGBUC3.New(readerNuclei, labelNuclei, UseBackground=True)



##########################
# wap
##########################

readerWap = itk.lsm(channel=2)
# mask the cytoplasm: there is too much noise
maskNWap = itk.MaskImageFilter.IUC3IUC3IUC3.New(readerWap, maskNuclei)
# again, remove some noise
medianWap = itk.MedianImageFilter.IUC3IUC3.New(maskNWap)
gaussianWap = itk.SmoothingRecursiveGaussianImageFilter.IUC3IUC3.New(medianWap, Sigma=0.1)
# keep the 4 more visible spots
maxtreeWap = itk.ImageToMaximumTreeFilter.IUC3CTUC3D.New(gaussianWap)
intensityMaxtreeWap = itk.IntensityComponentTreeFilter.CTUC3D.New(maxtreeWap)
keepMaxtreeWap = itk.KeepNLobesComponentTreeFilter.CTUC3D.New(intensityMaxtreeWap, NumberOfLobes=4)
leavesWap = itk.ComponentTreeLeavesToBinaryImageFilter.CTUC3DIUC3.New(keepMaxtreeWap)
maskWap = leavesWap
connectedWap = itk.ConnectedComponentImageFilter.IUC3IUC3.New(leavesWap, FullyConnected=True)
labelWap = connectedWap

overlayWap = itk.LabelOverlayImageFilter.IUC3IUC3IRGBUC3.New(readerWap, labelWap, UseBackground=True)


##########################
# casein
##########################

readerCas = itk.lsm(channel=0)
# again, remove some noise
medianCas = itk.MedianImageFilter.IUC3IUC3.New(readerCas)
gaussianCas = itk.SmoothingRecursiveGaussianImageFilter.IUC3IUC3.New(medianCas, Sigma=0.1)
# keep the 4 more visible spots
maxtreeCas = itk.ImageToMaximumTreeFilter.IUC3CTUC3D.New(gaussianCas)
intensityMaxtreeCas = itk.IntensityComponentTreeFilter.CTUC3D.New(maxtreeCas)
keepMaxtreeCas = itk.KeepNLobesComponentTreeFilter.CTUC3D.New(intensityMaxtreeCas, NumberOfLobes=4)
leavesCas = itk.ComponentTreeLeavesToBinaryImageFilter.CTUC3DIUC3.New(keepMaxtreeCas)
maskCas = leavesCas
connectedCas = itk.ConnectedComponentImageFilter.IUC3IUC3.New(leavesCas, FullyConnected=True)
labelCas = connectedCas

overlayCas = itk.LabelOverlayImageFilter.IUC3IUC3IRGBUC3.New(readerCas, labelCas, UseBackground=True)


labels = itk.NaryBinaryToLabelImageFilter.IUC3IUC3.New(maskNuclei, maskWap, maskCas)



def set_file_name( name ):
	readerNuclei.SetFileName( name )
	readerWap.SetFileName( name )
	readerCas.SetFileName( name )


set_file_name( "wap_cas_20070320_12m2.lsm" )
# set_file_name( "wap_cas_20070320_3m1.lsm" )
# set_file_name( "wap_cas_20070320_9m2.lsm" )
# set_file_name( "wap_cas_20070320_8m2.lsm" )

v = itk.show(labels, MaxOpacity=0.05)

