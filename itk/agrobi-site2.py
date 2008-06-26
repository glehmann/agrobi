import itk; itk.auto_progress()
itk.auto_not_in_place()

filename = "20080401_2.lsm"

readerNuclei = itk.lsm(filename)
medianNuclei = itk.MedianImageFilter.IUC3IUC3.New(readerNuclei)
gaussianNuclei = itk.SmoothingRecursiveGaussianImageFilter.IUC3IUC3.New(medianNuclei, Sigma=0.2)

fill2D = itk.GrayscaleFillholeImageFilter.IUC2IUC2.New(auto_progress=False)
fillNuclei = itk.SliceBySliceImageFilter.IUC3IUC3.New(gaussianNuclei, Filter=fill2D)
gradientNuclei = itk.GradientMagnitudeRecursiveGaussianImageFilter.IUC3IUC3.New(fillNuclei, Sigma=0.36)
robustNuclei = itk.RobustAutomaticThresholdImageFilter.IUC3IUC3IUC3.New(fillNuclei, gradientNuclei)
dilateNuclei = itk.BinaryDilateImageFilter.IUC3IUC3SE3.New(robustNuclei, Kernel=itk.strel(3, 3))
li2lmNuclei = itk.BinaryImageToShapeLabelMapFilter.IUC3LM3.New(dilateNuclei)
sizeNuclei = itk.ShapeOpeningLabelMapFilter.LM3.New(li2lmNuclei, Attribute="PhysicalSize", Lambda=200)
relabelNuclei = itk.ShapeRelabelLabelMapFilter.LM3.New(sizeNuclei)


readerCENP = itk.lsm(filename, channel=1)
medianCENP = itk.MedianImageFilter.IUC3IUC3.New(readerCENP)
gaussianCENP = itk.SmoothingRecursiveGaussianImageFilter.IUC3IUC3.New(medianCENP, Sigma=0.05)
maskCENP = itk.LabelMapMaskImageFilter.LM3IUC3.New(relabelNuclei, Label=1)
