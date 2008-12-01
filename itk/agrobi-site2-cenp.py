import itk; itk.auto_progress()
itk.auto_not_in_place()

import sys



nucleiFileName = "gml16_dapi_1_nuclei.nrrd"
inputFileName = "gml16_cenp-fitc_dapi_1.zvi"
nucleiFileName = sys.argv[1]
inputFileName = "gml16_cenp-fitc_dapi_%s.zvi" % nucleiFileName[11:-12]

# read the images
cenp = itk.bioformats(inputFileName, ImageType=itk.Image.US3)
nuclei = itk.ImageFileReader.IUC3.New(FileName=nucleiFileName)

# remove the background before masking the image, in case large background outside the nuclei leak inside the nuclei
median = itk.MedianImageFilter.IUS3IUS3.New(cenp, Radius=[1, 1, 0])
# add a gaussian filter after the median?
sopening = itk.PhysicalSizeOpeningImageFilter.IUS3IUS3.New(median, Lambda=0.2)
sub = itk.SubtractImageFilter.IUS3IUS3IUS3.New(median, sopening)

# translate the masks to label objects
li2lm = itk.LabelImageToStatisticsLabelMapFilter.IUC3IUS3LM3.New(nuclei, median) 
# and dilate the nuclei to be sure to not truncate the spots
dilate = itk.BinaryDilateImageFilter.IUC3IUC3SE3.New(Kernel=itk.strel(3,[3,3,0]))
obo = itk.ObjectImageLabelMapFilter.LM3.New(li2lm, Filter=dilate, PadSize=6)
lmNuclei = obo
lmNuclei.Update()

# keep only one nucleus in the pre-filtered cenp image
mask = itk.LabelMapMaskImageFilter.LM3IUS3.New(lmNuclei, sub, Crop=True, CropBorder=1, Label=3)
# and locate the spots
th = itk.BinaryThresholdImageFilter.IUS3IUC3.New(mask)
# convert them to a label map representation
bi2lm = itk.BinaryImageToStatisticsLabelMapFilter.IUC3IUS3LM3.New(th, median, FullyConnected=True) 

# visualisation
result = itk.LabelMap._3.New(Regions=itk.size(cenp))

# border = itk.LabelContourImageFilter.IUC3IUC3.New(auto_progress=False)#, FullyConnected=True)
# obo3 = itk.ObjectImageLabelMapFilter.LM3.New(li2lm, Filter=border, PadSize=6)
border = itk.LabelContourImageFilter.IUC2IUC2.New(auto_progress=False)#, FullyConnected=True)
slice = itk.SliceBySliceImageFilter.IUC3IUC3.New(nuclei, Filter=border)
obo3 = itk.LabelImageToLabelMapFilter.IUC3LM3.New(slice) 
obo3.Update()

dilate2 = itk.BinaryDilateImageFilter.IUC3IUC3SE3.New(Kernel=itk.strel(3,[3,3,1]), auto_progress=False)
dilate3 = itk.BinaryDilateImageFilter.IUC3IUC3SE3.New(dilate2, Kernel=itk.FlatStructuringElement._3.Cross(1), auto_progress=False)
border2 = itk.SubtractImageFilter.IUC3IUC3IUC3.New(dilate3, dilate2, auto_progress=False)
obo2 = itk.ObjectImageLabelMapFilter.LM3.New(bi2lm, InputFilter=dilate2, OutputFilter=border2, PadSize=6)

crop = itk.AutoCropLabelMapFilter.LM3.New(result, CropBorder=[30,30,10])
window = itk.IntensityWindowingImageFilter.IUS3IUS3.New(cenp, OutputMinimum=0, OutputMaximum=255)
rescale = itk.CastImageFilter.IUS3IUC3.New(window)
# extract = itk.ExtractImageFilter.IUS3IUS3.New(cenp)
# rescale = itk.RescaleIntensityImageFilter.IUS3IUC3.New(extract)
overlay = itk.LabelMapOverlayImageFilter.LM3IUC3IRGBUC3.New(crop, rescale, NumberOfThreads=1)#, Opacity=0.2)

# segmentation
result2 = itk.Image.UC3.New(Regions=itk.size(cenp), Spacing=itk.spacing(cenp))
result2.Allocate()
result2.FillBuffer(0)

for l in range(1, lmNuclei.GetOutput().GetNumberOfLabelObjects()+1):
  # select the right nucleus
  mask.SetLabel(l)
  # and set the threshold for the cenp spots
  maxValue = itk.range(mask)[1]
  th.SetLowerThreshold(maxValue/5)
  # auto crop filter is actually buggy - lets force it to see that the input has changed
  crop.Modified()
  #
  # visualisation
  obo2.UpdateLargestPossibleRegion()
  result.ClearLabels()
  window.SetWindowMinimum(int(li2lm.GetOutput().GetLabelObject(l).GetMinimum()))
  window.SetWindowMaximum(int(li2lm.GetOutput().GetLabelObject(l).GetMaximum()))
  #
  for l2 in range(1, bi2lm.GetOutput().GetNumberOfLabelObjects()+1):
    lo = bi2lm.GetOutput().GetLabelObject(l2)
    cog = lo.GetCenterOfGravity()
    idx = itk.physical_point_to_index(cenp, cog)
    result2.SetPixel(idx, l)
    print '"'+inputFileName+'"', l, " ".join(str(i) for i in idx), " ".join(str(i) for i in cog)
    #
    result.PushLabelObject(obo2.GetOutput().GetLabelObject(l2))
    #
  # push the nucleus shape to the result label map, to see the nucleus border in the image
  result.PushLabelObject(obo3.GetOutput().GetLabelObject(l))
  # and write the result
  itk.write(overlay, inputFileName + "-%s.tif" % l)

itk.write(result2, inputFileName + "-CENP.nrrd", True)


