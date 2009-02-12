import itk; itk.auto_progress()
itk.auto_not_in_place()

import sys

def mean(s):
  res = float(sum(s))
  return res / len(s)

def median(s):
  res = sorted(s)
  l = len(res)/2
  if len(res)%2:
    return res[l]
  else:
    return mean(res[l-1:l+1])


class ClosestLabelDilateImageFilter(itk.pipeline):
  def __init__(self, *args, **kargs):
    # call the constructor of the superclass but without args and kargs, because the attributes
    # are not all already there!
    # Set/GetRadius() is created in the constructor for example, with the expose() method
    itk.pipeline.__init__(self)
    
    # get the template parameters
    template_parameters = kargs["template_parameters"]
    
    # and store them in an easier way
    ImageType, DistanceMapType = template_parameters
    # the maximum value of the image type
    PixelType, dim = itk.template(ImageType)[1]
    maxValue = itk.NumericTraits[PixelType].max()
    # build the minipipeline
    # use a cast filter to dispatch the input image
    self.connect(itk.CastImageFilter[ImageType, ImageType].New(InPlace=False))
    # dilate the objects in the input image
    self.connect(itk.BinaryThresholdImageFilter[ImageType, ImageType].New(LowerThreshold=0, UpperThreshold=0, InsideValue=0, OutsideValue=maxValue))
    self.connect(itk.BinaryDilateImageFilter[ImageType, ImageType, itk.FlatStructuringElement[dim]].New())
    self.expose("Kernel")
    self.expose("Radius")
    # compute the voronoi map and cast it to a usable image type
    self.append(itk.DanielssonDistanceMapImageFilter[ImageType, DistanceMapType].New(self.filters[0], UseImageSpacing=True, SquaredDistance=False))
    self.append(itk.CastImageFilter[DistanceMapType, ImageType].New(self.filters[-1].GetVoronoiMap()))
    # and mask the voronoi map with the dilated objects
    self.connect(itk.MaskImageFilter[ImageType, ImageType, ImageType].New(Input2=self.filters[2]))
    
    # now we can parse the inputs
    itk.set_inputs(self, args, kargs)

ClosestLabelDilateImageFilter = itk.templated_class(ClosestLabelDilateImageFilter)
ClosestLabelDilateImageFilter.add_image_templates(itk.INTS, itk.REALS)


class ClosestDilateLabelMapFilter(itk.pipeline):
  def __init__(self, *args, **kargs):
    # call the constructor of the superclass but without args and kargs, because the attributes
    # are not all already there!
    # Set/GetRadius() is created in the constructor for example, with the expose() method
    itk.pipeline.__init__(self)
    
    # get the template parameters
    template_parameters = kargs["template_parameters"]
    
    # and store them in an easier way
    ImageType, DistanceMapType = template_parameters
    # the maximum value of the image type
    PixelType, dim = itk.template(ImageType)[1]
    LabelMapType = itk.LabelMap[itk.StatisticsLabelObject[itk.UL, dim]]
    # build the minipipeline
    self.connect(itk.LabelMapToLabelImageFilter[LabelMapType, ImageType].New())
    self.connect(ClosestLabelDilateImageFilter[ImageType, DistanceMapType].New())
    self.expose("Kernel")
    self.expose("Radius")
    self.connect(itk.LabelImageToLabelMapFilter[ImageType, LabelMapType].New())
    
    # now we can parse the inputs
    itk.set_inputs(self, args, kargs)

ClosestDilateLabelMapFilter = itk.templated_class(ClosestDilateLabelMapFilter)
ClosestDilateLabelMapFilter.add_image_templates(itk.INTS, itk.REALS)


#nucleiFileName = "gml16_dapi_1_nuclei.nrrd"
#inputFileName = "gml16_cenp-fitc_dapi_1.zvi"
nucleiFileName = sys.argv[1]
inputFileName = "gml16_cenp-fitc_dapi_%s.zvi" % nucleiFileName[11:-12]

# read the images
cenp = itk.bioformats(inputFileName, ImageType=itk.Image.US3)
nuclei = itk.ImageFileReader.IUC3.New(FileName=nucleiFileName)

# remove the background before masking the image, in case large background outside the nuclei leak inside the nuclei
cenpmedian = itk.MedianImageFilter.IUS3IUS3.New(cenp, Radius=[1, 1, 0])
# add a gaussian filter after the median?
sopening = itk.PhysicalSizeOpeningImageFilter.IUS3IUS3.New(cenpmedian, Lambda=0.2)
sub = itk.SubtractImageFilter.IUS3IUS3IUS3.New(cenpmedian, sopening)

# and dilate the nuclei to be sure to not truncate the spots
dilate = ClosestLabelDilateImageFilter.IUC3IF3.New(nuclei, Kernel=itk.strel(3,[3,3,0]))
# translate the masks to label objects
li2lm = itk.LabelImageToStatisticsLabelMapFilter.IUC3IUS3LM3.New(dilate, sub) 
lmNuclei = li2lm
lmNuclei()

# keep only one nucleus in the pre-filtered cenp image
mask = itk.LabelMapMaskImageFilter.LM3IUS3.New(lmNuclei, sub, Crop=True, CropBorder=1, Label=3)
# search a nice threshold - it will be based on the values of the brightests regions (spots)
rmax = itk.RegionalMaximaImageFilter.IUS3IUS3.New(mask, FullyConnected=True)
rmaxbi2lm = itk.BinaryImageToStatisticsLabelMapFilter.IUS3IUS3LM3.New(rmax, mask, FullyConnected=True)
rmaxRelabel = itk.StatisticsRelabelLabelMapFilter.LM3.New(rmaxbi2lm, Attribute="Maximum")
# and locate the spots
th = itk.BinaryThresholdImageFilter.IUS3IUC3.New(mask)
# convert them to a label map representation
bi2lm = itk.BinaryImageToStatisticsLabelMapFilter.IUC3IUS3LM3.New(th, sub, FullyConnected=True) 

# visualisation
result = itk.LabelMap._3.New(Regions=itk.size(cenp))

# border = itk.LabelContourImageFilter.IUC3IUC3.New(auto_progress=False)#, FullyConnected=True)
# obo3 = itk.ObjectByObjectLabelMapFilter.LM3.New(li2lm, Filter=border, PadSize=6)
border = itk.LabelContourImageFilter.IUC2IUC2.New(auto_progress=False)#, FullyConnected=True)
slice = itk.SliceBySliceImageFilter.IUC3IUC3.New(nuclei, Filter=border)
obo3 = itk.LabelImageToLabelMapFilter.IUC3LM3.New(slice) 
obo3()

dilate2 = itk.BinaryDilateImageFilter.IUC3IUC3SE3.New(Kernel=itk.strel(3,[3,3,1]), auto_progress=False)
dilate3 = itk.BinaryDilateImageFilter.IUC3IUC3SE3.New(dilate2, Kernel=itk.FlatStructuringElement._3.Cross(1), auto_progress=False)
border2 = itk.SubtractImageFilter.IUC3IUC3IUC3.New(dilate3, dilate2, auto_progress=False)
obo2 = itk.ObjectByObjectLabelMapFilter.LM3.New(bi2lm, InputFilter=dilate2, OutputFilter=border2, PadSize=6)

crop = itk.AutoCropLabelMapFilter.LM3.New(result, CropBorder=[30,30,10])
window = itk.IntensityWindowingImageFilter.IUS3IUS3.New(cenp, OutputMinimum=0, OutputMaximum=255)
rescale = itk.CastImageFilter.IUS3IUC3.New(window)
# extract = itk.ExtractImageFilter.IUS3IUS3.New(cenp)
# rescale = itk.RescaleIntensityImageFilter.IUS3IUC3.New(extract)
overlay = itk.LabelMapOverlayImageFilter.LM3IUC3IRGBUC3.New(crop, rescale, NumberOfThreads=1)#, Opacity=0.2)

# segmentation
result2 = itk.Image.UC3.New(Regions=itk.region(cenp), Spacing=itk.spacing(cenp))
result2.Allocate()
result2.FillBuffer(0)

for nucleus in lmNuclei[0]:
  l = nucleus.GetLabel()
  # select the right nucleus
  mask.SetLabel(l)
  # and set the threshold for the cenp spots
  rmaxRelabel()
#  oldValue = itk.range(mask)[1]/5
  v = int( rmaxRelabel[0].GetLabelObject(22).GetMaximum()/2 )
  values = [rmaxRelabel[0].GetLabelObject(i+1).GetMaximum() for i in range(11)]
#  vmean = mean(values) / 4
#  vmedian = median(values) / 4
  th.SetLowerThreshold( int(vmedian) )
  # auto crop filter is actually buggy - lets force it to see that the input has changed
  crop.Modified()
  #
  # visualisation
  obo2()
  result.ClearLabels()
  window.SetWindowMinimum(int(nucleus.GetMinimum()))
  window.SetWindowMaximum(int(nucleus.GetMaximum()))
  #
  for spot in bi2lm[0]:
    cog = spot.GetCenterOfGravity()
    idx = bi2lm[0].TransformPhysicalPointToIndex(cog)
    result2.SetPixel(idx, l)
    print '"'+inputFileName+'"', l, " ".join(str(i) for i in idx), " ".join(str(i) for i in cog)
    #
    result.PushLabelObject(obo2[0].GetLabelObject(spot.GetLabel()))
    #
  # push the nucleus shape to the result label map, to see the nucleus border in the image
  result.PushLabelObject(obo3.GetOutput().GetLabelObject(l))
  # and write the result
  itk.write(overlay, inputFileName + "-%s.tif" % l)

itk.write(result2, inputFileName + "-CENP.nrrd", True)


