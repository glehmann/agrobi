
import itk, sys, os
from optparse import OptionParser

parser = OptionParser(usage = '''"Usage: agrobi-site1.py image"
  image: the input image''')

parser.add_option("-v", "--verbose", action="store_true", dest="verbose", help="display infos about progress")
parser.add_option("-V", "--visual-validation", action="store_true", dest="visualValidation", help="write 3 images with the segmented zones overlayed on top of the input channels")
parser.add_option("-t", "--threads", type="int", dest="threads", default=0, help="number of threads to use. Defaults to the number of procs")
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
    nucleiFile.write( '"img" "nucleus" "size" "elongation" "nbSpots" "x" "y" "z"\n' )

if opts.genes == "-":
  genesFile = sys.stdout
else:
  genesFile = file( opts.genes, "a" )
  if os.path.getsize( opts.genes ) == 0:
    genesFile.write( '"img" "nucleus" "x" "y" "z" "px" "py" "pz"\n' )


itk.auto_not_in_place()

# 2 utility objects reused several times later
ni = itk.NearestNeighborInterpolateImageFunction.IUC3D.New()
fill = itk.GrayscaleFillholeImageFilter.IUC2IUC2.New(auto_progress=False)

def attribute_list( i, name ):
  """Returns a list of the specified attributes for the objects in the image.
  
  i: the input LabelImage
  name: the attribute name
  """
  i = itk.image(i)
  relabel = itk.StatisticsRelabelLabelMapFilter[i].New(i, Attribute=name, ReverseOrdering=True, InPlace=False)
  relabel.UpdateLargestPossibleRegion()
  r = relabel.GetOutput()
  l = []
  for i in range(1, r.GetNumberOfLabelObjects()+1):
    l.append( r.GetLabelObject(i).__getattribute__("Get"+name)() )
  return l


def attribute_dict( i, name ):
  """Returns a dict with the attribute values in keys and a list of the corresponding objects in value
  
  i: the input LabelImage
  name: the name of the attribute
  """
  i = itk.image(i)
  relabel = itk.StatisticsRelabelLabelMapFilter[i].New(i, Attribute=name, ReverseOrdering=True, InPlace=False)
  relabel.UpdateLargestPossibleRegion()
  r = relabel.GetOutput()
  d = {}
  for i in range(1, r.GetNumberOfLabelObjects()+1):
    lo = r.GetLabelObject(i)
    v = lo.__getattribute__("Get"+name)()
    l = d.get( v, [] )
    l.append( lo )
    d[v] = l
  return d


def number_of_objects( i ):
  """Returns the number of objets in the image.
  
  i: the input LabelImage
  """
  i.UpdateLargestPossibleRegion()
  i =  itk.image(i)
  return i.GetNumberOfLabelObjects()


def copyImage( f ):
  """Copy an itk.Image object.
  """
  f.UpdateLargestPossibleRegion()
  i = f.GetOutput()
  imgDuplicator = itk.ImageDuplicator[i].New(i)
  imgDuplicator.Update()
  return imgDuplicator.GetOutput()


readerNuclei = itk.lsm(channel=0)

# first reduce the number of pixels by 8, to be able to process the images
downSampleNuclei = itk.ResampleImageFilter.IUC3IUC3.New(readerNuclei, Interpolator=ni.GetPointer())
sMedianNuclei = itk.MedianImageFilter.IUC3IUC3.New(downSampleNuclei)
sSizeNuclei = itk.PhysicalSizeOpeningImageFilter.IUC3IUC3.New(sMedianNuclei, Lambda=100)
sGaussianNuclei = itk.SmoothingRecursiveGaussianImageFilter.IUC3IUC3.New(sSizeNuclei, Sigma=0.2)
sFillNuclei = itk.SliceBySliceImageFilter.IUC3IUC3.New(sGaussianNuclei, Filter=fill.GetPointer())
sGradientNuclei = itk.GradientMagnitudeRecursiveGaussianImageFilter.IUC3IUC3.New(sFillNuclei, Sigma=0.36)
sRatsNuclei = itk.RobustAutomaticThresholdImageFilter.IUC3IUC3IUC3.New(sFillNuclei, sGradientNuclei)
# sOpeningNuclei = itk.BinaryMorphologicalOpeningImageFilter.IUC3IUC3SE3.New(sRatsNuclei, Kernel=itk.strel(3, 5)) # TODO: remove it???
# sBorderNuclei = itk.BinaryShapeOpeningImageFilter.IUC3.New(sOpeningNuclei, Attribute="PhysicalSizeOnBorder", Lambda=50, ReverseOrdering=True)
sBorderNuclei = itk.BinaryShapeOpeningImageFilter.IUC3.New(sRatsNuclei, Attribute="PhysicalSizeOnBorder", Lambda=50, ReverseOrdering=True)
sBSizeNuclei = itk.BinaryShapeOpeningImageFilter.IUC3.New(sBorderNuclei, Attribute="PhysicalSize", Lambda=200)
upSampleNuclei = itk.ResampleImageFilter.IUC3IUC3.New(sBSizeNuclei, UseReferenceImage=True, ReferenceImage=readerNuclei.GetOutput(), Interpolator=ni.GetPointer())
upShapeNuclei = itk.BinaryImageToShapeLabelMapFilter.IUC3LM3.New(upSampleNuclei)

# now segment again, but 1 nucleus by 1 nucleus, to get a mask at the full resolution
singleNuclei = itk.ExtractImageFilter.IUC3IUC3.New(readerNuclei)
singleMedianNuclei = itk.MedianImageFilter.IUC3IUC3.New(singleNuclei)
singleSizeNuclei = itk.PhysicalSizeOpeningImageFilter.IUC3IUC3.New(singleMedianNuclei, Lambda=100)
singleGaussianNuclei = itk.SmoothingRecursiveGaussianImageFilter.IUC3IUC3.New(singleSizeNuclei, Sigma=0.2)
singleFillNuclei = itk.SliceBySliceImageFilter.IUC3IUC3.New(singleGaussianNuclei, Filter=fill.GetPointer())
singleGradientNuclei = itk.GradientMagnitudeRecursiveGaussianImageFilter.IUC3IUC3.New(singleFillNuclei, Sigma=0.36)
singleRatsNuclei = itk.RobustAutomaticThresholdImageFilter.IUC3IUC3IUC3.New(singleFillNuclei, singleGradientNuclei)
# keep only the object which contains the center of the nuclei in the lowres mask
singleShapeNuclei = itk.BinaryImageToShapeLabelMapFilter.IUC3LM3.New(singleRatsNuclei)
# singleLabelNuclei = itk.LabelSelectionLabelMapFilter.LM3.New(singleShapeNuclei, Label=1)
singleRelabelNuclei = itk.ShapeRelabelLabelMapFilter.LM3.New(singleShapeNuclei)
singleLabelNuclei = itk.ShapeKeepNObjectsLabelMapFilter.LM3.New(singleRelabelNuclei, NumberOfObjects=1)
singleChangeRegionNuclei = itk.RegionFromReferenceLabelMapFilter.LM3.New(singleLabelNuclei, readerNuclei, InPlace=False)
singleFullLabelNuclei = itk.LabelMapToLabelImageFilter.LM3IUC3.New(singleChangeRegionNuclei)

naryRelabelNuclei = itk.NaryRelabelImageFilter.IUC3IUC3.New()
if opts.visualValidation:
  overlayNuclei = itk.LabelOverlayImageFilter.IUC3IUC3IRGBUC3.New(readerNuclei, naryRelabelNuclei, UseBackground=True)


readerCENP = itk.lsm(channel=1)
maskCENP = itk.LabelMapMaskImageFilter.LM3IUC3.New(singleLabelNuclei, readerCENP, Crop=True, Label=1, CropBorder=4)
medianCENP = itk.MedianImageFilter.IUC3IUC3.New(maskCENP)
gaussianCENP = itk.SmoothingRecursiveGaussianImageFilter.IUC3IUC3.New(medianCENP, Sigma=0.05)
sizeCENP = itk.PhysicalSizeOpeningImageFilter.IUC3IUC3.New(gaussianCENP, Lambda=2.0)
subCENP = itk.SubtractImageFilter.IUC3IUC3IUC3.New(gaussianCENP, sizeCENP)
thCENP = itk.BinaryThresholdImageFilter.IUC3IUC3.New(subCENP, LowerThreshold=35)
bsizeCENP = itk.BinaryShapeOpeningImageFilter.IUC3.New(thCENP, Attribute="PhysicalSize", Lambda=0.01)

statsCENP = itk.BinaryImageToStatisticsLabelMapFilter.IUC3IUC3LM3.New(bsizeCENP, subCENP)

changeRegionCENP = itk.RegionFromReferenceLabelMapFilter.LM3.New(statsCENP, readerCENP, InPlace=False)
fullLabelCENP = itk.LabelMapToLabelImageFilter.LM3IUC3.New(changeRegionCENP)

outsideMask = itk.BinaryThresholdImageFilter.IUC3IUC3.New(naryRelabelNuclei, UpperThreshold=0)
naryRelabelCENP = itk.NaryRelabelImageFilter.IUC3IUC3.New(outsideMask)
if opts.visualValidation:
  overlayCENP = itk.LabelOverlayImageFilter.IUC3IUC3IRGBUC3.New(readerCENP, naryRelabelCENP, UseBackground=True)

readerNuclei.SetFileName(inputImageName)
readerCENP.SetFileName(inputImageName)

if 'blastocyste' in inputImageName:
  if opts.verbose:
    print >> sys.stderr, "swapping channels"
  readerCENP.SetChannel(0)
  readerNuclei.SetChannel(1)

print itk.size(readerNuclei)

# now set the values which can't be set without having loaded the image
newSpacing = [s*2 for s in itk.spacing(readerNuclei)]
newSize = [s/2 for s in itk.size(readerNuclei)]
downSampleNuclei.SetOutputSpacing(newSpacing)
downSampleNuclei.SetSize(newSize)
singleChangeRegionNuclei.SetRegion( readerNuclei.GetOutput().GetLargestPossibleRegion() )

# create a new image to store the CENP centers, so they can be easily reused to check the distribution
cenpSpotsImg = itk.Image.UC3.New(Regions=itk.size(readerNuclei), Spacing=itk.spacing(readerNuclei))
cenpSpotsImg.Allocate()
cenpSpotsImg.FillBuffer(0)

nucleiLabels = []
CENPLabels = []

upShapeNuclei.UpdateLargestPossibleRegion()


for l in range(5, upShapeNuclei.GetOutput().GetNumberOfLabelObjects()+1, 2):
  if opts.verbose:
    print >> sys.stderr, "processing nucleus", l
  # prepare the region for the extract filter
  region = upShapeNuclei.GetOutput().GetLabelObject(l).GetRegion()
  
  region.PadByRadius(4)
  region.Crop(readerNuclei.GetOutput().GetLargestPossibleRegion())
  print >> sys.stderr, region.GetIndex(), region.GetSize()
  
  singleNuclei.SetExtractionRegion(region)
  # keep only the object which contains the center of the nuclei in the lowres mask
#   singleShapeNuclei.UpdateLargestPossibleRegion()
#   print singleShapeNuclei.GetOutput().GetPixel(itk.physical_point_to_index(upShapeNuclei, upShapeNuclei.GetOutput().GetLabelObject(1).GetCentroid()))
#   singleLabelNuclei.SetLabel( singleShapeNuclei.GetOutput().GetPixel(itk.physical_point_to_index(upShapeNuclei, upShapeNuclei.GetOutput().GetLabelObject(1).GetCentroid())) )
  
  singleShapeNuclei.UpdateLargestPossibleRegion()
  singleLabelNuclei.UpdateLargestPossibleRegion()
#   print "singleShapeNuclei.GetOutput().GetNumberOfLabelObjects():", singleShapeNuclei.GetOutput().GetNumberOfLabelObjects()
#   print "singleLabelNuclei.GetOutput().GetNumberOfLabelObjects():", singleLabelNuclei.GetOutput().GetNumberOfLabelObjects()
#   itk.write(singleRatsNuclei, readerNuclei.GetFileName()+"-nuclei"+str(l)+".tif")
  nucleiLabels.append( copyImage( singleFullLabelNuclei ) )
  itk.write(singleFullLabelNuclei, readerNuclei.GetFileName()+"-nuclei-"+str(l)+".nrrd", True)
  naryRelabelNuclei.SetInput( l-1, nucleiLabels[-1] )
  
  nucleusLabelObject = singleShapeNuclei.GetOutput().GetLabelObject(1)
  # nulceusCentroid = itk.GetCentroid()
#   print l, nucleusLabelObject.GetPhysicalSize(), number_of_objects(statsCENP)
  
  medianCENP.UpdateLargestPossibleRegion()
#   print >> sys.stderr, "size: ", medianCENP.GetOutput().GetLargestPossibleRegion().GetSize()

  CENPLabels.append( copyImage( fullLabelCENP ) )
  naryRelabelCENP.SetInput( l, CENPLabels[-1] )

  for l2 in range(1, statsCENP.GetOutput().GetNumberOfLabelObjects()+1):
    cog = statsCENP.GetOutput().GetLabelObject(l2).GetCenterOfGravity()
    idx = itk.physical_point_to_index(cenpSpotsImg, cog)
    cenpSpotsImg.SetPixel(idx, 255)
    print >> genesFile, '"'+readerNuclei.GetFileName()+'"', l, " ".join(str(i) for i in idx), " ".join(str(i) for i in cog)
    
  # some basic outputs about nuclei
#   print nucleusLabelObject
  print >> nucleiFile, '"'+readerNuclei.GetFileName()+'"', l, nucleusLabelObject.GetPhysicalSize(), nucleusLabelObject.GetRegionElongation(), statsCENP.GetOutput().GetNumberOfLabelObjects(), ' '.join([str(i) for i in itk.physical_point_to_index(upShapeNuclei, nucleusLabelObject.GetCentroid())])
  
if opts.saveSegmentation:
  itk.write(naryRelabelNuclei, readerNuclei.GetFileName()+"-nuclei.nrrd", True)
  itk.write(cenpSpotsImg, readerNuclei.GetFileName()+"-CENP.nrrd", True)

if opts.visualValidation:
  itk.write(overlayNuclei, readerNuclei.GetFileName()+"-nuclei.tif") #, True)
  itk.write(overlayCENP, readerNuclei.GetFileName()+"-CENP.tif") #, True)

  
