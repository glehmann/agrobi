
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


def attributes_list( i, names ):
  """Returns a list of the specified attributes for the objects in the image.
  
  i: the input LabelImage
  name: the attribute name
  """
  i = itk.image(i)
  relabel = itk.StatisticsRelabelLabelMapFilter[i].New(i, Attribute=names[0], ReverseOrdering=True, InPlace=False)
  relabel.UpdateLargestPossibleRegion()
  r = relabel.GetOutput()
  l = []
  for i in range(1, r.GetNumberOfLabelObjects()+1):
    attrs = []
    for name in names :
      attrs.append( r.GetLabelObject(i).__getattribute__("Get"+name)() )
    l.append( tuple( attrs ) )
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

medianNuclei = itk.MedianImageFilter.IUC3IUC3.New(readerNuclei)
gaussianNuclei = itk.SmoothingRecursiveGaussianImageFilter.IUC3IUC3.New(medianNuclei, Sigma=0.2)
fillNuclei = itk.SliceBySliceImageFilter.IUC3IUC3.New(gaussianNuclei, Filter=fill.GetPointer())
sizeNuclei = itk.PhysicalSizeOpeningImageFilter.IUC3IUC3.New(fillNuclei, Lambda=200)
gradientNuclei = itk.GradientMagnitudeRecursiveGaussianImageFilter.IUC3IUC3.New(sizeNuclei, Sigma=0.36)
ratsNuclei = itk.RobustAutomaticThresholdImageFilter.IUC3IUC3IUC3.New(sizeNuclei, gradientNuclei)
# split and labelize the nuclei
maurerNuclei = itk.SignedMaurerDistanceMapImageFilter.IUC3IF3.New(ratsNuclei, UseImageSpacing=True)
watershedNuclei = itk.MorphologicalWatershedImageFilter.IF3IUC3.New(maurerNuclei, Level=2.0, MarkWatershedLine=False) #, FullyConnected=True)
maskWatershedNuclei = itk.MaskImageFilter.IUC3IUC3IUC3.New(watershedNuclei, ratsNuclei)
shapeNuclei = itk.LabelImageToShapeLabelMapFilter.IUC3LM3.New(maskWatershedNuclei)
# remove (again) the objects too small to be a nucleus
sizeOpeningNuclei = itk.ShapeOpeningLabelMapFilter.LM3.New(shapeNuclei, Attribute="PhysicalSize", Lambda=200)
# remove the nucleus on the border - note that they can touch "a little" the border
borderNuclei = itk.ShapeOpeningLabelMapFilter.LM3.New(sizeOpeningNuclei, Attribute="PhysicalSizeOnBorder", Lambda=50, ReverseOrdering=True)
relabelNuclei = itk.ShapeRelabelLabelMapFilter.LM3.New(borderNuclei)
lm2iNuclei = itk.LabelMapToLabelImageFilter.LM3IUC3.New(relabelNuclei)

if opts.visualValidation:
  overlayNuclei = itk.LabelOverlayImageFilter.IUC3IUC3IRGBUC3.New(readerNuclei, lm2iNuclei)


readerCENP = itk.lsm(channel=1)
medianCENP = itk.MedianImageFilter.IUC3IUC3.New(readerCENP)
gaussianCENP = itk.SmoothingRecursiveGaussianImageFilter.IUC3IUC3.New(medianCENP, Sigma=0.05)
sizeCENP = itk.PhysicalSizeOpeningImageFilter.IUC3IUC3.New(gaussianCENP, Lambda=0.36)
subCENP = itk.SubtractImageFilter.IUC3IUC3IUC3.New(gaussianCENP, sizeCENP)
thCENP = itk.BinaryThresholdImageFilter.IUC3IUC3.New(subCENP, LowerThreshold=35)
maskCENP = itk.LabelMapMaskImageFilter.LM3IUC3.New(relabelNuclei, thCENP, Label=0, Negated=True)
statsCENP = itk.BinaryImageToStatisticsLabelMapFilter.IUC3IUC3LM3.New(maskCENP, subCENP)
bsizeCENP = itk.ShapeOpeningLabelMapFilter.LM3.New(statsCENP, Attribute="PhysicalSize", Lambda=0.02)
relabelCENP = itk.ShapeRelabelLabelMapFilter.LM3.New(bsizeCENP)

fullCENP = itk.LabelMapToBinaryImageFilter.LM3IUC3.New(relabelCENP)

outsideMask = itk.BinaryThresholdImageFilter.IUC3IUC3.New(lm2iNuclei, UpperThreshold=0)
naryRelabelCENP = itk.NaryRelabelImageFilter.IUC3IUC3.New(outsideMask, fullCENP)
if opts.visualValidation:
  overlayCENP = itk.LabelOverlayImageFilter.IUC3IUC3IRGBUC3.New(readerCENP, naryRelabelCENP)

readerNuclei.SetFileName(inputImageName)
readerCENP.SetFileName(inputImageName)

if 'blasto' in inputImageName:
  if opts.verbose:
    print >> sys.stderr, "swapping channels"
  readerCENP.SetChannel(0)
  readerNuclei.SetChannel(1)

# create a new image to store the CENP centers, so they can be easily reused to check the distribution
cenpSpotsImg = itk.Image.UC3.New(Regions=itk.size(readerNuclei), Spacing=itk.spacing(readerNuclei))
cenpSpotsImg.Allocate()
cenpSpotsImg.FillBuffer(0)

relabelNuclei.UpdateLargestPossibleRegion()


if opts.saveSegmentation:
  itk.write(lm2iNuclei, readerNuclei.GetFileName()+"-nuclei.nrrd", True)

if opts.visualValidation:
  itk.write(overlayNuclei, readerNuclei.GetFileName()+"-nuclei.tif") #, True)
  itk.write(overlayCENP, readerNuclei.GetFileName()+"-CENP.tif") #, True)

nbOfSpots = {}
for l in range(1, relabelCENP.GetOutput().GetNumberOfLabelObjects()+1):
  cog = relabelCENP.GetOutput().GetLabelObject(l).GetCenterOfGravity()
  idx = itk.physical_point_to_index(cenpSpotsImg, cog)
  cenpSpotsImg.SetPixel(idx, 255)
  nucleus = lm2iNuclei.GetOutput().GetPixel(idx)
  if nucleus == 0:
    print >> sys.stderr, "spot outside a nucleus at position ", idx
  print >> genesFile, '"'+readerNuclei.GetFileName()+'"', nucleus, " ".join(str(i) for i in idx), " ".join(str(i) for i in cog)
  nbOfSpots[nucleus] = nbOfSpots.get(nucleus, 0) + 1
    
for l in range(1, relabelNuclei.GetOutput().GetNumberOfLabelObjects()+1):
  # some basic outputs about nuclei
  nucleusLabelObject = relabelNuclei.GetOutput().GetLabelObject(l)
  print >> nucleiFile, '"'+readerNuclei.GetFileName()+'"', l, nucleusLabelObject.GetPhysicalSize(), nucleusLabelObject.GetRegionElongation(), nbOfSpots.get(l, 0), ' '.join([str(i) for i in itk.physical_point_to_index(relabelNuclei, nucleusLabelObject.GetCentroid())])
  
if opts.saveSegmentation:
  itk.write(cenpSpotsImg, readerNuclei.GetFileName()+"-CENP.nrrd", True)


  
