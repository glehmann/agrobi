import itk, math, csv

itk.auto_progress(True)

toAnalyze = range(1,26)
imgDir = 'images/apotome/'
csvFile = 'data/nucleus_stats.csv'

## Define pipeline
reader = itk.ImageFileReader.IUC3.New()
labmap = itk.LabelImageToShapeLabelMapFilter.IUC3LM3.New(reader,ComputePerimeter=True)
## End of pipeline definition

toAnalyze = [str(nb) for nb in toAnalyze]
dataFile = open(csvFile,'w')
dataWriter = csv.writer(dataFile)
dataWriter.writerow(['img.nb','label','x','y','z','volume','diameter','surface','shape','elongation','height','a.x','a.y','a.z'])
for nb in toAnalyze:
   filename = imgDir+'gml16_dapi_'+nb+'_nuclei.nrrd'
   print 'Image filename',filename
   reader.SetFileName(filename)
   labmap.UpdateLargestPossibleRegion()
   lmo = labmap.GetOutput()
   n = lmo.GetNumberOfLabelObjects()
   print '\t',n,'nuclei'
   for i in range(n):
       row = [nb,i+1]
       nucleus = lmo.GetNthLabelObject(i)
       centroid = nucleus.GetCentroid()
       for j in range(3):
           row.append(centroid[j])
       volume = nucleus.GetPhysicalSize()
       row.append(volume)
       diameter = 2*math.pow(3*volume/(4*math.pi),1./3)
       row.append(diameter)
       surface = nucleus.GetPerimeter()
       row.append(surface)
       shape = surface/math.pow(volume,2./3)
       row.append(shape)
       elongation = nucleus.GetBinaryElongation()
       row.append(elongation)
       row.append(nucleus.GetRegion().GetSize()[2]*itk.spacing(reader)[2])
       axes = nucleus.GetBinaryPrincipalAxes().GetVnlMatrix()
       for j in range(3):
           row.append(axes.get(0,j))
       dataWriter.writerow(row)
dataFile.close()

