

import itkExtras

class leica( itkExtras.pipeline ):
	
	class LeicaDescriptor:
		name = ""
		seriesName = ""
		dim = 0
		channels = 0
		pixelSize = 0
		minValue = 0
		maxValue = 0
		size = []
		spacing = []
		origin = []
		
		def __init__( self, name="" ):
			if name:
				self.SetName( name )
		
		def SetName( self, name ):
			
			dimConv = {
				0: 0,
				1: 1,
				2: 3,
				3: 2,
				}
			
			self.name = name.replace(".txt", "")
			f = file( name )
			
			while not f.readline().startswith( "DIMENSION DESCRIPTION #0" ):
				pass
			
			currentDim = -1
			
			l = f.readline()
			while l and not l.startswith("*************************************** NEXT IMAGE *********************************"):
				l = l.replace('\x00', '')
				
				if l.startswith("Pixel Size in Byte:"):
					size = len("Pixel Size in Byte:")
					self.pixelSize = int(l[size:].strip())
					
				elif l.startswith("Max Value:"):
					size = len("Max Value:")
					self.maxValue = float(l[size:].strip())
					
				elif l.startswith("Min Value:"):
					size = len("Min Value:")
					self.minValue = float(l[size:].strip())
					
				elif l.startswith("Number of Dimensions:"):
					size = len("Number of Dimensions:")
					self.dim = int(l[size:].strip())
					self.size = [0] * self.dim
					self.spacing = [0] * self.dim
					self.origin = [0] * self.dim
	
				elif l.startswith("Channels"):
					size = len("Channels")
					self.channels = int(l[size:].strip())
	
				elif l.startswith("Dimension_"):
					size = len("Dimension_")
					currentDim = int(l[size])
					
				elif l.startswith("Logical Size:"):
					size = len("Logical Size:")
					self.size[ dimConv[ currentDim ] ] = int(l[size:].strip())
					
				elif l.startswith("Physical Length:"):
					size = len("Physical Length:")
					self.spacing[ dimConv[ currentDim ] ] = abs(float(l.replace('m', '')[size:].strip())) / self.size[ dimConv[ currentDim ] ] * 1000000
					
				elif l.startswith("Physical Origin:"):
					size = len("Physical Origin:")
					self.origin[ dimConv[ currentDim ] ] = float(l.replace('m', '')[size:].strip())
					
				elif l.startswith("Series Name:"):
					size = len("Series Name:")
					self.seriesName = l[size:].strip()
					
				l = f.readline()
					
					
			while self.spacing[-1] == 0.0:
				self.dim -= 1
				self.spacing.pop()
				self.size.pop()
				self.origin.pop()

	
	def __init__(self, fileName=None, channel=0, ImageType=None ):
		import itk
		itk.pipeline.__init__(self)
		# if ImageType is None, give it a default value
		# this is useful to avoid loading Base while loading this module
		if ImageType == None:
			ImageType = itk.Image.UC3
			
		self.descriptor = leica.LeicaDescriptor()
		self.names = itk.NumericSeriesFileNames.New()
		self.names.SetStartIndex( 0 )
		
		# remove useless SetInput() method created by the constructor of the pipeline class
		#     del self.SetInput
		# set up the pipeline
		self.connect( itk.ImageSeriesReader[ ImageType ].New() )
		self.connect( itk.ChangeInformationImageFilter[ImageType].New( ChangeSpacing=True, ChangeOrigin=True ) )
		# and configure the pipeline
		self.SetChannel( channel )
		if fileName:
			self.SetFileName( fileName )
	
	def SetFileName( self, fileName ):
		self.descriptor.SetName( fileName )
		self.__Update__()
		
	def SetChannel( self, channel ):
		self.channel = channel
		self.__Update__()
	
	def __Update__( self ):
		if self.descriptor.dim == 3:
			self.names.SetEndIndex( self.descriptor.size[2] - 1 )
			nameFormat = "%s_%s_z%%03d_ch%s.tif" % (self.descriptor.name, self.descriptor.seriesName, str(self.channel).zfill(2))
			self.names.SetSeriesFormat( nameFormat )
			self[0].SetFileNames( self.names.GetFileNames() )
			self[-1].SetOutputSpacing( self.descriptor.spacing )
			self[-1].SetOutputOrigin( self.descriptor.origin )
	

del itkExtras
