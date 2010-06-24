#ifndef __itkErodedVolumeFractionMapImageFilter_h
#define __itkErodedVolumeFractionMapImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkBarrier.h"

namespace itk {

/** \class ErodedVolumeFractionMapImageFilter
 * \brief Compute the eroded volume fraction map from a distance map
 *
 * 
 *
 * \author Gaetan Lehmann. Biologie du Developpement et de la Reproduction, INRA de Jouy-en-Josas, France.
 *
 * \sa SignedMaurerDistanceMapImageFilter, SignedDanielssonDistanceMapImageFilter
 * \ingroup ImageEnhancement  MathematicalMorphologyImageFilters
 */
template<class TInputImage, class TOutputImage>
class ITK_EXPORT ErodedVolumeFractionMapImageFilter : 
    public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Standard class typedefs. */
  typedef ErodedVolumeFractionMapImageFilter Self;
  typedef ImageToImageFilter<TInputImage, TOutputImage>
  Superclass;
  typedef SmartPointer<Self>        Pointer;
  typedef SmartPointer<const Self>  ConstPointer;

  /** Some convenient typedefs. */
  typedef TInputImage InputImageType;
  typedef TOutputImage OutputImageType;

  typedef typename InputImageType::Pointer        InputImagePointer;
  typedef typename InputImageType::ConstPointer   InputImageConstPointer;
  typedef typename InputImageType::RegionType     InputImageRegionType;
  typedef typename InputImageType::PixelType      InputImagePixelType;

  typedef typename OutputImageType::Pointer        OutputImagePointer;
  typedef typename OutputImageType::ConstPointer   OutputImageConstPointer;
  typedef typename OutputImageType::RegionType     OutputImageRegionType;
  typedef typename OutputImageType::PixelType      OutputImagePixelType;
  
  typedef typename OutputImageType::IndexType      IndexType;
  
  /** ImageDimension constants */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      TInputImage::ImageDimension);

  /** Standard New method. */
  itkNewMacro(Self);  

  /** Runtime information support. */
  itkTypeMacro(ErodedVolumeFractionMapImageFilter, 
               ImageToImageFilter);
  

protected:
  ErodedVolumeFractionMapImageFilter();
  ~ErodedVolumeFractionMapImageFilter() {};
  void PrintSelf(std::ostream& os, Indent indent) const;

  /** ErodedVolumeFractionMapImageFilter needs to request enough of the
   * marker image to account for the elementary structuring element.
   * The mask image does not need to be padded. Depending on whether
   * the filter is configured to run a single iteration or until
   * convergence, this method may request all of the marker and mask
   * image be provided. */
  void GenerateInputRequestedRegion();

  /** This filter will enlarge the output requested region to produce
   * all of the output if the filter is configured to run to
   * convergence.
   * \sa ProcessObject::EnlargeOutputRequestedRegion() */
  void EnlargeOutputRequestedRegion(DataObject *itkNotUsed(output));

  void BeforeThreadedGenerateData(void);
  void ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                            int threadId );
  void AfterThreadedGenerateData(void);
  
private:
  ErodedVolumeFractionMapImageFilter(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented

  typedef std::map< InputImagePixelType, unsigned long, std::greater< InputImagePixelType > > InMapType;
  std::vector<InMapType>          m_InMaps;
  // the map which will store the output histogram
  typedef std::map< InputImagePixelType, OutputImagePixelType > OutMapType;
  OutMapType m_OutMap;
  typename Barrier::Pointer m_Barrier;


} ; // end of class

} // end namespace itk
  
#ifndef ITK_MANUAL_INSTANTIATION
#include "itkErodedVolumeFractionMapImageFilter.txx"
#endif

#endif


