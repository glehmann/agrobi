/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkScalarImageToHistogramGenerator.h,v $
  Language:  C++
  Date:      $Date: 2006/03/14 22:01:52 $
  Version:   $Revision: 1.8 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkScalarImageToHistogramGenerator_h
#define __itkScalarImageToHistogramGenerator_h


#include "itkScalarImageToListAdaptor.h"
#include "itkImageToListGenerator.h"
#include "itkListSampleToHistogramGenerator.h"
#include "itkObject.h"


namespace itk {
namespace Statistics {

template< class TImageType, class TMaskImage = Image< unsigned char, TImageType::ImageDimension > >
class ScalarImageToHistogramGenerator : public Object
{
public:
  /** Standard typedefs */
  typedef ScalarImageToHistogramGenerator  Self ;
  typedef Object Superclass;
  typedef SmartPointer<Self> Pointer;
  typedef SmartPointer<const Self> ConstPointer;

  /** Run-time type information (and related methods). */
  itkTypeMacro(ScalarImageToHistogramGenerator, Object) ;

  /** standard New() method support */
  itkNewMacro(Self) ;

  typedef TImageType                                      ImageType;

  /** Mask Image typedefs */
  typedef TMaskImage                           MaskImageType;
  typedef typename MaskImageType::Pointer      MaskImagePointer ;
  typedef typename MaskImageType::ConstPointer MaskImageConstPointer ;
  typedef typename MaskImageType::PixelType    MaskPixelType ;

  typedef itk::Statistics::ScalarImageToListAdaptor< 
                                              ImageType 
                                                      >   AdaptorType;
  typedef typename AdaptorType::Pointer                   AdaptorPointer;
  typedef typename ImageType::PixelType                   PixelType;
  typedef typename NumericTraits< PixelType >::RealType   RealPixelType;

  typedef ImageToListGenerator< ImageType, MaskImageType > ListGeneratorType;
  typedef typename ListGeneratorType::Pointer ListGeneratorPointer;
  
  typedef itk::Statistics::ListSampleToHistogramGenerator< 
                                  AdaptorType, 
                                  RealPixelType,
                                  DenseFrequencyContainer
                                                          > GeneratorType;

  typedef itk::Statistics::ListSampleToHistogramGenerator< 
                                  typename ListGeneratorType::ListSampleType, 
                                  RealPixelType,
                                  DenseFrequencyContainer
                                                          > GeneratorType2;

  typedef typename GeneratorType::Pointer                   GeneratorPointer;
  typedef typename GeneratorType2::Pointer                   GeneratorPointer2;

  typedef typename GeneratorType::HistogramType             HistogramType;
  typedef typename HistogramType::Pointer                   HistogramPointer;
  typedef typename HistogramType::ConstPointer              HistogramConstPointer;
  
public:

  /** Triggers the Computation of the histogram */
  void Compute( void );

  /** Connects the input image for which the histogram is going to be computed */
  void SetInput( const ImageType * );
  
  /** Connects the input image for which the histogram is going to be computed */
  void SetMaskImage( const MaskImageType * );
  
  /** Return the histogram. o
   \warning This output is only valid after the Compute() method has been invoked 
   \sa Compute */
  const HistogramType * GetOutput() const;
  
  /** Set number of histogram bins */
  void SetNumberOfBins( unsigned int numberOfBins );
 
  /** Set marginal scale value to be passed to the histogram generator */
  void SetMarginalScale( double marginalScale );

  /** Set the minimum value from which the bins will be computed */
  void SetHistogramMin( RealPixelType minimumValue );

  /** Set the maximum value from which the bins will be computed */
  void SetHistogramMax( RealPixelType maximumValue );

  /** Set the pixel value treated as on in the mask. If a mask has been 
   * specified, only pixels with this value will be added to the list sample, if
   * no mask has been specified all pixels will be added as measurement vectors
   * to the list sample. */
  itkSetMacro( MaskValue, MaskPixelType );
  itkGetMacro( MaskValue, MaskPixelType );
  
protected:
  ScalarImageToHistogramGenerator();
  virtual ~ScalarImageToHistogramGenerator() {};
  void PrintSelf(std::ostream& os, Indent indent) const;
  MaskPixelType m_MaskValue;


private:

  AdaptorPointer      m_ImageToListAdaptor;

  ListGeneratorPointer m_ListGenerator;
  
  GeneratorPointer    m_HistogramGenerator;
  GeneratorPointer2    m_HistogramGenerator2;
  
  ScalarImageToHistogramGenerator(const Self&); //purposely not implemented
  void operator=(const Self&); //purposely not implemented
};


} // end of namespace Statistics 
} // end of namespace itk 

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkScalarImageToHistogramGenerator.txx"
#endif

#endif
