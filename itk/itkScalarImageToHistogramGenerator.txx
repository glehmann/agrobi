/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkScalarImageToHistogramGenerator.txx,v $
  Language:  C++
  Date:      $Date: 2006/03/14 22:01:52 $
  Version:   $Revision: 1.2 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef _itkScalarImageToHistogramGenerator_txx
#define _itkScalarImageToHistogramGenerator_txx

#include "itkScalarImageToHistogramGenerator.h"


namespace itk { 
namespace Statistics {


template < class TImage, class TMaskImage >
ScalarImageToHistogramGenerator< TImage, TMaskImage >
::ScalarImageToHistogramGenerator() 
{
  m_ImageToListAdaptor = AdaptorType::New();
  m_ListGenerator = ListGeneratorType::New();
  m_HistogramGenerator = GeneratorType::New();
  m_HistogramGenerator->SetListSample( m_ImageToListAdaptor );
  m_HistogramGenerator2 = GeneratorType2::New();
  m_HistogramGenerator2->SetListSample( m_ListGenerator->GetListSample() );
  m_MaskValue = NumericTraits<MaskPixelType>::max();
}



template < class TImage, class TMaskImage >
void
ScalarImageToHistogramGenerator< TImage, TMaskImage >
::SetInput( const ImageType * image ) 
{
  m_ImageToListAdaptor->SetImage( image );
  m_ListGenerator->SetInput( image );
}


template < class TImage, class TMaskImage >
void
ScalarImageToHistogramGenerator< TImage, TMaskImage >
::SetMaskImage( const MaskImageType * image ) 
{
  m_ListGenerator->SetMaskImage( image );
}


template < class TImage, class TMaskImage >
const typename ScalarImageToHistogramGenerator< TImage, TMaskImage >::HistogramType *
ScalarImageToHistogramGenerator< TImage, TMaskImage >
::GetOutput() const
{
  if( m_ListGenerator->GetMaskImage() )
    {
// std::cout << "Generator " << m_HistogramGenerator2->GetOutput()->GetTotalFrequency() << std::endl;
    return m_HistogramGenerator2->GetOutput();
    }
  else
    {
// std::cout << "Adaptor " << m_HistogramGenerator->GetOutput()->GetTotalFrequency() << std::endl;
    return m_HistogramGenerator->GetOutput();
    }
}



template < class TImage, class TMaskImage >
void
ScalarImageToHistogramGenerator< TImage, TMaskImage >
::Compute() 
{
  if( m_ListGenerator->GetMaskImage() )
    {
// std::cout << "Generator" << std::endl;
    m_ListGenerator->SetMaskValue( m_MaskValue );
    m_ListGenerator->Update();
    m_HistogramGenerator2->Update();
    // m_HistogramGenerator->SetListSample( m_ListGenerator->GetListSample() );
    }
  else
    {
// std::cout << "Adaptor" << std::endl;
    // m_HistogramGenerator->SetListSample( m_ImageToListAdaptor );
    m_HistogramGenerator->Update();
    }
  // m_HistogramGenerator->Update();
}



template < class TImage, class TMaskImage >
void
ScalarImageToHistogramGenerator< TImage, TMaskImage >
::SetNumberOfBins( unsigned int numberOfBins ) 
{
  typename HistogramType::SizeType size;
  size.Fill( numberOfBins );
  m_HistogramGenerator->SetNumberOfBins( size );
  m_HistogramGenerator2->SetNumberOfBins( size );
}


template < class TImage, class TMaskImage >
void
ScalarImageToHistogramGenerator< TImage, TMaskImage >
::SetHistogramMin( RealPixelType minimumValue ) 
{
  typedef typename GeneratorType::MeasurementVectorType     MeasurementVectorType;
  MeasurementVectorType minVector;
  minVector[0] = minimumValue;
  m_HistogramGenerator->SetHistogramMin( minVector );
  m_HistogramGenerator2->SetHistogramMin( minVector );
}


template < class TImage, class TMaskImage >
void
ScalarImageToHistogramGenerator< TImage, TMaskImage >
::SetHistogramMax( RealPixelType maximumValue ) 
{
  typedef typename GeneratorType::MeasurementVectorType     MeasurementVectorType;
  MeasurementVectorType maxVector;
  maxVector[0] = maximumValue;
  m_HistogramGenerator->SetHistogramMax( maxVector );
  m_HistogramGenerator2->SetHistogramMax( maxVector );
}



template < class TImage, class TMaskImage >
void
ScalarImageToHistogramGenerator< TImage, TMaskImage >
::SetMarginalScale( double marginalScale )
{
  m_HistogramGenerator->SetMarginalScale( marginalScale );
  m_HistogramGenerator2->SetMarginalScale( marginalScale );
}




template < class TImage, class TMaskImage >
void
ScalarImageToHistogramGenerator< TImage, TMaskImage >
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  os << "ImageToListSample adaptor = " << m_ImageToListAdaptor << std::endl;
  os << "ListGenerator = " << m_ListGenerator << std::endl;
  os << "HistogramGenerator = " << m_HistogramGenerator << std::endl;
  os << "HistogramGenerator2 = " << m_HistogramGenerator2 << std::endl;
}



} // end of namespace Statistics 
} // end of namespace itk

#endif


