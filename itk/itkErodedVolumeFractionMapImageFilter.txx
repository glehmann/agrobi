/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkErodedVolumeFractionMapImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2004/12/21 22:47:30 $
  Version:   $Revision: 1.12 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

    This software is distributed WITHOUT ANY WARRANTY; without even 
    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
    PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkErodedVolumeFractionMapImageFilter_txx
#define __itkErodedVolumeFractionMapImageFilter_txx

#include "itkErodedVolumeFractionMapImageFilter.h"
#include "itkProgressReporter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

namespace itk {

template <class TInputImage, class TOutputImage>
ErodedVolumeFractionMapImageFilter<TInputImage, TOutputImage>
::ErodedVolumeFractionMapImageFilter()
{
}


template <class TInputImage, class TOutputImage>
void 
ErodedVolumeFractionMapImageFilter<TInputImage, TOutputImage>
::GenerateInputRequestedRegion()
{
  // call the superclass' implementation of this method
  Superclass::GenerateInputRequestedRegion();
  
  // get pointers to the inputs
  InputImagePointer  inputPtr = 
    const_cast< InputImageType * >( this->GetInput() );
  
  if ( !inputPtr )
    { return; }

  // We need to
  // configure the inputs such that all the data is available.
  //
  inputPtr->SetRequestedRegion(inputPtr->GetLargestPossibleRegion());
}


template <class TInputImage, class TOutputImage>
void 
ErodedVolumeFractionMapImageFilter<TInputImage, TOutputImage>
::EnlargeOutputRequestedRegion(DataObject *)
{
  this->GetOutput()->SetRequestedRegion( this->GetOutput()->GetLargestPossibleRegion() );
}


template<class TInputImage, class TOutputImage>
void
ErodedVolumeFractionMapImageFilter<TInputImage, TOutputImage>
::BeforeThreadedGenerateData(void)
{
  typename TOutputImage::Pointer output = this->GetOutput();
  typename TInputImage::ConstPointer input = this->GetInput();

  long nbOfThreads = this->GetNumberOfThreads();
  if( itk::MultiThreader::GetGlobalMaximumNumberOfThreads() != 0 )
    {
    nbOfThreads = vnl_math_min( this->GetNumberOfThreads(), itk::MultiThreader::GetGlobalMaximumNumberOfThreads() );
    }
  // number of threads can be constrained by the region size, so call the SplitRequestedRegion
  // to get the real number of threads which will be used
  typename TOutputImage::RegionType splitRegion;  // dummy region - just to call the following method
  nbOfThreads = this->SplitRequestedRegion(0, nbOfThreads, splitRegion);

  m_Barrier = Barrier::New();
  m_Barrier->Initialize( nbOfThreads );
  m_InMaps.clear();
  m_InMaps.resize( nbOfThreads );
}

template<class TInputImage, class TOutputImage>
void
ErodedVolumeFractionMapImageFilter<TInputImage, TOutputImage>
::ThreadedGenerateData(const OutputImageRegionType& outputRegionForThread,
                       int threadId)
{
  // do not allocate the output now, to decrease the memory usage
  ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels()*2);

  InMapType & inMap = m_InMaps[threadId];

  // fill the map
  long inside = 0;
  ImageRegionConstIterator< InputImageType > iIt( this->GetInput(), outputRegionForThread );
  for( iIt.GoToBegin(); !iIt.IsAtEnd(); ++iIt )
    {
    const InputImagePixelType & v = iIt.Get();
    inMap[ v ]++;
    if( v > 0 )
      {
      inside++;
      }
    progress.CompletedPixel();
    }

  m_Barrier->Wait();

  if( threadId == 0 )
    {
    m_OutMap.clear();
    typename InMapType::const_iterator mapIt;

    // group the histograms
    for( unsigned int i=1; i<m_InMaps.size(); i++ )
      {
      InMapType & inMapT = m_InMaps[i];
      for( mapIt = inMapT.begin(); mapIt != inMapT.end(); mapIt++ )
        {
        inMap[ mapIt->first ] += mapIt->second;
        if( mapIt->first > 0 )
          {
          inside += mapIt->second;
          }
        }
      }
    // compute the density of the histogram
    long nb = 0;
    for( mapIt = inMap.begin(); mapIt != inMap.end(); mapIt++ )
      {
      m_OutMap[ mapIt->first ] = ( inside - nb ) / static_cast< OutputImagePixelType >( inside );
      nb += mapIt->second;
      }
    this->AllocateOutputs();
    }

  m_Barrier->Wait();

  // now clear the input map
  inMap.clear();

  // and fill the output image
  ImageRegionIterator< OutputImageType > oIt( this->GetOutput(), outputRegionForThread );
  for( iIt.GoToBegin(), oIt.GoToBegin();
    !iIt.IsAtEnd();
    ++iIt, ++oIt )
    {
    oIt.Set( m_OutMap[ iIt.Get() ] );
    progress.CompletedPixel();
    }
}


template<class TInputImage, class TOutputImage>
void
ErodedVolumeFractionMapImageFilter<TInputImage, TOutputImage>
::AfterThreadedGenerateData(void)
{
  m_InMaps.clear();
  m_OutMap.clear();
}


template<class TInputImage, class TOutputImage>
void
ErodedVolumeFractionMapImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  
//   os << indent << "FullyConnected: "  << m_FullyConnected << std::endl;
}
  
}// end namespace itk
#endif
