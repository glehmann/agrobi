/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkCentralIndexMapImageFilter.txx,v $
  Language:  C++
  Date:      $Date: 2004/12/21 22:47:30 $
  Version:   $Revision: 1.12 $

  Copyright (c) Insight Software Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

    This software is distributed WITHOUT ANY WARRANTY; without even 
    the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
    PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkCentralIndexMapImageFilter_txx
#define __itkCentralIndexMapImageFilter_txx

#include "itkCentralIndexMapImageFilter.h"
#include "itkProgressReporter.h"
#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"

namespace itk {

template <class TInputImage, class TOutputImage>
CentralIndexMapImageFilter<TInputImage, TOutputImage>
::CentralIndexMapImageFilter()
{
}


template <class TInputImage, class TOutputImage>
void 
CentralIndexMapImageFilter<TInputImage, TOutputImage>
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
CentralIndexMapImageFilter<TInputImage, TOutputImage>
::EnlargeOutputRequestedRegion(DataObject *)
{
  this->GetOutput()->SetRequestedRegion( this->GetOutput()->GetLargestPossibleRegion() );
}


template<class TInputImage, class TOutputImage>
void
CentralIndexMapImageFilter<TInputImage, TOutputImage>
::GenerateData()
{
  // do not allocate the output now, to decrease the memory usage
  ProgressReporter progress(this, 0, this->GetOutput()->GetRequestedRegion().GetNumberOfPixels()*2);

  // the map which will store the histogram of the distance map, with a high precision
  typedef std::map< InputImagePixelType, unsigned long, std::greater< InputImagePixelType > > InMapType;
  InMapType inMap;

  // fill the map
  long inside = 0;
  ImageRegionConstIterator< InputImageType > iIt( this->GetInput(), this->GetOutput()->GetRequestedRegion() );
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

  // the map which will store the output histogram
  typedef std::map< InputImagePixelType, OutputImagePixelType > OutMapType;
  OutMapType outMap;

  // compute the density of the histogram
  long nb = 0;
  typename InMapType::const_iterator mapIt;
  for( mapIt = inMap.begin(); mapIt != inMap.end(); mapIt++ )
    {
    outMap[ mapIt->first ] = ( inside - nb ) / static_cast< OutputImagePixelType >( inside );
    nb += mapIt->second;
    }

  // now clear the input map, and allocate the output
  inMap.clear();
  this->AllocateOutputs();

  // and fill the output image
  ImageRegionIterator< OutputImageType > oIt( this->GetOutput(), this->GetOutput()->GetRequestedRegion() );
  for( iIt.GoToBegin(), oIt.GoToBegin();
    !iIt.IsAtEnd();
    ++iIt, ++oIt )
    {
    oIt.Set( outMap[ iIt.Get() ] );
    progress.CompletedPixel();
    }
}



template<class TInputImage, class TOutputImage>
void
CentralIndexMapImageFilter<TInputImage, TOutputImage>
::PrintSelf(std::ostream &os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  
//   os << indent << "FullyConnected: "  << m_FullyConnected << std::endl;
}
  
}// end namespace itk
#endif
