#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCommand.h"
#include "itkSimpleFilterWatcher.h"

#include "itkShapeLabelObject.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkScalarImageToHistogramGenerator.h"
#include <algorithm>
#include "itkMersenneTwisterRandomVariateGenerator.h"
#include "itkLabelMapMaskImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include <fstream>
#include <iostream>

const int dim = 3;
const int hist_size = 1024;
const double hist_step = 0.04;
const int sim = 100;
typedef unsigned char PType;
typedef itk::Image< PType, dim > IType;
typedef itk::Index<dim> IndexType;
typedef itk::ImageFileReader< IType > ReaderType;
typedef itk::ShapeLabelObject< unsigned char, dim > LabelObjectType;
typedef itk::LabelMap< LabelObjectType > LabelMapType;
typedef itk::LabelImageToShapeLabelMapFilter< IType, LabelMapType > LI2LMType;
typedef itk::LabelMapMaskImageFilter< LabelMapType, IType > MaskType;
typedef itk::Image< float, dim > DistImgType;
typedef itk::SignedMaurerDistanceMapImageFilter< IType, DistImgType > DistType;
typedef itk::Statistics::ScalarImageToHistogramGenerator< DistImgType, IType > HistGenType;
typedef HistGenType::HistogramType HistogramType;
typedef HistogramType::Pointer HistogramPointerType;
typedef HistogramType::ConstPointer HistogramConstPointerType;
typedef std::vector<IndexType> IndexListType;
typedef std::vector<double> FListType;
typedef itk::ConnectedComponentImageFilter< IType, IType > ConnType;

static itk::Statistics::MersenneTwisterRandomVariateGenerator::Pointer randomGenerator;

long random( long lowest, long highest )
{
  if( randomGenerator.IsNull() )
    {
    randomGenerator = itk::Statistics::MersenneTwisterRandomVariateGenerator::New();
    randomGenerator->Initialize();
    }
  return randomGenerator->GetIntegerVariate( highest - lowest ) + lowest;
//   return lowest + static_cast< long >( ( ( highest - lowest ) + 1 ) *( rand() / ( RAND_MAX + 1.0 ) ) );
/*  return lowest+long(((highest-lowest)+1)*rand()/(RAND_MAX + 1.0));*/
}

double dist( const FListType & F, const FListType & Fmean )
{
double d = 0;
for( int i=0; i<hist_size; i++ )
  {
  d = std::max( d, std::abs( F[i] - Fmean[i] ) );
  }
  return d;
}

FListType festimator( const IType * input, const IType * spots, PType label, int nbOfRandIdx, long & realNbOfIdx )
{
  // convert it to a label map representation
  LI2LMType::Pointer li2lm = LI2LMType::New();
  li2lm->SetInput( input );
  li2lm->Update();

  // get the input objet
  const LabelObjectType * lo = li2lm->GetOutput()->GetLabelObject( label );
//   std::cout << "lo->GetSize(): " << lo->GetSize() << std::endl;
  
  // create a new image with the size of the bounding box of the label object where we will draw the points.
  IType::Pointer img;
  if( spots )
    {
    // keep only the interesting part in the image - the one inside the nucleus
    MaskType::Pointer mask = MaskType::New();
    mask->SetInput( li2lm->GetOutput() );
    mask->SetFeatureImage( spots );
    mask->SetLabel( label );
    mask->SetCrop( true );
    mask->Update();
    img = const_cast<IType*>( mask->GetOutput() );
    }
  else
    {
    img = IType::New();
    img->SetRegions( lo->GetRegion() );
    img->SetSpacing( input->GetSpacing() );
    img->Allocate();
    img->FillBuffer( 0 );
    }

  // random generated pixels
  for( int i=0; i<nbOfRandIdx; i++ )
    {
    IndexType idx = lo->GetIndex( random( 0, lo->GetSize() ) );
    if( img->GetPixel( idx ) == 255 )
      {
      std::cerr << "Warning, pixel already set." << std::endl;
      }
    img->SetPixel( idx, 255 );
    }

  // count the spots
  // std::cerr << "Warning, pixel already set." << std::endl;
  ConnType::Pointer connected = ConnType::New();
  connected->SetInput( img );
  connected->Update();
  realNbOfIdx = connected->GetObjectCount();
   
  // generate the distance map
  // std::cerr << "generate the distance map" << std::endl;
  DistType::Pointer dist = DistType::New();
  dist->SetInput( img );
  dist->SetUseImageSpacing( true );
  dist->SetSquaredDistance( false );
//   itk::SimpleFilterWatcher watcher(dist, "filter");
  dist->Update();
  
 
//   // crop the input label image to the bounding box of the object of interest so we can use it to mask the image
//   typedef itk::RegionFromReferenceLabelMapFilter< LabelMapType > CropType;
//   CropType::Pointer crop = CropType::New();
//   crop->SetInput( li2lm->GetOutput() );
//   crop->SetInPlace( false );  // to not alter the input image
//   crop->SetReference( img );

  // compute histogram
  // std::cerr << "compute histogram" << std::endl;
  HistGenType::Pointer histGen = HistGenType::New();
  histGen->SetInput( dist->GetOutput() );
  histGen->SetMaskImage( input );
  histGen->SetMaskValue( label );
  histGen->SetHistogramMin( 0 );
  histGen->SetHistogramMax( hist_size*hist_step );
  histGen->SetNumberOfBins( hist_size );
  histGen->Compute();
//   histGen->Print(std::cerr);

  HistogramConstPointerType histogram = histGen->GetOutput();

  // compute the ratio for all the distances
  // std::cerr << "compute the ratio for all the distances" << std::endl;
  FListType ratios;
  ratios.reserve( hist_size );
  long current_size = 0;
  for( int i=0; i<hist_size; i++ )
    {
//     std::cout << histogram->GetMeasurement( i, 0 ) << ": " << histogram->GetFrequency( i ) << std::endl;
    current_size += static_cast<long>( histogram->GetFrequency( i ) );
    ratios.push_back( current_size / (double)lo->GetSize() );
    }
  return ratios;
}





int main(int argc, char * argv[])
{

  if( argc != 3 )
    {
    std::cerr << "usage: " << argv[0] << " input label" << std::endl;
    // std::cerr << "  : " << std::endl;
    exit(1);
    }
  // a random seed: the time
//   srand((unsigned)time(0));

  std::cerr << argv[1] << "  nucleus: " << argv[2] << std::endl;
  // store the filename to reuse it later
  std::string basename = argv[1];
  // and the label to analyse
  std::string labelStr = argv[2];
  int label = atoi( labelStr.c_str() );

  // lets load the label image of the nuclei
  std::string nucleiFile = basename + "-nuclei.nrrd";
  ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName( nucleiFile.c_str() );
  reader->Update();

  // and the image of the spots
  std::string spotsFile = basename + "-CENP.nrrd";
  ReaderType::Pointer reader2 = ReaderType::New();
  reader2->SetFileName( spotsFile.c_str() );
  reader2->Update();

  // compute estimation of F for the indexes in the nucleus
  long nbOfSpots = 0;
  FListType F = festimator( reader->GetOutput(), reader2->GetOutput(), label, 0, nbOfSpots );

  // compute Fmean
  FListType Fmean( hist_size, 0 );
  std::vector<FListType> Fsims;

  for( int iter=0; iter<sim; iter++ )
    {
    std::cerr << iter << " "; // << std::endl;
    long dummy;
    FListType Fsim = festimator( reader->GetOutput(), NULL, label, nbOfSpots, dummy );
    Fsims.push_back( Fsim );
  //   std::cout << std::endl;
    for( int i=0; i<hist_size; i++ )
      {
      Fmean[i] += Fsim[i] / sim;
      }
    }
  std::cerr << std::endl;

  // to see what the F functions looks like
   std::string outFName = basename+"-"+labelStr+"-F.txt";
   std::ofstream outF(outFName.c_str(), std::ios::out);
   std::string outFMeanName = basename+"-"+labelStr+"-FMean.txt";
   std::ofstream outFMean(outFMeanName.c_str(), std::ios::out);
  for( int i=0; i<hist_size; i++ )
    {
    outF << F[i] << std::endl;
    outFMean << Fmean[i] << std::endl;
    }
  
  double d = dist( F, Fmean);
//   std::cout << "d: " << d << std::endl;

  // compute the distances for all the simulated F
  std::vector<double> dists;
  for( int i=0; i<sim; i++ )
    {
    dists.push_back( dist( Fsims[i], Fmean ) );
//     std::cout << dists[i] << std::endl;
    }
  // sort them
  std::sort( dists.begin(), dists.end() );

  int i;
  for( i=0; i<sim && d > dists[i]; i++ ); // search the index
//   std::cout << "i: " << i << std::endl;
  std::cout << "\"" << basename << "\" "
            << "\"" << labelStr << "\" "
            << ( i + 1.0 ) / ( sim + 1.0 ) << std::endl;

  return 0;
}

