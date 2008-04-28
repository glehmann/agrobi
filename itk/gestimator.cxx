#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkCommand.h"
#include "itkSimpleFilterWatcher.h"

#include "itkShapeLabelObject.h"
#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkBinaryImageToShapeLabelMapFilter.h"
#include "itkSignedMaurerDistanceMapImageFilter.h"
#include "itkScalarImageToHistogramGenerator.h"
#include <algorithm>
#include "itkMersenneTwisterRandomVariateGenerator.h"
#include "itkLabelMapMaskImageFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include <fstream>
#include <iostream>

const int dim = 3;
const int sim = 999;
typedef unsigned char PType;
typedef itk::Image< PType, dim > IType;
typedef itk::Index<dim> IndexType;
typedef itk::Point<double, dim> CIndexType;
typedef itk::ImageFileReader< IType > ReaderType;
typedef itk::ShapeLabelObject< unsigned char, dim > LabelObjectType;
typedef itk::LabelMap< LabelObjectType > LabelMapType;
typedef itk::LabelImageToShapeLabelMapFilter< IType, LabelMapType > LI2LMType;
typedef itk::BinaryImageToShapeLabelMapFilter< IType, LabelMapType > BI2LMType;
typedef itk::LabelMapMaskImageFilter< LabelMapType, IType > MaskType;
typedef itk::Image< float, dim > DistImgType;
typedef std::vector<double> GListType;
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

double dist( const GListType & G, const GListType & Gmean )
{
  double d = 0;
  int i = 0;
  int j = 0;
  while( i<G.size() && j<Gmean.size() )
    {
    if( G[i] == Gmean[j] )
      {
      i++;
      j++;
      }
    else if( G[i] > Gmean[j] )
      {
      while( G[i] > Gmean[j] && j<Gmean.size() )
        {
        j++;
        }
      }
    else if( G[i] < Gmean[j] )
      {
      while( G[i] < Gmean[j] && i<G.size() )
        {
        i++;
        }
      }
    d = std::max( std::abs( ( i + 1 ) / static_cast<double>( G.size() ) - ( j + 1 ) /static_cast<double>( Gmean.size() ) ), d );
    }

  return d;
}

GListType gestimator( const LabelMapType * input, const IType * spots, PType label, int nbOfRandIdx, long & realNbOfIdx )
{
  // get the input objet
  const LabelObjectType * lo = input->GetLabelObject( label );
//   std::cout << "lo->GetSize(): " << lo->GetSize() << std::endl;
  
  // get the spots from the image
  std::vector< CIndexType > spotList;
  if( spots )
    {
    // keep only the interesting part in the image - the one inside the nucleus
    MaskType::Pointer mask = MaskType::New();
    mask->SetInput( input );
    mask->SetFeatureImage( spots );
    mask->SetLabel( label );
    mask->SetCrop( true );
    BI2LMType::Pointer bi2lm = BI2LMType::New();
    bi2lm->SetInput( mask->GetOutput() );
    bi2lm->Update();
    for( int i=1; i<=bi2lm->GetOutput()->GetNumberOfLabelObjects(); i++ )
      {
      spotList.push_back( bi2lm->GetOutput()->GetLabelObject( i )->GetCentroid() );
      }
    }

  // random generated pixels
  for( int i=0; i<nbOfRandIdx; i++ )
    {
    IndexType idx = lo->GetIndex( random( 0, lo->GetSize() ) );
    CIndexType p;
    input->TransformIndexToPhysicalPoint( idx, p );
    spotList.push_back( p );
    }

  // count the spots
  realNbOfIdx = spotList.size();

  // measure the distances between spots
  std::vector< std::vector< double > > distances;
  distances.resize( spotList.size() );
  for( int i=0; i<spotList.size(); i++ )
    {
    distances[i].resize( spotList.size() );
    }
  for( int i=0; i<spotList.size(); i++ )
    {
    distances[i][i] = 0;
    for( int j=i+1; j<spotList.size(); j++ )
      {
      double dist = 0.0;
      for( int d=0; d<3; d++ )
        {
        dist += std::pow( spotList[i][d] - spotList[j][d], 2 );
        }
      dist = std::sqrt( dist );
      distances[i][j] = dist;
      }
    }

  GListType G;
  for( int i=0; i<spotList.size(); i++ )
    {
    double dist = itk::NumericTraits<double>::max();
    for( int j=i+1; j<spotList.size(); j++ )
      {
      if( distances[i][j] > 0 )
        {
        dist = std::min(dist, distances[i][j]);
        }
      }
      G.push_back( dist );
    }
  return G;
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

  // convert it to a label map representation
  LI2LMType::Pointer li2lm = LI2LMType::New();
  li2lm->SetInput( reader->GetOutput() );
  li2lm->Update();

  // compute estimation of G for the indexes in the nucleus
  long nbOfSpots = 0;
  GListType G = gestimator( li2lm->GetOutput(), reader2->GetOutput(), label, 0, nbOfSpots );
//   std::cerr << "  nbOfSpots: " << nbOfSpots << std::endl;

  // compute Gmean
  GListType Gmean;
  std::vector<GListType> Gsims;

  for( int iter=0; iter<sim; iter++ )
    {
//     std::cerr << iter << " "; // << std::endl;
    long dummy;
    GListType Gsim = gestimator( li2lm->GetOutput(), NULL, label, nbOfSpots, dummy );
    std::sort( Gsim.begin(), Gsim.end() );
    Gsims.push_back( Gsim );
  //   std::cout << std::endl;
    for( int i=0; i<Gsim.size(); i++ )
      {
      Gmean.push_back( Gsim[i] );
      }
    }
//   std::cerr << std::endl;
  std::sort( Gmean.begin(), Gmean.end() );

/*  // to see what the F functions looks like
   std::string outFName = basename+"-"+labelStr+"-F.txt";
   std::ofstream outF(outFName.c_str(), std::ios::out);
   std::string outFMeanName = basename+"-"+labelStr+"-FMean.txt";
   std::ofstream outFMean(outFMeanName.c_str(), std::ios::out);
  for( int i=0; i<hist_size; i++ )
    {
    outF << F[i] << std::endl;
    outFMean << Fmean[i] << std::endl;
    }
*/
  double d = dist( G, Gmean);
//   std::cout << "d: " << d << std::endl;

  // compute the distances for all the simulated F
  std::vector<double> dists;
  for( int i=0; i<sim; i++ )
    {
    dists.push_back( dist( Gsims[i], Gmean ) );
//     std::cout << dists[i] << std::endl;
    }
  // sort them
  std::sort( dists.begin(), dists.end() );

  int i;
  for( i=0; i<sim && d > dists[i]; i++ ); // search the index
//   std::cout << "i: " << i << std::endl;
  std::cout << "\"" << basename << "\" "
            << "\"" << labelStr << "\" "
            << ( sim + 1.0 - i ) / ( sim + 1.0 ) << std::endl;

  return 0;
}

