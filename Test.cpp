#include "WalshHadamardTransform.h"
#include "InverseWalshHadamardTransform.h"

#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkRescaleIntensityImageFilter.h"

#include <Eigen/Dense>

typedef itk::Image<unsigned char, 2> ImageType;
typedef itk::Image<float, 2> FloatImageType;

static void CreateImage(ImageType* const image, const itk::Index<2>& cornerOfSquare);

template <typename TImage>
static void WriteImage(TImage* const image, const std::string& filename);

static void TestWHT();
static void TestIWHT();

int main(int argc, char* argv[])
{
  WalshHadamardTransform<ImageType> whTransform;
  Eigen::MatrixXf hadamardMatrix = whTransform.ConstructMatrix(2);
  std::cout << "Hadamard: " << std::endl << hadamardMatrix << std::endl;
  whTransform.WriteBasisImages(hadamardMatrix, "BasisImage_");

  Eigen::MatrixXf reversedMatrix = whTransform.ReverseColumns(hadamardMatrix);
  std::cout << "Reversed: " << std::endl << reversedMatrix << std::endl;

  Eigen::MatrixXf sequencyMatrix = whTransform.SequencyOrdering(hadamardMatrix);
  whTransform.WriteBasisImages(sequencyMatrix, "SequencyOrdered_");
  std::cout << "sequencyMatrix: " << std::endl << sequencyMatrix << std::endl;

//   TestWHT();
//   TestIWHT();
  return EXIT_SUCCESS;
}

void CreateImage(ImageType* const image, const itk::Index<2>& cornerOfSquare)
{
  ImageType::IndexType start;
  start.Fill(0);

  ImageType::SizeType size;
  size.Fill(51);

  ImageType::RegionType region(start,size);

  image->SetRegions(region);
  image->Allocate();
  image->FillBuffer(0);

  itk::ImageRegionIterator<ImageType> imageIterator(image,region);

  unsigned int squareSize = 8;

  while(!imageIterator.IsAtEnd())
    {
    if(imageIterator.GetIndex()[0] > cornerOfSquare[0] &&
       imageIterator.GetIndex()[0] < cornerOfSquare[0] + static_cast<int>(squareSize) &&
       imageIterator.GetIndex()[1] > cornerOfSquare[1] &&
       imageIterator.GetIndex()[1] < cornerOfSquare[1] + static_cast<int>(squareSize))
      {
      imageIterator.Set(255);
      }

    ++imageIterator;
    }
}

template <typename TImage>
void WriteImage(TImage* const image, const std::string& filename)
{
  typedef  itk::ImageFileWriter<TImage> WriterType;
  typename WriterType::Pointer writer = WriterType::New();
  writer->SetFileName(filename);
  writer->SetInput(image);
  writer->Update();
}

void TestWHT()
{
  itk::Index<2> cornerOfSquare = {{10,20}};

  ImageType::Pointer image = ImageType::New();
  CreateImage(image, cornerOfSquare);

  typedef WalshHadamardTransform<ImageType> WHTFilterType;
  WHTFilterType whtFilter;
  whtFilter.SetImage(image);
  whtFilter.Compute();

  WHTFilterType::OutputImageType* outputImage = whtFilter.GetOutput();

  WriteImage(outputImage, "dct.mha");

}

void TestIWHT()
{
  itk::Index<2> cornerOfSquare = {{10,20}};

  ImageType::Pointer image = ImageType::New();
  CreateImage(image, cornerOfSquare);

  typedef WalshHadamardTransform<ImageType> WHTFilterType;
  WHTFilterType whtFilter;
  whtFilter.SetImage(image);
  whtFilter.Compute();

  typedef InverseWalshHadamardTransform<FloatImageType> IWHTFilterType;
  IWHTFilterType iwhtFilter;
  iwhtFilter.SetImage(whtFilter.GetOutput());
  iwhtFilter.Compute();

  typedef itk::RescaleIntensityImageFilter<FloatImageType, ImageType> RescaleFilterType;
  RescaleFilterType::Pointer rescaleFilter = RescaleFilterType::New();
  rescaleFilter->SetInput(iwhtFilter.GetOutput());
  rescaleFilter->SetOutputMinimum(0);
  rescaleFilter->SetOutputMaximum(255);
  rescaleFilter->Update();

  WriteImage(rescaleFilter->GetOutput(), "idct.png");

}
