#ifndef WalshHadamardTransform_HPP
#define WalshHadamardTransform_HPP

#include "WalshHadamardTransform.h"

#include "itkImageFileWriter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"

#include <iomanip>

template <typename TImage>
WalshHadamardTransform<TImage>::WalshHadamardTransform()
{
  this->Image = NULL;
  this->OutputImage = NULL;
}

template <typename TImage>
Eigen::MatrixXf WalshHadamardTransform<TImage>::ReverseColumns(const Eigen::MatrixXf& matrix)
{
  Eigen::MatrixXf reversed(matrix.rows(), matrix.cols());
  for(int i = 0; i < matrix.cols(); ++i)
  {
    reversed.col(i) = matrix.col(matrix.cols() - 1 - i);
  }

  return reversed;
}

template <typename TImage>
boost::dynamic_bitset<> WalshHadamardTransform<TImage>::BinaryToGray(const boost::dynamic_bitset<>& binary)
{
  return (binary >> 1) ^ binary;
}

template <typename TImage>
void WalshHadamardTransform<TImage>::WriteBasisImages(const Eigen::MatrixXf& hadamardMatrix, const std::string& prefix)
{
  // Compute the outerproduct of all pairs of column in the Hadamard matrix
  for(int i = 0; i < hadamardMatrix.cols(); ++i)
  {
    for(int j = 0; j < hadamardMatrix.cols(); ++j)
    {
      Eigen::MatrixXf imageMatrix = hadamardMatrix.col(i) * hadamardMatrix.col(j).transpose();
      typedef itk::Image<unsigned char, 2> ImageType;
      ImageType::Pointer image = ImageType::New();
      MatrixToImage(imageMatrix,image);

      std::stringstream ss;

      unsigned int basisNumber = i*hadamardMatrix.cols() + j;
      ss << prefix << std::setfill('0') << std::setw(3) << basisNumber << ".png";

      typedef itk::ImageFileWriter<ImageType> WriterType;
      WriterType::Pointer writer = WriterType::New();
      writer->SetFileName(ss.str());
      writer->SetInput(image);
      writer->Update();
    }
  }
}

template <typename TImage>
void WalshHadamardTransform<TImage>::MatrixToImage(const Eigen::MatrixXf& matrix,
                                                   itk::Image<unsigned char, 2>* const image)
{
  itk::Size<2> size = {{matrix.rows(), matrix.cols()}};
  itk::Index<2> corner = {{0,0}};
  itk::ImageRegion<2> region(corner,size);

  image->SetRegions(region);
  image->Allocate();

  itk::ImageRegionIterator<itk::Image<unsigned char, 2> > imageIterator(image, image->GetLargestPossibleRegion());

  while(!imageIterator.IsAtEnd())
  {
    float value = matrix(imageIterator.GetIndex()[0], imageIterator.GetIndex()[1]);
    if(value == 1)
    {
      imageIterator.Set(255);
    }
    else if(value == -1)
    {
      imageIterator.Set(0);
    }
    else
    {
      std::stringstream ss;
      ss << value << " is not a valid Hadamard value!";
      throw std::runtime_error(ss.str());
    }
    ++imageIterator;
  }
}

template <typename TImage>
Eigen::MatrixXf WalshHadamardTransform<TImage>::SequencyOrdering(const Eigen::MatrixXf& matrix)
{
  /* The sequency ordering of the rows of the Walsh matrix can be derived from the
   * ordering of the Hadamard matrix by first applying the bit-reversal permutation
   * and then the Gray code permutation.*/
  Eigen::MatrixXf sequencyMatrix = ReverseColumns(matrix);

  for(int r = 0; r < matrix.rows(); ++r)
  {
    std::stringstream ss;
    for(int c = 0; c < matrix.cols(); ++c)
    {
      if(sequencyMatrix(r,c) == 1)
      {
        ss << "1";
      }
      else if(sequencyMatrix(r,c) == -1)
      {
        ss << "0";
      }
    }
    std::cout << "Binary string " << ss.str() << std::endl;
    boost::dynamic_bitset<> binary(ss.str());
    boost::dynamic_bitset<> graycode = BinaryToGray(binary);
    // Replace the elements in reverse order (since the binary strings were constructed fromt he matrix elements in forward order,
    // but the graycode elements are read with [.] which is reverse order (i.e. least significant bit first)
    for(int c = 0; c < matrix.cols(); ++c)
    {
      if(graycode[c] == 1)
      {
        sequencyMatrix(r,matrix.cols() - 1 - c) = 1;
      }
      else if(graycode[c] == 0)
      {
        sequencyMatrix(r,matrix.cols() - 1 - c) = -1;
      }
    }
  }

  return sequencyMatrix;
}

template <typename TImage>
Eigen::MatrixXf WalshHadamardTransform<TImage>::ConstructMatrix(const unsigned int k)
{
  // This a recursive function that constructs k levels of the matrix
  // H_i = [H_{i-1} H_{i-1}]
  //       [H_{i-1} -H_{i-1}]
  // The size of the matrix returned will be 2^k x 2^k

  if(k == 0)
  {
    Eigen::MatrixXf matrix(1,1);
    matrix(0,0) = 1;
    return matrix;
  }

  Eigen::MatrixXf previousMatrix = ConstructMatrix(k - 1);
  Eigen::MatrixXf matrix(previousMatrix.rows() * 2, previousMatrix.cols() * 2);
  matrix << previousMatrix, previousMatrix, previousMatrix, -1.0f*previousMatrix;
  return matrix;
}

template <typename TImage>
void WalshHadamardTransform<TImage>::Compute()
{
//   if(!this->Image)
//   {
//     throw std::runtime_error("Must specify input image!");
//   }
//
//   if(!this->OutputImage)
//   {
//     this->OutputImage = OutputImageType::New();
//   }
//
//   this->OutputImage->SetRegions(this->Image->GetLargestPossibleRegion());
//   this->OutputImage->Allocate();
//
//   itk::ImageRegionIterator<OutputImageType> outputIterator(this->OutputImage,
//                                                            this->OutputImage->GetLargestPossibleRegion());
//
//   itk::Size<2> fullSize = this->Image->GetLargestPossibleRegion().GetSize();
//
//   while(!outputIterator.IsAtEnd())
//     {
//     float value = 0.0f; // Initialize
//
//     itk::ImageRegionConstIteratorWithIndex<TImage> imageIterator(this->Image,
//                                                                  this->Image->GetLargestPossibleRegion());
//     while(!imageIterator.IsAtEnd())
//       {
//       itk::Index<2> index = imageIterator.GetIndex();
//       // TODO: Write WHT
// //       value += imageIterator.Get() *
// //                cos(itk::Math::pi / fullSize[1] * (index[1] + 0.5f) * outputIterator.GetIndex()[1]) *
// //                cos(itk::Math::pi / fullSize[0] * (index[0] + 0.5f) * outputIterator.GetIndex()[0]);
//       ++imageIterator;
//       }
//
//     outputIterator.Set(value);
//
//     ++outputIterator;
//     }
}

template <typename TImage>
void WalshHadamardTransform<TImage>::SetImage(TImage* const image)
{
  this->Image = image;
}

template <typename TImage>
typename WalshHadamardTransform<TImage>::OutputImageType* WalshHadamardTransform<TImage>::GetOutput()
{
  return this->OutputImage;
}

#endif
