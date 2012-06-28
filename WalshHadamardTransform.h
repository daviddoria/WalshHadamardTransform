#ifndef WalshHadamardTransform_H
#define WalshHadamardTransform_H

#include "itkImage.h"

#include <Eigen/Dense>

#include "boost/dynamic_bitset.hpp"

template <typename TImage>
class WalshHadamardTransform
{
public:
  typedef itk::Image<float, 2> OutputImageType;

  WalshHadamardTransform();

  void SetImage(TImage* const image);

  Eigen::MatrixXf ConstructMatrix(const unsigned int k);

  Eigen::MatrixXf SequencyOrdering(const Eigen::MatrixXf& matrix);

  Eigen::MatrixXf ReverseColumns(const Eigen::MatrixXf& matrix);

  void WriteBasisImages(const Eigen::MatrixXf& matrix, const std::string& prefix);

  void MatrixToImage(const Eigen::MatrixXf& matrix, itk::Image<unsigned char, 2>* const image);

  boost::dynamic_bitset<> BinaryToGray(const boost::dynamic_bitset<>& binary);

  void Compute();

  OutputImageType* GetOutput();

private:
  TImage* Image;
  OutputImageType::Pointer OutputImage;
};

#include "WalshHadamardTransform.hpp"

#endif
