#ifndef InverseWalshHadamardTransform_H
#define InverseWalshHadamardTransform_H

#include "itkImage.h"

template <typename TImage>
class InverseWalshHadamardTransform
{
public:
  typedef itk::Image<float, 2> OutputImageType;

  InverseWalshHadamardTransform();

  void SetImage(TImage* const image);

  void Compute();

  OutputImageType* GetOutput();

private:
  TImage* Image;
  OutputImageType::Pointer OutputImage;
};

#include "InverseWalshHadamardTransform.hpp"

#endif
