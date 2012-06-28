#ifndef InverseWalshHadamardTransform_HPP
#define InverseWalshHadamardTransform_HPP

#include "InverseWalshHadamardTransform.h"

#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIteratorWithIndex.h"

template <typename TImage>
InverseWalshHadamardTransform<TImage>::InverseWalshHadamardTransform()
{
  this->Image = NULL;
  this->OutputImage = NULL;
}


template <typename TImage>
void InverseWalshHadamardTransform<TImage>::Compute()
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
//       // TODO: Write IWHT
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
void InverseWalshHadamardTransform<TImage>::SetImage(TImage* const image)
{
  this->Image = image;
}

template <typename TImage>
typename InverseWalshHadamardTransform<TImage>::OutputImageType*
InverseWalshHadamardTransform<TImage>::GetOutput()
{
  return this->OutputImage;
}

#endif
