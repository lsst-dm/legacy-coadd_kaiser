// -*- LSST-C++ -*-
/**
* @brief define addToMaskedImage
*
* @file
*
* @author Russell Owen.
*
* @todo: move to a better location
*/
#include "lsst/coadd/kaiser/addToMaskedImage.h"

namespace pexExcept = lsst::pex::exceptions;

/**
* @brief add good pixels from one MaskedImage to another MaskedImage
*
* @todo: move to a better location (afw?)
*
* @throw pexExcept::InvalidParameterException if the image sizes do not match.
*/
template <typename ImagePixelT, typename MaskPixelT, typename VariancePixelT>
void lsst::coadd::kaiser::addToMaskedImage(
    lsst::afw::image::MaskedImage<ImagePixelT, MaskPixelT, VariancePixelT> &outMaskedImage, ///< masked image to be added to
    lsst::afw::image::MaskedImage<ImagePixelT, MaskPixelT, VariancePixelT> const &inMaskedImage,    ///< masked image to add
    MaskPixelT const badPixelMask   ///< skip input pixel if input mask | badPixelMask != 0
) {
    typedef lsst::afw::image::MaskedImage<ImagePixelT, MaskPixelT, VariancePixelT> MaskedImageT;

    if (inMaskedImage.getDimensions() != outMaskedImage.getDimensions()) {
        throw LSST_EXCEPT(pexExcept::InvalidParameterException, "MaskedImage sizes do not match");
    }

    // Set the pixels row by row, to avoid repeated checks for end-of-row
    for (int y = 0, endY = inMaskedImage.getHeight(); y != endY; ++y) {
        typename MaskedImageT::const_x_iterator inPtr = inMaskedImage.row_begin(y);
        typename MaskedImageT::const_x_iterator inEndPtr = inMaskedImage.row_end(y);
        typename MaskedImageT::x_iterator outPtr = outMaskedImage.row_begin(y);
        for (; inPtr != inEndPtr; ++inPtr, ++outPtr) {
            if ((inPtr.mask() & badPixelMask) == 0) {
                *outPtr += *inPtr;
            }
        }
    }
}

//
// Explicit instantiations
//
#define INSTANTIATE(ImagePixelT) \
    template void lsst::coadd::kaiser::addToMaskedImage<ImagePixelT, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>( \
        lsst::afw::image::MaskedImage<ImagePixelT, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel> &outMaskedImage, \
        lsst::afw::image::MaskedImage<ImagePixelT, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel> const &inMaskedImage, \
        lsst::afw::image::MaskPixel const badPixelMask \
    );

INSTANTIATE(boost::uint16_t);
INSTANTIATE(int);
INSTANTIATE(float);
INSTANTIATE(double);
