// -*- LSST-C++ -*-
/**
* @brief define addToCoadd
*
* @file
*
* @author Russell Owen.
*
* @todo: move to a better location
*/
#include "lsst/pex/exceptions.h"
#include "lsst/coadd/kaiser/addToCoadd.h"

namespace pexExcept = lsst::pex::exceptions;

/**
* @brief add good pixels from an image to a coadd and associated depth map
*
* Depth map: the value of each pixel is the number of good image pixels
* that contributed to the associated pixel of the coadd.
*
* @todo: move to a better location (afw?)
*
* @throw pexExcept::InvalidParameterException if the image sizes do not match.
*/
template <typename ImagePixelT, typename MaskPixelT, typename VariancePixelT>
void lsst::coadd::kaiser::addToCoadd(
    lsst::afw::image::MaskedImage<ImagePixelT, MaskPixelT, VariancePixelT> &coadd, ///< coadd
    lsst::afw::image::Image<boost::uint16_t> &depthMap, ///< depth map
    lsst::afw::image::MaskedImage<ImagePixelT, MaskPixelT, VariancePixelT> const &image,    ///< image to add to coadd
    MaskPixelT const badPixelMask   ///< skip input pixel if input mask | badPixelMask != 0
) {
    typedef lsst::afw::image::MaskedImage<ImagePixelT, MaskPixelT, VariancePixelT> MaskedImageT;
    typedef lsst::afw::image::Image<boost::uint16_t> DepthMapT;

    if (image.getDimensions() != coadd.getDimensions()) {
        throw LSST_EXCEPT(pexExcept::InvalidParameterException, "MaskedImage sizes do not match");
    }

    // Set the pixels row by row, to avoid repeated checks for end-of-row
    for (int y = 0, endY = image.getHeight(); y != endY; ++y) {
        typename MaskedImageT::const_x_iterator imagePtr = image.row_begin(y);
        typename MaskedImageT::const_x_iterator imageEndPtr = image.row_end(y);
        typename MaskedImageT::x_iterator coaddPtr = coadd.row_begin(y);
        typename DepthMapT::x_iterator depthMapPtr = depthMap.row_begin(y);
        for (; imagePtr != imageEndPtr; ++imagePtr, ++coaddPtr) {
            if ((imagePtr.mask() & badPixelMask) == 0) {
                *coaddPtr += *imagePtr;
                ++(*depthMapPtr);
            }
        }
    }
}

//
// Explicit instantiations
//
#define INSTANTIATE(ImagePixelT) \
    template void lsst::coadd::kaiser::addToCoadd<ImagePixelT, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>( \
        lsst::afw::image::MaskedImage<ImagePixelT, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel> &coadd, \
        lsst::afw::image::Image<boost::uint16_t> &depthMap, \
        lsst::afw::image::MaskedImage<ImagePixelT, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel> const &image, \
        lsst::afw::image::MaskPixel const badPixelMask \
    );

INSTANTIATE(boost::uint16_t);
INSTANTIATE(int);
INSTANTIATE(float);
INSTANTIATE(double);
