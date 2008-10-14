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
#include "boost/cstdint.hpp"

#include "lsst/coadd/kaiser.h"
#include "lsst/coadd/kaiser/addToMaskedImage.h"

namespace pexExcept = lsst::pex::exceptions;

/**
* @brief add good pixels from one MaskedImage to another MaskedImage
*
* @todo: move to a better location (afw?)
*
* @throw pexExcept::InvalidParameter if the image sizes do not match.
*/
template <typename ImagePixelT, typename MaskPixelT>
void lsst::coadd::kaiser::addToMaskedImage(
    lsst::afw::image::MaskedImage<ImagePixelT, MaskPixelT> &outMaskedImage, ///< image to be added to
    lsst::afw::image::MaskedImage<ImagePixelT, MaskPixelT> const &inMaskedImage,    ///< image to be added from
    MaskPixelT const badPixelMask   ///< skip input pixel if mask | badPixelMask != 0
) {
//    typedef lsst::afw::image::MaskedPixelAccessor<ImagePixelT, MaskPixelT> MaskedPixelAccessor;

    const unsigned int nCols = inMaskedImage.getCols();
    const unsigned int nRows = inMaskedImage.getRows();
    if ((nCols != outMaskedImage.getCols()) ||
        (nRows != outMaskedImage.getRows())) {
        throw pexExcept::InvalidParameter("MaskedImage sizes do not match");
    }

    lsst::afw::image::MaskedPixelAccessor<ImagePixelT, MaskPixelT> inRowAcc(inMaskedImage);
    lsst::afw::image::MaskedPixelAccessor<ImagePixelT, MaskPixelT> outRowAcc(outMaskedImage);
    for (unsigned int row = 0; row < nRows; ++row, inRowAcc.nextRow(), outRowAcc.nextRow()) {
        lsst::afw::image::MaskedPixelAccessor<ImagePixelT, MaskPixelT> inColAcc(inRowAcc);
        lsst::afw::image::MaskedPixelAccessor<ImagePixelT, MaskPixelT> outColAcc(outRowAcc);
        for (unsigned int col = 0; col < nCols; ++col, inColAcc.nextCol(), outColAcc.nextCol()) {
            if (*(outColAcc.mask) | badPixelMask == 0) {
                *(outColAcc.image) += *(inColAcc.image);
                *(outColAcc.variance) += *(inColAcc.variance);
                *(outColAcc.mask) |= *(inColAcc.mask);
            }
        }
    }    
}

//
// Explicit instantiations
//
template void lsst::coadd::kaiser::addToMaskedImage<double, lsst::afw::image::maskPixelType>(
    lsst::afw::image::MaskedImage<double, lsst::afw::image::maskPixelType> &outMaskedImage,
    lsst::afw::image::MaskedImage<double, lsst::afw::image::maskPixelType> const &inMaskedImage,
    lsst::afw::image::maskPixelType const badPixelMask
);
template void lsst::coadd::kaiser::addToMaskedImage<float, lsst::afw::image::maskPixelType>(
    lsst::afw::image::MaskedImage<float, lsst::afw::image::maskPixelType> &outMaskedImage,
    lsst::afw::image::MaskedImage<float, lsst::afw::image::maskPixelType> const &inMaskedImage,
    lsst::afw::image::maskPixelType const badPixelMask
);
template void lsst::coadd::kaiser::addToMaskedImage<int, lsst::afw::image::maskPixelType>(
    lsst::afw::image::MaskedImage<int, lsst::afw::image::maskPixelType> &outMaskedImage,
    lsst::afw::image::MaskedImage<int, lsst::afw::image::maskPixelType> const &inMaskedImage,
    lsst::afw::image::maskPixelType const badPixelMask
);
template void lsst::coadd::kaiser::addToMaskedImage<boost::uint16_t, lsst::afw::image::maskPixelType>(
    lsst::afw::image::MaskedImage<boost::uint16_t, lsst::afw::image::maskPixelType> &outMaskedImage,
    lsst::afw::image::MaskedImage<boost::uint16_t, lsst::afw::image::maskPixelType> const &inMaskedImage,
    lsst::afw::image::maskPixelType const badPixelMask
);
