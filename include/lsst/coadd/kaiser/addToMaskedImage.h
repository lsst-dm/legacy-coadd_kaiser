// -*- LSST-C++ -*-
#ifndef LSST_COADD_KAISER_ADDTOMASKEDIMAGE_H
#define LSST_COADD_KAISER_ADDTOMASKEDIMAGE_H
/**
* @brief declare addToMaskedImage
*
* @file
*
* @author Russell Owen.
*
* @todo: move to a better location
*/
#include "lsst/afw/image.h"

namespace lsst {
namespace coadd {
namespace kaiser {

    template<typename ImagePixelT, typename MaskPixelT, typename VariancePixelT>
    void addToMaskedImage(
        lsst::afw::image::MaskedImage<ImagePixelT, MaskPixelT, VariancePixelT> &outMaskedImage,
        lsst::afw::image::MaskedImage<ImagePixelT, MaskPixelT, VariancePixelT> const &inMaskedImage,
        MaskPixelT const badPixelMask
    );

}}} // lsst::coadd::kaiser

#endif // !defined(LSST_COADD_KAISER_ADDTOMASKEDIMAGE_H)
