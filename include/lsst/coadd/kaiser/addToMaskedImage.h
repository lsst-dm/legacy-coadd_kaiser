// -*- LSST-C++ -*-
#ifndef LSST_COADD_KAISER_COPYMASKEDIMAGE_H
#define LSST_COADD_KAISER_COPYMASKEDIMAGE_H
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

template<typename ImagePixelT, typename MaskPixelT> 
void addToMaskedImage(
    lsst::afw::image::MaskedImage<ImagePixelT, MaskPixelT> &outMaskedImage,
    lsst::afw::image::MaskedImage<ImagePixelT, MaskPixelT> const &inMaskedImage,
    MaskPixelT const badPixelMask
);

}}} // lsst::coadd::kaiser

#endif // !defined(LSST_COADD_KAISER_COPYMASKEDIMAGE_H)
