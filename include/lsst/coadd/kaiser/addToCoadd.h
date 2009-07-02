// -*- LSST-C++ -*-
#ifndef LSST_COADD_KAISER_ADDTOCOADD_H
#define LSST_COADD_KAISER_ADDTOCOADD_H
/**
* @brief declare addToCoadd
*
* @file
*
* @author Russell Owen.
*
* @todo: move to a better location
*/
#include "boost/cstdint.hpp"

#include "lsst/afw/image.h"

namespace lsst {
namespace coadd {
namespace kaiser {

    template<typename ImagePixelT, typename MaskPixelT, typename VariancePixelT>
    void addToCoadd(
        lsst::afw::image::MaskedImage<ImagePixelT, MaskPixelT, VariancePixelT> &coadd,
        lsst::afw::image::Image<boost::uint16_t> &depthMap,
        lsst::afw::image::MaskedImage<ImagePixelT, MaskPixelT, VariancePixelT> const &image,
        MaskPixelT const badPixelMask
    );

}}} // lsst::coadd::kaiser

#endif // !defined(LSST_COADD_KAISER_ADDTOCOADD_H)
