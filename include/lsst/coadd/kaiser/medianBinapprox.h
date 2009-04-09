// -*- LSST-C++ -*-
#ifndef LSST_COADD_KAISER_MEDIANBINAPPROX_H
#define LSST_COADD_KAISER_MEDIANBINAPPROX_H
/**
* @brief define medianBinapprox
*
* @file
*
* @author Ryan J. Tibshirani, adapted from C to C++ by Russell Owen.
*
* @todo:
* * Remove this once afwMath::Statistics supports vectors and masked images.
*/
#include <iterator>

#include "lsst/afw/image.h"

namespace lsst {
namespace coadd {
namespace kaiser {

    template <class ForwardIterator>
    typename std::iterator_traits<ForwardIterator>::value_type medianBinapprox(
        ForwardIterator first,
        ForwardIterator last,
        int nBins = 1000
    );
    
    template <typename T>
    T medianBinapproxImage(lsst::afw::image::Image<T> const &image, int nBins = 1000);

}}} // lsst::coadd::kaiser

#ifndef SWIG // don't bother SWIG with .cc files
#include "lsst/coadd/kaiser/medianBinapprox.cc"
#endif

#endif // !defined(LSST_COADD_KAISER_MEDIANBINAPPROX_H)
