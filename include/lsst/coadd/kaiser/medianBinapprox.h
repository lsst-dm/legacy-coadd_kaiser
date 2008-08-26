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
* - get permission to use publicly (pending)
* - if this algorithm is of general interest them move to lsst::afw::math
* - check if RHL has a better version of this algorithm since it sounds similar
*   to something he's done in the past
* 
*/
#include <iterator>

#include "lsst/afw/image.h"

namespace lsst {
namespace coadd {
namespace kaiser {

#ifndef SWIG // don't wrap STL version in Python
template <class ForwardIterator>
typename std::iterator_traits<ForwardIterator>::value_type medianBinapprox(
    ForwardIterator first,
    ForwardIterator last,
    int nBins = 1000
);
#endif

template <typename T>
T medianBinapprox(
    lsst::afw::image::Image<T> const &image,
    int nBins = 1000
);

}}} // lsst::coadd::kaiser

#ifndef SWIG // don't bother SWIG with .cc files
#include "lsst/coadd/kaiser/medianBinapprox.cc"
#endif

#endif // !defined(LSST_COADD_KAISER_MEDIANBINAPPROX_H)
