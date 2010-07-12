// -*- LSST-C++ -*-

/* 
 * LSST Data Management System
 * Copyright 2008, 2009, 2010 LSST Corporation.
 * 
 * This product includes software developed by the
 * LSST Project (http://www.lsst.org/).
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the LSST License Statement and 
 * the GNU General Public License along with this program.  If not, 
 * see <http://www.lsstcorp.org/LegalNotices/>.
 */
 
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
