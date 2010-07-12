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
 
/**
* @brief Define medianBinapprox
*
* @file
*
* @author Ryan J. Tibshirani, adapted from C to C++ by Russell Owen.
*/
#include <cmath>
#include <valarray>
#include <stdexcept>

#include "lsst/pex/exceptions.h"

namespace pexExcept = lsst::pex::exceptions;

/**
* @brief Compute the median using the binapprox algorithm.
*
* The accuracy is 1/nBins of a standard deviation.
*
* Uses "binapprox", a fast algorithm that does not copy or rearrange the input data.
* The binapprox algorithm is described in Ryan J. Tibshirani's paper:
* "Fast computation of the median by successive binning", Jue 23, 2008
* <http://stat.stanford.edu/~ryantibs/median/medianpaper.pdf>
*
* @throw pexExcept::RangeErrorException if last <= first or nBins < 2
*
* @return approximate median
*/
template <class ForwardIterator>
typename std::iterator_traits<ForwardIterator>::value_type lsst::coadd::kaiser::medianBinapprox(
    ForwardIterator first,  ///< iterator to first element of array
    ForwardIterator last,   ///< iterator to last+1 element of array
    int nBins               ///< number of bins to use; 1000 is a typical value
) {
    if (first >= last) {
        throw LSST_EXCEPT(pexExcept::RangeErrorException, "last <= first");
    } else if (nBins < 2) {
        throw LSST_EXCEPT(pexExcept::RangeErrorException, "nBins < 2");
    }

    typedef typename std::iterator_traits<ForwardIterator>::value_type ValueType;

    // Compute the number of elements (n), mean (mu) and standard deviation (sigma)
    long int n = 0;
    double sum = 0;
    for (ForwardIterator it = first; it != last; ++it) {
        n++;
        sum += static_cast<double>(*it);
    }
    double mu = sum/static_cast<double>(n);
    if (n < 3) {
        return static_cast<ValueType>(mu);
    }

    double sumSq = 0;
    for (ForwardIterator it = first; it != last; ++it) {
        double val = static_cast<double>(*it) - mu;
        sumSq += val * val;
    }
    double sigma = std::sqrt(sumSq/static_cast<double>(n));

    // Bin data across the interval [mu-sigma, mu+sigma]
    int bottomcount = 0;
    std::valarray<int> bincounts(nBins);

    double scalefactor = static_cast<double>(nBins - 1)/(2.0 * sigma);
    if (std::isinf(sigma)) {
        // data are too closely spaced, just return mean
        return static_cast<ValueType>(mu);
    }
    double leftend =  mu-sigma;
    int bin;

    for (ForwardIterator it = first; it != last; ++it) {
        double val = static_cast<double>(*it);
        if (val < leftend) {
            bottomcount++;
        } else {
            bin = static_cast<int>((val -leftend) * scalefactor);
            if (bin < nBins) {
                bincounts[bin]++;
            }
        }
    }

    // Find the bin that contains the median
    if (n & 1) {
        // n is odd
        long int k = (n+1)/2;
        long int count = bottomcount;

        for (int i = 0; i < nBins; i++) {
            count += bincounts[i];

            if (count >= k) {
                return static_cast<ValueType>((static_cast<double>(i) + 0.5)/scalefactor + leftend);
            }
        }
    } else {
        // n is even
        long int k = n / 2;
        long int count = bottomcount;
        
        for (int i = 0; i < nBins; i++) {
            count += bincounts[i];
            
            if (count >= k) {
                int j = i;
                while (count == k) {
                    j++;
                    count += bincounts[j];
                }
                return static_cast<ValueType>(static_cast<double>(i + j + 1)/(2.0 * scalefactor) + leftend);
            }
        }
    }
    throw LSST_EXCEPT(pexExcept::RuntimeErrorException, "Unexpectedly failed to return a value");
} 


/**
* @brief Compute the median of an lsst::afw::image::Image using the binapprox algorithm.
*
* The accuracy is 1/nBins of a standard deviation.
*
* Uses "binapprox", a fast algorithm that does not copy or rearrange the input data.
* The binapprox algorithm is described in Ryan J. Tibshirani's paper:
* "Fast computation of the median by successive binning", Jue 23, 2008
* <http://stat.stanford.edu/~ryantibs/median/medianpaper.pdf>
*
* @note the name is different than medianBinapprox to avoid confusing SWIG;
* See the comments in kaiserLib.i for details.
*
* @throw pexExcept::RangeErrorException if no pixels or nBins < 2
*
* @return approximate median
*/
template <typename T>
T lsst::coadd::kaiser::medianBinapproxImage(
    lsst::afw::image::Image<T> const &image,   ///< image for which to compute median
    int nBins       ///< number of bins to use; 1000 is a typical value
) {
    return medianBinapprox(image.begin(), image.end(), nBins);
}
