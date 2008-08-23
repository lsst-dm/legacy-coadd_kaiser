/**
* @brief Define medianBinapprox
* @file
*
* @author Ryan J. Tibshirani, adapted from C to C++ by Russell Owen.
*/
#include <cmath>
#include <valarray>
#include <stdexcept>

#include "lsst/coadd/kaiser/medianBinapprox.h"

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
* @throw range_error if last <= first or nBins < 2
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
        throw std::range_error("last <= first");
    } else if (nBins < 2) {
        throw std::range_error("nBins < 2");
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
    double leftend =  mu-sigma;
    double rightend = mu+sigma;
    int bin;

    for (ForwardIterator it = first; it != last; ++it) {
        double val = static_cast<double>(*it);
        if (val < leftend) {
            bottomcount++;
        } else if (val < rightend) {
            bin = static_cast<int>((val -leftend) * scalefactor);
            bincounts[bin]++;
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
} 
