/**
* @brief Define medianBinapprox
* @file
*
* @author Ryan J. Tibshirani, adapted from C to C++ by Russell Owen.
*/
#include <cmath>
#include <valarray>
#include <stdexcept>

#include "medianBinapprox.h"

const int NBins = 1001;

/**
* @brief Compute the median using the binapprox algorithm.
*
* The accuracy is 1/1000th of a standard deviation (since the code uses 1000 bins).
*
* Uses "binapprox", a fast algorithm that does not copy or rearrange the input data.
*
* The binapprox algorithm is described in Ryan J. Tibshirani's paper:
* "Fast computation of the median by successive binning", Jue 23, 2008
* <http://stat.stanford.edu/~ryantibs/median/medianpaper.pdf>
*
* @throw range_error if last <= first
*
* @return approximate median
*/
template <class T>
T lsst::coadd::kaiser::medianBinapprox(
    T const * const first,  ///< iterator to first element of array
    T const * const last    ///< iterator to last+1 element of array
) {
    if (first >= last) {
        throw std::range_error("last <= first");
    }

    // Compute the number of elements (n), mean (mu) and standard deviation (sigma)
    long int n = 0;
    T sum = 0;
    for (T *it = first; it != last; ++it) {
        n++;
        sum += *it;
    }
    T mu = sum/static_cast<T>(n);

    sum = 0;
    for (T *it = first; it != last; ++it) {
        sum += ((*it)-mu)*((*it)-mu);
    }
    T sigma = std::sqrt(sum/static_cast<T>(n));

    // Bin data across the interval [mu-sigma, mu+sigma]
    int bottomcount = 0;
    std::valarray<int> bincounts(NBins);

    T scalefactor = (NBins-1)/(2*sigma);
    T leftend =  mu-sigma;
    T rightend = mu+sigma;
    int bin;

    for (T *it = first; it != last; ++it) {
        if ((*it) < leftend) {
            bottomcount++;
        }
        else if ((*it) < rightend) {
            bin = static_cast<int>(((*it)-leftend) * scalefactor);
            bincounts[bin]++;
        }
    }

    // Find the bin that contains the median
    if (n & 1) {
        // n is odd
        long int k = (n+1)/2;
        long int count = bottomcount;

        for (int i = 0; i < NBins; i++) {
            count += bincounts[i];

            if (count >= k) {
                return (i+0.5)/scalefactor + leftend;
            }
        }
    } else {
        // n is even
        long int k = n/2;
        long int count = bottomcount;
        
        for (int i = 0; i < NBins; i++) {
            count += bincounts[i];
            
            if (count >= k) {
                int j = i;
                while (count == k) {
                    j++;
                    count += bincounts[j];
                }
                return static_cast<T>(i+j+1)/(2*scalefactor) + leftend;
            }
        }
    } 
} 
