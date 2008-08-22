/**
* @brief define medianBinapprox
*
* @file
*
* @author Ryan J. Tibshirani, adapted from C to C++ by Russell Owen.
*/
#include <valarray>
#include <cmath>

namespace lsst {
namespace coadd {
namespace kaiser {

template <class T>
T medianBinapprox(
    T const * const first,  ///< iterator to first element of array
    T const * const last    ///< iterator to last+1 element of array
);

}}} // lsst::coadd::kaiser