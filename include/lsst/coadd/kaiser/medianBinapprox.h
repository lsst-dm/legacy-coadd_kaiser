/**
* @brief define medianBinapprox
*
* @file
*
* @author Ryan J. Tibshirani, adapted from C to C++ by Russell Owen.
*/
#include <iterator>

namespace lsst {
namespace coadd {
namespace kaiser {

template <class ForwardIterator>
typename std::iterator_traits<ForwardIterator>::value_type medianBinapprox(
    ForwardIterator first,
    ForwardIterator last,
    int nBins = 1000
);

}}} // lsst::coadd::kaiser