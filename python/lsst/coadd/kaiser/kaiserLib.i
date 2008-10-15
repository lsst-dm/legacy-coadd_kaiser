// -*- lsst-c++ -*-
%define kaiserLib_DOCSTRING
"
Python interface to lsst::coadd::kaiser functions and classes
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.coadd.kaiser",docstring=kaiserLib_DOCSTRING) kaiserLib

// Suppress swig complaints; see afw/image/imageLib.i for more
#pragma SWIG nowarn=362                 // operator=  ignored

// Everything we will need in the _wrap.cc file
// the lsst/afw files are actually included by kaiser.h but I'd rather be explicit
%{
#include "boost/cstdint.hpp"
#include "lsst/afw/image.h"
#include "lsst/afw/math.h"
#include "lsst/coadd/kaiser.h"
%}

%init %{
%}

// Everything whose bindings we will have to know about
%import "lsst/daf/data/LsstBase.h"  // avoid warning: Nothing known about base class 'lsst::daf::data::LsstBase'
%import "lsst/afw/image/Mask.h" // needed so SWIG knows lsst::afw::image::maskPixelType = boost::uint16_t
%include "lsst/p_lsstSwig.i"    // this needs to go first otherwise I do not know about e.g. boost
%include "lsst/afw/image/lsstImageTypes.i"  // vw and Image/Mask types and typedefs

// handle C++ arguments that should be outputs in python
%apply int& OUTPUT { int& };
%apply float& OUTPUT { float& };
%apply double& OUTPUT { double& };

%pythoncode %{
def version(HeadURL = r"$HeadURL: svn+ssh://svn.lsstcorp.org/DMS/afw/trunk/python/lsst/coadd/kaiser/kaiserLib.i $"):
    """Return a version given a HeadURL string. If a different version is setup, return that too"""

    """Return a version given a HeadURL string"""
    return guessSvnVersion(HeadURL)

%}

%ignore lsst::coadd::kaiser::medianBinapprox;
%include "lsst/coadd/kaiser/medianBinapprox.h"
%template(medianBinapproxImage)  lsst::coadd::kaiser::medianBinapproxImage<float>;

/*
// If both medianBinapprox functions have the same name then the following fails with:
//   /usr/include/c++/4.0.0/bits/stl_iterator_base_types.h:129: error: 'float' is not a class, struct, or union type
//   /usr/include/c++/4.0.0/bits/stl_iterator_base_types.h:130: error: 'float' is not a class, struct, or union type
%include "lsst/coadd/kaiser/medianBinapprox.h"
%template(medianBinapprox)  lsst::coadd::kaiser::medianBinapprox<float>;

// and this alternate approach fails in exactly the same way:
namespace lsst {
namespace coadd {
namespace kaiser {
    template <typename T>
    T medianBinapprox(
        lsst::afw::image::Image<T> const &image,
        int nBins = 1000
    );
}}} // lsst::coadd::kaiser
%template(medianBinapprox)  lsst::coadd::kaiser::medianBinapprox<float>;
*/

%include "lsst/coadd/kaiser/addToMaskedImage.h"
%template(addToMaskedImage) lsst::coadd::kaiser::addToMaskedImage<double, lsst::afw::image::maskPixelType>;
%template(addToMaskedImage) lsst::coadd::kaiser::addToMaskedImage<float, lsst::afw::image::maskPixelType>;
%template(addToMaskedImage) lsst::coadd::kaiser::addToMaskedImage<int, lsst::afw::image::maskPixelType>;
%template(addToMaskedImage) lsst::coadd::kaiser::addToMaskedImage<boost::uint16_t, lsst::afw::image::maskPixelType>;

%include "lsst/coadd/kaiser/CoaddComponent.h"
