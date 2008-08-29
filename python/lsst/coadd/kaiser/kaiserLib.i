// -*- lsst-c++ -*-
%define kaiserLib_DOCSTRING
"
Python interface to lsst::coadd::kaiser functions and classes
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.afw.math",docstring=kaiserLib_DOCSTRING) kaiserLib

// Suppress swig complaints
// copied from afw/image/imageLib.i
//#pragma SWIG nowarn=314                 // print is a python keyword (--> _print)
//#pragma SWIG nowarn=317                 // specialization of non-template
#pragma SWIG nowarn=362                 // operator=  ignored
//#pragma SWIG nowarn=389                 // operator[] ignored
//#pragma SWIG nowarn=503                 // Can't wrap 'operator unspecified_bool_type'

// Everything we will need in the _wrap.cc file
%{
#include "lsst/coadd/kaiser.h"
%}

%init %{
%}

//namespace boost {
//    class bad_any_cast; // remove warning: Nothing known about 'boost::bad_any_cast'
//}

// Everything whose bindings we will have to know about
%include "lsst/p_lsstSwig.i"    // this needs to go first otherwise i do not know about e.g. boost
%include "lsst/afw/image/lsstImageTypes.i"  // vw and Image/Mask types and typedefs

// handle C++ arguments that should be outputs in python
%apply int& OUTPUT { int& };
%apply float& OUTPUT { float& };
%apply double& OUTPUT { double& };

%pythoncode %{
#import lsst.daf.data
import lsst.utils

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

%include "lsst/coadd/kaiser/CoaddComponent.h"
