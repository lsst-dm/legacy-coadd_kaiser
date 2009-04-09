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

%include "lsst/p_lsstSwig.i"
%import  "lsst/afw/image/imageLib.i" 
%import  "lsst/afw/math/mathLib.i" 

%lsst_exceptions()

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

%include "lsst/coadd/kaiser/addToMaskedImage.h"
// %template(addToMaskedImage) lsst::coadd::kaiser::addToMaskedImage<double, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>;
%template(addToMaskedImage) lsst::coadd::kaiser::addToMaskedImage<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>;
// %template(addToMaskedImage) lsst::coadd::kaiser::addToMaskedImage<int, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>;
// %template(addToMaskedImage) lsst::coadd::kaiser::addToMaskedImage<boost::uint16_t, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>;

%include "lsst/coadd/kaiser/CoaddComponent.h"
