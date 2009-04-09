// -*- lsst-c++ -*-
%define kaiserLib_DOCSTRING
"
Python interface to lsst::coadd::kaiser functions and classes
"
%enddef

%feature("autodoc", "1");
%module(package="lsst.coadd.kaiser", docstring=kaiserLib_DOCSTRING) kaiserLib

// Everything we will need in the _wrap.cc file
%{
#include "boost/cstdint.hpp"
#include "lsst/coadd/kaiser.h"
%}

%include "lsst/p_lsstSwig.i"
%import  "lsst/afw/image/imageLib.i" 
%import  "lsst/afw/math/mathLib.i" 

%lsst_exceptions()

// it is not convenient to call the C++-iterator-based medianBinapprox version from Python
%ignore lsst::coadd::kaiser::medianBinapprox;
%include "lsst/coadd/kaiser/medianBinapprox.h"
%template(medianBinapproxImage)  lsst::coadd::kaiser::medianBinapproxImage<float>;

%include "lsst/coadd/kaiser/addToMaskedImage.h"
%template(addToMaskedImage) lsst::coadd::kaiser::addToMaskedImage<double, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>;
%template(addToMaskedImage) lsst::coadd::kaiser::addToMaskedImage<float, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>;
%template(addToMaskedImage) lsst::coadd::kaiser::addToMaskedImage<int, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>;
%template(addToMaskedImage) lsst::coadd::kaiser::addToMaskedImage<boost::uint16_t, lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>;

%include "lsst/coadd/kaiser/CoaddComponent.h"
