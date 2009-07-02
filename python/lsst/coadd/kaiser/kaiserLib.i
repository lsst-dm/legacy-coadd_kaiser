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

%include "lsst/coadd/kaiser/addToCoadd.h"
%define %ADDTOMASKEDIMAGE(TYPE)
    %template(addToCoadd) lsst::coadd::kaiser::addToCoadd<TYPE,
        lsst::afw::image::MaskPixel, lsst::afw::image::VariancePixel>;
%enddef
%ADDTOMASKEDIMAGE(double);
%ADDTOMASKEDIMAGE(float);
%ADDTOMASKEDIMAGE(int);
%ADDTOMASKEDIMAGE(boost::uint16_t);

SWIG_SHARED_PTR_DERIVED(CoaddComponent, lsst::daf::base::Citizen, lsst::coadd::kaiser::CoaddComponent)
%include "lsst/coadd/kaiser/CoaddComponent.h"
