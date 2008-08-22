// -*- lsst-c++ -*-
%{
#include "lsst/coadd/kaiser/medianBinapprox.h"
%}

%include "lsst/coadd/kaiser/medianBinapprox.h"

%template(medianBinapprox)  lsst::coadd::kaiser::medianBinapprox<int>;
%template(medianBinapprox)  lsst::coadd::kaiser::medianBinapprox<float>;
%template(medianBinapprox)  lsst::coadd::kaiser::medianBinapprox<double>;
