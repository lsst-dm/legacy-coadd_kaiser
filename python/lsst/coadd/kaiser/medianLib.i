// -*- lsst-c++ -*-
%{
#include "lsst/coadd/kaiser/medianBinapprox.h"
%}

%include "lsst/afw/image/Image.h"
%include "lsst/coadd/kaiser/medianBinapprox.h"

/* the following fails with:
/usr/include/c++/4.0.0/bits/stl_iterator_base_types.h: In instantiation of 'std::iterator_traits<float>':
python/lsst/coadd/kaiser/kaiserLib_wrap.cc:14536:   instantiated from here
/usr/include/c++/4.0.0/bits/stl_iterator_base_types.h:129: error: 'float' is not a class, struct, or union type
/usr/include/c++/4.0.0/bits/stl_iterator_base_types.h:130: error: 'float' is not a class, struct, or union type
/usr/include/c++/4.0.0/bits/stl_iterator_base_types.h:131: error: 'float' is not a class, struct, or union type
/usr/include/c++/4.0.0/bits/stl_iterator_base_types.h:132: error: 'float' is not a class, struct, or union type
/usr/include/c++/4.0.0/bits/stl_iterator_base_types.h:133: error: 'float' is not a class, struct, or union type
*/
// %template(medianBinapprox)  lsst::coadd::kaiser::medianBinapprox<float>;
