#include <iostream>

#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/image.h"
#include "lsst/afw/math.h"
#include "lsst/coadd/kaiser.h"

int main(int argc, char **argv) {
    typedef float Pixel;
    const double DefFwhm = 3.0;
    
    lsst::pex::logging::Trace::setDestination(std::cout);
    lsst::pex::logging::Trace::setVerbosity("lsst.coadd", 5);

    if (argc < 2) {
        std::cerr << "Usage: coaddComponent fitsFile [sigma]]" << std::endl;
        std::cerr << "fitsFile excludes the \"_img.fits\" suffix" << std::endl;
        return 1;
    }

    { // block in which to allocate and deallocate memory

        double fwhm = DefFwhm;
        if (argc > 2) {
            std::istringstream(argv[2]) >> fwhm;
        }
        
        // read in fits file
        lsst::afw::image::Exposure<Pixel> scienceExposure(argv[1]);
        
        // create psf kernel
        double sigma = fwhm / 2.35;
        lsst::afw::math::DoubleGaussianFunction2<double> psfFunc(sigma, sigma*10.0, 0.1);
        int kSize = 2 * static_cast<int>(fwhm + 0.5) + 1;
        lsst::afw::math::AnalyticKernel psfKernel(kSize, kSize, psfFunc);
        lsst::coadd::kaiser::CoaddComponent(scienceExposure, psfKernel);
        
        std::cout << "Got CoaddComponent" << std::endl;
    }

    //
    // Check for memory leaks
    //
    if (lsst::daf::base::Citizen::census(0) != 0) {
        std::cerr << "Leaked memory blocks:" << std::endl;
        lsst::daf::base::Citizen::census(std::cerr);
    }
}
