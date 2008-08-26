#include <iostream>

#include "lsst/pex/logging/Trace.h"
#include "lsst/afw/image.h"
#include "lsst/coadd/kaiser.h"

int main(int argc, char **argv) {
    typedef float pixelType;
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
        lsst::afw::image::MaskedImage<pixelType, lsst::afw::image::maskPixelType> mImage;
        mImage.readFits(argv[1]);
        
        // create psf kernel
        double sigma = fwhm / 2.35;
        int kSize = static_cast<int>(2.0 * fwhm);
        lsst::coadd::kaiser::DoubleGaussianFunction2<double> psfFunc(sigma, sigma*10.0, 0.1);
        lsst::afw::math::AnalyticKernel psfKernel(psfFunc, kSize, kSize);
        
        lsst::coadd::kaiser::CoaddComponent(mImage, psfKernel);
    }

     //
     // Check for memory leaks
     //
     if (lsst::daf::base::Citizen::census(0) != 0) {
         std::cerr << "Leaked memory blocks:" << std::endl;
         lsst::daf::base::Citizen::census(std::cerr);
     }
    
}
